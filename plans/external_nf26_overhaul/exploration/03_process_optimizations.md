# 03 - Portable Process Optimizations (INTERNAL -> EXTERNAL)

Scope: performance/efficiency improvements found in INTERNAL (`/nssd2/humebc/internal/rnaseq`)
that can be applied to EXTERNAL (`/nssd2/humebc/external/csgenetics_scrnaseq`) processes
WITHOUT changing what EXTERNAL computes or its published outputs.

Conventions:
- "Output-equiv risk = NO" means I assess the change cannot alter EXTERNAL outputs (same tool, same args, only threading/streaming/temp-file changes).
- "FLAG" means the change touches algorithm/tool/args and could alter outputs; only adopt with snapshot/regression testing or do not adopt.
- External process locations are `modules/processes.nf` line numbers. Internal are `modules/<mod>/main.nf`.

---

## A. Drop-in optimizations (Output-equiv risk = NO) -- recommended

| External process (file:line) | Internal counterpart | Specific optimization | Expected effect | Output-equiv risk |
|---|---|---|---|---|
| `star` (processes.nf:191) `--runThreadN 8` | `alignment/main.nf:8` `--runThreadN ${task.cpus}` | Hardcoded `8` while `conf/base.config:51` allocates `cpus = 16`. Use `--runThreadN ${task.cpus}`. STAR args otherwise identical. | ~2x faster alignment by using the cores already reserved (currently 8 idle). | NO (STAR output BAM is independent of thread count). |
| `initial_feature_count` (processes.nf:310) `samtools sort $bam` | `alignment/main.nf:29` `samtools sort -@ ${task.cpus} -m 2G` | External pre-sorts the STAR bam with single-threaded default `samtools sort`. Add `-@ ${task.cpus} -m <N>G`. Process has `cpus = 4` (base.config:67). | Faster coordinate sort of the largest BAM in the run. | NO (sort order deterministic; threading does not change it). |
| `sort_index_bam` (processes.nf:586-602) | `alignment/main.nf:29` pattern | `samtools sort` and `samtools index` both run single-threaded; `cpus = 1` (base.config:83). Bump cpus and add `-@`. `samtools view -c -f 16` and `index` also accept `-@`. | Faster sort+index of the annotated BAM. | NO. |
| `multimapper_transcript_assignment` (processes.nf:410) `samtools sort -n $bam \| ... gawk` | `gene_annotation/main.nf:213` `samtools sort -n -@ ${task.cpus} ... \| bam-assigner` | The `samtools sort -n` here is single-threaded (`cpus = 1`, base.config:168). Adding threads needs a cpus bump but the streamed pipe is already temp-file-free. | Faster name-sort feeding gawk. | NO (gawk logic unchanged; sort order with `-n` is deterministic). |
| `multimapper_exon_assignment` (processes.nf:450-451) | `gene_annotation/main.nf:239-245` | `featureCounts -T 4` but `cpus = 1` (base.config:172) -> over-subscribed/false threading. Align cpus to `-T`. `samtools sort -n` also single-threaded. | Correct core allocation; faster. | NO. |
| `dedup` (processes.nf:625) `umi_tools dedup --in-sam` | (internal replaced umi_tools entirely; see Section C) | `cpus = 4` (base.config:87) but `umi_tools dedup` is single-threaded -> 3 cores wasted. Right-size cpus to 1 (frees scheduler capacity, no speed change). Keep `--in-sam`/`--per-cell` exactly. | Better cluster packing; no output change. | NO. |

Notes on threading-only changes: `samtools` `-@` and STAR `--runThreadN` are documented to be
output-deterministic; only wall-time changes. featureCounts `-T` likewise does not change assignments.

---

## B. Streaming / temp-file / container improvements (Output-equiv risk = NO unless noted)

| External process (file:line) | Internal counterpart | Optimization | Expected effect | Output-equiv risk |
|---|---|---|---|---|
| `initial_feature_count` (processes.nf:307-310) | `gene_annotation/main.nf:16-18` | INTERNAL does NOT pre-sort before featureCounts; it runs featureCounts directly on the STAR Unsorted bam (its `star` process already produced a sorted bam but feeds the sorted one). EXTERNAL adds an extra `samtools sort` (line 307) purely to produce the `*_Aligned.sortedByCoord.out.bam` artifact published to `featureCounts/`. featureCounts does not require coordinate sort. **If the sorted BAM is a required published output, keep the sort; if it is only an intermediate, the sort can be dropped.** | Removing a full coordinate sort of the largest BAM if the sorted artifact is not consumed downstream. | FLAG: must confirm nothing downstream/no customer deliverable depends on `*_Aligned.sortedByCoord.out.bam`. featureCounts assignment output is identical regardless of input sort order. |
| `star` (processes.nf:209) `--outReadsUnmapped Fastx` | `alignment/main.nf` (no `--outReadsUnmapped`) | INTERNAL omits `--outReadsUnmapped Fastx`. EXTERNAL writes unmapped reads to disk. **Only drop if the unmapped FASTX files are not consumed/published by EXTERNAL.** | Less disk I/O during alignment. | FLAG: verify external workflow/QC does not use `*_Unmapped.out.mate1`. If used, do not change. |
| `io_count` (processes.nf:661) `samtools view \| awk` | n/a (internal uses Rust `count_table_builder`) | Pure stream already; `cpus = 4` (base.config:91) but `samtools view` single-threaded and awk single-threaded -> over-allocated. Add `-@` to `samtools view` or right-size cpus. | Better packing; minor speed. | NO (awk parse unchanged). |
| `count_matrix` container `scanpy_anndata:0.0.4` (images.config:32) | internal uses `scanpy_anndata:0.0.7` | Newer scanpy/anndata image (0.0.7) used internally; potentially faster anndata write. | Possible speedup, newer deps. | FLAG: must verify h5ad/mtx byte/format equivalence; library version bump can change outputs. Do NOT adopt without regression diff. |
| `dedup` container `umi-tools-csgx:0.9` (images.config:12) | internal `umi-tools-csgx:1.0` | Internal pins umi-tools image `1.0` vs external `0.9`. | Possible bugfixes/speed. | FLAG: umi_tools version change can alter dedup results. Keep `0.9` unless output-verified. |

---

## C. Algorithmic / process-fusion improvements (HIGH VALUE but FLAG for output verification)

INTERNAL has been substantially re-architected with Rust binaries that fuse several
EXTERNAL processes. These are the largest potential wins but each must be validated to
produce byte/semantically identical EXTERNAL outputs before adoption.

| External processes (file:line) | Internal counterpart | Optimization | Expected effect | Output-equiv risk |
|---|---|---|---|---|
| `filter_for_UMRs_mismatch` (319), `filter_for_multimappers_mismatch` (377), `count_high_conf_annotated_umr_multimap` (512) + several `samtools view -c` | `gene_annotation/main.nf:28` `split_and_count_feature_bam` -> Rust `bam-splitter` (`bam-tools:0.2`) | Single Rust pass over the featureCounts BAM splits UMR/multimapper, applies the `[NH]`/`[nM]` mismatch filters, and emits all count metrics at once -- replacing multiple `samtools view -e '...'` filter passes and separate counting processes. | Large: one BAM pass instead of N; eliminates intermediate BAMs and process-launch overhead. | FLAG: must prove the Rust filter expressions exactly reproduce EXTERNAL's `[NH]==1 && ([nM]==0..3)` etc. and that emitted counts match EXTERNAL semantics. Note: EXTERNAL counts/metrics differ from INTERNAL (mixed-species split differs), so the binary likely needs EXTERNAL-specific behavior. Do not adopt blindly. |
| `multimapper_transcript_assignment` (392) + `multimapper_exon_assignment` (435) gawk script | `gene_annotation/main.nf:201` `bam-assigner` (`bam-tools:0.2`) | Rust `bam-assigner` replaces the `gawk -f assign_multi_mappers.gawk` two-pass SAM-body approach, streamed from `samtools sort -n`. | Faster multimapper assignment, no `.sam_body` temp files. | FLAG: gawk logic is the assignment algorithm; Rust reimplementation must be proven output-identical to EXTERNAL gawk. EXTERNAL gawk script differs from INTERNAL phase5 gawk. High risk of output change. |
| `dedup` (607) + `io_count` (642) + `count_matrix` (674) | `make_counts/main.nf:11` `make_raw_count_matrices` -> Rust `count_table_builder` (`count-table-builder:0.2`) | INTERNAL fuses UMI grouping (replacing `umi_tools dedup`), the `bcGeneSummary` extraction (`io_count`), and count-matrix construction into one Rust binary, deriving dedup stats analytically (no `group.bam`/`group.tsv`/`dedup.log` files). | Very large: removes the umi_tools dedup bottleneck and two downstream passes. | FLAG: this CHANGES the dedup engine (umi_tools -> custom Rust grouping) and the count path. Almost certainly alters dedup edges and h5ad contents. EXTERNAL must keep umi_tools-equivalent output. Do NOT adopt as-is; only a like-for-like umi_tools acceleration would be output-safe. |

---

## D. Things that are NOT portable / would change outputs (do not adopt)

- INTERNAL `qc` (`qc/main.nf:155`) adds `--sss-trim ${sss_mer_length}` and extra filtered-barcode/QC outputs and emits cascade-format counts. EXTERNAL `qc` (processes.nf:158) runs `--no-filtered` and no `--sss-trim`. Both use container `qc:0.1`, so the binary is identical, but the **arguments differ and would change trimming/QC outputs.** Keep EXTERNAL invocation as-is.
- INTERNAL `categorize_reads` takes `--annotated_h5ad` + `--cell_caller_method`; EXTERNAL takes `--raw_count_matrix_h5ad` + `--fastp_json`. Different inputs/semantics -> not portable.
- INTERNAL `make_raw_count_matrices` supports downsampling and a different metrics set than EXTERNAL `count_matrix`. Not output-equivalent.
- INTERNAL `features_file` adds a `sed 's/mm10___/mm10_/g'` step (resource_preparation/main.nf:19) that EXTERNAL omits. This CHANGES the features file content -> do not port.

---

## E. Summary of safe, recommended actions (no output change)

1. STAR: `--runThreadN 8` -> `--runThreadN ${task.cpus}` (currently 8 cores idle of 16). [processes.nf:205]
2. Add `-@ ${task.cpus}` (and `-m`) to every `samtools sort` (initial_feature_count:307, sort_index_bam:599, multimapper_*:410/451) and align cpus to actual threads. [output-deterministic]
3. Add `-@` to heavy `samtools view`/`index` calls; right-size cpus where tools are single-threaded (`dedup` cpus=4 but umi_tools is single-threaded; `io_count` cpus=4 but samtools/awk single-threaded).
4. Investigate dropping the redundant `samtools sort` in `initial_feature_count` and `--outReadsUnmapped Fastx` in `star` -- only if those artifacts are confirmed unused (Section B, FLAGGED).

## F. High-value but MUST-verify (Section C)

The Rust-binary process fusions (`bam-splitter`, `bam-assigner`, `count_table_builder`) are the
biggest speedups but each changes the engine/algorithm and INTERNAL's outputs differ from
EXTERNAL's by design. Treat as a separate workstream requiring exact output-equivalence proof
(or EXTERNAL-specific binary builds) before adoption.
