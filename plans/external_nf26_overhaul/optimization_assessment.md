# External pipeline — process optimization assessment

Source: full read of `modules/processes.nf` (905 lines, 33 processes) + internal counterparts
(`internal/rnaseq/modules/{alignment,gene_annotation,make_counts,qc,...}`). Cross-ref exploration
doc `exploration/03_process_optimizations.md`.

Every rewrite below is classified by **output-risk**:
- **R0** = output cannot change (threading/cpu/IO only; tool + args identical). Verify by determinism probe.
- **R1** = structurally different commands, same tool semantics; output expected identical but MUST be proven per-stage.
- **R2** = different engine/algorithm; output WILL or MAY change; gated, may require a blessed re-baseline.

## Runtime hotspots (where the wall-clock actually goes, real data)
Ranked from the DAG + tool characteristics (to be confirmed by the Phase 0 per-process trace):
1. `star` — alignment, CPU-dominant. **Currently `--runThreadN 8` while `cpus=16` → ~half the box idle.**
2. `dedup` (`umi_tools dedup`) — single-threaded, scales with aligned reads; the single biggest per-process sink on large samples.
3. UMR + multimapper annotation branch (`initial_feature_count` → `filter_*` → `umr/multimapper_*_assignment` → merges) — several `featureCounts`, `samtools sort -n`, and `gawk` passes over large BAMs, largely serial.
4. `initial_feature_count` — `samtools sort` (single-threaded) of the largest BAM + `featureCounts -T 4`.
5. `sort_index_bam` — single-threaded `samtools sort` then separate `samtools index`.
6. `count_matrix` / `filter_count_matrix` / `cell_caller` — scanpy/anndata Python, moderate.
7. `single_sample_multiqc` / `multi_sample_multiqc` — moderate, little headroom.
8. Everything else (downloads, merges, `io_count`, `qc`, RSeQC, report render) — IO-bound or already fast.

## Tier 1 — R0 mechanical wins (do unconditionally, guaranteed identical output)

| Process (line) | Change | Effect |
|---|---|---|
| ~~`star` (205)~~ — MOVED TO R1, see below | `--runThreadN 8` → `--runThreadN ${task.cpus}` | ~2x alignment | **NOT R0.** STAR record ORDER is thread-dependent, and that order feeds the order-sensitive `gawk` multimapper assignment (`assign_multi_mappers.gawk:72,94`, last-write-wins on the gene tag) and `umi_tools dedup` representative selection → can change published counts. Gate behind the order-perturbation test (regression §5). |
| `initial_feature_count` (307,310) | `samtools sort` → `samtools sort -@ ${task.cpus} -m 2G`; `featureCounts -T 4` → `-T ${task.cpus}`; align `cpus` in base.config | Faster sort of largest BAM + correctly-threaded featureCounts. |
| `umr_exon_assignment` (372), `multimapper_exon_assignment` (450) | `featureCounts -T 4` → `-T ${task.cpus}` + bump `cpus` (currently 1, so `-T 4` is over-subscribed) | Correct core use. |
| `multimapper_transcript_assignment` (410), `multimapper_exon_assignment` (451) | `samtools sort -n` → `samtools sort -n -@ ${task.cpus}` + bump `cpus` | Faster name-sort feeding gawk. |
| `sort_index_bam` (599-600) | `samtools sort -@ ${task.cpus} -m 2G --write-index` (fuses index into the sort, drops the separate `samtools index` call); `samtools view -c -f 16 -@ ${task.cpus}` | One pass instead of sort+index; threaded. `--write-index` yields the same `.bai`. |
| `merge_*` (476,491,508) | `samtools merge -@ ${task.cpus} ...` + bump `cpus` | Threaded merge. |
| `count_high_conf_annotated_umr_multimap` (523), `io_count` (661) | add `-@ ${task.cpus}` to `samtools view`; right-size over-allocated `cpus` (io_count has cpus=4 but samtools view/awk are single-threaded) | Better scheduler packing; minor speed. |
| `dedup` (607) | right-size `cpus` 4→1 (umi_tools dedup is single-threaded) | No speed change, frees 3 cores per dedup task for parallel samples → higher throughput across many samples. |

VERIFICATION — IMPORTANT SCOPE CORRECTION (post-review B2/B3): because `umi_tools dedup` and the gawk
multimapper assignment select a representative by record-traversal order, ANY change here that alters BAM
record order UPSTREAM of them is NOT automatically R0. That implicates every `samtools sort`/`sort -n`/`merge`
threading change above (rows for `initial_feature_count`, `multimapper_*` name-sort, `sort_index_bam`
coord-sort feeding dedup, threaded `merge`). These ship only after the **order-perturbation test**
(regression §5) shows zero count drift on the covering profiles; if drift appears, pin ordering (stable
secondary sort key / `--outMultimapperOrder Random --runRNGseed` / order-invariant dedup) before adopting.
Genuinely R0 without the gate: `featureCounts -T` (assignment is thread-count-independent), `samtools view -c`
counts (order-independent), and pure `cpus` right-sizing (scheduler only, command unchanged) — e.g. the
`dedup` cpus 4→1 and `io_count` cpus changes are safe immediately.

## Tier 2 — R1 structural rewrites (same tool semantics, prove per-stage)

| Process (line) | Change | Why it helps | Risk note |
|---|---|---|---|
| `umr_exon_assignment` (369) | Replace the `samtools view -h \| sed 's/\tXS:Z:[^\t]*//' \| samtools view -h -b` SAM round-trip with `samtools view --remove-tag XS -b` | Avoids decoding the whole BAM to SAM text and back | Must confirm `--remove-tag XS` removes exactly the same tag the sed targets, and the installed samtools supports it. nf-test old-vs-new on a fixture. |
| `initial_feature_count` (307) | Drop the pre-`featureCounts` `samtools sort` IF the published `*_Aligned.sortedByCoord.out.bam` artifact is not consumed downstream/by the customer (featureCounts does not need coordinate sort) | Removes a full coord-sort of the largest BAM | FLAG: requires a published-output-consumer audit. If the sorted BAM is a deliverable, keep it (but still thread it via Tier 1). |
| `star` (209) | Drop `--outReadsUnmapped Fastx` IF the `*_Unmapped.out.mate*` files are unused by external QC/deliverables | Less alignment IO | FLAG: audit consumers first. |
| download_* (22-110) | Parallelize R1/R2 downloads / use `s5cmd` or `aws s3 cp` with higher concurrency; combine the 6 tiny metadata downloads | Faster cold-start data staging | R0 in practice (same bytes) but treat as R1 since tooling changes; trivial to verify by checksum. |
| `merge_lanes` (145) | Keep `cat *.gz` (already optimal — gzip member concatenation is valid) | n/a | No change recommended; flagged only to record it was reviewed. |

## Tier 3 — R2 engine rewrites (internal Rust fusions; gated, may change outputs)

These are the largest theoretical speedups and the source of the user's "rewrite to make faster" intent.
Each REPLACES an algorithm, so each needs an explicit equivalence verdict before adoption.

| External processes replaced | Internal binary (image) | What it fuses | Output-risk verdict |
|---|---|---|---|
| `filter_for_UMRs_mismatch` (319), `filter_for_multimappers_mismatch` (377), `count_high_conf_annotated_umr_multimap` (512), several `samtools view -c` | `bam-splitter` (`bam-tools:0.2`) | One Rust pass splits UMR/multimapper, applies the `[NH]`/`[nM]` filters, emits all counts | **R2-recoverable.** The filters are exact predicates (`[NH]==1 && [nM]<=3`). A Rust reimplementation CAN be byte/record-equivalent. Build an old-vs-new nf-test on the featureCounts BAM; adopt only on EXACT record-set match. Internal's binary differs (mixed-species split) so likely needs an external-mode build. |
| `umr/multimapper_transcript_assignment` + `*_exon_assignment` gawk (392-462) | `bam-assigner` (`bam-tools:0.2`) | Replaces the two-pass `gawk` SAM-body assignment, streamed from `samtools sort -n`, no `.sam_body` temp files | **R2.** The gawk script IS the assignment algorithm; a Rust rewrite must reproduce its tie-breaking exactly. External gawk differs from internal's. High proof burden; per-read assignment diff on a fixture required. |
| `dedup` (607) + `io_count` (642) + `count_matrix` (674) | `count_table_builder` (`count-table-builder:0.2`) | Replaces `umi_tools dedup` with a custom UMI grouper, derives dedup stats analytically, builds the count matrix in one pass | **R2-breaking (expected).** A from-scratch UMI grouper will NOT reproduce `umi_tools` directional-adjacency dedup edges exactly → h5ad counts, called cells, and dedup.log all shift. Treat as a deliberate re-baseline: quantify the delta (regression_strategy.md §8) and get user sign-off. Do NOT adopt as a transparent speedup. |

## Tier 4 — R2 dependency bumps (do NOT adopt without output proof)
- `count_matrix` container `scanpy_anndata:0.0.4` → `0.0.7`: anndata/scanpy version change can alter h5ad serialization/values. Only with h5ad equivalence diff.
- `dedup` container `umi-tools-csgx:0.9` → `1.0`: umi_tools version change can alter dedup results. Keep `0.9` unless output-verified.

## Recommended sequencing
1. Land all **Tier 1 (R0)** in Phase 3 — pure win, near-zero verification cost (determinism probe + end-to-end equivalence).
2. Land **Tier 2 (R1)** items that pass their consumer audit + per-stage nf-test.
3. Treat **Tier 3 (R2)** as a separate, explicitly-gated workstream:
   - `bam-splitter` (R2-recoverable) first — best risk/reward, provable equivalence.
   - `bam-assigner` next — heavier proof.
   - `count_table_builder`/dedup last — bring the user the quantified output delta and decide re-baseline vs keep `umi_tools`.
4. **Tier 4** only if a Tier 3 binary forces a specific dependency and the diff is clean.

## Things confirmed NOT portable (would change outputs — do not touch)
- `qc` external runs `--no-filtered` / no `--sss-trim`; internal adds `--sss-trim` + filtered outputs. Different args → different QC outputs. Keep external invocation.
- `categorize_reads` external takes `--raw_count_matrix_h5ad` + `--fastp_json`; internal takes `--annotated_h5ad` + `--cell_caller_method`. Different semantics.
- internal `features_file` adds `sed 's/mm10___/mm10_/g'`; external omits it — changes the features file. Do not port.
