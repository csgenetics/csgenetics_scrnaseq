# External pipeline speedup campaign

Objective (user, 2026-06-12): make the pipeline run FASTER with byte-identical output, by any
means. Safety gate: pipeline is byte-deterministic (PYTHONHASHSEED=0 + STAR runThreadN pinned 8),
so "same output" is machine-checkable. Every change is gated byte-exact.

Method: cheap LOCAL unit-equivalence on real fixtures (a real dedup BAM etc.) for fast iteration;
Seqera runs only to confirm speed of winners + full-pipeline byte-exact. Fixtures live in
/nssd2/humebc/scratch/external_speedup/ (BAMs not committed). Winners cherry-picked onto a speedup
branch; nothing merges; reviewed before landing.

Measured bottlenecks (8-sample human + mixed traces, total realtime): io_count (105-139 min/sample
on mixed, I/O-bound), multimapper_transcript_assignment (243 GB I/O, single-thread), dedup
(umi_tools single-thread). STAR is NOT a bottleneck (memory-bound, ~4-7 cores).

---

## DONE: resource right-sizing (commit 45a843b, local on epic)
star cpus 16->8, io_count 4->1, initial_feature_count 4->2, umr/multimapper_exon_assignment 4->2.
STAR runThreadN restored to original 8 + pinned (a prior change to =task.cpus=16 had silently
altered per-barcode counts). Throughput change; output byte-unchanged. CHANGELOG corrected.

---

## STRATEGY 1: io_count rewrite  -- WIN (pending pipeline-level validation)

io_count = `samtools view <bam> | awk '/XT:/{...}' > bcGeneSummary.txt`. KEY FINDING: the
production container (quay.io/biocontainers/samtools:1.17) ships **BusyBox awk**, one of the
slowest awks. samtools BAM decode is only ~1s (-@4); the awk text stage is the whole cost.

Fixture: MOR036_Sample2.dedup.bam (158 MB, 5,398,818 records). Golden md5 ee3772655a6f00eb195d471794a5eab1.
Benchmarks (byte-exact unless noted):
| strategy | time | note |
|---|---|---|
| busybox awk (PRODUCTION, in container) | 39 s | baseline |
| pysam BAM-native | 11.4 s | MISMATCH (dropped) |
| samtools -d XT prefilter + awk | 23.9 s | no gain |
| gawk + LC_ALL=C (beast) | 7.85 s | byte-exact, but container has no gawk |
| mawk (beast) | 6.3 s | byte-exact, but container has no mawk |
| **Rust static musl binary** (beast) | **2.0 s** | byte-exact |
| **Rust static binary IN CONTAINER** | **2 s** | byte-exact, ~19.5x vs busybox |

Decision: ship a static musl Rust binary `io_count_extract` (reads `samtools view` text on stdin,
replicates the awk byte-for-byte; streaming, low-memory for the 1 GB cap). Static => runs in the
existing public samtools container, NO image change, sidesteps the Wave-pull constraint entirely.
Source: bin/rust/io_count_extract (or scratch during dev). New io_count command:
  `samtools view -@ ${task.cpus} ${f} | io_count_extract > ${sample_id}_bcGeneSummary.txt`

PENDING: byte-exact on edge cases (empty/intergenic) + the 1.4 GB mixed monster (running);
then integrate + full-pipeline Seqera byte-exact gate. Open question for user: commit the
prebuilt binary to the repo vs add a build step (binary is x86_64 musl; deployment is AWS x86_64).

---

## NEXT STRATEGIES (queued)
2. multimapper_transcript_assignment (243 GB I/O, single-thread gawk): thread samtools sort -n,
   kill the BAM->SAM-text->gawk->SAM-text->BAM round-trip, Rust/BAM-native reimplementation.
3. Flow-level fusion of the UMR/multimapper filter->assign->merge chain (cut container-start + BAM I/O).
4. dedup (umi_tools single-thread) -- stretch, heavily gated.
