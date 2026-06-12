# External pipeline speedup campaign

Objective (user, 2026-06-12): make the pipeline run FASTER with byte-identical output, by any
means. Safety gate: the count data is run-to-run deterministic, so "same output" is machine-checkable.
Every change is gated; winners on `speedup/*` branches off epic; nothing merges; reviewed Monday.

Method: cheap LOCAL unit-equivalence on real fixtures for fast iteration; local test-profile
double-runs for end-to-end count-data byte-exactness; Seqera only to confirm speed of winners.
Fixtures in /nssd2/humebc/scratch/external_speedup/ (BAMs not committed).

Measured bottlenecks (8-sample human + mixed traces, total realtime): multimapper_transcript_assignment
(#1, ~200 min/sample, single-thread name sort + gawk, 243 GB I/O), io_count (#2, 105-139 min on
mixed, I/O-bound busybox-awk), dedup (umi_tools single-thread). STAR is NOT a bottleneck (memory-bound).

GATE NOTES learned:
- The count data (h5ad) is deterministic run-to-run; compare h5ad + decompressed mtx (the .mtx.gz md5
  differs only by gzip-header mtime). barcodes/features also byte-exact.
- SECOND ambient non-determinism source (separate from the umi_tools one PYTHONHASHSEED fixed): the
  multimapper gawk iterates `for key in array` in HASH order, so the annotated BAM byte order (and its
  RSeQC) wobbles by +/-1 read run-to-run. Does NOT affect counts. CHANGELOG "byte-reproducible" is
  therefore overstated for the annotated RSeQC; soften or fix.

---

## DONE: resource right-sizing (commit 45a843b on epic)
star cpus 16->8, io_count 4->1 (later 2), initial_feature_count 4->2, umr/multimapper_exon_assignment 4->2.
STAR runThreadN restored to original 8 + pinned (a prior change to =task.cpus=16 had silently altered
per-barcode counts -- a real output regression now fixed). Throughput change; output byte-unchanged.

## STRATEGY 1: io_count rewrite -- WIN, VALIDATED (commit 1b2a248, branch speedup/io-count-rust)
Prod container ships BUSYBOX awk (slowest). samtools decode ~1-2s; the awk text pass was the whole cost.
Replaced with static musl Rust binary `bin/io_count_extract` (src tools/io_count_extract) -- byte-for-byte,
streaming, runs in the existing public container (no image change, dodges Wave constraint).
Byte-exact (md5) on human 158MB (39->2s), empty+intergenic edges, 1.4GB mixed (409->20s, 52.5M lines). ~20x.
End-to-end gate PASSED: count matrices (h5ad) byte-identical. io_count cpus->2 (now samtools-decode-bound).

## STRATEGY 2: multimapper sort threading -- TRIED, REVERTED (output-coupled). KEY FINDING.
Lever A = thread the single-threaded `samtools sort -n` in multimapper_transcript_assignment +
multimapper_exon_assignment. SPEED was great (proxy 1.4GB BAM: -@1=349s, -@4=102s, -@8=81s, -@8 -m2G=65s,
~3.4-5x). But the COUNT-DATA GATE FAILED: all h5ad differed (e.g. mean_genes_detected_per_cell
66.2298->66.2330; Sample2 Total Tags +4).
ROOT CAUSE: the gawk keeps one alignment per gene by LAST-WRITE-WINS, so which alignment survives
depends on the PROCESSING order = the sort's tie order. Single-thread `samtools sort -n` has a
deterministic tie order (hence stable counts; gate 1 proved it). Threading reorders equal-name ties ->
a different (equally valid) multimapper alignment is kept -> it dedups differently -> different counts.
=> The output is INTRINSICALLY COUPLED to the single-thread sort order. Threading it, removing it, or
reordering it ANY way changes counts. The change is tiny (~0.005%, genuinely-ambiguous multimappers,
same envelope class as the PYTHONHASHSEED one-time change) but it is NOT byte-identical. Reverted to
honor "same output". Branch reset to 1b2a248.

### DECISION FOR USER (blocks a big win on the #1 bottleneck)
The largest single cost (multimapper sort, ~200 min/sample) cannot be sped up while keeping
byte-identical counts -- the result depends on the serial sort order. Options:
  (a) STRICT (default, current): leave it; the #1 bottleneck stays. io_count + right-sizing + flow
      fusion are the wins.
  (b) ACCEPT a one-time envelope re-baseline (like PYTHONHASHSEED): thread the sort for ~3-5x on the
      #1 bottleneck; counts shift by ~0.005% on genuinely-ambiguous multimappers (net-tiny, cell calls
      unaffected). Would also let us make the gawk deterministic (fix the ambient RSeQC wobble).
Awaiting user call. Until then, pursue only output-preserving strategies.

## STRATEGY 4: dedup contig-split parallelism -- COUNT-PRESERVING (verified), flagged for user
dedup = `umi_tools dedup --per-cell` (single-thread, ~80 min, #3 cost). Hypothesis: dedup is
position-local so splitting by contig -> parallel dedup -> merge is equivalent. VERIFIED on real
fixture (Sample2_sorted.bam, chr1+MT, 218,681 reads out): whole vs split-merge gave IDENTICAL
read-count-out AND IDENTICAL (barcode,gene) multiset (md5 6257d01686...) => COUNT MATRIX UNCHANGED.
Structural reason: a UMI group's reads are PCR dups of one molecule (same barcode/UMI/pos) -> same
gene, so the representative choice can't change counts.
CAVEAT: umi_tools representative selection is stateful (advances per group), so split keeps a
DIFFERENT (but equivalent, same-molecule) representative read for ~19% of groups -> the published
dedup.bam is NOT byte-identical, though every count/metric is. Same equivalence class as the ambient
RSeQC wobble.
=> FLAGGED option (c): ~Nx on the #3 bottleneck (limited by largest contig), counts identical,
dedup.bam representatives differ. More attractive than option (b) (which changed counts).
Build = split process + parallel dedup + merge + sum the dedup.log stats; gate on count-data.

## NEXT (output-preserving)
3. Flow-level fusion of the UMR/multimapper filter->assign->merge chain (cut container-start + BAM
   stage I/O) -- does not change computation, byte-exact safe.
4. multimapper post-sort round-trip: a Rust tool reading the (single-thread) sorted BAM directly to
   write the 2 output BAMs, skipping the samtools view->gawk->samtools view SAM-text round-trip. Keeps
   the sort (output-preserving); modest gain since the sort is the floor.
5. dedup (umi_tools single-thread) -- stretch, heavily gated.
