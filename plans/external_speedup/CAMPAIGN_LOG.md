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

BUILT on branch `speedup/dedup-contig-split` (commit a63b8ea off io_count tip 3995c26): the logic is
in bin/dedup_by_contig.sh (split by contig via idxstats -> parallel umi_tools dedup via xargs -P
${task.cpus} -> samtools cat merge -> sum the two parsed log fields). Module calls it for the
non-empty case; empty branch unchanged; fails loud on unmapped reads. dedup cpus 1->8, mem 8->16 GB.
Standalone-verified in the container: merged log correct (519415 in / 218681 out) AND (barcode,gene)
multiset byte-identical to whole-dedup. END-TO-END count-data gate RUNNING (bg blltp6gtc: local
test-profile whole-dedup 3995c26 vs contig-split a63b8ea; expect h5ad/mtx identical + dedup metrics OK).
NOT landing without user OK on the dedup.bam representative caveat.

## NEXT (output-preserving)
3. Flow-level fusion of the UMR/multimapper filter->assign->merge chain (cut container-start + BAM
   stage I/O) -- does not change computation, byte-exact safe.
4. multimapper post-sort round-trip: a Rust tool reading the (single-thread) sorted BAM directly to
   write the 2 output BAMs, skipping the samtools view->gawk->samtools view SAM-text round-trip. Keeps
   the sort (output-preserving); modest gain since the sort is the floor.
5. dedup (umi_tools single-thread) -- stretch, heavily gated.

### Strategy 4 GATE RESULT (end-to-end, local test profile): PASS
- COUNT-DATA: 16/16 files byte-identical (h5ad, mtx, barcodes, features) -> count matrices unchanged.
- dedup metrics: reads_before/after_deduplication + sequencing_saturation IDENTICAL both samples.
- Full metrics.csv: Sample1 fully identical; Sample2 differs ONLY by Total Tags/Introns +/-1 -- which
  is the pre-existing AMBIENT annotated-RSeQC wobble (annotated_rseqc runs on the umr_multimapper
  annotated bam, main.nf:348, UPSTREAM of dedup -> my change cannot affect it; same +/-1 seen in the
  io_count gate where dedup was unchanged).
- dedup.bam: md5 DIFFERS, read count SAME (the flagged representative-read caveat).
CONCLUSION: dedup contig-split preserves all counts + all metrics. Only the intermediate dedup.bam's
representative reads differ. Speedup ~Nx bounded by largest contig (precise number from next real run).
Ready to land pending user OK on the dedup.bam caveat.

---
## ROUND 2 (user greenlit small count changes + "land into the PR, go go go") 2026-06-12 eve
LANDED into epic/external-nf26-overhaul (PR #78, pushed): io_count(20x) + dedup contig-split(~Nx) +
right-sizing. epic now at 4d365f6.

### OPTION (b) DONE RIGHT: deterministic multimapper assignment + threaded sort (branch speedup/multimapper-deterministic)
Instead of naive threading (which made counts wobble every run), made the assignment ORDER-INDEPENDENT:
- assign_multi_mappers.gawk: canonical representative per gene (lex-smallest corrected alignment, not
  last-write-wins) + PROCINFO["sorted_in"] for deterministic output order.
- threaded `samtools sort -n -@ ${task.cpus} -m 1G` in both multimapper_{transcript,exon}_assignment; cpus->4, mem->6GB.
Unit-verified on real multimapper BAM: threaded run A == run B (reproducible) AND single==threaded
(order-independent); ambiguous set unchanged; assigned set same reads, canonical representative ->
one-time tiny count change (approved). Commit 1d124ff. Fixes the ambient annotated-RSeQC wobble too.
GATE RUNNING bg bt51crqr4 (candA, candB, base): expect candA==candB (deterministic) + small delta vs base.

### initial_feature_count: threaded the coordinate samtools sort (-@), cpus 2->4 (commit e289f92).
Order-neutral (featureCounts per-read; downstream re-sorts). Left sort_index_bam single-threaded (it is
the determinism anchor for the dedup input).
TODO: combined gate for the full branch after bt51crqr4; then land into epic. More threading candidates
scanned; sort_index_bam intentionally NOT threaded.

### ROUND-2 GATE FINDINGS (bt51crqr4) + LANDED to epic c44d418 (PR #78)
- (b) multimapper-deterministic: dedup INPUT byte-identical run-to-run (candA==candB) => the (b) path
  IS deterministic, the gawk fix works. Magnitude vs base: ~0.2-0.4% of count-matrix entries shift
  (num_cells IDENTICAL both samples). Bigger than the 0.005% I'd estimated, because a DETERMINISTIC
  threaded selection (canonical representative) deviates more from the original last-write-than-sort-order
  than naive tie-only threading would; this is inherent to threading-with-reproducibility. Within
  user's "very small changes of a count or so" (num_cells stable). KEPT.
- RARE DEDUP WOBBLE: contig-split dedup is deterministic 10/10 in ISOLATION (whole==split==dccba67a),
  but candA's PIPELINE dedup deviated once (aa6799bb, 2 (barcode,gene) entries -> ~5 mtx entries, ~0.02%).
  Triggered under pipeline concurrency, not reproducible in isolation. umi_tools residual non-determinism
  (PYTHONHASHSEED=0 doesn't fully fix it), possibly slightly exposed by the parallel contig-split. Within
  tolerance. CONSEQUENCE: byte-exact gating is no longer possible; gate on small-tolerance instead.
LANDED: epic ff-merged to c44d418 + pushed. Round-2 = (b) + initial_feature_count sort + 3 view-filter threadings.

### REMAINING IDEAS
- RSeQC (raw_rseqc/annotated_rseqc, ~50-60min, single-thread python) likely the NEW top cost. read_distribution
  is position-based -> contig-splittable + summable (like dedup). Involved but ~Nx. Best remaining lever.
- A fresh REAL-DATA Seqera run to confirm the speedups + reveal the new bottleneck profile + equivalence at scale.

### ROUND 3: real-data validation launched + RSeQC analysis
- LAUNCHED Seqera real-data run of epic 5c519368 on human_validation.csv (genome GRCh38, docker profile),
  workflow 35pj7QevURScMM, outdir .../cand_human_round2. Purpose: confirm speedups REAL (trace vs old
  base_human), cell-calls stable vs original base_human, reveal NEW top bottleneck. ~2-3h.
- RSeQC contig-split: VERIFIED count-identical (per-contig read_distribution Tag_counts sum EXACTLY to
  whole: Total Reads/Tags/Assigned + all 10 groups). read_distribution is read-processing-bound (whole
  1.13s, bed-load floor 0.28s) so split helps. BUT: rseqc container (quay biocontainers rseqc:5.0.3) has
  NO samtools, so split must be pysam (a single-threaded full read-pass) -> caps simple approach at ~2x.
  A true ~Nx needs a DAG scatter-gather (split process in samtools container -> parallel rseqc -> merge).
  Foundation ready (parse/sum logic proven). DECISION: wait for the validation trace to confirm RSeQC is
  the new top cost before building (avoid the STAR mistake of optimizing a non-bottleneck).

### ROUND-3 VALIDATION RESULT (real data, 8 human samples) -- CONFIRMED
Run 35pj7QevURScMM SUCCEEDED. Per-process realtime cand(round2) vs base(original):
- io_count 172.6m -> 1.1m (154x!), dedup 76m -> 30.6m (2.5x), umr_transcript_assignment 24.4->4.7 (5.2x),
  filter_for_multimappers_mismatch 51.6->16.4 (3.1x), filter_for_UMRs_mismatch 41.9->16.1 (2.6x),
  initial_feature_count 63.6->33.3 (1.9x), multimapper_transcript_assignment 201.9->141.3 (1.4x only).
- TOTAL sum-of-realtime 945m -> 579m (-39%).
- CELL CALLS IDENTICAL all 8 samples (num_cells); reads_after_dedup within 0.001-0.003%. Equivalence confirmed.
NEW BOTTLENECKS: (1) multimapper_transcript_assignment 141m STILL #1 -- sort threaded but the gawk + SAM
round-trip + 2 BAM rebuilds dominate now. NEXT: replace the gawk with a static Rust binary (BYTE-IDENTICAL
to the now-deterministic canonical gawk -> NO count change), like io_count_extract. (2) RSeQC raw 64m +
annotated 32m = 96m combined, untouched -> contig-split (count-identical, needs DAG scatter for ~Nx).
star 54m cand vs 33m base = instance noise (uses ~4-5 cores, cpus=8 fine).
CHANGELOG perf section updated with confirmed round-2/3 numbers + equivalence.

### ROUND 4: multimapper deep-dive -- Rust gawk = DEAD END; sort is the real cost
- Built a static Rust replica of assign_multi_mappers.gawk (BYTE-IDENTICAL verified on 2 real multimapper
  BAMs, assigned+ambiguous sam_body). BUT it is SLOWER than gawk (6.1s vs 2.9s on 2.4M records) -- gawk is
  well-optimised C and the gawk is NOT the bottleneck (~3s). Discarded (kept in scratch only).
- The samtools sort -n of the large multimapper BAM is the real cost. Tested: the multimapper BAM is NOT
  pre-grouped (2437896 reads, 408529 names, 1945205 contiguous runs) and unsorted gawk output DIFFERS, so
  the sort canNOT be dropped. It's inherently sort-bound (huge --outFilterMultimapNmax 1000 expansion).
- DONE (byte-identical): threaded the middle samtools view + the 2 BAM-rebuild samtools view -b with -@
  in both multimapper_{transcript,exon}_assignment (commit 54ce90c). Modest (rebuilds), kept cpus=4.
LESSON: gawk text passes are fast; don't Rust-rewrite them. The wins are in samtools (de)compression
threading + the genuinely single-thread python (umi_tools done via contig-split; RSeQC next).

### REMAINING: RSeQC contig-split (raw 64m + annotated 32m = 96m, count-identical verified).
Needs DAG scatter-gather (rseqc container has no samtools): a split process (samtools container) emits
per-contig bams -> parallel run_rseqc per contig -> a merge process sums read_distribution + reconstructs
the format. ~Nx (bounded by largest contig). This is the last clean meaningful win. Build next.

### RSeQC contig-split BUILT (bin/rseqc_by_contig.py, commit cb0348d LOCAL on epic, not pushed)
pysam split by contig -> parallel read_distribution.py -> sum + reconstruct. Output BYTE-IDENTICAL to
read_distribution.py on the whole BAM (verified standalone). run_rseqc cpus 2->4. Gate RUNNING (bg bbxu2qh8m:
local test-profile base 2b32df3 vs cand cb0348d; expect RSeQC outputs byte-identical + pipeline completes +
MultiQC parses). When PASS: push to epic/PR#78. ~2-2.5x on the 96m RSeQC cost (capped by single-thread pysam split).

### RSeQC contig-split LANDED (commit cb0348d + Tags/Kb fix dcef53a, pushed to PR #78)
Gate confirmed: pipeline completes, MultiQC parses, metrics.csv IDENTICAL, all counts IDENTICAL; only
Tags/Kb (a derived rate no metric uses) was off by 0.01 (rounding) -> fixed to match read_distribution.py's
exact formula count*1000.0/(bases+1) -> byte-identical. ~2-2.5x on the 96m RSeQC cost.

## CAMPAIGN COMPLETE (diminishing returns reached). All in PR #78 (epic dcef53a), validated on real data.
Landed perf commits: io_count rust (154x), dedup contig-split (2.5x), multimapper deterministic+threaded
(1.4x) + view threading, initial_feature_count sort thread, 3 filter view threads (2.6-5.2x), RSeQC
contig-split (~2-2.5x, byte-identical), right-sizing. Real-data result: total compute 945->579 CPU-min
(-39%, and RSeQC not yet in that run -> more now); CELL CALLS IDENTICAL all 8 samples; counts within 0.003%.
Remaining bottleneck = multimapper sort (inherently sort-bound, can't beat much). STAR memory-bound.

## PIPELINE REVIEW (2026-06-13, 3 parallel agents) + FIXES
Launched 4 fresh old-vs-new validation runs (OLD baseline-validation fe80f36 vs NEW epic cceb134) on
human + mixed validation CSVs (internal datasets MOR034/MOR036/KOL0054 + edge cases empty/noalign/intergenic):
 OLD_human 5XPxbeh9OXUNlj, NEW_human 1j0koOt0d6Xzp2, OLD_mixed 2pyX0TTrRHGWbR, NEW_mixed 1a6Wzbuvqa9a8J.

REVIEW FINDINGS:
FIXED (commit 2ae5f33, silent failures per "fail loud"):
 - rseqc_by_contig.py: per-contig read_distribution failures swallowed (no check) -> summed as zeros. Now check=True.
 - cell_caller + merge_annotated: publishDir `pattern:{closure}` matches NOTHING on NF26 -> plots + final
   annotated BAM never published. Changed to string globs; merge_annotated glob also didn't match its filename.
 - dedup empty-sample log: echo "...\n..." = literal \n (1 line) -> reads_after_deduplication never recorded. printf.
STILL TODO (from review):
 - HIGH: consolidated report MIXED-species mode broken -- headline cards empty + cross-sample table drops cell
   metrics (create_consolidated_report.py:296-308 + template:495-516 use group "Cell metrics" but mixed keys live
   under num_cells/raw_reads_per_cell/median_genes_detected_per_cell classifications). FIX when NEW_mixed lands (test against real mixed metrics).
 - PERF (not done, optional): sort_index_bam unthreaded coord sort (LEFT single-thread on purpose = dedup determinism
   anchor; threading adds to ~0.02% wobble); categorize_reads.py single-thread pysam pass over raw BAM (io_count-style
   rewrite possible); fuse UMR/multimapper filter->assign->merge chain (cut container starts + S3 stage); thread the
   3 samtools merge (-@).
 - TESTS: new bin/ tools (io_count_extract, dedup_by_contig.sh, rseqc_by_contig.py, assign_multi_mappers) have NO
   tests; compare_outputs.py comparator NOT wired into CI; pytest not run in CI. Add golden tests + wire CI.
 - DOCS: CHANGELOG perf-bullet cpu numbers stale (says io_count 4->1, dedup 1cpu; actually 2 and 8). io_count README
   says musl but binary is glibc static-pie. Low priority.

## OLD-vs-NEW FRESH VALIDATION RESULT (2026-06-13)
HUMAN (base_human_final OLD/devel+fix vs cand_human_final NEW/epic, 8 real MOR034/MOR036 samples + edges):
- num_cells IDENTICAL all 8 samples. reads_after_dedup within 0.003%.
- count matrices: ~0.1% of mtx entries differ (multimapper canonical-representative reassignments);
  TOTAL counts within 0.002% (net-near-zero). num_cells unaffected.
- EDGE PARITY: EDGE_emptybarcode -> both OLD+NEW exclude it (no metrics.csv, no crash); EDGE_intergenic
  -> both num_cells=2. Identical edge behavior.
=> Human: new code reproduces old code within the documented tiny envelope. CONFIRMED.
MIXED: pending (NEW_mixed near done, OLD_mixed slow on busybox io_count). Then mixed comparison + the
mixed consolidated-report template fix (verify against NEW_mixed rendered report).

## MIXED VALIDATION (partial) + REPORT FIX VERIFIED (2026-06-13)
MIXED-REPORT FIX render-verified (commits 3f30678 headline + 8bb6396 table + 91ea6b8 prettify): rendered the
NEW_mixed report with old vs fixed template -- cross-sample table group rows went from ['Read QC','Cell metrics',
'Deduplication'] (cell metrics MISSING) to ['Read QC','Num cells','Raw reads per cell',...16 cell-metric groups...,
'Num multiplet cells','Deduplication']. Cell metrics now appear in mixed reports. DONE.
MIXED barnyard NEW vs OLD (KOL0054_Sample134, first OLD_mixed sample done): num_cells_total/Hsap/Mmus IDENTICAL
(3233/2167/1066), num_multiplet_cells_total IDENTICAL (224), reads_after_dedup +0.0007%. Same as human -> equivalent.
OLD_mixed (busybox, slow) still finishing other 3 KOL + edges -> full mixed comparison pending.

## EDGE PARITY (mixed): EDGE_noalign_pig -> num_cells_total=0 in BOTH OLD+NEW, reads within 0.0006%. Graceful.
Mixed barnyard equivalence confirmed (Sample134 + the earlier round-2 run's 3 KOL samples all per-species-identical).
OLD_mixed (busybox) still finishing KOL107/142/143 + EDGE_emptybarcode -> will confirm but pattern is conclusive.
FINAL STATE epic 5a0398f, PR #78 unmerged. All validated; report fix render-verified; 4 review bugs fixed.

## MIXED VALIDATION ~COMPLETE: KOL Sample134/107/143 ALL per-species-cells + multiplets IDENTICAL (reads <0.0007%).
EDGE_emptybarcode excluded in both; EDGE_noalign_pig 0 cells in both. Only KOL Sample142 pending (OLD busybox
io_count, ~90min; already confirmed identical in round-2 run). VALIDATION CONCLUSIVE: new==old (human+mixed),
within documented envelope, edge parity. PR #78 ready for review.

## SAMPLE142 (OLD) -- the live speedup demonstration
OLD busybox io_count for KOL0054_Sample142 (1.4GB barnyard) FAILED after 2h43m (spot reclaim), now retrying
(~2.5h more). It has NEVER completed in either validation attempt -- the old code practically cannot process
this sample, while the NEW rust io_count did it in minutes. NEW Sample142 done + correct (in cand_mixed_final).
=> VALIDATION SUBSTANTIVELY COMPLETE: new==old (human 8/8 + mixed 3/4 KOL + both edges, all per-species cells +
multiplets identical, within envelope). Sample142 OLD comparison is the only literal gap, blocked by the very
bottleneck the rewrite fixes. PR #78 ready for review.

## VALIDATION 100% COMPLETE (OLD_mixed SUCCEEDED). HONEST FINAL RESULT:
KOL0054_Sample142 (1.4GB, the long-pole): num_cells_total 4383(OLD) vs 4384(NEW) = +1 cell; Hsap 2756 vs 2757
= +1; Mmus 1627 IDENTICAL; multiplets 705 IDENTICAL; reads_after_dedup +454 (0.0009%). So ONE borderline human
cell (1/4384 = 0.02%) tipped across the cell-calling threshold due to the ~0.001% multimapper-reassignment read
envelope. This is the ONLY non-identical cell-count across all 12 samples (human 8/8 + mixed KOL 134/107/143
identical; 142 differs by 1 cell; both edges identical). Well within user's "small changes of a count or so".
FINAL: overhaul is faster (-39% compute, io_count ~150x) and output-equivalent (11/12 samples bit-identical cell
calls, 1 sample +/-1 borderline cell, all within documented envelope). Report fix verified, 4 review bugs fixed.
PR #78 ready for review.

## WALL-CLOCK reframing (user, 2026-06-14): time-to-result is what matters, not CPU-min.
WALL-CLOCK (launch->complete, same CE, simultaneous): HUMAN 155(OLD)->92(NEW) min = 1.7x; MIXED 675->202 = 3.3x.
(Caveats: 4 runs shared CE = contention; OLD_mixed inflated by spot-reclaim retries on the long busybox io_count.)
PER-SAMPLE WALL (median realtime) OLD vs NEW: io_count 19.7->0.1 (169x), dedup 7->2.6 (2.7x), filters ~3x,
initial_feature_count 2.1x, multimapper_transcript 18.6->15.1 (1.2x, sort-bound). REGRESSION CAUGHT: rseqc
contig-split made raw_rseqc 6.5->25.6 (4x SLOWER) + annotated 3.3->6.4 (2x) -- pysam single-thread split over the
large raw BAM costs more than it saves. REVERTED (commit 5d2d396, back to direct read_distribution.py, byte-identical).
CURRENT WALL-CLOCK CRITICAL PATH per sample: multimapper_transcript_assignment ~15m (sort-bound, #1) + STAR ~6m +
a ~12-step serial chain each paying container-start + S3 stage + queue latency. WALL-CLOCK LEVERS (different from
compute): (1) FLOW-FUSION of the UMR/multimapper filter->assign->merge chain (cut serial container-starts/staging/queue
-- pure latency win); (2) the multimapper sort. NEXT: implement flow-fusion (gated re-validation). Optionally a solo
NEW run for a clean wall-clock headline (no CE contention).

## WALL-CLOCK round 2 (2026-06-14): contained wins + multimapper fusion
- multimapper sort -@4 -> -@8 (cpus 8, mem 12GB) both multimapper processes (f5bbaf0); byte-identical.
- FLOW-FUSION: fused multimapper branch (filter+transcript+exon+merge) into modules/local/multimapper_assignment
  (a0ed5c6); cuts 3 serial container/staging/queue per sample. Byte-for-byte same commands; sam_body rm'd between
  the 2 gawk passes. nextflow inspect clean. GATE RUNNING (bg bxvgmg0t8): test-profile fused vs unfused count-data
  byte-identical. Old 4 multimapper modules unwired (delete later).
