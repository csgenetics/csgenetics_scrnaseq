# CRITICAL finding: current pipeline is non-deterministic at the count-matrix level

Date: 2026-06-11. Source: Phase 0 determinism probe — two identical Seqera runs of CURRENT `devel`
(commit `81c094b`), `test` profile, `NXF_VER=25.10.2`, on CE `ed0zjbIeIvKoUPuGsy8fA`.
Outputs: `s3://csg-nextflow/external_nf26_overhaul/baseline/test_run{1,2}`.

## What the probe showed (run1 vs run2, same code, same input)
- **Published file SET: identical** (incl. cell-caller-threshold-derived h5ad names `Sample1.166`, `Sample2.100`).
- **All metric CSVs identical:** every `report/**/*.metrics.csv` and `report/multisample_out.csv` byte-identical.
- **RSeQC identical.** `dedup.log` differs ONLY in timestamps/pid/instance-UUID; the integer counts
  (Input Reads, Number of reads out, positions deduplicated) are identical.
- **Count matrix DIFFERS for Sample2** (raw AND filtered): same dims, same 51,079 triplets, but 2 entries
  shifted — barcode 5100: gene 164 `32↔31`, gene 166 `81↔82`. A **net-count-preserving reassignment** of one
  read between two genes. Sample1 matrix identical.

## Mechanism (matches review findings B2/B3)
STAR runs `--runThreadN 8` (hardcoded). Even at a FIXED thread count, STAR's multithreaded Unsorted-BAM
record ORDER varies run-to-run. That order feeds the order-sensitive multimapper assignment
(`assign_multi_mappers.gawk`, last-write-wins on the gene tag) and/or `umi_tools dedup` representative
selection. For a read that multimaps/ambiguously assigns between genes 164 and 166, the winner flips between
runs → the per-gene-per-cell count shifts by ±1, net-preserving. Summary metrics (totals, #genes, #cells)
are unaffected because the reassignment conserves totals — which is why `metrics.csv` is stable.

## Why this is critical
1. **Byte-exact output-equivalence is impossible even for baseline-vs-baseline.** The count matrices (h5ad/mtx
   — customer deliverables) are intrinsically non-deterministic. A naive EXACT comparator (regression §2 Class C)
   would FALSE-FAIL comparing the baseline to itself.
2. **B2/B3 is not theoretical — it is already happening.** Any optimization that perturbs record order (STAR
   threading, samtools sort threading) sits on top of pre-existing non-determinism. The equivalence question
   becomes "does the change move outputs BEYOND the intrinsic run-to-run envelope," not "are they byte-equal."
3. The probe-first methodology paid off: this surfaced before a single line of code changed.

## Two ways forward (decision needed)
- **(A) Characterize and accept the envelope.** Run the baseline N times, measure run-to-run variance per output
  class, and define equivalence as "candidate delta vs baseline is within the baseline-vs-baseline envelope."
  Comparator becomes statistical; baseline is a distribution, not a point. Keeps current outputs as-is.
- **(B) Make the pipeline DETERMINISTIC first, then byte-exact thereafter (recommended).** Pin the order
  sources: STAR `--outMultimapperOrder Random --runRNGseed <fixed>` (+ `--outSAMmultNmax` as needed) and a
  stable secondary sort key on the name-sorts feeding the gawk; verify the count matrix is then reproducible
  run-to-run. This is a one-time, blessed output change (the deterministic values differ slightly from any
  single current run) that ALSO improves the product — reproducible results are valuable for a customer-facing
  pipeline — and lets every subsequent optimization be held to a clean byte-exact bar.

## Recommendation
Path B. Establish determinism as the first substantive change (a blessed re-baseline of the count matrices,
metrics already stable), then the entire optimization programme is verifiable byte-exact. Quantify the B-vs-current
delta (expected tiny: a handful of net-preserving per-gene reassignments) and surface it for sign-off, same as
the dedup re-baseline.

## RESOLUTION (2026-06-11): Path A chosen after investigation

Path B was attempted and proved invasive for an irreducible residual:
- STAR `--outMultimapperOrder Random --runRNGseed` DID make the pipeline byte-deterministic (run-pair EQUAL),
  but it replaces STAR's best-scoring primary multimapper alignment with a random (seeded) one, which
  **corrupted the RSeQC read-distribution metrics** (e.g. Sample1 CDS 29.65%->25.73%, Intergenic 16.59%->20.22%) —
  a real change to customer-reported numbers. Rejected.
- Order-normalizing sorts alone (LC_ALL=C sort before the multimapper gawk; name-sort before the dedup
  coordinate sort) did NOT achieve full determinism. A local test proved `samtools sort` is stable (preserves
  input order for equal coordinates) and that name-sort-then-coordinate-sort yields a deterministic order — so
  the dedup *input* is deterministic, yet ~2 count entries still flipped. The residual is therefore a
  **genuinely ambiguous multimapper read** (a molecule mapping to two genes); the final gene attribution is
  irreducibly arbitrary. Forcing it deterministic is itself an arbitrary re-baseline of a customer deliverable.

Decision (user, 2026-06-11): **Path A.** Revert all Phase 1d determinism edits (commit `f251c59`); keep only the
NF26 parse fix. Characterize the inherent non-determinism as a **micro-envelope** and gate future changes against it:
- customer metrics CSVs: **byte-exact** (already stable run-to-run);
- count matrix: **per-barcode column sums byte-exact** (net-preserving — confirmed: column sums identical, only
  ~2 of 51k entries flip between two genes for one barcode), differing entries bounded by a measured N.
Implemented in `tests/regression/compare_outputs.py --envelope-max-flips N` (commit `7afed86`). This gives high
no-regression confidence with ZERO output change to the customer pipeline, and avoids re-baselining a deliverable.
