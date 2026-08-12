# Validation of pipeline version 2.0.0

Version 2.0.0 changed a great deal of how the pipeline runs — a Nextflow 26 migration, several
processes rewritten or parallelised, and a rebuilt report — while intending to change almost
nothing about what it computes. This document records how that was checked, so that the decision
to upgrade can rest on evidence rather than assurance.

**Summary:** across 20 samples in four datasets spanning both supported genome types, cell calling
was identical or differed by a single borderline cell, and read and count totals agreed to within
0.02%. Total compute fell by 39%, and time-to-result improved 1.7x on human data and 3.3x on
mixed-species data.

## Contents

- [What was compared](#what-was-compared)
- [Why the baseline was determinism-controlled](#why-the-baseline-was-determinism-controlled)
- [Results: cell calling](#results-cell-calling)
- [Results: reads and counts](#results-reads-and-counts)
- [Results: performance](#results-performance)
- [Known differences](#known-differences)
- [Multi-lane and heavy-sample testing](#multi-lane-and-heavy-sample-testing)
- [Checking equivalence yourself](#checking-equivalence-yourself)

## What was compared

Four internal CS Genetics datasets were processed through both the 1.x pipeline and 2.0.0 on the
same compute environment, and the outputs compared.

| Dataset | Genome | Samples | Purpose |
|---------|--------|---------|---------|
| Human validation set | GRCh38 | 8 | Primary human equivalence |
| Mixed-species validation set | `mouse_human_mix` | 4 | Barnyard / per-species cell calling |
| Human sign-off set | GRCh38 | 4 | Final confirmation on data not used during development |
| Mixed-species sign-off set | `mouse_human_mix` | 4 | As above, barnyard |

The two sign-off datasets were deliberately chosen to be ones that had played no part in
developing or tuning the changes, so that the final result could not be an artefact of having
optimised against the test data.

Comparisons covered per-sample cell counts (including per-species counts and multiplet counts for
barnyard data), read and count totals, median genes per cell, cell-caller thresholds, the count
matrices themselves, and the data behind the report figures.

## Why the baseline was determinism-controlled

The 1.x pipeline was not fully reproducible run to run. Two independent sources of ambiguity meant
that running the same data twice could produce slightly different count-matrix entries:

- `umi_tools` chose between equally-ranked reads using Python's hash-seed-randomised set
  iteration, which varies per process.
- The multimapper assignment step iterated an associative array in hash order and kept the last
  alignment written per gene, so which of several equally-valid alignments survived depended on
  processing order.

This matters for validation: comparing 2.0.0 against a raw 1.x run would mix *the effect of the
changes* with *1.x's own run-to-run noise*, and there would be no way to attribute a difference to
either.

The final sign-off therefore compared 2.0.0 against a **determinism-controlled 1.x baseline** —
1.x with `PYTHONHASHSEED=0` pinned, plus two metric-neutral crash fixes needed to get it through
the data at all. Both sides are then deterministic in the same way, so any remaining divergence is
attributable to the 2.0.0 changes rather than to noise.

One residual asymmetry is worth stating plainly: `PYTHONHASHSEED` fixes the `umi_tools` source of
non-determinism but not the multimapper one, which only 2.0.0's canonical-assignment change
addresses. A small amount of run-to-run wobble therefore remains on the baseline side. This makes
the comparison *conservative* — it can only overstate the difference, not hide it.

## Results: cell calling

Cell calling is the output most likely to matter to downstream analysis, and it is the most
stable result in the validation.

**Human sign-off set (4 samples):** cell counts agreed to within a single cell, a maximum relative
difference of 0.06%.

**Mixed-species sign-off set (4 samples):** `num_cells_total`, `num_cells_Hsap`, `num_cells_Mmus`
and `num_multiplet_cells_total` were **identical — 0.000% difference — across all four samples**.
Per-species barnyard calling is exactly preserved.

**Earlier validation set (12 samples):** cell calls were bit-identical on 11 of 12 samples. On the
twelfth, one borderline human cell tipped across the calling threshold: `num_cells_total` 4383 vs
4384, `num_cells_Hsap` 2756 vs 2757, with `num_cells_Mmus` (1627) and multiplets (705) both
identical. That is one cell in 4384, or 0.02%.

That single cell is the honest worst case across the whole campaign, and it is the expected
behaviour of a threshold applied to a value that moved by ~0.001%: a barcode sitting on the
boundary can fall either side of it.

## Results: reads and counts

| Measure | Agreement |
|---------|-----------|
| Reads and counts, human sign-off set | within 0.02% |
| Reads and counts, mixed-species sign-off set | within 0.002% |
| `reads_after_deduplication`, earlier validation set | within 0.001–0.003% |
| Median genes per cell | within 0.02% |
| Cell-caller thresholds | unchanged |

The `.csv` metric outputs are unchanged in format, and the count matrices are produced in the same
`.h5ad` and tripartite (`barcodes` / `features` / `matrix`) forms as before.

## Results: performance

Measured on real CS Genetics samples on the same compute environment.

| Measure | 1.x | 2.0.0 | Change |
|---------|-----|-------|--------|
| Total compute (sum of process runtime) | 945 min | 579 min | **-39%** |
| Wall-clock, human dataset | 155 min | 92 min | **1.7x faster** |
| Wall-clock, mixed-species dataset | 675 min | 202 min | **3.3x faster** |

The largest single contributions:

| Process | Speed-up | Why |
|---------|----------|-----|
| `io_count` | ~154x | The counting pass ran on BusyBox awk, the production container's awk and by far the slowest; replaced with a compiled binary that was byte-identical on the standard validation references while intentionally correcting punctuated custom identifiers (known difference 5) |
| `umr_transcript_assignment` | 5.2x | Threaded BAM compression |
| `filter_for_multimappers_mismatch` | 3.1x | Threaded BAM compression |
| `dedup` | 2.5x | `umi_tools dedup` is single-threaded but position-local, so it is now split by contig and run in parallel |
| `initial_feature_count` | 1.9x | Threaded coordinate sort |

Two caveats on the wall-clock figures, in the interest of not overstating them: the four
comparison runs shared a compute environment and so competed for resources, and the 1.x
mixed-species run was further inflated by spot-instance reclaims during its very long `io_count`
stage. The compute-time reduction of 39% is the more robust number; treat the wall-clock ratios as
indicative of the improvement rather than as precise measurements.

## Known differences

Being explicit about what *did* change is more useful than claiming nothing did.

**1. A one-time shift in ambiguous multimapper assignment.** Where a read maps equally well to more
than one gene, 1.x kept whichever alignment happened to be processed last, and 2.0.0 keeps a
canonical one chosen deterministically. This affects roughly 0.3% of count-matrix entries, once,
on the transition from 1.x to 2.0.0. It is not a correction of an error — both assignments are
equally valid — but it is the reason the counts are not byte-identical. From 2.0.0 onward the
choice is stable, so this does not recur.

**2. Two RSeQC read-distribution bins.** `TSS_up_10kb` and `TES_down_10kb` differ by more than 2%
in both sign-off datasets. These are the outermost bins of the read-distribution plot and hold
only tens to hundreds of tags, so a small absolute change is a large relative one. Because the
baseline was determinism-controlled, these are *proven* to follow from the multimapper change
above rather than being run-to-run noise. No other metric in either dataset diverged by more
than 2%.

**3. Intermediate BAM files are not byte-identical.** The deduplicated BAM keeps a different — but
equivalent — representative read for some UMI groups, because parallel deduplication advances
`umi_tools`' internal state differently. Every count and metric derived from it is unchanged; only
the choice of which PCR duplicate represents the molecule differs. If you consume `dedup.bam`
directly rather than the count matrices, this is worth knowing.

**4. Exact sparse count arithmetic can correct low-order decimal digits.** During release
hardening, count-statistics reductions were changed from dense `float32` accumulation to validated
exact-integer sparse arithmetic. Metric definitions, CSV schemas, integer totals, and the values
shown in customer reports retain their two-decimal formatting. The raw decimal spelling of a mean
or percentage in `*.metrics.csv` can differ where the former `float32` accumulation rounded (for
example, an approximation of 30% becomes exactly `30.0`). On the representative validation
fixtures the correction was below report display precision; sufficiently high valid counts can
also correct the last displayed decimal digits. These are approved arithmetic corrections, so
cross-version metric files remain diagnostic `TEXT_EXACT` differences rather than being claimed
byte-identical. The validation and overflow contract is documented in
[Count-statistics arithmetic](count-statistics.md).

**5. Punctuated custom-reference gene identifiers are corrected.** The 1.x `io_count` awk
character class accepted only letters, digits, and underscores in an `XT:Z:` gene assignment, so
it silently truncated otherwise valid identifiers containing hyphens, colons, dots, spaces, or
UTF-8 bytes. The 2.0 extractor reads the complete SAM optional field and removes only a terminal
numeric version suffix, using the same normalization as GTF feature extraction. Standard
validation references were byte-identical for this step; affected custom references intentionally
receive the complete corrected identifier. Malformed or duplicate `XT` fields and distinct GTF
identifiers that would collide after version normalization fail loudly.

## Multi-lane and heavy-sample testing

Beyond equivalence, 2.0.0 was tested on deliberately demanding input: four lane-merged samples,
each combining all lanes of a full sequencing run.

This found two genuine bugs, both of which were **already present in 1.x**:

- A channel race in `filter_count_matrix` that surfaced on cloud filesystems and could put a
  threshold value where a file path was expected.
- A fan-out bug where the cell-caller threshold channel was parsed once per input CSV *row*, so a
  four-lane sample produced four duplicate entries and the run collided downstream.

Both caused the run to **fail loudly**; neither could produce silently incorrect results. The
practical consequence is that 1.x could not complete a heavy multi-lane run — an attempt to
process this dataset on 1.x for comparison purposes failed partway through on the first of these
bugs, producing no metrics at all. 2.0.0 completed the same data end to end.

Single-lane testing would not have found either problem, which is the main argument for having
done this.

## Checking equivalence yourself

The comparator used for this work ships with the pipeline:

```bash
# Strict: requires output equivalence with no count-matrix envelope.
# Use this to compare two 2.0.0 runs.
python tests/regression/compare_outputs.py <outdir_a> <outdir_b>

# Optional cross-version diagnostic: stage matching relative paths into two
# non-empty directories and retain the per-file details for review.
python tests/regression/compare_outputs.py <old_common> <new_common> \
  --allow-subset --envelope-max-flips 200 --json cross-version.json
```

Do not point the cross-version command at both complete output directories.
2.0 intentionally removes and replaces report filenames and changes execution
metadata, so whole-tree path equality across the major-version boundary is
neither expected nor what was validated here. For this release validation,
`old_common` and `new_common` preserved matching relative paths for metrics,
RSeQC, dedup and raw/filtered matrix outputs. This was a diagnostic inventory,
not a command expected to exit zero: the strict `TEXT_EXACT` metric/RSeQC
classes and strict dedup-count contract correctly reported the measured
cross-version changes described above. Each `DIFFER` was reconciled with those
results, and matrix envelope verdicts were assessed separately for the
net-preserving ambiguity. `--allow-subset` disables the 2.0 complete-run
`pipeline_info` manifest and the 2.0-only deterministic gzip-byte assertion;
it still validates and compares the legacy archives' decompressed payloads,
does not relax any other selected file, and does not turn a known value change
into a pass.

Strict mode is the default and applies each output class's documented equivalence contract with no
count-matrix ambiguity allowance. Envelope mode relaxes *only* the count-matrix comparison, and
only on a specific condition: per-barcode column sums must still be identical, and the number of
differing entries must fall within the budget you give it. Every other class of output stays
strict regardless. That condition is what makes it a meaningful test rather than a loosened one —
it permits a read to move between two genes, but not to appear, disappear, or move between
barcodes.

Choose the budget to suit your data; `--json` writes the full result set if you want to inspect
what differed. Within 2.0.0, strict mode should pass — the pipeline is now reproducible run to
run under those per-class semantic/output contracts. This does not assert whole-directory byte
identity: for example, the deduplicated-BAM contract deliberately permits a different representative
read for the same stable molecule, while STAR, initial featureCounts and high-confidence annotated
BAMs retain their full alignment content; run-specific Nextflow task IDs/timings/resources are
normalised.
MultiQC timestamps, work/temp directories and equivalent provenance paths are likewise
normalised while its report data remains part of the comparison.
The raw and filtered tripartite count-matrix archives have the stronger guarantee: identical
inputs produce byte-identical `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` files.
The default complete-tree comparator validates each decompressed payload and then enforces those
archive bytes; the explicit cross-version `--allow-subset` diagnostic compares legacy decompressed
semantics because 1.x did not have the deterministic gzip-header contract.
