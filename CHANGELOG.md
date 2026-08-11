# Changelog

## 2.0.0 - 2026-08-11

Major overhaul: Nextflow 26 migration, performance, per-process module structure, a
consolidated report, and run-to-run reproducibility. Summary metrics and per-barcode total
counts are byte-identical to before; the changes below are runtime, structural, presentational,
and a one-time deterministic resolution of a pre-existing count-matrix ambiguity (see
Reproducibility).

Equivalence with 1.x was verified on 20 samples across four datasets covering both supported
genome types: cell calling was identical or differed by a single borderline cell, and read and
count totals agreed to within 0.02%. See [docs/validation.md](docs/validation.md) for the method,
the full results, and the known differences.

### Breaking

- **Requires Nextflow `>= 26.04.0`** (manifest `nextflowVersion = '!>=26.04.0'`). Nextflow 26
  evaluates process directive strings eagerly at compile time and uses the v2 config parser;
  the pipeline relies on these and will not run on older Nextflow.
- **The HTML report is consolidated into a single file.** The separate per-sample and
  multi-sample reports are replaced by one `report/consolidated_report.html`:
  - removed: `report/<sample>/<sample>_report.html`, `report/multisample_report.html`,
    `report/multisample_summary_plots.html`
  - added: `report/consolidated_report.html`
  - unchanged: per-sample `report/<sample>/<sample>.metrics.csv` and
    `report/<sample>/<sample>.qc_cascade.html`, `report/multisample_out.csv`,
    `report/multisample_qc_cascade.html`, `plots/*.html`, and all MultiQC outputs.

  If you consume the old report HTML filenames, switch to `consolidated_report.html`. The
  metric VALUES and the `.csv` outputs are identical to before.

### Changed

- **Nextflow 26 migration.** `publishDir`/`pattern` directives that interpolate input
  variables are closure-wrapped; the redundant `nextflow.enable.strict` flag is removed.
- **Performance (cluster throughput, not single-task speed).** Per-process cpu reservations were
  right-sized to *measured* core usage from profiling real CS Genetics samples on Seqera (human
  GRCh38 and the 60 GB mouse_human_mix barnyard index). The finding: STAR (~3-7 cores, memory-bound),
  featureCounts (~1 core), and the `io_count` awk pass (~0.4 cores, I/O-bound) do not saturate the
  cores previously reserved for them, and adding threads gives no measurable speedup at our read
  counts. Reservations were reduced accordingly (e.g. `star` 16->8, `io_count` 4->1, featureCounts
  steps 4->2) so more samples pack onto each instance, and `dedup` is set to 1 cpu (`umi_tools` is
  single-threaded). STAR's `--runThreadN` is kept at the original value of 8 and pinned (decoupled
  from the cpu reservation) because the thread count affects multimapper output order and therefore
  per-barcode counts; only the reservation changed. These are throughput/packing changes; pipeline
  outputs are byte-unchanged (verified by the output-equivalence comparator).
- **Performance (process-level speedups).** Profiling real samples on Seqera identified the dominant
  single-threaded steps and rewrote/parallelised them. Verified on 8 real human samples (cell calls
  identical to before, read counts within ~0.003%); total compute dropped ~39% (945 -> 579 CPU-min):
  - `io_count`: the awk pass ran on BusyBox awk (the production container's awk) and was the single
    largest cost. Replaced with a static Rust binary (`bin/io_count_extract`, byte-identical) -> ~150x.
  - `dedup`: `umi_tools dedup` is single-threaded but position-local, so it is now split by reference
    contig and run in parallel (`bin/dedup_by_contig.sh`), counts identical -> ~2.5x.
  - `multimapper_transcript_assignment`/`multimapper_exon_assignment`: the multimapper assignment is
    made deterministic and order-independent (canonical representative per gene in
    `bin/assign_multi_mappers.gawk`), which allows the previously single-threaded `samtools sort -n`
    to be threaded. One-time effect: a few genuinely-ambiguous multimapped reads (~0.3% of count-matrix
    entries; cell calls unchanged) resolve to a deterministic canonical alignment.
  - `samtools sort`/`view` threading in `initial_feature_count` and the filter steps (byte-identical).
  Remaining target: the multimapper assignment gawk and RSeQC `read_distribution` (both single-threaded).
- **Report.** A single, offline-safe (no CDN) consolidated report with CS Genetics branding,
  a searchable per-sample selector, and a cross-sample metrics table. Independently published
  Cell Caller and single-/multi-sample QC cascade plots are offline-safe too, and Seqera report
  mappings now follow the 2.0 output names with the cross-sample CSV available as a lightweight
  fallback for large HTML reports.
- **Structure.** The monolithic `modules/processes.nf` is split into per-process modules at
  `modules/local/<name>/main.nf` (nf-core local-module layout).

### Fixed

Three pre-existing bugs, all present in 1.x, found by heavy multi-lane testing. **All three caused
the run to fail loudly; none could produce silently incorrect results.**

- **Multi-lane runs failed.** The cell-caller threshold channel was parsed once per input-CSV
  *row*, but a multi-lane sample occupies one row per lane. A four-lane sample therefore produced
  four duplicate threshold entries, fanning `cell_caller` out four times and colliding downstream
  in `qc_cascade_plot_multi` with an input file name collision. Fixed by deduplicating the
  threshold channel; a no-op for single-lane input, which is why single-lane testing never
  surfaced it.
- **`filter_count_matrix` channel race.** On cloud filesystems (Seqera/Fusion), arrival-order
  dependent `groupTuple` ordering could place the cell-caller threshold in the slot expected to
  hold a file path, failing with `Not a valid path value`. Replaced with a deterministic
  `join(by: 0)`.
- **Empty-sample deduplication log.** Wrote a literal `\n` rather than a newline.

### Added

- An nf-test harness (`nf-test.config` + `tests/nf/`) with a `stub:` block on every process
  and a whole-pipeline `-stub` DAG test that exercises the full workflow wiring.
- An output-equivalence regression comparator (`tests/regression/compare_outputs.py`) with a
  strict mode and an envelope mode (see note below).
- [`docs/validation.md`](docs/validation.md) — how 2.0.0 was validated against 1.x, with the
  datasets, results, and known differences.

### Reproducibility

This release makes the pipeline **run-to-run byte-reproducible** by pinning `PYTHONHASHSEED=0`.
The prior non-determinism was `umi_tools` choosing between equally-ranked reads via Python's
hash-seed-randomized set iteration; pinning the seed makes that choice deterministic.

One-time consequence: a few genuinely-ambiguous multimapped reads (which previously landed on
either of two genes at random, run-to-run) now resolve deterministically. **Summary metrics and
per-barcode total counts are byte-identical to before**; only those few per-gene count-matrix
entries change, once. From this release, count matrices are reproducible — the output-equivalence
comparator (`tests/regression/compare_outputs.py`) can therefore gate future changes byte-exact.

**What this means in practice.** Under 1.x, processing the same FASTQs twice could give slightly
different count-matrix entries. If you are comparing samples processed at different times — a
longitudinal study, a re-analysis of archived data, or a batch split across several runs — that
variation was a floor on how precisely results could be compared. From 2.0.0 it is gone: the same
input produces the same output. Note that the floor applies to the 1.x-to-2.0.0 transition itself,
so a mixed-version comparison still carries the one-time shift described above; re-processing the
older data under 2.0.0 removes it.
