# Changelog

## 2.0.0

Major overhaul: Nextflow 26 migration, performance, per-process module structure, a
consolidated report, and run-to-run reproducibility. Summary metrics and per-barcode total
counts are byte-identical to before; the changes below are runtime, structural, presentational,
and a one-time deterministic resolution of a pre-existing count-matrix ambiguity (see
Reproducibility).

### Breaking

- **Requires Nextflow `>= 26.04.0`** (manifest `nextflowVersion = '!>=26.04.0'`). Nextflow 26
  evaluates process directive strings eagerly at compile time and uses the v2 config parser;
  the pipeline relies on these and will not run on older Nextflow.
- **The HTML report is consolidated into a single file.** The separate per-sample and
  multi-sample reports are replaced by one `report/consolidated_report.html`:
  - removed: `report/<sample>/<sample>_report.html`, `report/multisample_report.html`,
    `report/multisample_summary_plots.html`, `report/multisample_qc_cascade.html`
  - added: `report/consolidated_report.html`
  - unchanged: per-sample `report/<sample>/<sample>.metrics.csv`, `report/multisample_out.csv`,
    `plots/*.html`, and all MultiQC outputs.

  If you consume the old report HTML filenames, switch to `consolidated_report.html`. The
  metric VALUES and the `.csv` outputs are identical to before.

### Changed

- **Nextflow 26 migration.** `publishDir`/`pattern` directives that interpolate input
  variables are closure-wrapped; the redundant `nextflow.enable.strict` flag is removed.
- **Performance.** STAR (`--runThreadN`) and featureCounts (`-T`) now use the cores allocated
  to the task (`${task.cpus}`) instead of a hardcoded count; `dedup` is right-sized to 1 cpu
  (`umi_tools` is single-threaded). This gives roughly 2x faster alignment on production-scale
  samples (no measurable change on tiny, index-load-bound inputs). Pipeline outputs are
  unchanged (verified by an output-equivalence comparator).
- **Report.** A single, offline-safe (no CDN) consolidated report with CS Genetics branding,
  a searchable per-sample selector, and a cross-sample metrics table.
- **Structure.** The monolithic `modules/processes.nf` is split into per-process modules at
  `modules/local/<name>/main.nf` (nf-core local-module layout).

### Added

- An nf-test harness (`nf-test.config` + `tests/nf/`) with a `stub:` block on every process
  and a whole-pipeline `-stub` DAG test that exercises the full workflow wiring.
- An output-equivalence regression comparator (`tests/regression/compare_outputs.py`) with a
  strict mode and an envelope mode (see note below).

### Reproducibility

This release makes the pipeline **run-to-run byte-reproducible** by pinning `PYTHONHASHSEED=0`.
The prior non-determinism was `umi_tools` choosing between equally-ranked reads via Python's
hash-seed-randomized set iteration; pinning the seed makes that choice deterministic.

One-time consequence: a few genuinely-ambiguous multimapped reads (which previously landed on
either of two genes at random, run-to-run) now resolve deterministically. **Summary metrics and
per-barcode total counts are byte-identical to before**; only those few per-gene count-matrix
entries change, once. From this release, count matrices are reproducible — the output-equivalence
comparator (`tests/regression/compare_outputs.py`) can therefore gate future changes byte-exact.
