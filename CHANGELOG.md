# Changelog

## 2.0.0

Major overhaul: Nextflow 26 migration, performance, per-process module structure, and a
consolidated report. Pipeline metric values and `.csv` outputs are unchanged; the changes
below are runtime, structural, and presentational.

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

### Reproducibility note

The pipeline has a small, pre-existing, inherent run-to-run non-determinism in the count
matrix: a few genuinely-ambiguous multimapped reads can be attributed to either of two genes.
Each barcode's TOTAL counts and all summary metrics are unaffected (the difference is
net-preserving). This behaviour is unchanged by this release; the regression comparator's
envelope mode accounts for it while holding metrics and per-barcode totals byte-exact.
