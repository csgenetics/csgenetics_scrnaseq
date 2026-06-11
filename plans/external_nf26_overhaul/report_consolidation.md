# Report consolidation (Full Option C) - first cut

Goal: merge the separate per-sample reports and the multi-sample report family
into ONE navigable, styled, offline-safe HTML report, WITHOUT changing any
metric VALUES. Engine stays Jinja2. First cut for human visual review.

## 1. Design

### Wiring choice: replace the RENDERING layer only

Metric/plot producers untouched in computation; only rendering is replaced.

- KEPT (computation frozen): `summary_statistics.py`, `categorize_reads.py`, the
  `${sample_id}.metrics.csv` schema/values; the `cell_caller`,
  `summary_statistics`, `qc_cascade_plot_single`, `qc_cascade_plot_multi`
  processes.
- REMOVED: `single_summary_report`, `multi_sample_report`, and their templates.
- ADDED: one new process `consolidated_report` consuming ALL per-sample
  `metrics.csv` + all per-sample plot fragments + the multi-sample qc-cascade
  fragment, emitting ONE `consolidated_report.html`, and re-emitting
  `multisample_out.csv` (byte-identical to the old output).

`consolidated_report` takes inputs flat (collected file lists), discovers
samples from `*.metrics.csv` filenames, and locates each sample's plot fragments
by their conventional names in the same staging directory.

### Offline-safe / self-contained

Old reports loaded Bootstrap + icons from jsdelivr and Plotly.js from CDN. The
new report embeds everything inline:

- Plotly.js embedded ONCE via `plotly.offline.get_plotlyjs()`. Each plot fragment
  is now produced WITHOUT its own Plotly.js (`include_plotlyjs=False`).
- Bootstrap 5.3.0 CSS + JS bundle vendored locally and inlined.
- bootstrap-icons 1.5.0 + Lexend font vendored, with the font files inlined as
  base64 inside the CSS (`*.embedded.css`) - no external font loads.

Verified: rendered report has no external `<link>`, `<script src>`, `<img src>`
loads, no jsdelivr, no Google Fonts. Remaining `http(s)://` strings are SVG/XML
namespaces and inert doc/license/map-tile URLs baked inside the Plotly.js and
Bootstrap sources (never fetched; this report renders no maps).

NOTE: standalone published cell-caller plots
(`plots/${sample_id}*_pdf_with_cutoff.html`) are unchanged - they keep
`include_plotlyjs='cdn'` so they remain self-viewable on their own. Only the
report-bound fragments switched to `include_plotlyjs=False`.

### Layout / branding

CS green `#36BA00`, Lexend font, Bootstrap 5 cards. Sticky top navbar with the CS
logo + anchor links (Overview / Per-sample / Experiment summary / Cross-sample
table). Overview header: sample count + headline cards. Per-sample view: one
Bootstrap tab per sample reproducing exactly the old per-sample content (headline
cards, Read QC / Deduplication tables, Cell metrics flat-or-accordion,
cell-caller plot, qc-cascade plot, barnyard for mixed). Experiment summary:
across-sample violin/box plots + experiment-wide qc-cascade. Cross-sample table:
same content as `multisample_out.csv`, sticky header/first column, column headers
link back to that sample's tab.

## 2. OLD -> NEW published-filename mapping (CHANGELOG / customer deprecation)

Under `${params.outdir}/report/`:

| OLD file (removed)                                 | NEW file                          |
|----------------------------------------------------|-----------------------------------|
| `${sample_id}/${sample_id}_report.html`            | `consolidated_report.html`        |
| `multisample_report.html`                          | `consolidated_report.html`        |
| `multisample_summary_plots.html`                   | `consolidated_report.html`        |
| `multisample_qc_cascade.html` (standalone publish) | folded into `consolidated_report.html` |

PRESERVED (unchanged location, schema and values):

| File                                               | Status                            |
|----------------------------------------------------|-----------------------------------|
| `report/${sample_id}/${sample_id}.metrics.csv`     | unchanged (published by `summary_statistics`) |
| `report/multisample_out.csv`                       | unchanged schema + values (re-emitted by `consolidated_report`) |
| `report/${sample_id}/${sample_id}.qc_cascade.html` | still published (now a Plotly-free fragment) |
| `plots/${sample_id}*_pdf_with_cutoff.html`         | unchanged (standalone, still CDN Plotly) |
| `multiqc/...` MultiQC outputs                       | untouched                         |

Customer note: per-sample `${sample_id}_report.html` and the `multisample_*`
report HTMLs are replaced by a single `consolidated_report.html`. The
machine-readable `multisample_out.csv` and per-sample `metrics.csv` are unchanged.

## 3. Files changed / added / removed

Added:
- `templates/consolidated_report_template.html.jinja2`
- `bin/create_consolidated_report.py` (reuses `get_cell_stat_cat_dict_obj` and
  `MultipleSampleSummaries`; no metric recomputation)
- `modules/local/consolidated_report/main.nf`
- `assets/vendor/`: `bootstrap.min.css`, `bootstrap.bundle.min.js`,
  `bootstrap-icons.embedded.css`, `lexend.embedded.css`, `cs_logo_data_uri.txt`

Changed (computation untouched, only plot-bundling / wiring):
- `bin/cell_caller.py` - `output_plot_to_html` gained `include_plotlyjs` param
  (default `'cdn'`); report-bound cell-caller + barnyard fragments pass `False`.
- `bin/qc_cascade_plot.py` - single + multi qc-cascade write Plotly-free
  fragments via new `write_plotly_fragment`.
- `main.nf` - dropped the two old report includes/wiring; added
  `consolidated_report`; `qc_cascade_plot_multi` now consumes
  `summary_statistics.out.metrics_csv`; added consolidated template + vendor path.
- `conf/images.config`, `conf/conda_envs.config`, `conf/base.config` - replaced
  the two old report labels with `consolidated_report` (html_build container/conda).

Removed:
- `modules/local/single_summary_report/main.nf`
- `modules/local/multi_sample_report/main.nf`
- `templates/single_sample_report_template.html.jinja2`
- `templates/multi_sample_report_template.html.jinja2`

Retained but no longer wired as entrypoints (imported for shared helpers):
- `bin/create_single_sample_report.py`, `bin/create_multi_sample_report.py`

## 4. Verification

- `NXF_VER=26.04.1 nextflow inspect main.nf -profile test` -> exit 0.
- `NXF_VER=26.04.1 nf-test test tests/nf/pipeline_stub.nf.test` -> PASSED.
- Mock render harness (real Sample1/Sample2 metrics.csv + dummy fragments):
  `multisample_out.csv` byte-identical to pre-change output; examples at
  `/tmp/consolidated_report_example.html` (single) and
  `/tmp/consolidated_report_example_mixed.html` (mixed); no external
  CDN/jsdelivr/Google-Font/resource-load references in output.
- `bin/summary_statistics.py` and `bin/categorize_reads.py` unchanged.

## 5. Open design questions for the human

1. Standalone cell-caller plots still use CDN Plotly - leave, or make offline too?
2. Cross-sample column header links jump to the sample tab via small JS helper - ok?
3. Per-sample selector is Bootstrap tabs; for many samples a searchable dropdown
   may scale better - now or later?
4. Report is ~5.8 MB (one inlined Plotly.js copy) - acceptable, or slim the bundle?
5. The disabled alignment cards remain disabled (matching today) - re-enable now?
