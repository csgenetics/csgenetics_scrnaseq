# 05 - Report Generation: External vs Internal

Scope: make the EXTERNAL report "look much better" WITHOUT changing what metrics/sections/outputs are reported.

## 1. EXTERNAL: how reports are currently built

### 1.1 Process chain (modules/processes.nf, wired in main.nf)
Per-sample, then experiment-wide:

1. `single_sample_multiqc` (processes.nf:532) - MultiQC (`-m unified_qc -m rseqc`) -> `${sample_id}_multiqc.html` + `multiqc_data.json`. Container `quay.io/csgenetics/multiqc:0.1` (images.config:30). The JSON feeds `summary_statistics`; the HTML is published standalone.
2. `multi_sample_multiqc` (processes.nf:559) - same, experiment-wide -> `multisample_multiqc.html`. Container `multiqc:0.1` (images.config:31).
3. `cell_caller` (processes.nf:705) - `cell_caller.py` emits Plotly HTML fragments: `${sample_id}_counts_pdf_with_threshold.html`, `${sample_id}_barnyard_plot.html`. Container `scanpy_anndata:0.0.4`.
4. `summary_statistics` (processes.nf:790) - `summary_statistics.py` -> `${sample_id}.metrics.csv` (THE metrics contract). Container `scanpy_anndata:0.0.4`.
5. `qc_cascade_plot_single` / `qc_cascade_plot_multi` (processes.nf:820, 843) - `qc_cascade_plot.py` -> `${sample_id}.qc_cascade.html` / `multisample_qc_cascade.html`. Container `html_build:0.1.0` (images.config:37-38).
6. `single_summary_report` (processes.nf:863) - `create_single_sample_report.py` -> `${sample_id}_report.html`. Container `html_build:0.1.0`.
7. `multi_sample_report` (processes.nf:885) - `create_multi_sample_report.py` -> `multisample_report.html`, `multisample_out.csv`, `multisample_summary_plots.html`, `multisample_qc_cascade.html`. Container `html_build:0.1.0`.

Wiring: main.nf:337-470. Single report inputs = metrics_csv + cell-caller pdf html + barnyard html + qc_cascade html + `single_sample_report_template.html.jinja2`. Multi report consumes the collected per-sample `.metrics.csv` files.

### 1.2 Rendering engine (the customer report)
- Pure **Jinja2 `Template`** string render in Python (`create_single_sample_report.py:57,119`, `create_multi_sample_report.py:28,196`). No build step, no bundler, no JS modules.
- Templates: `templates/single_sample_report_template.html.jinja2` (233 lines) and `templates/multi_sample_report_template.html.jinja2` (~72KB). The large size is from inline Plotly HTML fragments embedded via `{{ plot | safe }}` at render time, not hand-written.
- Styling: **Bootstrap 5.3.0 + bootstrap-icons loaded from jsdelivr CDN** (`<link>` in template head) plus ~10 lines of inline `<style>`. Tooltips/accordions use native Bootstrap `data-bs-*` (no custom JS file).
- Plots are **standalone Plotly HTML fragments** generated independently in Python (`qc_cascade_plot.py`, `cell_caller.py`, `create_multi_sample_report.py`) with `pio.to_html(full_html=False, include_plotlyjs='cdn')` / `fig.write_html(include_plotlyjs='cdn')` (qc_cascade_plot.py:216,522; cell_caller.py:73). They are read back as strings and string-substituted into the Jinja template (`read_html_plot()`).
- Consequence: **Plotly.js is loaded from CDN** (`include_plotlyjs='cdn'`), so the report is NOT self-contained / not offline-safe. Each embedded fragment may pull its own plotly.

### 1.3 Sections / layout of the single-sample report
(template single_sample_report_template.html.jinja2)
- Header: sample_id (line 39-40).
- **Headline stat cards** (line 49-69): top row of key metrics with info-circle tooltips.
- **Card per group**: `Read QC` and `Deduplication` cards iterate `metrics_dict[group]` (line 84-101).
- **Cell Caller plot** card: `{{ pdf_plot | safe }}` (line 117), hidden if plot file size 0 (`show_cell_caller_plot`).
- **QC Cascade** card: `{{ qc_cascade_plot | safe }}` (line 134), gated on `show_qc_cascade_plot`.
- **Cell metrics** card (line 149): flat list for single-species; for mixed-species an **accordion** (`accordion-flush`, line 155-191) per sub-category driven by `cell_stat_cat_dict` (header/collapse IDs from `get_cell_stat_cat_dict_obj`, create_single_sample_report.py:9-29).
- **Barnyard plot** card (line 213): mixed-species only, gated on `show_barnyard_plot`.
- Alignment cards (`Post read QC alignment`, `Annotated reads alignment`) are coded but **temporarily disabled** by removing HTML; the Python `alignment_cat_dict` is still passed (create_single_sample_report.py:131-136).

### 1.4 The metrics contract (MUST be preserved unchanged)
`summary_statistics.py` writes `${sample_id}.metrics.csv` with header
`variable_name,value,human_readable_name,description,classification` (summary_statistics.py:170).
Five `classification` groups (the report sections): **Read QC, Cell metrics, Deduplication, Post read QC alignment, Annotated reads alignment** (summary_statistics.py:36-40). Each row carries its own human-readable name + tooltip description. Number formatting: floats -> 2 d.p., ints unchanged (`format_number_to_string`).

Multi-sample report: `create_multi_sample_report.py` reads every `.metrics.csv`, transposes to metric->list-over-samples, writes `multisample_out.csv`, and renders violin/box summary plots for 4 key metrics (reads_pre_qc, num_cells, raw_reads_per_cell, median_genes_detected_per_cell; mixed-species variants at lines 49-58) via `plotly.graph_objects` with CSG colors.

### 1.5 Published outputs (the output-equivalence baseline)
Under `${params.outdir}/`:
- `report/${sample_id}/${sample_id}_report.html`, `${sample_id}.metrics.csv`, `${sample_id}.qc_cascade.html`
- `report/multisample_report.html`, `multisample_out.csv`, `multisample_summary_plots.html`, `multisample_qc_cascade.html`
- `multiqc/single_sample_multiqc/${sample_id}/...`, `multiqc/multisample_multiqc.html` + data
- `plots/${sample_id}*_pdf_with_cutoff.html`
Any rework must keep the same set of metric values, the `.csv` schemas, and ideally the same published filenames.

## 2. INTERNAL: how the report is built (architecture only)

Module `modules/pipeline_report/` (main.nf + README.md + skills). It deliberately **replaces** MultiQC + scattered plot HTMLs with **one self-contained interactive report** (`*.pipeline_report.html`).

### 2.1 Three-layer architecture (skill: pipeline-report-ui)
**Python (raw data extraction) -> Jinja2 (HTML structure + data embedding) -> JavaScript (Plotly trace building + render).**
- Python (`report_generator/` package: `sections/`, `plotters/`, `domain/`, `utils/`, `presentation/`) extracts **raw native-type dicts only** - it must NOT build Plotly figure objects (avoids the numpy typed-array JSON serialization bug). Entry point `pipeline_report_generator.py`.
- One large Jinja2 template `templates/pipeline_report_template.html.jinja2` (~208KB). Its `<script>` block only embeds `{{ data|tojson|safe }}` + orchestrator registration - **no logic inline**.
- JavaScript lives as **ES modules** in `resources/usr/bin/js/src/` (66 files across ~17 section dirs: qc-cascade, mapping, headline-metrics, performance-summary, read-characterization, cell-caller, interp-extrap, barcode-diversity, sss-distance, pbmc-cell-types, plus `shared/`, `utils/`, `orchestration/`, `tab-navigation/`, `metric-visibility/`, `experiment-layout/`).

### 2.2 Build system
- **esbuild** bundles `src/index.js` -> single minified IIFE `report-bundle.js` (globalName `PipelineReport`), built **at container runtime** inside the process (main.nf:70-72 copies `js/` out of staging, `npm install --omit=dev --ignore-scripts` then `npm run build`). Python injects the bundle into the template.
- **Two-directory split** (skill pipeline-report-js-dev): production source in `resources/usr/bin/js/` (staged by Nextflow); dev/test tooling (vitest + jsdom) in `modules/pipeline_report/js/` so `node_modules/` never enters the Nextflow staging zone (Nextflow 1MB moduleBinaries limit).
- Tests: Vitest, `cd modules/pipeline_report/js && npm test`, importing source via `@src/` alias.

### 2.3 Shared infra worth knowing
- `js/src/shared/`: `trace-builders.js` (createBoxTrace/LineTrace/BarTrace/SEMBandTrace/LegendOnlyTrace), `layout-builders.js` (createBaseLayout/BoxPlotLayout/... + addHorizontal/VerticalLine), `color-utils.js` (CS_GREEN `#36BA00`, semantic colors, ALPHABET_COLORS, getColorPalette, colorWithAlpha), `selection-tree.js`, `toggle-utils.js`.
- Orchestrator pattern: `orchestrator.registerModule(name, initFn, {priority})` (10 tooltips, 50 section plots, 100 main); no standalone DOMContentLoaded.
- Manager-class-per-section pattern with lazy rendering, condition/units toggles, `Plotly.react()` updates.
- Self-contained output: Plotly.js + bundle embedded, no CDN/external deps; Bootstrap 5, Lexend font, `simple_white` template, fixed constants shared Python<->JS.

### 2.4 Key differences from external
| Aspect | External | Internal |
|---|---|---|
| Engine | Jinja2 string render only | Python raw-data -> Jinja2 -> JS ES modules + esbuild |
| Plots | standalone Plotly HTML fragments string-pasted | JS builds Plotly from embedded raw data |
| Plotly.js | CDN (`include_plotlyjs='cdn'`) | embedded, self-contained |
| Interactivity | Bootstrap accordions/tooltips only | toggles, dropdowns, sample trees, tabs, lazy render |
| Build/test | none | esbuild bundle + Vitest |
| Output files | many separate HTMLs (per-sample, multi, multiqc) | one consolidated `pipeline_report.html` |
| Scope of metrics | fixed 5-group metrics.csv | targets/interp-extrap/much richer (DIFFERENT functionality) |

NOTE: internal reports DIFFERENT metrics (targets, interp/extrap, PBMC, barcode diversity, etc.). Adopting its *engine* is allowed; adopting its *metrics/sections* is NOT (would change what external computes).

## 3. SCOPED OPTIONS to improve EXTERNAL report appearance

Output-equivalence anchor: same metric values, same `.metrics.csv`/`multisample_out.csv` schemas, same published `.html` filenames. Appearance/markup may change.

### Option A - CSS / template polish, same engine (LOW effort, LOW risk)
- Edit only the two `.jinja2` templates + inline `<style>`: spacing, typography, card styling, section headers, sticky nav, color scheme (apply CS green `#36BA00`, Lexend font), tidy headline cards, re-enable/restyle the disabled alignment cards if desired.
- No Python change, no container change, metric values untouched.
- Risk: cosmetic only. Verify Jinja variable references unchanged. Output HTML differs textually but reports identical numbers.
- Effort: ~0.5-1 day.

### Option B - Self-host assets + Plotly theming (LOW-MED effort, LOW-MED risk)
- Switch `include_plotlyjs='cdn'` -> `'inline'` (or a pinned local) in `qc_cascade_plot.py`, `cell_caller.py`, `create_multi_sample_report.py` so reports are offline/self-contained (matches internal). Vendor Bootstrap/icons locally instead of jsdelivr CDN.
- Apply a shared Plotly layout theme (Lexend, simple_white, CS colors, consistent margins/legends) across all three plot scripts - tightens look without changing data.
- Risk: must keep file sizes acceptable; verify plots still render. Container `html_build:0.1.0` / `scanpy_anndata` may need plotly assets - check before committing.
- Effort: ~1-2 days.

### Option C - Consolidate to a single styled report, keep Jinja-only engine (MED effort, MED risk)
- Merge per-sample cards + plots + a navigable layout into one cleaner template with sticky section nav and Bootstrap tabs/accordions (borrow internal's *visual* section-header pattern, not its JS). Possibly fold multi-sample summary + qc-cascade into a consolidated experiment page.
- Still pure Jinja2 string render; plots still pre-generated fragments.
- Risk: must preserve all published output files (some consumers may expect `${sample_id}_report.html` etc.). If consolidating, keep old files OR confirm no downstream dependency. Touches `create_*_report.py` orchestration.
- Effort: ~3-5 days.

### Option D - Adopt internal's Python->Jinja2->JS + esbuild engine (HIGH effort, MED-HIGH risk)
- Port external's 5-group metrics + cell-caller/qc-cascade/barnyard/violin plots onto internal's architecture: raw-data extraction in Python, ES-module trace builders reusing `shared/` (trace/layout/color utils), esbuild bundle, Vitest tests, single self-contained `report.html`.
- Pros: best-looking, interactive, offline, testable, shared with internal patterns; reuses mature infra.
- Cons: large rebuild; new container/build step (npm+esbuild at runtime) and two-dir JS layout; must re-derive identical metric values and decide output-file compatibility. Highest chance of accidental output drift.
- Risk to output-equivalence: MED-HIGH unless the `.metrics.csv` producers (`summary_statistics.py`) stay untouched and only presentation changes. Keep `summary_statistics.py` / `categorize_reads.py` as the single source of metric truth.
- Effort: ~2-4 weeks.

### Recommendation framing
A+B together are the cheap, safe win (modern look, self-contained, CS branding) with essentially zero output-equivalence risk. C is the natural mid-point if a single consolidated report is desired. D is only justified if interactivity/maintainability parity with internal is a hard requirement; gate it behind a strict "metrics.csv is frozen" contract test.

## 4. Output-equivalence guardrails for any option
- Freeze `summary_statistics.py` metric computation and the `variable_name,value,human_readable_name,description,classification` CSV schema; add a golden-file test on `.metrics.csv` and `multisample_out.csv`.
- Keep published filenames under `${params.outdir}/report/...` unless a downstream-consumer audit clears renaming.
- Treat MultiQC outputs (`single_sample_multiqc`, `multi_sample_multiqc`) as separate; they feed metrics via `multiqc_data.json` and are also published HTML - changing the customer report does not require touching MultiQC.
