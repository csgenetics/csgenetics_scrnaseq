#!/usr/bin/env python

"""
Create ONE consolidated, self-contained HTML report for a CS Genetics
scRNA-seq experiment.

This replaces the previous per-sample reports (``${sample_id}_report.html``)
and the multi-sample report family (``multisample_report.html``,
``multisample_summary_plots.html``) with a single navigable, styled report:

  * an experiment overview header (sample count + headline metrics),
  * a per-sample view (selectable via Bootstrap tabs) reproducing everything
    the old per-sample report showed: headline cards, the Read QC /
    Deduplication / Cell metrics groups, the cell-caller plot, the qc-cascade
    plot, and (mixed-species) the barnyard plot,
  * the experiment-wide multi-sample summary plots (violin/box) and the
    experiment-wide qc-cascade plot,
  * the cross-sample metrics table (same content as multisample_out.csv).

HARD CONSTRAINT: this script does NOT change any metric VALUES. It reads the
SAME ``${sample_id}.metrics.csv`` files produced by summary_statistics.py and
renders them. Number formatting matches the old reports exactly
(floats -> 2 d.p., ints unchanged). It also re-emits ``multisample_out.csv``
with the identical schema and values as before.

OFFLINE-SAFE: Plotly.js, Bootstrap CSS/JS, Bootstrap-icons and the Lexend font
are embedded inline from locally vendored assets. The rendered report contains
no external CDN URLs and renders with no network.

Usage:
  create_consolidated_report.py <template> <mixed_species> <vendor_dir> \\
      <multi_qc_cascade_html>

All per-sample inputs are staged FLAT into the working directory by Nextflow.
Samples are discovered from the ``*.metrics.csv`` files, and each sample's plot
fragments are located by their conventional filenames in the same directory:
  - <sample_id>.metrics.csv
  - <sample_id>_counts_pdf_with_threshold.html   (cell-caller fragment)
  - <sample_id>_barnyard_plot.html               (barnyard fragment; mixed only)
  - <sample_id>.qc_cascade.html                  (qc-cascade fragment)
"""

import sys
import os
from pathlib import Path
from collections import defaultdict

from jinja2 import Template
from plotly.offline import get_plotlyjs

from create_single_sample_report import get_cell_stat_cat_dict_obj


def format_number_to_string(number_str):
    """
    Format numbers exactly as the previous reports did:
    floats -> 2 d.p.; ints unchanged. Frozen formatting contract.
    Used ONLY for the multisample_out.csv data file (must stay byte-identical;
    thousands separators would corrupt the comma-delimited file).
    """
    if "." in number_str:
        return f"{float(number_str):.2f}"
    else:
        return number_str


def format_number_for_display(number_str):
    """
    Human-readable DISPLAY formatting for the on-screen tables and headline cards:
    add thousands separators to the integer part and keep the existing decimals
    (floats stay at 2 d.p., matching the value contract -- no rounding, so
    percentages/rates/fractions are unchanged; e.g. 181004808 -> "181,004,808",
    12106.54 -> "12,106.54", 94.30 -> "94.30"). DISPLAY ONLY -- never used for the
    CSV. Non-numeric sentinels (e.g. "nan", "nan_nan") pass through unchanged.
    """
    s = number_str.strip()
    try:
        if "." in s:
            return f"{float(s):,.2f}"
        return f"{int(s):,}"
    except ValueError:
        return number_str


def read_fragment(path):
    """Read a plotly HTML fragment from disk, or return None if absent/empty."""
    if path is None or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    with open(path, "r", encoding="utf-8") as f:
        return f.read()


def read_text_asset(vendor_dir, name):
    with open(os.path.join(vendor_dir, name), "r", encoding="utf-8") as f:
        return f.read()


class ConsolidatedReport:
    def __init__(self):
        self.template_path = sys.argv[1]
        self.mixed = sys.argv[2].upper() == "TRUE"
        self.vendor_dir = sys.argv[3]
        self.multi_qc_cascade_path = sys.argv[4]
        # All per-sample inputs are staged flat into the cwd.
        self.work_dir = "."

        with open(self.template_path) as fh:
            self.jinja_template = Template(fh.read())

        # Discover samples by their metrics csv. Sort for a stable, logical order.
        self.sample_ids = self._discover_sample_ids()

        # Build per-sample data structures (cards + plot fragments).
        self.samples = [self._build_sample(sid) for sid in self.sample_ids]

        # Build experiment-wide structures (cross-sample table + summary plots),
        # re-emitting multisample_out.csv with the unchanged schema/values.
        self.multi_metrics_dict, self.cell_metric_tooltip_dict = self._build_multi_metrics_and_csv()
        self.summary_plot_fragment = self._build_summary_plot_fragment()
        self.multi_qc_cascade_fragment = read_fragment(self.multi_qc_cascade_path)

        # Headline experiment metrics for the overview header.
        self.overview = self._build_overview()

        # Vendored, inlined assets (offline-safe).
        self.assets = self._load_assets()

        self.render_and_write()

    # ------------------------------------------------------------------ #
    # Sample discovery / per-sample data
    # ------------------------------------------------------------------ #
    def _metrics_csv_for(self, sample_id):
        return os.path.join(self.work_dir, f"{sample_id}.metrics.csv")

    def _discover_sample_ids(self):
        # Discover from the flat-staged metrics csvs. Sort for a stable order.
        csvs = sorted(Path(self.work_dir).glob("*.metrics.csv"))
        return [c.name[: -len(".metrics.csv")] for c in csvs]

    def _build_sample(self, sample_id):
        # Per-sample nested metrics dict, mirroring create_single_sample_report.py:
        # metrics_dict[classification][variable_name] = (human_readable, value, tooltip)
        metrics_dict = defaultdict(dict)
        with open(self._metrics_csv_for(sample_id)) as fh:
            next(fh)  # header
            for line in fh:
                var_name, var_value, human, tooltip, group = line.strip().split(",")
                metrics_dict[group][var_name] = (human, format_number_for_display(var_value), tooltip)

        cell_plot = read_fragment(os.path.join(self.work_dir, f"{sample_id}_counts_pdf_with_threshold.html"))
        qc_cascade = read_fragment(os.path.join(self.work_dir, f"{sample_id}.qc_cascade.html"))
        if self.mixed:
            barnyard = read_fragment(os.path.join(self.work_dir, f"{sample_id}_barnyard_plot.html"))
        else:
            barnyard = None

        return {
            "sample_id": sample_id,
            # plain dict so Jinja's `in` / attribute access behaves predictably
            "metrics_dict": {k: v for k, v in metrics_dict.items()},
            "cell_plot": cell_plot,
            "show_cell_plot": cell_plot is not None,
            "qc_cascade": qc_cascade,
            "show_qc_cascade": qc_cascade is not None,
            "barnyard": barnyard,
            "show_barnyard": barnyard is not None,
        }

    # ------------------------------------------------------------------ #
    # Experiment-wide data (reuses the frozen multi-sample logic)
    # ------------------------------------------------------------------ #
    def _build_multi_metrics_and_csv(self):
        """
        Parse the per-sample metrics csvs (no recomputation) to build both the
        re-emitted multisample_out.csv (byte-identical schema/values) and the
        cross-sample table. The CSV keeps the frozen format_number_to_string
        contract; the table gets thousands separators (see the return).
        """
        # metrics_dict[classification][(variable_name, human, description, classification)]
        #   = [value_per_sample, ...] in sample order
        metrics_dict = defaultdict(lambda: defaultdict(list))

        # Build from the SAME csvs, in the SAME sorted order as the table columns.
        csv_paths = [self._metrics_csv_for(sid) for sid in self.sample_ids]

        # Store the RAW per-sample values; format per-purpose below. The CSV keeps
        # the frozen format_number_to_string contract; the on-screen table gets
        # thousands separators. The two must not be conflated -- separators would
        # corrupt the comma-delimited CSV.
        for csv_path in csv_paths:
            with open(csv_path) as fh:
                next(fh)  # header
                for line in fh:
                    variable_name, value, human, description, classification = line.strip().split(",")
                    metrics_dict[classification][(variable_name, human, description, classification)].append(value)

        # Re-emit multisample_out.csv with the unchanged schema and values
        # (format_number_to_string: floats -> 2 d.p., ints unchanged; NO separators).
        with open("multisample_out.csv", "w") as csv_out:
            csv_out.write("variable_name,human_readable_name,description,classification,")
            csv_out.write(",".join(self.sample_ids) + "\n")
            for classification, subdict in metrics_dict.items():
                for (variable_name, human, description, _cls), vals in subdict.items():
                    csv_out.write(f"{','.join([variable_name, human, description, classification])},")
                    csv_out.write(",".join(format_number_to_string(v) for v in vals) + "\n")

        # Cell-metric tooltips for the cross-sample table (mixed only).
        cell_tooltip = get_cell_stat_cat_dict_obj(self.mixed)
        cell_tooltip = {k: v[0] for k, v in cell_tooltip.items()}

        # The cross-sample TABLE shows display-formatted values (thousands
        # separators); the CSV above kept the raw frozen contract.
        table_dict = {
            cls: {key: [format_number_for_display(v) for v in vals] for key, vals in sub.items()}
            for cls, sub in metrics_dict.items()
        }
        return table_dict, cell_tooltip

    def _build_summary_plot_fragment(self):
        """
        DISABLED (2026-06, by decision): the across-sample "Key metric distributions"
        violin/box summary plot is intentionally not generated. A violin/KDE drawn
        over the small number of samples these reports carry implies a continuous
        distribution that does not exist. Returning None hides the section cleanly
        via the template's `{% if summary_plot %}` guard (and the template block has
        also been removed). Reinstate by restoring this method from git history if a
        suitable small-N representation is ever chosen.
        """
        return None

    # ------------------------------------------------------------------ #
    # Overview header
    # ------------------------------------------------------------------ #
    def _build_overview(self):
        """A small set of headline experiment metrics for the top of the report."""
        # Per-sample headline metric keys (same as the old per-sample headline cards).
        if self.mixed:
            # In mixed-species mode summary_statistics.py classifies each cell metric under its own
            # name (e.g. num_cells), NOT under "Cell metrics" (which only holds counts_in/out_of_cells).
            headline_keys = [
                ("num_cells_total", "num_cells"),
                ("raw_reads_per_cell_total", "raw_reads_per_cell"),
                ("median_genes_detected_per_cell_total", "median_genes_detected_per_cell"),
            ]
        else:
            headline_keys = [
                ("num_cells", "Cell metrics"),
                ("raw_reads_per_cell", "Cell metrics"),
                ("median_genes_detected_per_cell", "Cell metrics"),
            ]
        return {"n_samples": len(self.sample_ids), "headline_keys": headline_keys}

    # ------------------------------------------------------------------ #
    # Vendored assets (inlined for offline-safety)
    # ------------------------------------------------------------------ #
    def _load_assets(self):
        plotlyjs = get_plotlyjs()  # full Plotly.js source, embedded once
        return {
            "bootstrap_css": read_text_asset(self.vendor_dir, "bootstrap.min.css"),
            "bootstrap_icons_css": read_text_asset(self.vendor_dir, "bootstrap-icons.embedded.css"),
            "lexend_css": read_text_asset(self.vendor_dir, "lexend.embedded.css"),
            "bootstrap_js": read_text_asset(self.vendor_dir, "bootstrap.bundle.min.js"),
            "plotly_js": plotlyjs,
            "logo_data_uri": read_text_asset(self.vendor_dir, "cs_logo_data_uri.txt").strip(),
        }

    # ------------------------------------------------------------------ #
    # Render
    # ------------------------------------------------------------------ #
    def render_and_write(self):
        cell_stat_cat_dict = get_cell_stat_cat_dict_obj(self.mixed)

        html = self.jinja_template.render(
            mixed=self.mixed,
            assets=self.assets,
            overview=self.overview,
            samples=self.samples,
            sample_name_list=self.sample_ids,
            multi_metrics_dict=self.multi_metrics_dict,
            cell_metric_tooltip_dict=self.cell_metric_tooltip_dict,
            cell_stat_cat_dict=cell_stat_cat_dict,
            summary_plot=self.summary_plot_fragment,
            multi_qc_cascade=self.multi_qc_cascade_fragment,
        )

        with open("consolidated_report.html", "w") as fh:
            fh.write(html)


if __name__ == "__main__":
    ConsolidatedReport()
