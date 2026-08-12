#!/usr/bin/env python

"""
Generate QC cascade plots for customer reports.

Single-sample mode: Stacked bar chart showing reads retained/lost at each QC step
Multi-sample mode: Box plots showing distribution of reads retained across samples at each step
"""

import argparse
import csv
import math
from pathlib import Path
import sys

import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots


PLOTLY_CONFIG = {"responsive": True, "displaylogo": False}
METRICS_HEADER = [
    "variable_name",
    "value",
    "human_readable_name",
    "description",
    "classification",
]
REQUIRED_READ_QC_METRICS = (
    "reads_pre_qc",
    "barcode_exact_match",
    "barcode_corrected",
    "barcode_failed",
    "barcode_valid_total",
    "sss_retained",
    "sss_too_short",
    "polyx_retained",
    "polyx_too_short",
    "quality_retained",
    "quality_failed",
    "n_base_retained",
    "n_base_failed",
    "polya_retained",
    "polya_too_short",
)


def safe_proportion(value, denominator):
    """Return a finite proportion, defining every zero/zero QC step as zero."""
    return 0.0 if denominator == 0 else value / denominator


def write_plotly_outputs(fig, fragment_filename, standalone_filename):
    """
    Write distinct consolidated-report and standalone Plotly outputs.

    ``fragment_filename`` deliberately omits Plotly.js because the consolidated
    report embeds one shared copy. ``standalone_filename`` embeds Plotly.js so
    the separately published plot remains usable without internet access.

    Keeping the two files distinct is important: publishing the fragment under
    the customer-facing filename creates a page that is blank when opened on
    its own, while feeding the standalone file into the consolidated report
    needlessly duplicates several megabytes of JavaScript for every plot.
    """
    fragment = pio.to_html(
        fig,
        full_html=False,
        include_plotlyjs=False,
        config=PLOTLY_CONFIG,
    )
    standalone_plot = pio.to_html(
        fig,
        full_html=False,
        include_plotlyjs=True,
        config=PLOTLY_CONFIG,
    )
    standalone = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src 'unsafe-inline'; style-src 'unsafe-inline'; img-src data: blob:; font-src data:; connect-src 'none'; frame-src 'none'; object-src 'none'; media-src 'none'; base-uri 'none'; form-action 'none'">
  <title>QC cascade plot</title>
  <style>
    body {{ font-family: Lexend, system-ui, -apple-system, "Segoe UI", sans-serif; margin: 0; }}
  </style>
</head>
<body>
{standalone_plot}
</body>
</html>
"""
    with open(fragment_filename, "w", encoding="utf-8") as f:
        f.write(fragment)
    with open(standalone_filename, "w", encoding="utf-8") as f:
        f.write(standalone)


class QCCascadePlotter:
    def __init__(self, mode, sample_id=None):
        self.mode = mode  # 'single' or 'multi'
        self.sample_id = sample_id

        # CSGenetics color palette
        self.cs_green = "#36BA00"
        self.cs_green_light = "#6CD34D"
        self.cs_blue = "#257FC1"
        self.cs_teal = "#34BBCF"
        self.loss_1 = "#A9A9A9"  # Dark gray
        self.loss_2 = "#808080"  # Medium gray
        self.loss_3 = "#D3D3D3"  # Light gray

    def read_metrics_csv(self, csv_path):
        """Read and validate the complete Read-QC schema used by the cascade."""
        metrics = {}
        with open(csv_path, 'r', newline='', encoding='utf-8') as f:
            reader = csv.reader(f)
            try:
                header = next(reader)
            except StopIteration as exc:
                raise ValueError(f"Metrics CSV is empty: {csv_path}") from exc
            if header != METRICS_HEADER:
                raise ValueError(
                    f"Metrics CSV {csv_path} has invalid header {header!r}; "
                    f"expected {METRICS_HEADER!r}"
                )

            for line_number, row in enumerate(reader, start=2):
                if len(row) != len(METRICS_HEADER):
                    raise ValueError(
                        f"Metrics CSV {csv_path} row {line_number} has "
                        f"{len(row)} columns; expected {len(METRICS_HEADER)}"
                    )
                var_name, value, _human_name, _description, classification = row
                if var_name not in REQUIRED_READ_QC_METRICS:
                    continue
                if classification != "Read QC":
                    raise ValueError(
                        f"Required metric {var_name!r} in {csv_path} row "
                        f"{line_number} has classification {classification!r}; "
                        "expected 'Read QC'"
                    )
                if var_name in metrics:
                    raise ValueError(
                        f"Required metric {var_name!r} occurs more than once "
                        f"in {csv_path} (duplicate at row {line_number})"
                    )
                try:
                    numeric = int(value)
                except ValueError:
                    try:
                        numeric = float(value)
                    except ValueError as exc:
                        raise ValueError(
                            f"Required metric {var_name!r} in {csv_path} row "
                            f"{line_number} is not numeric: {value!r}"
                        ) from exc
                if isinstance(numeric, float) and not math.isfinite(numeric):
                    raise ValueError(
                        f"Required metric {var_name!r} in {csv_path} row "
                        f"{line_number} is not finite: {value!r}"
                    )
                if isinstance(numeric, float) and not numeric.is_integer():
                    raise ValueError(
                        f"Required metric {var_name!r} in {csv_path} row "
                        f"{line_number} is not an integer count: {value!r}"
                    )
                if numeric < 0:
                    raise ValueError(
                        f"Required metric {var_name!r} in {csv_path} row "
                        f"{line_number} is negative: {value!r}"
                    )
                metrics[var_name] = int(numeric)

        missing = [name for name in REQUIRED_READ_QC_METRICS if name not in metrics]
        if missing:
            raise ValueError(
                f"Metrics CSV {csv_path} is missing required Read-QC metric(s): "
                + ", ".join(missing)
            )
        if metrics["reads_pre_qc"] == 0:
            nonzero = [
                name
                for name in REQUIRED_READ_QC_METRICS
                if name != "reads_pre_qc" and metrics[name] != 0
            ]
            if nonzero:
                raise ValueError(
                    f"Metrics CSV {csv_path} has zero reads_pre_qc but nonzero "
                    "downstream Read-QC metric(s): " + ", ".join(nonzero)
                )
        return metrics

    def create_single_sample_plot(self, metrics):
        """
        Create stacked bar chart for single sample showing QC cascade.
        Similar to internal pipeline qc_plotting.py but customer-friendly.
        """
        raw_reads = metrics['reads_pre_qc']

        # Build plotting data structure
        plotting_data = []

        # Step 0: Input reads (pre-QC)
        plotting_data.append({
            'step': 'Input',
            'metric': 'Input reads',
            'value': raw_reads,
            'color': self.cs_green,
            'type': 'input'
        })

        # Step 1: Barcode Validation - show stacked (exact + corrected) and failed
        plotting_data.append({
            'step': 'Barcode\nValidation',
            'metric': 'Exact match',
            'value': metrics['barcode_exact_match'],
            'color': self.cs_green,
            'type': 'retained'
        })
        plotting_data.append({
            'step': 'Barcode\nValidation',
            'metric': 'Corrected',
            'value': metrics['barcode_corrected'],
            'color': self.cs_green_light,
            'type': 'retained'
        })
        plotting_data.append({
            'step': 'Barcode\nValidation',
            'metric': 'Failed',
            'value': metrics['barcode_failed'],
            'color': self.loss_1,
            'type': 'lost'
        })

        # Step 3: SSS Trimming
        plotting_data.append({
            'step': 'SSS\nTrimming',
            'metric': 'Retained',
            'value': metrics['sss_retained'],
            'color': self.cs_green,
            'type': 'retained'
        })
        plotting_data.append({
            'step': 'SSS\nTrimming',
            'metric': 'Too short',
            'value': metrics['sss_too_short'],
            'color': self.loss_2,
            'type': 'lost'
        })

        # Step 4: PolyX Trimming
        plotting_data.append({
            'step': 'PolyX\nTrimming',
            'metric': 'Retained',
            'value': metrics['polyx_retained'],
            'color': self.cs_green,
            'type': 'retained'
        })
        plotting_data.append({
            'step': 'PolyX\nTrimming',
            'metric': 'Too short',
            'value': metrics['polyx_too_short'],
            'color': self.loss_2,
            'type': 'lost'
        })

        # Step 5: Quality Filtering
        plotting_data.append({
            'step': 'Quality\nFiltering',
            'metric': 'Retained',
            'value': metrics['quality_retained'],
            'color': self.cs_green,
            'type': 'retained'
        })
        plotting_data.append({
            'step': 'Quality\nFiltering',
            'metric': 'Failed',
            'value': metrics['quality_failed'],
            'color': self.loss_1,
            'type': 'lost'
        })

        # Step 6: N-Base Filtering
        plotting_data.append({
            'step': 'N-Base\nFiltering',
            'metric': 'Retained',
            'value': metrics['n_base_retained'],
            'color': self.cs_green,
            'type': 'retained'
        })
        plotting_data.append({
            'step': 'N-Base\nFiltering',
            'metric': 'Excessive Ns',
            'value': metrics['n_base_failed'],
            'color': self.loss_1,
            'type': 'lost'
        })

        # Step 7: PolyA Trimming
        plotting_data.append({
            'step': 'PolyA\nTrimming',
            'metric': 'Retained',
            'value': metrics['polya_retained'],
            'color': self.cs_green,
            'type': 'retained'
        })
        plotting_data.append({
            'step': 'PolyA\nTrimming',
            'metric': 'Too short',
            'value': metrics['polya_too_short'],
            'color': self.loss_2,
            'type': 'lost'
        })

        # Convert to DataFrame
        df = pd.DataFrame(plotting_data)

        # Create figure
        fig = go.Figure()

        # Add bars for each metric
        for _, row in df.iterrows():
            fig.add_trace(go.Bar(
                name=row['metric'],
                x=[row['step']],
                y=[row['value']],
                marker_color=row['color'],
                hovertemplate=f"<b>{row['step']}</b><br>" +
                              f"{row['metric']}: {row['value']:,}<br>" +
                              f"{safe_proportion(row['value'], raw_reads)*100:.2f}% of input<extra></extra>",
                showlegend=False
            ))

        # Update layout
        fig.update_layout(
            title=dict(
                text=f"QC Cascade - {self.sample_id}",
                x=0.5,
                xanchor='center',
                font=dict(size=18, family='Lexend, sans-serif')
            ),
            xaxis=dict(
                title="QC Step",
                tickfont=dict(size=11)
            ),
            yaxis=dict(
                title="Number of Reads",
                tickformat=',d'
            ),
            barmode='stack',
            font=dict(family='Lexend, sans-serif', color='black'),
            plot_bgcolor='white',
            paper_bgcolor='white',
            height=500,
            margin=dict(l=80, r=40, t=80, b=80)
        )

        # The fragment retains the historical task-local filename because the
        # consolidated-report generator discovers it by sample id. Nextflow
        # publishes only the standalone file, renaming it back to the stable
        # customer-facing ``<sample>.qc_cascade.html`` filename.
        fragment_filename = f"{self.sample_id}.qc_cascade.html"
        standalone_filename = f"{self.sample_id}.qc_cascade.standalone.html"
        write_plotly_outputs(fig, fragment_filename, standalone_filename)
        print(
            "Created single-sample QC cascade outputs: "
            f"{fragment_filename}, {standalone_filename}"
        )
        return fig

    def create_multi_sample_plot(self, csv_files):
        """
        Create box plots and line plots showing QC step retention across all samples.
        Toggle button to switch between Proportional (default) and Absolute views.
        """
        # Collect data from all samples
        all_samples_data = []

        for csv_file in csv_files:
            csv_name = Path(csv_file).name
            sample_id = csv_name.removesuffix('.metrics.csv')
            metrics = self.read_metrics_csv(csv_file)

            sample_data = {
                'sample_id': sample_id,
                'Input': metrics['reads_pre_qc'],
                'Barcode\nValidation': metrics['barcode_valid_total'],
                'SSS\nTrimming': metrics['sss_retained'],
                'PolyX\nTrimming': metrics['polyx_retained'],
                'Quality\nFiltering': metrics['quality_retained'],
                'N-Base\nFiltering': metrics['n_base_retained'],
                'PolyA\nTrimming': metrics['polya_retained']
            }
            all_samples_data.append(sample_data)

        df = pd.DataFrame(all_samples_data)

        # Define QC steps in order (removed redundant Final Output)
        steps = ['Input', 'Barcode\nValidation', 'SSS\nTrimming', 'PolyX\nTrimming',
                 'Quality\nFiltering', 'N-Base\nFiltering', 'PolyA\nTrimming']

        # Calculate proportional data (normalized to Input = 1.0)
        df_prop = df.copy()
        for step in steps:
            df_prop[step] = [
                safe_proportion(value, input_reads)
                for value, input_reads in zip(df[step], df['Input'])
            ]

        # Create color palette for samples
        import plotly.express as px
        colors = px.colors.qualitative.Set2
        if len(df) > len(colors):
            # Repeat colors if more samples than colors
            colors = colors * (len(df) // len(colors) + 1)

        # Create figure with 2 subplots (box plot on top, line plot below)
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=("Box Plot View", "Line Plot View"),
            vertical_spacing=0.20,  # Increased spacing to avoid overlap
            row_heights=[0.5, 0.5]
        )

        # === PROPORTIONAL TRACES (visible by default) ===

        # Box plots (row 1) - Proportional
        for i, step in enumerate(steps):
            step_data_prop = df_prop[step].dropna()
            fig.add_trace(
                go.Box(
                    y=step_data_prop,
                    x=[i] * len(step_data_prop),
                    name=step,
                    marker_color=self.cs_green,
                    boxmean='sd',
                    showlegend=False,
                    boxpoints=False,
                    visible=True,
                    hovertemplate='<b>' + step + '</b><br>' +
                                  'Median: %{median:.2%}<br>' +
                                  'Mean: %{mean:.2%}<br>' +
                                  'Q1: %{q1:.2%}<br>' +
                                  'Q3: %{q3:.2%}<br>' +
                                  'Min: %{lowerfence:.2%}<br>' +
                                  'Max: %{upperfence:.2%}<br>' +
                                  f'n = {len(step_data_prop)}' +
                                  '<extra></extra>'
                ),
                row=1, col=1
            )

        # Scatter points for box plots (row 1) - Proportional
        for i, step in enumerate(steps):
            step_data_prop = df_prop[step].dropna()
            fig.add_trace(
                go.Scatter(
                    y=step_data_prop,
                    x=[i + 0.35] * len(step_data_prop),
                    mode='markers',
                    marker=dict(
                        color='rgba(54, 186, 0, 0.6)',
                        size=6,
                        line=dict(width=1, color='white')
                    ),
                    customdata=df['sample_id'][step_data_prop.index],
                    hovertemplate='<b>%{customdata}</b><br>' +
                                  f'{step}<br>' +
                                  'Proportion: %{y:.2%}<extra></extra>',
                    showlegend=False,
                    visible=True
                ),
                row=1, col=1
            )

        # Line plots (row 2) - Proportional - one line per sample
        for sample_idx, sample_row in df_prop.iterrows():
            sample_id = sample_row['sample_id']
            y_values = [sample_row[step] for step in steps]

            fig.add_trace(
                go.Scatter(
                    x=list(range(len(steps))),
                    y=y_values,
                    mode='lines+markers',
                    name=sample_id,
                    line=dict(color=colors[sample_idx % len(colors)], width=2),
                    marker=dict(size=8),
                    customdata=steps,
                    hovertemplate=f'<b>{sample_id}</b><br>' +
                                  '%{customdata}<br>' +
                                  'Proportion: %{y:.2%}<extra></extra>',
                    showlegend=False,
                    visible=True
                ),
                row=2, col=1
            )

        # === ABSOLUTE TRACES (hidden by default) ===

        # Box plots (row 1) - Absolute
        for i, step in enumerate(steps):
            step_data_abs = df[step].dropna()
            fig.add_trace(
                go.Box(
                    y=step_data_abs,
                    x=[i] * len(step_data_abs),
                    name=step,
                    marker_color=self.cs_green,
                    boxmean='sd',
                    showlegend=False,
                    boxpoints=False,
                    visible=False,
                    hovertemplate='<b>' + step + '</b><br>' +
                                  'Median: %{median:,}<br>' +
                                  'Mean: %{mean:,}<br>' +
                                  'Q1: %{q1:,}<br>' +
                                  'Q3: %{q3:,}<br>' +
                                  'Min: %{lowerfence:,}<br>' +
                                  'Max: %{upperfence:,}<br>' +
                                  f'n = {len(step_data_abs)}' +
                                  '<extra></extra>'
                ),
                row=1, col=1
            )

        # Scatter points for box plots (row 1) - Absolute
        for i, step in enumerate(steps):
            step_data_abs = df[step].dropna()
            fig.add_trace(
                go.Scatter(
                    y=step_data_abs,
                    x=[i + 0.35] * len(step_data_abs),
                    mode='markers',
                    marker=dict(
                        color='rgba(54, 186, 0, 0.6)',
                        size=6,
                        line=dict(width=1, color='white')
                    ),
                    customdata=df['sample_id'][step_data_abs.index],
                    hovertemplate='<b>%{customdata}</b><br>' +
                                  f'{step}<br>' +
                                  'Reads: %{y:,}<extra></extra>',
                    showlegend=False,
                    visible=False
                ),
                row=1, col=1
            )

        # Line plots (row 2) - Absolute - one line per sample
        for sample_idx, sample_row in df.iterrows():
            sample_id = sample_row['sample_id']
            y_values = [sample_row[step] for step in steps]

            fig.add_trace(
                go.Scatter(
                    x=list(range(len(steps))),
                    y=y_values,
                    mode='lines+markers',
                    name=sample_id,
                    line=dict(color=colors[sample_idx % len(colors)], width=2),
                    marker=dict(size=8),
                    customdata=steps,
                    hovertemplate=f'<b>{sample_id}</b><br>' +
                                  '%{customdata}<br>' +
                                  'Reads: %{y:,}<extra></extra>',
                    showlegend=False,
                    visible=False
                ),
                row=2, col=1
            )

        # Create toggle buttons
        # Trace counts: box (len(steps)) + scatter (len(steps)) + lines (len(df)) per view
        box_scatter_count = len(steps) * 2
        line_count = len(df)
        traces_per_view = box_scatter_count + line_count

        # Proportional: first traces_per_view visible, rest hidden
        proportional_visible = [True] * traces_per_view + [False] * traces_per_view

        # Absolute: first traces_per_view hidden, rest visible
        absolute_visible = [False] * traces_per_view + [True] * traces_per_view

        # Update layout with toggle buttons
        fig.update_layout(
            updatemenus=[
                dict(
                    type="buttons",
                    direction="left",
                    buttons=list([
                        dict(
                            args=[{"visible": proportional_visible},
                                  {"yaxis.title": "Proportion of Input Reads",
                                   "yaxis.tickformat": ".0%",
                                   "yaxis2.title": "Proportion of Input Reads",
                                   "yaxis2.tickformat": ".0%"}],
                            label="Proportional",
                            method="update"
                        ),
                        dict(
                            args=[{"visible": absolute_visible},
                                  {"yaxis.title": "Number of Reads",
                                   "yaxis.tickformat": ",d",
                                   "yaxis2.title": "Number of Reads",
                                   "yaxis2.tickformat": ",d"}],
                            label="Absolute",
                            method="update"
                        )
                    ]),
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x=0.0,
                    xanchor="left",
                    y=1.08,  # Moved higher to be above plot area
                    yanchor="top",
                    bgcolor="white",
                    bordercolor="#36BA00",
                    borderwidth=2,
                    font=dict(size=12)
                ),
            ]
        )

        # Update layout
        fig.update_layout(
            title=dict(
                text="QC Cascade Across All Samples",
                x=0.5,
                xanchor='center',
                font=dict(size=18, family='Lexend, sans-serif')
            ),
            font=dict(family='Lexend, sans-serif', color='black'),
            plot_bgcolor='white',
            paper_bgcolor='white',
            height=1000,  # Taller to accommodate both plots
            showlegend=False,
            margin=dict(l=80, r=40, t=120, b=120)
        )

        # Update x-axes for both subplots
        # Top subplot: show labels (needed for hover text)
        fig.update_xaxes(
            tickmode='array',
            tickvals=list(range(len(steps))),
            ticktext=steps,  # Show step names
            tickangle=-90,
            tickfont=dict(size=10),
            range=[-0.5, len(steps) - 0.5],
            row=1, col=1
        )
        # Bottom subplot: show labels
        fig.update_xaxes(
            title="QC Step",
            tickangle=-90,
            tickfont=dict(size=10),
            tickmode='array',
            tickvals=list(range(len(steps))),
            ticktext=steps,
            range=[-0.5, len(steps) - 0.5],
            row=2, col=1
        )

        # Update y-axes (default to proportional)
        fig.update_yaxes(
            title="Proportion of Input Reads",
            tickformat='.0%',
            row=1, col=1
        )
        fig.update_yaxes(
            title="Proportion of Input Reads",
            tickformat='.0%',
            row=2, col=1
        )

        fragment_filename = "multisample_qc_cascade.fragment.html"
        standalone_filename = "multisample_qc_cascade.html"
        write_plotly_outputs(fig, fragment_filename, standalone_filename)
        print(
            "Created multi-sample QC cascade outputs: "
            f"{fragment_filename}, {standalone_filename}"
        )
        return fig


def main():
    parser = argparse.ArgumentParser(description='Generate QC cascade plots for customer reports')
    parser.add_argument('--mode', required=True, choices=['single', 'multi'],
                       help='Plot mode: single or multi sample')
    parser.add_argument('--sample-id', required=False,
                       help='Sample ID (required for single mode)')
    parser.add_argument('--metrics-csv', required=False,
                       help='Metrics CSV file (required for single mode)')
    parser.add_argument('--csv-files', required=False, nargs='+',
                       help='List of metrics CSV files (required for multi mode)')

    args = parser.parse_args()

    # Validate arguments
    if args.mode == 'single':
        if not args.sample_id or not args.metrics_csv:
            print("Error: --sample-id and --metrics-csv required for single mode")
            sys.exit(1)
    elif args.mode == 'multi':
        if not args.csv_files:
            print("Error: --csv-files required for multi mode")
            sys.exit(1)

    # Create plotter
    plotter = QCCascadePlotter(mode=args.mode, sample_id=args.sample_id)

    # Generate plot
    if args.mode == 'single':
        metrics = plotter.read_metrics_csv(args.metrics_csv)
        plotter.create_single_sample_plot(metrics)
    else:
        plotter.create_multi_sample_plot(args.csv_files)


if __name__ == "__main__":
    main()
