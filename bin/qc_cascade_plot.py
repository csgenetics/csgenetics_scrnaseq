#!/usr/bin/env python

"""
Generate QC cascade plots for customer reports.

Single-sample mode: Stacked bar chart showing reads retained/lost at each QC step
Multi-sample mode: Box plots showing distribution of reads retained across samples at each step
"""

import argparse
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys

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
        """Read metrics from CSV and return as dictionary"""
        metrics = {}
        with open(csv_path, 'r') as f:
            header = next(f)  # Skip header
            for line in f:
                parts = line.strip().split(',', 4)
                if len(parts) < 5:
                    continue
                var_name, value, human_name, description, classification = parts
                if classification == "Read QC":
                    try:
                        metrics[var_name] = int(value)
                    except ValueError:
                        try:
                            metrics[var_name] = float(value)
                        except ValueError:
                            pass  # Skip non-numeric values
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
                              f"{row['value']/raw_reads*100:.2f}% of input<extra></extra>",
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

        # Save as HTML
        fig.write_html(f"{self.sample_id}.qc_cascade.html", include_plotlyjs='cdn')
        print(f"Created single-sample QC cascade plot: {self.sample_id}.qc_cascade.html")

    def create_multi_sample_plot(self, csv_files):
        """
        Create box plots showing QC step retention across all samples.
        One box plot per QC step showing distribution of reads retained.
        """
        # Collect data from all samples
        all_samples_data = []

        for csv_file in csv_files:
            sample_id = csv_file.replace('.metrics.csv', '')
            metrics = self.read_metrics_csv(csv_file)

            sample_data = {
                'sample_id': sample_id,
                'Input': metrics.get('reads_pre_qc', 0),
                'Barcode\nValidation': metrics.get('barcode_valid_total', 0),
                'SSS\nTrimming': metrics.get('sss_retained', 0),
                'PolyX\nTrimming': metrics.get('polyx_retained', 0),
                'Quality\nFiltering': metrics.get('quality_retained', 0),
                'N-Base\nFiltering': metrics.get('n_base_retained', 0),
                'PolyA\nTrimming': metrics.get('polya_retained', 0),
                'Final\nOutput': metrics.get('reads_post_qc', 0)
            }
            all_samples_data.append(sample_data)

        df = pd.DataFrame(all_samples_data)

        # Create figure with subplots
        steps = ['Input', 'Barcode\nValidation', 'SSS\nTrimming', 'PolyX\nTrimming',
                 'Quality\nFiltering', 'N-Base\nFiltering', 'PolyA\nTrimming', 'Final\nOutput']

        n_steps = len(steps)

        # Calculate global y-axis range: max from first step, min from last step
        y_max = df[steps[0]].max()
        y_min = df[steps[-1]].min()
        y_range = [y_min * 0.95, y_max * 1.05]  # Add 5% padding

        fig = make_subplots(
            rows=1, cols=n_steps,
            subplot_titles=[f"{step}<br>(n={df[step].dropna().shape[0]})" for step in steps],
            horizontal_spacing=0.02
        )

        # Add box plot and scatter points for each step
        for i, step in enumerate(steps, 1):
            # Add box plot
            fig.add_trace(
                go.Box(
                    y=df[step],
                    name=step,
                    marker_color=self.cs_green,
                    boxmean='sd',
                    showlegend=False,
                    hoverinfo='skip'  # Disable box hover, we'll use scatter points instead
                ),
                row=1, col=i
            )

            # Add individual data points as scatter
            fig.add_trace(
                go.Scatter(
                    y=df[step],
                    x=[0] * len(df[step]),  # All points at x=0 (center of box)
                    mode='markers',
                    marker=dict(
                        color='rgba(54, 186, 0, 0.6)',  # Semi-transparent green
                        size=6,
                        line=dict(width=1, color='white')
                    ),
                    customdata=df['sample_id'],
                    hovertemplate='<b>%{customdata}</b><br>' +
                                  'Reads: %{y:,}<extra></extra>',
                    showlegend=False
                ),
                row=1, col=i
            )

            # Update y-axis for each subplot with fixed range
            fig.update_yaxes(
                title_text="Reads" if i == 1 else "",
                tickformat=',d',
                range=y_range,
                row=1, col=i
            )

            # Hide x-axis labels
            fig.update_xaxes(
                showticklabels=False,
                row=1, col=i
            )

        # Update overall layout
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
            height=600,
            showlegend=False,
            margin=dict(l=80, r=40, t=100, b=40)
        )

        # Save as HTML
        fig.write_html("multisample_qc_cascade.html", include_plotlyjs='cdn')
        print(f"Created multi-sample QC cascade plot: multisample_qc_cascade.html")


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
