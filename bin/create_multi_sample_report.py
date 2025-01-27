#!/usr/bin/env python

"""
Creates a multi sample report that summarizes the single summary report metrics
for each sample.

It reads in a single csv per sample that contains the metrics.

The metrics are then converted into a dictionary where k=metric name
and v=a list of metric values, where the metric values are reported in
order of a sample_names list that is also used to render the template.
"""

import glob
from jinja2 import Template
from pathlib import Path
from collections import defaultdict
import sys
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots 
from create_single_sample_report import get_cell_stat_cat_dict_obj
from cell_caller import output_plot_to_html

class MultipleSampleSummaries:
    def __init__(self):
        with open(sys.argv[1], "r") as html:
            self.template = Template(html.read())

        if sys.argv[2].upper() == "TRUE":
            self.mixed = True
        else:
            self.mixed = False

        self.sample_name_list, self.metrics_dict = self.make_metrics_dict_and_sample_list()
        self.render_and_write_report()
        self.write_out_csv_and_make_plotting_df()
        self.generate_summary_plots()

    def write_out_csv_and_make_plotting_df(self):
        """
        Write out all metrics as a csv and creates a DataFrame of select key metrics for plotting
        """
        # List of metrics to plot in the summary plots, and associated human-readable labels
        if self.mixed:
            self.metrics_to_plot = {"reads_pre_qc": "Number of reads pre-QC",
                                    "num_cells_total": "Number of single-cells (Hsap and Mmus)",
                                    "raw_reads_per_cell_total": "Raw reads per single-cell",
                                    "median_genes_detected_per_cell_total": "Median genes detected per single-cell"}
        else:
            self.metrics_to_plot = {"reads_pre_qc": "Number of reads pre-QC",
                                    "num_cells": "Number of cells",
                                    "raw_reads_per_cell": "Raw reads per cell",
                                    "median_genes_detected_per_cell": "Median genes detected per cell"
                                    }

        # Initialize a list to store the data for the plotting df
        data_for_plotting = {}

        # open a stream to write to
        with open("multisample_out.csv", "w") as csv_out:
            # Put in an empty field to allow for the variable names
            csv_out.write("variable_name,human_readable_name,description,classification,")
            csv_out.write(",".join(self.sample_name_list) + "\n")
            for subcategory, subdict in self.metrics_dict.items():
                for (variable_name, human_readable_name, description, classification), metric_val_list in subdict.items():
                    csv_out.write(f"{','.join([variable_name, human_readable_name, description, classification])},")
                    csv_out.write(",".join(metric_val_list) + "\n")

                    if variable_name in self.metrics_to_plot.keys():
                        data_for_plotting[variable_name] = metric_val_list

        # Create a DataFrame for plotting
        self.plotting_df = pd.DataFrame.from_dict(data_for_plotting)
        self.plotting_df.set_index(pd.Index(self.sample_name_list), inplace=True)
        for metric in self.plotting_df.columns:
            self.plotting_df[metric] = self.plotting_df[metric].astype(float)

    def generate_summary_plots(self):
        """
        Plots violin plots across all samples for a set of key metrics
        """
        # Calculate y-axis ranges for each metric
        y_ranges = self.calculate_y_ranges()

        # Generate subplot titles which include number of datapoints available for each metric
        subplot_titles = self.generate_subplot_titles()

        csgx_colors = [
                            "rgb(54,186,0)",  # CSG Green
                            "rgb(37,127,193)",  # CSG Blue
                            "rgb(52,187,207)"  # CSG Teal
                        ]

        n_metrics = len(self.metrics_to_plot)

        fig = make_subplots(
            rows=1, cols=n_metrics,
            subplot_titles=subplot_titles,
            horizontal_spacing=0.05  # Adjust this value to reduce the white space
        )

        for i, (metric, label) in enumerate(self.metrics_to_plot.items()):
            fig.add_trace(go.Violin(
                y=self.plotting_df[metric], 
                name=label, 
                box_visible=True, 
                points="all",
                line_color=csgx_colors[i % 3],  
                fillcolor=csgx_colors[i % 3],  
                opacity=0.6,
                box=dict(visible=True, line_color='black'),  # Ensure box plot is visible
                customdata=self.plotting_df.index,  # Add sample names as custom data
                hovertemplate=f'<b>%{{customdata}}</b><br>{label}: %{{text}}<extra></extra>',  # Customize hover text
                text=[self.abbreviate_number(y) for y in self.plotting_df[metric]]  # Use abbreviated values for hover text
            ),
            row=1, col=i+1
            )
            fig.update_yaxes(title_text=label, 
                             range=y_ranges[metric],
                             row=1, col=i+1)
            fig.update_xaxes(showticklabels=False, 
                             row=1, col=i+1)

        fig.update_layout(
            showlegend=False,
            font=dict(family = 'Lexend, sans-serif', color="black"),
            title_text="Multisample summary plots",
            margin = dict(l=45,r=45,t=150,b=150),
            autosize=True
        )

        output_plot_to_html({"multisample_summary_plots":fig}, "multisample_summary_plots.html")

    def generate_subplot_titles(self):
        # Create subplot titles with number of datapoints
        subplot_titles = []
        for metric in self.metrics_to_plot.keys():
            num_datapoints = self.plotting_df[metric].dropna().shape[0]
            label = self.metrics_to_plot[metric]
            title = f"{label} (n={num_datapoints})"
            subplot_titles.append(title)
        return subplot_titles

    def calculate_y_ranges(self):
        # Calculate y-axis ranges for each key metric
        y_ranges = {}
        for metric in self.plotting_df.columns:
            all_values = self.plotting_df[metric].dropna()
            y_min = np.max([all_values.min() - 2 * all_values.std(), 0])
            y_max = all_values.max() + 2 * all_values.std()
            y_ranges[metric] = [y_min, y_max]
        return y_ranges
    
    def abbreviate_number(self, value):
        # Custom function to abbreviate numbers for purposes of plot labels
        if value >= 1_000_000:
            return f"{value / 1_000_000:.1f} M"
        elif value >= 1_000:
            return f"{value / 1_000:.1f} k"
        else:
            return f"{value:.1f}"

    def render_and_write_report(self):
        """
        Render the template.
        We pass the alignmnet_tooltip_dict so that we can provide the tooltips
        for each of the alignment subcategories rather than for each
        individual alignment metric.
        NOTE the alignment card has temporarily been disabled by removing the
        associated html. I am leaving the associated code in place so that it
        can be reinstated.
        """
        
        # Populate cell_metric_tooltip_dict here if a mixed sample.
        # The primary keys will need to match the base metric primary keys that have been used in summary_statistics.py
        # e.g. "num_cells", "raw_reads_per_cell" etc. etc. Similar to the single sample summary report you may want
        # to implement tool tips for each of the individual metrics as well. This will require a small modificaiton to the
        # HTML.
        cell_metric_tooltip_dict_obj = get_cell_stat_cat_dict_obj(self.mixed)
        # slim down to just the tool tips (i.e. get rid of the accordian IDs)
        cell_metric_tooltip_dict_obj = {k: v[0] for k, v in cell_metric_tooltip_dict_obj.items()}

        print(cell_metric_tooltip_dict_obj)
        print(self.metrics_dict)
        final_report = self.template.render(
            sample_name_list=self.sample_name_list, metrics_dict=self.metrics_dict,
            mixed=self.mixed,
            alignment_tooltip_dict = {
                "Post read QC alignment": "Mapping of the post QC reads i.e. after trimming (polyX end and internal polyA) and barcode verification.",
                "Annotated reads alignment": "High confidence reads annotated with a gene ID (XT bam tag)."
                }, cell_metric_tooltip_dict = cell_metric_tooltip_dict_obj)

        # Write out the rendered template.
        with open(f"multisample_report.html", "w") as f_output:
            f_output.write(final_report)

    def make_metrics_dict_and_sample_list(self):
        """
        Create the two data structures that will be used to render the template.
        self.sample_list is simply a list that contains the names of the samples
        sorted in some sort of order.
        self.metrics_dict is a defaultdict(defaultdict(list)). key is a metric classification,
        value is a dict, where the key is a tuple of (variable_name, human_readable_name, description, classification)
        and the value is a list
        containing the values for that metric for each of the samples in the
        same sample order as self.sample_list.
        """

        self.metrics_dict = defaultdict(lambda: defaultdict(list))
        self.sample_name_list = []

        # Get all the csv paths and sort them
        # so that the samples will be reported in
        # a logical order
        csvs = [_ for _ in Path(".").glob('./*.metrics.csv')]
        csvs.sort()

        # For each sample csv in the sorted paths
        for csv_path in csvs:
            self.sample_name_list.append(csv_path.stem.replace(".metrics", ""))
            with open(csv_path, "r") as csv_f:
                header_line = next(csv_f).rstrip()
                for line in csv_f:
                    # read in the metrics and add them to the defaultdict
                    # Each line is in format: variable_name,value,human_readable_name,description,classification
                    # We will carry all of this information into the metrics dict.
                    # We only use the var_human_readable_name and var_toolstip, but we want to ouput all
                    # of the info in the csv we output.
                    # The classifications are:
                    #   Read QC
                    #   Cell metrics
                    #   Deduplication
                    #   Post read QC alignment
                    #   Annotated reads alignment
                    variable_name, value, human_readable_name, description, classification = line.strip().split(",")
                    self.metrics_dict[classification][(variable_name, human_readable_name, description, classification)].append(self.format_number_to_string(value))

        return self.sample_name_list, self.metrics_dict

    @staticmethod
    def format_number_to_string(number_str):
        """
        Ensure that the proper formatting of numbers
        for the output html.

        Floats should be formatted to 2 d.p.

        Integers should not have any decimal places.
        """
        if "." in number_str:
            # If float then return with 2 d.p.
            return f"{float(number_str):.2f}"
        else:
            # Else return original int
            return number_str

if __name__ == "__main__":
    MultipleSampleSummaries()