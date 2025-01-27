#!/usr/bin/env python

import sys
from jinja2 import Template
import base64
from collections import defaultdict
import os

def get_cell_stat_cat_dict_obj(mixed):
    if mixed:
        return {
                                "num_cells": ("Estimated number of cells; Number of barcodes passing the counts thresholds.", "num_cells_accord_header", "num_cells_collapse"),
                                "raw_reads_per_cell": ("Number of reads pre-QC / Number of cells", "raw_reads_per_cell_accord_header", "raw_reads_per_cell_collapse"), 
                                "mean_total_counts_per_cell": ("Mean sum of counts per cell.", "mean_total_counts_accord_header", "mean_total_counts_collapse"), 
                                "median_total_reads_per_cell": ("Median sum of counts per cell.", "median_total_reads_accord_header", "median_total_reads_collapse"), 
                                "mean_genes_detected_per_cell": ("Mean number of genes detected for each cell (including nuclear and mitochondrial genes).", "mean_genes_detected_accord_header", "mean_genes_detected_collapse"), 
                                "median_genes_detected_per_cell": ("Median number of genes detected for the cells (including nuclear and mitochondrial genes).", "median_genes_detected_accord_header", "median_genes_detected_collapse"), 
                                "mean_nuclear_genes_detected_per_cell": ("Mean nuclear (non-mitochondrial) genes detected per cell", "mean_nuclear_genes_detected_accord_header", "mean_nuclear_genes_detected_collapse"), 
                                "median_nuclear_genes_detected_per_cell": ("Median nuclear (non-mitochondrial) genes detected per cell", "median_nuclear_genes_detected_accord_header", "median_nuclear_genes_detected_collapse"), 
                                "mean_mito_genes_detected_per_cell": ("Mean number of mitochondrial genes detected for each cell.", "mean_mito_genes_detected_accord_header", "mean_mito_genes_detected_collapse"), 
                                "median_mito_genes_detected_per_cell": ("Median number of mitochondrial genes detected for each cell.", "median_mito_genes_detected_accord_header", "median_mito_genes_detected_collapse"), 
                                "percentage_counts_from_mito": ("(Mitochonrial counts / all counts) * 100", "percentage_counts_mito_accord_header", "percentage_counts_mito_collapse"), 
                                "num_unique_genes_detected_across_sample": ("Number of unique genes detected across the sample (each gene can be counted only once even if found in multiple cells).", "num_unique_genes_detected_accord_header", "num_unique_genes_detected_collapse"), 
                                "total_genes_detected_across_sample": ("Total number of genes detected across the sample (each gene can be counted more than once if detected in more than one cell).", "total_genes_detected_accord_header", "total_genes_detected_collapse"),
                                "num_raw_cells": ("Number of called cells (i.e. barcodes meeting either of the species count thresholds).", "num_raw_cells_accord_header", "num_raw_cells_collapse"),      
                                "num_multiplet_cells": ("Total number of multiplet cells (barcodes where Hsap counts exceed the Hsap threshold AND Mmus counts exceed the Mmus threshold).", "num_multiplet_cells_accord_header", "num_multiplet_cells_collapse")                                                          
                                }
    else:
        return {}

class SingleSampleHTMLReport:
    """
        This script takes a .csv input from the summary_statistics.py script
        and inserts the metrics into the HTML for each sample.
        We also inserts a premade png for cell caller.
    """
    def __init__(self):
        
        self.sample_id = sys.argv[1]
        self.pdf_plot_path = sys.argv[2]
        self.barnyard_plot_path = sys.argv[3]

        # Create a nested dict grouping by table
        self.metrics_dict = defaultdict(dict)
        with open(sys.argv[4]) as metrics_handle:
            header_line = next(metrics_handle)
            for line in metrics_handle:
                var_name, var_value, var_human_readable_name, var_tooltip, var_group = line.strip().split(",")
                self.metrics_dict[var_group][var_name] = (var_human_readable_name, self.format_number_to_string(var_value), var_tooltip)

        with open(sys.argv[5]) as html_template:
            self.jinja_template = Template(html_template.read())

        if sys.argv[6].upper() == "TRUE":
            self.mixed = True
        else:
            self.mixed = False

        # If size of cell caller plot file is 0, then set self.show_cell_caller_plot to False.
        # This will cause the Cell Caller plot to be hidden.
        if os.path.getsize(self.pdf_plot_path) == 0:
            self.show_cell_caller_plot=False
            self.pdf_plot = None
            print("Cell Caller plot size is 0. Hiding Cell Caller plot card.")
        else:
            self.show_cell_caller_plot=True
            self.pdf_plot = self.read_html_plot(self.pdf_plot_path)

        # If mixed species and the size of barnyard plot file is 0, 
        # then set self.show_barnyard_plot to False.
        # If single species, then set self.show_barnyard_plot to False always.
        if self.mixed:
            if os.path.getsize(self.barnyard_plot_path) == 0:
                self.show_barnyard_plot=False
                self.barnyard_plot = None
                print("Barnyard plot size is 0. Hiding Barnyard plot card.")
            else:
                self.show_barnyard_plot=True
                self.barnyard_plot = self.read_html_plot(self.barnyard_plot_path)
        else:
            self.show_barnyard_plot=False
            self.barnyard_plot = None

        self.render_and_write_report()

    def read_html_plot(self, path):
        with open(f"{path}", "r", encoding="utf-8") as html_file:
            return html_file.read()

    def render_and_write_report(self):
        # For mixed species cell metrics we need to create 'cell_stat_cat_dict'.
        # This should take exactly the same for as the alignment_cat_dict below.
        # The key should be the various cell base metric names i.e "num_cells", "raw_reads_per_cell"
        # (all of the primary keys that we used for the cell metrics in the summary_statistics.py).
        
        # The values should be a tuple where the first element is the tool tip,
        # is the name that will be used for the accordion for that stat e.g. "num_cells_accord_header", "raw_reads_per_cell_accord_header" etc.
        # and the third the name of the collapse e.g. "num_cells_collapse", "raw_reads_per_cell_accord_header" etc.

        # I have implemented the jinja logic in both of the templates to accept the 'cell_stat_cat_dict" based on what we had for
        # the alignments (which we deactivated by removing the associated HTML).

        # NOTE, unlike the alignment stats, where there was only a tool tip associated with the alignement type header and not with each of the individual metrics,
        # we want to implement a tool tip for each of the individual metrics. The tooltips will already
        # have been created by you in summary_statistics. This would require a small modification to the HTML for the single and multisample template.

        cell_stat_cat_dict_obj = get_cell_stat_cat_dict_obj(self.mixed)

        final_report = self.jinja_template.render(metrics_dict=self.metrics_dict,
                                            show_cell_caller_plot=self.show_cell_caller_plot,
                                            pdf_plot=self.pdf_plot,
                                            show_barnyard_plot=self.show_barnyard_plot,
                                            barnyard_plot=self.barnyard_plot,
                                            sample_id=self.sample_id,
                                            mixed = self.mixed,
                                            # These dicts are required to supply the tooltips, the accordion header ID and the collapse ID for each of the categories of alignment statistics
                                            # The keys match the alignment categories of the metrics dict, and the values (tuples) give the tooltips and IDs.
                                            # We rely on the order of the dictionary to populate the accordions in the output html.
                                            # NOTE the alignment section has temporarily been disabled by removing the corresponding html.
                                            # I have left the code in place so that it can be reestablished.
                                            alignment_cat_dict = {
                                                "Post read QC alignment": ("Mapping of the post QC reads i.e. after trimming (polyX end and internal polyA) and barcode verification.", "qc_accord_header", "qc_collapse"),
                                                "Annotated reads alignment": ("High confidence reads annotated with a gene ID (XT bam tag).", "ann_accord_header", "ann_collapse")
                                            },  
                                            cell_stat_cat_dict = cell_stat_cat_dict_obj)

        # Write out the rendered template.
        with open(f"{self.sample_id}_report.html", "w") as f_output:
            f_output.write(final_report)
    
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
    SingleSampleHTMLReport()