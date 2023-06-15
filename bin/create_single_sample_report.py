#!/usr/bin/python3
import sys
from jinja2 import Template
import base64

class SingleSampleHTMLReport:
    """
        This script will take CSV outputs for:
        - Cell statistics 
        - Mapping Statistics
        - Sequencing Statistics
        output from the web_sumamry.R script and insert them into the HTML for each sample.
        We also inserts a premade png for cell caller
    """
    def __init__(self):
        
        self.sample_id = sys.argv[1]
        self.plot_path = sys.argv[2]

        self.metric_name_to_jinja_var_name_dict = self.get_metric_name_to_jinja_var_name_dict()

        with open(sys.argv[3]) as metrics_handle:
            self.metrics_dict = {line.split("\t")[0]: line.split("\t")[1].rstrip() for line in metrics_handle}

        self.metrics_dict = {self.metric_name_to_jinja_var_name_dict[k]: v for k, v in self.metrics_dict.items()}

        self.metrics_dict["sample_id"] = self.sample_id

        with open(sys.argv[4]) as html_template:
            self.jinja_template = Template(html_template.read())

        self.encoded_png_str = self.generate_encoded_png_str()

        self.metrics_dict["encoded_png_str"] = self.encoded_png_str

        self.report = self.jinja_template.render(self.metrics_dict)

        self.write_report()

    def generate_encoded_png_str(self):
        with open(f"{self.plot_path}", "rb") as image_file:
            return base64.b64encode(image_file.read()).decode("utf-8")
        
    def insert_sample_id(self):
        self.insert_id("sample_id", self.sample_id, "h2")

    def insert_mapping_data(self):
        self.loop_dat(self.mapping_data, "td")

    def insert_cell_data(self):
        self.loop_dat(self.cell_data, "td")

    def insert_sequencing_data(self):
        self.loop_dat(self.seq_data, "td")

    def insert_cell_banner_data(self):
        self.loop_dat(self.cell_banner_data, "div")

    def write_report(self):
        with open(f"{self.sample_id}_report.html", "w") as report_out:
            report_out.write(self.report)

    def get_metric_name_to_jinja_var_name_dict(self):
        """
        Return a dictionary that is used to convert the metric
        names in the .metrics.csv input file to the 
        Jinja variable names in the report template.
        """
        return {
            "Reads Mapped to Genome": "mapped_reads",
            "Reads Mapped Confidently to Genome": "confidently_mapped_reads",
            "Reads Mapped Confidently to Intergenic Regions": "confidently_mapped_reads_intergenic",
            "Reads Mapped Confidently to Intronic Regions": "confidently_mapped_reads_intronic",
            "Reads Mapped Confidently to Exonic Regions": "confidently_mapped_reads_exonic",
            "Reads Mapped Antisense to Gene": "mapped_reads_antisense",
            "Estimated Number of Cells": "estimated_number_of_cells",
            "Estimated_Number_of_Cells": "estimated_number_of_cells",
            "Mean Reads per Cell": "mean_reads_per_cell",
            "Mean Genes per Cell": "mean_genes_per_cell",
            "Median Genes per Cell": "median_genes_per_cell",
            "Total Genes Detected": "total_genes_detected",
            "Mean Counts per Cell": "mean_counts_per_cell",
            "Median Counts per Cell": "median_counts_per_cell",
            "Mitochondrial Ratio": "mito_ratio",
            "Mean_Reads_per_Cell": "mean_reads_per_cell",
            "Total_genes_detected": "total_genes_detected",
            "Total Number of Reads": "total_num_reads",
            "Valid Barcodes / CSGX Cell Barcodes": "valid_barcodes",
            "Sequencing Saturation": "sequencing_saturation",
            "Q30 Bases in Barcode": "q30_barcode",
            "Q30 Bases in RNA Read": "q30_rna"
            }

if __name__ == "__main__":
    SingleSampleHTMLReport()