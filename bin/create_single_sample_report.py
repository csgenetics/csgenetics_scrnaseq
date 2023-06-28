#!/usr/bin/python3
import sys
from jinja2 import Template
import base64
from collections import defaultdict

class SingleSampleHTMLReport:
    """
        This script will take CSV outputs for:
        - Sample.metrics.csv file
        output from the summary_statistics.py script and insert them into the HTML for each sample.
        We also inserts a premade png for cell caller.
    """
    def __init__(self):
        
        self.sample_id = sys.argv[1]
        self.plot_path = sys.argv[2]

        # Create a nested dict grouping by table
        self.metrics_dict = defaultdict(dict)
        with open(sys.argv[3]) as metrics_handle:
            header_line = next(metrics_handle)
            for line in metrics_handle:
                var_name, var_value, var_human_readable_name, var_tooltip, var_group = line.strip().split(",")
                self.metrics_dict[var_group][var_name] = (var_human_readable_name, self.format_number_to_string(var_value), var_tooltip)

        with open(sys.argv[4]) as html_template:
            self.jinja_template = Template(html_template.read())

        self.encoded_png_str = self.generate_encoded_png_str()

        self.render_and_write_report()

    def generate_encoded_png_str(self):
        with open(f"{self.plot_path}", "rb") as image_file:
            return base64.b64encode(image_file.read()).decode("utf-8")
        
    def render_and_write_report(self):
        # Render the template.        
        final_report = self.jinja_template.render(metrics_dict=self.metrics_dict,
                                            encoded_png_str=self.encoded_png_str,
                                            sample_id=self.sample_id,
                                            # This dict is required to supply the tooltips, the accordion header ID and the collapse ID for each of the categories of alignment statistics
                                            # The keys match the alignment categories of the metrics dict, and the values (tuples) give the tooltips and IDs.
                                            # We rely on the order of the dictionary to populate the accordions in the output html.
                                            alignment_cat_dict = {
                                                "Post read QC alignment": ("Mapping of the post QC reads i.e. after trimming (polyX end and internal polyA) and barcode verification.", "qc_accord_header", "qc_collapse"),
                                                "High confidence read alignment": ("Reads with a single alignment and a maximum of 3 bp mismatch.", "hc_accord_header", "hc_collapse"),
                                                "Annotated reads alignment": ("High confidence reads annotated with a gene ID (XT bam tag).", "ann_accord_header", "ann_collapse")
                                                })

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