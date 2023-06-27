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
        # final_report = self.jinja_template.render(metrics_dict=self.metrics_dict["cell_stats"],
        #                                           duplication_stats_dict=self.metrics_dict["duplication_stats"],
        #                                     sequencing_metrics_dict=self.metrics_dict["sequence_stats"],
        #                                     mapping_metrics_dict=self.metrics_dict["alignment_stats"],
        #                                     encoded_png_str=self.encoded_png_str,
        #                                     sample_id=self.sample_id,
        #                                     estimated_number_of_cells=self.metrics_dict["cell_stats"]["num_cells"],
        #                                     mean_genes_detected_per_cell=self.metrics_dict["cell_stats"]["mean_genes_detected_per_cell"],
        #                                     raw_reads_aligned=self.metrics_dict["alignment_stats"]["raw_reads_aligned"])
        
        final_report = self.jinja_template.render(metrics_dict=self.metrics_dict,
                                            encoded_png_str=self.encoded_png_str,
                                            sample_id=self.sample_id)

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