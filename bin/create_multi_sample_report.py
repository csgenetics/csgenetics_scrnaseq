#!/usr/bin/python3

"""
Creates a multi sample report that summarises the single summary report metrics
for each sample.

It reads in a single csv per sample that contains the metrics.

The metrics are then converted into a dictionary where k=metric name
and v=a list of metric values, where the metric values are reported in
order of a sample_names list that is also used to render the template.
"""

from glob import glob
from jinja2 import Template
from pathlib import Path
from collections import defaultdict
import sys
from create_single_sample_report import get_cell_stat_cat_dict_obj

class MultipleSampleHTMLReport:
    def __init__(self):
        with open(sys.argv[1], "r") as html:
            self.template = Template(html.read())

        if sys.argv[2].upper() == "TRUE":
            self.mixed = True
        else:
            self.mixed = False

        self.sample_name_list, self.metrics_dict = self.make_metrics_dict_and_sample_list()
        self.render_and_write_report()
        self.write_out_csv()

    def write_out_csv(self):
        """
        write out the metrics as a csv
        """
        # open a stream to write to
        with open("multisample_out.csv", "w") as csv_out:
            # Put in an empty field to allow for the variable names
            csv_out.write("variable_name,human_readable_name,description,classification,")
            csv_out.write(",".join(self.sample_name_list) + "\n")
            for subcategory, subdict in self.metrics_dict.items():
                for (variable_name, human_readable_name, description, classification), metric_val_list in subdict.items():
                    csv_out.write(f"{','.join([variable_name, human_readable_name, description, classification])},")
                    csv_out.write(",".join(metric_val_list) + "\n")

    def render_and_write_report(self):
        # Render the template.
        # We pass the alignmnet_tooltip_dict so that we can provide the tooltips
        # for each of the alignment subcategories rather than for each
        # individual alignment metric.
        # NOTE the alignment card has temporarily been disabled by removing the
        # associated html. I am leaving the associated code in place so that it
        # can be reinstated.
        
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
                "High confidence read alignment": "Reads with a single alignment and a maximum of 3 bp mismatch.",
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
                    #   High confidence read alignment
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
    MultipleSampleHTMLReport()