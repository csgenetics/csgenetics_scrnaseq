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

class MultipleSampleHTMLReport:
    def __init__(self):
        with open(sys.argv[1], "r") as html:
            self.template = Template(html.read())

        self.sample_name_list, self.metrics_dict = self.make_metrics_dict_and_sample_list()
        self.render_and_write_report()
        self.write_out_tsv()

    def write_out_tsv(self):
        """
        write out the metrics as a tsv
        """
        # open a stream to write to
        with open("multisample_out.tsv", "w") as tsv_out:
            tsv_out.write("\t")
            tsv_out.write("\t".join(self.sample_name_list) + "\n")
            for metric_name, metric_val_list in self.metrics_dict.items():
                tsv_out.write(f"{metric_name}\t")
                tsv_out.write("\t".join(metric_val_list) + "\n")

    def render_and_write_report(self):
        # Render the template.
        final_report = self.template.render(sample_name_list=self.sample_name_list, metrics_dict=self.metrics_dict)

        # Write out the rendered template.
        with open(f"multisample_report.html", "w") as f_output:
            f_output.write(final_report)

    def make_metrics_dict_and_sample_list(self):
        """
        Create the two data structures that will be used to render the template.
        self.sample_list is simply a list that contains the names of the samples
        sorted in some sort of order.
        self.metrics_dict is a defaultdict(list) key is metric name, value is list
        containing the values for that metrics for each of the samples in the
        same sample order as self.sample_list.
        """

        self.metrics_dict = defaultdict(list)
        self.sample_name_list = []

        # Get all the csv paths and sort them
        # so that the samples will be reported in
        # a logical order
        csvs = [_ for _ in Path(".").glob('./*.csv')]
        csvs.sort()

        # For each sample csv in the sorted paths
        for csv_path in csvs:
            self.sample_name_list.append(csv_path.stem.replace(".metrics", ""))
            with open(csv_path, "r") as csv_f:
                for line in [_.rstrip() for _ in csv_f]:
                    # read in the metrics and add them to the defaultdict
                    metric_name, metric_val,table_group = line.split(",")
                    self.metrics_dict[metric_name].append(metric_val.strip())

        return self.sample_name_list, self.metrics_dict

if __name__ == "__main__":
    MultipleSampleHTMLReport()