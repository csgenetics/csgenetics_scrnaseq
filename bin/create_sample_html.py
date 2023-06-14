#!/usr/bin/python3
from bs4 import BeautifulSoup as bs
import re
import pandas as pd
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

        self.mapping_data = pd.read_csv(f"{sys.argv[3]}",
                                sep="\t",
                                header=None)
        
        self.cell_data = pd.read_csv(f"{sys.argv[4]}",
                                sep="\t",
                                header=None)

        self.cell_banner_data = pd.read_csv(f"{sys.argv[5]}",
                                sep="\t",
                                header=None)
        
        self.seq_data = pd.read_csv(f"{sys.argv[6]}",
                                sep="\t",
                                header=None)

        with open(sys.argv[7]) as html_template_path:
            self.soup = bs(html_template_path, 'html.parser')

        self.jinja_template = Template(str(self.soup))

        self.cell_caller_html = self.generate_cell_caller_plot_html()

    def generate_cell_caller_plot_html(self):
        with open(f"{self.plot_path}", "rb") as image_file:
            encoded_png_str = base64.b64encode(image_file.read()).decode("utf-8")
            return f"<img alt='...' class='img-fluid' src='data:image/png;base64,{encoded_png_str}'/>"

    def generate_html_report(self):
        self.insert_sample_id()
        self.insert_cell_caller_plot()
        self.insert_stats()
        self.write_report()

    def insert_sample_id(self):
        self.insert_id("sample_id", self.sample_id, "h2")

    def insert_cell_caller_plot(self):
        self.report = self.jinja_template.render(cell_caller_img=self.cell_caller_html)

    def insert_stats(self):
        self.insert_mapping_data()
        self.insert_cell_data()
        self.insert_sequencing_data()
        self.insert_cell_banner_data()

    def insert_mapping_data(self):
        self.loop_dat(self.mapping_dat, "td")

    def insert_cell_data(self):
        self.loop_dat(self.cell_data, "td")

    def insert_sequencing_data(self):
        self.loop_dat(self.seq_data, "td")

    def insert_cell_banner_data(self):
        self.loop_dat(self.cell_banner_data, "div")

    def write_report(self):
        with open(f"{self.sample_id}_report.html", "w") as report_out:
            report_out.write(self.report)
        
    def insert_id(self, id, insert, tag1):
        """
        Function will look for the td ID and insert a string
        """
        old_text = self.soup.find(f"{tag1}", {"id": f"{id}"})
        # Replace all existing values with insert
        new_text = old_text.find(text=re.compile(
        '.*.')).replace_with(f"{insert}")

    # --- This function will loop through the dataframe.
    # The td ids are identical to the dataframe row indices.
    # Therefore row indices are used as id, and corresponsing string is inserted
    def loop_dat(self, dat: pd.DataFrame, tag: str):
        """
        TODO document what form the df is in.
        """
        for index, row in dat.iterrows():
            print(row[0]) # rowname / td id
            print(row[1]) # string to insert
            print(f"{tag}")
            self.insert_id(row[0], row[1], f"{tag}")

if __name__ == "__main__":
    SingleSampleHTMLReport().generate_html_report()