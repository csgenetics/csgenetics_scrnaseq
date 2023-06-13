#!/usr/bin/python3

from bs4 import BeautifulSoup as bs
import re
import pandas as pd
import sys
from jinja2 import Template
import base64

# This script will take CSV outputs for:

# - Cell statistics 
# - Mapping Statistics
# - Sequencing Statistics
# output from the web_sumamry.R script and insert them into the HTML for each sample.
# We also inserts a premade png for cell caller

sample_id = sys.argv[1]
plot_path = sys.argv[2]

# ---  Encode cell caller PNG to insert into HTML.

with open(f"{plot_path}", "rb") as image_file:
    encoded_string = base64.b64encode(image_file.read())
    decoded_string =encoded_string.decode("ascii")

image_html = f"<img alt='...' class='img-fluid' src='data:image/png;base64,{decoded_string}'/>"

# --- Read template HTML
html = open('template_report.html')
# Parse HTML file in Beautiful Soup
soup = bs(html, 'html.parser')


# --- Function will look for the td ID and insert a string
def insert_id(id, insert, tag1):
    old_text = soup.find(f"{tag1}", {"id": f"{id}"})
    # Replace all existing values with insert
    new_text = old_text.find(text=re.compile(
    '.*.')).replace_with(f"{insert}")

# --- This function will loop through the dataframe.

# The td ids are identical to the dataframe row indices.
# Therefore row indices are used as id, and corresponsing string is inserted
def loop_dat(dat,tag):
    for index,row in dat.iterrows():
        print(row[0]) # rowname / td id
        print(row[1]) # string to insert
        print(f"{tag}")
        insert_id(row[0],row[1],f"{tag}")


# --- Read data tables
mapping_dat = pd.read_csv(f"{sample_id}_map_stats.csv",
                          sep="\t",
                          header=None)

cell_dat = pd.read_csv(f"{sample_id}_cell_stats.csv",
                          sep="\t",
                          header=None)

cell_banner_dat = pd.read_csv(f"{sample_id}_cell_stats_banner.tmp",
                          sep="\t",
                          header=None)

seq_dat = pd.read_csv(f"{sample_id}_seq_stats.csv",
                          sep="\t",
                          header=None)

# --- Apply function to insert strings
loop_dat(mapping_dat, "td")
loop_dat(cell_dat,"td")
loop_dat(seq_dat, "td")
loop_dat(cell_banner_dat, "div")

# --- Insert sample_id and plot path
insert_id("sample_id", sample_id, "h2")


# --- Input encoded image for cell caller into html
html_table_jinj = Template(str(soup))
final_report = html_table_jinj.render(cell_caller_img=image_html)

# --- Wrtie
with open(f"{sample_id}_report.html", "w") as f_output:
    f_output.write(final_report)


