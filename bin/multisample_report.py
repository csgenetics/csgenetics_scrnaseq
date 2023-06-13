#!/usr/bin/python3
import pandas as pd
import sys
from glob import glob
from pathlib import Path
import json
from jinja2 import Template


# ------------------ Load CSVs and concat
path = "."
csv_files = Path(path).glob('*.csv') 
dfs = list()

for file in csv_files:
    dat = pd.read_csv(file, sep = "\t", index_col=None, header=None)
    dat['Path'] = file.stem.split("_")[0]
    dfs.append(dat)


concat_dat = pd.concat(dfs,
                       ignore_index=True)

concat_dat = concat_dat.rename({0:"Statistics", 1:"value"}, axis=1)

# Pivot table
pivot_table = concat_dat.pivot(index = "Statistics",
                               columns="Path",
                               values="value")

# Remove pivot indices
rename_pivot_table = pivot_table.rename_axis(None, axis=0) 


# Make Sample ID the first row
rename_pivot_table.loc['Sample ID'] = rename_pivot_table.columns
index_no_sample = list(rename_pivot_table.index)
index_no_sample.remove('Sample ID')
new_index = ['Sample ID'] + index_no_sample
final_df = rename_pivot_table.loc[new_index]
# Write CSV file
final_df.to_csv("multisample_out.csv", sep=",")

# Write JSON file
final_df.to_json("multisample_out.json",orient="index")
json_file = open("multisample_out.json")
json_data = json.load(json_file)

# ------------------ HTML manipulation

# --- create jinja template for html table

# This loops through the values in the merged table and creates 
# table tags for values

table_template = Template("""

{% for key,value in users.items() %}
{% if key == "Sample ID" %}
<tr>
    <th >
        
     </th>

    {% for id,value2 in value.items() %}

     <td style="font-weight: bold">
      {{value2}}
     </td>

    {% endfor %}

{% else %}
<tr>
    <th >
        {{key}}
     </th>

{% for id,value2 in value.items() %}

     <td>
      {{value2}}
     </td>

    {% endfor %}
{% endif %}
</tr>

{% endfor %}
""")

# --- insert json data into table template
insert_html = table_template.render(users = json_data)

# ------- Open CS html to insert table and save

html_txt = open("multisample_template.html", "r").read()
CS_template_html = Template(html_txt)

# --- Insert table
final_report = CS_template_html.render(table_insert=insert_html)

#  --- write
with open(f"multisample_report.html", "w") as f_output:
    f_output.write(final_report)
