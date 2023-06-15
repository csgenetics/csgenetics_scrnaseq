#!/usr/bin/python3

"""
Create a multi sample report html.
Takes as input a set of .csv files (one per sample)
to be found in the current directory.

TODO check if we can simply stack the existing csv files.
"""

import pandas as pd
from glob import glob
import json
from jinja2 import Template

class SingleSampleHTMLReport:
    def __init__(self):
        self.table_template = self.return_table_template()
        self.make_multisample_df()
        self.json_data = self.create_write_final_df_json()

    def create_html_report(self):
        """
        Inject the json data into the html template
        and output the report
        """

        # Create the html to be inserted by injecting the json into
        # the table_template
        insert_html = self.table_template.render(users = self.json_data)

        # Generate a Jinja2 Template from the outer template html
        with open("multisample_template.html", "r") as html_txt:
            CS_template_html = Template(html_txt.read())

        # Insert table
        # TODO rewrite so that we are not doing a nested render.
        final_report = CS_template_html.render(table_insert=insert_html)

        # Write out the rendered template.
        with open(f"multisample_report.html", "w") as f_output:
            f_output.write(final_report)

    def make_multisample_df(self):
        """
        Create a single csv from the input csvs.
        """
        self.concatenated_df = self.make_concatenated_df()
        self.final_df = self.curate_concat_df()
        self.wite_final_df_csv()
        
    def wite_final_df_csv(self):
        self.final_df.to_csv("multisample_out.csv", sep=",")

    @staticmethod
    def make_concatenated_df():
        """
        Load in a df for each input csv
        and concatenate them together into a single dataframe
        """
        dfs = list()
        for file in glob('./*.csv'):
            dat = pd.read_csv(file, sep = "\t", index_col=None, header=None)
            dat['Path'] = file.stem.split("_")[0]
            dfs.append(dat)

        return pd.concat(dfs, ignore_index=True)

    def curate_concat_df(self):
        self.concatenated_df = self.concatenated_df.rename({0:"Statistic", 1:"value"}, axis=1)

        # Pivot table
        self.concatenated_df = self.concatenated_df.pivot(index = "Statistic",
                                    columns="Path",
                                    values="value")

        # Remove pivot indices
        self.concatenated_df = self.concatenated_df.rename_axis(None, axis=0) 

        # Make Sample ID the first row
        self.concatenated_df.loc['Sample ID'] = self.concatenated_df.columns
        index_no_sample = list(self.concatenated_df.index)
        index_no_sample.remove('Sample ID')
        new_index = ['Sample ID'] + index_no_sample
        return self.concatenated_df.loc[new_index]

    def create_write_final_df_json(self):
        # Write JSON file
        self.final_df.to_json("multisample_out.json", orient="index")
        with open("multisample_out.json", "r") as jsonhandle:
            return json.load(jsonhandle)
        
    @staticmethod
    def return_table_template():
        return Template("""

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
