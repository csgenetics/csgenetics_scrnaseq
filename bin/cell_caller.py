#!/usr/bin/env python

import anndata as ad
import numpy as np
import pandas as pd
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde
import sys, argparse
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio

"""
This is the cell caller function.
It operates on the distribution of total counts for cellular barcodes.

This function reads in an anndata object (h5ad), gets the total counts per cellular barcode,
transforms the counts according to log10(counts + 1) and evaluates the density distribution of those transformed counts.
The transformation is necessary to condense the data enough to get a noise peak and a single-cell peak.
When values are untransformed, there's too much spread to identify the populations.

The cannonical profile that the cell caller assumes is a bimodal distribution - with the rightmost peak representing the cells.

This function aims to find the local minima of the probabilty density (i.e. the lowest point) between the 
two peaks, and use that as the counts threshold to call cells. That is to say that any cellular barcode
with counts below this threshold would be called as noise, and any barcode with counts above this threshold
would be called as a cell.

If there are no local minima above the minimum_count_threshold, then the function will 
next look for inflection points above the minimum_count_threshold where the gradient is 
close to 0. If there are none of these either, then the function will default to the minimum_count_threshold. 
This is generally set to 100.

If working with mixed species it identifies two sub-populations of barcodes: those with Mmus_counts > Hsap_counts 
(i.e. candidate mouse cells) and those with Hsap_counts > Mmus_counts (i.e. candidate human cells).

The total human counts per barcode in the majority human barcodes are transformed according to: log10(Hsap_counts + 1), 
and the local minmia of the probability density function is found.  The value of Hsap_counts at the minima is our Hsap_threshold.

The same process is repeated for the pool of candidate mouse cells, this time looking at the log-transformed distribution of mouse counts
to identify our Mmus_threshold.

In both cases, if there are no local minima above the minimum_count_threshold, then the function will next look for a 
suitable inflection point to use instead, and if it fails to find one it will finally default to the minimum_count_threshold. 

The combination of the two thresholds is used to classify barcodes according to the following rules:

Hsap_counts >= Hsap_threshold & Mmus_counts <= Mmus_threshold is a single-cell (human).
Mmus_counts >= Mmus_threshold & Hsap_counts <= Hsap_threshold is a single-cell (mouse).
Mmus_counts < Mmus_threshold & Hsap_counts < Hsap_threshold is a noisy barcode.
Mmus_counts > Mmus_threshold & Hsap_counts > Hsap_threshold is a multiplet.
"""

# The plot to html method is defined outside the class as it will also be used later in the multi-sample summary report process. 
def output_plot_to_html(dict_of_figs_and_names, html_filename):
   """
   Output a suite of plots to an html file, incorporating the Lexend font.
   Expects as input a dictionary of plotly figures and their names.
   Figure names are used to name the .svg files which can be downloaded from the html.
   """
   # Add custom CSS to embed the Lexend font
   html_content = """
   <head>
      <link href="https://fonts.googleapis.com/css2?family=Lexend:wght@400;700&display=swap" rel="stylesheet">
   </head>
   <body>
   """

   # Add each figure to the HTML content
   for fig_name in dict_of_figs_and_names.keys():
      fig = dict_of_figs_and_names[fig_name]
      html_content += pio.to_html(fig, full_html=False, include_plotlyjs='cdn', config={"responsive": True, 
                                                                                          'displaylogo': False, 
                                                                                          'toImageButtonOptions': {'format': 'svg',
                                                                                                                  'filename': fig_name,
                                                                                                                  'scale': 1 
                                                                                                                  }
                                                                                       })

   # Close the HTML tags
   html_content += """
   </body>
   """

   # Write the HTML content to the file
   with open(html_filename, "w") as f:
      f.write(html_content)
class CellCaller:

   def __init__(self):
      self.parse_arguments()
      self.read_in_anndata_and_handle_error()
      self.derive_count_threshold()
      self.generate_pdf_plot()
      self.generate_barnyard_plot()

   def parse_arguments(self):
      parser = argparse.ArgumentParser(description = "Arguments for cell caller script to calculate count threshold for calling cells.")
      parser.add_argument("--sample_name", help="Sample ID", required=True)
      parser.add_argument("--single_species", help="Whether this is a single species run.", required=True)
      parser.add_argument("--minimum_count_threshold", default=100, type=float, help="Default minimum number of counts to call a barcode a cell in the absence of a dynamically determined threshold.")
      parser.add_argument("--count_matrix", help="Path to the h5ad count matrix.", required=True)
      parser.add_argument("--manual_threshold_str", help="A string of manual thresholds as specified by the user in the optional threshold input csv.", required=True)
      
      args_dict = vars(parser.parse_args())
      for key in args_dict:
         setattr(self, key, args_dict[key])

      if self.single_species.lower() == "true":
         self.single_species = True
         self.mixed_species = False
      else:
         self.single_species = False
         self.mixed_species = True

      # Extract any user specified cell caller thresholds
      if self.mixed_species:
         # Split the input thresholds string on _, the first value is the human threshold, the second is the mouse threshold
         self.manual_thresholds = self.manual_threshold_str.split("_")
         if self.manual_thresholds[0] == "nan":
            self.mixed_species_specified_threshold_hsap = None
         else:
            self.mixed_species_specified_threshold_hsap = float(self.manual_thresholds[0])
         if self.manual_thresholds[1] == "nan":
            self.mixed_species_specified_threshold_mmus = None
         else:
            self.mixed_species_specified_threshold_mmus = float(self.manual_thresholds[1])
      else:
         if self.manual_threshold_str == "nan":
            self.single_species_specified_threshold = None
         else:
            self.single_species_specified_threshold = float(self.manual_threshold_str)

   def derive_count_threshold(self):
      """
      Derive the count threshold above which
      a barcode will be considered a cell
      """
      if self.single_species:
         self.log10_counts = np.log10(self.adata.X.sum(axis=1).A1 + 1)

         if len(set(self.log10_counts)) <= 1:
            self.clean_exit_on_error()
         else:
            self.pdf_df = self.get_prob_dens_data(self.log10_counts)
            if self.single_species_specified_threshold:
               self.log_cutoff = self.single_species_specified_threshold
            else:
               self.log_cutoff = self.get_cutoff(self.pdf_df)
               
            # Transform back, and round to the nearest integer
            self.single_species_thres = round(10 ** self.log_cutoff -1)
            print(f"{self.single_species_thres}", end="") 
      else:
         # Identify the sub-populations of barcodes which are: 1) majority human counts and 2) majority mouse counts
         self.hsap_majority = self.adata.obs[self.adata.obs["hsap_counts"] > self.adata.obs["mmus_counts"]]
         self.mmus_majority = self.adata.obs[self.adata.obs["mmus_counts"] > self.adata.obs["hsap_counts"]]

         # Log transform the total human counts (+1) in the majority human population, and the total mouse counts (+1) in the majority mouse population  
         self.log10_hsap_counts_hsap_majority = np.log10(self.hsap_majority.hsap_counts+1)
         self.log10_mmus_counts_mmus_majority = np.log10(self.mmus_majority.mmus_counts+1)

         # Compute the human and mouse thresholds on a log10 scale, and output the probability density plots
         if (len(set(self.log10_hsap_counts_hsap_majority)) <= 1) or (len(set(self.log10_mmus_counts_mmus_majority)) <= 1):
            self.clean_exit_on_error()
         else:
            self.hsap_pdf_df = self.get_prob_dens_data(self.log10_hsap_counts_hsap_majority)
            self.mmus_pdf_df = self.get_prob_dens_data(self.log10_mmus_counts_mmus_majority)
            if self.mixed_species_specified_threshold_hsap:
               self.hsap_log_cutoff = self.mixed_species_specified_threshold_hsap
            else: 
               self.hsap_log_cutoff = self.get_cutoff(self.hsap_pdf_df)
            if self.mixed_species_specified_threshold_mmus:
               self.mmus_log_cutoff = self.mixed_species_specified_threshold_mmus
            else:
               self.mmus_log_cutoff = self.get_cutoff(self.mmus_pdf_df)

            # Transform back, and round to the nearest integer
            self.hsap_thres = round(10 ** self.hsap_log_cutoff -1)
            self.mmus_thres = round(10 ** self.mmus_log_cutoff -1)
            
            print(f"{self.hsap_thres}_{self.mmus_thres}", end="")

   def clean_exit_on_error(self):
      """
      If we encounter an error, we output empty figures and return the
      default minimum_count_threshold in either single or mixed species format
      """
      open(f"{self.sample_name}_counts_pdf_with_threshold.html", "w").close()
      open(f"{self.sample_name}_barnyard_plot.html", "w").close()
      if self.single_species:
         print(int(self.minimum_count_threshold), end="")
      else:
         print(f"{int(self.minimum_count_threshold)}_{int(self.minimum_count_threshold)}", end="")
      sys.exit(0)
      
   def read_in_anndata_and_handle_error(self):
      try:
         self.adata = ad.read_h5ad(self.count_matrix)
      except OSError:
         # If we encounter an empty h5ad then we output an empty figure and return the
         # default minimum_count_threshold
         self.clean_exit_on_error()

   def get_prob_dens_data(self, counts):
      """
      Make a pdf and calculate probability density
      for the distribution of species-specific counts (log10 scale)
      """
      pdf = gaussian_kde(counts)

      # return evenly spaced numbers over the number of species-specific counts range
      data_space = np.linspace(counts.min(), counts.max(), 200)

      # get probability values based on the pdf and the counts range
      evaluated = pdf.evaluate(data_space)
      
      # make a df
      pdf_df=pd.DataFrame({'data_space':data_space, 'evaluated':evaluated}) 
      return pdf_df

   def get_cutoff(self, pdf_df):
      """
      Calculate and return cutoff
      """
      # ----------------------
      # MINIMA IDENTIFICATION
      # ----------------------
      # Find the local minima of the probability density curve
      local_minima = argrelextrema(pdf_df['evaluated'].values, np.less, order=1)
      # Find the number of species-specific counts at which the minima occur
      minima_locations = pdf_df['data_space'].values[local_minima]
      # Filter to retain only the minima above the minimum_count_threshold
      potential_minima_cutoffs = minima_locations[minima_locations >= np.log10(self.minimum_count_threshold + 1)]

      # --------------------------------
      # INFLECTION POINT IDENTIFICATION
      # --------------------------------
      # Calculate the first and second derivatives of the probability density curve to use in identifying suitable inflection points
      first_derivatives = np.diff(pdf_df['evaluated'].values)/np.diff(pdf_df['data_space'].values)
      second_derivatives = np.diff(first_derivatives)/np.diff(pdf_df['data_space'].values[:-1])
      derivatives_df = pd.DataFrame({"data_space":pdf_df['data_space'].values[1:-1], "first_derivative":first_derivatives[:-1], "second_derivative":second_derivatives})
      # Add a binary flag to the derivatives df to indicate whether the sign of the second derivative changes from positive to negative between consecutive points (i.e. an inflection point)
      derivatives_df["sign_change"] = np.where(derivatives_df["second_derivative"].shift(1) * derivatives_df["second_derivative"] < 0, 1, 0)
      
      # Find the inflection points which are above the minimum_count_threshold, which also have a gradient close to 0 
      inflection_points_df = derivatives_df[(derivatives_df["sign_change"] == 1) & (derivatives_df["data_space"] >= np.log10(self.minimum_count_threshold + 1)) & (derivatives_df["first_derivative"].abs() < 0.2)]
      # Find the number of species-specific counts at which the inflection point(s) occur
      potential_inflection_point_cutoffs = inflection_points_df["data_space"].values

      # --------------------
      # THRESHOLD SELECTION
      # --------------------

      if len(potential_minima_cutoffs) == 0:
         # If there are no suitable minima but there is a back-up inflection point, select the inflection point
         # Otherwise, default to the minimum_count_threshold
         if len(potential_inflection_point_cutoffs) == 0:
            log_cutoff = np.log10(self.minimum_count_threshold + 1)
         else:
            log_cutoff = min(potential_inflection_point_cutoffs)
      else:
         # If there are any minima that are above the minimum threshold, find the smallest one.
         log_cutoff = min(potential_minima_cutoffs)
      return log_cutoff
         
   def generate_pdf_plot(self):
      pdf_html_filename = f"{self.sample_name}_counts_pdf_with_threshold.html"
      if self.single_species:
         total_counts_pdf_fig = self.pdf_plotter(self.pdf_df, self.log_cutoff, "total")
         output_plot_to_html({f"{self.sample_name}_cellcaller_plot":total_counts_pdf_fig}, pdf_html_filename)
      else:
         human_counts_pdf_fig = self.pdf_plotter(self.hsap_pdf_df, self.hsap_log_cutoff, "human")
         mouse_counts_pdf_fig = self.pdf_plotter(self.mmus_pdf_df, self.mmus_log_cutoff, "mouse")
         output_plot_to_html({f"{self.sample_name}_hsap_cellcaller_plot":human_counts_pdf_fig, f"{self.sample_name}_mmus_cellcaller_plot":mouse_counts_pdf_fig}, 
                                   pdf_html_filename)

   def pdf_plotter(self, input_pdf_df, input_log_cutoff, count_type_str):
      """
      Plot probability density with the cell caller cutoff, 
      and save the plot as an html (two plots for mixed species)
      """

      title_str = f"Cell caller minimum {count_type_str}-count threshold derived from <br> distribution of {count_type_str} counts"
      if self.mixed_species:
         title_str += f" in majority {count_type_str}-count barcodes"

      pdf_fig = px.line(input_pdf_df, 
                        x="data_space", 
                        y="evaluated",
                        labels={"data_space":f"log10({count_type_str}_counts + 1)", "evaluated":"Density"},
                        title=title_str
                        )
      
      pdf_fig.update_traces(line_color='#0081D4', 
                            line=dict(width=3),
                            customdata = [10**x -1 for x in input_pdf_df["data_space"]],
                            hovertemplate=f"log10({count_type_str}_counts + 1): %{{x:.3f}}<br>{count_type_str}_counts: %{{customdata:.0f}}<br>Density: %{{y:.3f}}"
                            )
      
      pdf_fig.add_vline(x=input_log_cutoff, 
                        line_color="#36BA00", 
                        line_width=3
                        )
      
      pdf_fig.add_trace(go.Scatter(x=[input_log_cutoff]*200,
                                   y=np.linspace(input_pdf_df["evaluated"].min(), input_pdf_df["evaluated"].max(), 200),
                                   mode='lines',
                                   name=f"Cell caller threshold =  <br>{str(round(10**input_log_cutoff -1))} {count_type_str} counts",
                                   line=dict(color = "#36BA00", width=3),  
                                   hovertemplate=f"Cell caller threshold = {round(input_log_cutoff,3)} (log10 scale)<extra></extra>"
                                   )
                        )
            
      pdf_fig.update_layout(font=dict(family = 'Lexend, sans-serif', color="black"),
                            title_x=0.5,
                            plot_bgcolor="white",
                            autosize=True,
                            xaxis = dict(
                               mirror=True,
                               ticks='outside',
                               showline=True,
                               linecolor='black',
                               linewidth=3,
                               tickwidth=3,
                               tickformat=".1f"
                               ),
                            yaxis = dict(
                               mirror=True,
                               ticks='outside',
                               showline=True,
                               linecolor='black',
                               linewidth=3,
                               tickwidth=3,
                               tickformat=".1f"
                               ),
                            legend=dict(
                               yanchor="top", y=0.99, 
                               xanchor="right",x=0.99, 
                               bordercolor="Black",
                               borderwidth=2,
                               title=""
                               )
                           )

      return pdf_fig

   def assign_barcode_type(self, hsap_counts, mmus_counts):
      """
      In a mixed-species experiment, barcode-types are assigned based on the following rules:

      Hsap_counts >= Hsap_threshold & Mmus_counts <= Mmus_threshold is a single-cell (human).
      Mmus_counts >= Mmus_threshold & Hsap_counts <= Hsap_threshold is a single-cell (mouse).
      Mmus_counts < Mmus_threshold & Hsap_counts < Hsap_threshold is a noisy barcode.
      Mmus_counts > Mmus_threshold & Hsap_counts > Hsap_threshold is a multiplet.
      """

      if hsap_counts > self.hsap_thres and mmus_counts > self.mmus_thres:
         barcode_type = "multiplet" 
      elif hsap_counts < self.hsap_thres and mmus_counts < self.mmus_thres:
         barcode_type = "noise"
      else:
         barcode_type = "single-cell"
      return barcode_type
      
   def generate_barnyard_plot(self):
      """
      Create barnyard plot for mixed species samples. If single species then output an empty html. 
      """
      if self.single_species:
         open(f"{self.sample_name}_barnyard_plot.html", "w").close()
      else:
         # Create a new column in the adata object to assign barcode types
         self.adata.obs["barcode_type"] = self.adata.obs.apply(lambda x:self.assign_barcode_type(x["hsap_counts"], x["mmus_counts"]), axis=1)
         self.adata.obs["cell_name"] = self.adata.obs.index

         # Set axis upper limit to just higher than the max count of either species 
         axis_upper_limit = max([self.adata.obs["hsap_counts"].max(), self.adata.obs["mmus_counts"].max()])*1.01
         # Set axis lower limit to a fixed % below 0 - allows some separation from the axis for better readability
         axis_lower_limit = 0-(axis_upper_limit*0.02)

         barnyard_fig = px.scatter(self.adata.obs, 
                                 x="mmus_counts", 
                                 y="hsap_counts", 
                                 color="barcode_type", 
                                 labels={"hsap_counts":"Number of counts from human genes", "mmus_counts":"Number of counts from mouse genes", "barcode_type":""},
                                 title="Barnyard plot of human and mouse counts per barcode",
                                 opacity=0.7,
                                 color_discrete_map={"noise": "#9AA3A8", "single-cell":"#36BA00", "multiplet":"#0081D4"},
                                 custom_data=['cell_name', 'barcode_type']
                                 )

         barnyard_fig.update_traces(hovertemplate=(
                                       "Cell Name: %{customdata[0]}<br>"
                                       "Barcode Type: %{customdata[1]}<br>"
                                       "Mouse Counts: %{x:.0f}<br>"
                                       "Human Counts: %{y:.0f}"
                                       )
                                    )
      
         barnyard_fig.add_hline(y=self.hsap_thres, 
                              line_dash="dot", 
                              line_color="black",
                              line_width = 1.5,
                              annotation_position = "bottom right",
                              annotation_text = f"Hsap threshold = {self.hsap_thres}"
                              )

         # Add the threshold legend to the graph
         barnyard_fig.add_trace(go.Scatter(y=[self.hsap_thres, self.hsap_thres],
                                          x=[0, 0.1],
                                          mode='lines',
                                          name=f"Hsap threshold = {self.hsap_thres}",
                                          line=dict(dash='dot', color = "black", width=1.5))
                                 )
         
         barnyard_fig.add_vline(x=self.mmus_thres, 
                           line_dash="dash", 
                           line_color="black",
                           line_width = 1.5,
                           annotation_position = "top right",
                           annotation_text = f"Mmus threshold = {self.mmus_thres}"
                           )
         
         barnyard_fig.add_trace(go.Scatter(x=[self.mmus_thres, self.mmus_thres],
                                          y=[0, 0.1],
                                          mode='lines',
                                          name=f"Mmus threshold = {self.mmus_thres}",
                                          line=dict(dash='dash', color = "black", width=1.5))
                                 )
      
         barnyard_fig.update_layout(font=dict(family = 'Lexend, sans-serif', color="black"),
                                    plot_bgcolor="white",
                                    margin = dict(l=45,r=45,t=50,b=45),
                                    autosize=True,
                                    xaxis = dict(
                                       range=[axis_lower_limit, axis_upper_limit],
                                       mirror=True,
                                       ticks='outside',
                                       showline=True,
                                       linecolor='black',
                                       linewidth=2,
                                       tickwidth=2,
                                       gridcolor='lightgrey'
                                       ),
                                    yaxis = dict(
                                       range=[axis_lower_limit, axis_upper_limit],
                                       mirror=True,
                                       ticks='outside',
                                       showline=True,
                                       linecolor='black',
                                       linewidth=2,
                                       tickwidth=2,
                                       gridcolor='lightgrey'
                                       )
                                    )

         barnyard_html_filename = f"{self.sample_name}_barnyard_plot.html"

         output_plot_to_html({f"{self.sample_name}_barnyard_plot":barnyard_fig}, barnyard_html_filename)
      
if __name__ == "__main__":
   CellCaller()