#!/usr/bin/env python

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde
import sys, argparse

"""
This is the cell caller function.
It operates on the distribution of total counts for cellular barcodes.

This function reads in an anndata object (h5ad), gets the total counts per cellular barcode,  
log10 + 1 transforms the counts and evaluates the density distribution of those transformed counts.
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
This is generally set to 100 i.e. 2 on the log10 scale.

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

class CellCaller:

   def __init__(self):
      self.parse_arguments()
      self.read_in_anndata_and_handle_error()
      self.derive_count_threshold()

   def parse_arguments(self):
      parser = argparse.ArgumentParser(description = "Arguments for cell caller script to calculate count threshold for calling cells")
      parser.add_argument("--sample_name", help="Sample ID", required=True)
      parser.add_argument("--single_species", help="Whether this is a single species run.", required=True)
      parser.add_argument("--minimum_count_threshold", default=100, type=float, help="minimal number of counts to call cell")
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

         if len(set(self.log10_counts)) == 1:
            self.clean_exit_on_error()
         else:
            self.pdf_df = self.get_prob_dens_data(self.log10_counts)
            if self.single_species_specified_threshold:
               self.log_cutoff = self.single_species_specified_threshold
            else:
               self.log_cutoff = self.get_cutoff(self.pdf_df)

            self.make_pd_plots()
            # Transform back, and round to the nearest integer
            print(f"{round(10 ** self.log_cutoff)}", end="")
      else:
         # Identify the sub-populations of barcodes which are: 1) majority human counts and 2) majority mouse counts
         self.hsap_majority = self.adata.obs[self.adata.obs["hsap_counts"] > self.adata.obs["mmus_counts"]]
         self.mmus_majority = self.adata.obs[self.adata.obs["mmus_counts"] > self.adata.obs["hsap_counts"]]

         # Log transform (+1) the total human counts in the majority human population, and the total mouse counts in the majority mouse population  
         self.log10_hsap_counts_hsap_majority = np.log10(self.hsap_majority.hsap_counts+1)
         self.log10_mmus_counts_mmus_majority = np.log10(self.mmus_majority.mmus_counts+1)

         # Compute the human and mouse thresholds on a log10 scale, and output the probability density plots
         if (len(set(self.log10_hsap_counts_hsap_majority)) == 1) or (len(set(self.log10_mmus_counts_mmus_majority)) == 1):
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

            self.make_pd_plots()
            hsap_thres = round(10 ** self.hsap_log_cutoff)
            mmus_thres = round(10 ** self.mmus_log_cutoff)
            # Transform back, and round to the nearest integer
            print(f"{hsap_thres}_{mmus_thres}", end="")

   def clean_exit_on_error(self):
      """
      If we encounter an error we output an empty figure and return the
      default minimum_count_threshold in either single or mixed species format
      """
      open(f"{self.sample_name}_pdf_with_cutoff.png", "w").close()
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
      # **********************
      # MINIMA IDENTIFICATION
      # **********************
      # Find the local minima of the probability density curve
      local_minima = argrelextrema(pdf_df['evaluated'].values, np.less, order=1)
      # Find the number of species-specific counts at which the minima occur
      minima_locations = pdf_df['data_space'].values[local_minima]
      # Filter to retain only the minima above the minimum_count_threshold
      potential_minima_cutoffs = minima_locations[minima_locations >= np.log10(self.minimum_count_threshold)]

      # ********************************
      # INFLECTION POINT IDENTIFICATION
      # ********************************
      # Calculate the first and second derivatives of the probability density curve to use in identifying suitable inflection points
      first_derivatives = np.diff(pdf_df['evaluated'].values)/np.diff(pdf_df['data_space'].values)
      second_derivatives = np.diff(first_derivatives)/np.diff(pdf_df['data_space'].values[:-1])
      derivatives_df = pd.DataFrame({"data_space":pdf_df['data_space'].values[:-2], "first_derivative":first_derivatives[:-1], "second_derivative":second_derivatives})
      # Add a binary flag to the derivatives df to indicate whether the sign of the second derivative changes from positive to negative between consecutive points (i.e. an inflection point)
      derivatives_df["sign_change"] = np.where(derivatives_df["second_derivative"].shift(1) * derivatives_df["second_derivative"] < 0, 1, 0)
      
      # Find the inflection points which are above the minimum_count_threshold, which also have a gradient close to 0 
      inflection_points_df = derivatives_df[(derivatives_df["sign_change"] == 1) & (derivatives_df["data_space"] >= np.log10(self.minimum_count_threshold)) & (derivatives_df["first_derivative"].abs() < 0.2)]
      # Find the number of species-specific counts at which the inflection point(s) occur
      potential_inflection_point_cutoffs = inflection_points_df["data_space"].values

      # **********************
      # THRESHOLD SELECTION
      # **********************

      if len(potential_minima_cutoffs) == 0:
         # If there are no suitable minima but there is a back-up inflection point, select the inflection point
         # Otherwise, default to the minimum_count_threshold
         if len(potential_inflection_point_cutoffs) == 0:
            log_cutoff = np.log10(self.minimum_count_threshold)
         else:
            log_cutoff = min(potential_inflection_point_cutoffs)
      else:
         # If there are any minima that are above 2, find the smallest one.
         log_cutoff = min(potential_minima_cutoffs)
      return log_cutoff

   def make_pd_plots(self):
      """
      # TODO convert to plotly.
      Plot probability density with the default cutoffs
      and cell caller cutoff, and save the plot as a png
      """
      if self.single_species:
         plt.plot(self.pdf_df['data_space'], self.pdf_df['evaluated'])
         plt.axvline(self.log_cutoff, color = 'red', label = f"Cell Caller\n {str(round(10 ** self.log_cutoff))}")
         plt.axvline(np.log10(self.minimum_count_threshold), color = 'black', label = f"Default\n{str(round(self.minimum_count_threshold))}")
         plt.title("Cell Caller minimum count\nthreshold calculation", fontdict = {'family':'sans-serif','color':'black','size':12,'fontweight':'bold'})
         plt.xlabel("log10(total_counts + 1)")
         plt.ylabel("Density")
         plt.legend(bbox_to_anchor = (1.0, 1), loc = 'upper right', title="Method", fontsize=7)
         plt.savefig(f"{self.sample_name}_pdf_with_cutoff.png", dpi=900) 
         plt.close()    # close the figure window
      else:
         # Create a new figure
         plt.figure(figsize=(7, 10))

         # First subplot for hsap
         plt.subplot(2, 1, 1)  # (rows, columns, panel number)
         plt.plot(self.hsap_pdf_df['data_space'], self.hsap_pdf_df['evaluated'])
         plt.axvline(self.hsap_log_cutoff, color='red', label=f"Cell Caller\n {str(round(10 ** self.hsap_log_cutoff))}")
         plt.axvline(np.log10(self.minimum_count_threshold), color='black', label=f"Default\n{str(round(self.minimum_count_threshold))}")
         plt.title("Cell Caller minimum human count threshold calculation\nin majority human-count barcodes", fontdict={'family': 'sans-serif', 'color': 'black', 'size': 12, 'fontweight': 'bold'})
         plt.xlabel("log10(hsap_counts + 1)")
         plt.ylabel("Density")
         plt.legend(bbox_to_anchor=(1.0, 1), loc='upper right', title="Method", fontsize=7)

         # Second subplot for mmus
         plt.subplot(2, 1, 2)  # (rows, columns, panel number)
         plt.plot(self.mmus_pdf_df['data_space'], self.mmus_pdf_df['evaluated'])
         plt.axvline(self.mmus_log_cutoff, color='red', label=f"Cell Caller\n {str(round(10 ** self.mmus_log_cutoff))}")
         plt.axvline(np.log10(self.minimum_count_threshold), color='black', label=f"Default\n{str(round(self.minimum_count_threshold))}")
         plt.title("Cell Caller minimum mouse count threshold calculation\nin majority mouse-count barcodes", fontdict={'family': 'sans-serif', 'color': 'black', 'size': 12, 'fontweight': 'bold'})
         plt.xlabel("log10(mmus_counts + 1)")
         plt.ylabel("Density")
         plt.legend(bbox_to_anchor=(1.0, 1), loc='upper right', title="Method", fontsize=7)

         # Adjust layout to prevent overlapping content
         plt.tight_layout()

         # Save the figure as a single .png image
         plt.savefig(f"{self.sample_name}_combined_pdf_with_cutoff.png", dpi=900)

         # Close the figure window to free up memory
         plt.close()

if __name__ == "__main__":
   CellCaller()
