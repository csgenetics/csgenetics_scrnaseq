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

This function reads in an anndata object (h5ad), get the total counts per cell, and 
log10 + 1 transforms the counts to find multiple troughs in the distribution of number
of counts per cell. The transformation is necessary to get the second peak.
When values are untransformed as there's too much spread.

The distribution is usually bi or trimodal with the last peak representing cells.

This function aims to find the local minima of the probabilty density
(i.e. the lowest point) between the last two peaks, and use that as the
counts threshold to call cells.

If there are no local minima above the minimum_count_threshold,
then the function will default to the minimum_count_threshold. This is generally set to 100 or 2
on the log10 scale.

If working with mixed species it identifies two sub-populations of barcodes: those with Mmus_counts > Hsap_counts 
(i.e. candidate mouse cells) and those with Hsap_counts > Mmus_counts (i.e. candidate human cells).

The total human counts per barcode in the majority human barcodes are transformed according to: log10(Hsap_counts + 1), 
and the local minmia of the probability density function is found.  The value of Hsap_counts at the minima is our Hsap_threshold.

The same process is repeated for the pool of candidate mouse cells, this time looking at the log-transformed distribution of mouse counts
to identify our Mmus_threshold.

In both cases, if there are no local minima above the minimum_count_threshold, then the function will default to the minimum_count_threshold. 
This is generally set to 100 or 2 on the log10 scale.

The combination of the two thresholds is used to classify barcodes according to the following rules in the annotate_cells_on_count_matrix.py script:

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
      
      args_dict = vars(parser.parse_args())
      for key in args_dict:
         setattr(self, key, args_dict[key])

      if self.single_species.lower() == "true":
         self.single_species = True
         self.mixed_species = False
      else:
         self.single_species = False
         self.mixed_species = True

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
            #self.hsap_log_cutoff = self.get_cutoff(self.hsap_pdf_df)
            #self.mmus_log_cutoff = self.get_cutoff(self.mmus_pdf_df)
            self.hsap_log_cutoff = 3.351725014652899
            self.mmus_log_cutoff = 3.2493156433793557
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
      # find the "dips" in the probability density - local minima
      local_minima = argrelextrema(pdf_df['evaluated'].values, np.less)
      # find the number of species-specific counts at which the minima occur
      potential_cutoffs = pdf_df['data_space'].values[local_minima]

      if len(np.where(potential_cutoffs > np.log10(self.minimum_count_threshold))[0]) == 0:
         log_cutoff = np.log10(self.minimum_count_threshold)
      else:
      # If there are any that are above 2, find the smallest one.
         log_cutoff = min(potential_cutoffs[np.where(potential_cutoffs>=np.log10(self.minimum_count_threshold))])
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
