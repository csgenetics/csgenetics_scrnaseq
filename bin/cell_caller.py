#!/usr/bin/env python

import scanpy as sc 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde
import sys, argparse

"""
This is the cell caller function. It is counts based. It operates on both all counts (i.e. including mito).

This function reads in an anndata object (h5ad), get the total counts per cell, and 
log10 + 1 transforms the counts to find multiple troughs in the distribution of number
of counts per cell. The transformation is necessary to get the  second peak.
When values are untransformed as there's too much spread.

The distribution is usually bi or trimodal with the last peak representing cells.

This function aims to find the local minima of the probabilty density
(i.e. the lowest point) between the last two peaks, and use that as the
counts threshold to call cells.

If there are no local minima above the minimum_count_threshold,
then the function will default to the minimum_count_threshold. This is generally set to 100 or 2
on the log10 scale."""

class CellCaller:

   def __init__(self):
      self.parse_arguments()
      self.read_in_anndata_and_handle_error()
      self.derive_count_threshold()

   def parse_arguments(self):
      parser = argparse.ArgumentParser(description = "Arguments for cell caller script to calculate count threshold for calling cells")
      parser.add_argument("--sample_name", help="Sample ID", required=True)
      parser.add_argument("--minimum_count_threshold", default=100, type=float, help="minimal number of counts to call cell")
      parser.add_argument("--count_matrix", help="Path to the h5ad count matrix.", required=True)
      
      args_dict = vars(parser.parse_args())
      for key in args_dict:
         setattr(self, key, args_dict[key])

   def derive_count_threshold(self):
      """
      Derive the count threshold above which
      a barcode will be considered a cell
      """
      self.log10_counts = np.log10(self.adata.X.sum(axis=1).A1 + 1)

      if len(set(self.log10_counts)) == 1:
         self.clean_exit_on_error()
      else:
         self.pdf_df = self.get_prob_dens_data()
         self.log_cutoff = self.get_cutoff()
         self.make_pd_plots()

      # Transform back, and round to the nearest integer
      print(round(10 ** self.log_cutoff), end="")

   def clean_exit_on_error(self):
      """
      If we encounter an error we output an empty figure and return the
      default minimum_count_threshold
      """
      open(f"{self.sample_name}_pdf_with_cutoff.png", "w").close()
      print(int(self.minimum_count_threshold), end="")
      sys.exit(0)
      
   def read_in_anndata_and_handle_error(self):
      try:
         self.adata = sc.read_h5ad(self.count_matrix)
      except OSError:
         # If we encounter an empty h5ad then we output an empty figure and return the
         # default minimum_count_threshold
         self.clean_exit_on_error()

   def get_prob_dens_data(self):
      """
      Make a pdf and calculate probability density
      for the distribution of nuclear genes (log10 scale)
      """
      pdf = gaussian_kde(self.log10_counts)

      # return evenly spaced numbers over the number of nuclear genes range
      data_space = np.linspace(self.log10_counts.min(), self.log10_counts.max(), 200)

      # get probability values based on the pdf and the gene range
      evaluated = pdf.evaluate(data_space)
      
      # make a df
      pdf_df=pd.DataFrame({'data_space':data_space, 'evaluated':evaluated}) 
      return pdf_df

   def get_cutoff(self):
      """
      Calculate and return cutoff
      """
      # find the "dips" in the probability density - local minima
      local_minima = argrelextrema(self.pdf_df['evaluated'].values, np.less)
      # find the number of genes at which the minima occur
      potential_cutoffs = self.pdf_df['data_space'].values[local_minima]

      if len(np.where(potential_cutoffs > np.log10(self.minimum_count_threshold))[0]) == 0:
         log_cutoff = np.log10(self.minimum_count_threshold)
      else:
      # If there are any that are above 2, find the smallest one.
      # Transform back by 10^ to the value, and round to the nearest integer.
         log_cutoff = min(potential_cutoffs[np.where(potential_cutoffs>=np.log10(self.minimum_count_threshold))])
      return log_cutoff

   def make_pd_plots(self):
      """
      Plot probability density with the default cutoffs
      and cell caller cutoff, and save the plot as a png
      """
      plt.plot(self.pdf_df['data_space'], self.pdf_df['evaluated'])
      plt.axvline(self.log_cutoff, color = 'red', label = f"Cell Caller\n {str(round(10 ** self.log_cutoff))}")
      plt.axvline(np.log10(self.minimum_count_threshold), color = 'black', label = f"Default\n{str(round(self.minimum_count_threshold))}")
      plt.title("Cell Caller minimum count\nthreshold calculation", fontdict = {'family':'sans-serif','color':'black','size':12,'fontweight':'bold'})
      plt.xlabel("log10(total_counts + 1)")
      plt.ylabel("Density")
      plt.legend(bbox_to_anchor = (1.0, 1), loc = 'upper right', title="Method", fontsize=7)
      plt.savefig(f"{self.sample_name}_pdf_with_cutoff.png", dpi=900) 
      plt.close()    # close the figure window

if __name__ == "__main__":
   CellCaller()
