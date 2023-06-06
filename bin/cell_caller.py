#!/usr/bin/env python

import scanpy as sc 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde
import re
import sys, argparse

"""This is the cell caller function. This function reads in an anndata object (h5ad), counts total number of genes per cells, and 
log10 transforms the counts to see multiple peaks in the distribution of number of genes per cell (no second peak when values are
untransformed as there's too much spread).  From looking at the density plots of the number of genes (on a log 10 scale), it is
generally fairly obvious that there are multiple peaks. In the more recent versions of the assay, the distribution is usually
bimodal - one peak for "noise", one peak for "cells". This function aims to find the local minima of the probabilty density
(i.e. the lowest point) between the two peaks, and use that as the nuclear gene threshold to call cells. This will be called for
every sample, rather than having a fixed 100 gene threshold for everything. If there are no local minima above 100 nuclear genes,
then the function will default to 100 - i.e. 2 on the log10 scale, as log10(100) = 2."""

def parse_arguments(args):
    parser = argparse.ArgumentParser(description = "Arguments for cell caller script to calculate number of genes threshold")
    parser.add_argument("--sample", help="Sample ID")
    parser.add_argument("--min_nucGene", default=100, type=int, help="minimal number of nuclear gene to call single cell")

    return parser.parse_args()


def getlog10NucGenes(sample):
   '''function to read in sample and get relevant values into adata.obs'''
   #read in the counts matrix
   adata = sc.read_h5ad('{sample}_tmp.h5ad'.format(sample=sample))
   #calculate QC metrics (total number of genes per cell, total counts per cell)
   sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
   # calculate the number of NUCLEAR genes per cell
   adata.obs['nNuc_genes'] = adata.X[:,~adata.var_names.str.contains("mt-", flags=re.IGNORECASE, regex=True)].toarray().astype(bool).sum(axis=1)
   # get a log 10 of the number of genes  -  need +1 as some values in nNuc_genes are 0 
   adata.obs['log10_Nuc_genes'] = np.log10(adata.obs['nNuc_genes'] +1)
   log10_Nuc_genes = adata.obs['log10_Nuc_genes'].to_numpy()
   return log10_Nuc_genes

def get_prob_dens_data(log10_Nuc_genes):
   '''function to make a pdf and calculate probability density for the distribution of nuclear genes (log10 scale)'''
   pdf = gaussian_kde(log10_Nuc_genes)
   # return evenly spaced numbers over the number of nuclear genes range
   data_space = np.linspace(log10_Nuc_genes.min(), log10_Nuc_genes.max(), 200)
   # get probability values based on the pdf and the gene range
   evaluated = pdf.evaluate(data_space)
   # make a df
   pdf_df=pd.DataFrame({'data_space':data_space, 'evaluated':evaluated}) 
   return pdf_df

def get_cutoff(pdf_df, min_nucGene):
   '''function to return a cutoff'''
   # find the "dips" in the probability density - local minima
   local_minima = argrelextrema(pdf_df['evaluated'].values, np.less)
   # find the number of genes at which the minima occur
   potential_cutoffs = pdf_df['data_space'].values[local_minima]
   if len(np.where(potential_cutoffs > np.log10(min_nucGene))[0]) == 0:
      cutoff = np.log10(min_nucGene)
   # if there are any that are above 2, find the smallest one. because this is on a log10 scale, I transform back by 10^ (** in the script) to the value, and round to the nearest integer. 
   else: 
      cutoff = min(potential_cutoffs[np.where(potential_cutoffs>=np.log10(min_nucGene))])
   return cutoff

def make_pd_plots(pdf_df, min_nucGene, cutoff, sample):
   '''function to plot probability density with the default cutoffs and cell caller cutoff, and save the plot as a png'''
   plt.plot(pdf_df['data_space'], pdf_df['evaluated'])
   plt.axvline(cutoff, color = 'red',label = 'Nuclear genes'+'\n'+str(round(10 ** cutoff)))
   plt.axvline(np.log10(min_nucGene), color = 'black',label = 'Default'+'\n'+str(round(min_nucGene)))
   plt.title("Cell Caller Plot", fontdict = {'family':'sans-serif','color':'black','size':20,'fontweight':'bold'})
   plt.xlabel("log10(nNuc_genes+1)")
   plt.ylabel("Density")
   plt.legend(bbox_to_anchor = (1.0, 1), loc = 'upper right', title="Cells", fontsize=7)
   plt.savefig('{0}_pdf_with_cutoff.png'.format(sample), dpi=900)   # save the figure to file
   plt.close()    # close the figure window


def cell_caller(args):
   '''launch function'''
   sample = args.sample
   min_nucGene = float(args.min_nucGene)
   log10_Nuc_genes = getlog10NucGenes(sample)
   if len(set(log10_Nuc_genes)) == 1:
      cutoff = 2
   else:
      pdf_df = get_prob_dens_data(log10_Nuc_genes)
      cutoff = get_cutoff(pdf_df, min_nucGene)   
      make_pd_plots(pdf_df, min_nucGene, cutoff, sample)
   # because cutoff is on a log10 scale, I transform back by 10^ (** in the script) to the value, and round to the nearest integer
   print(round(10 ** cutoff), end="")

   # instead of return, write to txt file and don't forget tests

if __name__ == "__main__":
    args = parse_arguments(sys.argv[1:])
    # not sure what the best way to return this back to nextflow is, so leaving it like this for now (until I learn more about nextflow). 
    cell_caller(args)
