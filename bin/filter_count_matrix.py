#!/usr/bin/env python

"""
Takes in an h5ad raw count matrix and a nuclear gene threshold to filter to
and outputs a filtered h5ad containing only those barcodes meeting the threshold.
Also outputs the filtered matrix in the tripartite format with:
    barcodes.tsv.gz
    features.tsv.gz
    matrix.mtx.gz
"""

import sys
import gzip
from scipy.io import mmwrite
import os
import anndata
import numpy as np
import re

from empty_h5ad import is_empty_h5ad_sentinel

class FilterCountMatrix:
    def __init__(self):
        # Whether we are working with a Human Mouse mixed samples
        # that has been mapped to a synthetic genome.
        if sys.argv[4] == "TRUE":
            self.mixed = True
            # The counts required to consider a barcode a cell
            # hsap_counts >= hsap_count_threshold & mmus_counts <= mmus_count_threshold is a single-cell (human).
            # mmus_counts >= mmus_count_threshold & hsap_counts <= hsap_count_threshold is a single-cell (mouse).
            # mmus_counts < mmus_count_threshold & hsap_counts < hsap_count_threshold is a noisy barcode.
            # mmus_counts > mmus_count_threshold & hsap_counts > hsap_count_threshold is a multiplet.
            self.hsap_count_threshold = int(sys.argv[1].split("_")[0])
            self.mmus_count_threshold = int(sys.argv[1].split("_")[1])

        else:
            self.mixed = False
            # The counts required to consider a barcode a cell
            self.single_cell_count_threshold = int(sys.argv[1])
       
        self.sample_name = sys.argv[3]

        # Only the pipeline's explicit, zero-byte *.empty.h5ad sentinel may
        # produce downstream empty sentinels. Corrupt or generic zero-byte
        # inputs must fail rather than being relabelled as intentional emptiness.
        if is_empty_h5ad_sentinel(sys.argv[2]):
            open(f"{self.sample_name}.{sys.argv[1]}.filtered_feature_bc_matrix.empty.h5ad", "w").close()
            open(f"{self.sample_name}.{sys.argv[1]}.raw_feature_bc_matrix.empty.h5ad", "w").close()
            sys.exit(0)
        self.anndata_obj = anndata.read_h5ad(sys.argv[2])

        if self.mixed:
            # For mixed we specifically annotate Hsap and Mmus cells as well as the more generic is_called_cell and is_single_cell
            self.anndata_obj.obs["is_called_cell"] = ((self.anndata_obj.obs['hsap_counts'] >= self.hsap_count_threshold) | (self.anndata_obj.obs['mmus_counts'] >= self.mmus_count_threshold))
            self.anndata_obj.obs["is_single_cell"] = (
                ((self.anndata_obj.obs['hsap_counts'] >= self.hsap_count_threshold) & (self.anndata_obj.obs['mmus_counts'] <= self.mmus_count_threshold)) | 
                ((self.anndata_obj.obs['hsap_counts'] <= self.hsap_count_threshold) & (self.anndata_obj.obs['mmus_counts'] >= self.mmus_count_threshold))
                )
            self.anndata_obj.obs["is_hsap_cell"] = (
                (self.anndata_obj.obs['hsap_counts'] >= self.hsap_count_threshold) & (self.anndata_obj.obs['mmus_counts'] <= self.mmus_count_threshold)
                )
            # To prevent cases where the Hsap and Mmus counts both exactly meet the thresholds leading to either a cell
            # classifed as both Hsap and Mmus or as neither, we will classify mmus_cell as is_single_cell that are not hsap_cell
            self.anndata_obj.obs["is_mmus_cell"] = (
                (self.anndata_obj.obs['is_single_cell']) & (self.anndata_obj.obs['is_hsap_cell'] == False)
                )
            
            assert self.anndata_obj.obs["is_hsap_cell"].sum() + self.anndata_obj.obs["is_mmus_cell"].sum() == self.anndata_obj.obs["is_single_cell"].sum()
            
        else:
            # Filter for only those barcodes that have >= to the count threshold
            self.anndata_obj.obs["is_single_cell"] = np.where(
                self.anndata_obj.obs['total_counts']>= self.single_cell_count_threshold,
                True, False
                )

        # Write out the raw count matrix
        self.anndata_obj.write(f"{self.sample_name}.{sys.argv[1]}.raw_feature_bc_matrix.h5ad", compression="gzip", compression_opts=9)

        # Filter down to single cells and write out
        self.anndata_obj_filtered = self.anndata_obj[self.anndata_obj.obs['is_single_cell'] == True]
        self.anndata_obj_filtered.write(f"{self.sample_name}.{sys.argv[1]}.filtered_feature_bc_matrix.h5ad", compression="gzip", compression_opts=9)

        # Write out the barcodes, features, matrix files
        self.write_out_tripartite_filtered_matrix_files()

    def write_out_tripartite_filtered_matrix_files(self):
        """
        write out the
            barcodes.tsv.gz
            features.tsv.gz
            matrix.mtx.gz
        """
        # Write filtered AnnData object into matrix file
        # Use a fixed gzip header timestamp so identical filtered matrices are
        # byte-identical across runs.  The filename and decompressed payload stay
        # unchanged.
        with gzip.GzipFile(filename='matrix.mtx.gz', mode='wb', mtime=0) as mtx_file:
            mmwrite(mtx_file, a = self.anndata_obj_filtered.X.T, comment='', field='integer', precision=None, symmetry='general')

        # Write barcode table
        self.anndata_obj_filtered.obs.index.to_frame().to_csv('barcodes.tsv.gz', sep='\t', compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1}, header=None, index=False)

        # Write features table
        ft_table = self.anndata_obj_filtered.var.loc[:, ['gene_id', 'gene_name']]
        ft_table['feature_type'] = 'Gene Expression'
        ft_table.to_csv('features.tsv.gz', sep="\t", compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1}, header=None, index=False)

if __name__ == "__main__":
    FilterCountMatrix()
