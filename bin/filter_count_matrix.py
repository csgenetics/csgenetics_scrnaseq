#!/usr/bin/env python

"""
Takes in an h5ad raw count matrix and a nuclear gene threshold to filter to
and outputs a filtered h5ad containing only those barcodes meeting the threshold.
Also outputs the filtered matrix in the tripartite format with:
    cell_only.barcodes.tsv.gz
    cell_only.features.tsv.gz
    cell_only.matrix.mtx.gz
"""

import sys
import gzip
from scipy.io import mmwrite
import os
import anndata
import numpy as np

class FilterCountMatrix:
    def __init__(self):
        # The number of nuclear genes that must have
        # been covered in order to consider a barcode a cell
        self.single_cell_nuc_gene_threshold = int(sys.argv[1])

        # Read in the h5ad matrix to an anndata object
        self.anndata_obj = anndata.read_h5ad(sys.argv[2])

        self.sample_name = sys.argv[3]

        # Compute the number of nuclear genes covered per barcode
        self.anndata_obj.obs['num_genes_covered_per_barcode'] = self.anndata_obj.X.toarray().astype(bool).sum(axis=1)
        self.anndata_obj.obs['num_mt_genes_covered_per_barcode'] = self.anndata_obj.X[:,self.anndata_obj.var_names.str.startswith('MT-')].toarray().astype(bool).sum(axis=1)
        self.anndata_obj.obs['num_nuc_genes_covered_per_barcode'] = self.anndata_obj.obs['num_genes_covered_per_barcode'] - self.anndata_obj.obs['num_mt_genes_covered_per_barcode']

        # Filter for only those barcodes that have >= to the threshold
        self.anndata_obj.obs["is_single_cell"] = np.where((self.anndata_obj.obs['num_nuc_genes_covered_per_barcode']>= self.single_cell_nuc_gene_threshold), True, False)
        self.anndata_obj_filtered = self.anndata_obj[self.anndata_obj.obs['is_single_cell'] == True]

        # Write out the filtered h5ad
        self.anndata_obj_filtered.write(f"{self.sample_name}.{self.single_cell_nuc_gene_threshold}.cell_only.count_matrix.h5ad")

        # Write out the barcodes, features, matrix files
        self.write_out_tripartite_filtered_matrix_files()

    def write_out_tripartite_filtered_matrix_files(self):
        """
        write out the
            cell_only.barcodes.tsv.gz
            cell_only.features.tsv.gz
            cell_only.matrix.mtx.gz
        """
        # Write filtered AnnData object into matrix file
        with gzip.open(f'{self.sample_name}.{self.single_cell_nuc_gene_threshold}.cell_only.matrix.mtx.gz', 'w') as mtx_file:
            mmwrite(mtx_file, a = self.anndata_obj_filtered.X.T, comment='', field='integer', precision=None, symmetry='general')

        # Write barcode table
        self.anndata_obj_filtered.obs.index.to_frame().to_csv(f'{self.sample_name}.{self.single_cell_nuc_gene_threshold}.cell_only.barcodes.tsv.gz', sep='\t', compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1}, header=None, index=False)

        # Write features table
        ft_table = self.anndata_obj_filtered.var.loc[:, ['ensID', 'geneSym']]
        ft_table['feature_type'] = 'Gene Expression'
        ft_table.to_csv(f'{self.sample_name}.{self.single_cell_nuc_gene_threshold}.cell_only.features.tsv.gz', sep="\t", compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1}, header=None, index=False)

if __name__ == "__main__":
    FilterCountMatrix()