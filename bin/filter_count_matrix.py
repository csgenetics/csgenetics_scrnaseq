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
import re

class FilterCountMatrix:
    def __init__(self):
        # The number of nuclear genes that must have
        # been covered in order to consider a barcode a cell
        self.single_cell_nuc_gene_threshold = int(sys.argv[1])

        self.sample_name = sys.argv[3]

        # Read in the h5ad matrix to an anndata object
        # If the h5ad matrix is an empty file, simply write out
        # another empty file and exit
        try:
            self.anndata_obj = anndata.read_h5ad(sys.argv[2])
        except OSError:
            open(f"{self.sample_name}.{self.single_cell_nuc_gene_threshold}.filtered_feature_bc_matrix.empty.h5ad", "w").close()
            open(f"{self.sample_name}.{self.single_cell_nuc_gene_threshold}.raw_feature_bc_matrix.empty.h5ad", "w").close()
            sys.exit(0)

        # Whether we are working with a Human Mouse mixed samples
        # that has been mapped to a synthetic genome.
        if sys.argv[4] == "TRUE":
            self.mixed = True
            # We have a separate gene prefix string for each
            # of the species.
            self.hsap_mito_chr = sys.argv[5]
            self.mmus_mito_chr = sys.argv[6]

            # We have a separate gene prefix string for each
            # of the species.
            self.hsap_gene_prefix = sys.argv[7]
            self.mmus_gene_prefix = sys.argv[8]



            # The purity threshold used to classify a barcode as a cell
            # (in addition to the num nuc genes detected threshold)
            # if we are working with a mixed species sample
            self.purity = float(sys.argv[9])
        else:
            self.mixed = False
            self.mito_chr = sys.argv[5]

        if self.mixed:
            # Compute the number of nuclear genes covered per barcode
            # in order to be able to filter according to num_nuc_genes_covered_per_barcode

            self.anndata_obj.var['is_mito'] = np.where((self.anndata_obj.var['chromosome'] == self.hsap_mito_chr) | (self.anndata_obj.var['chromosome'] == self.mmus_mito_chr), True, False)
            self.anndata_obj.var['is_mito_hsap'] = np.where(self.anndata_obj.var['chromosome'] == self.hsap_mito_chr , True, False)
            self.anndata_obj.var['is_mito_mmus'] = np.where(self.anndata_obj.var['chromosome'] == self.mmus_mito_chr, True, False)
            self.anndata_obj.obs['num_mt_genes_covered_per_barcode'] = self.anndata_obj.X[:,self.anndata_obj.var['is_mito']].toarray().astype(bool).sum(axis=1)


            self.anndata_obj.obs['num_genes_covered_per_barcode'] = self.anndata_obj.X.toarray().astype(bool).sum(axis=1)
            self.anndata_obj.obs['num_nuc_genes_covered_per_barcode'] = self.anndata_obj.obs['num_genes_covered_per_barcode'] - self.anndata_obj.obs['num_mt_genes_covered_per_barcode']
           
            # For mixed species we additionally
            # need to filter according to a purity threshold.
            # The purity threshold is based on the number of gene detected by species.
            self.anndata_obj.obs['num_genes_detected_Hsap'] = self.anndata_obj.X[:,self.anndata_obj.var_names.str.startswith(self.hsap_gene_prefix)].toarray().astype(bool).sum(axis=1)
            self.anndata_obj.obs['num_genes_detected_Mmus'] = self.anndata_obj.X[:,self.anndata_obj.var_names.str.startswith(self.mmus_gene_prefix)].toarray().astype(bool).sum(axis=1)
            # To filter for the mitochondrial genes
            self.anndata_obj.obs['num_mito_genes_detected_Hsap'] = self.anndata_obj.X[:, self.anndata_obj.var['is_mito_hsap']].toarray().astype(bool).sum(axis=1)
            self.anndata_obj.obs['num_mito_genes_detected_Mmus'] = self.anndata_obj.X[:, self.anndata_obj.var['is_mito_mmus']].toarray().astype(bool).sum(axis=1)
            self.anndata_obj.obs['num_nuc_genes_detected_Hsap'] = self.anndata_obj.obs['num_genes_detected_Hsap'] - self.anndata_obj.obs['num_mito_genes_detected_Hsap']
            self.anndata_obj.obs['num_nuc_genes_detected_Mmus'] = self.anndata_obj.obs['num_genes_detected_Mmus'] - self.anndata_obj.obs['num_mito_genes_detected_Mmus']
            self.anndata_obj.obs['num_nuc_genes_detected_total'] = self.anndata_obj.obs['num_nuc_genes_detected_Hsap'] + self.anndata_obj.obs['num_nuc_genes_detected_Mmus']

            # Calculate purity metrics
            self.anndata_obj.obs['percent_nuc_genes_detected_Hsap'] = self.anndata_obj.obs['num_nuc_genes_detected_Hsap']/self.anndata_obj.obs['num_nuc_genes_detected_total']
            self.anndata_obj.obs['percent_nuc_genes_detected_Mmus'] = self.anndata_obj.obs['num_nuc_genes_detected_Mmus']/self.anndata_obj.obs['num_nuc_genes_detected_total']
            self.anndata_obj.obs['nuc_gene_purity'] = np.where((self.anndata_obj.obs['percent_nuc_genes_detected_Hsap']-self.anndata_obj.obs['percent_nuc_genes_detected_Mmus'])>0, self.anndata_obj.obs['percent_nuc_genes_detected_Hsap'], self.anndata_obj.obs['percent_nuc_genes_detected_Mmus'])
            self.anndata_obj.obs['species_based_on_nuc_gene_purity'] = np.where((self.anndata_obj.obs['percent_nuc_genes_detected_Hsap']-self.anndata_obj.obs['percent_nuc_genes_detected_Mmus'])>0, 'Hsap', 'Mmus')
        else:
            # Compute the number of nuclear genes covered per barcode
            self.anndata_obj.obs['num_genes_covered_per_barcode'] = self.anndata_obj.X.toarray().astype(bool).sum(axis=1)
            self.anndata_obj.var['is_mito'] = np.where(self.anndata_obj.var['chromosome'] == self.mito_chr, True, False)
            self.anndata_obj.obs['num_mt_genes_covered_per_barcode'] = self.anndata_obj.X[:,self.anndata_obj.var['is_mito']].toarray().astype(bool).sum(axis=1)
            self.anndata_obj.obs['num_nuc_genes_covered_per_barcode'] = self.anndata_obj.obs['num_genes_covered_per_barcode'] - self.anndata_obj.obs['num_mt_genes_covered_per_barcode']

        if self.mixed:
            # Filter for those barcodes that meet the num nuclear gene detected threshold AND
            # the purity threshold
            # We do this in two separate processes to collect the set of barcodes that meet
            # the num nuclear genes detected threshold independent of the purity (i.e. including
            # impure cells) as this designation will be used in internal metrics.
            self.anndata_obj.obs["is_called_cell"] = np.where(
                self.anndata_obj.obs['num_nuc_genes_covered_per_barcode']>= self.single_cell_nuc_gene_threshold, True, False
                )
            
            self.anndata_obj.obs["is_single_cell"] = np.where(
                self.anndata_obj.obs['is_called_cell'] & (self.anndata_obj.obs['nuc_gene_purity']>= self.purity), 
                True, False
                )
        else:
            # Filter for only those barcodes that have >= to the num nuclear gene detected threshold
            self.anndata_obj.obs["is_single_cell"] = np.where(
                self.anndata_obj.obs['num_nuc_genes_covered_per_barcode']>= self.single_cell_nuc_gene_threshold,
                True, False
                )
            
        # Write out the raw count matrix
        self.anndata_obj.write(f"{self.sample_name}.{self.single_cell_nuc_gene_threshold}.raw_feature_bc_matrix.h5ad")

        # Filter down to single cells and write out
        self.anndata_obj_filtered = self.anndata_obj[self.anndata_obj.obs['is_single_cell'] == True]
        self.anndata_obj_filtered.write(f"{self.sample_name}.{self.single_cell_nuc_gene_threshold}.filtered_feature_bc_matrix.h5ad")

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
        with gzip.open('matrix.mtx.gz', 'w') as mtx_file:
            mmwrite(mtx_file, a = self.anndata_obj_filtered.X.T, comment='', field='integer', precision=None, symmetry='general')

        # Write barcode table
        self.anndata_obj_filtered.obs.index.to_frame().to_csv('barcodes.tsv.gz', sep='\t', compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1}, header=None, index=False)

        # Write features table
        ft_table = self.anndata_obj_filtered.var.loc[:, ['gene_id', 'gene_name']]
        ft_table['feature_type'] = 'Gene Expression'
        ft_table.to_csv('features.tsv.gz', sep="\t", compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1}, header=None, index=False)

if __name__ == "__main__":
    FilterCountMatrix()