#!/usr/bin/env python

from pandas import read_csv, read_table, merge, DataFrame, errors
from pandas.api.types import CategoricalDtype
from scipy.sparse import csr_matrix, hstack
from scipy.io import mmwrite
from numpy import zeros, float32, where
from anndata import AnnData
import sys, argparse, gzip

class CountMatrix:
    def __init__(self):
        self.parse_args()
    # Parse the command line arguments
    def parse_args(self):
        parser = argparse.ArgumentParser(description = "Arguments for int_metrics script to calculate internal metrics")
        parser.add_argument("--barcode_list", help="path to allow_list.csv file")
        parser.add_argument("--count_table", help="path to count table (bcGeneSummary.txt) file ")
        parser.add_argument("--gene_list", help="path to features.tsv.gz file")
        parser.add_argument("--sample", help="Sample ID")
        parser.add_argument("--mixed_species", help="Flag to indicate if the sample is mixed species", required=True) # String 'True' or 'False'
        parser.add_argument("--mito_chr", help="Mitochondrial chromosome symbol", required=False) # For single species
        parser.add_argument("--hsap_mito_chr", help="Human mitochondrial chromosome symbol", required=False) # For mixed species only
        parser.add_argument("--mmus_mito_chr", help="Mouse mitochondrial chromosome symbol", required=False) # For mixed species only
        parser.add_argument("--hsap_gene_prefix", help="Human gene prefix", required=False) # For mixed species only
        parser.add_argument("--mmus_gene_prefix", help="Mouse gene prefix", required=False) # For mixed species only

        args_dict = vars(parser.parse_args())
        for key in args_dict:
            setattr(self, key, args_dict[key])

        if self.mixed_species.lower() == "true":
            self.single_species = False
            self.mixed_species = True
        else:
            self.single_species = True
            self.mixed_species = False

    def make_count_matrix(self):
        # Makes self.anndata_obj
        self.make_base_count_matrix()

        # Adds gene annotations and others
        self.annotate_count_matrix()
    
        # Write out the matrix
        self.write_count_matrix()
    
    def annotate_count_matrix(self):
        if self.mixed_species:
            # Annotate the mitochondrial genes to report the mito and nuclear genes in the summary statistics
            self.anndata_obj.var['is_mito'] = where((self.anndata_obj.var['chromosome'] == self.hsap_mito_chr) | (self.anndata_obj.var['chromosome'] == self.mmus_mito_chr), True, False)
            self.anndata_obj.var['is_mito_hsap'] = where(self.anndata_obj.var['chromosome'] == self.hsap_mito_chr , True, False)
            self.anndata_obj.var['is_mito_mmus'] = where(self.anndata_obj.var['chromosome'] == self.mmus_mito_chr, True, False)
           
            # Annotate the Hsap and Mmus genes to be able to report Hsap and Mmus specific stats in summary_statistics
            self.anndata_obj.var['is_hsap'] = where(self.anndata_obj.var_names.str.startswith(self.hsap_gene_prefix), True, False)
            self.anndata_obj.var['is_mmus'] = where(self.anndata_obj.var_names.str.startswith(self.mmus_gene_prefix), True, False)

            self.anndata_obj.obs["hsap_counts"] = self.anndata_obj.X[:, self.anndata_obj.var["is_hsap"]].toarray().sum(axis=1)
            self.anndata_obj.obs["mmus_counts"] = self.anndata_obj.X[:, self.anndata_obj.var["is_mmus"]].toarray().sum(axis=1)
            self.anndata_obj.obs["total_counts"] = self.anndata_obj.X.toarray().sum(axis=1)
            
        else:
            # Annotate mitochondrial genes to report the mito and nuclear genes in the summary statistics
            self.anndata_obj.var['is_mito'] = where(self.anndata_obj.var['chromosome'] == self.mito_chr, True, False)
            self.anndata_obj.obs['total_counts'] = self.anndata_obj.X.toarray().sum(axis=1)

    def make_base_count_matrix(self):
        # Load CS Genetics barcode_list of IOs
        bcl = read_csv(self.barcode_list)
        bcl.columns = ['cell', 'io','ioID']
        bcl = bcl.iloc[:,[0,1]]

        # Load count table from io_count
        try:
            counts = read_table(self.count_table, header=None)
        except errors.EmptyDataError:
            # The input file is empty and we simply write out an empty .h5ad
            # to be picked up by the process then exit
            open(f"{self.sample}.raw_feature_bc_matrix.empty.h5ad", "w").close()
            sys.exit(0)

        counts.columns = ['io', 'gene_id']

        # Load feature (gene) names from the genome GTF file
        genes_all = read_table(self.gene_list).drop_duplicates()

        # First merge by IOs
        counts = merge(bcl, counts)

        # Sum together counts corresponding to the same cell-gene_name pair
        counts = counts.groupby(['cell', 'gene_id']).size().reset_index(name='count')

        # Second merge by gene ID
        counts = merge(counts, genes_all, how='left')

        zero_genes = genes_all.loc[-genes_all.gene_id.isin(counts.gene_id)]

        # Transform a long counts DataFrame into a sparse matrix
        cell_c = CategoricalDtype(sorted(counts.cell.unique()), ordered=True)
        name_c = CategoricalDtype(sorted(counts.gene_id.unique()), ordered=True)
        row = counts.gene_id.astype(name_c).cat.codes
        col = counts.cell.astype(cell_c).cat.codes

        non_zero_matrix = csr_matrix((counts["count"], (col,row)), shape=(cell_c.categories.size,name_c.categories.size))
        # Create a zero matrix for all zero-genes
        zero_matrix = csr_matrix((non_zero_matrix.shape[0], zero_genes.shape[0]))

        # Sparse_matrix for annData needs to be in the format of n_obs x n_vars
        self.sparse_matrix = hstack([non_zero_matrix, zero_matrix])

        # Create an AnnData object for downstream analysis
        self.ft_names= DataFrame(name_c.categories.tolist() + zero_genes.gene_id.tolist(),columns=['gene_id'])
        self.ft_names = merge(genes_all, self.ft_names,how='right')

        # Have to specify float32 so that it is compatible with BPCells package in R for Seurat v5
        # and anndata package in Seurat v4.
        # Convert sparse matrix to float32, and create new AnnData object
        # If you try to convert adata.X directly on anndata with only 1 row, it will throw an error
        sparse_matrix_float32 = self.sparse_matrix.astype(float32)
        self.anndata_obj = AnnData(sparse_matrix_float32,var=self.ft_names)
        self.anndata_obj.var_names = self.ft_names.gene_name.tolist()
        self.anndata_obj.var_names_make_unique()
        self.anndata_obj.obs_names = cell_c.categories

        # It is important to add the sample name to make the barcode names unique
        # for use in Seurat v5.
        self.anndata_obj.obs_names = [self.sample + "_" + _ for _ in self.anndata_obj.obs_names]

    def write_count_matrix(self):
        self.anndata_obj.write(f"{self.sample}.raw_feature_bc_matrix.h5ad")

        # Write sparse matrix format
        mtx_file = gzip.open('matrix.mtx.gz', 'w')
        mmwrite(mtx_file, a = self.sparse_matrix.T, comment='', field='integer', precision=None, symmetry='general')
        # Write feature table
        self.ft_names['feature_type'] = 'Gene Expression'
        # Drop the chromosome column 
        self.ft_names.drop('chromosome', axis=1, inplace=True)

        self.ft_names.to_csv('features.tsv.gz', sep="\t",compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1},header=None, index=False)
        
        # Write barcode table
        self.anndata_obj.obs.index.to_frame().to_csv('barcodes.tsv.gz',sep='\t',compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1},header=None, index=False)


if __name__ == "__main__":
    CountMatrix().make_count_matrix()
