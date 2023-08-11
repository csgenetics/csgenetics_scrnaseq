#!/usr/bin/env python

from pandas import read_csv, read_table, merge, DataFrame, errors
from pandas.api.types import CategoricalDtype
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
from numpy import zeros,hstack
from anndata import AnnData
import sys, argparse, gzip

# Parse the command line arguments
def parse_arguments(args):
    parser = argparse.ArgumentParser(description = "Arguments for int_metrics script to calculate internal metrics")
    parser.add_argument("--white_list", help="path to white_list.csv file")
    parser.add_argument("--count_table", help="path to count table (bcGeneSummary.txt) file ")
    parser.add_argument("--gene_list", help="path to features.tsv.gz file")
    parser.add_argument("--sample", help="Sample ID")

    return parser.parse_args()

def make_count_matrix(args):

    # load CSGX whitelist of IOs
    wl = read_csv(args.white_list)
    wl.columns = ['cell', 'io','ioID']
    wl = wl.iloc[:,[0,1]]

    # load count table obtained using UMI-tools
    try:
        counts = read_table(args.count_table,header=None)
    except errors.EmptyDataError:
        # The input file is empty and we simply write out an empty .h5ad
        # to be picked up by the process then exit
        open(f"{args.sample}.raw_feature_bc_matrix.empty.h5ad", "w").close()
        sys.exit(0)

    counts.columns = ['io', 'gene_id']

     # load feature (gene) names obtained from the genome GTF file
    genes_all_unfiltered = read_table(args.gene_list).drop_duplicates()
    genes_all = genes_all_unfiltered[genes_all_unfiltered['gene_id'].isin(counts['gene_id'])]    

    # first merge by IOs
    counts = merge(wl, counts)

    # sum together counts corresponding to the same cell-gene_name pair
    counts = counts.groupby(['cell', 'gene_id']).size().reset_index(name='count')

    # second merge by gene ID
    counts = merge(counts, genes_all,how='left')

    zero_genes = genes_all.loc[-genes_all.gene_id.isin(counts.gene_id)]

    # transform a long counts DataFrame into a sparse matrix
    cell_c = CategoricalDtype(sorted(counts.cell.unique()), ordered=True)
    name_c = CategoricalDtype(sorted(counts.gene_id.unique()), ordered=True)
    row = counts.gene_id.astype(name_c).cat.codes
    col = counts.cell.astype(cell_c).cat.codes

    #sparse_matrix = csr_matrix((counts["count"], (col, row)), shape=(cell_c.categories.size, name_c.categories.size))
    non_zero_matrix = csr_matrix((counts["count"], (col,row)), shape=(cell_c.categories.size,name_c.categories.size)).toarray()
    # Create a zero matrix for all zero-genes
    zero_matrix = zeros((non_zero_matrix.shape[0],zero_genes.shape[0]))
    all_matrix = hstack([non_zero_matrix, zero_matrix])
    # sparse_matrix for annData needs to be in the format of n_obs x n_vars
    sparse_matrix = csr_matrix(all_matrix)

    # create an AnnData object for downstream analysis
    ft_names= DataFrame(name_c.categories.tolist() + zero_genes.gene_id.tolist(),columns=['gene_id'])
    ft_names = merge(genes_all,ft_names,how='right')

    adata = AnnData(sparse_matrix,var=ft_names)
    adata.var_names = ft_names.gene_name.tolist()
    adata.var_names_make_unique()
    adata.obs_names = cell_c.categories


    # write AnnData object into H5 file
    adata.write(f"{args.sample}.raw_feature_bc_matrix.h5ad")

    # Write sparse matrix format
    mtx_file = gzip.open('matrix.mtx.gz', 'w')
    mmwrite(mtx_file, a = sparse_matrix.T, comment='', field='integer', precision=None, symmetry='general')
    # Write feature table
    ft_names['feature_type'] = 'Gene Expression'
    ft_names.to_csv('features.tsv.gz', sep="\t",compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1},header=None, index=False)
    # Write barcode table
    adata.obs.index.to_frame().to_csv('barcodes.tsv.gz',sep='\t',compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1},header=None, index=False)


if __name__ == "__main__":
    args = parse_arguments(sys.argv[1:])
    make_count_matrix(args)
