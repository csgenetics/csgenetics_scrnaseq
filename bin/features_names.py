#!/usr/bin/env python
"""
Script to produce the _features_names.tsv used in the count_matrix process.

Takes in a .gtf file, filters down to XXX and XXX.
"""

import math
import sys
from gtfparse import read_gtf
# surpress the FutureWarning that is being output by the read_gtf function.
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
# supress the SettingWithCopyWarning warning
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

def isnan(value):
    try:
        return math.isnan(float(value))
    except:
        return False
    
def replace_na_gene_names_with_gene_id(gene_id, gene_name):
    """
    If gene_name is NA then replace with gene_id
    """
    if isnan(gene_name) or not gene_name:
        gene_name2 = gene_id
    else:
        gene_name2 = gene_name
    return(gene_name2)

def annotate_mito_genes_with_mt_prefix(chromosome, gene_name):
    """
    Annotate the MT Genes with a "MT_" prefix
    """
    if chromosome == "MT" and gene_name != "InterGenic":
        gene_name2 = chromosome + "_" + gene_name
    else:
        # do not annotate if not mitochondrial
        gene_name2 = gene_name
    return gene_name2

# Get in and out paths from 
gtf_in_path = sys.argv[1]
feature_names_out_path = sys.argv[2]

gtf_obj = read_gtf(gtf_in_path, usecols=['gene_id','seqname','gene_name'])

# Keep selected columns
feature_names_obj = gtf_obj.drop_duplicates()

# Split chromosome list and keep first value (they are all identical)
feature_names_obj['seqname'] = feature_names_obj.seqname.str.split(";", expand=True)[0]

feature_names_obj['gene_name1'] = feature_names_obj.apply(lambda x: replace_na_gene_names_with_gene_id(x['gene_id'], x['gene_name']), axis=1)
feature_names_obj['gene_name'] = feature_names_obj.apply(lambda x: annotate_mito_genes_with_mt_prefix(x['seqname'], x['gene_name1']), axis=1)

# Select and rename columns
feature_names_obj = feature_names_obj.loc[:,['gene_id','gene_name']]
feature_names_obj = feature_names_obj.rename(columns={"gene_id":"ensID", "gene_name":"geneSym"})

# write output
feature_names_obj.to_csv(feature_names_out_path,
                    sep="\t",
                    quoting=None,
                    header=True,
                    index=False
                    )
