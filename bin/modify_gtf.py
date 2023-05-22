#!/usr/bin/env python
import pandas as pd
import numpy as np
import math
import sys
from pathlib import Path
from gtfparse import read_gtf

def isnan(value):
    try:
        return math.isnan(float(value))
    except:
        return False

file1 = sys.argv[1]
file2 = Path(file1).stem

output_name = sys.argv[2]


dat = read_gtf(file1)

# Keep selected columns
dat_select = dat.loc[:,['gene_id','seqname','gene_name']].drop_duplicates()
# Split chromosome list and keep first value (they are all identical)
dat_select['seqname'] = dat_select.seqname.str.split(";",expand=True)[0]

# If gene_name is NA then replace with gene_id
def gene_NA(gene_id, gene_name):
    if isnan(gene_name) or not gene_name:
        gene_name2 = gene_id
    else:
        gene_name2 = gene_name
    return(gene_name2)

# Annotate MT Genes with a MT prefix
def label_MT(chromosome, gene_name):
    if chromosome == "MT" and gene_name != "InterGenic":
        gene_name2 = chromosome + "_" + gene_name
    else:
        # do not annotate if not mitochondrial
        gene_name2 = gene_name
    return gene_name2

dat_select['gene_name1'] = dat_select.apply(lambda x: gene_NA(x['gene_id'], x['gene_name']), axis=1)
dat_select['gene_name'] = dat_select.apply(lambda x: label_MT(x['seqname'], x['gene_name1']), axis=1)


# Select and rename columns
dat_final = dat_select.loc[:,['gene_id','gene_name']]
dat_final = dat_final.rename(columns={"gene_id":"ensID", "gene_name":"geneSym"})
# write output
dat_final.to_csv(output_name,
                    sep="\t",
                    quoting=None,
                    header=True,
                    index=False
                    )
