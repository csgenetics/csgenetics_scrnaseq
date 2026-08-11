#!/usr/bin/env python
"""
Script to produce the _features_names.tsv used in the count_matrix process.

Takes in a .gtf file, and creates a tsv with two columns:
    ensID and geneSym

These are the gene ID and gene Name
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


VERSION_SUFFIX_PATTERN = r"\.[0-9]+$"


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

# Get in and out paths from 
gtf_in_path = sys.argv[1]
feature_names_out_path = sys.argv[2]

# Keep selected columns
gtf_obj = read_gtf(gtf_in_path, usecols=['gene_id','seqname','gene_name'])

# Keep gene-ID normalization in lockstep with tools/io_count_extract: remove only a
# terminal numeric version suffix (e.g. ".1"), preserving all other punctuation and
# Unicode. Reject collisions rather than allowing the count-table merge to duplicate
# counts for two distinct source identifiers that normalize to the same value.
if gtf_obj["gene_id"].isna().any() or (gtf_obj["gene_id"].str.len() == 0).any():
    raise ValueError("GTF contains an empty gene_id")

original_gene_ids = gtf_obj["gene_id"].copy()
normalized_gene_ids = original_gene_ids.str.replace(
    VERSION_SUFFIX_PATTERN, "", regex=True
)
if (normalized_gene_ids.str.len() == 0).any():
    raise ValueError("GTF gene_id is empty after version normalization")

gene_id_mapping = pd.DataFrame({
    "original": original_gene_ids,
    "normalized": normalized_gene_ids,
}).drop_duplicates()
collision_counts = gene_id_mapping.groupby("normalized")["original"].nunique()
colliding_ids = collision_counts[collision_counts > 1].index.tolist()
if colliding_ids:
    collision_details = []
    for normalized_id in colliding_ids:
        originals = sorted(
            gene_id_mapping.loc[
                gene_id_mapping["normalized"] == normalized_id, "original"
            ].tolist()
        )
        collision_details.append(
            f"{normalized_id!r} <- {', '.join(repr(value) for value in originals)}"
        )
    raise ValueError(
        "GTF gene IDs collide after version normalization: "
        + "; ".join(collision_details)
    )

gtf_obj["gene_id"] = normalized_gene_ids

feature_names_obj = gtf_obj.drop_duplicates()

# Split chromosome list and keep first value (they are all identical)
feature_names_obj['seqname'] = feature_names_obj.seqname.str.split(";", expand=True)[0]

feature_names_obj['gene_name'] = feature_names_obj.apply(lambda x: replace_na_gene_names_with_gene_id(x['gene_id'], x['gene_name']), axis=1)

# Select and rename columns
feature_names_obj = feature_names_obj.loc[:,['gene_id','gene_name','seqname']]
feature_names_obj = feature_names_obj.rename(columns={"seqname":"chromosome"})

# write output
feature_names_obj.to_csv(feature_names_out_path,
                    sep="\t",
                    quoting=None,
                    header=True,
                    index=False
                    )
