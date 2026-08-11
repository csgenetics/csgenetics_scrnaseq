#!/usr/bin/env python3

"""
Categorize reads and counts from STAR BAM as in cells or out of cells.
Reads cell information from raw count matrix H5AD file with cell annotations.
Calculates both read-based and count-based metrics.

Customer pipeline version - handles single cell caller method only.
"""

import argparse
import json
from numbers import Integral, Real
from pathlib import Path
import re

import anndata as ad
import numpy as np
import pandas as pd
import pysam
from scipy import sparse

from empty_h5ad import is_empty_h5ad_sentinel


MAX_SIGNED_64_BIT = (1 << 63) - 1


def as_nonnegative_integer(value, label):
    """Validate and return one finite, non-negative, signed-64-bit integer."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{label} must be a numeric integer, got {value!r}")

    if isinstance(value, Integral):
        integer = int(value)
    else:
        numeric = float(value)
        if not np.isfinite(numeric) or not numeric.is_integer():
            raise ValueError(f"{label} must be a finite integer, got {value!r}")
        integer = int(numeric)

    if integer < 0:
        raise ValueError(f"{label} must be non-negative, got {integer}")
    if integer > MAX_SIGNED_64_BIT:
        raise ValueError(f"{label} exceeds the signed 64-bit limit: {integer}")
    return integer

def extract_barcode_from_read_name(read_name, barcode_length):
    """
    Extract barcode from read name. Barcode is surrounded by underscores.
    Example: VH01429:149:22253TGNX:1:2601:69579:24995_TGAGCCACATCGA_

    Args:
        read_name: The read name from the BAM file
        barcode_length: Expected length of the barcode sequence
    """
    # Look for barcode pattern - sequence of ACGTN surrounded by underscores
    pattern = rf'_([ACGTN]{{{barcode_length}}})_'
    match = re.search(pattern, read_name)
    if match:
        return match.group(1)
    return None

def get_raw_reads_from_fastp(fastp_json_path):
    """Extract total reads from fastp JSON file."""
    with open(fastp_json_path, 'r') as f:
        fastp_data = json.load(f)

    # fastp stores total reads in summary.before_filtering.total_reads
    return as_nonnegative_integer(
        fastp_data['summary']['before_filtering']['total_reads'],
        "fastp summary.before_filtering.total_reads",
    )


def load_count_matrix(h5ad_path):
    """Load an H5AD, accepting only the pipeline's named empty sentinel."""
    path = Path(h5ad_path)
    if is_empty_h5ad_sentinel(path):
        return None
    return ad.read_h5ad(path)


def _validate_count_matrix_values(adata):
    """Reject values that cannot represent raw non-negative integer counts."""
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    values = np.asarray(values)
    if not np.issubdtype(values.dtype, np.number) or not np.isrealobj(values):
        raise ValueError("H5AD count matrix must contain real numeric values")
    if not np.all(np.isfinite(values)):
        raise ValueError("H5AD count matrix contains non-finite values")
    if np.any(values < 0):
        raise ValueError("H5AD count matrix contains negative values")
    if np.any(values != np.floor(values)):
        raise ValueError("H5AD count matrix contains non-integer values")
    if np.any(values > MAX_SIGNED_64_BIT):
        raise ValueError("H5AD count matrix contains a value above the signed 64-bit limit")

def get_cells_from_h5ad(adata, mixed_species=False):
    """
    Extract cell barcodes from raw count matrix H5AD file.
    For customer pipeline, cells are marked with 'is_single_cell' column.

    Returns dict with cell sets for total and optionally by species.
    """
    cells = {}

    if mixed_species:
        # For mixed species, we have is_single_cell, is_hsap_cell, is_mmus_cell
        if 'is_single_cell' not in adata.obs.columns:
            raise ValueError("H5AD file does not have 'is_single_cell' column in obs")

        cell_mask = adata.obs['is_single_cell'].astype(bool)
        cells['total'] = set(adata.obs_names[cell_mask])

        if 'is_hsap_cell' in adata.obs.columns:
            hsap_mask = adata.obs['is_hsap_cell'].astype(bool)
            cells['hsap'] = set(adata.obs_names[hsap_mask])

        if 'is_mmus_cell' in adata.obs.columns:
            mmus_mask = adata.obs['is_mmus_cell'].astype(bool)
            cells['mmus'] = set(adata.obs_names[mmus_mask])
    else:
        # For single species, just use is_single_cell
        if 'is_single_cell' not in adata.obs.columns:
            raise ValueError("H5AD file does not have 'is_single_cell' column in obs")

        cell_mask = adata.obs['is_single_cell'].astype(bool)
        cells['total'] = set(adata.obs_names[cell_mask])

    return cells

def get_counts_by_cell_status(adata, sample_id):
    """
    Calculate counts for cells vs non-cells using dataframe operations.

    Args:
        adata: AnnData object with count matrix
        sample_id: Sample identifier
    """
    # Check if is_single_cell column exists
    if 'is_single_cell' not in adata.obs.columns:
        raise ValueError("H5AD file does not have 'is_single_cell' column in obs")

    _validate_count_matrix_values(adata)

    # Create boolean mask for cells
    cell_mask = adata.obs['is_single_cell'].astype(bool).to_numpy()

    # Get the count matrix as a dataframe
    # Sum counts across all genes for each barcode
    counts_per_barcode = np.array(adata.X.sum(axis=1)).flatten()

    # Calculate counts in cells vs out of cells using boolean indexing
    counts_in_cells = as_nonnegative_integer(
        counts_per_barcode[cell_mask].sum(), "counts in cells"
    )
    counts_out_of_cells = as_nonnegative_integer(
        counts_per_barcode[~cell_mask].sum(), "counts out of cells"
    )

    return counts_in_cells, counts_out_of_cells


def categorize_primary_reads(bam_path, cell_barcodes, barcode_length):
    """Count each primary BAM record once, excluding secondary/supplementary records."""
    reads_in_cells = 0
    reads_out_of_cells = 0
    total_reads_processed = 0
    skipped_secondary = 0
    skipped_supplementary = 0

    # The pipeline's intentional empty BAM carries @HD but no @SQ records.
    # check_sq=False accepts that valid sentinel while pysam still validates the
    # BAM container itself; ordinary STAR BAMs retain their normal sequence
    # dictionary and behavior.
    with pysam.AlignmentFile(bam_path, 'rb', check_sq=False) as bam:
        for read in bam:
            if read.is_secondary:
                skipped_secondary += 1
                continue
            if read.is_supplementary:
                skipped_supplementary += 1
                continue

            total_reads_processed += 1
            barcode = extract_barcode_from_read_name(read.query_name, barcode_length)
            if barcode is None:
                print(f"Warning: Failed to extract barcode from read: {read.query_name}")
                continue

            if barcode in cell_barcodes:
                reads_in_cells += 1
            else:
                reads_out_of_cells += 1

            if total_reads_processed % 1000000 == 0:
                print(f"Processed {total_reads_processed:,} reads...")

    return {
        "reads_in_cells": reads_in_cells,
        "reads_out_of_cells": reads_out_of_cells,
        "total_reads_processed": total_reads_processed,
        "skipped_secondary": skipped_secondary,
        "skipped_supplementary": skipped_supplementary,
    }

def main():
    parser = argparse.ArgumentParser(description='Categorize reads and counts by cell status')
    parser.add_argument('--sample_id', required=True, help='Sample identifier')
    parser.add_argument('--star_bam', required=True, help='STAR output BAM file')
    parser.add_argument('--raw_count_matrix_h5ad', required=True, help='Raw count matrix H5AD with cell annotations')
    parser.add_argument('--fastp_json', required=True, help='Fastp JSON file with read counts')
    parser.add_argument('--barcode_length', type=int, required=True, help='Expected barcode length')
    parser.add_argument('--mixed_species', action='store_true', help='Flag for mixed species samples')

    args = parser.parse_args()

    print(f"Processing {args.sample_id}")
    print(f"Mixed species: {args.mixed_species}")
    print(f"Expected barcode length: {args.barcode_length}")

    # Get raw reads count from fastp JSON
    raw_reads = get_raw_reads_from_fastp(args.fastp_json)
    print(f"Raw reads from fastp: {raw_reads:,}")

    # Load the raw count matrix, or recognize the pipeline's explicit
    # zero-byte *.empty.h5ad sentinel. Generic zero-byte/corrupt H5ADs are not
    # silently treated as empty data.
    adata = load_count_matrix(args.raw_count_matrix_h5ad)
    if adata is None:
        print("Loaded named empty H5AD sentinel")
        cell_barcodes = set()
        counts_in_cells = 0
        counts_out_of_cells = 0
    else:
        print(f"Loaded H5AD with {adata.n_obs} barcodes and {adata.n_vars} genes")
        cells_dict = get_cells_from_h5ad(adata, args.mixed_species)
        cell_barcodes = cells_dict['total']
        counts_in_cells, counts_out_of_cells = get_counts_by_cell_status(
            adata, args.sample_id
        )

    print(f"Found {len(cell_barcodes)} cells in H5AD")
    print(f"Counts in cells: {counts_in_cells:,}")
    print(f"Counts out of cells: {counts_out_of_cells:,}")

    # Strip sample prefix from cell barcodes for faster BAM processing
    prefix = f"{args.sample_id}_"
    cell_barcodes_stripped = {
        bc[len(prefix):] for bc in cell_barcodes if bc.startswith(prefix)
    }

    # Parse STAR BAM for read categorization
    print(f"Processing BAM file: {args.star_bam}")
    bam_metrics = categorize_primary_reads(
        args.star_bam, cell_barcodes_stripped, args.barcode_length
    )
    reads_in_cells = bam_metrics["reads_in_cells"]
    reads_out_of_cells = bam_metrics["reads_out_of_cells"]

    print(f"Total reads processed from BAM: {bam_metrics['total_reads_processed']:,}")
    print(f"Skipped secondary alignments: {bam_metrics['skipped_secondary']:,}")
    print(f"Skipped supplementary alignments: {bam_metrics['skipped_supplementary']:,}")
    print(f"Reads in cells: {reads_in_cells:,}")
    print(f"Reads out of cells: {reads_out_of_cells:,}")

    # Calculate unusable reads
    total_mapped = reads_in_cells + reads_out_of_cells
    if total_mapped > raw_reads:
        raise ValueError(
            "Primary alignments with extractable barcodes exceed raw reads: "
            f"mapped={total_mapped}, raw={raw_reads}"
        )
    if bam_metrics["total_reads_processed"] > raw_reads:
        raise ValueError(
            "Primary BAM records exceed raw reads: "
            f"primary={bam_metrics['total_reads_processed']}, raw={raw_reads}"
        )
    unusable_reads = raw_reads - total_mapped

    print(f"Unusable reads: {unusable_reads:,}")

    # Create output dataframe with all metrics
    metrics = {
        'sample_id': args.sample_id,
        'raw_reads': raw_reads,
        'reads_in_cells': reads_in_cells,
        'reads_out_of_cells': reads_out_of_cells,
        'unusable_reads': unusable_reads,
        'counts_in_cells': counts_in_cells,
        'counts_out_of_cells': counts_out_of_cells
    }

    # Add species-specific metrics if mixed
    if args.mixed_species:
        # For mixed species, we could add species-specific metrics here if needed
        # This would require additional BAM parsing with gene assignments
        pass

    # Write CSV output
    df = pd.DataFrame([metrics])
    output_file = f'{args.sample_id}.read_categorization.csv'
    df.to_csv(output_file, index=False)

    print(f"\nWrote metrics to {output_file}")

if __name__ == '__main__':
    main()
