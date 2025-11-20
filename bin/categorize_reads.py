#!/usr/bin/env python3

"""
Categorize reads and counts from STAR BAM as in cells or out of cells.
Reads cell information from raw count matrix H5AD file with cell annotations.
Calculates both read-based and count-based metrics.

Customer pipeline version - handles single cell caller method only.
"""

import pysam
import anndata as ad
import pandas as pd
import numpy as np
import argparse
import json
import re
import sys

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
    return fastp_data['summary']['before_filtering']['total_reads']

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

    # Create boolean mask for cells
    cell_mask = adata.obs['is_single_cell'].astype(bool)

    # Get the count matrix as a dataframe
    # Sum counts across all genes for each barcode
    counts_per_barcode = np.array(adata.X.sum(axis=1)).flatten()

    # Calculate counts in cells vs out of cells using boolean indexing
    counts_in_cells = counts_per_barcode[cell_mask].sum()
    counts_out_of_cells = counts_per_barcode[~cell_mask].sum()

    return int(counts_in_cells), int(counts_out_of_cells)

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

    # Load raw count matrix H5AD which has cell annotations
    adata = ad.read_h5ad(args.raw_count_matrix_h5ad)
    print(f"Loaded H5AD with {adata.n_obs} barcodes and {adata.n_vars} genes")

    # Extract cells from the H5AD
    cells_dict = get_cells_from_h5ad(adata, args.mixed_species)
    cell_barcodes = cells_dict['total']  # Will raise KeyError if not found
    print(f"Found {len(cell_barcodes)} cells in H5AD")

    # Get counts in cells vs out of cells
    counts_in_cells, counts_out_of_cells = get_counts_by_cell_status(adata, args.sample_id)
    print(f"Counts in cells: {counts_in_cells:,}")
    print(f"Counts out of cells: {counts_out_of_cells:,}")

    # Strip sample prefix from cell barcodes for faster BAM processing
    prefix = f"{args.sample_id}_"
    cell_barcodes_stripped = {bc.replace(prefix, '') for bc in cell_barcodes if bc.startswith(prefix)}

    # Parse STAR BAM for read categorization
    reads_in_cells = 0
    reads_out_of_cells = 0
    total_reads_processed = 0
    skipped_secondary = 0
    skipped_supplementary = 0

    print(f"Processing BAM file: {args.star_bam}")
    with pysam.AlignmentFile(args.star_bam, 'rb') as bam:
        for read in bam:
            # Skip secondary and supplementary alignments
            if read.is_secondary:
                skipped_secondary += 1
                continue
            if read.is_supplementary:
                skipped_supplementary += 1
                continue

            total_reads_processed += 1

            # Extract barcode from read name
            barcode = extract_barcode_from_read_name(read.query_name, args.barcode_length)

            if barcode is None:
                print(f"Warning: Failed to extract barcode from read: {read.query_name}")
                continue

            # Check if barcode is a cell
            if barcode in cell_barcodes_stripped:
                reads_in_cells += 1
            else:
                reads_out_of_cells += 1

            # Progress update
            if total_reads_processed % 1000000 == 0:
                print(f"Processed {total_reads_processed:,} reads...")

    print(f"Total reads processed from BAM: {total_reads_processed:,}")
    print(f"Skipped secondary alignments: {skipped_secondary:,}")
    print(f"Skipped supplementary alignments: {skipped_supplementary:,}")
    print(f"Reads in cells: {reads_in_cells:,}")
    print(f"Reads out of cells: {reads_out_of_cells:,}")

    # Calculate unusable reads
    total_mapped = reads_in_cells + reads_out_of_cells
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

    # Validation check
    if reads_in_cells + reads_out_of_cells + unusable_reads != raw_reads:
        print(f"WARNING: Read counts don't add up to raw_reads")
        print(f"  raw_reads: {raw_reads}")
        print(f"  sum: {reads_in_cells + reads_out_of_cells + unusable_reads}")

if __name__ == '__main__':
    main()