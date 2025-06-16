#!/usr/bin/env python3
"""
Create a curated barcode correction list by removing collision-prone 1-Hamming distance variants.

This script:
1. Reads the original barcode list CSV
2. Generates all 1-Hamming distance variants for each barcode
3. Identifies variants that appear for multiple original barcodes (collisions)
4. Creates a TSV where collision-prone variants are removed
5. Ensures barcodes with no remaining corrections are still included for exact matching

Output format is compatible with UMI-tools whitelist:
- Tab-separated values
- Column 1: Original barcode
- Column 2: Comma-separated list of 1-Hamming distance variants (or empty if no safe variants)
"""

import argparse
import csv
from collections import defaultdict
from typing import Dict, List, Set


def generate_1_hamming_variants(barcode: str) -> List[str]:
    """Generate all 1-Hamming distance variants of a barcode."""
    variants = []
    nucleotides = ['N', 'A', 'C', 'G', 'T']
    
    for i in range(len(barcode)):
        original_char = barcode[i]
        for replacement in nucleotides:
            if replacement != original_char:
                variant = barcode[:i] + replacement + barcode[i+1:]
                variants.append(variant)
    
    return variants


def main():
    parser = argparse.ArgumentParser(description='Create curated barcode correction list')
    parser.add_argument('input_csv', help='Input barcode CSV file')
    parser.add_argument('output_tsv', help='Output correction TSV file')
    parser.add_argument('--report', help='Optional collision report file')
    
    args = parser.parse_args()
    
    # Read original barcodes
    print("Reading original barcode list...")
    original_barcodes = []
    
    with open(args.input_csv, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) >= 2:
                barcode = row[1].strip()  # Second column contains barcode
                if barcode:  # Skip empty barcodes
                    original_barcodes.append(barcode)
    
    print(f"Found {len(original_barcodes)} original barcodes")
    
    # Generate all variants and track their frequency
    print("Generating 1-Hamming distance variants...")
    variant_to_origins = defaultdict(list)  # variant -> list of original barcodes that generate it
    barcode_to_variants = {}  # original barcode -> list of its variants
    
    for barcode in original_barcodes:
        variants = generate_1_hamming_variants(barcode)
        barcode_to_variants[barcode] = variants
        
        for variant in variants:
            variant_to_origins[variant].append(barcode)
    
    # Identify collision-prone variants (appear for multiple original barcodes)
    print("Identifying collision-prone variants...")
    collision_variants = set()
    collision_report = []
    
    for variant, origins in variant_to_origins.items():
        if len(origins) > 1:
            collision_variants.add(variant)
            collision_report.append({
                'variant': variant,
                'origins': origins,
                'count': len(origins)
            })
    
    print(f"Found {len(collision_variants)} collision-prone variants")
    
    # Create curated correction list
    print("Creating curated correction list...")
    curated_corrections = {}
    stats = {
        'total_barcodes': len(original_barcodes),
        'barcodes_with_corrections': 0,
        'barcodes_without_corrections': 0,
        'total_variants_removed': len(collision_variants),
        'total_safe_variants': 0
    }
    
    for barcode in original_barcodes:
        # Filter out collision-prone variants
        safe_variants = [v for v in barcode_to_variants[barcode] if v not in collision_variants]
        curated_corrections[barcode] = safe_variants
        
        if safe_variants:
            stats['barcodes_with_corrections'] += 1
            stats['total_safe_variants'] += len(safe_variants)
        else:
            stats['barcodes_without_corrections'] += 1
    
    # Write output TSV
    print(f"Writing output to {args.output_tsv}...")
    with open(args.output_tsv, 'w') as f:
        for barcode in original_barcodes:
            variants = curated_corrections[barcode]
            variants_str = ','.join(variants) if variants else ''
            f.write(f"{barcode}\t{variants_str}\n")
    
    # Write collision report if requested
    if args.report:
        print(f"Writing collision report to {args.report}...")
        with open(args.report, 'w') as f:
            f.write("Barcode Correction Collision Report\n")
            f.write("====================================\n\n")
            
            f.write(f"Total original barcodes: {stats['total_barcodes']}\n")
            f.write(f"Barcodes with safe corrections: {stats['barcodes_with_corrections']}\n")
            f.write(f"Barcodes without corrections: {stats['barcodes_without_corrections']}\n")
            f.write(f"Total collision-prone variants removed: {stats['total_variants_removed']}\n")
            f.write(f"Total safe variants retained: {stats['total_safe_variants']}\n\n")
            
            f.write("Collision Details:\n")
            f.write("-" * 50 + "\n")
            
            # Sort by collision count (most problematic first)
            collision_report.sort(key=lambda x: x['count'], reverse=True)
            
            for item in collision_report:
                f.write(f"Variant: {item['variant']} (appears in {item['count']} barcodes)\n")
                f.write(f"  Origins: {', '.join(item['origins'])}\n\n")
    
    # Print summary
    print("\nSummary:")
    print(f"  Total original barcodes: {stats['total_barcodes']}")
    print(f"  Barcodes with safe corrections: {stats['barcodes_with_corrections']}")
    print(f"  Barcodes without corrections: {stats['barcodes_without_corrections']}")
    print(f"  Collision-prone variants removed: {stats['total_variants_removed']}")
    print(f"  Safe variants retained: {stats['total_safe_variants']}")
    
    print(f"\nOutput written to: {args.output_tsv}")
    if args.report:
        print(f"Collision report written to: {args.report}")


if __name__ == '__main__':
    main()