#!/usr/bin/env python3
"""
SNP Classification Script
Classifies SNPs based on parent genotypes from HDF5 files.
"""

import h5py
import numpy as np
import os
import time
import multiprocessing as mp
from functools import partial
import argparse
import sys

# Classification mapping for cleaner code
CLASSIFICATION_RULES = {
    ((0, 0), (0, 0)): ("0/0", 0.0001),
    ((1, 1), (1, 1)): ("1/1", 0.9999),
    ((0, 0), (1, 1)): ("0/1", 0.5),
    ((1, 1), (0, 0)): ("1/0", 0.5)
}

def classify_genotype_pair(parent1_gt, parent2_gt):
    """Classify a pair of parent genotypes"""
    key = (tuple(sorted(parent1_gt)), tuple(sorted(parent2_gt)))
    return CLASSIFICATION_RULES.get(key, (None, None))

def write_file_header(file_path, parent1_id, parent2_id):
    """Write header to output file"""
    with open(file_path, 'w') as f:
        f.write(f"# SNP classification for parent pair: {parent1_id} and {parent2_id}\n")
        f.write("# Chromosome\tPosition\tClassification\tProbability\n")

def write_batch_results(file_path, batch_positions):
    """Append batch results to file"""
    with open(file_path, 'a') as f:
        for chrom, pos, classification, probability in batch_positions:
            f.write(f"{chrom}\t{pos}\t{classification}\t{probability}\n")

def process_parent_pair(parent_pair, h5file_path, output_dir, sample_to_idx,
                        variant_chrom, variant_pos, chrom_str, total_variants, batch_size=1000):
    """Process a single parent pair to classify SNPs based on parent genotypes"""
    parent1_id, parent2_id = parent_pair

    # Initialize counts
    counts = {"0/0": 0, "0/1": 0, "1/0": 0, "1/1": 0}

    # Skip if any parent is not found
    if parent1_id not in sample_to_idx or parent2_id not in sample_to_idx:
        message = f"Skipping pair: {parent1_id} and {parent2_id} - Not all parents found in HDF5"
        return parent_pair, [], message, counts

    # Get parent indices
    parent1_idx = sample_to_idx[parent1_id]
    parent2_idx = sample_to_idx[parent2_id]

    print(f"Processing parent pair: {parent1_id} and {parent2_id}")

    # Create output file
    pair_key = f"{parent1_id}_{parent2_id}"
    pair_file = os.path.join(output_dir, f"{pair_key}.txt")
    write_file_header(pair_file, parent1_id, parent2_id)

    positions = []

    # Process variants in batches
    with h5py.File(h5file_path, 'r') as h5f:
        for start_idx in range(0, total_variants, batch_size):
            end_idx = min(start_idx + batch_size, total_variants)

            # Progress reporting
            if start_idx % (batch_size * 10) == 0:
                progress = (start_idx / total_variants) * 100
                print(f"  {pair_key} - Progress: {progress:.1f}% ({start_idx}/{total_variants})")

            # Get genotypes for both parents
            parent1_gts = h5f['calldata/GT'][start_idx:end_idx, parent1_idx, :]
            parent2_gts = h5f['calldata/GT'][start_idx:end_idx, parent2_idx, :]

            # Process batch
            batch_positions = []
            for i in range(end_idx - start_idx):
                variant_idx = start_idx + i
                parent1_gt = parent1_gts[i]
                parent2_gt = parent2_gts[i]

                # Skip missing data
                if -1 in parent1_gt or -1 in parent2_gt:
                    continue

                # Classify genotype
                classification, probability = classify_genotype_pair(parent1_gt, parent2_gt)
                if classification is None:
                    continue

                # Get position info
                chrom = chrom_str[variant_idx]
                pos = variant_pos[variant_idx]

                batch_positions.append((chrom, pos, classification, probability))
                counts[classification] += 1

            # Write batch results
            if batch_positions:
                write_batch_results(pair_file, batch_positions)
                positions.extend(batch_positions)

    message = f"Classified {len(positions)} SNPs for {pair_key}: " + \
              ", ".join([f"{count} as {cls}" for cls, count in counts.items()])
    print(message)

    # Return examples
    examples = positions[:5] if positions else []
    return parent_pair, examples, message, counts

def load_parent_pairs(parent_file_path):
    """Load parent pairs from file"""
    parent_pairs = []
    with open(parent_file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                parent_pairs.append((parts[1], parts[2]))
            elif len(parts) == 2:
                parent_pairs.append((parts[1], parts[1]))

    # Remove duplicates
    return list(set(parent_pairs))

def write_summary_header(summary_file):
    """Write summary file header"""
    with open(summary_file, 'w') as f:
        f.write("# SNP classification based on parent genotypes\n")
        f.write("# 0/0: both parents have 0/0 genotype\n")
        f.write("# 0/1: proper heterozygous (0/0 maternal) and (1/1 progeny)\n")
        f.write("# 1/0: artificial heterozygous, maternal is alternative\n")
        f.write("# 1/1: both parents have 1/1 genotype\n")
        f.write(f"# Generated on {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

def write_summary_results(summary_file, results, output_dir, runtime):
    """Write results to summary file"""
    with open(summary_file, 'a') as f:
        for parent_pair, examples, message, counts in results:
            parent1_id, parent2_id = parent_pair
            pair_key = f"{parent1_id}_{parent2_id}"
            pair_file = os.path.join(output_dir, f"{pair_key}.txt")

            f.write(message + "\n")

            if examples:
                f.write("Example positions:\n")
                for chrom, pos, classification, probability in examples:
                    f.write(f"  {chrom}:{pos} - {classification}\t{probability}\n")
            f.write(f"Complete list saved to: {pair_file}\n\n")

        f.write(f"Total runtime: {runtime:.2f} seconds ({runtime/60:.2f} minutes)\n")

def find_informative_snps(h5file_path, parent_file_path, output_dir,
                         num_processes=4, batch_size=1000, limit=None):
    """Main function to classify SNPs for all parent pairs"""
    start_time = time.time()

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Setup summary file
    summary_file = os.path.join(output_dir, "summary.txt")
    write_summary_header(summary_file)

    # Load parent pairs
    parent_pairs = load_parent_pairs(parent_file_path)
    print(f"Found {len(parent_pairs)} unique parent pairs")

    with open(summary_file, 'a') as f:
        f.write(f"Found {len(parent_pairs)} unique parent pairs\n\n")

    # Load HDF5 metadata
    with h5py.File(h5file_path, 'r') as h5f:
        samples = h5f['samples'][:]
        samples_str = [s.decode() if isinstance(s, bytes) else s for s in samples]
        sample_to_idx = {sample: idx for idx, sample in enumerate(samples_str)}

        variant_chrom = h5f['variants/CHROM'][:]
        variant_pos = h5f['variants/POS'][:]
        chrom_str = [c.decode() if isinstance(c, bytes) else c for c in variant_chrom]

        total_variants = h5f['calldata/GT'].shape[0]
        if limit:
            total_variants = min(limit, total_variants)

        print(f"Processing {total_variants} variants")
        with open(summary_file, 'a') as f:
            f.write(f"Processing {total_variants} variants\n\n")

    # Process in parallel
    print(f"Using {num_processes} processes")
    with mp.Pool(processes=num_processes) as pool:
        process_func = partial(
            process_parent_pair,
            h5file_path=h5file_path,
            output_dir=output_dir,
            sample_to_idx=sample_to_idx,
            variant_chrom=variant_chrom,
            variant_pos=variant_pos,
            chrom_str=chrom_str,
            total_variants=total_variants,
            batch_size=batch_size
        )

        results = pool.map(process_func, parent_pairs)

    # Write summary
    runtime = time.time() - start_time
    write_summary_results(summary_file, results, output_dir, runtime)

    print(f"Total runtime: {runtime:.2f} seconds ({runtime/60:.2f} minutes)")
    return results

def main():
    """Main function with command line argument parsing"""
    parser = argparse.ArgumentParser(
        description="Classify SNPs based on parent genotypes from HDF5 files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "h5file",
        help="Path to HDF5 file containing genotype data"
    )

    parser.add_argument(
        "parent_file",
        help="Path to file containing parent relationships"
    )

    parser.add_argument(
        "output_dir",
        help="Output directory for results"
    )

    parser.add_argument(
        "--processes", "-p",
        type=int,
        default=mp.cpu_count() - 1,
        help="Number of processes to use"
    )

    parser.add_argument(
        "--batch-size", "-b",
        type=int,
        default=5000,
        help="Batch size for processing variants"
    )

    parser.add_argument(
        "--limit", "-l",
        type=int,
        default=None,
        help="Limit number of variants to process (for testing)"
    )

    args = parser.parse_args()

    # Validate input files
    if not os.path.exists(args.h5file):
        print(f"Error: HDF5 file not found: {args.h5file}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.parent_file):
        print(f"Error: Parent file not found: {args.parent_file}", file=sys.stderr)
        sys.exit(1)

    # Run analysis
    try:
        results = find_informative_snps(
            args.h5file,
            args.parent_file,
            args.output_dir,
            num_processes=args.processes,
            batch_size=args.batch_size,
            limit=args.limit
        )

        print(f"Results saved to {args.output_dir}")

    except Exception as e:
        print(f"Error during analysis: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
