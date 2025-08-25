#!/usr/bin/env python3
"""
Variant Processor - Extract and classify variants from HDF5 genotype data
"""

import h5py
import numpy as np
import os
import time
import multiprocessing as mp
from functools import partial
import argparse
import sys

def process_sample(entry, h5file_path, output_dir, sample_to_idx, batch_size=5000):
    """Process a single sample to report its variants with parent information"""
    sample_id, parent1_id, parent2_id = entry

    # Skip if sample is not found in HDF5
    if sample_id not in sample_to_idx:
        return entry, [], f"Skipping sample: {sample_id} - Not found in HDF5", {}

    sample_idx = sample_to_idx[sample_id]
    print(f"Processing sample: {sample_id} (parents: {parent1_id}, {parent2_id})")

    # Create output file for this sample
    output_file = os.path.join(output_dir, f"{sample_id}_{parent1_id}_{parent2_id}.txt")
    with open(output_file, 'w') as f:
        f.write(f"# Variant report for sample: {sample_id} (parents: {parent1_id}, {parent2_id})\n")
        f.write("# Sample_ID\tChromosome\tPosition\tParent1\tParent2\tGenotype\tDepth\n")

    variants_processed = 0
    examples = []

    with h5py.File(h5file_path, 'r') as h5f:
        total_variants = h5f['calldata/GT'].shape[0]
        variant_chrom = h5f['variants/CHROM'][:]
        variant_pos = h5f['variants/POS'][:]
        chrom_str = [c.decode() if isinstance(c, bytes) else c for c in variant_chrom]

        # Process variants in batches
        for start_idx in range(0, total_variants, batch_size):
            end_idx = min(start_idx + batch_size, total_variants)

            # Progress reporting
            if start_idx % (batch_size * 10) == 0:
                progress = (start_idx / total_variants) * 100
                print(f"  {sample_id} - Progress: {progress:.1f}% ({start_idx}/{total_variants})")

            # Get genotypes and depths
            sample_gts = h5f['calldata/GT'][start_idx:end_idx, sample_idx, :]
            try:
                sample_depths = h5f['calldata/DP'][start_idx:end_idx, sample_idx]
            except KeyError:
                sample_depths = np.full(end_idx - start_idx, -1)

            # Process batch
            batch_variants = []
            for i in range(end_idx - start_idx):
                variant_idx = start_idx + i
                sample_gt = sample_gts[i]
                depth = sample_depths[i] if sample_depths is not None else -1

                chrom = chrom_str[variant_idx]
                pos = variant_pos[variant_idx]

                # Format genotype
                gt_str = "./." if -1 in sample_gt else f"{sample_gt[0]}/{sample_gt[1]}"
                batch_variants.append((chrom, pos, gt_str, depth))

            # Write batch to file
            if batch_variants:
                with open(output_file, 'a') as f:
                    for chrom, pos, gt_str, depth in batch_variants:
                        f.write(f"{sample_id}\t{chrom}\t{pos}\t{parent1_id}\t{parent2_id}\t{gt_str}\t{depth}\n")
                variants_processed += len(batch_variants)

                # Store examples from first batch
                if not examples:
                    examples = batch_variants[:5]

    message = f"Processed {variants_processed} variants for sample {sample_id}"
    print(message)
    return entry, examples, message, {"variants_processed": variants_processed}


def load_parent_data(parent_file_path):
    """Load sample-parent relationships from file"""
    sample_entries = []
    try:
        with open(parent_file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split()
                if len(parts) >= 3:
                    sample_entries.append((parts[0], parts[1], parts[2]))
                elif len(parts) == 2:
                    sample_entries.append((parts[0], parts[1], parts[1]))
                else:
                    print(f"Warning: Invalid format at line {line_num}: {line}")
    except FileNotFoundError:
        print(f"Error: Parent file not found: {parent_file_path}")
        sys.exit(1)

    return sample_entries


def main():
    parser = argparse.ArgumentParser(
        description="Extract and classify variants from HDF5 genotype data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input HDF5 file containing genotype data'
    )

    parser.add_argument(
        '-p', '--parents',
        required=True,
        help='Parent file (format: sample_id parent1_id [parent2_id])'
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory for results'
    )

    parser.add_argument(
        '--processes',
        type=int,
        default=mp.cpu_count() - 1,
        help='Number of parallel processes'
    )

    parser.add_argument(
        '--batch-size',
        type=int,
        default=5000,
        help='Batch size for processing variants'
    )

    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)

    if not os.path.exists(args.parents):
        print(f"Error: Parent file not found: {args.parents}")
        sys.exit(1)
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    print(f"Input HDF5: {args.input}")
    print(f"Parent file: {args.parents}")
    print(f"Output directory: {args.output}")
    print(f"Using {args.processes} processes")
    start_time = time.time()
    # Load parent data
    sample_entries = load_parent_data(args.parents)
    print(f"Found {len(sample_entries)} sample entries")
    # Create summary file
    summary_file = os.path.join(args.output, "summary.txt")
    with open(summary_file, 'w') as f:
        f.write("# Variant Processing Summary\n")
        f.write(f"# Generated on {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Input: {args.input}\n")
        f.write(f"# Parent file: {args.parents}\n")
        f.write(f"# Found {len(sample_entries)} sample entries\n\n")
    # Load HDF5 sample information
    try:
        with h5py.File(args.input, 'r') as h5f:
            samples = h5f['samples'][:]
            samples_str = [s.decode() if isinstance(s, bytes) else s for s in samples]
            sample_to_idx = {sample: idx for idx, sample in enumerate(samples_str)}
            total_variants = h5f['calldata/GT'].shape[0]

            print(f"HDF5 contains {total_variants} variants and {len(samples_str)} samples")
    except Exception as e:
        print(f"Error reading HDF5 file: {e}")
        sys.exit(1)

    # Process samples in parallel
    print("Starting parallel processing...")
    with mp.Pool(processes=args.processes) as pool:
        process_func = partial(
            process_sample,
            h5file_path=args.input,
            output_dir=args.output,
            sample_to_idx=sample_to_idx,
            batch_size=args.batch_size
        )

        results = pool.map(process_func, sample_entries)

    # Write summary
    with open(summary_file, 'a') as f:
        total_processed = 0
        for entry, examples, message, stats in results:
            sample_id, parent1_id, parent2_id = entry
            output_file = f"{sample_id}_{parent1_id}_{parent2_id}.txt"

            f.write(f"{message}\n")
            f.write(f"Output file: {output_file}\n")

            if examples and args.verbose:
                f.write("Example variants:\n")
                for chrom, pos, gt_str, depth in examples:
                    f.write(f"  {chrom}:{pos} - {gt_str} (depth: {depth})\n")
            f.write("\n")

            total_processed += stats.get('variants_processed', 0)

        # Final summary
        runtime = time.time() - start_time
        f.write(f"Total variants processed: {total_processed}\n")
        f.write(f"Total runtime: {runtime:.2f} seconds ({runtime/60:.2f} minutes)\n")

    print(f"\nProcessing complete!")
    print(f"Total runtime: {runtime:.2f} seconds ({runtime/60:.2f} minutes)")
    print(f"Results saved to: {args.output}")
    print(f"Summary file: {summary_file}")


if __name__ == "__main__":
    main()
