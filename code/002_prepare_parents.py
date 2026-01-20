#!/usr/bin/env python3
"""
This code classifies parent genotypes and calculates theta (expected allele frequency)
"""

import h5py
import os
import time
import sys

# Read the command line arguments
if len(sys.argv) < 4:
    print("Usage: python script.py <hdf5_file> <parent_file> <output_folder>")
    sys.exit(1)

hdf5_file = sys.argv[1]
parent_file = sys.argv[2]
output_folder = sys.argv[3]

# Check that files exist
if not os.path.exists(hdf5_file):
    print(f"Error: HDF5 file not found: {hdf5_file}")
    sys.exit(1)
if not os.path.exists(parent_file):
    print(f"Error: Parent file not found: {parent_file}")
    sys.exit(1)

# Create output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

print("=" * 60)
print("PARENT CLASSIFIER")
print("=" * 60)
print(f"HDF5 file:     {hdf5_file}")
print(f"Parent file:   {parent_file}")
print(f"Output folder: {output_folder}")
print("=" * 60)

start_time = time.time()

# STEP 1: Read the parent file to get all unique parent pairs
print("\n[1/4] Reading parent file...")

parent_pairs = []  # Will store: [(parent1, parent2), ...]

with open(parent_file, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 3:
            parent1 = parts[1]
            parent2 = parts[2]
            # Add this pair if we haven't seen it before
            if (parent1, parent2) not in parent_pairs:
                parent_pairs.append((parent1, parent2))

print(f"      Found {len(parent_pairs)} unique parent pairs")

# STEP 2: Open the HDF5 file and get basic info
print("\n[2/4] Reading HDF5 file...")

h5f = h5py.File(hdf5_file, 'r')

# Get the list of all samples in the HDF5 file
all_samples = h5f['samples'][:]
all_samples = [s.decode() if isinstance(s, bytes) else s for s in all_samples]

# Create a dictionary to quickly find sample index
sample_index = {}
for i, sample_name in enumerate(all_samples):
    sample_index[sample_name] = i

total_variants = h5f['calldata/GT'].shape[0]
print(f"      {len(all_samples)} samples in HDF5")
print(f"      {total_variants:,} variants in HDF5")

# Get chromosome and position data
chromosomes = h5f['variants/CHROM'][:]
positions = h5f['variants/POS'][:]
chromosomes = [c.decode() if isinstance(c, bytes) else c for c in chromosomes]

# STEP 3: Process each parent pair
print(f"\n[3/4] Processing parent pairs...")

pair_results = []  # Store results for summary

for parent1, parent2 in parent_pairs:
    # Check if both parents exist in the HDF5 file
    if parent1 not in sample_index or parent2 not in sample_index:
        print(f"  WARNING: {parent1} or {parent2} not found in HDF5, skipping...")
        continue

    print(f"  Processing: {parent1} + {parent2}")

    # Get the column indices for these parents
    parent1_col = sample_index[parent1]
    parent2_col = sample_index[parent2]

    # Create output file
    output_file = os.path.join(output_folder, f"{parent1}_{parent2}.txt")
    f_out = open(output_file, 'w')
    f_out.write("Chromosome\tPosition\tClassification\tTheta\n")

    # Get all genotype data for both parents
    parent1_genotypes = h5f['calldata/GT'][:, parent1_col, :]
    parent2_genotypes = h5f['calldata/GT'][:, parent2_col, :]

    # Counters for each classification type
    count_0_0 = 0
    count_0_1 = 0
    count_1_1 = 0
    total_classified = 0

    # Process each variant
    for i in range(total_variants):
        chrom = chromosomes[i]
        pos = positions[i]

        # Get genotypes for both parents
        p1_gt = parent1_genotypes[i]  # e.g., [0, 0] or [1, 1] or [0, 1]
        p2_gt = parent2_genotypes[i]

        # Skip if either parent has missing data
        if -1 in p1_gt or -1 in p2_gt:
            continue

        # Sort genotypes to make comparison easier
        p1_sorted = tuple(sorted(p1_gt))
        p2_sorted = tuple(sorted(p2_gt))

        # Classify based on parent genotypes
        classification = None
        theta = None

        # Case 1: Both parents are 0/0 → progeny should be 0/0
        if p1_sorted == (0, 0) and p2_sorted == (0, 0):
            classification = "0/0"
            theta = 0.0001  # Very low chance of alt allele
            count_0_0 += 1

        # Case 2: Both parents are 1/1 → progeny should be 1/1
        elif p1_sorted == (1, 1) and p2_sorted == (1, 1):
            classification = "1/1"
            theta = 0.9999  # Very high chance of alt allele
            count_1_1 += 1

        # Case 3: One parent is 0/0, other is 1/1 → progeny should be 0/1
        elif (p1_sorted == (0, 0) and p2_sorted == (1, 1)) or \
             (p1_sorted == (1, 1) and p2_sorted == (0, 0)):
            classification = "0/1"
            theta = 0.5  # 50% chance of alt allele
            count_0_1 += 1

        # Skip other cases (e.g., one or both parents are heterozygous)
        else:
            continue

        # Write to file
        f_out.write(f"{chrom}\t{pos}\t{classification}\t{theta}\n")
        total_classified += 1

    f_out.close()

    print(f"    Done! ({total_classified} variants: {count_0_0} as 0/0, {count_0_1} as 0/1, {count_1_1} as 1/1)")
    pair_results.append((parent1, parent2, total_classified, count_0_0, count_0_1, count_1_1))

# Close the HDF5 file
h5f.close()

# STEP 4: Create summary file
print("\n[4/4] Creating summary file...")

summary_file = os.path.join(output_folder, "summary.txt")
f_summary = open(summary_file, 'w')
f_summary.write("# Parent1\tParent2\tTotal_Variants\tCount_0/0\tCount_0/1\tCount_1/1\n")

for parent1, parent2, total, c00, c01, c11 in pair_results:
    f_summary.write(f"{parent1}\t{parent2}\t{total}\t{c00}\t{c01}\t{c11}\n")

f_summary.close()

# Print summary
runtime = time.time() - start_time
print("\n" + "=" * 60)
print("DONE!")
print("=" * 60)
print(f"Processed {len(pair_results)} parent pairs")
print(f"Runtime: {runtime:.1f} seconds ({runtime/60:.1f} minutes)")
print(f"Results saved in: {output_folder}")
print(f"Summary file: {summary_file}")
print("=" * 60)
