#!/usr/bin/env python3
"""
This code extracts genetic variants and creates txt files for each sample with minimal depth 5
"""

import h5py
import os
import time
import sys

if len(sys.argv) < 4:
    print("Usage: python script.py <hdf5_file> <parent_file> <output_folder>") # Check that arguments are ok
    sys.exit(1)
hdf5_file = sys.argv[1]
parent_file = sys.argv[2]
output_folder = sys.argv[3]
if not os.path.exists(hdf5_file): # Check that file exist
    print(f"Error: HDF5 file not found: {hdf5_file}")
    sys.exit(1)
if not os.path.exists(parent_file):
    print(f"Error: Parent file not found: {parent_file}")
    sys.exit(1)
# Create output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

print("=" * 60)
print("VARIANT PROCESSOR")
print("=" * 60)
print(f"HDF5 file:     {hdf5_file}")
print(f"Parent file:   {parent_file}")
print(f"Output folder: {output_folder}")
print("=" * 60)

start_time = time.time()
# Read the parent file
print("\n[1/4] Reading parent file...")

sample_list = []  # Will store: [(sample_name, parent1, parent2)]

with open(parent_file, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'): # Skip comment lines
            continue
        parts = line.split() # Split the line into parts
        # Create at least 3 columns sample, parent1, parent2
        if len(parts) >= 3:
            sample_name = parts[0]
            parent1 = parts[1]
            parent2 = parts[2]
            sample_list.append((sample_name, parent1, parent2))
print(f"      Found {len(sample_list)} samples") # How many samples found

# Open the HDF5 file and get basic info
print("\n[2/4] Reading HDF5 file...")
h5f = h5py.File(hdf5_file, 'r')
# Get the list of all samples in the HDF5 file
all_samples = h5f['samples'][:]
all_samples = [s.decode() if isinstance(s, bytes) else s for s in all_samples] # Convert from bytes to string
sample_index = {} # Create sample dictionary
for i, sample_name in enumerate(all_samples):
    sample_index[sample_name] = i
total_variants = h5f['calldata/GT'].shape[0] #Get total number of variants
print(f"      {len(all_samples)} samples in HDF5")
print(f"      {total_variants:,} variants in HDF5")
# Get chromosome and position data
chromosomes = h5f['variants/CHROM'][:]
positions = h5f['variants/POS'][:]
chromosomes = [c.decode() if isinstance(c, bytes) else c for c in chromosomes]# Convert from bytes to string

#  Process each sample one by one
print(f"\n[3/4] Processing samples...")
total_processed = 0
sample_results = []  # Store results for summary

for sample_name, parent1, parent2 in sample_list:
    # Check if this sample exists in the HDF5 file
    if sample_name not in sample_index:
        print(f"  WARNING: {sample_name} not found in HDF5 file, skipping...")
        continue
    print(f"  Processing: {sample_name}")
    # Get the column index for this sample
    sample_col = sample_index[sample_name]
    output_file = os.path.join(output_folder, f"{sample_name}_{parent1}_{parent2}.txt")
    f_out = open(output_file, 'w')
    f_out.write("Sample_ID\tChromosome\tPosition\tParent1\tParent2\tGenotype\tAllelic_depth\tTotal_depth\n")
    # Get all data for this sample
    genotypes = h5f['calldata/GT'][:, sample_col, :]  # Get all genotype data
    # Get allelic depth data if it exists
    try:
        allelic_depths = h5f['calldata/AD'][:, sample_col, :]  # Ref/Alt allelic depth
    except:
        # If allelic depth doesn't exist, use -1 for all
        allelic_depths = [[-1, -1]] * total_variants
    try:
        depths = h5f['calldata/DP'][:, sample_col] # Total depth data
    except:
        # If depth data doesn't exist, use -1 for all
        depths = [-1] * total_variants

    # Process each variant
    variant_count = 0

    for i in range(total_variants):
        chrom = chromosomes[i]
        pos = positions[i]
        genotype = genotypes[i]
        depth = depths[i]
        ad = allelic_depths[i]  # This is an array like [ref_depth, alt_depth]
        if depth < 5: # Filter on depth >= 5
            continue
        # Format the genotype as "0/1" or "./." if missing
        if -1 in genotype:
            genotype_string = "./."
        else:
            genotype_string = f"{genotype[0]}/{genotype[1]}"

        # Format allelic depth as "ref,alt" (e.g., "10,5")
        if len(ad) >= 2:
            ad_string = f"{ad[0]}/{ad[1]}"
        else:
            ad_string = "./."

        # Write this variant to the file
        f_out.write(f"{sample_name}\t{chrom}\t{pos}\t{parent1}\t{parent2}\t{genotype_string}\t{ad_string}\t{depth}\n")
        variant_count += 1

    f_out.close()
    print(f"    Done! ({variant_count} variants)")
    total_processed += variant_count
    sample_results.append((sample_name, parent1, parent2, variant_count))

# Close the HDF5 file
h5f.close()

# Create summary file
print("\n[4/4] Creating summary file...")
summary_file = os.path.join(output_folder, "summary.txt")
f_summary = open(summary_file, 'w')
f_summary.write("# Sample\tParent1\tParent2\tVariants_Processed\n")

for sample_name, parent1, parent2, variant_count in sample_results:
    f_summary.write(f"{sample_name}\t{parent1}\t{parent2}\t{variant_count}\n")

f_summary.close()

# Print summary
runtime = time.time() - start_time
print("\n" + "=" * 60)
print("DONE!")
print("=" * 60)
print(f"Total variants processed: {total_processed:,}")
print(f"Runtime: {runtime:.1f} seconds ({runtime/60:.1f} minutes)")
print(f"Results saved in: {output_folder}")
print(f"Summary file: {summary_file}")
print("=" * 60)
