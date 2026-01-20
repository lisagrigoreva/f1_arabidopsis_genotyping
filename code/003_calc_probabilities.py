#!/usr/bin/env python3
"""
Simple Parent Finder - With Probability and Log-Likelihood
For each progeny, test all parent pairs and output top 20 matches with probabilities
"""

import sys
import os
from scipy.stats import binom
from scipy.special import logsumexp
import glob
import math
import numpy as np


if len(sys.argv) < 4:
    print("Usage: python simple_intersection.py <progeny_file> <parent_folder> <output_file>")
    print("\nExample: python simple_intersection.py 299528_6909_9938.txt parent_classifications/ results_299528.txt")
    sys.exit(1)

progeny_file = sys.argv[1]
parent_folder = sys.argv[2]
output_file = sys.argv[3]

print("=" * 70)
print("SIMPLE PARENT FINDER - PROBABILITY METHOD")
print("=" * 70)
print(f"Progeny file:     {progeny_file}")
print(f"Parent folder:    {parent_folder}")
print(f"Output:           {output_file}")
print("=" * 70)

# ============================================================================
# READ PROGENY FILE
# ============================================================================

print("\n[1/4] Reading progeny file...")

progeny_variants = {}

with open(progeny_file, 'r') as f:
    header = f.readline()  # Skip header

    for line in f:
        if line.startswith('#'):
            continue

        parts = line.strip().split('\t')
        if len(parts) < 8:
            continue

        try:
            chrom = parts[1]
            pos = parts[2]
            allelic_depth = parts[6]  # Format: "ref/alt"
            total_depth = int(parts[7])

            # Parse allelic depth
            alt_count = int(allelic_depth.split('/')[1])

            if total_depth > 0:
                key = f"{chrom}:{pos}"
                progeny_variants[key] = (alt_count, total_depth)
        except:
            continue

progeny_id = os.path.basename(progeny_file).replace('.txt', '')
print(f"      Progeny ID: {progeny_id}")
print(f"      Loaded {len(progeny_variants)} variants")

# ============================================================================
# TEST ALL PARENT PAIRS
# ============================================================================

print("\n[2/4] Testing parent pairs...")

parent_files = glob.glob(os.path.join(parent_folder, '*.txt'))
parent_files = [f for f in parent_files if 'summary' not in f]

print(f"      Found {len(parent_files)} parent pairs to test")

results = []

for parent_file in parent_files:
    parent_basename = os.path.basename(parent_file).replace('.txt', '')

    # Extract parent IDs (format: Parent1_Parent2.txt)
    parts = parent_basename.split('_')
    if len(parts) != 2:
        continue

    parent1, parent2 = parts

    # Read parent classifications
    parent_classifications = {}

    with open(parent_file, 'r') as f:
        header = f.readline()  # Skip header

        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue

            try:
                chrom = parts[0]
                pos = parts[1]
                theta = float(parts[3])

                key = f"{chrom}:{pos}"
                parent_classifications[key] = theta
            except:
                continue

    # Calculate log-likelihood
    log_likelihood = 0.0
    num_snps = 0

    for position, (alt_count, total_depth) in progeny_variants.items():
        if position not in parent_classifications:
            continue

        theta = parent_classifications[position]

        # Calculate PMF
        pmf = binom.pmf(alt_count, total_depth, theta)

        # Add to log-likelihood
        if pmf > 0:
            log_likelihood += math.log(pmf)
        else:
            log_likelihood += -1000  # Very negative value for PMF = 0

        num_snps += 1

    if num_snps > 0:
        # Normalize by number of SNPs
        norm_log_likelihood = log_likelihood / num_snps

        results.append({
            'Parent1': parent1,
            'Parent2': parent2,
            'Log_Likelihood': log_likelihood,
            'Norm_Log_Likelihood': norm_log_likelihood,
            'Num_SNPs': num_snps
        })

print(f"      Tested {len(results)} parent pairs")


# Calculate probabilities 
print("\n[3/4] Calculating probabilities...")

if results:
    # Extract TOTAL log-likelihoods (NOT normalized!)
    log_likelihoods = np.array([r['Log_Likelihood'] for r in results])
    # Convert to probabilities using logsumexp for numerical stability
    log_sum = logsumexp(log_likelihoods)
    probabilities = np.exp(log_likelihoods - log_sum)
    # Add probabilities to results
    for i, r in enumerate(results):
        r['Probability'] = probabilities[i]
    # Sort by probability (descending)
    results.sort(key=lambda x: -x['Probability'])

    print(f"      Sum of probabilities: {sum(probabilities):.6f}")
# Top 20 results
print("\n[4/4] Writing results...")

# Write top 20
with open(output_file, 'w') as f:
    f.write(f"Progeny_ID\tParent1\tParent2\tProbability_Percent\tLog_Likelihood\tNorm_Log_Likelihood\tNum_SNPs\tRank\n")

    for rank, r in enumerate(results[:20], 1):
        prob_percent = r['Probability'] * 100  # Convert to percentage
        f.write(f"{progeny_id}\t{r['Parent1']}\t{r['Parent2']}\t"
                f"{prob_percent:.4f}\t{r['Log_Likelihood']:.6f}\t"
                f"{r['Norm_Log_Likelihood']:.6f}\t{r['Num_SNPs']}\t{rank}\n")

## Summary 
print("\n" + "=" * 70)
print("DONE!")
print("=" * 70)
print(f"Progeny:         {progeny_id}")
print(f"Parent pairs:    {len(results)}")
print(f"Top 20 saved:    {output_file}")
print("=" * 70)

if results:
    print("\nTop 5 Matches:")
    print(f"\n{'Rank':<6} {'Parent1':<10} {'Parent2':<10} {'Probability':<15} {'Norm_LL':<12} {'SNPs':<8}")
    print("-" * 70)

    for rank, r in enumerate(results[:5], 1):
        prob_percent = r['Probability'] * 100  # Convert to percentage
        print(f"{rank:<6} {r['Parent1']:<10} {r['Parent2']:<10} "
              f"{prob_percent:>14.4f}% {r['Norm_Log_Likelihood']:>11.4f} {r['Num_SNPs']:>7}")

print("\n" + "=" * 70)
