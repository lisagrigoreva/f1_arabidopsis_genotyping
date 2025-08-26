#!/usr/bin/env python3
"""
Create Likelihood Table from Intersection Files with PMF
Processes intersection files that already contain PMF values and creates likelihood tables.
"""

import os
import sys
import glob
import numpy as np
import math
import pandas as pd
from scipy.special import logsumexp
import argparse
from io import StringIO

def process_intersection_file(file_path):
    """
    Process an intersection file with PMF column and return the parents, log-likelihood, and number of SNPs.
    """
    try:
        # Extract parent IDs from filename
        basename = os.path.basename(file_path)
        if basename.startswith('intersect_'):
            parts = basename.replace('intersect_', '').replace('.txt', '').split('_')
            if len(parts) >= 2:
                parent1, parent2 = parts[0], parts[1]
            else:
                # Try to extract from file header
                with open(file_path, 'r') as f:
                    header = f.readline().strip()
                    if '(parents:' in header:
                        sample_info = header.split('(parents: ')[1].split(')')[0]
                        parent1, parent2 = sample_info.split(', ')
                    else:
                        return None, None, 0, 0

        # Read the file - skip the first two comment lines, then read with header
        with open(file_path, 'r') as f:
            lines = f.readlines()

        # Find the header line (starts with # and contains column names)
        header_line = None
        data_start = 0
        for i, line in enumerate(lines):
            if line.startswith('# Sample_ID'):
                header_line = line.strip().replace('# ', '')
                data_start = i + 1
                break

        if header_line is None:
            print(f"Warning: Header line not found in {file_path}")
            return None, None, 0, 0

        # Create a temporary file content with proper header
        temp_content = header_line + '\n' + ''.join(lines[data_start:])

        # Read with pandas using StringIO
        df = pd.read_csv(StringIO(temp_content), sep='\t')

        if 'PMF' not in df.columns:
            print(f"Warning: PMF column not found in {file_path}")
            print(f"Available columns: {df.columns.tolist()}")
            return None, None, 0, 0

        # Calculate log-likelihood
        log_likelihood = 0
        num_snps = 0

        for _, row in df.iterrows():
            try:
                pmf = float(row['PMF'])
                if pmf > 0:  # Avoid log(0)
                    log_likelihood += math.log(pmf)
                else:
                    # Handle PMF = 0 by using a very small value
                    log_likelihood += -1000  # A very negative value
                num_snps += 1
            except (ValueError, TypeError):
                pass

        return parent1, parent2, log_likelihood, num_snps

    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return None, None, 0, 0

def process_intersection_directory(dir_path, progeny_id=None):
    """
    Process all intersection files in a specific directory.
    """
    if not os.path.exists(dir_path):
        print(f"Directory not found: {dir_path}")
        return False

    # Extract progeny ID from directory name if not provided
    if progeny_id is None:
        dir_name = os.path.basename(dir_path)
        if dir_name.endswith('_intersections'):
            progeny_id = dir_name.replace('_intersections', '')
        else:
            progeny_id = dir_name

    print(f"Processing directory: {dir_path}")
    print(f"Progeny ID: {progeny_id}")

    # Find all intersection files (intersect_*.txt)
    intersection_files = glob.glob(os.path.join(dir_path, 'intersect_*.txt'))

    if not intersection_files:
        print(f"No intersection files found in {dir_path}")
        return True

    print(f"Found {len(intersection_files)} intersection files to process")

    # Create a list to store results
    results = []

    # Process each file
    for file_path in intersection_files:
        parent1, parent2, log_likelihood, num_snps = process_intersection_file(file_path)
        if parent1 and parent2:
            # Normalize log-likelihood by number of SNPs
            norm_log_likelihood = log_likelihood / num_snps if num_snps > 0 else 0
            results.append({
                'Progeny': progeny_id,
                'Parent1': parent1,
                'Parent2': parent2,
                'Log_Likelihood': log_likelihood,
                'Num_SNPs': num_snps,
                'Norm_Log_Likelihood': norm_log_likelihood,
                'File': os.path.basename(file_path)
            })

    # Convert to probabilities
    if results:
        df = pd.DataFrame(results)

        # Extract normalized log-likelihoods
        norm_log_likelihoods = df['Norm_Log_Likelihood'].values

        # Convert to probabilities using logsumexp for numerical stability
        # P(parent_pair | data) = exp(log_likelihood) / sum(exp(all_log_likelihoods))
        log_sum = logsumexp(norm_log_likelihoods)
        probabilities = np.exp(norm_log_likelihoods - log_sum)

        # Add probability column
        df['Probability'] = probabilities

        # Sort by probability (descending)
        df = df.sort_values('Probability', ascending=False)

        # Save to CSV
        output_file = os.path.join(dir_path, f"{progeny_id}_likelihood_table.csv")
        df.to_csv(output_file, index=False)
        print(f"Saved likelihood table to {output_file}")

        # Print top 5 results
        print("\nTop 5 most likely parent pairs:")
        print(df.head(5)[['Parent1', 'Parent2', 'Probability', 'Norm_Log_Likelihood', 'Num_SNPs']].to_string(index=False))

        # Print sum of probabilities (should be 1.0)
        print(f"\nSum of all probabilities: {df['Probability'].sum():.6f}")

        return df
    else:
        print("No valid results found")
        return None

def process_intersection_directory(dir_path, progeny_id=None):
    """
    Process all intersection files in a specific directory.
    """
    if not os.path.exists(dir_path):
        print(f"Directory not found: {dir_path}")
        return None

    # Extract progeny ID from directory name if not provided
    if progeny_id is None:
        dir_name = os.path.basename(dir_path)
        if dir_name.endswith('_intersections'):
            progeny_id = dir_name.replace('_intersections', '')
        else:
            progeny_id = dir_name

    print(f"Processing directory: {dir_path}")
    print(f"Progeny ID: {progeny_id}")

    # Find all intersection files (intersect_*.txt)
    intersection_files = glob.glob(os.path.join(dir_path, 'intersect_*.txt'))

    if not intersection_files:
        print(f"No intersection files found in {dir_path} - skipping")
        return None  # Changed from True to None

    print(f"Found {len(intersection_files)} intersection files to process")

    # Create a list to store results
    results = []

    # Process each file
    for file_path in intersection_files:
        parent1, parent2, log_likelihood, num_snps = process_intersection_file(file_path)
        if parent1 and parent2:
            # Normalize log-likelihood by number of SNPs
            norm_log_likelihood = log_likelihood / num_snps if num_snps > 0 else 0
            results.append({
                'Progeny': progeny_id,
                'Parent1': parent1,
                'Parent2': parent2,
                'Log_Likelihood': log_likelihood,
                'Num_SNPs': num_snps,
                'Norm_Log_Likelihood': norm_log_likelihood,
                'File': os.path.basename(file_path)
            })

    # Convert to probabilities
    if results:
        df = pd.DataFrame(results)

        # Extract normalized log-likelihoods
        norm_log_likelihoods = df['Norm_Log_Likelihood'].values

        # Convert to probabilities using logsumexp for numerical stability
        # P(parent_pair | data) = exp(log_likelihood) / sum(exp(all_log_likelihoods))
        log_sum = logsumexp(norm_log_likelihoods)
        probabilities = np.exp(norm_log_likelihoods - log_sum)

        # Add probability column
        df['Probability'] = probabilities

        # Sort by probability (descending)
        df = df.sort_values('Probability', ascending=False)

        # Save to CSV
        output_file = os.path.join(dir_path, f"{progeny_id}_likelihood_table.csv")
        df.to_csv(output_file, index=False)
        print(f"Saved likelihood table to {output_file}")

        # Print top 5 results
        print("\nTop 5 most likely parent pairs:")
        print(df.head(5)[['Parent1', 'Parent2', 'Probability', 'Norm_Log_Likelihood', 'Num_SNPs']].to_string(index=False))

        # Print sum of probabilities (should be 1.0)
        print(f"\nSum of all probabilities: {df['Probability'].sum():.6f}")

        return df
    else:
        print("No valid results found")
        return None

def process_all_directories(base_dir):
    """
    Process all intersection directories in the base directory.
    """
    # Get all intersection directories
    intersection_dirs = [d for d in os.listdir(base_dir)
                        if d.endswith('_intersections') and os.path.isdir(os.path.join(base_dir, d))]

    if not intersection_dirs:
        print(f"No intersection directories found in {base_dir}")
        return

    # Sort directories to ensure consistent ordering
    intersection_dirs.sort()

    print(f"Found {len(intersection_dirs)} intersection directories to process")

    all_results = []
    processed_count = 0
    skipped_count = 0

    # Process each directory
    for i, dir_name in enumerate(intersection_dirs, 1):
        dir_path = os.path.join(base_dir, dir_name)
        print(f"\n{'='*60}")
        print(f"Processing directory {i} of {len(intersection_dirs)}: {dir_name}")
        print(f"{'='*60}")

        df = process_intersection_directory(dir_path)
        if df is not None:
            # Keep top result for summary
            top_result = df.iloc[0].copy()
            all_results.append(top_result)
            processed_count += 1
        else:
            # Directory was empty or had no valid results
            skipped_count += 1

    # Create summary table
    if all_results:
        summary_df = pd.DataFrame(all_results)
        summary_file = os.path.join(base_dir, "likelihood_summary.csv")
        summary_df.to_csv(summary_file, index=False)
        print(f"\n{'='*60}")
        print(f"SUMMARY - Top matches for all progeny:")
        print(f"{'='*60}")
        print(summary_df[['Progeny', 'Parent1', 'Parent2', 'Probability', 'Num_SNPs']].to_string(index=False))
        print(f"\nSummary saved to: {summary_file}")
        print(f"Processed: {processed_count} directories")
        print(f"Skipped: {skipped_count} empty directories")
    else:
        print(f"\nNo valid results found in any directory")
        print(f"Skipped: {skipped_count} empty directories")

def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Create likelihood tables from intersection files with PMF values",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "input_path",
        help="Path to base directory containing intersection directories, or path to single intersection directory"
    )

    parser.add_argument(
        "--single", "-s",
        action="store_true",
        help="Process single intersection directory instead of all directories"
    )

    parser.add_argument(
        "--progeny-id", "-p",
        help="Progeny ID (only used with --single option)"
    )

    args = parser.parse_args()

    # Validate input path
    if not os.path.exists(args.input_path):
        print(f"Error: Path not found: {args.input_path}", file=sys.stderr)
        sys.exit(1)

    try:
        if args.single:
            # Process single directory
            result = process_intersection_directory(args.input_path, args.progeny_id)
            if result is None:  # Changed from 'if not result:' to 'if result is None:'
                sys.exit(1)
        else:
            # Process all directories
            process_all_directories(args.input_path)

        print("\nLikelihood table creation completed successfully")

    except Exception as e:
        print(f"Error during analysis: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
