#!/usr/bin/env python3
"""
All Intersections Script - Improved
Creates intersections for each progeny ID with all possible parental combinations.
Uses pandas for robust data handling and column management.
"""

import os
import sys
import pandas as pd
import time
import argparse
import glob
from multiprocessing import Pool, cpu_count
from functools import partial
from scipy.stats import binom

def extract_progeny_id(filename):
    """Extract progeny ID from the filename."""
    parts = os.path.basename(filename).split('_')
    if len(parts) >= 1:
        return parts[0]
    return None

def extract_parent_ids(filename):
    """Extract parent IDs from the parent filename."""
    basename = os.path.basename(filename)
    # Remove .txt extension and split by underscore
    name_without_ext = basename.replace('.txt', '')
    parts = name_without_ext.split('_')
    
    if len(parts) >= 2:
        parent1 = parts[0]
        parent2 = parts[1]
        return parent1, parent2
    return None, None

def read_progeny_file(filepath):
    """Read progeny file with proper column names using pandas"""
    try:
        df = pd.read_csv(filepath, sep='\t', comment='#',
                        names=["Sample_ID", "Chromosome", "Position", "Genotype",
                              "Ref_Count", "Alt_Count", "Total_Depth"])
        return df
    except Exception as e:
        print(f"Error reading progeny file {filepath}: {e}")
        return None

def read_parent_file(filepath):
    """Read parent classification file with proper column names using pandas"""
    try:
        df = pd.read_csv(filepath, sep='\t', comment='#',
                        names=["Chromosome", "Position", "Classification", "Probability"])
        return df
    except Exception as e:
        print(f"Error reading parent file {filepath}: {e}")
        return None

def create_intersection_key(df, chrom_col, pos_col):
    """Create a key for joining dataframes"""
    return df[chrom_col].astype(str) + ':' + df[pos_col].astype(str)

def write_intersection_file(output_file, progeny_id, parent1, parent2, merged_df):
    """Write the intersection results to file with proper headers"""
    with open(output_file, 'w') as f:
        f.write(f"# Intersection for sample: {progeny_id} (parents: {parent1}, {parent2})\n")
        f.write("# Sample_ID\tChromosome\tPosition\tGenotype\tRef_Count\tAlt_Count\tTotal_Depth\tParentClassification\tProbability\tPMF\n")

    # Rename the Classification column to ParentClassification for clarity
    output_df = merged_df.copy()
    output_df = output_df.rename(columns={'Classification': 'ParentClassification'})

    # Reorder columns to match header
    column_order = ["Sample_ID", "Chromosome", "Position", "Genotype",
                   "Ref_Count", "Alt_Count", "Total_Depth", "ParentClassification", "Probability", "PMF"]
    output_df = output_df[column_order]

    output_df.to_csv(output_file, sep='\t', index=False, mode='a', header=False)

def calculate_pmf(alt_count, total_depth, theta):
    """Calculate binomial PMF for given parameters."""
    try:
        alt_count = int(alt_count) if pd.notna(alt_count) else 0
        total_depth = int(total_depth) if pd.notna(total_depth) else 0
        theta = float(theta) if pd.notna(theta) else None
        
        if theta is not None and total_depth > 0:
            pmf = binom.pmf(alt_count, total_depth, theta)
        else:
            pmf = 0
        return pmf
    except (ValueError, TypeError):
        return 0

def intersect_progeny_with_parents(progeny_df, parent_df):
    """Intersect progeny and parent dataframes based on chromosome and position."""
    try:
        # Create keys for joining
        progeny_df = progeny_df.copy()
        parent_df = parent_df.copy()
        
        progeny_df['key'] = create_intersection_key(progeny_df, 'Chromosome', 'Position')
        parent_df['key'] = create_intersection_key(parent_df, 'Chromosome', 'Position')

        # Find intersections using pandas merge
        merged = pd.merge(progeny_df, parent_df[['key', 'Classification', 'Probability']], 
                         on='key', how='inner')
        merged = merged.drop('key', axis=1)

        # Calculate PMF for each row
        merged['PMF'] = merged.apply(
            lambda row: calculate_pmf(row['Alt_Count'], row['Total_Depth'], row['Probability']), 
            axis=1
        )

        return merged
    except Exception as e:
        print(f"Error intersecting dataframes: {e}")
        return pd.DataFrame()

def process_single_progeny(args):
    """Process a single progeny file with all parent files - for multiprocessing."""
    progeny_file, progeny_dir, parent_dir, output_base_dir, parent_files = args
    
    start_time = time.time()
    
    try:
        # Extract progeny ID and create output directory
        progeny_id = extract_progeny_id(progeny_file)
        if not progeny_id:
            print(f"Could not extract progeny ID from {progeny_file}")
            return progeny_file, 0, 0

        progeny_output_dir = os.path.join(output_base_dir, f"{progeny_id}_intersections")
        os.makedirs(progeny_output_dir, exist_ok=True)

        # Read the progeny file
        progeny_path = os.path.join(progeny_dir, progeny_file)
        progeny_df = read_progeny_file(progeny_path)
        if progeny_df is None or progeny_df.empty:
            print(f"Failed to read progeny file {progeny_file}")
            return progeny_file, 0, 0

        print(f"Processing {progeny_file} with {len(progeny_df)} rows against {len(parent_files)} parent combinations")

        # Process each parent file
        successful_intersections = 0
        total_intersections = 0
        
        for parent_file in parent_files:
            parent_path = os.path.join(parent_dir, parent_file)
            parent1, parent2 = extract_parent_ids(parent_file)

            if not parent1 or not parent2:
                continue

            # Read the parent file
            parent_df = read_parent_file(parent_path)
            if parent_df is None or parent_df.empty:
                continue

            # Intersect the dataframes
            merged_df = intersect_progeny_with_parents(progeny_df, parent_df)
            
            if merged_df.empty:
                continue

            # Create output file path
            output_file = os.path.join(progeny_output_dir, f"intersect_{parent1}_{parent2}.txt")

            # Write the intersection results
            write_intersection_file(output_file, progeny_id, parent1, parent2, merged_df)

            successful_intersections += 1
            total_intersections += len(merged_df)

        end_time = time.time()
        elapsed_time = end_time - start_time
        
        print(f"Completed {progeny_file}: {successful_intersections} successful intersections, "
              f"{total_intersections} total matches in {elapsed_time:.2f}s")
        
        return progeny_file, successful_intersections, total_intersections

    except Exception as e:
        print(f"Error processing progeny file {progeny_file}: {e}")
        return progeny_file, 0, 0

def get_progeny_files(progeny_dir, progeny_pattern="*.txt"):
    """Get all progeny files matching the pattern."""
    progeny_pattern_full = os.path.join(progeny_dir, progeny_pattern)
    progeny_files = [os.path.basename(f) for f in glob.glob(progeny_pattern_full)]
    
    # Filter out summary files
    progeny_files = [f for f in progeny_files if "summary" not in f.lower()]
    
    return progeny_files

def get_parent_files(parent_dir):
    """Get all parent files."""
    parent_files = [f for f in os.listdir(parent_dir) 
                   if f.endswith('.txt') and '_' in f and f != 'summary.txt']
    # Filter out summary files and ensure it has the parent1_parent2.txt format
    parent_files = [f for f in parent_files if "summary" not in f.lower()]
    return parent_files

def process_all_intersections(progeny_dir, parent_dir, output_dir, num_processes=None, progeny_pattern="*.txt"):
    """Process all progeny files with all parent combinations."""
    
    if num_processes is None:
        num_processes = max(1, cpu_count() - 1)
    
    print(f"Starting intersection analysis with {num_processes} processes")
    
    # Get all progeny files
    progeny_files = get_progeny_files(progeny_dir, progeny_pattern)
    
    if not progeny_files:
        print(f"No progeny files found matching pattern: {progeny_pattern}")
        return
    
    # Get all parent files
    parent_files = get_parent_files(parent_dir)
    
    if not parent_files:
        print(f"No parent files found in {parent_dir}")
        return
    
    print(f"Found {len(progeny_files)} progeny files and {len(parent_files)} parent combinations")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Prepare arguments for multiprocessing
    args_list = [(pf, progeny_dir, parent_dir, output_dir, parent_files) for pf in progeny_files]
    
    # Process files
    start_time = time.time()
    
    if num_processes == 1:
        # Single process
        results = [process_single_progeny(args) for args in args_list]
    else:
        # Multiprocessing
        with Pool(processes=num_processes) as pool:
            results = pool.map(process_single_progeny, args_list)
    
    # Summary
    total_successful = sum(r[1] for r in results)
    total_intersections = sum(r[2] for r in results)
    end_time = time.time()
    
    print(f"\n=== SUMMARY ===")
    print(f"Processed {len(progeny_files)} progeny files")
    print(f"Total successful intersections: {total_successful}")
    print(f"Total intersection matches: {total_intersections}")
    print(f"Total processing time: {(end_time - start_time):.2f} seconds")
    
    # Write summary file
    write_summary_file(output_dir, progeny_files, parent_files, results, end_time - start_time)
    
def write_summary_file(output_dir, progeny_files, parent_files, results, total_time):
    """Write a comprehensive summary file."""
    summary_file = os.path.join(output_dir, "intersection_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("# Intersection Analysis Summary\n")
        f.write(f"# Generated on {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Progeny files processed: {len(progeny_files)}\n")
        f.write(f"# Parent combinations: {len(parent_files)}\n")
        f.write(f"# Total processing time: {total_time:.2f} seconds\n\n")
        f.write("Progeny_File\tSuccessful_Intersections\tTotal_Matches\n")
        
        for progeny_file, successful, total in results:
            f.write(f"{progeny_file}\t{successful}\t{total}\n")
    
    print(f"Summary saved to: {summary_file}")

def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description="Create intersections for each progeny ID with all possible parental combinations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "progeny_dir",
        help="Directory containing progeny files"
    )
    
    parser.add_argument(
        "parent_dir",
        help="Directory containing parent files (prob_*.txt)"
    )
    
    parser.add_argument(
        "output_dir",
        help="Output directory for intersection results"
    )
    
    parser.add_argument(
        "--processes", "-p",
        type=int,
        default=None,
        help="Number of processes to use (default: CPU count - 1)"
    )
    
    parser.add_argument(
        "--progeny-pattern", "-pp",
        default="*.txt",
        help="Pattern to match progeny files"
    )
    
    parser.add_argument(
        "--single-progeny", "-s",
        help="Process only a single progeny file (filename only, not full path)"
    )
    
    args = parser.parse_args()
    
    # Validate input directories
    if not os.path.exists(args.progeny_dir):
        print(f"Error: Progeny directory not found: {args.progeny_dir}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.parent_dir):
        print(f"Error: Parent directory not found: {args.parent_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Process single progeny file or all files
    if args.single_progeny:
        # Process single progeny file
        parent_files = get_parent_files(args.parent_dir)
        single_args = (args.single_progeny, args.progeny_dir, args.parent_dir, args.output_dir, parent_files)
        result = process_single_progeny(single_args)
        print(f"Single progeny processing complete: {result}")
    else:
        # Process all progeny files
        try:
            process_all_intersections(
                args.progeny_dir,
                args.parent_dir, 
                args.output_dir,
                num_processes=args.processes,
                progeny_pattern=args.progeny_pattern
            )
            print("All intersections completed successfully")
            
        except Exception as e:
            print(f"Error during analysis: {e}", file=sys.stderr)
            sys.exit(1)

if __name__ == "__main__":
    main()
