#!/usr/bin/env python

# get_summary_metrics.py

# This script combines the FastQC metrics from raw and cleaned fastq files and calculates additional metrics.

import os
import pandas as pd
import argparse

def combine_fastqc_summaries(sample_names, wd, project_name):
    df_raw_list = []
    df_filtered_list = []

    # Loop through each sample and read raw and filtered data
    for sample in sample_names:
        print(f"Processing sample: {sample}")
        # Read in raw summary file
        path_raw = os.path.join(wd, 'fastqc_raw', sample, f'{sample}_summary_metrics.tsv')
        df_raw = pd.read_csv(path_raw, sep='\t')
        df_raw_list.append(df_raw)

        # Read in filtered summary file
        path_filtered = os.path.join(wd, 'fastqc_clean', sample, f'{sample}_summary_metrics.tsv')
        df_filtered = pd.read_csv(path_filtered, sep='\t')
        df_filtered_list.append(df_filtered)

    # Concatenate dataframes and set the index
    raw_df = pd.concat(df_raw_list).reset_index(drop=True)
    raw_df = raw_df.set_index('sample_name')

    filtered_df = pd.concat(df_filtered_list).reset_index(drop=True)
    filtered_df = filtered_df.set_index('sample_name')

    # Join the raw and filtered dataframes
    combined_df = raw_df.join(filtered_df, how='left', lsuffix='_raw', rsuffix='_filtered')
    combined_df = combined_df.reset_index()

    # Calculate total reads difference and paired reads
    combined_df['total_reads_diff'] = combined_df['r1_total_reads_raw'] - combined_df['r1_total_reads_filtered']
    combined_df['raw_total_reads_paired'] = combined_df['r1_total_reads_raw']
    combined_df['filtered_total_reads_paired'] = combined_df['r1_total_reads_filtered']

    # Define the order of columns
    col_order = ['sample_name', 'raw_total_reads_paired', 'filtered_total_reads_paired',
                 'total_reads_diff', 'r1_total_reads_raw', 'r1_flagged_reads_as_poor_quality_raw',
                 'r1_read_len_raw', 'r2_total_reads_raw', 'r2_flagged_reads_as_poor_quality_raw',
                 'r2_read_len_raw', 'r1_total_reads_filtered', 'r1_flagged_reads_as_poor_quality_filtered',
                 'r1_read_len_filtered', 'r2_total_reads_filtered',
                 'r2_flagged_reads_as_poor_quality_filtered', 'r2_read_len_filtered']

    combined_df = combined_df[col_order]
    
    # Add project name column
    combined_df['project_name'] = project_name

    # Save to CSV
    outfile = os.path.join(wd, f'{project_name}_reads_QC_summary.csv')
    combined_df.to_csv(outfile, index=False)
    print(f"Summary metrics are stored in {outfile}")

    return combined_df

def main():
    # argument parser
    parser = argparse.ArgumentParser(description="Combine FastQC summary for raw and filtered data.")
    parser.add_argument("sample_names", nargs="+", help="List of sample names")
    parser.add_argument("wd", help="Working directory containing FastQC data")
    parser.add_argument("project_name", help="Name of the project")

    args = parser.parse_args()

    # Call the function
    combined_df = combine_fastqc_summaries(args.sample_names, args.wd, args.project_name)

    # example arguments:
    # sample_names = 2407110180 SG20240710
    # wd = /home/irin_paul/h5/test_h5_0005_nextseq/avrl_h5n1_250bp
    # project_name = avrl_h5n1


if __name__ == "__main__":
    main()
    