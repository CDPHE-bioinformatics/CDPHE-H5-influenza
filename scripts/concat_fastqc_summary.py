#!/usr/bin/env python

# concat_fastqc_summary.py

# This script combines the FastQC metrics from raw and cleaned fastq files and calculates additional metrics.

import os
import re
import pandas as pd
import argparse

def get_sample_name(filename):
    sample_name = re.sub(r'_clean_summary_metrics.tsv$|_raw_summary_metrics.tsv$', '', filename)
    return sample_name

def combine_fastqc_summaries(summarized_fastqcs, project_name):
    df_raw_list = []
    df_filtered_list = []

    for file in summarized_fastqcs:
        sample_name = get_sample_name(os.path.basename(file))

        # read clean data
        if file.endswith('_clean_summary_metrics.tsv'):
            df_filtered = pd.read_csv(file, sep = '\t')
            df_filtered['sample_name'] = sample_name
            df_filtered_list.append(df_filtered)

        elif file.endswith('_raw_summary_metrics.tsv'):
            df_raw = pd.read_csv(file, sep = '\t')
            df_raw['sample_name'] = sample_name
            df_raw_list.append(df_raw)
        else:
            print(f"Skipping unrecognized file: {file}")

    
    # concatenate dataframes and join
    raw_df = pd.concat(df_raw_list).reset_index(drop = True)
    raw_df = raw_df.set_index('sample_name')

    filtered_df = pd.concat(df_filtered_list).reset_index(drop = True)
    filtered_df = filtered_df.set_index('sample_name')


    combined_df = raw_df.join(filtered_df, how = 'left', lsuffix = '_raw', rsuffix = '_filtered')
    combined_df = combined_df.reset_index()

    # calculate total reads diff
    combined_df['total_reads_diff'] = combined_df['r1_total_reads_raw'] - combined_df['r1_total_reads_filtered']
    combined_df['raw_total_reads_paired'] = combined_df['r1_total_reads_raw']
    combined_df['filtered_total_reads_paired'] = combined_df['r1_total_reads_filtered']

    # add primer set
    for row in range(combined_df.shape[0]):
        sample = combined_df.sample_name[row]

    # print(j.columns)

    col_order = ['sample_name', 'raw_total_reads_paired', 'filtered_total_reads_paired',
                'total_reads_diff', 'r1_total_reads_raw', 'r1_flagged_reads_as_poor_quality_raw',
                'r1_read_len_raw', 'r2_total_reads_raw', 'r2_flagged_reads_as_poor_quality_raw',
                'r2_read_len_raw', 'r1_total_reads_filtered', 'r1_flagged_reads_as_poor_quality_filtered',
                'r1_read_len_filtered', 'r2_total_reads_filtered',
                'r2_flagged_reads_as_poor_quality_filtered', 'r2_read_len_filtered']

    combined_df = combined_df[col_order]
    combined_df['project_name'] = project_name

    # Save to CSV
    outfile = f'{project_name}_reads_QC_summary.csv'
    combined_df.to_csv(outfile, index=False)
    print(f"Summary metrics are stored in {outfile}")

    return combined_df

def main():
    # argument parser
    parser = argparse.ArgumentParser(description="Combine FastQC summary for raw and filtered data.")
    parser.add_argument("--summarized_fastqcs", nargs='+', help="Array of Fastqc summary files")
    parser.add_argument("--project_name", help="Name of the project")

    args = parser.parse_args()

    # Call the function
    combined_df = combine_fastqc_summaries(args.summarized_fastqcs, args.project_name)

if __name__ == "__main__":
    main()
    