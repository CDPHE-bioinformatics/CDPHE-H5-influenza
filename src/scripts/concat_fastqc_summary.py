#!/usr/bin/env python

# concat_fastqc_summary.py

# This script combines the FastQC metrics from raw and cleaned fastq files and calculates additional metrics.

import sys
import argparse
import pandas as pd


def get_args(argv = sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Combine FastQC summary for raw and filtered data.")
    parser.add_argument("--raw_summarized_fastqcs", nargs='+', help="Array of raw fastqc summary files")
    parser.add_argument("--clean_summarized_fastqcs", nargs='+', help="Array of clean fastqc summary files")
    parser.add_argument("--out_fn", help="Name for output file")

    args = parser.parse_args(argv)
    
    return args


def combine_fastqc_summaries(raw_summarized_fastqcs,
                             clean_summarized_fastqcs,
                             out_fn):
    
    clean_dfs = []
    raw_dfs = []
    
    for file in clean_summarized_fastqcs:
        clean_dfs.append(pd.read_csv(file))
    for file in raw_summarized_fastqcs:
        raw_dfs.append(pd.read_csv(file))

    df_clean = pd.concat(clean_dfs, ignore_index = True)
    df_raw = pd.concat(raw_dfs, ignore_index = True)

    df = df_raw.merge(df_clean, 
                      on = ['sample_name', 'project_name', 'primer_name'], 
                      how = 'outer', 
                      suffixes = ('_raw', '_filtered'))

    # calculate total reads diff
    df['total_reads_diff'] = df['r1_total_reads_raw'] - df['r1_total_reads_filtered']
    df['raw_total_reads_paired'] = df['r1_total_reads_raw']
    df['filtered_total_reads_paired'] = df['r1_total_reads_filtered']

    col_order = ['sample_name', 'project_name', 'primer_name', 
                 'raw_total_reads_paired', 'filtered_total_reads_paired',
                 'total_reads_diff', 'r1_total_reads_raw', 'r1_flagged_reads_as_poor_quality_raw',
                 'r1_read_len_raw', 'r2_total_reads_raw', 'r2_flagged_reads_as_poor_quality_raw',
                 'r2_read_len_raw', 'r1_total_reads_filtered', 'r1_flagged_reads_as_poor_quality_filtered',
                 'r1_read_len_filtered', 'r2_total_reads_filtered',
                 'r2_flagged_reads_as_poor_quality_filtered', 'r2_read_len_filtered']

    df = df[col_order]

    df.to_csv(out_fn, index=False)
    print(f"Summary metrics are stored in {out_fn}")

    return 


def main():
    args = get_args()
    combine_fastqc_summaries(args.raw_summarized_fastqcs, 
                             args.clean_summarized_fastqcs, 
                             args.out_fn)

if __name__ == "__main__":
    main()
    