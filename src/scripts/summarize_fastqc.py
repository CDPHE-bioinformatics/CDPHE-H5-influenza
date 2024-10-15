#!/usr/bin/env python3

# summarize_fastqc.py

# This script processes FastQC output for raw/cleaned FASTQ files and extracts quality metrics.

import subprocess
import pandas as pd
import argparse

def get_fastqc_metrics(fastqc_file):
    metrics = {}
    grep_commands = {
        "Total Sequences": "r_total_reads",
        "Sequences flagged as poor quality": "flagged_reads_as_poor_quality",
        "Sequence length": "read_len"
    }
    for pattern, metric_name in grep_commands.items():
        command = f'grep "{pattern}" {fastqc_file} | cut -f 2'
        metrics[metric_name] = subprocess.check_output(command, shell=True, text=True).strip()
    return metrics

def summarize_reads(sample_names, fastqc1_data_array, fastqc2_data_array, fastqc_type):
    for i in range(len(sample_names)):
        sample_name = sample_names[i]
        fastqc1_data = fastqc1_data_array[i]
        fastqc2_data = fastqc2_data_array[i]

        print(f"Processing sample: {sample_name}")

        # Extract metrics for R1 and R2
        metrics_r1 = get_fastqc_metrics(fastqc1_data)
        metrics_r2 = get_fastqc_metrics(fastqc2_data)

        # Combine metrics in a dataframe
        df = pd.DataFrame({
            'sample_name': [sample_name],
            'r1_total_reads': [metrics_r1['r_total_reads']],
            'r1_flagged_reads_as_poor_quality': [metrics_r1['flagged_reads_as_poor_quality']],
            'r1_read_len': [metrics_r1['read_len']],
            'r2_total_reads': [metrics_r2['r_total_reads']],
            'r2_flagged_reads_as_poor_quality': [metrics_r2['flagged_reads_as_poor_quality']],
            'r2_read_len': [metrics_r2['read_len']]
        })

        # Save metrics to a TSV file

        outfile = f'{sample_name}_{fastqc_type}_summary_metrics.tsv'
        df.to_csv(outfile, sep='\t', index=False)
        print(f"Summary metrics written to {outfile}")

def main():
    parser = argparse.ArgumentParser(description="Summarize FastQC reads for samples.")
    parser.add_argument("--sample_names", nargs="+", help="List of sample names")
    parser.add_argument("--fastqc1_data_array", nargs="+", help="List of FastQC data files for R1")
    parser.add_argument("--fastqc2_data_array", nargs="+", help="List of FastQC data files for R2")
    parser.add_argument("--fastqc_type", help="type of fastqc file, clean/raw")

    args = parser.parse_args()
    summarize_reads(args.sample_names, args.fastqc1_data_array, args.fastqc2_data_array, args.fastqc_type)

if __name__ == "__main__":
    main()
