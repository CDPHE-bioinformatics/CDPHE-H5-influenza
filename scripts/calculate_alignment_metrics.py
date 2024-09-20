#!/usr/bin/env python

# calculate_alignment_metrics.py

# This script calculates alignment metrics and percent coverage.
import os
import re
import pandas as pd
from Bio import SeqIO
import argparse
import json

# Function to calculate percent coverage
def calculate_percent_coverage(sample, consensus, total_reference_length, primer_insert_length, aln_metrics_dir):
    record = SeqIO.read(consensus, 'fasta')
    seq_len = len(record.seq)
    number_ns = record.seq.count('N')
    total_seq_len = seq_len - number_ns

    # Count N's in front and back of sequence
    seq_str = str(record.seq)
    front_ns = len(re.findall('^N+', seq_str)[0]) if re.search('^N+', seq_str) else 0
    back_ns = len(re.findall('N+$', seq_str)[0]) if re.search('N+$', seq_str) else 0

    percent_coverage_total_length = round(total_seq_len / total_reference_length * 100, 1)
    percent_coverage_primer_insert_length = round(total_seq_len / primer_insert_length * 100, 1)

    # Generate reference and subdir paths
    consensus_path_split = consensus.split('/')
    if len(consensus_path_split) == 4:
        reference = f'{consensus_path_split[1]}: {consensus_path_split[2]}'
        subdir = os.path.join(consensus_path_split[1], consensus_path_split[2])
    else:
        reference = consensus_path_split[1]
        subdir = consensus_path_split[1]

    # Create DataFrame for percent coverage results
    df = pd.DataFrame({
        'sample_name': [sample],
        'reference': [reference],
        'percent_coverage_total_length': [percent_coverage_total_length],
        'percent_coverage_primer_insert_length': [percent_coverage_primer_insert_length],
        'number_aln_bases': [seq_len],
        'ambiguous_bases': [number_ns],
        'number_aln_nonambiguous_bases': [total_seq_len],
        'front_ns': [front_ns],
        'back_ns': [back_ns]
    })

    # Save output to CSV
    os.makedirs(os.path.join(aln_metrics_dir, subdir), exist_ok=True)
    fname = f'{consensus_path_split[-1].removesuffix("_consensus.fa")}_percent_coverage.csv'
    outfile = os.path.join(aln_metrics_dir, subdir, fname)
    df.to_csv(outfile, index=False)

    return df

# Function to summarize coverage
def summarize_coverage(sample, stats_dir, wd, aln_metrics_dir):
    stats_dir_split = stats_dir.split('/')
    ref_str = f'{stats_dir_split[1]}_{stats_dir_split[2]}' if len(stats_dir_split) == 3 else f'{stats_dir_split[1]}'

    bam_coverage = os.path.join(stats_dir, f'{sample}_{ref_str}_coverage.txt')
    bam_cov_df = pd.read_csv(bam_coverage, sep='\t')

    # Parse BAM stats for reads mapped
    bam_stats = os.path.join(stats_dir, f'{sample}_{ref_str}_stats.txt')
    reads_mapped = 0
    with open(bam_stats, 'r') as file:
        for line in file:
            if re.search('reads mapped:', line):
                reads_mapped = float(line.split('\t')[2].strip())

    # Retrieve total raw reads from FastQC
    fastqc_metrics_file = os.path.join(wd, "fastqc_clean", sample, f'{sample}_summary_metrics.tsv')
    fastqc_df = pd.read_csv(fastqc_metrics_file, sep='\t')
    total_raw_reads = fastqc_df['r1_total_reads'][0] + fastqc_df['r2_total_reads'][0]

    # Calculate reads unmapped and percent mapped
    percent_reads_mapped = round(reads_mapped / total_raw_reads * 100, 1) if total_raw_reads != 0 else 0
    reads_unmapped = total_raw_reads - reads_mapped

    # Retrieve percent coverage data
    percent_coverage_file = os.path.join(stats_dir, f'{sample}_{ref_str}_percent_coverage.csv')
    percent_df = pd.read_csv(percent_coverage_file)

    # Create DataFrame for summarized coverage
    df = pd.DataFrame({
        'sample_name': [sample],
        'reference': [percent_df['reference'].iloc[0]],
        'percent_reads_mapped': [percent_reads_mapped],
        'total_raw_reads': [total_raw_reads],
        'reads_mapped': [reads_mapped],
        'reads_unmapped': [reads_unmapped],
        'av_depth': [bam_cov_df['meandepth'].iloc[0]],
        'meanmapq': [bam_cov_df['meanmapq'].iloc[0]],
        'percent_coverage_total_length': [percent_df['percent_coverage_total_length'].iloc[0]],
        'percent_coverage_primer_insert_length': [percent_df['percent_coverage_primer_insert_length'].iloc[0]],
        'front_ns': [percent_df['front_ns'].iloc[0]],
        'back_ns': [percent_df['back_ns'].iloc[0]]
    })

    # Save output to CSV
    outfile = os.path.join(stats_dir, f'{sample}_{ref_str}_aln_metrics_summary.csv')
    df.to_csv(outfile, index=False)

    return df


def main():
    parser = argparse.ArgumentParser(description='Calculate alignment metrics and percent coverage.')
    parser.add_argument('sample_names', nargs='+', help='List of sample names')
    parser.add_argument('reference_lengths', help='Path to a file with reference lengths')
    parser.add_argument('consensus_dir', help='Directory containing consensus files')
    parser.add_argument('aln_metrics_dir', help='Directory to save alignment metrics')
    parser.add_argument('project_name', help='Name of the project')
    parser.add_argument('reference_sequence_paths', help='Path to reference sequences')

    args = parser.parse_args()

    # sample_names = ['2407110180', 'SG20240710']
    # consensus_dir = /home/irin_paul/h5/test_h5_0005_nextseq/avrl_h5n1_250bp/consensus_sequences
    # aln_metrics_dir = /home/irin_paul/h5/test_h5_0005_nextseq/avrl_h5n1_250bp/metrics
    # project_name = avrl_h5n1

    # Load reference lengths and sequence paths
    with open(args.reference_lengths) as f:
        reference_lengths = json.load(f)
    
    reference_sequence_paths = args.reference_sequence_paths.split(',')

    percent_coverage_dfs = []
    for sample in args.sample_names:
        for ref, segments in reference_lengths.items():
            if isinstance(segments, dict):
                for segment in segments:
                    consensus = os.path.join(args.consensus_dir, ref, segment, f'{sample}{ref}{segment}_consensus.fa')
                    percent_coverage_df = calculate_percent_coverage(
                        sample,
                        consensus,
                        segments[segment][0],
                        segments[segment][1],
                        args.aln_metrics_dir
                    )
                    percent_coverage_dfs.append(percent_coverage_df)
            else:
                consensus = os.path.join(args.consensus_dir, ref, f'{sample}_{ref}_consensus.fa')
                percent_coverage_df = calculate_percent_coverage(
                    sample,
                    consensus,
                    segments[0],
                    segments[1],
                    args.aln_metrics_dir
                )
                percent_coverage_dfs.append(percent_coverage_df)

    # Concatenate percent coverage dataframes and output the summary
    concat_percent_coverage_df = pd.concat(percent_coverage_dfs).reset_index(drop=True)
    concat_percent_coverage_outfile = os.path.join(args.aln_metrics_dir, f'{args.project_name}_percent_coverage.csv')
    concat_percent_coverage_df.to_csv(concat_percent_coverage_outfile, index=False)

    # Summarize coverage for each sample
    coverage_summary_dfs = []
    for sample in args.sample_names:
        for ref_path in reference_sequence_paths:
            records = SeqIO.parse(ref_path, format='fasta')
            num_records = len(list(records))
            if num_records > 1:
                for segment in SeqIO.parse(ref_path, format='fasta'):
                    stats_dir = os.path.join(args.aln_metrics_dir, ref_path)
                    coverage_summary_df = summarize_coverage(sample, stats_dir, args.aln_metrics_dir, args.aln_metrics_dir)
                    coverage_summary_dfs.append(coverage_summary_df)
            else:
                stats_dir = os.path.join(args.aln_metrics_dir, ref_path)
                coverage_summary_df = summarize_coverage(sample, stats_dir, args.aln_metrics_dir, args.aln_metrics_dir)
                coverage_summary_dfs.append(coverage_summary_df)

    # Concatenate coverage summary dataframes and output the summary
    concat_coverage_summary_df = pd.concat(coverage_summary_dfs).reset_index(drop=True)
    concat_coverage_summary_outfile = os.path.join(args.aln_metrics_dir, f'{args.project_name}_aln_metrics_summary.csv')
    concat_coverage_summary_df.to_csv(concat_coverage_summary_outfile, index=False)


if __name__ == "__main__":
    main()
