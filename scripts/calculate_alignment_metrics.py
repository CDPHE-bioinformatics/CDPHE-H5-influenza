#!/usr/bin/env python

# calculate_alignment_metrics.py

# This script calculates alignment metrics and percent coverage.
import os
import re
import pandas as pd
from Bio import SeqIO
import argparse

# Function to calculate percent coverage
def calculate_percent_coverage(sample, consensus, total_reference_length, primer_insert_length):
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

    # Generate reference from the consensus path
    consensus_path_split = consensus.split('/')
    reference = f'{consensus_path_split[-2]}: {consensus_path_split[-1]}' if len(consensus_path_split) ==4 else consensus_path_split[0]

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

    # Save output as CSV
    fname = f'{sample}_percent_coverage.csv'
    df.to_csv(fname, index=False)

    return df

# Function to summarize coverage
def summarize_coverage(sample, coverage_file, stats_file, fastqc_metrics_file):
    # Read coverage file
    bam_cov_df = pd.read_csv(coverage_file, sep='\t')

    # Parse stats for reads mapped
    reads_mapped = 0
    with open(stats_file, 'r') as file:
        for line in file:
            if re.search('reads mapped:', line):
                reads_mapped = float(line.split('\t')[2].strip())

    # Retrieve total raw reads from FastQC
    fastqc_df = pd.read_csv(fastqc_metrics_file, sep='\t')
    total_raw_reads = fastqc_df['r1_total_reads'][0] + fastqc_df['r2_total_reads'][0]

    # Calculate reads unmapped and percent mapped
    percent_reads_mapped = round(reads_mapped / total_raw_reads * 100, 1) if total_raw_reads != 0 else 0
    reads_unmapped = total_raw_reads - reads_mapped

    # Create DataFrame for summarized coverage
    df = pd.DataFrame({
        'sample_name': [sample],
        'percent_reads_mapped': [percent_reads_mapped],
        'total_raw_reads': [total_raw_reads],
        'reads_mapped': [reads_mapped],
        'reads_unmapped': [reads_unmapped],
        'av_depth': [bam_cov_df['meandepth'].iloc[0]],
        'meanmapq': [bam_cov_df['meanmapq'].iloc[0]]
    })

    # Save output as CSV
    outfile = f'{sample}_aln_metrics_summary.csv'
    df.to_csv(outfile, index=False)

    return df

# Function to calculate reference lengths with and without primers
def calculate_reference_lengths(reference_sequence_path, bed_df):
    records = list(SeqIO.parse(reference_sequence_path, 'fasta'))
    total_length = len(records[0].seq)
    segment_bed_df = bed_df[bed_df[0] == records[0].id]
    
    if segment_bed_df.empty:
        raise ValueError(f'Reference {records[0].id} not found in BED')
        
    primer_insert_length = segment_bed_df[1].max() - segment_bed_df[2].min()
    return total_length, primer_insert_length

# Function to gather percent coverage dataframes for all samples and references
def gather_percent_coverage(args, reference_lengths):
    percent_coverage_dfs = []
    for sample in args.sample_names:
        for i, consensus in enumerate(args.consensus_fastas):
            total_length, primer_insert_length = reference_lengths[i]
            percent_coverage_df = calculate_percent_coverage(sample, consensus, total_length, primer_insert_length)
            percent_coverage_dfs.append(percent_coverage_df)
    return pd.concat(percent_coverage_dfs).reset_index(drop=True)

# Function to gather coverage summaries for all samples
def gather_coverage_summaries(args):
    coverage_summary_dfs = []
    for sample in args.sample_names:
        bam_coverage_file = f'{sample}_coverage.txt'
        bam_stats_file = f'{sample}_stats.txt'
        fastqc_metrics_file = f'{sample}_clean_summary_metrics.tsv'
        
        coverage_summary_df = summarize_coverage(sample, bam_coverage_file, bam_stats_file, fastqc_metrics_file)
        coverage_summary_dfs.append(coverage_summary_df)
    return pd.concat(coverage_summary_dfs).reset_index(drop=True)

def main():
    parser = argparse.ArgumentParser(description='Calculate alignment metrics and percent coverage.')
    parser.add_argument('sample_names', nargs='+', help='List of sample names')
    parser.add_argument('consensus_fastas', nargs='+', help='List of consensus FASTA files')
    parser.add_argument('project_name', help='Name of the project')
    parser.add_argument('reference_fasta', help='reference FASTA file')
    parser.add_argument('primer_bed', help='BED file with primer information')

    args = parser.parse_args()

    # Load BED file into DataFrame
    bed_df = pd.read_csv(args.primer_bed, sep='\t', header=None)

    # Load reference lengths using the BED file and reference sequences
    reference_lengths = calculate_reference_lengths(args.reference_fasta, bed_df)

    # Calculate and save percent coverage
    concat_percent_coverage_df = gather_percent_coverage(args, [reference_lengths])
    concat_percent_coverage_outfile = f'{args.project_name}_percent_coverage.csv'
    concat_percent_coverage_df.to_csv(concat_percent_coverage_outfile, index=False)

    # Calculate and save coverage summaries
    concat_coverage_summary_df = gather_coverage_summaries(args)
    concat_coverage_summary_outfile = f'{args.project_name}_aln_metrics_summary.csv'
    concat_coverage_summary_df.to_csv(concat_coverage_summary_outfile, index=False)

if __name__ == "__main__":
    main()