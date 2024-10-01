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
    reference = f'{consensus_path_split[1]}: {consensus_path_split[2]}' if len(consensus_path_split) == 4 else consensus_path_split[1]

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
    fname = f'{consensus_path_split[-1].removesuffix("_consensus.fa")}_percent_coverage.csv'
    df.to_csv(fname, index=False)

    return df

# Function to summarize coverage
def summarize_coverage(sample, stats_dir, wd):
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

    # Save output as CSV
    outfile = f'{sample}_{ref_str}_aln_metrics_summary.csv'
    df.to_csv(outfile, index=False)

    return df


# Function to calculate reference lengths with and without primers
def calculate_reference_lengths(reference_sequence_paths, bed_df):
    reference_lengths = {}
    for ref_path in reference_sequence_paths:
        records = list(SeqIO.parse(ref_path, 'fasta'))
        num_records = len(records)
        if num_records > 1:
            segment_lengths = {}
            for segment in records:
                total_length = len(segment.seq)
                segment_bed_df = bed_df[bed_df[0] == segment.id]

                # Calculate primer insert length
                primer_insert_length = segment_bed_df[1].max() - segment_bed_df[2].min()
                segment_lengths[segment.id] = (total_length, primer_insert_length)
            reference_lengths[ref_path] = segment_lengths
        else:
            record = records[0]
            total_length = len(record.seq)
            segment_bed_df = bed_df[bed_df[0] == record.id]
            if segment_bed_df.empty:
                raise ValueError(f'Reference {record.id} not found in BED')
            primer_insert_length = segment_bed_df[1].max() - segment_bed_df[2].min()
            reference_lengths[ref_path] = (total_length, primer_insert_length)

    return reference_lengths

# Function to process segments for a given sample and reference
def process_segments(sample, ref, segments, args):
    percent_coverage_dfs = []
    if isinstance(segments, dict):
        for segment in segments:
            consensus = os.path.join(args.consensus_dir, ref, segment, f'{sample}_{ref}_{segment}_consensus.fa')
            percent_coverage_df = calculate_percent_coverage(
                sample,
                consensus,
                segments[segment][0],
                segments[segment][1]
            )
            percent_coverage_dfs.append(percent_coverage_df)
    else:
        consensus = os.path.join(args.consensus_dir, ref, f'{sample}_{ref}_consensus.fa')
        percent_coverage_df = calculate_percent_coverage(
            sample,
            consensus,
            segments[0],
            segments[1]
        )
        percent_coverage_dfs.append(percent_coverage_df)
    
    return percent_coverage_dfs

# Function to gather percent coverage dataframes for all samples and references
def gather_percent_coverage(args, reference_lengths):
    percent_coverage_dfs = [
        df for sample in args.sample_names for ref, segments in reference_lengths.items()
        for df in process_segments(sample, ref, segments, args)
    ]
    return pd.concat(percent_coverage_dfs).reset_index(drop=True)

# Function to gather coverage summaries for all samples
def gather_coverage_summaries(args, reference_sequence_paths):
    coverage_summary_dfs = [
        summarize_coverage(sample, os.path.dirname(ref_path), os.getcwd())
        for sample in args.sample_names for ref_path in reference_sequence_paths
    ]
    return pd.concat(coverage_summary_dfs).reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(description='Calculate alignment metrics and percent coverage.')
    parser.add_argument('sample_names', nargs='+', help='List of sample names')
    parser.add_argument('consensus_dir', help='Directory containing consensus files')
    parser.add_argument('project_name', help='Name of the project')
    parser.add_argument('reference_dir', help='Directory containing reference sequence files')
    parser.add_argument('bed_file', help='Path to BED file with primer information')

    args = parser.parse_args()

    # Load BED file into DataFrame
    bed_df = pd.read_csv(args.bed_file, sep='\t', header=None)

    # Load reference lengths using the BED file and reference sequences
    reference_sequence_paths = [os.path.join(args.reference_dir, f) for f in os.listdir(args.reference_dir) if f.endswith('.fasta')]
    reference_lengths = calculate_reference_lengths(reference_sequence_paths, bed_df)

    # Calculate and save percent coverage
    concat_percent_coverage_df = gather_percent_coverage(args, reference_lengths)
    concat_percent_coverage_outfile = f'{args.project_name}_percent_coverage.csv'
    concat_percent_coverage_df.to_csv(concat_percent_coverage_outfile, index=False)

    # Calculate and save coverage summaries
    concat_coverage_summary_df = gather_coverage_summaries(args, reference_sequence_paths)
    concat_coverage_summary_outfile = f'{args.project_name}_aln_metrics_summary.csv'
    concat_coverage_summary_df.to_csv(concat_coverage_summary_outfile, index=False)

if __name__ == "__main__":
    main()
