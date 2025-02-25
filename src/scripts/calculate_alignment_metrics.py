#!/usr/bin/env python

import sys
import re
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(description='Calculate alignment metrics and percent coverage.')
    parser.add_argument('--sample_names', nargs='+', help='List of sample names')
    parser.add_argument('--consensus_fastas', nargs='+', help='List of consensus FASTA files')
    parser.add_argument('--samtools_coverages', nargs='+', help='List of samtools coverage.txt files')
    parser.add_argument('--mapped_reads', nargs='+', help='List of samtools mapped reads strings')
    parser.add_argument('--reads_qc_summary', help='Concatenated fastQC summary file')
    parser.add_argument('--project_name', help='Name of the project')
    parser.add_argument('--reference_fasta', help='Reference FASTA file')
    parser.add_argument('--reference_name', help='Reference FASTA name')
    parser.add_argument('--primer_bed', help='BED file with primer information')

    return parser.parse_args()


def calculate_reference_lengths(reference_sequence_path, bed_df):
    """Calculate reference segment lengths with and without primers"""

    records = list(SeqIO.parse(reference_sequence_path, 'fasta'))
    reference_lengths  = {}
    for segment in records:
        total_length = len(segment.seq)
        segment_bed_df = bed_df[bed_df[0] == segment.id]
        if segment_bed_df.empty:
            raise ValueError(f'Reference {segment.id} not found in BED')
        primer_insert_length = segment_bed_df[1].max() - segment_bed_df[2].min()
        reference_lengths [segment.id] = [total_length, primer_insert_length]
    
    total_length = sum(i for i, j in reference_lengths.values())
    primer_insert_length = sum(j for i, j in reference_lengths.values())

    return total_length, primer_insert_length 


def parse_consensus(row):
    """Parse consensus fasta for various metrics."""

    record = SeqIO.read(row.consensus, 'fasta')
    seq_len = len(record.seq)
    number_ns = record.seq.count('N')

    # Count N's in front and back of sequence
    seq_str = str(record.seq)
    front_ns = len(re.findall('^N+', seq_str)[0]) if re.search('^N+', seq_str) else 0
    back_ns = len(re.findall('N+$', seq_str)[0]) if re.search('N+$', seq_str) else 0

    return [seq_len, number_ns, front_ns, back_ns]


# Function to gather percent coverage dataframes for all samples and references
def calculate_percent_coverage(args, total_reference_length, primer_insert_length):
    """Calculate percent coverages and other metrics from consensus fastas."""

    df = pd.DataFrame(data = {'sample_name': args.sample_names, 
                              'consensus': args.consensus_fastas, 
                              'reference': args.reference_name})
    df[['aligned_bases', 'ambiguous_bases', 'front_ns', 'back_ns', ]] = \
        df.apply(lambda x: parse_consensus(x), axis = 1, result_type = 'expand')
    df['aligned_nonambiguous_bases'] = df.aligned_bases - df.ambiguous_bases
    df['percent_coverage_total_length'] = round(df.aligned_nonambiguous_bases / total_reference_length * 100, 1)
    df['percent_coverage_primer_insert_length'] = round(df.aligned_nonambiguous_bases / primer_insert_length * 100, 1)

    df_qc_reads = pd.read_csv(args.reads_qc_summary,
                              usecols = ['sample_name', 
                                         'project_name', 
                                         'primer_name']) 
    df = df.merge(df_qc_reads, on = 'sample_name')

    cols = ['sample_name', 
            'project_name', 
            'primer_name',
            'reference', 
            'percent_coverage_total_length',
            'percent_coverage_primer_insert_length',
            'aligned_bases',
            'ambiguous_bases',
            'aligned_nonambiguous_bases',
            'front_ns',
            'back_ns']
    
    return df[cols]


def parse_samtools_coverages(fn):
    """Return mean depth and mean map quality, weighted by segment length."""
    df = pd.read_csv(fn, sep = '\t')
    
    if df.shape[0] > 1:
        df['len_x_meandepth'] = (df.endpos + 1) * df.meandepth
        total_len = df.endpos.sum() + df.shape[0]
        total_len_x_meandepth = df.len_x_meandepth.sum()
        weighted_meandepth = round(total_len_x_meandepth / total_len, 2)

        df['len_x_meanmapq'] = (df.endpos + 1) * df.meanmapq

        total_len_x_meanmapq = df.len_x_meanmapq.sum()
        weighted_meanmapq = round(total_len_x_meanmapq / total_len, 2)

        return weighted_meandepth, weighted_meanmapq
    else:
        return df.meandepth[0], df.meanmapq[0]
    

def summarize_coverage(args):
    """Read in information from QC reads file and samtools to summarize percent coverage."""
    df = pd.DataFrame(data = {'sample_name': args.sample_names, 
                              'bam_coverage_file': args.samtools_coverages, 
                              'reads_mapped': [int(x) for x in args.mapped_reads]})
    df[['mean_depth', 'mean_mapq']] = df.apply(lambda x: parse_samtools_coverages(fn = x.bam_coverage_file),
                                               result_type = 'expand',
                                               axis = 1)
    
    df_qc_reads = pd.read_csv(args.reads_qc_summary,
                              usecols = ['sample_name', 
                                         'project_name', 
                                         'primer_name',
                                         'r1_total_reads_filtered', 
                                         'r2_total_reads_filtered']) 
    df = df.merge(df_qc_reads, on = 'sample_name')
    
    df['total_filtered_reads'] = df.r1_total_reads_filtered + df.r2_total_reads_filtered
    df['percent_reads_mapped'] = np.where(df.total_filtered_reads != 0, 
                                          round(df.reads_mapped / df.total_filtered_reads * 100, 1), 
                                          0)
    df['reads_unmapped'] = df.total_filtered_reads - df.reads_mapped

    cols = ['sample_name',
            'project_name', 
            'primer_name',
            'percent_reads_mapped',
            'total_filtered_reads',
            'reads_mapped',
            'reads_unmapped',
            'mean_depth',
            'mean_mapq'
            ]

    return df[cols]


def main():
    args = get_args()

    bed_df = pd.read_csv(args.primer_bed, sep='\t', header=None)

    # Load reference lengths using the BED file and reference sequences
    total_reference_length, primer_insert_length = calculate_reference_lengths(args.reference_fasta, bed_df)

    # Calculate and save percent coverage
    concat_percent_coverage_df = calculate_percent_coverage(args, total_reference_length, primer_insert_length)
    concat_percent_coverage_df.to_csv(f'{args.project_name}_percent_coverage.csv', index=False)

    # Calculate and save coverage summaries
    concat_coverage_summary_df = summarize_coverage(args)
    concat_coverage_summary_df.to_csv(f'{args.project_name}_aligned_metrics_summary.csv', index=False)

if __name__ == "__main__":
    main()