#!/usr/bin/env python

import sys
import re
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(description='Calculate alignment metrics and percent coverage.')
    parser.add_argument('--consensus_fastas', nargs='+', help='List of consensus FASTA files')
    parser.add_argument('--out_fn_prefix', help = 'Out filename prefix')
    parser.add_argument('--primer_bed', help='BED file with primer information')
    parser.add_argument('--project_name', help='Name of the project')
    parser.add_argument('--reads_qc_summary', help='Concatenated fastQC summary file')
    parser.add_argument('--reference_fasta', help='Reference FASTA file')
    parser.add_argument('--reference_name', help='Reference FASTA name')
    parser.add_argument('--sample_names', nargs='+', help='List of sample names')
    parser.add_argument('--samtools_coverages', nargs='+', help='List of samtools coverage.txt files')
    return parser.parse_args()



def calculate_reference_lengths(reference_sequence_path, bed_df):
    """Calculate reference segment lengths with and without primers"""

    records = list(SeqIO.parse(reference_sequence_path, 'fasta'))
    segment_lengths  = {}
    for segment in records:
        total_length = len(segment.seq)
        segment_bed_df = bed_df[bed_df[0] == segment.id]
        if segment_bed_df.empty:
            raise ValueError(f'Reference {segment.id} not found in BED')
        primer_insert_span = segment_bed_df[1].max() - segment_bed_df[2].min()
        segment_lengths[segment.id] = [total_length, primer_insert_span]

    return segment_lengths 


def parse_consensus(row):
    """Parse consensus fasta for various metrics."""

    records = list(SeqIO.parse(row.consensus, 'fasta'))
    record_ids = [r.id for r in records]
    if len(record_ids) == 1:
        record_idx = 0
    else:
        try:
            record_idx = [row.record_id in r for r in record_ids].index(True)
        except ValueError:
            print(f'{row.record_id} could not be found in the consensus fasta record ids:')
            print(record_ids, sep = '\n')
            sys.exit(1)
    
    record = records[record_idx]
    seq_len = len(record.seq)
    number_ns = record.seq.count('N')

    seq_str = str(record.seq)
    front_ns = len(seq_str) - len(seq_str.lstrip('N'))
    back_ns = len(seq_str) - len(seq_str.rstrip('N'))

    return [seq_len, number_ns, front_ns, back_ns]


def calculate_segment_stats(args, segment_lengths):
    """Calculate percent coverages and other metrics for each record in sample fasta."""

    df = pd.DataFrame(data = {'sample_name': args.sample_names, 
                              'consensus': args.consensus_fastas, 
                              'reference': args.reference_name,
                              'record_id': [segment_lengths.keys()] * len(args.sample_names),
                              'record_lengths': [segment_lengths.values()] * len(args.sample_names)})
    df = df.explode(['record_id', 'record_lengths'])
    df['record_length'], df['record_primer_span'] = zip(*list(df['record_lengths'].values))

    df[['aligned_bases', 'ambiguous_bases', 'front_ns', 'back_ns', ]] = \
        df.apply(lambda x: parse_consensus(x), axis = 1, result_type = 'expand')
    df['aligned_nonambiguous_bases'] = df.aligned_bases - df.ambiguous_bases
    df['percent_coverage_total_length'] = round(df.aligned_nonambiguous_bases / df.record_length * 100, 1)
    df['percent_coverage_primer_span'] = round(df.aligned_nonambiguous_bases / df.record_primer_span * 100, 1)

    df_bam_cov = parse_samtools_coverages(args.sample_names, args.samtools_coverages)
    df = df.merge(df_bam_cov, on = ['sample_name', 'record_id'], how = 'outer')

    df_qc_reads = pd.read_csv(args.reads_qc_summary,
                              usecols = ['sample_name', 
                                         'project_name', 
                                         'primer_name']) 
    df = df.merge(df_qc_reads, on = 'sample_name')

    cols = ['sample_name', 
            'project_name', 
            'primer_name',
            'reference', 
            'record_id',
            'record_length',
            'record_primer_span',
            'percent_coverage_total_length',
            'percent_coverage_primer_span',
            'aligned_bases',
            'ambiguous_bases',
            'aligned_nonambiguous_bases',
            'front_ns',
            'back_ns',
            'num_reads',
            'num_covered_bases',
            'samtools_coverage',
            'mean_depth',
            'mean_base_quality',
            'mean_map_quality']
    
    return df[cols]


def average_samtools_coverages(df):
    """Return mean percent coverage, depth, and map quality, weighted by segment length."""
    
    if df.shape[0] > 1:
        total_len = df.record_length.sum()

        df['len_x_percent_cov'] = df.record_length * df.percent_coverage_total_length
        total_len_x_percent_cov = df.len_x_percent_cov.sum()
        weighted_percent_cov = round(total_len_x_percent_cov / total_len, 2)
        
        df['len_x_mean_depth'] = df.record_length * df.mean_depth
        total_len_x_mean_depth = df.len_x_mean_depth.sum()
        weighted_mean_depth = round(total_len_x_mean_depth / total_len, 2)

        df['len_x_mean_map_quality'] = df.record_length * df.mean_map_quality
        total_len_x_mean_map_quality = df.len_x_mean_map_quality.sum()
        weighted_mean_map_quality = round(total_len_x_mean_map_quality / total_len, 2)

        ser = pd.Series(data = {'avg_percent_coverage': weighted_percent_cov,
                                'avg_depth': weighted_mean_depth,
                                'avg_map_quality': weighted_mean_map_quality})
    else:
        ser = (df[['percent_coverage_total_length', 'mean_depth', 'mean_map_quality']]
               .rename(columns = {'percent_coverage_total_length': 'avg_percent_coverage',
                                  'mean_depth': 'avg_depth', 
                                  'mean_map_quality': 'avg_map_quality'})
               .squeeze(axis = 0)
        )
    
    return ser
    

def parse_samtools_coverages(sample_names, fns):
    """Read in samtools coverage files."""
    dfs = []
    for i in range(len(sample_names)):
        df = pd.read_csv(fns[i], sep = '\t')
        df['sample_name'] = sample_names[i]
        dfs.append(df)
    concat_df = (pd
                 .concat(dfs, ignore_index = True)
                 .rename(columns = {'#rname': 'record_id',
                                    'coverage': 'samtools_coverage',
                                    'covbases': 'num_covered_bases',
                                    'numreads': 'num_reads',
                                    'meandepth': 'mean_depth',
                                    'meanbaseq': 'mean_base_quality',
                                    'meanmapq': 'mean_map_quality'})
                 .drop(columns = ['startpos', 'endpos'])
    )

    return concat_df


def calculate_sample_stats(segmented_df, args):
    """Calculate sample level stats."""
    segmented_df = segmented_df[['sample_name', 
                                 'project_name',
                                 'primer_name',
                                 'reference', 
                                 'record_id',
                                 'record_length',
                                 'num_reads',
                                 'percent_coverage_total_length',
                                 'mean_depth',
                                 'mean_map_quality']]

    df_meta = (segmented_df[['sample_name', 'project_name', 'primer_name', 'reference']]
               .copy()
               .drop_duplicates()
    )
    df_sum = (segmented_df[['sample_name', 'num_reads']]
              .copy()
              .groupby(['sample_name'])
              .sum()
              .reset_index(names = ['sample_name'])
              .rename(columns = {'num_reads': 'reads_mapped'})
              )
    df_avg = (segmented_df[['sample_name', 'record_id', 'record_length', 
                            'percent_coverage_total_length', 'mean_depth', 'mean_map_quality']]
              .copy()
              .groupby('sample_name')
              .apply(lambda x: average_samtools_coverages(x))
              .reset_index()
    )

    df_qc_reads = pd.read_csv(args.reads_qc_summary,
                              usecols = ['sample_name',
                                         'r1_total_reads_filtered', 
                                         'r2_total_reads_filtered']) 
    
    df = (df_meta
          .merge(df_sum, on = 'sample_name')
          .merge(df_avg, on = 'sample_name')
          .merge(df_qc_reads, on = 'sample_name')
    )

    df['total_filtered_reads'] = df.r1_total_reads_filtered + df.r2_total_reads_filtered
    df['percent_reads_mapped'] = np.where(df.total_filtered_reads != 0, 
                                          round(df.reads_mapped / df.total_filtered_reads * 100, 1), 
                                          0)
    df['reads_unmapped'] = df.total_filtered_reads - df.reads_mapped

    cols = ['sample_name',
            'project_name', 
            'primer_name',
            'reference',
            'avg_percent_coverage',
            'percent_reads_mapped',
            'total_filtered_reads',
            'reads_mapped',
            'reads_unmapped',
            'avg_depth',
            'avg_map_quality'
            ]

    
    return df[cols]
    

def main():
    args = get_args()

    bed_df = pd.read_csv(args.primer_bed, sep='\t', header=None)

    # Load reference lengths using the BED file and reference sequences
    segment_lengths = calculate_reference_lengths(args.reference_fasta, bed_df)

    # Calcuate segmented metrics
    segmented_df = calculate_segment_stats(args, segment_lengths)
    segmented_df.to_csv(f'{args.out_fn_prefix}_segment_metrics.csv', index=False)

    # Calculate overall sample metrics
    sample_df = calculate_sample_stats(segmented_df, args)
    sample_df.to_csv(f'{args.out_fn_prefix}_sample_metrics.csv', index=False)

if __name__ == "__main__":
    main()