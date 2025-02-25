#!/usr/bin/env python

import argparse
import re
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description='Sort bed file for samtools ampliconclip')
    parser.add_argument('--bed_file', help = 'Primer bed file', required = True)
    parser.add_argument('--out_fn', help = 'Sorted bed file name', required = True)

    return parser.parse_args()


def sort_bed(bed_file):
    """Sort bed file"""
    df = pd.read_csv(bed_file, sep = '\t', header = None)
    df['amplicon'] = df[3].str.replace('[_-](LEFT|RIGHT)', '', regex = True)
    df['amplicon_number'] = df.amplicon.str.extract(r'(\d+)(?!.*\d)', expand = False).astype(int)

    df = (df
          .sort_values([0, 'amplicon_number', 5])
          .drop(columns = ['amplicon', 'amplicon_number'])
    )
    
    return df


def main():
    args = get_args()
    df = sort_bed(args.bed_file)
    df.to_csv(args.out_fn, sep = '\t', index = False, header = False)

if __name__ == "__main__":
    main()