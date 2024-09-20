#!/usr/bin/env python3

# summarize_fastqc.py

# This script processes the FastQC output for raw/cleaned FASTQ files and extracts quality metrics.

import os
import glob
import subprocess
import pandas as pd
import argparse

def summarize_reads(sample_names, fastqc_dir):
    for sample in sample_names:
        print(f"Processing sample: {sample}")
        
        # set outpath

        outpath = os.path.join(fastqc_dir, sample)
        os.makedirs(outpath, exist_ok=True)

        # create dataframe
        df = pd.DataFrame()
        df['sample_name'] = [sample]
        

        # R1
        path = glob.glob(os.path.join(fastqc_dir, f'{sample}*1_fastqc_data.txt'))[0]
        

        command = f'cat {path} | grep "Total Sequences" | cut -f 2'
        total_seqs = subprocess.check_output(command, shell=True, text=True).strip()
        df['r1_total_reads'] = [total_seqs]

        command = f'cat {path} | grep "Sequences flagged as poor quality" | cut -f 2'
        flagged_seqs = subprocess.check_output(command, shell=True, text=True).strip()
        df['r1_flagged_reads_as_poor_quality'] = [flagged_seqs]

        command = f'cat {path} | grep "Sequence length" | cut -f 2'
        seq_len = subprocess.check_output(command, shell=True, text=True).strip()
        df['r1_read_len'] = [seq_len]
        
        # R2
        path = glob.glob(os.path.join(fastqc_dir, f'{sample}*2_fastqc_data.txt'))[0]
        
        command = f'cat {path} | grep "Total Sequences" | cut -f 2'
        total_seqs = subprocess.check_output(command, shell=True, text=True).strip()
        df['r2_total_reads'] = [total_seqs]

        command = f'cat {path} | grep "Sequences flagged as poor quality" | cut -f 2'
        flagged_seqs = subprocess.check_output(command, shell=True, text=True).strip()
        df['r2_flagged_reads_as_poor_quality'] = [flagged_seqs]

        command = f'cat {path} | grep "Sequence length" | cut -f 2'
        seq_len = subprocess.check_output(command, shell=True, text=True).strip()
        df['r2_read_len'] = [seq_len]
        
        # Write summary metrics to file
        outfile = os.path.join(outpath, f'{sample}_summary_metrics.tsv')
        df.to_csv(outfile, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description="Summarize FastQC reads for samples.")
    parser.add_argument("sample_names", nargs="+", help="List of sample names")
    parser.add_argument("fastqc_dir", help="Directory containing unzipped FastQC results")

    args = parser.parse_args()
    
    # Example arguments:
    #sample_names = 2407110180 SG20240710
    #fastqc_dir = /home/irin_paul/h5/test_h5_0005_nextseq/avrl_h5n1_250bp/fastqc_raw
    # OR
    #fastqc_dir = /home/irin_paul/h5/test_h5_0005_nextseq/avrl_h5n1_250bp/fastqc_clean


    summarize_reads(args.sample_names, args.fastqc_dir)

if __name__ == "__main__":
    main()
