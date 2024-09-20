#!/usr/bin/env python 

# main.py

import subprocess

def run_script(script_name, *args):
    command = f'python {script_name} {" ".join(args)}'
    subprocess.run(command, shell=True)

if __name__ == "__main__":
    # paths and parameters (to be edited)
    sample_names = 'list_of_samples'
    ref_seq_path_dict = {'ref1': '/path/to/ref1.fasta'}
    aln_metrics_dir = '/path/to/alignment_metrics'
    consen_dir = '/path/to/consensus_sequences'
    wd = '/path/to/working_directory'

    # running each script:
    run_script('alignment_metrics.py', *sample_names, ref_seq_path_dict, aln_metrics_dir, consen_dir, wd)
    
