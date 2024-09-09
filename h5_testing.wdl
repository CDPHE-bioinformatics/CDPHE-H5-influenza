version 1.0

workflow h5 {
    input {
        Array[String] primers
        Array[String] samples
        Array[File] fastq1s
        Array[File] fastq2s
        String project_name
        String gs_dir
        File contaminants_fasta
        File human_h5_250_bed
        File houston_bed
        File AVRL_H5N1_250bp_bed

        # All of these reference fastas are HA gene segments
        File bovine_texas_029328_01_UtoT_ha_fasta
        File darwin_9_2021_h3n2_ha_h3_fasta
        File victoria_4897_2022_h1n1_ha_h1_fasta
        File vietnam_1203_2024_h5n1_ha_v2_fasta    
    
        # Whole genome references
        File bovine_texas_029328_01_UtoT_fasta

        # temporary because seqyclean is slow
        Array[File] clean_fastq1s
        Array[File] clean_fastq2s
    }

    # private declarations
    String fastqc_docker = 'staphb/fastqc:0.12.1'
    String seqyclean_docker = 'staphb/seqyclean:1.10.09'
    String bwa_docker = 'staphb/bwa:0.7.17'
    String samtools_docker = 'staphb/samtools:1.10'
    String ivar_docker = 'staphb/ivar:1.4.2'
    String python_docker = 'ariannaesmith/py3.10.9-bio'
    String viral_core_docker = 'quay.io/broadinstitute/viral-core:2.2.3'
    String multiqc_docker = 'staphb/multiqc:1.8'
    String jammy_docker = 'ubuntu:jammy-20240627.1'
    String utility_docker = 'theiagen/utility:1.0'

    Array[Int] indexes = range(length(samples))

    String project_outdir = gs_dir + "/" +  project_name + "/"
}