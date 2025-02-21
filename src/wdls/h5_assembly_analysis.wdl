version 1.0

import "structs.wdl" as initializations
import "version_capture_tasks.wdl" as vc
import "reference_tasks.wdl" as rt
import "primer_tasks.wdl" as pt
import "other_tasks.wdl" as ot

workflow h5_assembly_analysis {
    input {
        Array[String] primers
        Array[String] samples
        Array[File] fastq1s
        Array[File] fastq2s
        String project_name
        String gs_dir
        File contaminants_fasta
        File version_capture_py
    }

    meta {
        allowNestedInputs: true
    }

    # private declarations
    String fastqc_docker = 'staphb/fastqc:0.12.1'
    String seqyclean_docker = 'staphb/seqyclean:1.10.09'
    String ivar_docker = 'staphb/ivar:1.4.2'
    String python_docker = 'ariannaesmith/py3.10.9-bio'
    String multiqc_docker = 'multiqc/multiqc:v1.24'
    String jammy_docker = 'ubuntu:jammy-20240627.1'
    String utility_docker = 'theiagen/utility:1.0'
    String h5_docker = 'ariannaesmith/cdphe_h5_influenza:latest'
    String version_capture_docker = 'ariannaesmith/cdphe_wdl_version_capture:latest'
    
    String workflow_name = 'h5_assembly_analysis'
    String workflow_version = 'v0.1.0-alpha'
    String workflow_version_und = 'v0_1_0-alpha'

    Array[Int] indexes = range(length(samples))

    call vc.workflow_metadata as w_meta { 
        input: 
            docker = jammy_docker,
            workflow_name = workflow_name,
            workflow_version = workflow_version
    }
    String project_outdir = gs_dir + "/" +  project_name + "/terra_outputs/" + workflow_version_und + "/"

    # Struct initilizations (subworkflow)
    call initializations.declare_structs as ini { input: h5_docker = h5_docker}

    # Scatter samples to create structs
    scatter (idx in indexes) {
        Sample sample = object {
            name: samples[idx],
            primer: primers[idx],
            fastq1: fastq1s[idx],
            fastq2: fastq2s[idx],
            i: idx
        }
    }

    Array[Sample] all_samples = sample

    # Group samples by primer
    scatter (ps in ini.primer_schemes) {
        scatter (all_samp in all_samples) {
            if (all_samp.primer == ps.name) {
                # Only add to list if fastqs are not empty
                Float fastqs_size = size([all_samp.fastq1, all_samp.fastq2], "MiB")
                if (fastqs_size > 1) {
                    Sample primer_sample = all_samp
                }
            }
        }
        Array[Sample] primer_samples = select_all(primer_sample)
        
        # Only call downstream tasks if primer was used
        if (length(primer_samples) > 0) {
            String p_name = ps.name
            String primer_outdir = project_outdir + p_name + "/"
            # Call primer level tasks (subworkflow)
            call pt.primer_level_tasks as p_sub {
                input:
                    primer_samples = primer_samples,
                    primer_outdir = primer_outdir,
                    project_name = project_name,
                    contaminants_fasta = contaminants_fasta,
                    fastqc_docker = fastqc_docker,
                    seqyclean_docker = seqyclean_docker,
                    multiqc_docker = multiqc_docker,
                    utility_docker = utility_docker,
                    h5_docker = h5_docker
            }

            Array[File] seqyclean_output = flatten([p_sub.cleaned_PE1, p_sub.cleaned_PE2])

            # Call reference level tasks (subworkflow)
            Array[Int] num_samples = range(length(primer_samples)) 
            String ref_name = ps.reference_name           
            
            call rt.reference_level_tasks as r_sub {
                input: 
                    reference_name = ps.reference_name,
                    reference_fasta = ps.reference_fasta,
                    project_name = project_name,
                    reference_outdir = primer_outdir,
                    num_samples = num_samples,
                    primer_samples = primer_samples,
                    cleaned_PE1 = p_sub.cleaned_PE1,
                    cleaned_PE2 = p_sub.cleaned_PE2,
                    fastqc_clean_summary_metrics = p_sub.fastqc_clean_summary_metrics,
                    primer_bed = ps.bed,
                    ivar_docker = ivar_docker,
                    python_docker = python_docker,
                    multiqc_docker = multiqc_docker,
                    utility_docker = utility_docker,
                    h5_docker = h5_docker
            }   
        }
    }

    # Collect various version information
    VersionInfo fastqc_version = select_first(p_sub.fastqc_version)
    VersionInfo seqyclean_version = select_first(p_sub.seqyclean_version)
    VersionInfo multiqc_version = select_first(p_sub.multiqc_version)
    VersionInfo h5_docker_version = select_first(p_sub.h5_docker_version)
    VersionInfo samtools_version = select_first(r_sub.samtools_version)
    VersionInfo bwa_version = select_first(r_sub.bwa_version)
    VersionInfo ivar_version = select_first(r_sub.ivar_version)
    Array[VersionInfo] version_array = [w_meta.version_info, fastqc_version, seqyclean_version, 
                                        multiqc_version, h5_docker_version, samtools_version, 
                                        bwa_version, ivar_version]

    call vc.capture_versions as version_cap {
        input:
            version_array = version_array,
            workflow_name = workflow_name,
            workflow_version = workflow_version_und,
            project_name = project_name,
            analysis_date = w_meta.analysis_date,
            docker = version_capture_docker
    }

    call ot.transfer as transfer_vc {
        input:
            out_dir = project_outdir,
            task_dir = 'summary_results',
            task_files = [version_cap.output_file],
            docker = utility_docker
    }

    output { 
        Array[String] primers_used = select_all(p_name)
        Array[Array[File]] primers_fastqc_raw_outputs = select_all(p_sub.fastqc_raw_outputs)
        Array[Array[File]] primers_fastqc_clean_outputs = select_all(p_sub.fastqc_clean_outputs)
        Array[Array[File]] primers_seqyclean_outputs = select_all(seqyclean_output)
        Array[Array[File]] primers_summary_outputs = select_all(p_sub.p_summary_outputs)
        Array[Array[File]] primers_alignment_outputs = select_all(r_sub.alignment_outputs)
        Array[Array[File]] primers_consensus_outputs = select_all(r_sub.consensus_outputs)
        Array[Array[File]] primers_ref_summary_outputs = select_all(r_sub.ref_summary_outputs)
        Array[VersionInfo] version_capture = version_array
    }    
}



