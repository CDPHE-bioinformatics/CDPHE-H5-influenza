version development

import "h5_structs.wdl" as sub

struct Sample {
    String name
    String primer
    File fastq1
    File fastq2 
    Int i
}

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

    # Struct initilizations
    call sub.declare_structs as s {
        input: 
            bovine_texas_029328_01_UtoT_ha_fasta = bovine_texas_029328_01_UtoT_ha_fasta,
            darwin_9_2021_h3n2_ha_h3_fasta = darwin_9_2021_h3n2_ha_h3_fasta,
            victoria_4897_2022_h1n1_ha_h1_fasta = victoria_4897_2022_h1n1_ha_h1_fasta,
            vietnam_1203_2024_h5n1_ha_v2_fasta = vietnam_1203_2024_h5n1_ha_v2_fasta,
            bovine_texas_029328_01_UtoT_fasta = bovine_texas_029328_01_UtoT_fasta,
            human_h5_250_bed = human_h5_250_bed,
            houston_bed = houston_bed,
            AVRL_H5N1_250bp_bed = AVRL_H5N1_250bp_bed
    }

    # Scatter samples to create structs
    scatter (idx in indexes) {
        Sample sample = Sample {
            name: samples[idx],
            primer: primers[idx],
            fastq1: fastq1s[idx],
            fastq2: fastq2s[idx],
            i: idx
        }
    }

    Array[Sample] all_samples = sample
    # Group samples by primer
    scatter (ps in s.primer_schemes) {
        
        scatter (all_samp in all_samples) {
            if (all_samp.primer == ps.name) {
                Sample primer_sample = all_samp
                Int match_index = all_samp.i
            }
        }

        Boolean primer_sample_matches = defined(primer_sample)
        if (!primer_sample_matches) {
            call exit_wdl {input: exit_reason = "primer not used"}
        }        

        # This first line is only included because the following fails without it?
        Array[Int] index_matches = select_all(match_index)
        Array[Sample] primer_samples = select_all(primer_sample)

        String primer_outdir = project_outdir + ps.name +"/"
        # Call primer level tasks
        scatter (samp in primer_samples) {
            # Call sample level tasks
            call fastqc as fastqc_raw {
                input: 
                    sample = samp, 
                    fastq1 = samp.fastq1, 
                    fastq2 = samp.fastq2, 
                    docker = fastqc_docker
            }

            # call sample_qc_file as sample_qc_file_raw {
            #     input: 
            #         sample_name = samp.name,
            #         fastqc1_data = fastqc_raw.fastqc1_data,
            #         fastqc2_data = fastqc_raw.fastqc2_data,
            #         docker = python_docker
            # }
            
            call seqyclean {
                input: 
                    sample = samp,
                    contaminants_fasta = contaminants_fasta,
                    docker = seqyclean_docker
            }

            call fastqc as fastqc_clean {
                input: 
                    sample = samp, 
                    fastq1 = seqyclean.PE1, 
                    fastq2 = seqyclean.PE2, 
                    docker = fastqc_docker
            }

            # call sample_qc_file as sample_qc_file_clean {
            #     input: 
            #         sample_name = samp.name,
            #         fastqc1_data = fastqc_clean.fastqc1_data,
            #         fastqc2_data = fastqc_clean.fastqc2_data,
            #         docker = python_docker
            # }
        }

        Array[String] task_dirs = ["fastqc_raw", "fastqc_clean", "seqyclean"]
        # Array[Array[File]] task_files = [flatten([fastqc_raw.fastqc1_data, fastqc_raw.fastqc2_data, sample_qc_file_raw.summary_metrics]),
        #                     flatten([fastqc_clean.fastqc1_data, fastqc_clean.fastqc2_data, sample_qc_file_clean.summary_metrics]),
        #                     flatten([seqyclean.PE1, seqyclean.PE2])]
        Array[Array[File]] task_files = [flatten([fastqc_raw.fastqc1_data, fastqc_raw.fastqc2_data]),
                            flatten([fastqc_clean.fastqc1_data, fastqc_clean.fastqc2_data]),
                            flatten([seqyclean.PE1, seqyclean.PE2])]
        call transfer as transfer_primer_tasks {
            input:
                out_dir = primer_outdir,
                task_dirs = task_dirs,
                task_files = task_files,
                docker = utility_docker
        }
        
    }

}

task transfer {
    input {
        String out_dir
        Array[String] task_dirs
        Array[Array[File]] task_files
        String docker
    }

    # Have to re-declare variables in bash due to syntax clash
    command <<<
        task_dirs_bash=(~{sep(' ', task_dirs)})
        task_files_bash=~{task_files}
        for i in "${!task_dirs_bash[@]}"; 
        do
            files="${task_files_bash[$i]//,/}"
            gsutil -m cp "${files}" "~{out_dir}${task_dirs_bash[$i]}/"; 
        done;
    >>>
    runtime {
        #cpu: ,
        #memory: ,
        docker: docker
    }
}

task multiqc {
    input {
        Array[File] fastq1s
        Array[File] fastq2s
        String fastq_type
        String docker
    }
    String fastq_dir = "~{fastq_type}_fastqs/"
    String out_dir = "multiqc_~{fastq_type}/"

    command <<<
        mkdir ~{fastq_dir}
        cp ~{sep(' ', fastq1s)} ~{fastq_dir}
        cp ~{sep(' ', fastq2s)} ~{fastq_dir}
        multiqc ~{fastq_dir} --outdir
    >>>

    output {
        Array[File] fastqs_dir = glob("~{out_dir}*")
    }
    runtime {
        #cpu: ,
        #memory: ,
        docker: docker
    }
}

task fastqc {
    input {
        Sample sample
        File fastq1
        File fastq2
        String docker
    }

    String fastq1_name = basename(fastq1, ".fastq.gz")
    String fastq2_name = basename(fastq2, ".fastq.gz")

    command <<<
        fastqc --outdir $PWD --extract --delete ~{sample.fastq1} ~{sample.fastq2}
        fastqc --version | awk '/FastQC/ {print $2}' | tee VERSION  
        cp "~{fastq1_name}_fastqc/fastqc_data.txt" "~{fastq1_name}_fastqc_data.txt"
        cp "~{fastq2_name}_fastqc/fastqc_data.txt" "~{fastq2_name}_fastqc_data.txt"  
    >>>

    output {
        File fastqc1_data = "~{fastq1_name}_fastqc_data.txt"
        File fastqc2_data = "~{fastq2_name}_fastqc_data.txt"
        String version = read_string('VERSION')
    }

    runtime {
        #cpu: ,
        #memory: ,
        docker: docker
    }
}

task sample_qc_file {
    input {
        String sample_name
        File fastqc1_data
        File fastqc2_data
        String docker
    }

    command <<<
        # python things
    >>>

    output {
        File summary_metrics = "${sample_name}_summary_metrics.tsv"
    }

    runtime {
        #cpu: ,
        #memory: ,
        docker: docker
    }
}

task seqyclean {
    input {
        Sample sample
        File contaminants_fasta
        String docker
    }

    String out_name = "seqyclean/~{sample.name}_clean"

    command <<<
        seqyclean -minlen 25 -qual 30 30 -gz -1 ~{sample.fastq1} -2 ~{sample.fastq2} -c ~{contaminants_fasta} -o ~{out_name}
    >>>

    output {
        File PE1 = "seqyclean/~{sample.name}_clean_PE1.fastq.gz"
        File PE2 = "seqyclean/~{sample.name}_clean_PE2.fastq.gz"
    }

    runtime {
        #cpu: ,
        #memory: ,
        docker: docker
    }
}

task align_bwa {
    input {
        String sample_name
        File fastq1
        File fastq2
        String reference_name
        File reference_fasta
        String docker
    }

    String sam_fn = sample_name + ".sam"
    String bam_fn = sample_name + "aln.sorted.bam"

    command <<<
        bwa index -p ~{reference_name} -a is ~{reference_fasta}
        bwa mem -t 6 ~{reference_name} ~{fastq1} ~{fastq2} -f ~{sam_fn}
        samtools view -bS ~{sam_fn} | samtools sort -o {bam_fn}
    >>>

    output {
        File bam = bam_fn
    }

    runtime {
        #cpu: ,
        #memory: ,
        docker: docker
    }
}

task trim_primers_ivar {
    input {
        String sample_name
        File bam
        File primer_bed
        String docker
    }

    String trim_fn = sample_name + "_trimmed.bam"
    String trim_sort_bam_fn = sample_name + "_trimmed.sorted.bam"
    String trim_sort_bai_fn = sample_name + "_trimmed.sorted.bai"
    
    command <<<
        ivar trim -e -i ~{bam} -b ~{primer_bed} -p ~{trim_fn}
        samtools sort ~{trim_fn} -o ~{trim_sort_bam_fn}
        samtools index ~{trim_sort_bam_fn}
    >>>

    output {
        File trimmed_sorted_bam = trim_sort_bam_fn
        File trimmed_sorted_bai = trim_sort_bai_fn
    }

    runtime {
        #cpu: ,
        #memory: ,
        docker: docker
    }
}

task generate_consensus_ivar {
    input {
        String sample_name
        File trimmed_sorted_bam
        File reference_fasta
       String docker
    }

    String pileup_fn = sample_name + "_pileup.txt"
    String consensus_fn_prefix = sample_name + "_consensus"

    command <<<
        samtools faidx ~{reference_fasta}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ~{reference_fasta} ~{trimmed_sorted_bam} -o ~{pileup_fn}
        cat ~{pileup_fn} | ivar consensus -p ~{consensus_fn_prefix} -q 20 -t 0.6 -m 10
    >>>

    output {
        File pileup = pileup_fn
        File consensus_fasta = consensus_fn_prefix + ".fa"
    }

    runtime {
        #cpu: ,
        #memory: ,
        docker: docker
    }
}

task alignment_metrics {
    input {
        String sample_name
        File trimmed_sorted_bam
        File reference_fasta
        String docker
    }

    String coverage_fn = sample_name + "_coverage.txt"
    String stats_fn = sample_name + "_stats.txt"

    command <<<
        samtools coverage -o ~{coverage_fn} ~{trimmed_sorted_bam}
        samtools stats ~{trimmed_sorted_bam} > ~{stats_fn}
    >>>

    output {
        File coverage = coverage_fn
        File stats = stats_fn
    }

    runtime {
        #cpu: ,
        #memory: ,
        docker: docker
    }
}

task calculate_coverage_stats {
    input {
        Int ref_length
        Int seq_ref_length
        String sample_name
        File consensus_fasta
        String docker
    }

    command <<<
        # python stuff
    >>>

    output {
        File coverage_stats = "~{sample_name}_coverage_stats.csv"
    }

    runtime {
        #cpu: ,
        #memory: ,
        docker: docker
    }
}

task concat_sample_reference_metrics {
    input {
        String sample_name
        File samtools_coverage
        File samtools_stats
        File coverage_stats
        String docker
    }

    command <<<
        # python things
    >>>

    output {
        File sample_reference_metrics = "~{sample_name}_reference_metrics.csv"
    }

    runtime {
        #cpu: ,
        #memory: ,
        docker: docker
    }
}

task concat_sample_metrics {
    input {
        String sample_name
        Array[File] sample_reference_metrics
        String docker
    }

    command <<<
        # python things
    >>>

    output {
        File sample_metrics = "~{sample_name}_metrics.csv"
    }

    runtime {
        #cpu: ,
        #memory: ,
        docker: docker
    }
}

task concat_all_samples_metrics {
    input {
        Array[File] sample_metrics_files
        String docker
    }

    command <<<
        # python things
    >>>

    output {
        File samples_metrics = "summary_metrics.csv"
    }

    runtime {
        #cpu: ,
        #memory: ,
        docker: docker
    }
}

task exit_wdl {
    input {
        String exit_reason
    }
    
    command <<<
        echo "~{exit_reason}"
        exit 1
    >>>

    runtime {
        returnCodes: 0
        docker: "ubuntu:latest"
    }
}