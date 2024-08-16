version development

import "h5_structs.wdl" as sub

workflow h5 {
    input {
        Array[String] primers
        Array[String] samples
        Array[File] fastq1s
        Array[File] fastq2s
        String project_name
        String gs_dir
        File contaminants_fasta
    }

    meta {
        allowNestedInputs: true
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

    # Struct initilizations (subworkflow)
    call sub.declare_structs as s {}

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

        # Only call downstream tasks if primer was used
        Array[Sample] primer_samples = select_all(primer_sample)
        if (length(primer_samples) > 0) {          

            # Call primer level tasks
            scatter (p_samp in primer_samples) {
                # Call sample level tasks
                call fastqc as fastqc_raw {
                    input: 
                        fastq1 = p_samp.fastq1, 
                        fastq2 = p_samp.fastq2, 
                        docker = fastqc_docker
                }

                # Anything calling/using output from a python script task is commented out right now
                # call sample_qc_file as sample_qc_file_raw {input: 
                #         sample_name = p_samp.name,
                #         fastqc1_data = fastqc_raw.fastqc1_data,
                #         fastqc2_data = fastqc_raw.fastqc2_data,
                #         docker = python_docker
                # }
                
                call seqyclean {
                    input: 
                        sample = p_samp,
                        contaminants_fasta = contaminants_fasta,
                        docker = seqyclean_docker
                }

                call fastqc as fastqc_clean {
                    input: 
                        fastq1 = seqyclean.PE1, 
                        fastq2 = seqyclean.PE2, 
                        docker = fastqc_docker
                }

                # call sample_qc_file as sample_qc_file_clean {
                #     input: 
                #         sample_name = p_samp.name,
                #         fastqc1_data = fastqc_clean.fastqc1_data,
                #         fastqc2_data = fastqc_clean.fastqc2_data,
                #         docker = python_docker
                # }
    
            }

            # Transfer primer level files
            String primer_outdir = project_outdir + ps.name + "/"
            Array[String] primer_task_dirs = ["fastqc_raw", "fastqc_clean", "seqyclean"]
            Array[Array[File]] primer_task_files = [flatten([fastqc_raw.fastqc1_data, fastqc_raw.fastqc2_data]),
                                flatten([fastqc_clean.fastqc1_data, fastqc_clean.fastqc2_data]),
                                flatten([seqyclean.PE1, seqyclean.PE2])]       

            scatter (dir_files in zip(primer_task_dirs, primer_task_files)) {       
                call transfer as transfer_primer_tasks {
                    input:
                        out_dir = primer_outdir,
                        task_dir = dir_files.left,
                        task_files = dir_files.right,
                        docker = utility_docker
                }
            }

            # Call reference level tasks
            Array[Int] num_samples = range(length(primer_samples))            
            scatter (p_ref in ps.references) {
                String reference_outdir = primer_outdir + p_ref.name + "/"

                scatter (n in num_samples) {
                    Sample r_samp = primer_samples[n]
                    File PE1 = seqyclean.PE1[n]
                    File PE2 = seqyclean.PE2[n]

                    call align_bwa {
                        input:
                            sample_name = r_samp.name,
                            fastq1 = PE1,
                            fastq2 = PE2,
                            reference_name = p_ref.name,
                            reference_fasta = p_ref.fasta,
                            docker = viral_core_docker
                    }  

                    call trim_primers_ivar {
                        input: 
                            sample_name = r_samp.name,
                            bam = align_bwa.bam,
                            primer_bed = ps.bed,
                            docker = ivar_docker
                    }

                    call generate_consensus_ivar {
                        input: 
                            sample_name = r_samp.name,
                            trim_sort_bam = trim_primers_ivar.trim_sort_bam,
                            reference_fasta = p_ref.fasta,
                            docker = samtools_docker
                    }   

                    call alignment_metrics {
                        input:
                            sample_name = r_samp.name,
                            trim_sort_bam = trim_primers_ivar.trim_sort_bam,
                            docker = samtools_docker
                    }

                    # call calculate_coverage_stats {
                    #     input:
                    #         sample_name = r_samp.name,
                    #         ref_length = p_ref.length,
                    #         seq_ref_length = p_ref.captured_length,
                    #         consensus_fasta = generate_consensus_ivar.consensus_fasta,
                    #         docker = python_docker
                    # }

                    # call concat_sample_reference_metrics {
                    #     input: 
                    #         sample_name = r_samp.name,
                    #         samtools_coverage = alignment_metrics.coverage,
                    #         samtools_stats = alignment_metrics.stats,
                    #         coverage_stats = calculate_coverage_stats.coverage_stats,
                    #         docker = python_docker
                    # }

                } 

                Array[String] reference_task_dirs = ["alignments", "consensus_sequences", "metrics"]
                Array[Array[File]] reference_task_files = [flatten([trim_primers_ivar.trim_sort_bam, trim_primers_ivar.trim_sort_bai,]),
                                                        generate_consensus_ivar.consensus_fasta,
                                                        flatten([alignment_metrics.coverage, alignment_metrics.stats])]

                scatter (dir_files in zip(reference_task_dirs, reference_task_files)) {       
                    call transfer as transfer_reference_tasks {
                        input:
                            out_dir = reference_outdir,
                            task_dir = dir_files.left,
                            task_files = dir_files.right,
                            docker = utility_docker
                    }
                }
            }
        }


    }
}

task transfer {
    input {
        String out_dir
        String task_dir
        Array[File] task_files
        String docker
    }

    command <<<
        cat "~{write_lines(task_files)}" | gsutil -m cp -I "~{out_dir}~{task_dir}/"
    >>>

    runtime {
        #cpu: 
        #memory: 
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
        #cpu: 
        #memory: 
        docker: docker
    }
}

task fastqc {
    input {
        File fastq1
        File fastq2
        String docker
    }

    String fastq1_name = basename(fastq1, ".fastq.gz")
    String fastq2_name = basename(fastq2, ".fastq.gz")

    command <<<
        fastqc --outdir $PWD --extract --delete ~{fastq1} ~{fastq2}
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
        #cpu: 
        #memory: 
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
        #cpu: 
        #memory: 
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
        #cpu: 
        #memory: 
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

    String sam_fn = "~{sample_name}.sam"
    String bam_fn = "~{sample_name}_aln.sorted.bam"

    command <<<
        samtools --version-only | tee samtools_version
        bwa index -p ~{reference_name} -a is ~{reference_fasta}
        bwa mem -t 6 ~{reference_name} ~{fastq1} ~{fastq2} > ~{sam_fn}
        samtools view -b -@ 6 ~{sam_fn} | samtools sort -m 2G -@ 6 -o ~{bam_fn}
    >>>

    output {
        String samtools_version = read_string('samtools_version')
        File bam = bam_fn
    }

    runtime {
        cpu: 3
        memory: "12 GB"
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

    String trim_fn = "~{sample_name}_trimmed.bam"
    String trim_sort_bam_fn = "~{sample_name}_trimmed.sorted.bam"
    String trim_sort_bai_fn = "~{sample_name}_trimmed.sorted.bai"
    
    command <<<
        ivar trim -e -i ~{bam} -b ~{primer_bed} -p ~{trim_fn}
        samtools sort -@ 6 -o ~{trim_sort_bam_fn} ~{trim_fn}
        samtools index -@ 6 ~{trim_sort_bam_fn} -o ~{trim_sort_bai_fn}
    >>>

    output {
        File trim_sort_bam = trim_sort_bam_fn
        File trim_sort_bai = trim_sort_bai_fn
    }

    runtime {
        #cpu: 
        #memory: 
        docker: docker
    }
}

task generate_consensus_ivar {
    input {
        String sample_name
        File trim_sort_bam
        File reference_fasta
       String docker
    }

    String pileup_fn = "~{sample_name}_pileup.txt"
    String consensus_fn_prefix = "~{sample_name}_consensus"

    command <<<
        samtools faidx ~{reference_fasta}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ~{reference_fasta} ~{trim_sort_bam} -o ~{pileup_fn}
        cat ~{pileup_fn} | ivar consensus -p ~{consensus_fn_prefix} -q 20 -t 0.6 -m 10
    >>>

    output {
        File pileup = pileup_fn
        File consensus_fasta = consensus_fn_prefix + ".fa"
    }

    runtime {
        #cpu: 
        #memory: 
        docker: docker
    }
}

task alignment_metrics {
    input {
        String sample_name
        File trim_sort_bam
        String docker
    }

    String coverage_fn = "~{sample_name}_coverage.txt"
    String stats_fn = "~{sample_name}_stats.txt"

    command <<<
        samtools coverage -o ~{coverage_fn} ~{trim_sort_bam}
        samtools stats ~{trim_sort_bam} > ~{stats_fn}
    >>>

    output {
        File coverage = coverage_fn
        File stats = stats_fn
    }

    runtime {
        #cpu: 
        #memory: 
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
        #cpu: 
        #memory: 
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
        # python stuff
    >>>

    output {
        File sample_reference_metrics = "~{sample_name}_reference_metrics.csv"
    }

    runtime {
        #cpu: 
        #memory: 
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
        # python stuff
    >>>

    output {
        File sample_metrics = "~{sample_name}_metrics.csv"
    }

    runtime {
        #cpu: 
        #memory: 
        docker: docker
    }
}

task concat_all_samples_metrics {
    input {
        Array[File] sample_metrics_files
        String docker
    }

    command <<<
        # python stuff
    >>>

    output {
        File samples_metrics = "summary_metrics.csv"
    }

    runtime {
        #cpu: 
        #memory: 
        docker: docker
    }
}

