version development

import "h5_structs.wdl" as sub
import "version_capture_task.wdl" as vc

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
    #String samtools_docker = 'staphb/samtools:1.10'
    String ivar_docker = 'staphb/ivar:1.4.2'
    String python_docker = 'ariannaesmith/py3.10.9-bio'
    String viral_core_docker = 'quay.io/broadinstitute/viral-core:2.2.3'
    String multiqc_docker = 'multiqc/multiqc:v1.24'
    String jammy_docker = 'ubuntu:jammy-20240627.1'
    String utility_docker = 'theiagen/utility:1.0'

    Array[Int] indexes = range(length(samples))

    call vc.workflow_metadata as w_meta {}
    String project_outdir = gs_dir + "/" +  project_name + "/terra_outputs/" + w_meta.version + "/"

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
            String p_name = ps.name
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

            # Call multiqc
            String fastqc_cl_config = "sp: { fastqc/data: { fn: '*_fastqc_data.txt' } }" 
            call multiqc as multiqc_raw {
                input:
                    files = flatten([fastqc_raw.fastqc1_data, fastqc_raw.fastqc2_data]),
                    module = "fastqc",
                    task_name = "raw",
                    cl_config = fastqc_cl_config,
                    docker = multiqc_docker
            }

            call multiqc as multiqc_clean {
                input:
                    files = flatten([fastqc_clean.fastqc1_data, fastqc_clean.fastqc2_data]),
                    module = "fastqc",
                    task_name = "clean",
                    cl_config = fastqc_cl_config,
                    docker = multiqc_docker
            }

            call multiqc as multiqc_seqyclean {
                input:
                    files = seqyclean.summary_stats,
                    module = "seqyclean",
                    docker = multiqc_docker
            }

            
            # Transfer primer level files
            String primer_outdir = project_outdir + p_name + "/"
            Array[File] fastqc_raw_output = flatten([fastqc_raw.fastqc1_data, fastqc_raw.fastqc2_data])
            Array[File] fastqc_clean_output = flatten([fastqc_clean.fastqc1_data, fastqc_clean.fastqc2_data])
            Array[File] seqyclean_output = flatten([seqyclean.PE1, seqyclean.PE2])
            Array[File] p_summary_output = [multiqc_raw.html_report, multiqc_clean.html_report, multiqc_seqyclean.html_report]

            Array[String] primer_task_dirs = ["fastqc_raw", "fastqc_clean", "seqyclean", "summary_files"]
            Array[Array[File]] primer_task_files = [fastqc_raw_output, fastqc_clean_output, seqyclean_output, p_summary_output]       

            scatter (dir_files in zip(primer_task_dirs, primer_task_files)) {       
                call transfer as transfer_primer_tasks {
                    input:
                        out_dir = primer_outdir,
                        task_dir = dir_files.left,
                        task_files = dir_files.right,
                        docker = utility_docker
                }
            }

            # Grab primer level task versions
            # VersionInfo fastqc_version = select_first(fastqc_raw.version_info)
            # VersionInfo seqyclean_version = select_first(seqyclean.version_info)
            # VersionInfo multiqc_version = multiqc_raw.version_info

            # Call reference level tasks
            Array[Int] num_samples = range(length(primer_samples))            
            scatter (p_ref in ps.references) {
                String ref_name = p_ref.name
                String reference_outdir = primer_outdir + ref_name + "/"

                scatter (n in num_samples) {
                    Sample r_samp = primer_samples[n]
                    File PE1 = seqyclean.PE1[n]
                    File PE2 = seqyclean.PE2[n]

                    call align_bwa {
                        input:
                            sample_name = r_samp.name,
                            fastq1 = PE1,
                            fastq2 = PE2,
                            reference_name = ref_name,
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
                            docker = ivar_docker
                    }   

                    call alignment_metrics {
                        input:
                            sample_name = r_samp.name,
                            trim_sort_bam = trim_primers_ivar.trim_sort_bam,
                            docker = ivar_docker
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

                String samtools_cl_config = "extra_fn_clean_exts: [ '_coverage', '_stats']"
                call multiqc as multiqc_samtools {
                    input:
                        files = flatten([trim_primers_ivar.idxstats, alignment_metrics.coverage, alignment_metrics.stats]),
                        module = "samtools",
                        task_name = "alignment",
                        cl_config = samtools_cl_config,
                        docker = multiqc_docker
                }
                
                # Transfer reference level files
                Array[File] alignment_output = flatten([trim_primers_ivar.trim_sort_bam, trim_primers_ivar.trim_sort_bai])
                Array[File] consensus_output = generate_consensus_ivar.consensus_fasta
                Array[File] ref_summary_output = flatten([alignment_metrics.coverage, alignment_metrics.stats, [multiqc_samtools.html_report]])
                
                Array[String] reference_task_dirs = ["alignments", "consensus_sequences", "metrics"]
                Array[Array[File]] reference_task_files = [alignment_output, consensus_output, ref_summary_output]

                scatter (dir_files in zip(reference_task_dirs, reference_task_files)) {       
                    call transfer as transfer_reference_tasks {
                        input:
                            out_dir = reference_outdir,
                            task_dir = dir_files.left,
                            task_files = dir_files.right,
                            docker = utility_docker
                    }
                }

                # Grab reference level task versions
                # VersionInfo viral_core_samtools_version = select_first(align_bwa.samtools_version_info)
                # VersionInfo ivar_samtools_version = select_first(trim_primers_ivar.samtools_version_info)
                # VersionInfo ivar_version = select_first(trim_primers_ivar.ivar_version_info)
            }

            Array[String] refs_used = ref_name
            Array[Array[File]] refs_alignment_outputs = alignment_output
            Array[Array[File]] refs_consensus_outputs = consensus_output
            Array[Array[File]] refs_summary_outputs = ref_summary_output
        }
    }

    VersionInfo fastqc_version = select_first(select_first(fastqc_raw.version_info))
    VersionInfo seqyclean_version = select_first(select_first(seqyclean.version_info))
    VersionInfo multiqc_version = select_first(multiqc_raw.version_info)
    VersionInfo viral_core_samtools_version = select_first(select_first(select_first(align_bwa.samtools_version_info)))
    VersionInfo ivar_samtools_version = select_first(select_first(select_first(trim_primers_ivar.samtools_version_info)))
    VersionInfo ivar_version = select_first(select_first(select_first(trim_primers_ivar.ivar_version_info)))
    Array[VersionInfo] version_array = [fastqc_version, seqyclean_version, multiqc_version, 
                                viral_core_samtools_version, ivar_samtools_version, ivar_version]


    # VersionInfo fastqc_v = select_first(fastqc_version)
    # VersionInfo seqyclean_v = select_first(seqyclean_version)
    # VersionInfo multiqc_v = select_first(multiqc_version)
    # VersionInfo viral_core_samtools_v = select_first(select_first(viral_core_samtools_version))
    # VersionInfo ivar_samtools_v = select_first(select_first(ivar_samtools_version))
    # VersionInfo ivar_v = select_first(select_first(ivar_version))
    # Array[VersionInfo] version_array = [fastqc_v, seqyclean_v, multiqc_v, viral_core_samtools_v, ivar_samtools_v, ivar_v]

    # call vc.capture_versions {
    #     input:
    #         version_array,
    #         workflow_name = 'h5_testing',
    #         workflow_version = w_meta.version,
    #         project_name,
    #         project_outdir,
    #         analysis_date = w_meta.analysis_date,
    #         version_capture_py,
    #         docker = 'lets make a new docker'
    # }

    output { 
        Array[String] primers_used = select_all(p_name)
        Array[Array[File]] p_fastqc_raw_outputs = select_all(fastqc_raw_output)
        Array[Array[File]] p_fastqc_clean_outputs = select_all(fastqc_clean_output)
        Array[Array[File]] p_seqyclean_outputs = select_all(seqyclean_output)
        Array[Array[File]] p_summary_outputs = select_all(p_summary_output)
        Array[Array[String]] p_refs_used = select_all(refs_used)
        Array[Array[Array[File]]] p_refs_alignment_outputs = select_all(refs_alignment_outputs)
        Array[Array[Array[File]]] p_refs_consensus_outputs = select_all(refs_consensus_outputs)
        Array[Array[Array[File]]] p_refs_summary_outputs = select_all(refs_summary_outputs)
        Array[VersionInfo] version_capture = version_array
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
        Array[File] files
        String module
        String? task_name
        String? cl_config
        String docker
    }

    String? task_prefix = if defined(task_name) then "~{module}_~{task_name}" else None
    String html_fn = "~{select_first([task_prefix, module])}_multiqc_report.html"

    command <<<
        multiqc -m "~{module}" -l ~{write_lines(files)} -n ~{html_fn} ~{if defined(cl_config) then '--cl-config "${cl_config}"' else ''} 
        multiqc --version | awk ' {print $3}' | tee VERSION
    >>>

    output {
        File html_report = html_fn
        VersionInfo version_info = VersionInfo {
            software: "multiqc",
            docker: docker,
            version: read_string('VERSION')
        }
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
        VersionInfo version_info = VersionInfo {
            software: "fastqc",
            docker: docker,
            version: read_string('VERSION')
        }
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

    String out_name = "~{sample.name}_clean"
    # seqyclean doesn't have a version command
    String VERSION = sub(docker, "staphb/seqyclean:", "")

    command <<<
        seqyclean -minlen 25 -qual 30 30 -gz -1 ~{sample.fastq1} -2 ~{sample.fastq2} -c ~{contaminants_fasta} -o ~{out_name}
    >>>

    output {
        File PE1 = "~{sample.name}_clean_PE1.fastq.gz"
        File PE2 = "~{sample.name}_clean_PE2.fastq.gz"
        File summary_stats = "~{sample.name}_clean_SummaryStatistics.tsv"
        VersionInfo version_info = VersionInfo {
            software: "seqyclean",
            docker: docker,
            version: VERSION
        }
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
        samtools --version-only | tee SAMTOOLS_VERSION
        bwa index -p ~{reference_name} -a is ~{reference_fasta}
        bwa mem -t 6 ~{reference_name} ~{fastq1} ~{fastq2} > ~{sam_fn}
        samtools view -b -@ 6 ~{sam_fn} | samtools sort -m 2G -@ 6 -o ~{bam_fn}
    >>>

    output {
        File bam = bam_fn
        VersionInfo samtools_version_info = VersionInfo {
            software: "samtools",
            docker: docker,
            version: read_string('SAMTOOLS_VERSION')
        }
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
    String idxstats_fn = "~{sample_name}_trimmed_sorted_idxstats.txt"
    
    command <<<
        ivar trim -e -i ~{bam} -b ~{primer_bed} -p ~{trim_fn}
        samtools sort -@ 6 -o ~{trim_sort_bam_fn} ~{trim_fn}
        samtools index -@ 6 ~{trim_sort_bam_fn} -o ~{trim_sort_bai_fn}
        samtools idxstats ~{trim_sort_bam_fn} > ~{idxstats_fn}
        ivar version | awk ' /iVar/ {print $3}' | tee IVAR_VERSION
        samtools --version-only | tee SAMTOOLS_VERSION
    >>>

    output {
        File trim_sort_bam = trim_sort_bam_fn
        File trim_sort_bai = trim_sort_bai_fn
        File idxstats = idxstats_fn        
        VersionInfo samtools_version_info = VersionInfo {
            software: "samtools",
            docker: docker,
            version: read_string('SAMTOOLS_VERSION')
        }
        VersionInfo ivar_version_info = VersionInfo {
            software: "ivar",
            docker: docker,
            version: read_string('IVAR_VERSION')
        }
    }

    runtime {
        cpu: 3
        memory: "12 GB"
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
        cpu: 2
        memory: "4 GB"
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

