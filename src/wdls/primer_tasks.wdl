version 1.0

import "other_tasks.wdl" as ot

workflow primer_level_tasks {
    input {
        Array[Sample] primer_samples
        String primer_outdir
        String project_name
        File contaminants_fasta
        String fastqc_docker
        String seqyclean_docker
        String multiqc_docker
        String utility_docker
        String h5_docker
    }

    scatter (sample in primer_samples) {
        String sample_name = sample.name
        # Call sample level tasks
        call fastqc as fastqc_raw {
            input: 
                fastq1 = sample.fastq1, 
                fastq2 = sample.fastq2, 
                docker = fastqc_docker
        }
        
        call seqyclean {
            input: 
                sample = sample,
                contaminants_fasta = contaminants_fasta,
                docker = seqyclean_docker
        }

        call fastqc as fastqc_clean {
            input: 
                fastq1 = seqyclean.PE1, 
                fastq2 = seqyclean.PE2, 
                docker = fastqc_docker
        }

    }

    call summarize_fastqc as summarize_fastqc_raw {
        input: 
            sample_names = sample_name,
            fastqc1_data_array = fastqc_raw.fastqc1_data,
            fastqc2_data_array = fastqc_raw.fastqc2_data,
            fastqc_type = "raw",
            docker = h5_docker
    }

    call summarize_fastqc as summarize_fastqc_clean {
        input: 
            sample_names = sample_name,
            fastqc1_data_array = fastqc_clean.fastqc1_data,
            fastqc2_data_array = fastqc_clean.fastqc2_data,
            fastqc_type = "clean",
            docker = h5_docker
    }

    call concat_fastqc_summary {
        input:
            sample_names = sample_name,
            summarized_fastqcs = flatten([summarize_fastqc_raw.summary_metrics, 
                                            summarize_fastqc_clean.summary_metrics]),
            project_name = project_name,
            docker = h5_docker
    }

    # Call multiqc
    call ot.multiqc as multiqc_fastqc {
        input:
            files = flatten([fastqc_raw.fastqc1_data, fastqc_raw.fastqc2_data,
                            fastqc_clean.fastqc1_data, fastqc_clean.fastqc2_data]),
            module = "fastqc",
            cl_config = "sp: { fastqc/data: { fn: '*_fastqc_data.txt' } }",
            docker = multiqc_docker
    }

    call ot.multiqc as multiqc_seqyclean {
        input:
            files = seqyclean.summary_stats,
            module = "seqyclean",
            cl_config = 'extra_fn_clean_exts: ["_clean"]',
            docker = multiqc_docker
    }
    
    # Transfer primer level files
    
    Array[File] fastqc_raw_output = flatten([fastqc_raw.fastqc1_data, fastqc_raw.fastqc2_data])
    Array[File] fastqc_clean_output = flatten([fastqc_clean.fastqc1_data, fastqc_clean.fastqc2_data])
    Array[File] seqyclean_output = flatten([seqyclean.PE1, seqyclean.PE2])
    Array[File] p_summary_output = [multiqc_fastqc.html_report, multiqc_seqyclean.html_report, concat_fastqc_summary.fastqc_summary]

    Array[String] primer_task_dirs = ["fastqc_raw", "fastqc_clean", "seqyclean", "summary_results"]
    Array[Array[File]] primer_task_files = [fastqc_raw_output, fastqc_clean_output, seqyclean_output, p_summary_output]       

    scatter (dir_files in zip(primer_task_dirs, primer_task_files)) {       
        call ot.transfer {
            input:
                out_dir = primer_outdir,
                task_dir = dir_files.left,
                task_files = dir_files.right,
                docker = utility_docker
        }
    }

    output {
        Array[File] fastqc_raw_outputs = fastqc_raw_output
        Array[File] fastqc_clean_outputs = fastqc_clean_output
        Array[File] cleaned_PE1 = seqyclean.PE1
        Array[File] cleaned_PE2 = seqyclean.PE2
        Array[File] p_summary_outputs = p_summary_output
        Array[File] fastqc_clean_summary_metrics = summarize_fastqc_clean.summary_metrics
        File fastqc_summary = concat_fastqc_summary.fastqc_summary
        VersionInfo fastqc_version = select_first(fastqc_raw.version_info)
        VersionInfo seqyclean_version = select_first(seqyclean.version_info)
        VersionInfo multiqc_version = multiqc_fastqc.version_info
        VersionInfo h5_docker_version = concat_fastqc_summary.version_info
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
        fastqc --version | awk "/FastQC/ {print $2}" | tee VERSION  
        cp "~{fastq1_name}_fastqc/fastqc_data.txt" "~{fastq1_name}_fastqc_data.txt"
        cp "~{fastq2_name}_fastqc/fastqc_data.txt" "~{fastq2_name}_fastqc_data.txt"  
    >>>

    output {
        File fastqc1_data = "~{fastq1_name}_fastqc_data.txt"
        File fastqc2_data = "~{fastq2_name}_fastqc_data.txt"
        String version = read_string('VERSION')
        VersionInfo version_info = {
            "software": "fastqc",
            "docker": docker,
            "version": read_string('VERSION')
        }
    }

    runtime {
        #cpu: 
        #memory: 
        docker: docker
    }
}

task summarize_fastqc {
    input {
        Array[String] sample_names
        Array[File] fastqc1_data_array
        Array[File] fastqc2_data_array
        String fastqc_type
        String docker
    }

    meta {
        volatile: true
    }
    
    command <<<
        summarize_fastqc.py --sample_names ~{sep=" " sample_names} \
            --fastqc1_data_array ~{sep=" " fastqc1_data_array} \
            --fastqc2_data_array ~{sep=" "  fastqc2_data_array} \
            --fastqc_type ~{fastqc_type}
    >>>

    output {
        Array[File] summary_metrics = glob("*_summary_metrics.tsv")
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
        VersionInfo version_info = {
            "software": "seqyclean",
            "docker": docker,
            "version": VERSION
        }
    }

    runtime {
        cpu: 8
        memory: "8G"
        docker: docker
    }
}

task concat_fastqc_summary {
    input {
        Array[String] sample_names
        Array[File] summarized_fastqcs
        String project_name
        String docker
    }

    meta {
        volatile: true
    }

    command <<<
        concat_fastqc_summary.py --summarized_fastqcs ~{sep=" "  summarized_fastqcs} \
            --project_name ~{project_name}
        echo $DOCKER_VERSION > VERSION
    >>>

    output {
        File fastqc_summary = "~{project_name}_reads_QC_summary.csv"
        VersionInfo version_info = {
            "software": "cdphe_h5_influenza docker",
            "docker": docker,
            "version": read_string("VERSION")
        }
    }

    runtime {
        #cpu: 
        #memory: 
        docker: docker
    }
}
