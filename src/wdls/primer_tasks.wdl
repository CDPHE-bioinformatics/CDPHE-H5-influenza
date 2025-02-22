version 1.0

import "other_tasks.wdl" as ot

workflow primer_level_tasks {
    input {
        Array[Sample] primer_samples
        String primer_name
        String primer_outdir
        String project_name
        String fastqc_docker
        String fastp_docker
        String multiqc_docker
        String utility_docker
        String h5_docker
    }

    scatter (sample in primer_samples) {
        String sample_name = sample.name
        String fastq_raw = "raw"
        String fastq_clean = "clean"
        # Call sample level tasks
        call fastqc as fastqc_raw {
            input: 
                sample_name = sample_name,
                project_name = project_name,
                primer_name = primer_name,
                fastq1 = sample.fastq1, 
                fastq2 = sample.fastq2, 
                fastq_type = fastq_raw,
                docker = fastqc_docker
        }
        

        call fastp {
            input:
                sample = sample,
                docker = fastp_docker
        }

        call fastqc as fastqc_clean {
            input: 
                sample_name = sample_name,
                project_name = project_name,
                primer_name = primer_name,
                fastq1 = fastp.fastq_1_cleaned, 
                fastq2 = fastp.fastq_2_cleaned, 
                fastq_type = fastq_clean,
                docker = fastqc_docker
        }

    }

    call concat_fastqc_summary {
        input:
            raw_summarized_fastqcs = fastqc_raw.summary_metrics, 
            clean_summarized_fastqcs = fastqc_clean.summary_metrics,
            primer_name = primer_name,
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

    call ot.multiqc as multiqc_fastp {
        input:
            files = fastp.fastp_json,
            module = "fastp",
            docker = multiqc_docker
    }
    
    # Transfer primer level files
    Array[File] fastqc_raw_output = flatten([fastqc_raw.fastqc1_data, fastqc_raw.fastqc2_data])
    Array[File] fastqc_clean_output = flatten([fastqc_clean.fastqc1_data, fastqc_clean.fastqc2_data])
    Array[File] fastp_output = flatten([fastp.fastq_1_cleaned, fastp.fastq_2_cleaned])
    Array[File] p_summary_output = [multiqc_fastqc.html_report, multiqc_fastp.html_report, concat_fastqc_summary.fastqc_summary]

    Array[String] primer_task_dirs = ["fastqc_raw", "fastqc_clean", "fastp", "summary_results"]
    Array[Array[File]] primer_task_files = [fastqc_raw_output, fastqc_clean_output, fastp_output, p_summary_output]       

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
        Array[File] cleaned_PE1 = fastp.fastq_1_cleaned
        Array[File] cleaned_PE2 = fastp.fastq_2_cleaned
        Array[File] p_summary_outputs = p_summary_output
        Array[File] fastqc_clean_summary_metrics = fastqc_clean.summary_metrics
        File fastqc_summary = concat_fastqc_summary.fastqc_summary
        VersionInfo fastqc_version = select_first(fastqc_raw.version_info)
        VersionInfo fastp_version = select_first(fastp.version_info)
        VersionInfo multiqc_version = multiqc_fastqc.version_info
        VersionInfo h5_docker_version = concat_fastqc_summary.version_info
    }
}

task fastqc {
    input {
        String sample_name
        String project_name
        String primer_name
        File fastq1
        File fastq2
        String fastq_type
        String docker
    }

    String fastq1_name = basename(fastq1, ".fastq.gz")
    String fastq2_name = basename(fastq2, ".fastq.gz")
    String fastq1_data_name = "~{fastq1_name}_fastqc_data.txt"
    String fastq2_data_name = "~{fastq2_name}_fastqc_data.txt"
    String summary_metrics_fn = "~{sample_name}_~{fastq_type}_summary_metrics.tsv"

    command <<<
        fastqc --outdir $PWD --extract --delete ~{fastq1} ~{fastq2}
        fastqc --version | awk "{print $NF}" | tee VERSION  
        cp "~{fastq1_name}_fastqc/fastqc_data.txt" ~{fastq1_data_name}
        cp "~{fastq2_name}_fastqc/fastqc_data.txt" ~{fastq2_data_name}

        # Summarize output to simpler csv file
        summarize_fastqc () {
            total_seqs=$(grep "Total Sequences" $1 | cut -f 2)
            flagged_reads=$(grep "Sequences flagged as poor quality" $1 | cut -f 2)
            sequence_length=$(grep "Sequence length" $1 | cut -f 2)
            echo $total_seqs,$flagged_reads,$sequence_length
        }

        fastq1_data_name="~{fastq1_data_name}"
        fastq2_data_name="~{fastq2_data_name}"
        r1_info=$(summarize_fastqc ${fastq1_data_name})
        r2_info=$(summarize_fastqc ${fastq2_data_name})
        echo "sample_name,project_name,primer_name,r1_total_reads,r1_flagged_reads_as_poor_quality,r1_read_len,r2_total_reads,r2_flagged_reads_as_poor_quality,r2_read_len" >> ~{summary_metrics_fn} 
        echo "~{sample_name},~{project_name},~{primer_name},${r1_info},${r2_info}" >> ~{summary_metrics_fn}
    >>>

    output {
        File fastqc1_data = fastq1_data_name
        File fastqc2_data = fastq2_data_name
        File summary_metrics = summary_metrics_fn
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

task fastp {
    input {
        Sample sample
        String docker
    }

    String fastq1_name = basename(sample.fastq1, ".fastq.gz")
    String fastq2_name = basename(sample.fastq2, ".fastq.gz")
    String cleaned_1_name = "~{fastq1_name}_clean.fastq.gz"
    String cleaned_2_name = "~{fastq2_name}_clean.fastq.gz"
    String unpaired_1_name = "~{fastq1_name}_unpaired.fastq.gz"
    String unpaired_2_name = "~{fastq2_name}_unpaired.fastq.gz"
    String fastp_html_name = "~{sample.name}_fastp.html"
    String fastp_json_name = "~{sample.name}_fastp.json"

    command <<<
        fastp --version | tee VERSION
        fastp \
            --in1 ~{sample.fastq1} --in2 ~{sample.fastq2} \
            --out1 ~{cleaned_1_name} --out2 ~{cleaned_2_name}  \
            --unpaired1 ~{unpaired_1_name} --unpaired2 ~{unpaired_2_name} \
            --cut_tail \
            --cut_tail_window_size 4 \
            --cut_tail_mean_quality 30 \
            --length_required 70 \
            --detect_adapter_for_pe \
            --trim_poly_g \
            --html ~{fastp_html_name} \
            --json ~{fastp_json_name}
    >>>

    output {
        File fastq_1_cleaned = cleaned_1_name
        File fastq_2_cleaned = cleaned_2_name
        File fastq_1_unpaired = unpaired_1_name
        File fastq_2_unpaired = unpaired_2_name
        File fastp_html = fastp_html_name
        File fastp_json = fastp_json_name

        VersionInfo version_info = {
            "software": "fastp",
            "docker": docker,
            "version": read_string("VERSION")
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
        Array[File] clean_summarized_fastqcs
        Array[File] raw_summarized_fastqcs
        String primer_name
        String project_name
        String docker
    }

    meta {
        volatile: true
    }

    String out_fn = "~{project_name}_~{primer_name}_reads_QC_summary.csv"

    command <<<
        concat_fastqc_summary.py \
            --raw_summarized_fastqcs ~{sep=" " raw_summarized_fastqcs} \
            --clean_summarized_fastqcs ~{sep=" " clean_summarized_fastqcs} \
            --out_fn ~{out_fn} 
        echo $DOCKER_VERSION > VERSION
    >>>

    output {
        File fastqc_summary = out_fn
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
