version 1.0

import "structs.wdl"

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

    String prefix = if defined(task_name) then "~{module}_~{task_name}" else "~{module}"
    String html_fn = "~{prefix}_multiqc_report.html"
    String base_cl_config = '--cl-config "show_analysis_paths: False'

    command <<<
        multiqc -m "~{module}" -l ~{write_lines(files)} -n ~{html_fn} \
        ~{if defined(cl_config) then '~{base_cl_config}\n~{cl_config}"' else '~{base_cl_config}"'}  

        multiqc --version | awk ' {print $3}' | tee VERSION
    >>>

    output {
        File html_report = html_fn
        VersionInfo version_info = {
            "software": "multiqc",
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

task concat_all_samples_metrics {
    input {
        String project_name
        Array[File] segment_metrics_files
        Array[File] sample_metrics_files
        String docker
    }

    String segment_concat_fn = "~{project_name}_segment_metrics_summary.csv"
    String sample_concat_fn = "~{project_name}_sample_metrics_summary.csv"

    command <<<
        awk 'NR==1||FNR>1' ~{sep=" " segment_metrics_files} > ~{segment_concat_fn}
        awk 'NR==1||FNR>1' ~{sep=" " sample_metrics_files}> ~{sample_concat_fn}
    >>>

    output {
        File segment_summary = segment_concat_fn
        File sample_summary = sample_concat_fn
    }

    runtime {
        docker: docker
    }
}
