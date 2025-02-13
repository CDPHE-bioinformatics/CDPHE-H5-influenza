version 1.0

import "h5_structs.wdl"

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
    String base_cl_config = '--cl-config "show_analysis_paths: False'

    command <<<
        multiqc -m "~{module}" -l ~{write_lines(files)} -n ~{html_fn} \
        ~{if defined(cl_config) then '~{base_cl_config}\n~{cl_config}"' else '~{base_cl_config}"'}  

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
