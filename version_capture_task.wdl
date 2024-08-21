version development

struct VersionInfo {
    String software
    String docker
    String version
}

# workaround cromwell bug with read_json of Array
# https://github.com/openwdl/wdl/issues/409
struct VersionInfoArray {
  Array[VersionInfo] versions
}

task workflow_metadata {

    meta {
        description: "capture repository version release"
    }

    command <<<
        date +"%Y-%m-%d" > TODAY
    >>>

    output {
        String analysis_date = read_string("TODAY")
        String version = "v0-0-0-alpha" 
    }
}

task capture_versions {
    input {
        Array[VersionInfo] version_array
        String workflow_name
        String workflow_version
        String project_name
        String project_outdir
        String analysis_date
        File version_capture_py
        String docker
    }

    VersionInfoArray versions = VersionInfoArray { versions: version_array }
    String out_fn = "version_capture_~{workflow_name}_~{project_name}_~{workflow_version}.csv"

    command <<<
        python ~{version_capture_py} \
        --versions_json ~{write_json(versions)} \
        --workflow_name ~{workflow_name} \
        --workflow_version ~{workflow_version} \
        --project_name ~{project_name} \
        --analysis_date ~{analysis_date}

        gsutil cp ~{out_fn} "~{project_outdir}summary_results/"
    >>>

    output {
        File output_file = out_fn
    }

    runtime {
        docker: docker
    }
}