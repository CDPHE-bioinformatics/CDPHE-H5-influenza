version 1.0

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
    input {
        String docker
        String workflow_name
        String workflow_version
    }
    meta {
        description: "capture GitHub repository version"
        volatile: true
    }

    command <<<
        date +"%Y-%m-%d" > TODAY
    >>>

    output {
        String analysis_date = read_string("TODAY")
        VersionInfo version_info = {
            "software": workflow_name,
            "docker": "",
            "version": workflow_version
        }
    }

    runtime {
        docker: docker
    }
}

task capture_versions {
    meta {
        volatile: true
    }

    input {
        Array[VersionInfo] version_array
        String workflow_name
        String workflow_version
        String project_name
        String analysis_date
        String docker
    }

    VersionInfoArray versions = object {versions: version_array}
    String out_fn = "version_capture_~{workflow_name}_~{project_name}_~{workflow_version}.csv"

    command <<<
        cp ../../usr/src/app/* .
        python3 version_capture.py \
        --versions_json ~{write_json(versions)} \
        --workflow_name ~{workflow_name} \
        --workflow_version ~{workflow_version} \
        --project_name ~{project_name} \
        --analysis_date ~{analysis_date} \
        --docker_name $NAME \
        --docker_host $HOST  \
        --docker_version $VERSION  \
    >>>

    output {
        File output_file = out_fn
    }

    runtime {
        docker: docker
    }
}