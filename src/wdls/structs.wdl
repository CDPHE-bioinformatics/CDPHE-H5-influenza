version 1.0

struct Sample {
    String name
    String primer
    File fastq1
    File fastq2 
    Int i
}

struct PrimerScheme {
    String name
    File bed
    String reference_name
    File reference_fasta
}

struct PrimerSchemeArray {
    Array[PrimerScheme] primer_schemes
}

struct VersionInfo {
  String software
  String docker
  String version
}

workflow declare_structs {
    input {
        String h5_docker
    }

    call download_references {
        input:
            docker = h5_docker
    }

    output {
        Array[PrimerScheme] primer_schemes = download_references.primer_schemes.primer_schemes
    }
}

task download_references {
    input {
        String docker
    }

    command <<<
        cp $APPDIR/references/* .
        echo "Outputting references..."
    >>>

    output {
        PrimerSchemeArray primer_schemes = read_json("primer_scheme_structs.json")
    }

    runtime {
        docker: docker
    }
}