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
        File A_Bovine_Texas_24_029328_01_2024_H5N1_multi_fasta = "A_Bovine_Texas_24-029328-01_2024_H5N1_multi.fasta"
        File A_Texas_37_2024_H5N1_HA_H5_fasta = "A_Texas_37_2024_H5N1_HA-H5.fasta"
        File houston_fluA_multi_fasta = "houston_fluA_multi.fasta"
        File houston_bed = "houston_fluA_primer.bed"
        File human_h5_200_bed = "human_h5_200bp_primer.bed"
        File human_h5_250_bed = "human_h5_250bp_primer.bed"
        File AVRL_H5N1_250bp_bed = "AVRL_H5N1_250bpAmpWGS_primer_v2.bed"
    }

    runtime {
        docker: docker
    }
}