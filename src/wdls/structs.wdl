version development

struct Sample {
    String name
    String primer
    File fastq1
    File fastq2 
    Int i
}

struct Reference {
    String name
    File fasta
}

struct PrimerScheme {
    String name
    Array[Reference] references
    File bed
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

    call download_references as refs{
        input:
            docker = h5_docker
    }

    Reference A_Bovine_Texas_24-029328-01_2024_H5N1_multi = Reference {
        name: "A_Bovine_Texas_24-029328-01_2024_H5N1_multi",
        fasta: refs.A_Bovine_Texas_24-029328-01_2024_H5N1_multi_fasta
    }
    Reference houston_fluA_multi = Reference {
        name: "houston_fluA_multi",
        fasta: refs.houston_fluA_multi_fasta
    }
    PrimerScheme human_h5_200 = PrimerScheme {
        name: "human_h5_200",
        references:   [A_Texas_37_2024_H5N1_HA-H5],
        bed: refs.human_h5_200_bed
    }
    PrimerScheme human_h5_250 = PrimerScheme {
        name: "human_h5_250",
        references:   [A_Texas_37_2024_H5N1_HA-H5],
        bed: refs.human_h5_250_bed
    }
    PrimerScheme houston = PrimerScheme {
        name: "houston",
        references:  [H1N1-2024-01-20-H5N1-2022-05-25-NA-consensus],
        bed: refs.houston_bed
    }
    PrimerScheme AVRL_H5N1_250bp = PrimerScheme {
        name: "AVRL_H5N1_250bp",
        references: [A_Bovine_Texas_24-029328-01_2024_H5N1_multi],
        bed: refs.AVRL_H5N1_250bp_bed
    }

    output {
        Array[PrimerScheme] primer_schemes = [human_h5_200, human_h5_250, houston, AVRL_H5N1_250bp] 
        Array[Reference] references = [A_Texas_37_2024_H5N1_HA-H5, 
                                        A_Bovine_Texas_24-029328-01_2024_H5N1_multi,
                                        houston_fluA_multi]
    }
}

task download_references {
    input {
        String docker
    }

    output {
        File A_Texas_37_2024_H5N1_HA-H5_fasta = "references/A_Texas_37_2024_H5N1_HA-H5.fasta"
        File A_Bovine_Texas_24-029328-01_2024_H5N1_multi_fasta = "references/A_Bovine_Texas_24-029328-01_2024_H5N1_multi.fasta"
        File houston_fluA_multi_fasta = "references/houston_fluA_multi.fasta"
        File human_h5_200_bed = "references/human_h5_200bp_primer.bed"
        File human_h5_250_bed = "references/human_h5_250bp_primer.bed"
        File houston_bed = "references/houston_fluA_primer.bed"
        File AVRL_H5N1_250bp_bed = "references/AVRL_H5N1_250bpAmpWGS_primer_v2"
    }

    runtime {
        docker: docker
    }
}