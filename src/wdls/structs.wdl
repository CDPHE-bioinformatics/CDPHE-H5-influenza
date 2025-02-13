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

    call download_references as refs {
        input:
            docker = h5_docker
    }

    Reference A_Bovine_Texas_24-029328-01_2024_H5N1_multi = Reference {
        name: "A_Bovine_Texas_24-029328-01_2024_H5N1_multi",
        fasta: refs.A_Bovine_Texas_24-029328-01_2024_H5N1_multi_fasta
    }
    Reference A_Texas_37_2024_H5N1_HA-H5 = Reference {
        name: "A_Texas_37_2024_H5N1_HA-H5_fasta",
        fasta: refs.A_Texas_37_2024_H5N1_HA-H5_fasta
    }
    Reference houston_fluA_multi = Reference {
        name: "houston_fluA_multi_fasta",
        fasta: refs.houston_fluA_multi_fasta
    } 
    PrimerScheme AVRL_H5N1_250bp = PrimerScheme {
        name: "AVRL_H5N1_250bp",
        reference: A_Bovine_Texas_24-029328-01_2024_H5N1_multi,
        bed: refs.AVRL_H5N1_250bp_bed
    }
    PrimerScheme houston = PrimerScheme {
        name: "houston",
        reference:  houston_fluA_multi,
        bed: refs.houston_bed
    }
    PrimerScheme human_h5_200 = PrimerScheme {
        name: "human_h5_200",
        reference:   A_Texas_37_2024_H5N1_HA-H5,
        bed: refs.human_h5_200_bed
    }
    PrimerScheme human_h5_250 = PrimerScheme {
        name: "human_h5_250",
        reference:   A_Texas_37_2024_H5N1_HA-H5,
        bed: refs.human_h5_250_bed
    }


    output {
        Array[PrimerScheme] primer_schemes = [AVRL_H5N1_250bp, houston, human_h5_200, human_h5_250] 
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
        File A_Bovine_Texas_24-029328-01_2024_H5N1_multi_fasta = "references/A_Bovine_Texas_24-029328-01_2024_H5N1_multi.fasta"
        File A_Texas_37_2024_H5N1_HA-H5_fasta = "references/A_Texas_37_2024_H5N1_HA-H5.fasta"
        File houston_fluA_multi_fasta = "references/houston_fluA_multi.fasta"
        File houston_bed = "references/houston_fluA_primer.bed"
        File human_h5_200_bed = "references/human_h5_200bp_primer.bed"
        File human_h5_250_bed = "references/human_h5_250bp_primer.bed"
        File AVRL_H5N1_250bp_bed = "references/AVRL_H5N1_250bpAmpWGS_primer_v2.bed"
    }

    runtime {
        docker: docker
    }
}