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
        File bovine_texas_029328_01_UtoT_ha_fasta
        File bovine_texas_029328_01_UtoT_fasta
        File darwin_9_2021_h3n2_ha_h3_fasta
        File victoria_4897_2022_h1n1_ha_h1_fasta
        File vietnam_1203_2024_h5n1_ha_v2_fasta  
        File human_h5_200_bed
        File human_h5_250_bed
        File houston_bed
        File AVRL_H5N1_250bp_bed
    }
    Reference vietnam_1203_2024_h5n1_ha_v2 = Reference {
        name: "vietnam_1203_2024_h5n1_ha_v2",
        fasta: vietnam_1203_2024_h5n1_ha_v2_fasta
    }
    Reference bovine_texas_029328_01_UtoT_ha = Reference {
        name: "bovine_texas_029328_01_UtoT_ha",
        fasta: bovine_texas_029328_01_UtoT_ha_fasta
    }
    Reference bovine_texas_029328_01_UtoT = Reference {
        name: "bovine_texas_029328_01_UtoT",
        fasta: bovine_texas_029328_01_UtoT_fasta
    }
    Reference darwin_9_2021_h3n2_ha_h3 = Reference {
        name: "darwin_9_2021_h3n2_ha_h3",
        fasta: darwin_9_2021_h3n2_ha_h3_fasta
    }
    Reference victoria_4897_2022_h1n1_ha_h1 = Reference {
        name: "victoria_4897_2022_h1n1_ha_h1",
        fasta: victoria_4897_2022_h1n1_ha_h1_fasta
    }

    PrimerScheme human_h5_200 = PrimerScheme {
        name: "human_h5_200",
        references:   [vietnam_1203_2024_h5n1_ha_v2, 
                        bovine_texas_029328_01_UtoT_ha],
        bed: human_h5_200_bed
    }
    PrimerScheme human_h5_250 = PrimerScheme {
        name: "human_h5_250",
        references:   [vietnam_1203_2024_h5n1_ha_v2, 
                            bovine_texas_029328_01_UtoT_ha],
        bed: human_h5_250_bed
    }
    PrimerScheme houston = PrimerScheme {
        name: "houston",
        references:  [darwin_9_2021_h3n2_ha_h3, 
                            victoria_4897_2022_h1n1_ha_h1, 
                            bovine_texas_029328_01_UtoT_ha],
        bed: houston_bed
    }
    PrimerScheme AVRL_H5N1_250bp = PrimerScheme {
        name: "AVRL_H5N1_250bp",
        references: [bovine_texas_029328_01_UtoT],
        bed: AVRL_H5N1_250bp_bed
    }

    output {
        Array[PrimerScheme] primer_schemes = [human_h5_200, human_h5_250, houston, AVRL_H5N1_250bp] 
        Array[Reference] references = [vietnam_1203_2024_h5n1_ha_v2, 
                                        bovine_texas_029328_01_UtoT_ha, 
                                        darwin_9_2021_h3n2_ha_h3, 
                                        victoria_4897_2022_h1n1_ha_h1,
                                        bovine_texas_029328_01_UtoT]
    }
}