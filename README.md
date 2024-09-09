# DISCLAIMER
This repository is in active development and is not yet ready for production use.

## CDPHE-H5-Influenza workflow

### Disclaimer
Next generation sequencing and bioinformatic and genomic analysis at the Colorado Department of Public Health and Environment (CDPHE) is not CLIA validated at this time. These workflows and their outputs are not to be used for diagnostic purposes and should only be used for public health action and surveillance purposes. CDPHE is not responsible for the incorrect or inappropriate use of these workflows or their results.

### Overview
The following documentation describes the Colorado Department of Public Health and Environment's workflows for the assembly and analysis of next genome sequencing data of H5 influenza on GCP's Terra.bio platform. 

The workflow currently allows for various primer schemes and references for alignment. Some are only for the HA gene segment while others are whole genome. See Workspace data below.

The workflow is split into multiple wdl files, but is all launched as one workflow.

### Files

#### Subworkflows

| Name | Description | Subworkflow/task calls |
| --- | --- | --- |
| `h5_testing` | The main workflow that calls all other subworkflows and tasks. | `version_capture_tasks.workflow_metadata`<br>`h5_structs.declare_structs`<br>`primer_tasks.primer_level_tasks`<br>`reference_tasks.reference_level_tasks`<br>`version_capture_tasks.capture_versions`<br>`other_tasks.transfer_vc` |
| `h5_structs` | Contains struct definitions and subworkflow to declare the primer schemes and references structs. | None |
| `primer_tasks` | Subworkflow and task declarations for primer-level tasks.| `fastqc_raw`<br>`seqyclean`<br>`fastqc_clean`<br>`sample_qc_file`<br>`multiqc_fastqc`<br>`multiqc_seqyclean`<br>`transfer` |
| `reference_tasks` | Subworkflow and task declarations for reference-level tasks.| `align_bwa`<br>`trim_primers_ivar`<br>`generate_consensus_ivar`<br>`alignment_metrics`<br>`calculate_percent_coverage`<br>`concat_sample_reference_metrics`<br>`multiqc_samtools`<br>`transfer` |
| `version_capture_tasks` | Task declarations for tool and workflow version capture.| `workflow_metadata`<br>`capture_versions` |
| `other_tasks` | Task declarations for all other tasks. | `transfer`<br>`multiqc`<br>`concat_all_samples_metrics` |

### Flow
Note- Commonly-used inputs such as `sample_name` are not noted below.

- Call `version_capture_tasks.workflow_metadata` task to get the workflow version and analysis date.
- Call `h5_structs.declare_structs` subworkflow to declare reference and primer scheme struct objects.
- Scatter samples to create `Sample` struct objects.
- Scatter `PrimerScheme` array.
  - Create primer sample list, excluding those with empty fastq files.
  - Call `primer_tasks.primer_level_tasks` subworkflow
    - Scatter samples
      - Call `fastqc_raw`. Input - raw fastq files. 
      - Call `seqyclean`. Input - raw fastq files.
      - Call `fastqc_clean`. Input - cleaned fastq files from `seqyclean`.
    - Call `sample_qc_file`. Input - `data.txt` files from `fastqc_raw` and `fastqc_clean`, `SummaryStatistics.txt` from `seqyclean`.
    - Call `multiqc_fastqc`. Input- `data.txt` files from `fastqc_raw` and `fastqc_clean`.
    - Call `multiqc_seqyclean`. Input- `SummaryStatistics.txt` files from `seqyclean`.
    - Call `other_tasks.transfer` task to transfer all primer-level task outputs to their respective directories.
  - Call `reference_tasks.reference_level_tasks` subworkflow
    - Scatter samples
      - Call `align_bwa`. Input - `reference` name and fasta file, cleaned fastq files from `seqyclean`.
      - Call `trim_primers_ivar`. Input - `bwa` from `align_bwa`, `primer_bed` file
      - Call `generate_consensus_ivar`. Input - `trimmed.sorted.bam` from `trim_primers_ivar`, `reference` fasta file.
      - Call `alignment_metrics`. Input - `trimmed.sorted.bam` from `trim_primers_ivar`
      - Call `calculate_percent_coverage`. Input - {to fill in}
      - Call `concat_sample_reference_metrics`.  Input - `stats.txt` and `coverage.txt` from `alignment_metrics`, `coverage_stats.csv` from `calculate_percent_coverage`.
    - Call `multiqc_samtools`. Input - `stats.txt` and `coverage.txt` from `alignment_metrics`.
    - Call `other_tasks.transfer` task to transfer all reference-level task outputs to their respective directories.
- Call `version_capture_tasks.capture_versions`. Input:`analysis_date` and `workflow_version` from `workflow_metadata`, `Array[VersionInfo]` objects from all tools used.
- Call `other_tasks.transfer` task to transfer `version_capture.csv` from `version_capture_tasks.capture_versions`.

### Workspace data

| Name | Description |
| -- | -- |
| `AVRL_H5N1_250bp_bed` | cattle-specific, tiled whole genome |
| `houston_bed` | Rice, H1, H3, H5, N1, N2 - Olivar (tiled) |
| `human_h5_200_bed` | cattle-H5-specific, 200bp, tiled HA gene |
| `human_h5_250_bed` | cattle-H5-specific, 250bp, tiled HA gene |
| `bovine_texas_029328_01_UtoT_fasta` | Cattle reference, whole genome |
| `bovine_texas_029328_01_UtoT_ha_fasta` | Cattle reference, HA gene |
| `darwin_9_2021_h3n2_ha_h3_fasta` | CDC vaccine strain reference for H3N2, HA gene |
| `victoria_4897_2022_h1n1_ha_h1_fasta` | CDC vaccine strain reference for H1N1, HA gene |
| `vietnam_1203_2024_h5n1_ha_v2_fasta` | CDC vaccine strain reference for H5N1, HA gene |
| `contaminants_fasta` | Adapters and contaminants fasta file |