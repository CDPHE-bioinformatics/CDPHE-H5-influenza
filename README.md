# DISCLAIMER
This repository is in active development and is not yet ready for production use.

# CDPHE-H5-Influenza workflow

## Disclaimer
Next generation sequencing and bioinformatic and genomic analysis at the Colorado Department of Public Health and Environment (CDPHE) is not CLIA validated at this time. These workflows and their outputs are not to be used for diagnostic purposes and should only be used for public health action and surveillance purposes. CDPHE is not responsible for the incorrect or inappropriate use of these workflows or their results.

## Overview
The following documentation describes the Colorado Department of Public Health and Environment's workflows for the assembly and analysis of next genome sequencing data of H5 influenza on GCP's Terra.bio platform. 

The workflow currently allows for various primer schemes and references for alignment. Some are only for certain gene segments while others are whole genome. See Workspace data below.

The workflow is split into multiple WDL files, but is all launched as one workflow.

## WDL files

### Main workflow

The main workflow, `h5_assembly_analysis`, is a set-level workflow that calls all other subworkflows and tasks.

#### Inputs

| Variable | Type | Description |
| -- | -- | -- |
| `fastq1s` | `Array[File]` | R1 fastq files |
| `fastq2s` | `Array[File]` | R2 fastq files |
| `out_dir` | `String` | Directory prefix to copy files to |
| `primers` | `Array[String]` | Primer schemes - names must match those in `structs.wdl` |
| `project_name` | `String` | Sequencing run name |
| `sample_names` | `Array[String]` | Sample names - must be unique |

### Subworkflows

| Name | Description | Task definitions |
| --- | --- | --- |
| `structs` | Contains struct definitions and subworkflow to download primer schemes structs. | `download_references` |
| `primer_tasks` | Subworkflow and task declarations for primer-level tasks.| `fastqc`<br>`fastp`<br>`concat_fastqc_summary` |
| `reference_tasks` | Subworkflow and task declarations for reference-level tasks.| `align_bwa`<br>`calculate_alignment_metrics`<br>`calculate_metrics_samtools`<br>`concat_sample_reference_metrics`<br>`generate_consensus_ivar`<br>`trim_primers_samtools` |
| `version_capture_tasks` | Task declarations for tool and workflow version capture.| `capture_versions`<br>`workflow_metadata` |
| `other_tasks` | Task declarations for tasks used by multiple workflows. | `concat_all_samples_metrics`<br>`multiqc`<br>`transfer` |

### High-level overview
- Call `version_capture_tasks.workflow_metadata`
- Call `structs.declare_structs` subworkflow 
- Create `Sample` structs.
- Scatter `PrimerScheme` array.
  - Create primer sample list, excluding those with empty fastq files.
  - Call `primer_tasks.primer_level_tasks` subworkflow
    - Scatter samples
      - Call `fastqc_raw`
      - Call `fastp`
      - Call `fastqc_clean`
    - Call `concat_fastqc_summary`
    - Call `multiqc_fastqc`
    - Call `multiqc_fastp`
    - Call `other_tasks.transfer`
  - Call `reference_tasks.reference_level_tasks` subworkflow
    - Scatter samples
      - Call `align_bwa`
      - Call `trim_primers_samtools`
      - Call `generate_consensus_ivar`
      - Call `calculate_metrics_samtools`
    - Call `calculate_alignment_metrics`
    - Call `multiqc_samtools`
    - Call `other_tasks.transfer`
- Call `other_tasks.concat_all_samples_metrics`
- Call `other_tasks.transfer_concat_metrics`
- Call `version_capture_tasks.capture_versions`
- Call `other_tasks.transfer_vc` 

## Docker container

Located at [ariannaesmith/cdphe_h5_influenza](https://hub.docker.com/repository/docker/ariannaesmith/cdphe_h5_influenza/general).

### Reference data

| Primer bed | Associated reference fasta | Description | Source |
| -- | -- | -- | -- 
| `AVRL_H5N1_250bpAmpWGS_primer_v2.bed` | `A_Bovine_Texas_24-029328-01_2024_H5N1_multi.fasta` | Tiled - whole genome - H5N1 | [Paper](dx.doi.org/10.17504/protocols.io.kqdg322kpv25/v1) |
| `houston_fluA_primer.bed` | `houston_fluA_multi.fasta` | Tiled - HA + NA genes from H5N1 and H3N2; HA gene from H1N1 | [Github repository](https://github.com/treangenlab/InfA-amplicon) |
| `human_h5_200bp_primer.bed` | `A_Texas_37_2024_H5N1_HA-H5.fasta` | Tiled - HA gene - H5N1 | made in-house  with [PrimalScheme](https://primalscheme.com/) |
| `human_h5_250bp_primer.bed` | `A_Texas_37_2024_H5N1_HA-H5.fasta` | Tiled - HA gene - H5N1 | made in-house with [PrimalScheme](https://primalscheme.com/) |
| `h5_PB2_primer.bed` | `h5_PB2_consensus.fasta` | Tiled - PB2 gene - H5N1 (B3.13 and D.1.1)| made in-house with [PrimalScheme](https://primalscheme.com/) |
| `olivar_A_HA-H1_primer.bed` | `A_Victoria_4897_2022_H1N1_HA-H1.fasta` | Tiled - HA gene - H1N1 | Made in-house with [Olivar](https://github.com/treangenlab/Olivar) |
| `olivar_A_HA-H3_primer.bed` | `A_Darwin_9_2021_H3N2_HA-H3.fasta` | Tiled - HA gene - H3N2 | Made in-house with [Olivar](https://github.com/treangenlab/Olivar) |
| `olivar_B_HA_primer.bed` | `B_Austria_1359417_2021_vic_HA.fasta` | Tiled - HA gene - B/Victoria | Made in-house with [Olivar](https://github.com/treangenlab/Olivar) |

### Python scripts

| Name | Associated workflow | Description |
| -- | -- | -- |
| `calculate_alignment_metrics.py` | `reference_tasks` | Calculate alignment metrics |
| `concat_fastqc_summary.py` | `primer_tasks` | Calculate QC metrics | 