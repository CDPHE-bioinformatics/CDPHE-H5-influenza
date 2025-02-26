version 1.0

import "other_tasks.wdl" as ot

workflow reference_level_tasks {
    input {
        String reference_name
        File reference_fasta
        String project_name
        String reference_outdir
        Array[Int] num_samples
        Array[Sample] primer_samples
        Array[File] cleaned_PE1
        Array[File] cleaned_PE2
        File reads_qc_summary
        File primer_bed
        String ivar_docker
        String multiqc_docker
        String utility_docker
        String h5_docker
    }

    call sort_bed {
        input: 
            primer_bed = primer_bed,
            docker = h5_docker
    }

    scatter (n in num_samples) {
        Sample sample = primer_samples[n]
        String sample_name = sample.name
        File PE1 = cleaned_PE1[n]
        File PE2 = cleaned_PE2[n]

        call align_bwa {
            input:
                sample_name = sample_name,
                fastq1 = PE1,
                fastq2 = PE2,
                reference_name = reference_name,
                reference_fasta = reference_fasta,
                docker = ivar_docker
        }  

        call trim_primers_samtools {
            input: 
                sample_name = sample_name,
                aligned_bam = align_bwa.bam,
                sorted_bed = sort_bed.sorted_bed,
                docker = ivar_docker
        }

        call generate_consensus_ivar {
            input: 
                sample_name = sample_name,
                trim_sort_bam = trim_primers_samtools.trim_sort_bam,
                reference_fasta = reference_fasta,
                docker = ivar_docker
        }   

        call calculate_metrics_samtools {
            input:
                sample_name = sample_name,
                trim_sort_bam = trim_primers_samtools.trim_sort_bam,
                docker = ivar_docker
        }

    } 

    call calculate_alignment_metrics {
        input:
            sample_names = sample_name,
            consensus_fastas = generate_consensus_ivar.consensus_fasta,
            samtools_coverages = calculate_metrics_samtools.coverage,
            mapped_reads = calculate_metrics_samtools.mapped_reads,
            reads_qc_summary = reads_qc_summary,
            project_name  = project_name,
            reference_fasta = reference_fasta,
            reference_name = reference_name,
            primer_bed = primer_bed,
            docker = h5_docker
    }

    call ot.multiqc as multiqc_samtools {
        input:
            files = flatten([trim_primers_samtools.idxstats, calculate_metrics_samtools.coverage, calculate_metrics_samtools.stats]),
            module = "samtools",
            task_name = "alignment",
            cl_config = "extra_fn_clean_exts: [ '_coverage', '_stats']",
            docker = multiqc_docker
    }
    
    # Transfer reference level files
    Array[File] alignment_output = flatten([trim_primers_samtools.trim_sort_bam, trim_primers_samtools.trim_sort_bai])
    Array[File] consensus_output = generate_consensus_ivar.consensus_fasta
    Array[File] samtools_output = flatten([calculate_metrics_samtools.coverage, 
                                            calculate_metrics_samtools.stats])
    Array[File] summary_output = [multiqc_samtools.html_report, 
                                    calculate_alignment_metrics.percent_coverage,
                                    calculate_alignment_metrics.aln_metrics]
    
    Array[String] reference_task_dirs = ["alignments", "assemblies", "samtools", "summary_results"]
    Array[Array[File]] reference_task_files = [alignment_output, consensus_output, samtools_output, summary_output]

    scatter (dir_files in zip(reference_task_dirs, reference_task_files)) {       
        call ot.transfer {
            input:
                out_dir = reference_outdir,
                task_dir = dir_files.left,
                task_files = dir_files.right,
                docker = utility_docker
        }
    }

    output {
        VersionInfo samtools_version = select_first(align_bwa.samtools_version_info)
        VersionInfo bwa_version = select_first(align_bwa.bwa_version_info)
        VersionInfo ivar_version = select_first(align_bwa.ivar_version_info)
        Array[File] alignment_outputs = alignment_output
        Array[File] consensus_outputs = consensus_output
        Array[File] summary_outputs = summary_output   
    }
}

task align_bwa {
    input {
        String sample_name
        File fastq1
        File fastq2
        String reference_name
        File reference_fasta
        String docker
    }

    Float input_size = size([fastq1, fastq2], "GiB")
    Int dynamic_disk_size = ceil((input_size) * 10) + ceil(input_size) + 10
    String sam_fn = "~{sample_name}.sam"
    String bam_fn = "~{sample_name}_aln.sorted.bam"

    command <<<
        samtools --version-only | tee SAMTOOLS_VERSION
        ivar version | awk ' /iVar/ {print $3}' | tee IVAR_VERSION
        bwa index -p ~{reference_name} -a is ~{reference_fasta}
        bwa mem -t 4 ~{reference_name} ~{fastq1} ~{fastq2} > ~{sam_fn}
        samtools view -b -@ 4 ~{sam_fn} | samtools sort -m 2G -@ 4 -o ~{bam_fn}
    >>>

    output {
        File bam = bam_fn
        VersionInfo bwa_version_info = {
            "software": "bwa",
            "docker": docker,
            "version": '0.7.17-r1188'
        }
        VersionInfo samtools_version_info = {
            "software": "samtools",
            "docker": docker,
            "version": read_string('SAMTOOLS_VERSION')
        }
        VersionInfo ivar_version_info = {
            "software": "ivar",
            "docker": docker,
            "version": read_string('IVAR_VERSION')
        }
    }

    runtime {
        cpu: 4
        memory: "15 GB"
        disks: "local-disk ~{dynamic_disk_size} HDD"
        docker: docker
        maxRetries: 2
    }
}

task sort_bed {
    input {
        File primer_bed
        String docker
    }

    String out_fn = basename(primer_bed, ".bed") + "_sorted.bed"

    command <<<
        sort_bed.py --bed_file ~{primer_bed} --out_fn ~{out_fn}
    >>>

    output {
        File sorted_bed = out_fn
    }

    runtime {
        docker: docker
    }
}

task trim_primers_samtools {
    input {
        String sample_name
        File sorted_bed
        File aligned_bam
        String docker
    }

    Int dynamic_disk_size = ceil(size(aligned_bam, "GiB")) * 2 + 10
    String trim_bam_fn = "~{sample_name}_trimmed.bam"
    String trim_sort_bam_fn = "~{sample_name}_trimmed.sorted.bam"
    String trim_sort_bai_fn = "~{sample_name}_trimmed.sorted.bai"
    String idxstats_fn = "~{sample_name}_trimmed_sorted_idxstats.txt"

    command <<<
        aligned_bam=~{aligned_bam}
        sorted_bed=~{sorted_bed}
        sn_segments=$(samtools view -H $aligned_bam | grep '^@SQ' | cut -f 2)
        segments=$(for s in ${sn_segments}; do echo "${s#SN:}"; done) # test this bash line again
        bed_segments=$(cat ${sorted_bed} | cut -f 1 | uniq)

        for elem in ${segments[@]}; do 
            if [[ " ${bed_segments[*]} " =~ [[:space:]]${elem}[[:space:]] ]]; then 
                echo "${elem} found in bed file"
            else
                echo "The reference ${elem} from the sample's BAM file count not be matched in the BED file"
                exit 1
            fi
        done

        samtools ampliconclip --both-ends -b ~{sorted_bed} -o ~{trim_bam_fn} ~{aligned_bam}
        samtools sort ~{trim_bam_fn} -o ~{trim_sort_bam_fn}
        samtools index ~{trim_sort_bam_fn}
        samtools idxstats ~{trim_sort_bam_fn} > ~{idxstats_fn}
    >>>

    output {
        File trim_sort_bam = trim_sort_bam_fn
        File trim_sort_bai = trim_sort_bai_fn
        File idxstats = idxstats_fn 
    }

    runtime {
        cpu: 4
        memory: "15 GB"
        disks: "local-disk ~{dynamic_disk_size} HDD"
        docker: docker
        maxRetries: 2
    }
}

task generate_consensus_ivar {
    input {
        String sample_name
        File trim_sort_bam
        File reference_fasta
        String docker
    }

    Int dynamic_disk_size = ceil(size(trim_sort_bam, "GiB")) * 2 + 10
    String pileup_fn = "~{sample_name}_pileup.txt"
    String consensus_fn_prefix = "~{sample_name}_consensus"
    String consensus_fn = consensus_fn_prefix + ".fa"

    command <<<
        samtools faidx ~{reference_fasta}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ~{reference_fasta} ~{trim_sort_bam} -o ~{pileup_fn}

        reference_fasta=~{reference_fasta}
        num_records=$(grep -c ">" $reference_fasta)
        if (( num_records > 1 )); then 
            awk -F '|' '/^>/ {close(F); ID=$1; gsub("^>", "", ID); gsub(" $", "", ID); F=ID".fasta"; array+="$ID"} {print > F}' ~{reference_fasta}
            rm ~{reference_fasta}
            files=(*.fasta)
            for filename in "${files[@]}"; do
                segment_id="${filename%'.fasta'}"
                prefix="~{sample_name}_${segment_id}_consensus"
                cat ~{pileup_fn} | grep ${segment_id} | ivar consensus -p ${prefix} -q 20 -t 0.6 -m 10
            done
            cat *consensus.fa > ~{consensus_fn}
        else
            cat ~{pileup_fn} | ivar consensus -p ~{consensus_fn_prefix} -q 20 -t 0.6 -m 10
        fi

    >>>

    output {
        File pileup = pileup_fn
        File consensus_fasta = consensus_fn
        Array[File]? segment_consensus_fastas = glob("~{sample_name}_*_consensus.fa")
    }

    runtime {
        cpu: 4
        memory: "15 GB"
        disks: "local-disk ~{dynamic_disk_size} HDD"
        docker: docker
        maxRetries: 2
    }
}

task calculate_metrics_samtools {
    input {
        String sample_name
        File trim_sort_bam
        String docker
    }

    String coverage_fn = "~{sample_name}_coverage.txt"
    String stats_fn = "~{sample_name}_stats.txt"

    command <<<
        samtools coverage -o ~{coverage_fn} ~{trim_sort_bam}
        samtools stats ~{trim_sort_bam} > ~{stats_fn}
        grep "reads mapped:" ~{stats_fn} | cut -f 3 >> MAPPED_READS
    >>>

    output {
        File coverage = coverage_fn
        File stats = stats_fn
        String mapped_reads = read_string('MAPPED_READS')
    }

    runtime {
        cpu: 2
        memory: "4 GB"
        docker: docker
    }
}

task calculate_alignment_metrics {
    input {
        Array[String] sample_names
        Array[File] consensus_fastas
        Array[File] samtools_coverages
        Array[String] mapped_reads
        File reads_qc_summary
        String project_name
        File reference_fasta
        String reference_name
        File primer_bed
        String docker
    }

    meta {
        volatile: true
    }

    command <<<
        calculate_alignment_metrics.py \
                --sample_names ~{sep=" " sample_names} \
                --consensus_fastas ~{sep=" " consensus_fastas} \
                --samtools_coverages ~{sep=" " samtools_coverages} \
                --mapped_reads ~{sep=" " mapped_reads} \
                --reads_qc_summary ~{reads_qc_summary} \
                --project_name ~{project_name} \
                --reference_fasta ~{reference_fasta} \
                --reference_name ~{reference_name} \
                --primer_bed ~{primer_bed} 
    >>>

    output {
        File aln_metrics = "~{project_name}_aligned_metrics_summary.csv"
        File percent_coverage = "~{project_name}_percent_coverage.csv"
    }

    runtime {
        #cpu: 
        #memory: 
        docker: docker
    }
}

