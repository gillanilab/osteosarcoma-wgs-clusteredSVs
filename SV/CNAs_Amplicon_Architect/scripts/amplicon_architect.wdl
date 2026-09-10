workflow AmpliconSuitePipeline {
    File        bam_file
    File        bam_index
    File        cnv_bed
    File        aa_data_repo_tar
    File        mosek_license
    String      sample_name
    String      reference_genome
    Int         threads
    Int         memory_gb
    Int         disk_gb
    Float       cn_gain

    call RunAmpliconSuite {
        input:
            bam_file           = bam_file,
            bam_index          = bam_index,
            cnv_bed            = cnv_bed,
            aa_data_repo_tar   = aa_data_repo_tar,
            mosek_license      = mosek_license,
            sample_name        = sample_name,
            reference_genome   = reference_genome,
            threads            = threads,
            memory_gb          = memory_gb,
            disk_gb            = disk_gb,
            cn_gain            = cn_gain
    }

    output {
            File  aa_results_tar               = RunAmpliconSuite.aa_results_tar
            File  stdout_log                   = RunAmpliconSuite.stdout_log
            File?  summary_txt                  = RunAmpliconSuite.summary_txt
            File?  result_table                 = RunAmpliconSuite.result_table
            File?  amplicon_classification      = RunAmpliconSuite.amplicon_classification
            File?  gene_list                    = RunAmpliconSuite.gene_list
        }
}


task RunAmpliconSuite {
    File        bam_file
    File        bam_index
    File        cnv_bed
    File        aa_data_repo_tar
    File        mosek_license
    String      sample_name
    String      reference_genome
    Int         threads   = 8
    Int         memory_gb = 32
    Int         disk_gb   = 250
    Float       cn_gain   = 4.5

    command {
        set -euo pipefail

        # Mosek license
        mkdir -p mosek
        cp ${mosek_license} mosek/mosek.lic
        export MOSEKLM_LICENSE_FILE=mosek/

        # Localize to AA Data Repo
        mkdir -p data_repo
        tar -xzf ${aa_data_repo_tar} -C data_repo/
        export AA_DATA_REPO=data_repo

        mkdir -p input
        cp ${bam_file}  input/
        cp ${bam_index} input/
        cp ${cnv_bed}   input/

        BAM_FILENAME=$(basename ${bam_file})
        BED_FILENAME=$(basename ${cnv_bed})
        # 09 add chr to bed filename 
        sed '/^chr/! s/^/chr/' ${cnv_bed} > input/$BED_FILENAME
        echo "Input files staged."

        export AA_SRC=/home/programs/AmpliconArchitect-master/src
        export AC_SRC=/home/programs/AmpliconClassifier-main
        export NCM_HOME=/home/programs/NGSCheckMate-master/

        # Run AmpliconSuite-pipeline
        echo "Launching AmpliconSuite-pipeline."

        python3 /home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py \
            -s ${sample_name} \
            -t ${threads} \
            --bam input/$BAM_FILENAME \
            --cnv_bed input/$BED_FILENAME \
            --ref ${reference_genome} \
            --run_AA \
            --run_AC \
            --cngain ${cn_gain} \
            -o . \
            2>&1 | tee AS-p_stdout.log

        echo "Pipeline complete."

        tar --exclude="${sample_name}_outputs.tar.gz" \
            --exclude="*.tar" \
            --exclude="*.tar.gz" \
            --exclude="./data_repo" \
            --exclude="./programs" \
            --exclude="./testdata" \
            --exclude="./input" \
            --exclude="*.bam" \
            --exclude="*.bai" \
            --exclude="*.fastq*" \
            --exclude="*.fq*" \
            -zcf ${sample_name}_outputs.tar.gz .

        echo "Done."
    }

    output {
        File  aa_results_tar               = "${sample_name}_outputs.tar.gz"
        File  stdout_log                   = "AS-p_stdout.log"
        File?  summary_txt                  = "${sample_name}_AA_results/${sample_name}_summary.txt"
        File?  result_table                 = "${sample_name}_classification/${sample_name}_result_table.tsv"
        File?  amplicon_classification      = "${sample_name}_classification/${sample_name}_amplicon_classification_profiles.tsv"
        File?  gene_list                    = "${sample_name}_classification/${sample_name}_gene_list.tsv"
    }

    runtime {
        docker:      "jluebeck/ampliconsuite-pipeline"
        memory:      "${memory_gb} GB"
        cpu:         threads
        disks:       "local-disk ${disk_gb} SSD"
        preemptible: 1
        maxRetries:  1
    }
}