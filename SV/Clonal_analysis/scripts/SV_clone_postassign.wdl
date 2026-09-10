version 1.0

workflow SVclone_PostAssign {
    input {
        String sample_id
        File svclone_results_tar
    }

    call run_postassign {
        input:
            sample_id = sample_id,
            svclone_results_tar = svclone_results_tar,
    }

    output {
        File postassign_cluster_results = run_postassign.results
    }
}

task run_postassign {
    input {
        String sample_id
        File svclone_results_tar
        
        # Lower memory and disk requirements since we aren't processing BAMs
        Int memory_gb = 64
        Int disk_gb = 100 
    }

    command <<<
        set -ex

        echo "1. Extracting previous SVclone results..."
        # This will extract the ~{sample_id}/ folder into the current working directory
        tar -xzvf ~{svclone_results_tar}

        echo "2. Running svclone PostAssign..."
        svclone postassign -s ~{sample_id} --joint

        echo "3. Re-zipping final results..."
        # Create a new tarball with a slightly different name to distinguish it
        tar -czvf ~{sample_id}_svclone_postassign_joint.tar.gz ~{sample_id}/
    >>>

    runtime {
        docker: "quay.io/biocontainers/svclone:1.1.4--pyr44hdfd78af_0"
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} SSD"
        cpu: 4
    }

    output {
        File results = "~{sample_id}_svclone_postassign_joint.tar.gz"
    }
}