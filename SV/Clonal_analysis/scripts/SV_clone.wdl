version 1.0

workflow SVclone_Clustering {
    input {
        String sample_id
        File tumor_bam
        File tumor_bai
        File sv_vcf
        File cnv_filex
        File snv_file
        File config_file
        Float purity
        Float ploidy
    }

    call run_svclone {
        input:
            sample_id = sample_id,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            sv_vcf = sv_vcf,
            cnv_file = cnv_file,
            snv_file = snv_file,
            config_file = config_file,
            purity = purity,
            ploidy = ploidy
    }

    output {
        File cluster_results = run_svclone.results
    }
}

task run_svclone {
    input {
        String sample_id
        File tumor_bam
        File tumor_bai
        File sv_vcf
        File cnv_file
        File snv_file
        File config_file
        Float purity
        Float ploidy
        Int memory_gb = 64
        Int disk_gb = 300
    }

    command <<<
        set -ex
        
        echo "1. Patching config file..."
        cp ~{config_file} custom_config.ini
        sed -i 's/insert_mean.*/insert_mean: 400/g' custom_config.ini
        sed -i 's/insert_std.*/insert_std: 50/g' custom_config.ini
        
        echo "1.5 Patching SVclone to allow mixed read lengths..."
        # Using the syntax-specific match to change the default argument to "mean"
        BAMTOOLS_PY="/usr/local/lib/python3.11/site-packages/SVclone/SVprocess/bamtools.py"
        if [ -f "$BAMTOOLS_PY" ]; then
            sed -i 's/multiple="error"):/multiple="mean"):/g' "$BAMTOOLS_PY"
            echo "Successfully patched bamtools.py!"
        else
            echo "Warning: bamtools.py not found at expected path!"
        fi

        echo "2. Creating purity ploidy file..."
        echo -e "sample\tpurity\tploidy" > purity_ploidy.txt
        echo -e "~{sample_id}\t~{purity}\t~{ploidy}" >> purity_ploidy.txt
        
echo "3. Reformatting SVs (Ensuring 'chr' prefix and syncing X/Y)..."
        awk -F'\t' -v OFS='\t' 'NR==1 {print; next} {
            # Add "chr" if it is missing
            if($1 !~ /^chr/) $1="chr"$1;
            if($4 !~ /^chr/) $4="chr"$4;
            # Handle X/Y (if your BAM uses chrX/chrY, keep them as chrX/chrY)
            gsub(/chr23/, "chrX", $1); gsub(/chr24/, "chrY", $1);
            gsub(/chr23/, "chrX", $4); gsub(/chr24/, "chrY", $4);
            print
        }' ~{sv_vcf} > formatted_svs.txt
        
        echo "4. Reformatting CNVs (Ensuring 'chr' prefix)..."
        # Adjusted to output standard tab-separated format with chr prefix
        awk -F',' -v OFS='\t' 'BEGIN {print "chr", "startpos", "endpos", "norm_total", "norm_minor", "tumour_total", "tumour_minor"} NR>1 {
            chrom=$2;
            if(chrom !~ /^chr/) chrom="chr"chrom;
            gsub(/chr23/, "chrX", chrom); gsub(/chr24/, "chrY", chrom);
            print chrom, $3, $4, $5, $6, $7, $8
        }' ~{cnv_file} > formatted_cnvs.txt

        echo "4.5 Reformatting SNV VCF (Ensuring 'chr' prefix)..."
        sed -E '/^[#][#]/ s/Description=([^"][^>]*)>/Description="\1">/' ~{snv_file} > patched_header_snvs.vcf
        
        awk -F'\t' -v OFS='\t' '{
            if ($0 ~ /^[#]/) {
                # Ensure contig headers match BAM (chr1, etc)
                print $0;
            } else {
                if($1 !~ /^chr/) $1="chr"$1;
                gsub(/chr23/, "chrX", $1); gsub(/chr24/, "chrY", $1);
                print $0;
            }
        }' patched_header_snvs.vcf > formatted_snvs.vcf
        
        echo "5. Running svclone Annotate..."
        svclone annotate -i formatted_svs.txt -b ~{tumor_bam} -s ~{sample_id} --sv_format simple --config custom_config.ini
        
        echo "6. Running svclone Count..."
        svclone count -i ~{sample_id}/~{sample_id}_svin.txt -b ~{tumor_bam} -s ~{sample_id} --config custom_config.ini
        
        echo "7. Running svclone Filter..."
        # We use the new 'formatted_cnvs.txt' here
        svclone filter -i ~{sample_id}/~{sample_id}_svinfo.txt \
         -s ~{sample_id} \
         -c formatted_cnvs.txt \
         -p purity_ploidy.txt \
         --snvs  formatted_snvs.vcf \
         --snv_format mutect \
         --config custom_config.ini
        
        echo "8. Running svclone Cluster..."
        svclone cluster -s ~{sample_id} --config custom_config.ini
        
        echo "9. Zipping results..."
        tar -czvf ~{sample_id}_svclone_results.tar.gz ~{sample_id}/
    >>>

    runtime {
        docker: "quay.io/biocontainers/svclone:1.1.4--pyr44hdfd78af_0"
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} SSD"
        cpu: 4
    }

    output {
        File results = "~{sample_id}_svclone_results.tar.gz"
    }
}