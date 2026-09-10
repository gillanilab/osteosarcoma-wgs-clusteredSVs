# Process multiple samples with ShatterSeek and combine results

library(ShatterSeek)
library(data.table)
library(dplyr)

# Configuration
CREATE_PLOTS <- FALSE  # Set to FALSE for full batch run to save space
PLOT_FIRST_N_SAMPLES <- 10  # Only create plots for first N samples
GENOME <- "hg38"

# Load formatted data files
sv_data <- read.csv("../data/shatterseek_sv_data.csv", stringsAsFactors = FALSE)
cnv_data <- read.csv("../data/shatterseek_cnv_data.csv", stringsAsFactors = FALSE)

# Get unique samples that have both SV and CNV data
sv_samples <- unique(sv_data$sample_id)
cnv_samples <- unique(cnv_data$sample_id)
samples_to_process <- intersect(sv_samples, cnv_samples)

cat(sprintf("Found %d samples with both SV and CNV data\n", length(samples_to_process)))

# select samples to process (REMOVE this line when running full batch)
#samples_to_process <- c("SJOS013_D-SJOS013_H", "SJOS001111_M1-SJOS001111_G1", "SJOS009_D-SJOS009_G")
#cat(sprintf("Processing %d samples for testing\n", length(samples_to_process)))

# Initialize master results dataframe
master_results <- data.frame()

# Process each sample
for (i in seq_along(samples_to_process)) {
  sample_id <- samples_to_process[i]
  
  cat(sprintf("\n=== Processing sample %d/%d: %s ===\n", i, length(samples_to_process), sample_id))
  
  tryCatch({
    # Filter data for this sample
    sample_sv <- sv_data[sv_data$sample_id == sample_id, ]
    sample_cnv <- cnv_data[cnv_data$sample_id == sample_id, ]
    
    cat(sprintf("Sample %s: %d SVs, %d CNV segments\n", 
                sample_id, nrow(sample_sv), nrow(sample_cnv)))
    
    # Skip if no data
    if (nrow(sample_sv) == 0 || nrow(sample_cnv) == 0) {
      cat("Skipping sample with no SV or CNV data\n")
      next
    }
        # Check and clean chromosome names
    cat("Checking chromosome names...\n")
    
    # Valid chromosomes for ShatterSeek
    valid_chroms <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                      "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                      "20", "21", "22", "X")
    
    # Clean SV chromosome names
    sample_sv$chrom1 <- as.character(sample_sv$chrom1)
    sample_sv$chrom2 <- as.character(sample_sv$chrom2)
    
    # Remove any remaining "chr" prefix and clean up
    sample_sv$chrom1 <- gsub("^chr", "", sample_sv$chrom1, ignore.case = TRUE)
    sample_sv$chrom2 <- gsub("^chr", "", sample_sv$chrom2, ignore.case = TRUE)
    
    # Clean CNV chromosome names
    sample_cnv$chromosome <- as.character(sample_cnv$chromosome)
    sample_cnv$chromosome <- gsub("^chr", "", sample_cnv$chromosome, ignore.case = TRUE)
    
    # Check for invalid chromosomes
    invalid_sv_chroms1 <- unique(sample_sv$chrom1[!sample_sv$chrom1 %in% valid_chroms])
    invalid_sv_chroms2 <- unique(sample_sv$chrom2[!sample_sv$chrom2 %in% valid_chroms])
    invalid_cnv_chroms <- unique(sample_cnv$chromosome[!sample_cnv$chromosome %in% valid_chroms])
    
    if (length(invalid_sv_chroms1) > 0) {
      cat("Invalid SV chrom1 found:", paste(invalid_sv_chroms1, collapse = ", "), "\n")
    }
    if (length(invalid_sv_chroms2) > 0) {
      cat("Invalid SV chrom2 found:", paste(invalid_sv_chroms2, collapse = ", "), "\n")
    }
    if (length(invalid_cnv_chroms) > 0) {
      cat("Invalid CNV chromosomes found:", paste(invalid_cnv_chroms, collapse = ", "), "\n")
    }
    
    # Filter to only valid chromosomes
    sample_sv <- sample_sv[sample_sv$chrom1 %in% valid_chroms & sample_sv$chrom2 %in% valid_chroms, ]
    sample_cnv <- sample_cnv[sample_cnv$chromosome %in% valid_chroms, ]
    
    cat(sprintf("After chromosome filtering: %d SVs, %d CNV segments\n", 
                nrow(sample_sv), nrow(sample_cnv)))
    
    # Skip if no valid data remaining
    if (nrow(sample_sv) == 0 || nrow(sample_cnv) == 0) {
      cat("Skipping sample - no data with valid chromosomes\n")
      next
    }
    
    # Convert to ShatterSeek objects
    cat("Converting to ShatterSeek objects...\n")
    
    SV_object <- SVs(
      chrom1 = as.character(sample_sv$chrom1), 
      pos1 = as.numeric(sample_sv$start1),
      chrom2 = as.character(sample_sv$chrom2), 
      pos2 = as.numeric(sample_sv$end2),
      SVtype = as.character(sample_sv$svclass), 
      strand1 = as.character(sample_sv$strand1),
      strand2 = as.character(sample_sv$strand2)
    )
    
    CN_object <- CNVsegs(
      chrom = as.character(sample_cnv$chromosome),
      start = sample_cnv$start,
      end = sample_cnv$end,
      total_cn = sample_cnv$total_cn
    )
    
    # Run ShatterSeek analysis
    cat("Running ShatterSeek analysis...\n")
    start_time <- Sys.time()
    chromothripsis <- shatterseek(
      SV.sample = SV_object, 
      seg.sample = CN_object, 
      genome = GENOME
    )
    end_time <- Sys.time()
    processing_time <- round(as.numeric(end_time - start_time), 2)
    cat(sprintf("Analysis completed in %.2f seconds\n", processing_time))
    
    # Extract results and add sample information
    sample_results <- chromothripsis@chromSummary
    sample_results$sample_id <- sample_id
    sample_results$processing_time_seconds <- processing_time
    sample_results$analysis_date <- as.character(Sys.time())
    
    # Add to master results
    if (nrow(master_results) == 0) {
      master_results <- sample_results
    } else {
      master_results <- rbind(master_results, sample_results)
    }
    
    # Create plots for first few samples only
    if (CREATE_PLOTS && i <= PLOT_FIRST_N_SAMPLES) {
      cat("Creating plots for significant chromosomes...\n")

        # Create plots directory if it doesn't exist
      plots_dir <- "../results/plots"
      if (!dir.exists(plots_dir)) {
        dir.create(plots_dir, recursive = TRUE)
        cat(sprintf("Created plots directory: %s\n", plots_dir))
      }
      
      # Find chromosomes with significant chromothripsis
      sig_chroms <- sample_results[
        !is.na(sample_results$chr_breakpoint_enrichment) & 
        sample_results$chr_breakpoint_enrichment < 0.05, "chrom"]
      
      if (length(sig_chroms) > 0) {
        # Limit to top 3 most significant chromosomes to save space
        if (length(sig_chroms) > 3) {
          sig_chroms <- sig_chroms[1:3]
        }
        
        for (chr in sig_chroms) {
          plot_file <- file.path(plots_dir, sprintf("%s_chr%s.png", sample_id, chr))
          
          tryCatch({
            plots <- plot_chromothripsis(
              ShatterSeek_output = chromothripsis,
              chr = paste0("chr", chr),
              sample_name = sample_id,
              genome = GENOME
            )
            
            png(plot_file, width = 14, height = 10, units = "in", res = 300)
            par(mar = c(5, 6, 4, 2))  # bottom, left, top, right margins
            
            if (is.list(plots) && length(plots) > 1) {
              if (require(gridExtra, quietly = TRUE)) {
                grid.arrange(grobs = plots, ncol = 1)
              } else {
                print(plots[[1]])
              }
            } else if (is.list(plots)) {
              print(plots[[1]])
            } else {
              print(plots)
            }
            
            dev.off()
            
            cat(sprintf("Plot saved: %s\n", plot_file))
            
          }, error = function(e) {
            cat(sprintf("Warning: Could not create plot for chr%s: %s\n", chr, e$message))
          })
        }
      } else {
        cat("No chromothripsis events found for plotting\n")
      }
    }
    
    # Save intermediate results (in case of crash)
    write.csv(master_results, "../results/shatterseek_results_intermediate.csv", row.names = FALSE)
    
    cat(sprintf("Sample %s completed successfully\n", sample_id))
    
  }, error = function(e) {
    cat(sprintf("ERROR processing sample %s: %s\n", sample_id, e$message))
    
    # Add error record to results
    error_row <- data.frame(
      chrom = NA,
      start = NA,
      end = NA,
      number_DEL = NA,
      number_DUP = NA,
      number_h2hINV = NA,
      number_t2tINV = NA,
      number_TRA = NA,
      clusterSize_including_TRA = NA,
      number_SVs_sample = nrow(sample_sv),
      number_CNV_segments = NA,
      pval_fragment_joins = NA,
      chr_breakpoint_enrichment = NA,
      pval_exp_chr = NA,
      pval_exp_cluster = NA,
      sample_id = sample_id,
      processing_time_seconds = NA,
      analysis_date = as.character(Sys.time()),
      error_message = as.character(e$message),
      stringsAsFactors = FALSE
    )
    
    # Add any missing columns to match existing results structure
    missing_cols <- setdiff(colnames(master_results), colnames(error_row))
    for (col in missing_cols) {
      error_row[[col]] <- NA
    }
    
    if (nrow(master_results) == 0) {
      master_results <- error_row
    } else {
      master_results <- rbind(master_results, error_row)
    }
  })
  
  # Memory cleanup
  gc()
}

# Final results summary
cat("\n=== FINAL SUMMARY ===\n")
cat(sprintf("Total samples processed: %d\n", length(samples_to_process)))
cat(sprintf("Total results rows: %d\n", nrow(master_results)))

# Count successful vs failed samples
if (nrow(master_results) > 0) {
  # Add error_message column if it doesn't exist
  if (!"error_message" %in% colnames(master_results)) {
    master_results$error_message <- NA
  }
  
  successful_samples <- length(unique(master_results$sample_id[is.na(master_results$error_message)]))
  failed_samples <- length(unique(master_results$sample_id[!is.na(master_results$error_message)]))
} else {
  successful_samples <- 0
  failed_samples <- 0
}

cat(sprintf("Successful: %d, Failed: %d\n", successful_samples, failed_samples))

# Count chromothripsis events
sig_events <- sum(!is.na(master_results$pval_exp_chr) & master_results$pval_exp_chr < 0.05)
cat(sprintf("Significant chromothripsis events detected: %d\n", sig_events))

# Save final results
final_output_file <- "../results/shatterseek_master_results.csv"
write.csv(master_results, final_output_file, row.names = FALSE)
cat(sprintf("Final results saved to: %s\n", final_output_file))

# Create summary table by sample
if (nrow(master_results) > 0) {
  sample_summary <- master_results %>%
    group_by(sample_id) %>%
    summarise(
      total_SVs = first(number_SVs_sample),
      processing_time = first(processing_time_seconds),
      has_error = any(!is.na(error_message)),
      .groups = 'drop'
    )
  write.csv(sample_summary, "../results/shatterseek_sample_summary.csv", row.names = FALSE)
  cat("Sample summary saved to: ../results/shatterseek_sample_summary.csv\n")
} else {
  cat("No results to summarize\n")
}

cat("\nAnalysis complete!\n")