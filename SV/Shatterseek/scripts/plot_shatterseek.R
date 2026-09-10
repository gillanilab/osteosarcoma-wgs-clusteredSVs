# Plot ShatterSeek results for CGR events from existing results file
# This script reads the filtered results and creates plots based on a target region without re-running ShatterSeek analysis

library(ShatterSeek)
library(data.table)
library(dplyr)

# Configuration
GENOME <- "hg38"
FILTER_CONDITION <- "is_cgr"  # Column to filter on
FILTER_VALUE <- "True"  # Value to filter for (use "True" if it's a string in CSV)
MAX_PLOTS_PER_SAMPLE <- 2 # Limit plots per sample to avoid too many files

TARGET_REGION <- list( chrom = "19", start = 5000000, end = 140000000)
TARGET_SAMPLE <- "SJOS012_D-SJOS012_G"
# TARGET_SAMPLE <- NULL
# TARGET_REGION <- NULL

# Function to check if an event overlaps with target region
overlaps_region <- function(event_chrom, event_start, event_end, target_region) {
  # Catch any errors and return FALSE
  tryCatch({
    if (is.null(target_region)) {
      return(TRUE)  # No filtering if target_region is NULL
    }
    
    # Check for NULL or length 0
    if (is.null(event_chrom) || length(event_chrom) == 0 ||
        is.null(event_start) || length(event_start) == 0 ||
        is.null(event_end) || length(event_end) == 0) {
      return(FALSE)
    }
    
    # Check for NA values - return FALSE if any are NA
    if (is.na(event_chrom) || is.na(event_start) || is.na(event_end)) {
      return(FALSE)
    }
    
    # Ensure chromosome matches (handle both with and without 'chr' prefix)
    event_chrom_clean <- gsub("^chr", "", as.character(event_chrom), ignore.case = TRUE)
    target_chrom_clean <- gsub("^chr", "", as.character(target_region$chrom), ignore.case = TRUE)
    
    if (event_chrom_clean != target_chrom_clean) {
      return(FALSE)
    }
    
    # Check for overlap: event overlaps if it's not completely before or after the region
    # Overlap exists if: event_end >= region_start AND event_start <= region_end
    overlap <- (as.numeric(event_end) >= as.numeric(target_region$start)) & 
               (as.numeric(event_start) <= as.numeric(target_region$end))
    
    if (is.na(overlap)) {
      return(FALSE)
    }
    return(overlap)
    
  }, error = function(e) {
    # If any error occurs, return FALSE
    return(FALSE)
  })
}

# Input files
results_file <- "../results/shatterseek_results_filtering_criteria.csv"
sv_data_file <- "../data/shatterseek_sv_data.csv"
cnv_data_file <- "../data/shatterseek_cnv_data.csv"
plots_dir <- "../results/plots"

# Create plots directory
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
  cat(sprintf("Created plots directory: %s\n", plots_dir))
}

# Read the results file
cat("Reading results file...\n")
results <- read.csv(results_file, stringsAsFactors = FALSE)
cat(sprintf("Total results rows: %d\n", nrow(results)))

# Filter for desired events (is_cgr == TRUE)
if (is.logical(FILTER_VALUE)) {
  filtered_results <- results[results[[FILTER_CONDITION]] == FILTER_VALUE | 
                               results[[FILTER_CONDITION]] == as.character(FILTER_VALUE), ]
} else {
  filtered_results <- results[results[[FILTER_CONDITION]] == FILTER_VALUE, ]
}

# Remove rows with NA chromosomes
filtered_results <- filtered_results[!is.na(filtered_results$chrom), ]

# Apply region filter if specified
if (!is.null(TARGET_REGION)) {
  cat(sprintf("Filtering for region: chr%s:%d-%d\n", 
              TARGET_REGION$chrom, TARGET_REGION$start, TARGET_REGION$end))
  
  # Apply overlap check to each row
  region_overlaps <- sapply(1:nrow(filtered_results), function(i) {
    overlaps_region(
      filtered_results$chrom[i],
      filtered_results$start[i],
      filtered_results$end[i],
      TARGET_REGION
    )
  })
  
  filtered_results <- filtered_results[region_overlaps, ]
  
  cat(sprintf("Events after region filtering: %d\n", nrow(filtered_results)))
}

# Apply sample filter if specified
if (!is.null(TARGET_SAMPLE)) {
  cat(sprintf("Filtering for sample: %s\n", TARGET_SAMPLE))
  filtered_results <- filtered_results[filtered_results$sample_id == TARGET_SAMPLE, ]
  cat(sprintf("Events after sample filtering: %d\n", nrow(filtered_results)))
}

cat(sprintf("Filtered events to plot: %d\n", nrow(filtered_results)))

if (nrow(filtered_results) == 0) {
  cat("No events found matching filter criteria. Exiting.\n")
  quit()
}

# Load original SV and CNV data
cat("Loading original SV and CNV data...\n")
sv_data <- read.csv(sv_data_file, stringsAsFactors = FALSE)
cnv_data <- read.csv(cnv_data_file, stringsAsFactors = FALSE)

# Get unique sample-chromosome combinations to plot
events_to_plot <- filtered_results %>%
  select(sample_id, chrom) %>%
  distinct() %>%
  arrange(sample_id, chrom)

cat(sprintf("Unique sample-chromosome combinations: %d\n", nrow(events_to_plot)))

# Valid chromosomes for ShatterSeek
valid_chroms <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                  "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                  "20", "21", "22", "X")

# Process each event
for (i in seq_len(nrow(events_to_plot))) {
  sample_id <- events_to_plot$sample_id[i]
  target_chr <- as.character(events_to_plot$chrom[i])
  
  cat(sprintf("\n=== Processing %d/%d: %s, chr%s ===\n", 
              i, nrow(events_to_plot), sample_id, target_chr))
  
  tryCatch({
    # Filter data for this sample
    sample_sv <- sv_data[sv_data$sample_id == sample_id, ]
    sample_cnv <- cnv_data[cnv_data$sample_id == sample_id, ]
    
    cat(sprintf("Sample %s: %d SVs, %d CNV segments\n", 
                sample_id, nrow(sample_sv), nrow(sample_cnv)))
    
    if (nrow(sample_sv) == 0 || nrow(sample_cnv) == 0) {
      cat("Skipping - no SV or CNV data\n")
      next
    }
    
    # Clean chromosome names
    sample_sv$chrom1 <- gsub("^chr", "", as.character(sample_sv$chrom1), ignore.case = TRUE)
    sample_sv$chrom2 <- gsub("^chr", "", as.character(sample_sv$chrom2), ignore.case = TRUE)
    sample_cnv$chromosome <- gsub("^chr", "", as.character(sample_cnv$chromosome), ignore.case = TRUE)
    
    # Filter to only valid chromosomes
    sample_sv <- sample_sv[sample_sv$chrom1 %in% valid_chroms & 
                           sample_sv$chrom2 %in% valid_chroms, ]
    sample_cnv <- sample_cnv[sample_cnv$chromosome %in% valid_chroms, ]
    
    if (nrow(sample_sv) == 0 || nrow(sample_cnv) == 0) {
      cat("Skipping - no data with valid chromosomes\n")
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
    
    # Run ShatterSeek analysis (needed to generate chromothripsis object for plotting)
    cat("Running ShatterSeek analysis...\n")
    chromothripsis <- shatterseek(
      SV.sample = SV_object, 
      seg.sample = CN_object, 
      genome = GENOME
    )
    
    # Create plot for the target chromosome
    cat(sprintf("Creating plot for chr%s...\n", target_chr))
    
    plot_file <- file.path(plots_dir, sprintf("%s_chr%s_CGR.png", sample_id, target_chr))

    png(plot_file, width = 14, height = 10, units = "in", res = 300)
    par(mar = c(5, 6, 4, 2))
    
    plots <- plot_chromothripsis(
      ShatterSeek_output = chromothripsis,
      chr = paste0("chr", target_chr),
      sample_name = sample_id,
      genome = GENOME
    )

    # Print the plots
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
    
    # Get event details from filtered results
    event_details <- filtered_results[filtered_results$sample_id == sample_id & 
                                      filtered_results$chrom == target_chr, ]
    
    if (nrow(event_details) > 0) {
      cat(sprintf("  pval_exp_cluster: %s\n", event_details$pval_exp_cluster[1]))
      cat(sprintf("  intrachrom_svs: %s\n", event_details$intrachrom_svs[1]))
      cat(sprintf("  is_chromothripsis: %s\n", event_details$is_chromothripsis[1]))
    }
    
  }, error = function(e) {
    cat(sprintf("ERROR creating plot for %s chr%s: %s\n", sample_id, target_chr, e$message))
  })
  
  # Memory cleanup
  gc()
}

# Create summary
cat("\n=== PLOTTING SUMMARY ===\n")
cat(sprintf("Total events plotted: %d\n", nrow(events_to_plot)))
cat(sprintf("Plots saved to: %s\n", plots_dir))

# Create a summary table of what was plotted
plotted_summary <- filtered_results %>%
  select(sample_id, chrom, pval_exp_cluster, intrachrom_svs, is_chromothripsis) %>%
  arrange(sample_id, chrom)

write.csv(plotted_summary, file.path(plots_dir, "plotted_events_summary.csv"), row.names = FALSE)
cat("Summary of plotted events saved to: plotted_events_summary.csv\n")

cat("\nPlotting complete!\n")