suppressPackageStartupMessages({
  library(fishHook)
  library(GenomicRanges)
  library(data.table)
  library(dplyr)
})

cat("Starting fishHook analysis...\n")

# Load data
covs <- read.table('../data/covariates_for_fishhook.txt', header = TRUE, stringsAsFactors = FALSE)
breaks <- fread('../data/breakpoints_for_fishhook.bed', header = FALSE)
names(breaks)[1:3] <- c("chrom", "start", "end")

cat("Loaded", nrow(covs), "covariate regions and", nrow(breaks), "breakpoints\n")

# Standardize chromosome names
if(!any(grepl("chr", covs$chrom))) covs$chrom <- paste0("chr", covs$chrom)
if(!any(grepl("chr", breaks$chrom))) breaks$chrom <- paste0("chr", breaks$chrom)

# Get common chromosomes
common_chrs <- intersect(unique(covs$chrom), unique(breaks$chrom))
cat("Analyzing", length(common_chrs), "chromosomes\n")

# Create GRanges objects
# Covariates
create_cov_gr <- function(data, score_col) {
  gr <- GRanges(seqnames = data$chrom, ranges = IRanges(start = data$start + 1, end = data$end))
  mcols(gr)$score <- as.numeric(data[[score_col]])
  return(gr)
}

# Create all covariate GRanges
gc_gr <- create_cov_gr(covs, "gc_content")
rt_gr <- create_cov_gr(covs, "replication_timing") 
map_gr <- create_cov_gr(covs, "mappability")
repressed_gr <- create_cov_gr(covs, "frac_heterochromatin")
sine_gr <- create_cov_gr(covs, "frac_sine")
fragile_site_gr <- create_cov_gr(covs, "frac_fragile_site")

# Combine all the covariates
covariates <- c(
  Cov(gc_gr, field = "score", name = 'GC'),
  Cov(rt_gr, field = "score", name = 'ReplicationTiming'),
  Cov(map_gr, field = "score", name = 'Mappability'),
  Cov(repressed_gr, field = "score", name = 'Heterochromatin'),
  Cov(sine_gr, field = "score", name = 'SINE'),
  Cov(fragile_site_gr, field = "score", name = 'FragileSite')
)

# Convert original breakpoints to point events  
original_events_gr <- GRanges(
  seqnames = breaks$chrom,
  ranges = IRanges(start = breaks$end, end = breaks$end)  # Point events
)

# Use the bins defined in the covariates file (where there is complete covariate data) as hypotheses (regions to test)
hypotheses_gr <- GRanges(
  seqnames = covs$chrom,
  ranges = IRanges(start = covs$start + 1, end = covs$end)
)

# eligible regions construction (not needed, but here in case)
eligible_gr <- GRanges()
for(chr in common_chrs) {
  chr_covs <- covs[covs$chrom == chr, ]
  if(nrow(chr_covs) > 0) {
    # Create individual ranges for each covariate region (exludes blacklisted & gap regions)
    chr_ranges <- GRanges(
      seqnames = chr_covs$chrom,
      ranges = IRanges(start = chr_covs$start + 1, end = chr_covs$end)
    )
    # Reduce to merge overlapping regions
    chr_ranges <- reduce(chr_ranges)
    eligible_gr <- c(eligible_gr, chr_ranges)
  }
}

cat("Created", length(original_events_gr), "point events,", length(hypotheses_gr), "hypothesis bins,", length(eligible_gr), "eligible regions\n")

# checks
cat("Checking overlaps between events and eligible regions...\n")
events_in_eligible <- countOverlaps(original_events_gr, eligible_gr) > 0
cat("Events within eligible regions:", sum(events_in_eligible), "out of", length(original_events_gr), "\n")

cat("Checking overlaps between hypotheses and eligible regions...\n")
hyp_in_eligible <- countOverlaps(hypotheses_gr, eligible_gr) > 0
cat("Hypotheses within eligible regions:", sum(hyp_in_eligible), "out of", length(hypotheses_gr), "\n")

# Run fishHook
cat("Running fishHook analysis...\n")

# Create Fish object with point events and bin hypotheses
fh <- Fish(
  events = original_events_gr,     # Individual breakpoints
  hypotheses = hypotheses_gr,      # Bins to test
  eligible = hypotheses_gr,        # may have to change this to hypotheses_gr?? 
  covariates = covariates
)

# Score with negative binomial
fh$score(nb = TRUE, verbose = FALSE)

# Check model type
model_type <- "Unknown"
if(!is.null(fh$model) && "glm" %in% class(fh$model) && !is.null(fh$model$family)) {
  if(grepl("Negative Binomial", fh$model$family$family)) {
    model_type <- "Negative Binomial"
  } else if(grepl("poisson", fh$model$family$family, ignore.case = TRUE)) {
    model_type <- "Poisson"
  }
}
cat("Model used:", model_type, "\n")

# Extract and display model coefficients
if(!is.null(fh$model)) {
  cat("\n=== MODEL COEFFICIENTS ===\n")
  cat("Model convergence:", fh$model$converged, "\n")
  if(!is.null(fh$model$boundary)) {
    cat("Hit boundary:", fh$model$boundary, "\n")
  }
  # Get coefficient summary
  model_summary <- summary(fh$model)
  coef_table <- model_summary$coefficients
  
  cat("Coefficient table:\n")
  print(coef_table)
  
  # Create coefficient table for saving
  coef_df <- data.frame(
    covariate = rownames(coef_table),
    estimate = coef_table[, "Estimate"],
    std_error = coef_table[, "Std. Error"],
    z_value = coef_table[, "z value"],
    p_value = coef_table[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  
  # Add interpretations
  coef_df$direction <- ifelse(coef_df$estimate > 0, "increases", "decreases")
  coef_df$percent_change <- round((exp(coef_df$estimate) - 1) * 100, 1)
  coef_df$significance <- ifelse(coef_df$p_value < 0.001, "***", 
                                ifelse(coef_df$p_value < 0.01, "**", 
                                      ifelse(coef_df$p_value < 0.05, "*", "ns")))
  
  write.table(coef_df, '../results/fishhook_coefficients.txt', sep = '\t', row.names = FALSE, quote = FALSE)


  # interpretations for coefficients
  cat("\nCoefficient interpretation:\n")
  coef_names <- rownames(coef_table)
  coef_values <- coef_table[, "Estimate"]
  coef_pvals <- coef_table[, "Pr(>|z|)"]
  
  for(i in 1:length(coef_names)) {
    name <- coef_names[i]
    value <- coef_values[i]
    pval <- coef_pvals[i]
    
    if(name == "(Intercept)") {
      cat("- Intercept:", round(value, 4), "(p =", format(pval, scientific = TRUE, digits = 3), ")\n")
      cat("  → Baseline log breakpoint rate\n")
    } else {
      direction <- ifelse(value > 0, "increases", "decreases")
      significance <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "ns")))
      
      cat("-", name, ":", round(value, 4), "(p =", format(pval, scientific = TRUE, digits = 3), ")", significance, "\n")
      cat("  → Higher", gsub("track", "", name), direction, "breakpoint rate by", round((exp(value) - 1) * 100, 1), "%\n")
    }
  }
  
  # Model fit stats 
  cat("\nModel fit statistics:\n")
  cat("- AIC:", round(fh$model$aic, 2), "\n")
  cat("- Deviance:", round(fh$model$deviance, 2), "\n")
  cat("- Null deviance:", round(fh$model$null.deviance, 2), "\n")
  cat("- Pseudo R-squared:", round(1 - fh$model$deviance/fh$model$null.deviance, 3), "\n")
  
  if(model_type == "Negative Binomial" && !is.null(fh$model$theta)) {
    cat("- Theta (dispersion):", round(fh$model$theta, 4), "\n")
    cat("- Log(theta) SE:", round(fh$model$SE.theta, 4), "\n")
  }
}

# Extract results
fish_results <- fh$res

if(class(fish_results)[1] == "GRanges" && length(fish_results) > 0) {
  # Debug: Check what columns are actually in the results
  cat("\nDEBUG: Available columns in fishHook results:\n")
  print(colnames(mcols(fish_results)))
  
  # Convert GRanges to data frame
  mcols_data <- mcols(fish_results)
  result_data <- data.frame(
    chr = as.character(seqnames(fish_results)),
    start = start(fish_results),
    end = end(fish_results),
    as.data.frame(mcols_data),
    stringsAsFactors = FALSE
  )
  
  # Debug: Check initial column names
  cat("Initial column names:", paste(names(result_data), collapse = ", "), "\n")
  

  if("p" %in% names(result_data)) {
    names(result_data)[names(result_data) == "p"] <- "pvalue"
  }
  if("fdr" %in% names(result_data)) {
    names(result_data)[names(result_data) == "fdr"] <- "qvalue"
  }
  if("count.pred" %in% names(result_data)) {
    names(result_data)[names(result_data) == "count.pred"] <- "expected"
  }
  
  if(!"count" %in% names(result_data)) {
    # If count is missing, we need to recalculate it
    cat("WARNING: 'count' column missing from results. Recalculating...\n")
    
    # count events in each hypothesis region
    overlaps <- findOverlaps(original_events_gr, fish_results)
    count_table <- table(subjectHits(overlaps))
    
    result_data$count <- 0
    result_data$count[as.numeric(names(count_table))] <- as.numeric(count_table)
  }
  
  # Debug: Check final column names and data types
  cat("Final column names:", paste(names(result_data), collapse = ", "), "\n")
  cat("Count column summary:\n")
  print(summary(result_data$count))
  
  # fold change
  if("expected" %in% names(result_data)) {
    result_data$fold_change <- result_data$count / pmax(result_data$expected, 1e-10)
  }
  
  # Sort by p-value
  if("pvalue" %in% names(result_data)) {
    result_data <- result_data[order(result_data$pvalue), ]
  }
  
  # Summary and diagnostics
  n_sig_p <- if("pvalue" %in% names(result_data)) sum(result_data$pvalue < 0.05, na.rm = TRUE) else 0
  n_sig_q <- if("qvalue" %in% names(result_data)) sum(result_data$qvalue < 0.05, na.rm = TRUE) else 0
  
  cat("Results: ", nrow(result_data), " regions, ", n_sig_p, " significant (p<0.05), ", n_sig_q, " significant (q<0.05)\n", sep="")
  
  # Check expected counts
  if("expected" %in% names(result_data)) {
    cat("\nDiagnostic - Expected count summary:\n")
    print(summary(result_data$expected))
    cat("Total expected:", round(sum(result_data$expected, na.rm = TRUE)), "\n")
  }
  
  # Check observed counts
  if("count" %in% names(result_data) && !all(is.na(result_data$count))) {
    cat("Total breakpoints:", sum(result_data$count, na.rm = TRUE), "\n")
    
    if("expected" %in% names(result_data)) {
      total_obs <- sum(result_data$count, na.rm = TRUE)
      total_exp <- sum(result_data$expected, na.rm = TRUE)
      if(total_exp > 0) {
        cat("Overall observed/expected ratio:", round(total_obs / total_exp, 2), "\n")
      }
      
      mean_expected <- mean(result_data$expected, na.rm = TRUE)
      expected_per_bin <- total_obs / nrow(result_data)
      cat("Mean expected per bin (fishHook):", round(mean_expected, 3), "\n")
      cat("Mean observed per bin (uniform):", round(expected_per_bin, 3), "\n")
      
      # NA checks
      if(!is.na(mean_expected) && !is.na(expected_per_bin) && mean_expected < expected_per_bin * 0.5) {
        cat("WARNING: fishHook expected counts seem too low!\n")
        cat("This suggests over-correction by covariates or model issues.\n")
      }
    }
  } else {
    cat("WARNING: Count data is missing or all NA\n")
  }
  
  # Save results
  write.table(result_data, '../results/fishhook_results.txt', sep = '\t', row.names = FALSE, quote = FALSE)

  # Save significant results
  if(n_sig_p > 0 && "pvalue" %in% names(result_data)) {
    sig_results <- result_data[result_data$pvalue < 0.05 & !is.na(result_data$pvalue), ]
    write.table(sig_results, '../results/fishhook_significant_results.txt', sep = '\t', row.names = FALSE, quote = FALSE)
    
    cat("\nTop significant regions:\n")
    display_cols <- intersect(c("chr", "start", "end", "count", "pvalue", "qvalue", "fold_change"), names(sig_results))
    print(head(sig_results[, display_cols], 10))
  } else {
    cat("No significant regions found\n")
  }
  
  cat("\nAnalysis complete!\n")
  
} else {
  cat("No results found or results object is not a GRanges\n")
  cat("Results class:", class(fh$res), "\n")
  if(!is.null(fh$res)) {
    cat("Results length:", length(fh$res), "\n")
  }
}