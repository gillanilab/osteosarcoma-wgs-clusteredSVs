################################################################################
## Figure 3D
## Pairwise cSV (CGR) co-occurrence / mutual exclusivity
## Firth penalized logistic regression, purity-adjusted
################################################################################


library(data.table)   
library(dplyr)        
library(ggplot2)     
library(logistf)    

metadata_file   <- "../../../metadata/metadata_v12.tsv"
cooc_plot_out    <- "../../../Images/Figures/Co-occurence_firth_purity_corrected.pdf"
cooc_table_out   <- "../../../Images/Figures/Tables/Co-occurence_firth_purity_corrected_cSV_vs_cSV.txt"
mutex_plot_out   <- "../../../Images/Figures/Mutually_exclusivity_firth_purity_corrected.pdf"
mutex_table_out  <- "../../../Images/Tables/Mutual_exclusivity_firth_purity_corrected_cSV_vs_cSV.txt"


# Load metadata

metadata     <- fread(metadata_file)
meta_heatmap <- copy(metadata)

# Remove multiple-tumor samples (matches the main heatmap script)
if ("filter status" %in% names(meta_heatmap)) {
  meta_heatmap <- meta_heatmap[`filter status` != "multiple-tumor"]
}

CGR_cols <- c("MYC_CGR", "CCND3_CGR", "CCNE1_CGR", "TP53_CGR", "CDK4-MDM2_CGR")

# Derive cSV Hotspots column
# A sample is a "hotspot carrier" if it has a clustered SV (CGR) event
# at ANY of the tracked hotspot loci (MYC, CCND3, CCNE1, TP53, CDK4/MDM2).

meta_heatmap[, `cSV Hotspots` := {
  # TRUE if any CGR column for this row indicates a positive call
  any_hit <- Reduce(`|`, lapply(.SD, function(col) col %in% c("Yes", "cSV", "CGR")))
  fifelse(any_hit, "Yes", "No")
}, .SDcols = CGR_cols]


#Analysis functions

# Pairwise purity-adjusted Firth regression across the given cSV columns.
run_cgr_firth <- function(df, cgr_cols) {
  
  cgr_pairs <- combn(cgr_cols, 2, simplify = FALSE)
  res <- data.frame()
  
  for (pair in cgr_pairs) {
    gene1 <- pair[1]
    gene2 <- pair[2]
    
    df_temp <- df[!is.na(df[[gene1]]) &
                    !is.na(df[[gene2]]) &
                    !is.na(df$purity), ]
    
    y          <- ifelse(df_temp[[gene1]] %in% c("Yes", "cSV", "CGR"), 1, 0)
    x          <- ifelse(df_temp[[gene2]] %in% c("Yes", "cSV", "CGR"), 1, 0)
    purity_val <- as.numeric(df_temp$purity)
    
    # Firth tolerates zero cells but still needs 2 levels in each predictor
    if (length(unique(y)) < 2 | length(unique(x)) < 2) next
    
    model_df <- data.frame(y = y, x = x, purity_val = purity_val)
    
    fit <- tryCatch(
      logistf(y ~ x + purity_val, data = model_df),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    
    idx_x <- which(names(coef(fit)) == "x")
    if (length(idx_x) == 0) next
    
    log_or <- coef(fit)[idx_x]   # natural-log OR (penalized)
    pval   <- fit$prob[idx_x]    # penalized likelihood-ratio p-value
    
    res <- rbind(res, data.frame(
      Gene1      = gene1,
      Gene2      = gene2,
      N          = nrow(model_df),
      N_co       = sum(x == 1 & y == 1),   # double-positive samples
      Odds_Ratio = exp(as.numeric(log_or)),
      P_Value    = pval,
      stringsAsFactors = FALSE
    ))
  }
  
  if (nrow(res) == 0) return(res)
  res$P_Adj <- p.adjust(res$P_Value, method = "BH")
  res
}

# Lower-triangle heatmap of log2(OR) with FDR significance stars.
plot_cgr <- function(results, cgr_cols, title) {
  plot_df_clean <- results %>%
    mutate(
      Gene1 = factor(Gene1, levels = cgr_cols),
      Gene2 = factor(Gene2, levels = cgr_cols),
      LogOR = log2(Odds_Ratio),
      star  = case_when(
        P_Adj < 0.001 ~ "***",
        P_Adj < 0.01  ~ "**",
        P_Adj < 0.05  ~ "*",
        TRUE          ~ ""
      )
    ) %>%
    filter(as.numeric(Gene1) < as.numeric(Gene2))   # lower triangle only
  
  ggplot(plot_df_clean, aes(x = Gene1, y = Gene2, fill = LogOR)) +
    geom_tile(color = "white", lwd = 0.8) +
    geom_text(aes(label = star), color = "white", size = 7, vjust = 0.8) +
    scale_fill_gradient2(
      low = "#1065AB", mid = "white", high = "#B31529",
      midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish,
      name   = "Odds Ratio",
      breaks = c(-2, -1, 0, 1, 2),
      labels = c("0.25", "0.5", "1", "2", "4")
    ) +
    theme_minimal() +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      panel.grid  = element_blank(),
      plot.title  = element_text(size = 16, face = "bold")
    )
}


# Analysis 1: full cohort

results_full <- run_cgr_firth(meta_heatmap, CGR_cols)
print(results_full)

p_full <- plot_cgr(results_full, CGR_cols,
                   "CGR co-occurrence and mutual exclusivity (Firth, purity-adjusted)")


p_full

# Save plot and table to local folder
ggsave(cooc_plot_out, p_full, width = 6, height = 6)
fwrite(results_full, cooc_table_out, sep = "\t")


# Analysis 2: cSV-Hotspot carriers only

CGR_hotspot_carriers <- meta_heatmap[meta_heatmap[["cSV Hotspots"]] == "Yes", ]

results_hotspot <- run_cgr_firth(CGR_hotspot_carriers, CGR_cols)
print(results_hotspot)

p_hotspot <- plot_cgr(results_hotspot, CGR_cols,
                      "CGR co-occurrence within hotspot carriers (Firth, purity-adjusted)")

p_hotspot

# Save plot and table to local folder
ggsave(mutex_plot_out, p_hotspot, width = 6, height = 6)
fwrite(results_hotspot, mutex_table_out, sep = "\t")