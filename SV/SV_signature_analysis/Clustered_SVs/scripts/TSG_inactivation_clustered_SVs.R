################################################################################
## Figure 3C
## Purity-corrected TSG inactivation vs cSV enrichment
## Firth's penalized logistic regression
################################################################################

library(data.table) 
library(dplyr)       
library(stringr)     
library(ggplot2)     
library(logistf)      

#metadata_file   <- "../metadata/2026_06_08_Metadata_v11.tsv"
metadata_file   <- "../../../metadata/metadata_v12.tsv" #" 2026_06_08_Metadata_v11.tsv #v12 is union and v13 is intersection
firth_table_out <- "../../../Images/Figures/Tables/Table_Firth_purity_corrected_gene_vs_cSV.txt"
firth_plot_out  <- "../../../Images/Figures/TSG_cSV_Heatmap_Firth_corrected.pdf"



# Load metadata

metadata <- fread(metadata_file)
meta_heatmap <- copy(metadata)


# Gene / ordering configuration
# Define TSG included in Figure
tsg_genes <- c("ATRX", "RB1", "PTEN", "CDKN2A", "NF2") 

# Define Heatmap ordering
target_order <- c("cSV overall", "TP53 cSV", "MYC cSV", "CCND3 cSV", "CCNE1 cSV",
                  "CDK4/MDM2 cSV")
# Gene order from top to bottom
gene_order <- rev(c("ATRX", "RB1", "PTEN", "CDKN2A", "NF2")) 


# Build required gene columns

# 1.) # Drop multi-tumor samples: including them would let one patient's
# recurrent mutation "vote" more than once in downstream enrichment stats.
if ("filter status" %in% names(meta_heatmap)) {
  meta_heatmap <- meta_heatmap[`filter status` != "multiple-tumor"]
}

# 2.) Recode mutation flags Yes/No -> Mutation/No Mutation
# Identify all _Mutation_yes_No columns
mut_yes_no_cols <- paste0(tsg_genes, "_Mutation_yes_no")
existing_mut_cols <- mut_yes_no_cols[mut_yes_no_cols %in% names(meta_heatmap)]
# Relabel yes and no into Mutation and No Mutation for better understanding in table
meta_heatmap[, (existing_mut_cols) := lapply(.SD, function(x) {
  fifelse(x == "Yes", "Mutation", "No Mutation")
}), .SDcols = existing_mut_cols]

# Derive gene names from all "<gene>_SV_gene_body" columns in metadata.
# For each gene, copy the gene-body SV call directly into a new
# "<gene>_SV_binary" column (assuming it already holds only "SV"/"No")

sv_genes <- gsub("_SV_gene_body", "",
                 grep("_SV_gene_body", names(metadata), value = TRUE))
for (gene in sv_genes) {
  gb_col  <- paste0(gene, "_SV_gene_body")
  bin_col <- paste0(gene, "_SV_binary")
  
  meta_heatmap[, (bin_col) := get(gb_col)]  
}


# Build cSV_mat 

# Does a sample has clustered SVs or not and encode it in column named 'cSV'
meta_heatmap <- meta_heatmap %>%
  mutate("cSV" = case_when(
    n_clustered_events == 0 ~ "No",
    n_clustered_events >= 1 ~ "Yes"
  ))


# Assemble the cSV annotation matrix (which sample has a cSV at which locus)
# and reformat labels (renaming of CGR to cSV)
cSV_mat <- meta_heatmap %>%
  select("cSV", grep("_CGR$", colnames(meta_heatmap), value = TRUE)) %>%
  mutate(across(ends_with("_CGR"), ~ str_replace_all(., "CGR", "cSV"))) %>%
  mutate(across(ends_with("_CGR"), ~ str_replace_all(., "noCGR", "no cSV"))) %>%
  rename(
    "cSV overall"    = "cSV",
    "MYC cSV"        = "MYC_CGR",
    "CCND3 cSV"      = "CCND3_CGR",
    "CCNE1 cSV"      = "CCNE1_CGR",
    "TP53 cSV"       = "TP53_CGR",
    "CDK4/MDM2 cSV"  = "CDK4-MDM2_CGR"
  )


# Transform into dataframe
cSV_mat <- as.data.frame(cSV_mat)
# Annotate Sample names again to cSV mat
rownames(cSV_mat) <- meta_heatmap$Tumor_Sample_Barcode

# cSV targets tested in the Firth models (saved into variables)
specific_csv_cols <- grep(" cSV$", colnames(cSV_mat), value = TRUE)
all_csv_targets   <- c("cSV overall", specific_csv_cols)


# Step 1: Classify Any_Hit vs No_Hit per gene
# A sample is "No_Hit" only with no mutation, no SV, and retained copy number
# (Total_CN >= 1.5). Everything else is "Any_Hit".

for (gene in tsg_genes) {
  
  mut_col <- paste0(gene, "_Mutation_yes_no")
  sv_col  <- paste0(gene, "_SV_gene_body") #_SV_binary
  tcn_col <- paste0(gene, "_total_cn")
  
  curr_mut <- if (mut_col %in% colnames(meta_heatmap)) meta_heatmap[[mut_col]] else rep("No", nrow(meta_heatmap))
  curr_sv  <- if (sv_col  %in% colnames(meta_heatmap)) meta_heatmap[[sv_col]]  else rep("No", nrow(meta_heatmap))
  curr_tcn <- if (tcn_col %in% colnames(meta_heatmap)) meta_heatmap[[tcn_col]] else rep(NA,   nrow(meta_heatmap))
  
  is_no_hit <- curr_mut != "Mutation" &
    curr_sv != "SV" &
    (is.na(curr_tcn) | curr_tcn >= 1.5)
  
  meta_heatmap[[paste0(gene, "_AnyHit_vs_NoHit")]] <- case_when(
    !is_no_hit ~ "Any_Hit",
    is_no_hit  ~ "No_Hit",
    TRUE       ~ NA_character_
  )
  
  cat("\n", gene, ":\n")
  print(table(meta_heatmap[[paste0(gene, "_AnyHit_vs_NoHit")]], useNA = "always"))
}

hit_summary <- sapply(tsg_genes, function(gene) {
  col <- paste0(gene, "_AnyHit_vs_NoHit")
  n_any  <- sum(meta_heatmap[[col]] == "Any_Hit", na.rm = TRUE)
  n_no   <- sum(meta_heatmap[[col]] == "No_Hit", na.rm = TRUE)
  n_na   <- sum(is.na(meta_heatmap[[col]]))
  n_total <- nrow(meta_heatmap)
  
  c(Any_Hit_n = n_any, Any_Hit_pct = round(100 * n_any / n_total, 1),
    No_Hit_n  = n_no,  No_Hit_pct  = round(100 * n_no  / n_total, 1),
    NA_n = n_na)
})

# Define a short hit summary to see if correctly classified as expected
hit_summary <- as.data.frame(t(hit_summary))
print(hit_summary)

# Step 2: Firth regression
# For every gene x cSV-target pair, fit an unadjusted and a purity-adjusted
# Firth model, extracting the csv_status effect as log2(OR).

results_firth <- list()

for (gene in tsg_genes) {
  
  combined_col <- paste0(gene, "_AnyHit_vs_NoHit")
  if (!combined_col %in% colnames(meta_heatmap)) next
  
  gene_status <- meta_heatmap[[combined_col]]
  purity_vec  <- as.numeric(meta_heatmap$purity)
  
  for (target in all_csv_targets) {
    if (!target %in% colnames(cSV_mat)) next
    
    csv_array <- ifelse(cSV_mat[[target]] %in% c("Yes", "cSV", 1, TRUE), 1, 0)
    
    model_df <- data.frame(
      outcome    = as.integer(gene_status == "Any_Hit"),
      csv_status = csv_array,
      purity     = purity_vec
    ) %>%
      filter(!is.na(outcome), !is.na(purity), outcome %in% c(0, 1))
    
    # Require variation in both predictor and outcome, and a minimum n
    if (length(unique(model_df$csv_status)) < 2) next
    if (length(unique(model_df$outcome))    < 2) next
    if (nrow(model_df) < 10) next
    
    # Firth model: purity-adjusted
    fit_adj <- tryCatch(
      logistf(outcome ~ csv_status + purity, data = model_df),
      error = function(e) NULL
    )
    
    if (is.null(fit_adj)) next
    
    idx_adj <- which(names(coef(fit_adj)) == "csv_status")
    idx_pur <- which(names(coef(fit_adj)) == "purity")
    if (length(idx_adj) == 0) next
    
    # Coefficients are log-odds; convert to log2(OR) for the heatmap
    beta_adj <- coef(fit_adj)[idx_adj]
    
    results_firth[[paste0(gene, "_vs_", target)]] <- data.frame(
      Gene        = gene,
      cSV_Target  = target,
      N_total     = nrow(model_df),
      N_AnyHit    = sum(model_df$outcome == 1),
      N_cSV       = sum(model_df$csv_status == 1),
      # Purity-adjusted (Firth)
      Log2_OR_adj = beta_adj / log(2),
      CI_low_adj  = fit_adj$ci.lower[idx_adj] / log(2),   # profile-penalized CI, log2
      CI_high_adj = fit_adj$ci.upper[idx_adj] / log(2),
      P_adj_model = fit_adj$prob[idx_adj],
      # Purity's own effect on the outcome
      P_purity    = if (length(idx_pur)) fit_adj$prob[idx_pur] else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}


# Step 3: FDR-correct within each model, flag changes

final_firth <- do.call(rbind, results_firth) %>%
  mutate(
    P_Adj_adj = p.adjust(P_adj_model, method = "BH"),
    star_adj  = case_when(P_Adj_adj < 0.01 ~ "**",
                          P_Adj_adj < 0.05 ~ "*", TRUE ~ "")
  ) %>%
  arrange(P_Adj_adj)

print(head(final_firth, 20))

# Make any Inf/-Inf Excel-safe before export (shouldn't occur with Firth)
final_firth_export <- final_firth %>%
  mutate(across(where(is.numeric), ~ ifelse(is.infinite(.), NA, .)))
# Export result table with exact values
fwrite(final_firth_export, firth_table_out, sep = "\t")


# Step 4: Heatmap of purity-adjusted log2(OR)
# Generates Heatmap from Firth regression results

# Isolate columns required for plotting from firth results
plot_firth <- final_firth %>%
  mutate(
    cSV_Target = factor(cSV_Target, levels = intersect(target_order, unique(cSV_Target))),
    Gene       = factor(Gene,       levels = intersect(gene_order,   unique(Gene)))
  )


# Plot 
p_firth <- ggplot(plot_firth, aes(x = cSV_Target, y = Gene)) +
  geom_tile(aes(fill = Log2_OR_adj), width = 1, height = 1,
            color = "white", linewidth = 0.8) +
  geom_text(aes(label = star_adj), color = "white", size = 5, vjust = 0.75) +
  scale_fill_gradient2(
    low = "#1065AB", mid = "white", high = "#B31529",
    midpoint = 0, limits = c(-4, 4), oob = scales::squish,
    name = "Log2(OR)\npurity-adjusted (Firth)", na.value = "grey90"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),
    plot.title  = element_text(face = "bold", size = 13)
  ) +
  labs(
    title    = "TSG Inactivation Enrichment in cSV Backgrounds",
    subtitle = "Firth penalized logistic regression, purity-adjusted | * FDR < 0.05, ** FDR < 0.01",
    x = "cSV Target", y = "Gene"
  )

# Print plot
p_firth

# Save Plot
ggsave(firth_plot_out, p_firth, width = 10, height = 3)



# ATRX vs cSV Overall: Contingency Table & Summary
# Testing if ATRX has 0 samples in N

library(dplyr)
library(tidyr)

# 1. Filter the dataset for valid ATRX and cSV overall status
atrx_df <- meta_heatmap %>%
  filter(
    !is.na(ATRX_AnyHit_vs_NoHit),
    !is.na(cSV)
  ) %>%
  mutate(
    # Clean up names for presentation
    ATRX_Status = factor(ATRX_AnyHit_vs_NoHit, levels = c("Any_Hit", "No_Hit")),
    cSV_Cohort  = factor(ifelse(cSV == "Yes", "cSV Cohort", "No cSV Cohort"), 
                         levels = c("No cSV Cohort", "cSV Cohort"))
  )

# 2. Generate raw counts and column percentages (what % of each cohort has an ATRX hit)
atrx_table_summary <- atrx_df %>%
  group_by(cSV_Cohort, ATRX_Status) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(
    Total_Cohort = sum(Count),
    Percentage   = round((Count / Total_Cohort) * 100, 1),
    `Count (%)`  = paste0(Count, " (", Percentage, "%)")
  )

# 3. Reshape into a clean wide 2x2 table
atrx_2x2_wide <- atrx_table_summary %>%
  select(ATRX_Status, cSV_Cohort, `Count (%)`) %>%
  pivot_wider(names_from = cSV_Cohort, values_from = `Count (%)`, values_fill = "0 (0%)")

cat("\n=======================================================\n")
cat(" ATRX Inactivation Status by cSV Cohort\n")
cat("=======================================================\n")
print(as.data.frame(atrx_2x2_wide), row.names = FALSE)

# 4. Run quick Fisher's Exact Test for raw numbers check
raw_matrix <- table(atrx_df$ATRX_Status, atrx_df$cSV_Cohort)
cat("\n--- Raw Counts Matrix ---\n")
print(raw_matrix)

cat("\n--- Unadjusted Fisher's Exact Test ---\n")
print(fisher.test(raw_matrix))

