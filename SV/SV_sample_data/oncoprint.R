################################################################################
## OncoPrint of Oncogenes and Tumor Suppressor Genes (Figure 1B)
##
## Generates:
##   - OncoPrint of oncogenes and TSGs across cohorts (CNV / LOH / SNV / SV)
##   - Demographic stacked bar plots (sex, disease stage, age group)
################################################################################

# Install Bioconductor package manager and ComplexHeatmap if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap")
}

# CRAN packages used in this script
cran_packages <- c("data.table", "dplyr", "forcats", "ggplot2",
                   "circlize", "patchwork")
missing <- cran_packages[!vapply(cran_packages, requireNamespace,
                                 logical(1), quietly = TRUE)]
if (length(missing) > 0) install.packages(missing)

library(ComplexHeatmap)
library(data.table)
library(dplyr)
library(forcats)
library(ggplot2)
library(circlize)
library(patchwork)


# 1. Set working directory to local path:
setwd('')

# Load files required for Oncoprint 
metadata_file <- "../../../metadata/metadata_v12.tsv"
sn_rna_file   <- "../../../metadata/2026_samples_with_sn_rna.tsv" # List of samples with snRNA data available
bulk_rna_file <- "../../../metadata/2026_samples_with_bulk_rna.tsv" # List of samples with bulk RNA available

# Save output path and names to variable
oncoprint_pdf     <- "../../../Images/Figures/OncoPrint_correct.pdf"
demographics_pdf  <- "../../../Images/Figures/Demographics_StackedBars.pdf"
demographics_svg  <- "../../../Images/Figures/Demographics_StackedBars.svg"


# Load metadata and annotation files

metadata <- fread(metadata_file)
sn_rna   <- fread(sn_rna_file,   header = FALSE)
bulk_rna <- fread(bulk_rna_file, header = FALSE)


# Prepare metadata

# Flag samples that have bulk RNA / snRNA data for the manuscript
metadata <- metadata %>%
  mutate(
    bulk_RNA_manuscript = Tumor_Sample_Barcode %in% bulk_rna$V1,
    snRNA_manuscript    = Tumor_Sample_Barcode %in% sn_rna$V1
  )


# Copy number alterations

# Collapse the detailed CNV status strings into broader biological categories
collapse_cnv <- function(x) {
  fct_collapse(x,
               "Deep Deletion" = c("Deep Deletion"),
               "Deletion"      = c("Het. D.", "Het. D. (LOH)",
                                   "Severe Post-WGD loss (LOH)", "Severe Post-WGD Loss"),
               "CN-LOH"        = c("CN-LOH"),
               "Neutral"       = c("Neutral"),
               "Amplification" = c("Focal Gain", "Amplification"),
               other_level     = "Neutral"   # handles NAs and any remaining levels
  )
}

meta_heatmap <- metadata %>%
  mutate(across(ends_with("_CNV_Status"), collapse_cnv, .names = "{.col}_Heatmap"))


# Mutations (SNV)

# Recode per-gene mutation flags from Yes/No to SNV / No SNV for the legend
all_genes <- c("TP53", "RB1", "ATRX", "PTEN", "CDKN2A",
               "MYC", "CCND3", "CCNE1", "CDK4", "MDM2", "NF2")

mut_yes_no_cols  <- paste0(all_genes, "_Mutation_yes_no")
existing_mut_cols <- mut_yes_no_cols[mut_yes_no_cols %in% names(meta_heatmap)]

meta_heatmap[, (existing_mut_cols) := lapply(.SD, function(x) {
  fifelse(x == "Yes", "SNV", "No SNV")
}), .SDcols = existing_mut_cols]


# Clinical annotation

# Complex genomic rearrangement (CGR) flag from number of clustered events
meta_heatmap <- meta_heatmap %>%
  mutate(CGR = case_when(
    n_clustered_events == 0 ~ "No",
    n_clustered_events >= 1 ~ "Yes"
  ))

# Age group (2 categories, NA handled)
meta_heatmap <- meta_heatmap %>%
  mutate("Age Group" = case_when(
    Age_at_diagnosis <= 25 ~ "<= 25",
    Age_at_diagnosis  > 25 ~ ">25",
    TRUE                   ~ "NA"
  ))


# Empty ecDNA strings -> NA
meta_heatmap <- meta_heatmap %>%
  rename(any_of(c("any_ecDNA_events" = "has_ecDNA"))) %>%
  mutate(any_ecDNA_events = na_if(any_ecDNA_events, ""),
         any_ecDNA_events = case_when(
           tolower(any_ecDNA_events) == "yes" ~ "Yes",
           tolower(any_ecDNA_events) == "no"  ~ "No",
           TRUE                               ~ NA_character_
         ))

# Assemble the clinical annotation table
clinic_annot <- meta_heatmap %>%
  mutate(
    Sample_type = recode(Sample_type, "Recurrent, Metastatic" = "Metastatic"),
    study       = recode(study,
                         "DFCI_from_CZ_Zhang" = "DFCI",
                         "DFCI_with_snRNAseq" = "DFCI")
  ) %>%
  select("study", "sex", "Age Group", "Sample_type", "purity",
         "Genome doublings", "any_ecDNA_events", "Multi_Tumor_Label",
         "bulk_RNA_manuscript", "snRNA_manuscript") %>%
  rename(
    "Cohort"        = "study",
    "Sex"           = "sex",
    "Disease stage" = "Sample_type",
    "Purity"        = "purity",
    "WGD"           = "Genome doublings",
    "ecDNA"         = "any_ecDNA_events",
    "Multi-sample"  = "Multi_Tumor_Label",
    "RNA"           = "bulk_RNA_manuscript",
    "snRNA"         = "snRNA_manuscript"
  )

clinic_annot <- as.data.frame(clinic_annot)
rownames(clinic_annot) <- meta_heatmap$Tumor_Sample_Barcode

# Standardise types and fill missing values with the literal "NA" label
clinic_annot$Purity <- suppressWarnings(as.numeric(as.character(clinic_annot$Purity)))
clinic_annot$WGD    <- as.character(clinic_annot$WGD)
clinic_annot[is.na(clinic_annot)] <- "NA"
clinic_annot[clinic_annot == ""]  <- "NA"
# Re-coerce Purity to numeric (the blanket NA fill above coerced it to character)
clinic_annot$Purity <- suppressWarnings(as.numeric(as.character(clinic_annot$Purity)))

# Annotation colour palettes
clin_cols <- list(
  Cohort = c(
    "dbGaP" = "#ee876cff", "DFCI" = "#3d85c6ff", "ICGC" = "#6aa84fff",
    "MD_Anderson" = "#674ea7ff", "StJude" = "#b22222ff",
    "TARGET" = "#ccccccff", "NA" = "#d9d9d9ff"
  ),
  Sex = c(
    "Male" = "#9fc5e8ff", "Female" = "#ee876cda",
    "Unknown" = "#d9d9d9ff", "NA" = "#d9d9d9ff"
  ),
  "Age Group" = c(
    "<= 25" = "#b6d7a8ff", ">25" = "#6aa84fff", "NA" = "#d9d9d9ff"
  ),
  "Disease stage" = c(
    "Primary" = "#ee876cff", "Metastatic" = "#a61c00ff",
    "Relapse" = "#b4a7d6ff", "NA" = "#d9d9d9ff"
  ),
  Purity = colorRamp2(c(0.15, 0.5, 0.98),
                      c("#f7fbff", "#9fc5e8ff", "#3d85c6ff")),
  WGD    = c("0" = "#f7f7f7", "1" = "#b4a7d6ff", "2" = "#674ea7ff", "NA" = "#d9d9d9ff"),
  ecDNA  = c("Yes" = "#674ea7ff", "No" = "#d9d9d9ff", "NA" = "#8c8c8c"),
  RNA    = c("TRUE" = "black", "FALSE" = "#d9d9d9ff", "NA" = "#8c8c8c"),
  snRNA  = c("TRUE" = "black", "FALSE" = "#d9d9d9ff", "NA" = "#8c8c8c")
)

clin_cols$`SV data` <- c("Available" = "#d9d9d9ff", "NA" = "#8c8c8c")

# Top annotation build

# Multi-sample text track: keep only patients that appear more than once
clean_multi_text <- trimws(as.character(clinic_annot$`Multi-sample`))
clean_multi_text[is.na(clean_multi_text)] <- ""
clean_multi_text[clean_multi_text %in% c("NA", "No", "NaN")] <- ""

standard_cols <- setdiff(names(clinic_annot), "Multi-sample")

# SV counts for the top barplot but label samples without sv sequencing data
sv_counts <- as.numeric(meta_heatmap$number_of_SVs)
sv_counts[is.na(sv_counts)]  <- 0

sv_cap <- 1500
sv_counts <- pmin(sv_counts, sv_cap)

# Build the top annotation in the exact display order
annotation_args <- c(
  # 1. SV barplot at the top
  list("SVs" = anno_barplot(
    sv_counts,
    height     = unit(1.5, "cm"),
    gp         = gpar(fill = "#a7a7a8", col = "transparent"),
    axis_param = list(gp = gpar(fontsize = 8)),
    border     = FALSE
  )),
  
  # 2. Standard clinical columns
  as.list(clinic_annot[, standard_cols]),
  
  # Spacer
  list("  " = anno_empty(height = unit(1, "mm"), border = FALSE)),
  
  # 3. Multi-sample text track
  list("Multi-sample" = anno_text(
    clean_multi_text,
    rot = 0, just = "center", location = 0.5,
    gp = gpar(fontsize = 8)
  )),
  
  # Spacer
  list(" " = anno_empty(height = unit(2, "mm"), border = FALSE)),
  
  # 4. Formatting parameters
  list(
    col                  = clin_cols,
    na_col               = "#d9d9d9ff",
    annotation_name_side = "left",
    annotation_name_gp   = gpar(fontsize = 11),
    simple_anno_size     = unit(5, "mm"),
    gp                   = gpar(col = "white", lwd = 1),
    annotation_legend_param = list(
      Cohort      = list(ncol = 2),
      "Age Group" = list(nrow = 1)
    )
  )
)

top_ann <- do.call(HeatmapAnnotation, annotation_args)


# Oncogene / TSG gene matrices

oncogenes <- c("MYC", "CCND3", "CCNE1", "CDK4", "MDM2")
tumor_suppressors <- c("TP53", "RB1", "ATRX", "CDKN2A", "PTEN")

# Build a per-gene x per-sample matrix combining CNV, LOH, SNV and SV events.
# TSGs and oncogenes follow different rules for which events are "hits".
process_gene_matrix_smart <- function(input_df, gene_list, is_oncogene = FALSE) {
  df <- as.data.frame(input_df)
  
  mat_rows <- lapply(gene_list, function(gene) {
    cnv_col <- paste0(gene, "_CNV_Status_Heatmap")
    loh_col <- paste0(gene, "_LOH")
    snv_col <- paste0(gene, "_Mutation_yes_no")
    sv_col  <- paste0(gene, "_SV_gene_body")
    
    sapply(seq_len(nrow(df)), function(i) {
      
      cnv <- if (cnv_col %in% names(df)) df[[cnv_col]][i] else NA
      loh <- if (loh_col %in% names(df)) df[[loh_col]][i] else NA
      snv <- if (snv_col %in% names(df)) df[[snv_col]][i] else NA
      sv  <- if (sv_col  %in% names(df)) df[[sv_col]][i]  else NA
      
      has_loh <- isTRUE(loh)
      
      cnv <- ifelse(is.na(cnv) | cnv == "", "Neutral", as.character(cnv))
      #loh <- ifelse(is.na(loh) | loh == "", "No",      as.character(loh))
      snv <- ifelse(is.na(snv) | snv == "", "No",      as.character(snv))
      sv  <- ifelse(is.na(sv)  | sv  == "", "No",      as.character(sv))
      
      has_snv <- snv == "SNV"
      has_sv  <- sv  == "SV"
      
      
      base_state <- character(0)
      
      if (cnv %in% c("Amplification", "Focal Gain")) {
        # Amplification/gain shown only for oncogenes; neutral for TSGs
        if (is_oncogene) base_state <- "Amplification"
        
      } else if (cnv == "Deep Deletion") {
        base_state <- "Deep Deletion"
        
      } else if (cnv == "Deletion") {
        base_state <- "Loss (CN < 2)/ LOH"
        
      } else if (!is_oncogene) {
        # TSG-specific rules
        if (cnv == "CN-LOH") {
          # CN-LOH counts only with an accompanying functional event
          if (has_snv | has_sv) base_state <- "CN-LOH"
        } else if (has_loh) {
          # Plain LOH counts only with an accompanying functional event
          if (has_snv | has_sv) base_state <- "Loss (CN < 2)/ LOH"
        }
      }
      # Oncogenes: CN-LOH / plain LOH are not shown (loss events)
      
      events <- c(base_state)
      if (has_snv) events <- c(events, "SNV")
      if (has_sv)  events <- c(events, "SV")
      
      if (length(events) == 0) return("")
      paste(unique(events), collapse = ";")
    })
  })
  
  res_mat <- do.call(rbind, mat_rows)
  rownames(res_mat) <- gene_list
  colnames(res_mat) <- df$Tumor_Sample_Barcode
  res_mat
}

tsg_complheat_mat  <- process_gene_matrix_smart(meta_heatmap, tumor_suppressors, is_oncogene = FALSE)
onco_complheat_mat <- process_gene_matrix_smart(meta_heatmap, oncogenes,         is_oncogene = TRUE)


# Complex genomic rearrangement (CGR) matrix

CGR_mat <- meta_heatmap %>%
  select("CGR", grep("_CGR$", colnames(meta_heatmap))) %>%
  rename(
    "any"        = "CGR",
    "MYC"        = "MYC_CGR",
    "CCND3"      = "CCND3_CGR",
    "CCNE1"      = "CCNE1_CGR",
    "TP53"       = "TP53_CGR",
    "CDK4/MDM2"  = "CDK4-MDM2_CGR"
  ) %>%
  mutate(across(everything(), ~ case_when(
    .x %in% c("Yes", "CGR")   ~ "cSV",
    .x %in% c("No", "noCGR")  ~ "no cSV",
    TRUE                      ~ NA_character_
  )))

CGR_mat <- as.data.frame(CGR_mat)
rownames(CGR_mat) <- meta_heatmap$Tumor_Sample_Barcode


# Combine gene matrices

final_mutation_matrix <- rbind(onco_complheat_mat, tsg_complheat_mat)

# Row split separating oncogenes from TSGs
gene_row_split <- factor(
  c(rep("Oncogenes", nrow(onco_complheat_mat)),
    rep("TSGs",      nrow(tsg_complheat_mat))),
  levels = c("Oncogenes", "TSGs")
)


# Right annotation: Number of samples with priorization

n_samples <- ncol(final_mutation_matrix)   # 236

# Priority order used ONLY to pick one "primary" event per sample per gene,
# so each sample contributes to exactly one stacked-bar segment (no over-counting).
# Order encodes clinical severity: deletions/CN-loss outrank point mutations,
# which outrank SVs, matching the same logic as score_tsg()/score_oncogene().
priority_order <- c("Deep Deletion", "Loss (CN < 2)/ LOH", "CN-LOH",
                    "Amplification", "SNV", "SV")

assign_primary_type <- function(cell) {
  if (cell == "" || is.na(cell)) return(NA_character_)
  events <- strsplit(cell, ";")[[1]]
  hit <- priority_order[priority_order %in% events]
  if (length(hit) == 0) return(NA_character_)
  hit[1]   # highest-priority event wins; sample counted once
}

# Build a genes x category count matrix (each row sums to <= 236)
stack_categories <- priority_order
count_mat <- t(apply(final_mutation_matrix, 1, function(gene_row) {
  primary <- vapply(gene_row, assign_primary_type, character(1))
  tab <- table(factor(primary, levels = stack_categories))
  as.numeric(tab)
}))
colnames(count_mat) <- stack_categories
rownames(count_mat) <- rownames(final_mutation_matrix)

# Sanity check while developing (remove once confirmed):
# stopifnot(all(rowSums(count_mat) <= n_samples))

right_ann <- rowAnnotation(
  "# Altered" = anno_barplot(
    count_mat,
    ylim       = c(0, n_samples),
    width      = unit(2.5, "cm"),
    gp         = gpar(fill = col[stack_categories], col = "transparent"),
    axis_param = list(
      at     = c(0, 50, 100, 150, 200, n_samples),
      labels = c("0", "50", "100", "150", "200", as.character(n_samples)),
      gp     = gpar(fontsize = 8)
    ),
    border = FALSE
  ),
  annotation_name_side = "bottom",
  annotation_name_gp   = gpar(fontsize = 7)
)


col <- c(
  "Amplification"      = "#C14E4E",
  "Deep Deletion"      = "#3d85c6ff",
  "Loss (CN < 2)/ LOH" = "#9fc5e8ff",
  "CN-LOH"             = "#9fc5e8ff",
  "SNV"                = "#ee876cff",
  "SV"                 = "#f5dc5d"
)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = "#EEEEEE", col = NA))
  },
  "Amplification" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = col["Amplification"], col = NA))
  },
  "Deep Deletion" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = col["Deep Deletion"], col = NA))
  },
  "Loss (CN < 2)/ LOH" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = col["Loss (CN < 2)/ LOH"], col = NA))
  },
  "CN-LOH" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = col["CN-LOH"], col = NA))
  },
  # SNV: thin bar slightly above centre
  "SNV" = function(x, y, w, h) {
    grid.rect(x, y + h * 0.15, w - unit(0.5, "mm"), h * 0.25,
              gp = gpar(fill = col["SNV"], col = NA))
  },
  # SV: thin bar slightly below centre
  "SV" = function(x, y, w, h) {
    grid.rect(x, y - h * 0.15, w - unit(0.5, "mm"), h * 0.25,
              gp = gpar(fill = col["SV"], col = NA))
  }
)


# Alteration priority scores (lower = higher priority for sorting)
score_tsg <- function(gene_row) {
  case_when(
    grepl("Deep Deletion",      gene_row)             ~ 1,
    grepl("Loss (CN < 2)/ LOH", gene_row, fixed = TRUE) ~ 2,
    grepl("CN-LOH",             gene_row)             ~ 3,
    grepl("SV",                 gene_row)             ~ 4,
    grepl("SNV",                gene_row)             ~ 5,
    TRUE                                              ~ 6
  )
}

score_oncogene <- function(gene_row) {
  case_when(
    grepl("Amplification", gene_row) ~ 1,
    TRUE                             ~ 2
  )
}

# cSV priority
csv_prio <- case_when(
  CGR_mat$any == "cSV"    ~ 1,
  CGR_mat$any == "no cSV" ~ 2,
  is.na(CGR_mat$any)      ~ 3
)

# Cohort rank for the outer sort (must match column_split ordering)
cohort_levels <- c("StJude", "TARGET", "ICGC", "dbGaP", "MD_Anderson", "DFCI")
cohort_rank   <- as.integer(factor(clinic_annot$Cohort,
                                   levels = cohort_levels, exclude = NULL))
cohort_rank[is.na(cohort_rank)] <- length(cohort_levels) + 1

# Per-gene priority scores
tp53_prio  <- score_tsg(final_mutation_matrix["TP53",  ])
rb1_prio   <- score_tsg(final_mutation_matrix["RB1",   ])
atrx_prio  <- score_tsg(final_mutation_matrix["ATRX",  ])
myc_prio   <- score_oncogene(final_mutation_matrix["MYC",   ])
ccnd3_prio <- score_oncogene(final_mutation_matrix["CCND3", ])
num_alts   <- colSums(final_mutation_matrix != "")

# Primary sort: cohort first (aligns with column_split), then alteration priorities
custom_order <- order(
  cohort_rank,
  csv_prio,
  tp53_prio,
  rb1_prio,
  atrx_prio,
  myc_prio,
  ccnd3_prio,
  -num_alts
)

# Keep multi-sample patients adjacent for readability
patient_id <- clean_multi_text
patient_id[patient_id == ""] <- NA
tab        <- table(patient_id)
shared_ids <- names(tab[tab > 1])

ord     <- custom_order
ord_pid <- patient_id[ord]

# Collapse shared-patient samples to their first appearance; singletons stay put
first_pos <- ave(seq_along(ord), ord_pid, FUN = function(i) min(i))
group_key <- ifelse(ord_pid %in% shared_ids, first_pos, seq_along(ord))
local_order          <- order(group_key, seq_along(ord))  # stable tie-break
custom_order_grouped <- ord[local_order]


# Plotting

op <- oncoPrint(
  final_mutation_matrix,
  alter_fun               = alter_fun,
  col                     = col,
  get_type                = function(x) strsplit(x, ";")[[1]],
  column_order            = custom_order_grouped,
  top_annotation          = top_ann,
  right_annotation        = right_ann,
  bottom_annotation       = NULL,
  row_split               = gene_row_split,
  alter_fun_is_vectorized = FALSE,
  column_split            = factor(clinic_annot$Cohort, levels = cohort_levels),
  row_order               = c(oncogenes, tumor_suppressors),
  column_title            = "OncoPrint of Oncogenes and TSGs",
  show_column_names       = FALSE,
  row_names_side          = "left",
  pct_side                = "right",
  column_gap              = unit(2, "mm"),
  row_gap                 = unit(2, "mm"),
  remove_empty_columns    = FALSE,
  remove_empty_rows       = FALSE
)

pdf(oncoprint_pdf, width = 30, height = 7)
draw(op,
     merge_legends          = TRUE,
     heatmap_legend_side    = "bottom",
     annotation_legend_side = "bottom")
dev.off()


# Demographic bar plots

# Order categories by frequency, forcing "NA" to the end
clinic_annot <- clinic_annot %>%
  mutate(
    Sex             = fct_relevel(fct_infreq(Sex), "NA", after = Inf),
    `Disease stage` = fct_relevel(fct_infreq(`Disease stage`), "NA", after = Inf),
    `Age Group`     = fct_relevel(fct_infreq(`Age Group`), "NA", after = Inf)
  )

# Shared theme for the demographic bars
demo_theme <- theme_classic() +
  theme(
    axis.text.y     = element_text(face = "bold", size = 12, color = "black"),
    axis.text.x     = element_text(size = 10, color = "black"),
    legend.position = "right",
    legend.title    = element_text(size = 8, face = "bold"),
    legend.text     = element_text(size = 7),
    legend.key.size = unit(4, "mm")
  )

p_sex <- ggplot(clinic_annot, aes(y = "Sex", fill = Sex)) +
  geom_bar(position = position_stack(reverse = TRUE),
           color = "black", linewidth = 0.2, width = 0.5) +
  scale_fill_manual(values = clin_cols$Sex) +
  labs(x = "Number of Samples", y = NULL) +
  demo_theme

p_stage <- ggplot(clinic_annot, aes(y = "Disease Stage", fill = `Disease stage`)) +
  geom_bar(position = position_stack(reverse = TRUE),
           color = "black", linewidth = 0.2, width = 0.5) +
  scale_fill_manual(values = clin_cols$`Disease stage`) +
  labs(x = "Number of Samples", y = NULL) +
  demo_theme

p_age <- ggplot(clinic_annot, aes(y = "Age Group", fill = `Age Group`)) +
  geom_bar(position = position_stack(reverse = TRUE),
           color = "black", linewidth = 0.2, width = 0.5) +
  scale_fill_manual(values = clin_cols$`Age Group`) +
  labs(x = "Number of Samples", y = NULL) +
  demo_theme

demographics_plot <- p_sex / p_stage / p_age

pdf(demographics_pdf, width = 4, height = 3)
print(demographics_plot)
dev.off()

svg(demographics_svg, width = 4, height = 2)
print(demographics_plot)
dev.off()

