##############################################################
#  OncoPrint – MutSig Driver Genes (TP53, ATRX, RB1, PTEN, NF2)
##############################################################

library(ComplexHeatmap)
library(data.table)
library(tidyverse)
library(circlize)

# Set working directory
setwd("...")
# Adapt path to metadata for plotting Mutsig significant gene results
# Preperation is as in done in Data_preperation_for_oncoprint.ipynb
metadata <- fread('../../../metadata/metadata_v12.tsv')

# 1. Driver genes (determined from MutSig2CV)
driver_genes <- c("TP53", "ATRX", "RB1", "PTEN", "NF2")

# 2. Mutation type annotation (oncogenic filter)
genes_with_mut   <- gsub("_Mutation_Status", "",
                         grep("_Mutation_Status", names(metadata), value = TRUE))
genes_to_process <- intersect(driver_genes, genes_with_mut)

meta_heatmap <- as.data.table(metadata)

for (gene in genes_to_process) {
  mut_col <- paste0(gene, "_Mutation_Status")
  okb_col <- paste0(gene, "_OncoKB_Status")
  new_col <- paste0(gene, "_Mutation_Heatmap")
  
  meta_heatmap[, (new_col) := {
    mapply(function(status, okb) {
      if (is.na(status) | status == "No Mutation") return("No Mutation")
      
      m_vec <- trimws(unlist(strsplit(as.character(status), ",")))
      o_vec <- trimws(unlist(strsplit(as.character(okb),    ",")))
      if (length(o_vec) < length(m_vec)) o_vec <- rep(o_vec, length(m_vec))
      
      results <- sapply(seq_along(m_vec), function(i) {
        m <- m_vec[i]; o <- o_vec[i]
        is_oncogenic <- o %in% c("Oncogenic", "Likely Oncogenic")
        if (grepl("Missense|Nonsense|Frame_Shift|Splice_Site", m) & is_oncogenic) return(m)
        if (grepl("3'UTR|5'UTR|Silent|Intron|Flank", m) & is_oncogenic) return(m)
        return(NA)
      })
      
      valid_muts <- results[!is.na(results)]
      if (length(valid_muts) == 0) return("No Mutation")
      return(paste(valid_muts, collapse = "; "))
      
    }, get(mut_col), get(okb_col))
  }]
}

# Clean "_Mutation" suffix from cell values (Frame_Shift_Del etc.)
mut_heatmap_cols <- grep("_Mutation_Heatmap", names(meta_heatmap), value = TRUE)
meta_heatmap[, (mut_heatmap_cols) := lapply(.SD, function(x) gsub("_Mutation", "", x)),
             .SDcols = mut_heatmap_cols]

# 3. Remove multi-tumor samples
meta_heatmap <- meta_heatmap[`filter status` != "multiple-tumor"]

# 4. Build mutation-only matrix (genes × samples)
driver_mat <- do.call(rbind, lapply(driver_genes, function(gene) {
  col_name <- paste0(gene, "_Mutation_Heatmap")
  if (!col_name %in% names(meta_heatmap)) return(rep("", nrow(meta_heatmap)))
  vals <- meta_heatmap[[col_name]]
  ifelse(is.na(vals) | vals == "No Mutation", "", vals)   # NA → blank; "No Mutation" kept as-is
}))
rownames(driver_mat) <- driver_genes
colnames(driver_mat) <- meta_heatmap$Tumor_Sample_Barcode

# 5. Sort columns: oncoprint-style gene priority
gene_priority <- c("TP53", "RB1", "ATRX", "PTEN", "NF2")

priority_score <- sapply(colnames(driver_mat), function(s) {
  for (i in seq_along(gene_priority)) {
    g <- gene_priority[i]
    val <- driver_mat[g, s]
    if (!is.na(val) && val != "") return(i)
  }
  return(length(gene_priority) + 1)
})

num_real_muts <- colSums(driver_mat != "")
col_order     <- order(priority_score, -num_real_muts)



# 8. Right Annotation (Frequency Barplot)
right_ann <- rowAnnotation(
  "Frequency" = anno_oncoprint_barplot(
    show_fraction = TRUE,  # Forces the 0.0 to 1.0 percentage scale
    axis_param = list(side = "bottom", labels_rot = 0)
  ),
  annotation_name_side = "bottom",
  annotation_name_gp   = gpar(fontsize = 10)
)


# 9. Alteration colours & alter_fun (mutations only)
col <- c(
  "Missense"        = "#6aa84fff",
  "Nonsense"        = "#ee876cff",
  "Splice_Site"     = "#3d85c6ff",
  "Frame_Shift_Del" = "#b22222ff",
  "Frame_Shift_Ins" = "#674ea7ff"
)

alter_fun <- list(
  # White background = no data / not in matrix
  background = function(x, y, w, h)
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = "#F5F5F5", col = NA)),
  
  # Point mutations: full-height coloured block
  Missense = function(x, y, w, h)
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = col["Missense"], col = NA)),
  Nonsense = function(x, y, w, h)
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = col["Nonsense"], col = NA)),
  Splice_Site = function(x, y, w, h)
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = col["Splice_Site"], col = NA)),
  Frame_Shift_Del = function(x, y, w, h)
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA)),
  Frame_Shift_Ins = function(x, y, w, h)
    grid.rect(x, y, w - unit(0.3, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
)

# 10. Build OncoPrint
op <- oncoPrint(
  driver_mat,
  alter_fun               = alter_fun,
  col                     = col,
  get_type                = function(x) strsplit(x, ";")[[1]],
  column_order            = col_order,
  row_order               = gene_priority,
  top_annotation          = NULL,          # Set to top_ann if you want Clinical data shown
  bottom_annotation       = NULL,    # cSV BOTTOM bottom_ann
  right_annotation        = right_ann,     # Frequency RIGHT
  alter_fun_is_vectorized = FALSE,
  column_title            = "OncoPrint \u2013 MutSig Driver Genes (Mutations Only)",
  show_column_names       = FALSE,
  
  # --- Percentage setup ---
  row_names_side          = "left",
  show_pct                = TRUE,          
  pct_side                = "right",
  
  column_gap              = unit(2, "mm"),
  row_gap                 = unit(4, "mm"),
  remove_empty_columns    = FALSE,
  remove_empty_rows       = FALSE
)

# 11. Save
pdf("../Images/Figures/OncoPrint_MutSig_Drivers_with_Clinical_cSV_Frequency.pdf", width = 14, height = 4.5)
draw(op,
     merge_legends          = TRUE,
     heatmap_legend_side    = "bottom",
     annotation_legend_side = "bottom")
dev.off()

op



