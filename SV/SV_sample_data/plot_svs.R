if (!requireNamespace("ReConPlot", quietly = TRUE)) {
  devtools::install_github("cortes-ciriano-lab/ReConPlot")
}

library(ggplot2)
library(ReConPlot)
library(data.table)

# patch
ReConPlot_patched <- ReConPlot
body_lines <- deparse(body(ReConPlot_patched))
body_lines[232] <- '            chr_selection$chr[i]])) {'
body_lines[234] <- '                chr_selection$chr[i]])'
body(ReConPlot_patched) <- parse(text = paste(body_lines, collapse = "\n"))


TARGET_SAMPLE <- "SJOS001107_M2-SJOS001107_G1"
chrom_map <- c("23" = "X", "24" = "Y")

# load SVs
sv_raw <- fread("data/sv_file_v2.csv")
sv_filtered <- sv_raw[sv_raw$name == TARGET_SAMPLE, ]

# reformat to ReConPlot input
sv_data <- data.frame(
  chr1    = paste0("chr", sv_filtered$chr1),
  pos1    = sv_filtered$pos1,
  chr2    = paste0("chr", sv_filtered$chr2),
  pos2    = sv_filtered$pos2,
  strands = paste0(ifelse(sv_filtered$str1 == 1, "+", "-"),
                   ifelse(sv_filtered$str2 == 1, "+", "-"))
)
sv_data$chr1 <- gsub("chr23", "chrX", gsub("chr24", "chrY", sv_data$chr1))
sv_data$chr2 <- gsub("chr23", "chrX", gsub("chr24", "chrY", sv_data$chr2))
print(head(sv_data))

# load CNVs
cnv_raw <- read.csv("data/cnv_segments.csv")
cnv_filtered <- cnv_raw[cnv_raw$sample == TARGET_SAMPLE, ]
cn_data <- data.frame(
  chr                   = paste0("chr", cnv_filtered$Chromosome),
  start                 = cnv_filtered$Start.bp,
  end                   = cnv_filtered$End.bp,
  copyNumber            = round(cnv_filtered$rescaled.cn.a1 + cnv_filtered$rescaled.cn.a2),
  minorAlleleCopyNumber = round(pmin(cnv_filtered$rescaled.cn.a1, cnv_filtered$rescaled.cn.a2))
)
cn_data$chr <- gsub("chr23", "chrX", gsub("chr24", "chrY", cn_data$chr))
print(head(cn_data))

#Chromosome selection
chrs=c("chr8")
chr_selection = data.frame(
  chr=chrs,
  start= 50000000,
  end= 145138636) 
print(chr_selection)

# plot
p = ReConPlot_patched(sv_data,
  cn_data,
  scale_ticks=15000000,
  chr_selection=chr_selection,
  genes = c("MYC", "PTEN"),
  legend_SV_types=T,
  scaling_cn_SVs = 1/4,
  max.cn = 14,
  karyotype_rel_size = 0.1,
  pos_SVtype_description=1000000,
  scale_separation_SV_type_labels=1/15,
  genome_version="hg38",
  curvature_intrachr_SVs = -0.08,
  npc_now = 0.00625 * 6
)


p$layers[[1]]$aes_params$fill <- "white"
p$layers[[2]]$aes_params$fill <- "white"

old_to_new <- c(
  "black"       = "#5588e0ff",
  "forestgreen" = "#e07a5fff",
  "darkblue"    = "#3d405bff",
  "orange"      = "#e05555ff"
)

curve_idx <- which(sapply(p$layers, function(l) class(l$geom)[1] == "GeomCurve"))

for (i in curve_idx) {
  old_col <- p$layers[[i]]$aes_params$colour
  if (!is.null(old_col) && old_col %in% names(old_to_new)) {
    p$layers[[i]]$aes_params$colour <- unname(old_to_new[old_col])
  }
}

for (i in curve_idx) {
    cat("\nLayer", i, "\n")
    cat("Color:", p$layers[[i]]$aes_params$colour, "\n")
    print(unique(p$layers[[i]]$data$strands))
}
#text_idx <- which(sapply(p$layers, function(l) class(l$geom)[1] == "GeomText"))
#if (length(text_idx) > 0) {
#  p$layers[text_idx] <- NULL
#}
#ls(envir = asNamespace("ReConPlot"))
#karyotype <- get("karyotype_data", envir = asNamespace("ReConPlot"))
#karyotype[karyotype$chr == "chr13", ]
#getAnywhere(ReConPlot)
# get the function source
ggsave(filename = paste0("plots/", TARGET_SAMPLE, "_ReConPlot.pdf"), plot = p, width = 20, height = 6, units = "cm")
