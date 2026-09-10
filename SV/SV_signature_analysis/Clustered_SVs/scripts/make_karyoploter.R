if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("karyoploteR", quietly = TRUE))
  BiocManager::install("karyoploteR")

library(karyoploteR)
library(GenomicRanges)

sv_density <- read.csv("../results/hotspot_windows.csv")
#chromothripsis_density <- read.csv("../../../Shatterseek/results/chromothripsis_windows.csv")
ec_density <- read.csv("../../../CNAs_Amplicon_Architect/results/ecdna_windows.csv")
hotspots <- read.csv("../results/hotspot_event_definitions.csv")
out_dir <- "../results"

recode_chr <- function(x) {
  x <- as.character(x)
  x[x == "23"] <- "X"
  x[x == "24"] <- "Y"
  paste0("chr", x)
}

sv_gr <- GRanges(
  seqnames = recode_chr(sv_density$chromosome),
  ranges   = IRanges(sv_density$window_start_bp + 1, sv_density$window_end_bp),
  score    = sv_density$clustered_breakpoints
)
if (FALSE) {
chromothripsis_gr <- GRanges(
  seqnames = chromothripsis_density$chrom,
  ranges   = IRanges(chromothripsis_density$window_start_bp + 1, chromothripsis_density$window_end_bp),
  score    = chromothripsis_density$n_samples
)
}
ec_gr <- GRanges(
  seqnames = ec_density$chrom,
  ranges   = IRanges(ec_density$window_start_bp + 1, ec_density$window_end_bp),
  score    = ec_density$n_samples
)

plot_params <- getDefaultPlotParams(plot.type = 2)
plot_params$data1height <- 200
#plot_params$data2height <- 100
plot_params$ideogramheight <- 30

sv_max <- max(sv_gr$score)
ec_max <- max(ec_gr$score)
# chromothripsis_max <- max(chromothripsis_gr$score)

# left: chr1-11
pdf(file.path(out_dir, "karyotype_left_svs_only_v2.pdf"), width = 8, height = 6)
kp_l <- plotKaryotype(genome = "hg38", plot.type = 2,
                      chromosomes = paste0("chr", 1:11),
                      plot.params = plot_params)

for (i in 1:nrow(hotspots)) {
  chrom <- paste0("chr", hotspots$chrom[i])
  kpRect(kp_l, chr = chrom,
         x0 = hotspots$start[i], x1 = hotspots$end[i],
         y0 = 0.25, y1 = 1,
         col = "#e07a5f65", border = NA,
         data.panel = "all")
}

kpArea(kp_l, data = sv_gr, y = sv_gr$score, ymax = sv_max, col = "#b22222ff", border = NA, data.panel = 1)
#kpArea(kp_l, data = chromothripsis_gr, y = chromothripsis_gr$score, ymax = chromothripsis_max, col = "#b22222ff", border = NA, data.panel = 1)
#kpArea(kp_l, data = ec_gr, y = ec_gr$score, ymax = ec_max, col = "#3d85c6ff", border = NA, data.panel = 2)

dev.off()

plot_params <- getDefaultPlotParams(plot.type = 2)
plot_params$data1height <- 200
#plot_params$data2height <- 100
plot_params$ideogramheight <- 30

# right: chr22 to chr12, then chrX
pdf(file.path(out_dir, "karyotype_right_svs_only_v2.pdf"), width = 10, height = 6)
kp_r <- plotKaryotype(genome = "hg38", plot.type = 2,
                      chromosomes = c(paste0("chr", 22:12), "chrX"),
                      plot.params = plot_params)

for (i in 1:nrow(hotspots)) {
  chrom <- paste0("chr", hotspots$chrom[i])
  kpRect(kp_r, chr = chrom,
         x0 = hotspots$start[i], x1 = hotspots$end[i],
         y0 = 0.25, y1 = 1,
         col = "#e07a5f65", border = NA,
         data.panel = "all")
}

kpArea(kp_r, data = sv_gr, y = sv_gr$score, ymax = sv_max, col = "#b22222ff", border = NA, data.panel = 1)
#kpArea(kp_r, data = chromothripsis_gr, y = chromothripsis_gr$score, ymax = chromothripsis_max, col = "#b22222ff", border = NA, data.panel = 1)
#kpArea(kp_r, data = ec_gr, y = ec_gr$score, ymax = ec_max, col = "#3d85c6ff", border = NA, data.panel = 2)

# add labels manually on the right side
kpAddChromosomeNames(kp_r,
                     chr.names = c(as.character(22:12), "X"),
                     xoffset = 1,       # push to right end (0=left, 1=right)
                     adj = c(0, 0.5))   # left-justify text at that right position

dev.off()

# chr17 ideogram only with TP53 hotspot highlight

if (FALSE) {
pdf(file.path(out_dir, "karyotype_chr17.pdf"), width = 8, height = 2)
kp_17 <- plotKaryotype(genome = "hg38", plot.type = 1, chromosomes = "chr17")

tp53 <- hotspots[hotspots$chrom == 17, ]
kpRect(kp_17, chr = "chr17",
       x0 = tp53$start, x1 = tp53$end,
       y0 = -0.3, y1 = 1.3,
       col = "#e07a5f65", border = NA,
       data.panel = "ideogram")

dev.off()
}