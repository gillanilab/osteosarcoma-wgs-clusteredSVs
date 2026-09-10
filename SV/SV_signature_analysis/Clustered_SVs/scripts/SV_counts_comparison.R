library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# Paths
metadata_dir <- "../../../metadata"

metadata_file <- file.path(metadata_dir, "metadata_v12.tsv")
sv_counts_file <- file.path(metadata_dir, "clustered_and_nonclustered_sv_counts_by_sample.csv")

output_dir <- "../../../Images/Figures"
clustered_svg    <- file.path(output_dir, "n_clustered_count_primary_metastatic_v2.svg")
nonclustered_svg <- file.path(output_dir, "Primary_Metastatic_clinical_association_non_cSV_v2.svg")

# Load and filter metadata

metadata <- fread(metadata_file)

df_clean <- metadata %>%
  mutate(Sample_type_clean = case_when(
    Sample_type == "Primary" ~ "Primary",
    Sample_type %in% c("Metastatic", "Recurrent, Metastatic") ~ "Metastatic",
    TRUE ~ NA_character_  # blanks and "Relapse" excluded, explicitly
  )) %>%
  filter(!is.na(Sample_type_clean)) %>%
  mutate(Sample_type = factor(Sample_type_clean, levels = c("Primary", "Metastatic"))) %>%
  select(-Sample_type_clean)

# Load SV counts and join

ncSVs <- fread(sv_counts_file, sep = ",")

df_extended <- left_join(df_clean, ncSVs, by = c("Tumor_Sample_Barcode" = "sample"))

# Plot 1: Clustered SVs, Primary vs Metastatic

clustered_summary <- df_extended %>%
  filter(!is.na(n_clustered)) %>%
  group_by(Sample_type) %>%
  summarise(
    N_Samples     = n(),
    Mean_Events   = round(mean(n_clustered), 2),
    Median_Events = median(n_clustered),
    Max_Events    = max(n_clustered),
    .groups = "drop"
  )

wilcox_clustered <- wilcox.test(n_clustered ~ Sample_type, data = df_extended)

cat("\n── Clustered SVs: Primary vs Metastatic (Summary) ──\n")
print(clustered_summary)
cat("\nWilcoxon Rank-Sum Test P-value (clustered):", wilcox_clustered$p.value, "\n")

p_boxplot <- ggplot(
  df_extended %>% filter(!is.na(n_clustered)),
  aes(x = Sample_type, y = n_clustered, fill = Sample_type)
) +
  geom_boxplot(alpha = 1, outlier.shape = NA, width = 0.75, color = "black",
               box.linewidth = 0.5, whisker.linewidth = 0.5, staple.linewidth = 0.5) +
  geom_jitter(width = 0.35, alpha = 1, size = 1.5, color = "black") +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("#ccccccff", "#b22222ff")) +
  labs(
    x = "Sample Type",
    y = "Number of Clustered Events"
  ) +
  theme(
    legend.position = "none",
    plot.title  = element_text(face = "bold", size = 14),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10)
  )

print(p_boxplot)
ggsave(clustered_svg, p_boxplot, width = 2.5, height = 3, bg = "white")


# Plot 2: Non-clustered SVs, Primary vs Metastatic

nonclustered_summary <- df_extended %>%
  filter(!is.na(n_non_clustered)) %>%
  group_by(Sample_type) %>%
  summarise(
    N_Samples     = n(),
    Mean_Events   = round(mean(n_non_clustered), 2),
    Median_Events = median(n_non_clustered),
    Max_Events    = max(n_non_clustered),
    .groups = "drop"
  )

wilcox_nonclustered <- wilcox.test(n_non_clustered ~ Sample_type, data = df_extended)

cat("\nNon-Clustered SVs: Primary vs Metastatic (Summary):\n")
print(nonclustered_summary)
cat("\nWilcoxon Rank-Sum Test P-value (non-clustered):", wilcox_nonclustered$p.value, "\n")

ncSVs_boxplot <- ggplot(
  df_extended %>% filter(!is.na(n_non_clustered)),
  aes(x = Sample_type, y = n_non_clustered, fill = Sample_type)
) +
  geom_boxplot(alpha = 1, outlier.shape = NA, width = 0.75, color = "black",
               box.linewidth = 0.5, whisker.linewidth = 0.5, staple.linewidth = 0.5) +
  geom_jitter(width = 0.35, alpha = 1, size = 1.5, color = "black") +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("#ccccccff", "#b22222ff")) +
  labs(
    x = "Sample Type",
    y = "Number of Non-Clustered Events"
  ) +

  theme(
    legend.position = "none",
    plot.title  = element_text(face = "bold", size = 10),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10)
  )

print(ncSVs_boxplot)
ggsave(nonclustered_svg, ncSVs_boxplot, width = 2.5, height = 3, bg = "white")


# Combined BH-adjusted p-values across the two tests

p_adj <- p.adjust(
  c(clustered = wilcox_clustered$p.value,
    non_clustered = wilcox_nonclustered$p.value),
  method = "BH"
)

cat("\nBH-adjusted p-values (clustered vs non-clustered):\n")
print(p_adj)

