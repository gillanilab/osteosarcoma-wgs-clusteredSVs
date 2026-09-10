library(ComplexHeatmap)
library(data.table)
library(dplyr)
library(forcats)
library(ggplot2)
library(circlize)
library(patchwork)


# Set working directory to local path:
setwd('')

# Load files required for Oncoprint 
metadata_file <- "../../../metadata/metadata_v12.tsv"
metadata <- fread(metadata_file)

# Generate a column which classifies if sample has a clustered SV
metadata <- metadata %>% 
  mutate("cSV" = case_when(n_clustered_events == 0 ~ "No",
                           n_clustered_events >= 1 ~ "Yes"))

# Generate Plot of Purity on y axis and clustered SV yes vs no on x axis
p_cSV_purity <- metadata %>%
  mutate(cSV = factor(cSV, levels = c("Yes", "No"))) %>%
  ggplot(aes(x = cSV, y = purity, fill = cSV)) +
  geom_boxplot(alpha = 1, position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.75), 
              size = 1, alpha = 0.8, color = "black") +
  theme_bw() +
  scale_fill_manual(values = c("Yes" = "#9fc5e8ff", "No" = "#d9d9d9ff")) + 
  labs(
    title = "Sample Purity by cSV", 
    y = "Sample Purity",
    fill = "cSV Status",
    x = "cSV Status"
  ) +
  theme(
    plot.title = element_text(size = 16),
    strip.background = element_rect(fill = "#EEEEEE"),   
    axis.text.x = element_text(angle = 0, size = 12), 
    legend.position = "right", 
    legend.title = element_text(face = "bold")
  )

p_cSV_purity

# Testing significance
wilcox.test(purity ~ cSV, data = metadata)

# Save the Plot
ggsave("../Images/Figures/Purity_of_cSV_samples_vs_no_cSV.svg",
       p_cSV_purity, width = 3, height = 5, device = "svg")