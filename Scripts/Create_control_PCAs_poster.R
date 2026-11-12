# Run Pipeline_compare_controls.R first: this script uses the `vsd` object it creates.
# Paths are relative to the repository root (open ZF-Temporal-RNAseq.Rproj, or setwd() to the repo root).
library(ggplot2)
library(ggrepel)
library(dplyr)
library(DESeq2)

# Extract PCA data
pca_df <- plotPCA(vsd, intgroup = c("timepoint_code", "vehicle_type"), returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))
pca_df$sample_id <- rownames(pca_df)

# Build convex hulls per timepoint
hulls <- pca_df %>%
  group_by(timepoint_code) %>%
  slice(chull(PC1, PC2))

# Define custom pastel palette
pastel_palette <- c("#C1A3CD", "#604A7B", "#B0E298", "#FC98A9", "#8596CC")

# === Base PCA Plot with large text and lines ===
p_pub <- ggplot(pca_df, aes(PC1, PC2, color = timepoint_code, shape = vehicle_type)) +
  stat_ellipse(aes(group = timepoint_code),
               level = 0.90, linewidth = 2.5, alpha = 0.3) +
  geom_point(size = 10, stroke = 2) +
  scale_color_manual(values = pastel_palette) +
  labs(
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Timepoint", shape = "Vehicle",
    title = "Temporal dynamics of neural tube closure (controls)"
  ) +
  theme_bw(base_size = 40) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", size = 42),
    axis.text = element_text(size = 36),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 38),
    legend.text = element_text(size = 34),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 50),
    axis.ticks = element_line(linewidth = 2),
    panel.border = element_rect(linewidth = 2)
  )

# === Labeled version ===
p_lab <- p_pub +
  ggrepel::geom_text_repel(
    aes(label = sample_id),
    size = 10,             # increased for poster size
    max.overlaps = 50,
    seed = 42
  )

# === Convex hull version ===
p_hull <- ggplot(pca_df, aes(PC1, PC2, color = timepoint_code, shape = vehicle_type)) +
  geom_polygon(data = hulls, aes(group = timepoint_code),
               fill = NA, linewidth = 2.5, alpha = 0.5) +
  geom_point(size = 10, stroke = 2) +
  scale_color_manual(values = pastel_palette) +
  labs(
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Timepoint", shape = "Vehicle",
    title = "Temporal dynamics of neural tube closure (controls)"
  ) +
  theme_bw(base_size = 40) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", size = 42),
    axis.text = element_text(size = 36),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 38),
    legend.text = element_text(size = 34),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 50),
    axis.ticks = element_line(linewidth = 2),
    panel.border = element_rect(linewidth = 2)
  )

# === Save plots in high resolution ===
ggsave("Output/PCA/PCA_controls_pub_poster.png", p_pub,
       width = 48, height = 16, dpi = 600)

ggsave("Output/PCA/PCA_controls_pub_labeled_poster.png", p_lab,
       width = 48, height = 16, dpi = 600)

ggsave("Output/PCA/PCA_controls_pub_hulls_poster.png", p_hull,
       width = 48, height = 16, dpi = 600)
