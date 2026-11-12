# Run Pipeline_compare_controls.R first: this script uses the `vsd` object it creates.
# Paths are relative to the repository root (open ZF-Temporal-RNAseq.Rproj, or setwd() to the repo root).
# --- Assumes dds and vsd already exist from your DESeq2 workflow ---

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

# === Plot 1: Publication-style PCA with ellipses ===
p_pub <- ggplot(pca_df, aes(PC1, PC2,
                            color = timepoint_code,
                            shape = vehicle_type)) +
  stat_ellipse(aes(group = timepoint_code), level = 0.90,
               linewidth = 0.6, alpha = 0.4) +
  geom_point(size = 3.5, stroke = 0.6) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       color = "Timepoint", shape = "Vehicle",
       title = "Controls: PCA (ellipses)") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5))

# === Plot 2: Same but labeled ===
p_lab <- p_pub +
  ggrepel::geom_text_repel(aes(label = sample_id),
                           size = 3, max.overlaps = 50, seed = 42)

# === Save both plots ===
ggsave("Output/PCA/PCA_controls_pub.png", p_pub,
       width = 12, height = 4, dpi = 300)
ggsave("Output/PCA/PCA_controls_pub_labeled.png", p_lab,
       width = 12, height = 4, dpi = 300)

# (Optional) === Plot 3: Convex hull version ===
p_hull <- ggplot(pca_df, aes(PC1, PC2, color = timepoint_code, shape = vehicle_type)) +
  geom_polygon(data = hulls, aes(group = timepoint_code), fill = NA, linewidth = 0.6) +
  geom_point(size = 3.5, stroke = 0.6) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       color = "Timepoint", shape = "Vehicle",
       title = "Controls: PCA (convex hulls)") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("Output/PCA/PCA_controls_pub_hulls.png", p_hull,
       width = 12, height = 4, dpi = 300)
