# --- Self-contained pastel control PCA (builds vsd itself; no prior script needed) ---
# Reproduces the control vsd the way Pipeline_compare_controls.R does, then plots
# the PCA with the named pastel timepoint palette (matches the GSVA heatmaps).
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)

# ---- paths (relative to the repository root: open ZF-Temporal-RNAseq.Rproj, or setwd() to the repo root) ----
counts_path <- "ZF_RNA_Raw_Count.txt"
meta_path   <- "metadata/20250520_ZF_RNAseq_metadata_JB.csv"
out         <- "Output/Figures/"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---- build vsd from control samples ----
counts <- read.table(counts_path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
meta   <- read.csv(meta_path, stringsAsFactors = FALSE)
colnames(counts) <- gsub("^X", "", colnames(counts))
colnames(counts) <- gsub("\\.", "-", colnames(counts))

meta_ctrl <- subset(meta, dose_code == "C" & sample_id %in% colnames(counts))
rownames(meta_ctrl) <- meta_ctrl$sample_id
cnt <- round(counts[, meta_ctrl$sample_id, drop = FALSE])
stopifnot(all(colnames(cnt) == rownames(meta_ctrl)))

dds <- DESeqDataSetFromMatrix(countData = cnt, colData = meta_ctrl, design = ~ vehicle_type)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)

# ---- pastel PCA (named mapping = identical to the GSVA heatmaps) ----
pastel <- c("T1" = "#C1A3CD", "T2" = "#604A7B", "T3" = "#B0E298",
            "T4" = "#FC98A9", "T5" = "#8596CC")

pca_df <- plotPCA(vsd, intgroup = c("timepoint_code", "vehicle_type"), returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))
pca_df$sample_id <- rownames(pca_df)

p_notitle <- ggplot(pca_df, aes(PC1, PC2, color = timepoint_code, shape = vehicle_type)) +
  stat_ellipse(aes(group = timepoint_code), level = 0.90, linewidth = 0.6, alpha = 0.4) +
  geom_point(size = 3.5, stroke = 0.6) +
  scale_color_manual(values = pastel) +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"),
       color = "Timepoint", shape = "Vehicle") +
  theme_bw(base_size = 14) + theme(panel.grid = element_blank())

p_labeled <- p_notitle +
  ggrepel::geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 50, seed = 42)

# PDF = vector (highest quality). If your R lacks cairo, drop `device = cairo_pdf`.
ggsave(paste0(out, "PCA_controls_pastel_notitle.pdf"), p_notitle, width = 12, height = 4, device = cairo_pdf)
ggsave(paste0(out, "PCA_controls_pastel_labeled.pdf"), p_labeled, width = 12, height = 4, device = cairo_pdf)
