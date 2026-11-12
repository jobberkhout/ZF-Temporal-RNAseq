# ================================
# Control-only DESeq2 pipeline (robust, uses ENSEMBL for GO)
# ================================
# INPUT:
# Paths are relative to the repository root (open ZF-Temporal-RNAseq.Rproj, or setwd() to the repo root).
#   - counts: norm_data_controls_ALL.txt  (genes x samples, rownames=ENSDARG...)
#   - metadata: metadata/20250520_ZF_RNAseq_metadata_JB.csv     (must include: sample_id, vehicle_type, timepoint_code)
# OUTPUT: Output/Control/
#   - PCA_controls_by_vehicle_time.png
#   - DE_vehicle_DMSO_vs_E3_<T>.csv         (all genes per timepoint)
#   - DE_vehicle_DMSO_vs_E3_<T>_sig.csv     (only significant genes)
#   - GO_BP_enrichment_<T>.csv              (if any sig genes)
#   - GO_BP_dotplot_<T>.pdf                 (if any enriched terms)
#   - All_significant_vehicle_effect_genes.csv (combined)
# ================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(clusterProfiler)
  library(org.Dr.eg.db)
})

out_dir <- "Output/Control"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load data ----
counts <- read.table("norm_data_controls_ALL.txt",
                     header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
meta <- read.csv("metadata/20250520_ZF_RNAseq_metadata_JB.csv", stringsAsFactors = FALSE)

# Fix column names in counts to match metadata sample_id
colnames(counts) <- gsub("^X", "", colnames(counts))
colnames(counts) <- gsub("\\.", "-", colnames(counts))

# Keep only samples present in counts
meta_ctrl <- subset(meta, sample_id %in% colnames(counts))
rownames(meta_ctrl) <- meta_ctrl$sample_id
counts <- counts[, rownames(meta_ctrl)]

stopifnot(all(colnames(counts) == rownames(meta_ctrl)))

# ---- PCA (QC) ----
dds_all <- DESeqDataSetFromMatrix(countData = round(counts),
                                  colData = meta_ctrl,
                                  design = ~ vehicle_type)
vsd <- vst(dds_all, blind = TRUE)
p <- plotPCA(vsd, intgroup = c("vehicle_type", "timepoint_code"))
ggsave(file.path(out_dir, "PCA_controls_by_vehicle_time.png"), p, width = 6, height = 5, dpi = 300)

# ---- Per-timepoint DE + GO (using ENSEMBL IDs) ----
all_res <- list()
tp_order <- unique(meta_ctrl$timepoint_code)

for (t in tp_order) {
  sub_meta <- subset(meta_ctrl, timepoint_code == t)
  if (length(unique(sub_meta$vehicle_type)) < 2) next  # need both DMSO and E3
  
  message("Processing timepoint: ", t)
  sub_counts <- counts[, sub_meta$sample_id, drop = FALSE]
  
  sub_dds <- DESeqDataSetFromMatrix(countData = round(sub_counts),
                                    colData = sub_meta,
                                    design = ~ vehicle_type)
  sub_dds <- DESeq(sub_dds)
  
  res <- results(sub_dds, contrast = c("vehicle_type", "DMSO", "E3"))
  res <- as.data.frame(res[order(res$padj), ])
  
  # Save full DE table
  write.csv(res, file.path(out_dir, paste0("DE_vehicle_DMSO_vs_E3_", t, ".csv")))
  
  # Significant genes
  sig <- res %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
  
  if (nrow(sig) > 0) {
    sig$timepoint <- t
    write.csv(sig, file.path(out_dir, paste0("DE_vehicle_DMSO_vs_E3_", t, "_sig.csv")))
    all_res[[t]] <- sig
    
    # Ensembl IDs are the rownames (e.g., ENSDARG...); tell clusterProfiler to use ENSEMBL
    ens_ids <- rownames(sig)
    
    ego <- tryCatch(
      enrichGO(gene          = ens_ids,
               OrgDb         = org.Dr.eg.db,
               keyType       = "ENSEMBL",
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.05),
      error = function(e) NULL
    )
    
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      ego_df <- as.data.frame(ego)
      write.csv(ego_df, file.path(out_dir, paste0("GO_BP_enrichment_", t, ".csv")), row.names = FALSE)
      
      pdf(file.path(out_dir, paste0("GO_BP_dotplot_", t, ".pdf")), width = 8, height = 6)
      print(dotplot(ego, showCategory = 20))
      dev.off()
    } else {
      # write an empty file to show no enrichment found
      write.csv(data.frame(), file.path(out_dir, paste0("GO_BP_enrichment_", t, ".csv")), row.names = FALSE)
    }
  } else {
    # save an empty sig file so you know this timepoint had no hits
    write.csv(sig, file.path(out_dir, paste0("DE_vehicle_DMSO_vs_E3_", t, "_sig.csv")))
    # and an empty GO file
    write.csv(data.frame(), file.path(out_dir, paste0("GO_BP_enrichment_", t, ".csv")), row.names = FALSE)
  }
}

# ---- Combined significant table ----
if (length(all_res)) {
  combined_sig <- dplyr::bind_rows(all_res)
  write.csv(combined_sig, file.path(out_dir, "All_significant_vehicle_effect_genes.csv"), row.names = TRUE)
} else {
  write.csv(data.frame(), file.path(out_dir, "All_significant_vehicle_effect_genes.csv"), row.names = FALSE)
}

message("Pipeline complete. Outputs saved in: ", out_dir)

