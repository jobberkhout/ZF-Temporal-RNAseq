# ===== Packages =====
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")

for (p in c("clusterProfiler","AnnotationDbi","org.Dr.eg.db","GO.db","enrichplot","DOSE")) {
  if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)
}
for (p in c("dplyr","stringr","readr","tibble","purrr","ggplot2","ggridges","ggraph","igraph")) {
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p, dependencies=TRUE)
}

# add to your CRAN installs at the top if missing
for (p in c("ggraph","igraph","ggridges","patchwork")) {
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(clusterProfiler); library(AnnotationDbi); library(org.Dr.eg.db); library(GO.db); library(enrichplot); library(DOSE)
  library(dplyr); library(stringr); library(readr); library(tibble); library(purrr); library(ggplot2)
})

# ---- single-threaded (prevents PSOCK serialize errors on Windows) ----
if (requireNamespace("BiocParallel", quietly = TRUE)) BiocParallel::register(BiocParallel::SerialParam())
options(mc.cores = 1L)
if (requireNamespace("future", quietly = TRUE)) future::plan(future::sequential)
Sys.setenv(R_FGSEA_NPROC = "1")
set.seed(42)
library(patchwork)
# ======= USER SETTINGS =======
in_dir   <- "C:\\Users\\jobbe\\RNAseq_analysis\\DEGs_NTC"
out_dir  <- "C:\\Users\\jobbe\\RNAseq_analysis\\Output\\GSEA_NTC"
file_pat <- "\\.(tsv|txt|csv)$"
minGS    <- 5
maxGS    <- 5000
include_offspring <- TRUE
sig_padj <- 0.05
gsea_padj_filter <- 0.25
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ======= NTC keywords/seeds =======
broad_pathway_terms <- c("\\bsignall?ing pathway\\b", "\\bpathway\\b", "\\bcascade\\b")
ntc_keywords <- c(
  "\\bneural tube closure\\b", "\\bneurulation\\b",
  "\\bneural plate\\b", "\\bneural fold\\b", "\\bneuropore\\b",
  "\\bneural groove\\b", "\\bneural ridge(?:s)?\\b",
  "\\bhinge point?\\b", "\\bmedian hinge point?\\b", "\\bdorsolateral hinge point?\\b",
  "\\bDLHPs?\\b", "\\bMHPs?\\b",
  "\\bfloor plate\\b", "\\broof plate\\b",
  "\\bneural crest\\b", "\\bnotochord\\b", "\\bparaxial mesoderm\\b",
  "\\bneuroectoderm\\b", "\\bsurface ectoderm\\b", "\\bectoderm\\b",
  "\\bapical constriction\\b", "\\binterkinetic nuclear migration\\b",
  "\\bnuclear migration\\b", "\\bconvergent extension\\b",
  "\\bplanar cell polarity\\b|\\bPCP\\b", "\\bfusion\\b",
  "\\bapico[- ]?basal polarity\\b|\\bapicobasal polarity\\b",
  "\\btight junctions?\\b", "\\badherens? junctions?\\b", "\\btricellular junctions?\\b",
  "\\bactin cable(?:s)?\\b", "\\bactomyosin\\b", "\\bRho GTPase\\b",
  "\\bRAC1\\b|\\bCDC42\\b|\\bRHOA\\b|\\bROCK\\b",
  "\\blamellipodia?\\b|\\bfilopodia?\\b|\\bpurse[- ]string\\b",
  "\\bcadherin\\b", "\\bE-?cadherin\\b|\\bCDH1\\b", "\\bN-?cadherin\\b|\\bCDH2\\b",
  "\\bgrainyhead[- ]like\\s*2\\b|\\bGRHL2\\b",
  "\\bnoggin\\b|\\bNOG\\b|\\bchordin\\b|\\bCHRD\\b|\\bfollistatin\\b|\\bFST\\b",
  "\\bZIC\\b|\\bZIC1\\b|\\bZIC3\\b", "\\bGLI2\\b", "\\bMSX1\\b",
  "\\bSonic hedgehog\\b|\\bSHH\\b", "\\bBMP\\b", "\\bWNT\\b|\\bWnt\\b",
  "\\bFGF\\b|\\bFGF8\\b", "\\bNotch\\b", "\\bTGF-?β\\b|\\bTGF-?b\\b|\\bSMAD\\b",
  "\\bEGF\\b|\\bERBB\\b|\\bHER\\b", "\\bVEGF\\b|\\bPDGF\\b",
  "\\bEphrin\\b|\\bEph\\b", "\\bSemaphorin\\b|\\bNetrin\\b|\\bSlit\\b|\\bRobo\\b",
  "\\bretinoic acid\\b|\\bretinoid\\b|\\bRALDH\\b|\\bCYP26\\b|\\bRDH10\\b|\\bDHRS3\\b",
  "\\bRAR\\b|\\bRXR\\b",
  "\\bdorsoventral pattern\\w*\\b|\\bdorsal[- ]ventral pattern\\w*\\b",
  "\\baxis elongation\\b|\\bbody axis extension\\b|\\bsomitogenesis\\b",
  # extras
  "\\bcili(a|um|ary)\\b", "\\bintraflagellar transport\\b|\\bIFT\\b",
  "\\baxoneme\\b", "\\bbasal body\\b", "\\btransition zone\\b",
  "\\bextracellular matrix\\b|\\bECM\\b", "\\bbasement membrane\\b",
  "\\blaminin\\b|\\bfibronectin\\b|\\bcollagen\\b",
  "\\bintegrin\\b|\\bfocal adhesion\\b",
  "\\bPar3\\b|\\bPar6\\b|\\baPKC\\b|\\bCrumbs\\b|\\bScribble\\b",
  "\\bnon[- ]?muscle myosin\\b|\\bMYH9\\b|\\bMYH10\\b|\\bMYH11\\b|\\bShroom\\b",
  "\\bERM\\b|\\bezrin\\b|\\bmoesin\\b|\\bradixin\\b",
  "\\bRhoGEF\\b|\\bARHGEF\\b|\\bRhoGAP\\b",
  "\\bVangl\\b|\\bCelsr\\b|\\bPrickle\\b|\\bPtk7\\b|\\bDchs\\b|\\bFat4\\b",
  "\\bHippo\\b|\\bYAP\\b|\\bTAZ\\b|\\bTEAD\\b"
)

ntc_goids_seed <- c(
  "GO:0001843","GO:0021993","GO:0021994","GO:0021995","GO:0021506","GO:0021507","GO:0061713","GO:0060608","GO:0001841",
  "GO:0001839","GO:0021990","GO:0021502","GO:0021503","GO:0021504","GO:0021505",
  "GO:0033504","GO:0033505","GO:0021508","GO:0021509",
  "GO:0030903","GO:0048339","GO:0007398",
  "GO:0014033","GO:0014037","GO:0001837",
  "GO:0090177","GO:0090178","GO:0060071",
  "GO:0003383","GO:0030036","GO:0045216","GO:0070830","GO:0034332","GO:0098609","GO:0030032","GO:0046847",
  "GO:0007266","GO:0051056",
  "GO:0045197",
  "GO:0022007","GO:0060029","GO:0022027","GO:0007097",
  "GO:0030509","GO:0030514","GO:0008543","GO:0007224","GO:0048384",
  "GO:0042573","GO:0002138","GO:0001523",
  # extras
  "GO:0021915","GO:0002009","GO:0060271","GO:0044782","GO:0042073",
  "GO:0030198","GO:0071711","GO:0030199","GO:0007160",
  "GO:0034330","GO:0000226","GO:0030031","GO:0035023","GO:0008360","GO:0001755"
)
ntc_goids_seed <- unique(ntc_goids_seed)

# ===== Helpers =====
read_deg <- function(path) {
  delim <- if (grepl("\\.csv$", path, ignore.case=TRUE)) "," else "\t"
  x <- suppressMessages(readr::read_delim(path, delim = delim, guess_max = 1e6, show_col_types = FALSE))
  if (!"ENSEMBL" %in% names(x)) names(x)[1] <- "ENSEMBL"
  x$ENSEMBL <- sub("\\..*$","", x$ENSEMBL)
  x
}

build_ntc_T2G <- function(all_gene_ids) {
  eg2go <- suppressMessages(
    AnnotationDbi::select(org.Dr.eg.db, keys=unique(all_gene_ids), keytype="ENSEMBL",
                          columns=c("ENSEMBL","GO","ONTOLOGY"))
  ) %>% dplyr::filter(!is.na(GO), ONTOLOGY=="BP")
  
  all_goids <- unique(eg2go$GO)
  term_map <- suppressMessages(AnnotationDbi::mapIds(GO.db::GO.db, keys=all_goids, keytype="GOID", column="TERM", multiVals="first"))
  syn_tbl <- tryCatch(
    suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys=all_goids, keytype="GOID", columns="SYNONYM")) %>%
      dplyr::group_by(GOID) %>% dplyr::summarise(SYN=paste(unique(na.omit(SYNONYM)), collapse=" | "), .groups="drop"),
    error = function(e) tibble(GOID=all_goids, SYN=NA_character_)
  )
  text_blob <- tibble(GOID = names(term_map), TERM = unname(term_map)) %>%
    dplyr::left_join(syn_tbl, by="GOID") %>% dplyr::mutate(TEXT = paste(TERM, ifelse(is.na(SYN),"",SYN)))
  
  kw_rx   <- stringr::regex(paste(ntc_keywords, collapse="|"), ignore_case=TRUE)
  path_rx <- stringr::regex(paste(broad_pathway_terms, collapse="|"), ignore_case=TRUE)
  axis_rx <- stringr::regex("(retinoic acid|retinoid|RALDH|CYP26|RDH10|DHRS3|RAR|RXR|hedgehog|SHH|BMP|WNT|FGF|NOTCH|TGF\\s*β|SMAD|ERBB|EGF|VEGF|PDGF|Eph|Ephrin|Semaphorin|Netrin|Slit|Robo)", ignore_case=TRUE)
  
  kw_hits    <- text_blob$GOID[ stringr::str_detect(text_blob$TEXT, kw_rx) ]
  broad_hits <- text_blob$GOID[ stringr::str_detect(text_blob$TERM, path_rx) & stringr::str_detect(text_blob$TERM, axis_rx) ]
  
  chosen <- unique(c(intersect(ntc_goids_seed, all_goids), kw_hits, broad_hits))
  if (include_offspring && length(chosen)) {
    kids <- unique(unlist(AnnotationDbi::mget(chosen, GO.db::GOBPOFFSPRING, ifnotfound=NA), use.names=FALSE))
    chosen <- unique(c(chosen, intersect(kids, all_goids)))
  }
  chosen <- intersect(chosen, all_goids)
  
  T2G <- eg2go %>% dplyr::filter(GO %in% chosen) %>% dplyr::distinct(GO, ENSEMBL) %>% dplyr::rename(term=GO, gene=ENSEMBL)
  T2N <- tibble(term=chosen, name=unname(term_map[chosen])) %>% dplyr::mutate(name = ifelse(is.na(name), term, name))
  list(T2G=T2G, T2N=T2N)
}

# ===== 1) Files & gene universe =====
files <- list.files(in_dir, pattern=file_pat, full.names = TRUE)
stopifnot(length(files) > 0)

message("Scanning files to collect gene universe …")
all_ids <- files %>% map(read_deg) %>% map(~ .x$ENSEMBL) %>% unlist(use.names = FALSE) %>% unique()

# ===== 2) Build custom TERM2GENE/TERM2NAME =====
message("Building NTC GO TERM2GENE …")
maps <- build_ntc_T2G(all_ids)
T2G <- maps$T2G; T2N <- maps$T2N
if (nrow(T2G) == 0) stop("No NTC GO terms found for your gene universe.")

readr::write_csv(T2G, file.path(out_dir, "NTC_GO_TERM2GENE.csv"))
readr::write_csv(T2N, file.path(out_dir, "NTC_GO_TERM2NAME.csv"))

sym_map <- suppressMessages(AnnotationDbi::mapIds(org.Dr.eg.db, keys=unique(T2G$gene),
                                                  keytype="ENSEMBL", column="SYMBOL", multiVals="first"))




# collect per-contrast barplots, grouped by chemical
barplots_sig05_by_chem <- list()

# helper: guess chemical from filename (edit this regex if your naming differs)
get_chem_from_bn <- function(bn) {
  # takes the part before first "_"  e.g. "THIR_S2_bud" -> "THIR"
  sub("_.*$", "", bn)
}

# ===== 3) Run GSEA per file =====
summary_list <- list()

for (f in files) {
  bn <- tools::file_path_sans_ext(basename(f))
  this_out <- file.path(out_dir, bn)
  dir.create(this_out, showWarnings = FALSE, recursive = TRUE)
  
  message("Processing: ", bn)
  de <- read_deg(f)
  
  # ranking vector
  rank_stat <- if ("stat" %in% names(de)) de$stat else {
    if (!all(c("log2FoldChange","pvalue") %in% names(de))) stop("Need 'stat' or both 'log2FoldChange' & 'pvalue' in: ", f)
    sign(de$log2FoldChange) * -log10(pmax(de$pvalue, 1e-300))
  }
  names(rank_stat) <- de$ENSEMBL
  
  gene_df <- tibble(gene = names(rank_stat), stat = as.numeric(rank_stat)) %>%
    dplyr::filter(!is.na(gene), !is.na(stat)) %>%
    dplyr::arrange(dplyr::desc(abs(stat))) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::arrange(dplyr::desc(stat), gene) %>%
    dplyr::group_by(stat) %>%
    dplyr::mutate(stat = stat + (dplyr::row_number() - 1L) * 1e-12) %>%
    dplyr::ungroup()
  
  geneList <- gene_df$stat; names(geneList) <- gene_df$gene
  geneList <- sort(geneList, decreasing = TRUE)
  
  sig_genes <- if ("padj" %in% names(de)) de$ENSEMBL[!is.na(de$padj) & de$padj < sig_padj] else character(0)
  
  egsea <- clusterProfiler::GSEA(geneList, TERM2GENE=T2G, TERM2NAME=T2N,
                                 minGSSize=minGS, maxGSSize=maxGS,
                                 pvalueCutoff=1, pAdjustMethod="BH", seed=TRUE, verbose=FALSE)
  
  res <- as.data.frame(egsea)
  if (nrow(res) == 0) { message("  No enriched NTC terms for: ", bn); next }
  
  term_members <- split(T2G$gene, T2G$term)
  res$n_deg_padj_lt_0_05 <- vapply(res$ID, function(id) length(intersect(term_members[[id]], sig_genes)), integer(1))
  res$deg_ENSEMBL <- vapply(res$ID, function(id) { g <- intersect(term_members[[id]], sig_genes); if (length(g)) paste(g, collapse=";") else "" }, character(1))
  res$deg_SYMBOLS <- vapply(res$ID, function(id) { g <- intersect(term_members[[id]], sig_genes); if (!length(g)) "" else paste(na.omit(unname(sym_map[g])), collapse=";") }, character(1))
  
  readr::write_csv(res %>% arrange(p.adjust), file.path(this_out, paste0(bn, "_NTC_GO_GSEA_FULL.csv")))
  readr::write_csv(res %>% filter(!is.na(p.adjust), p.adjust < gsea_padj_filter) %>% arrange(p.adjust),
                   file.path(this_out, paste0(bn, "_NTC_GO_GSEA_SIG.csv")))
  
  # ---- top-leading-edge plot ----
  res_toprow <- res %>% filter(!is.na(p.adjust)) %>% arrange(p.adjust) %>% slice(1)
  if (nrow(res_toprow) == 1) {
    top_id <- res_toprow$ID[1]
    p <- try(enrichplot::gseaplot2(egsea, geneSetID = top_id, title = res_toprow$Description[1]), silent = TRUE)
    if (!inherits(p, "try-error")) ggsave(file.path(this_out, paste0(bn, "_NTC_GO_GSEA_top.png")), p, width = 8, height = 5, dpi = 300)
  }
  
# ---- NES barplot (significant terms) FDR < 0.05 ----
sig_cut <- 0.05
res_sig <- res %>%
  dplyr::filter(!is.na(p.adjust), p.adjust < sig_cut) %>%
  dplyr::arrange(desc(abs(NES)))

if (nrow(res_sig)) {

  p_bar <- ggplot(res_sig, aes(x = reorder(Description, NES), y = NES)) +
    geom_col(aes(fill = NES > 0)) +
    scale_fill_manual(values = c(`TRUE` = "#d95f02", `FALSE` = "#1b9e77"),
                      breaks = c(FALSE, TRUE), labels = c("Down","Up"), name = "Direction") +
    coord_flip() +
    labs(x = NULL, y = "NES",
         title = bn) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 7))

  chem <- get_chem_from_bn(bn)
  if (is.null(barplots_sig05_by_chem[[chem]])) barplots_sig05_by_chem[[chem]] <- list()
  barplots_sig05_by_chem[[chem]][[bn]] <- p_bar

} else {
  message("  No terms pass FDR<", sig_cut, " for bar plot: ", bn)
}

# ---- NES barplot (significant terms) FDR < 0.05 ----
sig_cut <- 0.05
res_sig <- res %>%
  dplyr::filter(!is.na(p.adjust), p.adjust < sig_cut) %>%
  dplyr::arrange(desc(abs(NES)))

if (nrow(res_sig)) {
  
  p_bar <- ggplot(res_sig, aes(x = reorder(Description, NES), y = NES)) +
    geom_col(aes(fill = NES > 0)) +
    scale_fill_manual(values = c(`TRUE` = "#d95f02", `FALSE` = "#1b9e77"),
                      breaks = c(FALSE, TRUE), labels = c("Down","Up"), name = "Direction") +
    coord_flip() +
    labs(x = NULL, y = "NES",
         title = bn) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 7))
  
  chem <- get_chem_from_bn(bn)
  if (is.null(barplots_sig05_by_chem[[chem]])) barplots_sig05_by_chem[[chem]] <- list()
  barplots_sig05_by_chem[[chem]][[bn]] <- p_bar
  
} else {
  message("  No terms pass FDR<", sig_cut, " for bar plot: ", bn)
}

  
  # ---- NES barplot (significant terms) ----
  sig_cut <- gsea_padj_filter   # e.g., 0.25 is standard for GSEA
  res_sig <- res %>% filter(!is.na(p.adjust), p.adjust < sig_cut) %>% arrange(desc(abs(NES)))
  if (nrow(res_sig)) {
    p_bar <- ggplot(res_sig, aes(x = reorder(Description, NES), y = NES)) +
      geom_col(aes(fill = NES > 0)) +
      scale_fill_manual(values = c(`TRUE` = "#d95f02", `FALSE` = "#1b9e77"),
                        breaks = c(FALSE, TRUE), labels = c("Down","Up"), name = "Direction") +
      coord_flip() +
      labs(x = NULL, y = "NES", title = paste0(bn, " — Significant NTC GO terms (FDR < ", sig_cut, ")")) +
      theme_minimal() + theme(axis.text.y = element_text(size = 9))
    ggsave(file.path(this_out, paste0(bn, "_NTC_GO_GSEA_bar_sig25.png")),
           p_bar, width = 8, height = max(3, 0.3*nrow(res_sig)), dpi = 300)
  } else {
    message("  No terms pass FDR<", sig_cut, " for bar plot: ", bn)
  }
  
  # ---- Top-N dot plot (by FDR) ----
  topN <- 25L
  res_top <- dplyr::slice_min(res, order_by = p.adjust, n = min(topN, nrow(res)), with_ties = FALSE)
  if (nrow(res_top)) {
    p_dot <- ggplot(res_top, aes(x = NES, y = reorder(Description, NES))) +
      geom_point(aes(size = -log10(p.adjust), fill = NES), shape = 21, alpha = 0.9, stroke = 0.2) +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, name = "NES") +
      geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
      labs(y = NULL, x = "NES", size = "-log10(FDR)", title = paste0(bn, " — Top ", nrow(res_top), " NTC GO terms")) +
      theme_minimal()
    ggsave(file.path(this_out, paste0(bn, "_NTC_GO_GSEA_dot_top.png")),
           p_dot, width = 8, height = max(3, 0.35*nrow(res_top)), dpi = 300)
  }
  
  # ---- Top-K leading-edge plots (robust replacement for ridgeplot) ----
  k <- min(6L, nrow(res))  # how many sets to show
  if (k > 0) {
    ids <- res$ID[order(res$p.adjust)][seq_len(k)]
    # A multipage PDF (one page per pathway)
    pdf(file.path(this_out, paste0(bn, "_NTC_GO_GSEA_topK.pdf")), width = 8, height = 5)
    for (id in ids) {
      ttl <- res$Description[match(id, res$ID)]
      print(enrichplot::gseaplot2(egsea, geneSetID = id, title = ttl))
    }
    dev.off()
  }
  
  # ---- Enrichment map & cnet ----
  if (nrow(res) >= 2 && requireNamespace("ggraph", quietly = TRUE)) {
    egsea_sim <- try(enrichplot::pairwise_termsim(egsea), silent = TRUE)
    if (!inherits(egsea_sim, "try-error")) {
      # draw emap only if we have some similarity edges
      sim_mat <- try(egsea_sim@termsim, silent = TRUE)
      has_edges <- !inherits(sim_mat, "try-error") && any(sim_mat[upper.tri(sim_mat, diag=FALSE)] > 0, na.rm = TRUE)
      if (has_edges) {
        em <- try(enrichplot::emapplot(egsea_sim, showCategory = min(40L, nrow(res)), layout = "kk"), silent = TRUE)
        if (!inherits(em, "try-error")) {
          ggsave(file.path(this_out, paste0(bn, "_NTC_GO_GSEA_emap.png")), em, width = 9, height = 7, dpi = 300)
        } else message("  Skipping emapplot: ", attr(em, "condition")$message)
      } else message("  Skipping emapplot: no term-term similarity edges.")
    } else message("  Skipping term similarity: ", attr(egsea_sim, "condition")$message)
    
    fc <- setNames(de$log2FoldChange, de$ENSEMBL)
    cn <- try(enrichplot::cnetplot(egsea, foldChange = fc, showCategory = min(10L, nrow(res))), silent = TRUE)
    if (!inherits(cn, "try-error")) {
      ggsave(file.path(this_out, paste0(bn, "_NTC_GO_GSEA_cnet.png")), cn, width = 9, height = 7, dpi = 300)
    } else message("  Skipping cnetplot: ", attr(cn, "condition")$message)
  } else message("  Skipping emap/cnet: need >=2 terms and package 'ggraph'.")
  
  # add to summary
  summary_list[[bn]] <- res %>% mutate(file = bn)
}

# ===== 3b) Save one combined page per chemical (FDR < 0.05 barplots) =====
out_appendix <- file.path(out_dir, "APPENDIX_barplots_sig05_byChemical")
dir.create(out_appendix, showWarnings = FALSE, recursive = TRUE)

for (chem in names(barplots_sig05_by_chem)) {
  plots <- barplots_sig05_by_chem[[chem]]
  if (length(plots) == 0) next
  
  # layout: 2 columns works well for “many contrasts”
  ncol <- 2
  nrow <- ceiling(length(plots) / ncol)
  
  combined <- wrap_plots(plots, ncol = ncol) +
    plot_annotation(
      title = paste0(chem, " — NTC GO GSEA (FDR < 0.05)"),
      theme = theme(plot.title = element_text(face = "bold"))
    )
  
  # save as ONE-PAGE PDF per chemical (nice for appendix)
  pdf(file.path(out_appendix, paste0(chem, "_GSEA_bar_sig05_ALL.pdf")),
      width = 14, height = max(7, 3.5 * nrow), onefile = TRUE)
  print(combined)
  dev.off()
  
  # optional PNG too
  ggsave(file.path(out_appendix, paste0(chem, "_GSEA_bar_sig05_ALL.png")),
         combined, width = 14, height = max(7, 3.5 * nrow), dpi = 300)
}

message("Appendix pages saved in: ", out_appendix)

# ===== 4) Combined summary (then heatmap) =====
if (length(summary_list)) {
  combo <- bind_rows(summary_list)
  
  readr::write_csv(combo %>% arrange(file, p.adjust),
                   file.path(out_dir, "ALL_FILES_NTC_GO_GSEA_FULL.csv"))
  readr::write_csv(combo %>% filter(!is.na(p.adjust), p.adjust < gsea_padj_filter) %>% arrange(file, p.adjust),
                   file.path(out_dir, "ALL_FILES_NTC_GO_GSEA_SIG.csv"))
  
  # ---- Cross-contrast heatmap 0.05 ----
  library(tidyr)
  mat <- combo %>%
    filter(!is.na(p.adjust)) %>%
    mutate(score = sign(NES) * -log10(p.adjust)) %>%
    group_by(Description) %>%
    filter(any(p.adjust < 0.05, na.rm = TRUE)) %>%
    ungroup() %>%
    select(Description, file, score) %>%
    tidyr::pivot_wider(names_from = file, values_from = score, values_fill = 0)
  
  if (nrow(mat) > 0) {
    mat_long <- mat %>% tidyr::pivot_longer(-Description, names_to = "file", values_to = "score")
    hp <- ggplot(mat_long, aes(x = file, y = Description, fill = score)) +
      geom_tile() +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                           name = "sign(NES)·-log10(FDR)") +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.grid = element_blank()) +
      labs(x = NULL, y = NULL, title = "NTC GO — disrupted terms across contrasts")
    ggsave(file.path(out_dir, "ALL_FILES_NTC_GO_GSEA_heatmap05.png"),
           hp, width = max(8, 0.18 * (ncol(mat) - 1)), height = max(4, 0.25 * nrow(mat)), dpi = 300)
  } else {
    message("No terms passed FDR<0.05 across files — not drawing heatmap.")
  }
  
}

message("Done. Results in: ", out_dir)
