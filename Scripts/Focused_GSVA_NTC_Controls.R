# ============================================================
# Focused GSVA for ALL signaling/signalling pathways
#   - Works on raw counts (DESeq2 VST) or normalized data
#   - Picks all GO:BP terms whose TERM matches "signaling/signalling pathway"
#   - Summarizes ssGSEA scores by timepoint within each group
#   - Produces per-group heatmaps + ONE combined heatmap across groups
#   - Optional delta (treated - vehicle) heatmap if controls found
#   - Exports long-format CSV for custom visualizations
# ============================================================

# -----------------
# USER SETTINGS
# -----------------
counts_path <- "D:/RNAseq/Raw_data/Job_Raw_Count.txt"
meta_path   <- "D:/RNAseq/20250520_ZF_RNAseq_metadata_JB.csv"
out_dir     <- "D:/RNAseq/Focused_GSVA_Selected_GO_NTC"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# What column defines the analysis groups? (vehicle/compound/treatment/…)
group_col_preferred <- c("vehicle_type","vehicle","compound_code","treatment","group","Group","Label","label")

# Optionally limit to specific groups, e.g. c("DMSO","E3") or c("ChemX","ChemY")
groups_include <- NULL

# Timepoint column (factor T1..T10). If missing, inferred from sample_id "-T#".
timepoint_col_preferred <- c("timepoint_code","timepoint","tp")

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

## Final controls used by your script:
focus_keywords    <- c(broad_pathway_terms, ntc_keywords)
focus_go_ids      <- unique(ntc_goids_seed)
include_offspring <- TRUE      # pull in children of seed GO IDs (great for NTC coverage)
require_GOdb      <- TRUE      # keep TRUE for name mapping & offspring expansion


# GSVA gene set size limits
gsva_min_size <- 5
gsva_max_size <- 2000

# Heatmap settings
order_pathways_by <- "variance"   # "variance" or "input"
heatmap_max_height_px <- 12000
font_row_min <- 5
min_samples_per_group <- 2

# Optional: build a "treated - vehicle" delta heatmap if a vehicle exists
infer_controls_regex <- "(?i)^(E3|DMSO)$"  # infer controls by label; edit if needed
make_delta_heatmap <- TRUE

# --- Heatmap helper: vector PDF with dynamic sizing & long-rowname padding ---
save_heatmap_pdf <- function(mat, file, ann_row=NULL, ann_col=NULL,
                             show_rownames=TRUE, main=NULL,
                             fontsize_row=NULL, fontsize_col=10,
                             cluster_rows=FALSE, cluster_cols=FALSE,
                             angle_col=0) {
  nR <- nrow(mat); nC <- ncol(mat)
  
  # Extra width if row names are long (prevents left gutter from crushing the matrix)
  rn_chars <- if (show_rownames && !is.null(rownames(mat))) {
    suppressWarnings(max(nchar(rownames(mat)), na.rm = TRUE))
  } else 0
  left_extra <- if (is.finite(rn_chars) && rn_chars > 0) pmin(8, 0.085 * rn_chars) else 0
  
  # Dynamic size in inches (bounded)
  height <- max(5,  min(80, 2 + 0.12 * nR))
  width  <- max(6 + left_extra, min(80, 3 + 0.22 * nC))
  
  if (is.null(fontsize_row)) fontsize_row <- max(5, 12 - 0.03 * nR)
  
  pheatmap::pheatmap(
    mat,
    filename       = file,       # let pheatmap manage the device
    width          = width,
    height         = height,
    scale          = "none",
    show_rownames  = show_rownames,
    annotation_row = ann_row,
    annotation_col = ann_col,
    cluster_rows   = cluster_rows,
    cluster_cols   = cluster_cols,  # keep T1..T# order
    na_col         = "grey90",
    main           = main,
    fontsize_row   = fontsize_row,
    fontsize_col   = fontsize_col,
    border_color   = NA,
    angle_col      = angle_col,
    silent         = TRUE
  )
}

# -----------------
# PACKAGES
# -----------------
cran_needed <- c("ggplot2","pheatmap","dplyr","stringr","readr","tibble","matrixStats","tidyr","purrr")
bio_needed  <- c("DESeq2","GSVA","AnnotationDbi","org.Dr.eg.db","GO.db")

ipkg <- function(pkgs) for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) install.packages(p, dependencies=TRUE)
if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
ipkg(cran_needed)
for(p in bio_needed) if(!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)

suppressPackageStartupMessages({
  library(ggplot2); library(pheatmap)
  library(dplyr);   library(stringr);  library(readr); library(tibble)
  library(tidyr);   library(purrr)
  library(matrixStats)
  library(DESeq2)
  library(GSVA)
  library(AnnotationDbi); library(org.Dr.eg.db); library(GO.db)
})

# -----------------
# LOAD DATA
# -----------------
counts <- read.table(counts_path, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta   <- read.csv(meta_path, stringsAsFactors=FALSE)

# Fix sample-name quirks
colnames(counts) <- gsub("^X","", colnames(counts))
colnames(counts) <- gsub("\\.", "-", colnames(counts))

# Keep only samples present in both
meta <- meta %>% filter(sample_id %in% colnames(counts))

# Timepoint: ensure factor T1..T10
tp_col <- intersect(timepoint_col_preferred, names(meta))[1]
if (is.na(tp_col)) {
  tp <- stringr::str_match(meta$sample_id, "-T([1-9]|10)")[,2]
  meta$timepoint_code <- paste0("T", tp)
  tp_col <- "timepoint_code"
}
meta[[tp_col]] <- factor(meta[[tp_col]], levels = paste0("T", 1:10))

# Grouping column
grp_col <- intersect(group_col_preferred, names(meta))[1]
if (is.na(grp_col)) stop("Could not find a grouping column. Add one of: ", paste(group_col_preferred, collapse=", "))

# Standardize obvious control labels
if (grp_col %in% c("vehicle","vehicle_type","label","Label")) {
  v <- meta[[grp_col]]
  v2 <- ifelse(grepl("E3", v, ignore.case=TRUE), "E3",
               ifelse(grepl("DM|DMSO", v, ignore.case=TRUE), "DMSO", v))
  meta[[grp_col]] <- v2
}

# Subset groups if requested
if (!is.null(groups_include)) meta <- meta %>% filter(.data[[grp_col]] %in% groups_include)

# Reorder counts to match meta
counts <- counts[, meta$sample_id, drop=FALSE]
stopifnot(all(colnames(counts) == meta$sample_id))

# -----------------
# HELPERS
# -----------------
is_integer_matrix <- function(m) {
  if (!is.numeric(as.matrix(m))) return(FALSE)
  all(abs(m - round(m)) < 1e-6, na.rm=TRUE)
}

# Build GO:BP gene sets keyed by ENSEMBL (zebrafish)
build_gobp_gene_sets_ensembl <- function() {
  eg2go <- suppressMessages(AnnotationDbi::select(
    x      = org.Dr.eg.db,
    keys   = AnnotationDbi::keys(org.Dr.eg.db, keytype = "ENSEMBL"),
    columns= c("ENSEMBL","GO","ONTOLOGY"),
    keytype= "ENSEMBL"
  ))
  eg2go <- eg2go %>% filter(ONTOLOGY == "BP", !is.na(GO))
  split(eg2go$ENSEMBL, eg2go$GO)
}

# Pick GO terms by keywords and/or explicit IDs; optionally add offspring
pick_go_terms <- function(all_goids, focus_keywords, focus_go_ids, include_offspring, require_godb=TRUE) {
  chosen <- character(0)
  # Start with explicit IDs
  if (length(focus_go_ids)) chosen <- union(chosen, intersect(focus_go_ids, all_goids))
  # Keyword matches (TERM contains any keyword)
  if (length(focus_keywords)) {
    if (!requireNamespace("GO.db", quietly = TRUE)) {
      if (require_godb) stop("GO.db is required to search by keyword. Install with BiocManager::install('GO.db').")
    } else {
      term_map <- AnnotationDbi::mapIds(GO.db::GO.db, keys = all_goids,
                                        column = "TERM", keytype = "GOID", multiVals = "first")
      hit <- rep(FALSE, length(term_map))
      for (kw in focus_keywords) hit <- hit | grepl(kw, term_map, ignore.case = TRUE, perl = TRUE)
      chosen <- union(chosen, names(term_map)[hit])
    }
  }
  # Expand to offspring if requested (then filter back to keywords to keep only "... signaling pathway")
  if (include_offspring && length(chosen)) {
    kids_list <- AnnotationDbi::mget(chosen, GO.db::GOBPOFFSPRING, ifnotfound = NA)
    kids <- unique(unlist(kids_list, use.names = FALSE))
    chosen <- unique(c(chosen, intersect(kids, all_goids)))
    term_map <- AnnotationDbi::mapIds(GO.db::GO.db, keys = chosen,
                                      column = "TERM", keytype = "GOID", multiVals = "first")
    keep <- rep(FALSE, length(term_map))
    for (kw in focus_keywords) keep <- keep | grepl(kw, term_map, ignore.case = TRUE, perl = TRUE)
    chosen <- names(term_map)[keep]
  }
  intersect(chosen, all_goids)
}

# Label GO IDs with term names (ensuring uniqueness)
label_go_terms <- function(goids) {
  if (!requireNamespace("GO.db", quietly=TRUE)) return(goids)
  tm <- suppressMessages(AnnotationDbi::mapIds(
    GO.db::GO.db, keys = goids, column = "TERM", keytype = "GOID", multiVals = "first"
  ))
  labels <- ifelse(is.na(tm), goids, tm)
  dup <- duplicated(labels) | duplicated(labels, fromLast=TRUE)
  labels[dup] <- paste0(labels[dup], " (", goids[dup], ")")
  labels
}

# Robust GSVA (works across GSVA versions)
run_gsva_ssgsea <- function(expr_mat, gene_sets, min_size, max_size) {
  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    par <- GSVA::ssgseaParam(as.matrix(expr_mat), gene_sets, minSize = min_size, maxSize = max_size)
    GSVA::gsva(par, verbose = FALSE)
  } else {
    GSVA::gsva(expr=as.matrix(expr_mat), gset.idx.list=gene_sets,
               method="ssGSEA", min.sz=min_size, max.sz=max_size, ssgsea.norm=TRUE)
  }
}

# -----------------
# GLOBAL: build gene sets + select ALL "… signaling pathway" GO terms
# -----------------
message("Building GO:BP gene sets (org.Dr.eg.db)…")
gobp_sets_all <- build_gobp_gene_sets_ensembl()
all_goids <- names(gobp_sets_all)

message("Selecting GO terms that match: ", paste(focus_keywords, collapse=" | "))
chosen_goids_global <- pick_go_terms(
  all_goids        = all_goids,
  focus_keywords   = focus_keywords,
  focus_go_ids     = focus_go_ids,
  include_offspring= include_offspring,
  require_godb     = require_GOdb
)
if (!length(chosen_goids_global)) stop("No signaling-pathway GO terms found. Check GO.db installation.")

# Order by TERM alphabetically initially; we'll re-order by variance later
chosen_goids_global <- sort(chosen_goids_global)
chosen_labels_global <- label_go_terms(chosen_goids_global)

# -----------------
# MAIN: iterate groups, compute ssGSEA, collect timepoint means
# -----------------
all_groups <- unique(meta[[grp_col]])
tp_means_by_group <- list()
present_rows_union <- character(0)  # track union of pathways that pass size filter per group
sample_scores_long <- list()

for (g in all_groups) {
  message("\n==============================")
  message("Group: ", g)
  message("==============================")
  
  gdir <- file.path(out_dir, make.names(grp_col), make.names(g))
  dir.create(gdir, showWarnings = FALSE, recursive = TRUE)
  
  md  <- meta %>% filter(.data[[grp_col]] == g)
  if (nrow(md) < min_samples_per_group) {
    message("Skipping '", g, "' (only ", nrow(md), " samples; need >= ", min_samples_per_group, ").")
    next
  }
  cnt <- counts[, md$sample_id, drop=FALSE]
  
  # Filter low-count genes: expressed in >=20% of samples (min 2)
  keep <- rowSums(cnt > 0) >= max(2, ceiling(0.2*ncol(cnt)))
  cnt  <- cnt[keep, , drop=FALSE]
  if (nrow(cnt) == 0) { message("No genes after filtering. Skipping ", g); next }
  
  # Normalize for GSVA
  if (is_integer_matrix(cnt)) {
    dds <- DESeqDataSetFromMatrix(countData = round(cnt), colData = md, design = ~ 1)
    dds <- suppressMessages(estimateSizeFactors(dds))
    vsd <- suppressMessages(vst(dds, blind = TRUE))
    expr_mat <- assay(vsd)
  } else {
    expr_mat <- log2(cnt + 1)
  }
  rownames(expr_mat) <- sub("\\..*$", "", rownames(expr_mat))  # strip Ensembl versions
  
  # Intersect gene sets with available genes
  gsets_focus <- gobp_sets_all[chosen_goids_global]
  gsets_focus <- lapply(gsets_focus, function(s) intersect(s, rownames(expr_mat)))
  
  # Size filter per group
  keep_set <- vapply(gsets_focus, function(x) {
    lx <- length(x); !is.na(lx) && lx >= gsva_min_size && lx <= gsva_max_size
  }, logical(1))
  gsets_focus <- gsets_focus[keep_set]
  
  if (!length(gsets_focus)) { message("No pathway sets passed size filter for ", g); next }
  
  # Save chosen sets for this group
  sizes <- vapply(gsets_focus, length, integer(1))
  term_names <- label_go_terms(names(gsets_focus))
  write.csv(
    data.frame(GOID = names(gsets_focus), TERM = term_names, SIZE = sizes),
    file.path(gdir, "Focused_GO_sets.csv"), row.names = FALSE
  )
  
  present_rows_union <- union(present_rows_union, names(gsets_focus))
  
  # Run GSVA
  message("Running GSVA (ssGSEA) on ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples...")
  gsva_scores <- run_gsva_ssgsea(expr_mat, gsets_focus, gsva_min_size, gsva_max_size)
  
  # Save sample-level scores
  write.csv(gsva_scores, file.path(gdir, "GSVA_scores_per_sample.csv"))
  
  # Build a mapping table from meta and give it deterministic names
  md_map <- md %>%
    dplyr::select(sample_id, dplyr::all_of(c(tp_col, grp_col))) %>%
    as.data.frame(stringsAsFactors = FALSE)
  names(md_map)[names(md_map) == tp_col]  <- "Timepoint"
  names(md_map)[names(md_map) == grp_col] <- "Group"
  
  # Long table (sample-level) for downstream plots
  sample_scores_long[[g]] <- gsva_scores %>%
    as.data.frame() %>%
    tibble::rownames_to_column("GOID") %>%
    dplyr::mutate(TERM = label_go_terms(GOID)) %>%
    tidyr::pivot_longer(-c(GOID, TERM), names_to = "sample_id", values_to = "ssGSEA") %>%
    dplyr::left_join(md_map, by = "sample_id")
  
  # Summarize by timepoint (replicate means)
  tp_levels <- levels(md[[tp_col]])
  tp_levels <- tp_levels[tp_levels %in% md[[tp_col]]]
  tp_means <- sapply(tp_levels, function(tp) {
    cols <- md$sample_id[md[[tp_col]] == tp]
    if (length(cols) == 1) gsva_scores[, cols, drop = FALSE] else rowMeans(gsva_scores[, cols, drop = FALSE])
  })
  if (is.null(dim(tp_means))) tp_means <- matrix(tp_means, ncol = length(tp_levels))
  rownames(tp_means) <- rownames(gsva_scores); colnames(tp_means) <- tp_levels
  
  # Reorder pathways (variance across timepoints)
  if (order_pathways_by == "variance") {
    v_path <- matrixStats::rowVars(as.matrix(tp_means), na.rm = TRUE)
    tp_means <- tp_means[order(v_path, decreasing = TRUE), , drop = FALSE]
  } else {
    tp_means <- tp_means[names(gsets_focus), , drop = FALSE]
  }
  
  # Label rows with TERM
  rownames(tp_means) <- label_go_terms(rownames(tp_means))
  
  # Save timepoint means + per-group heatmap (PDF, vector)
  write.csv(tp_means, file.path(gdir, "GSVA_scores_timepoint_means.csv"))
  
  n_paths <- nrow(tp_means)
  fs_row  <- max(font_row_min, 12 - 0.03 * n_paths)
  
  ann_col <- data.frame(Timepoint = factor(colnames(tp_means), levels = colnames(tp_means)))
  row.names(ann_col) <- colnames(tp_means)
  
  save_heatmap_pdf(
    mat = tp_means,
    file = file.path(gdir, "GSVA_focus_heatmap.pdf"),
    ann_col = ann_col,
    main = paste0("Signaling pathways - ", grp_col, ": ", g, " (", n_paths, " pathways)"),
    show_rownames = TRUE,
    fontsize_row = fs_row,
    fontsize_col = 10,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    angle_col = 0
  )
  
  tp_means_by_group[[g]] <- tp_means
  message("Done group: ", g, " -> ", gdir)
}

# -----------------
# COMBINE across groups -> one grand heatmap
# -----------------
if (length(tp_means_by_group)) {
  # master row order = union that passed size filter anywhere, labeled as TERM
  master_rows <- label_go_terms(sort(present_rows_union))
  
  # build combined matrix with columns "Group|T#"; fill missing with NA
  combo_list <- lapply(names(tp_means_by_group), function(g) {
    x <- tp_means_by_group[[g]]
    x2 <- matrix(NA_real_, nrow=length(master_rows), ncol=ncol(x),
                 dimnames = list(master_rows, colnames(x)))
    common <- intersect(rownames(x), master_rows)
    x2[common, ] <- x[common, , drop=FALSE]
    colnames(x2) <- paste0(g, "|", colnames(x2))
    x2
  })
  combo <- do.call(cbind, combo_list)
  
  # reorder rows by variance across ALL columns
  v_all <- matrixStats::rowVars(as.matrix(combo), na.rm = TRUE)
  ord <- order(v_all, decreasing = TRUE, na.last = TRUE)
  combo <- combo[ord, , drop = FALSE]
  
  write.csv(combo, file.path(out_dir, "GSVA_signaling_pathways_ALL_groups_timepoints.csv"))
  
  # build column annotations by parsing "Group|T#"
  parts <- tidyr::separate(
    data.frame(col = colnames(combo)),
    col, into = c("Group","Timepoint"), sep = "\\|", remove = TRUE, fill = "right"
  )
  tp_order_all <- paste0("T", 1:10)
  ann_col_all <- data.frame(
    Group     = factor(parts$Group),
    Timepoint = factor(parts$Timepoint, levels = tp_order_all[tp_order_all %in% parts$Timepoint])
  )
  row.names(ann_col_all) <- colnames(combo)
  
  n_paths_all <- nrow(combo)
  fs_row_all  <- max(font_row_min, 12 - 0.03 * n_paths_all)
  
  save_heatmap_pdf(
    mat = combo,
    file = file.path(out_dir, "GSVA_signaling_pathways_COMBINED.pdf"),
    ann_col = ann_col_all,
    main = paste0("Signaling pathways – combined (", n_paths_all, " pathways; columns = Group|Timepoint)"),
    show_rownames = TRUE,
    fontsize_row = fs_row_all,
    fontsize_col = 8,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    angle_col = 0
  )
}

# -----------------
# OPTIONAL: Delta (treated - vehicle) heatmap
# -----------------
if (make_delta_heatmap && length(tp_means_by_group)) {
  controls <- names(tp_means_by_group)[grepl(infer_controls_regex, names(tp_means_by_group))]
  if (length(controls)) {
    # average control per timepoint across all control groups present
    ctrl_mats <- tp_means_by_group[controls]
    common_rows <- Reduce(intersect, lapply(ctrl_mats, rownames))
    tps_all <- paste0("T", 1:10)
    tps_all <- tps_all[tps_all %in% unique(unlist(lapply(ctrl_mats, colnames)))]
    ctrl_avg <- matrix(NA_real_, nrow=length(common_rows), ncol=length(tps_all),
                       dimnames=list(common_rows, tps_all))
    for (tp in tps_all) {
      cols <- lapply(ctrl_mats, function(m) if (tp %in% colnames(m)) m[common_rows, tp, drop=FALSE] else NULL)
      cols <- cols[!sapply(cols, is.null)]
      if (length(cols)) ctrl_avg[, tp] <- Reduce("+", cols) / length(cols)
    }
    
    # deltas for each non-control group
    treated_groups <- setdiff(names(tp_means_by_group), controls)
    if (length(treated_groups)) {
      delta_list <- lapply(treated_groups, function(g) {
        M <- tp_means_by_group[[g]]
        rows <- intersect(rownames(M), rownames(ctrl_avg))
        tps  <- intersect(colnames(M), colnames(ctrl_avg))
        if (!length(rows) || !length(tps)) return(NULL)
        D <- M[rows, tps, drop=FALSE] - ctrl_avg[rows, tps, drop=FALSE]
        colnames(D) <- paste0(g, "|ΔvsCtrl|", colnames(D))
        D
      })
      delta_list <- delta_list[!sapply(delta_list, is.null)]
      if (length(delta_list)) {
        delta_combo <- do.call(cbind, delta_list)
        v_delta <- matrixStats::rowVars(as.matrix(delta_combo), na.rm = TRUE)
        delta_combo <- delta_combo[order(v_delta, decreasing=TRUE), , drop=FALSE]
        write.csv(delta_combo, file.path(out_dir, "GSVA_signaling_pathways_DELTA_vs_CONTROL.csv"))
        
        # Build annotations "Group|ΔvsCtrl|T#"
        parts_d <- tidyr::separate(
          data.frame(col = colnames(delta_combo)),
          col, into = c("Group","Delta","Timepoint"), sep = "\\|", fill = "right"
        )
        tp_order_all <- paste0("T", 1:10)
        ann_col_delta <- data.frame(
          Group     = factor(parts_d$Group),
          Timepoint = factor(parts_d$Timepoint, levels = tp_order_all[tp_order_all %in% parts_d$Timepoint])
        )
        row.names(ann_col_delta) <- colnames(delta_combo)
        
        fs_row_delta <- max(font_row_min, 12 - 0.03 * nrow(delta_combo))
        
        save_heatmap_pdf(
          mat = delta_combo,
          file = file.path(out_dir, "GSVA_signaling_pathways_DELTA_vs_CONTROL.pdf"),
          ann_col = ann_col_delta,
          main = "Δ ssGSEA (treated - control) – signaling pathways",
          show_rownames = TRUE,
          fontsize_row = fs_row_delta,
          fontsize_col = 8,
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          angle_col = 0
        )
      }
    }
  }
}

# -----------------
# EXPORT: tidy long CSV across all groups (sample-level)
# -----------------
if (length(sample_scores_long)) {
  long_all <- dplyr::bind_rows(sample_scores_long)
  write.csv(long_all, file.path(out_dir, "GSVA_signaling_pathways_scores_long.csv"), row.names = FALSE)
}

message("\nAll done. Outputs in: ", out_dir)
