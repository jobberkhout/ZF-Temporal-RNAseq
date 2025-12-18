# ============================================================
# Focused GSVA (Conservative NTC)
#   • Robust group parsing for CHEM|S# (incl. '2-4D')
#   • Safe folder names (no pipes, whitespace trimmed, length-capped)
#   • Per-group + combined heatmaps + Δ vs matched vehicle + paired Δ
#   • Timepoint heatmaps per sample
#   • NTC significance (CONSERVATIVE): joint Chem×TP scoring + limma eBayes + broader FDR + effect-size floor
#   • Expanded NTC GO set builder: seeds + offspring + keyword hits (+ guarded pathway hits)
# ============================================================

# -----------------
# USER SETTINGS
# -----------------
counts_path <- "C:\\Users\\jobbe\\RNAseq_analysis\\ZF_RNA_Raw_Count.txt"
meta_path   <- "C:\\Users\\jobbe\\RNAseq_analysis\\metadata\\20250520_ZF_RNAseq_metadata_JB.csv"
out_dir     <- "C:\\Users\\jobbe\\RNAseq_analysis\\Output\\GSVA_NTC"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Optionally limit to specific groups, e.g. c("E3","DMSO") or c("CAFF|S1","CYCL|S2")
groups_include <- NULL

# Timepoint column preferences (used only if present)
timepoint_col_preferred <- c("timepoint_code","timepoint","tp","TP","time")

# Δ output toggles
make_delta_heatmap_avg    <- TRUE   # matched-vehicle Δ (treated − avg E3/DMSO or its vehicle)
make_delta_heatmap_paired <- TRUE   # paired Δ (treated − its matched control sample ID)

# GSVA gene set size limits (slightly stricter)
gsva_min_size <- 5
gsva_max_size <- 20000

# Heatmap settings
order_pathways_by <- "variance"
font_row_min <- 5
min_samples_per_group <- 2

# Identify vehicle groups by name
infer_controls_regex <- "(?i)^(E3|DMSO)$"

# ----- Conservative NTC testing toggles -----
ntc_use_limma           <- TRUE      # use limma (eBayes)
ntc_joint_scoring       <- TRUE      # score treated + vehicle jointly within Chem×TP
ntc_effect_size_floor   <- 0.02      # require |Δ ssGSEA| >= this
ntc_min_samples_per_arm <- 2         # need >=2 treated and >=2 controls at a timepoint
ntc_fdr_family          <- "by_group_by_time" #"by_chemical"
# options: "by_group_by_time", "by_group_across_time", "by_chemical"
ntc_pool_S1_S2          <- TRUE   # <— NEW: also test pooled S1+S2 vs Ctrl

# -----------------
# NTC-focused GO selection (keywords + seeds)
# -----------------
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
# -----------------
# PACKAGES
# -----------------
cran_needed <- c("ggplot2","pheatmap","dplyr","stringr","readr","tibble","matrixStats","tidyr","purrr")
bio_needed  <- c("DESeq2","GSVA","AnnotationDbi","org.Dr.eg.db","GO.db","limma")
ipkg <- function(pkgs) for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) install.packages(p, dependencies=TRUE)
if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
ipkg(cran_needed)
for(p in bio_needed) if(!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)
suppressPackageStartupMessages({
  library(ggplot2); library(pheatmap)
  library(dplyr);   library(stringr);  library(readr); library(tibble)
  library(tidyr);   library(purrr);    library(matrixStats)
  library(DESeq2);  library(GSVA);     library(limma)
  library(AnnotationDbi); library(org.Dr.eg.db); library(GO.db)
})

# -----------------
# LOAD DATA
# -----------------
counts <- read.table(counts_path, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta   <- read.csv(meta_path, stringsAsFactors=FALSE)

# Fix sample-name quirks + trim whitespace
colnames(counts) <- gsub("^X","", colnames(counts))
colnames(counts) <- gsub("\\.", "-", colnames(counts))
colnames(counts) <- trimws(colnames(counts))

# If metadata had no headers (V1, V2, ...), assign positional names per your example
if (all(grepl("^V\\d+$", names(meta)))) {
  pos_names <- c("sample_id","idx","chem4","chem_name","dose_code","dose_value","timepoint","stage",
                 "date1","blank1","date2","replicate","vehicle","matched_control_id")
  names(meta)[seq_along(pos_names)] <- pos_names
}

# Basic trims/normalization on key fields (handles stray spaces)
if (!"sample_id" %in% names(meta)) stop("Metadata must contain 'sample_id' or have it as the first column.")
meta$sample_id <- trimws(meta$sample_id)

# Keep only samples present in both
meta <- meta %>% filter(sample_id %in% colnames(counts))

# -----------------
# MAP COLUMNS (flexible) + build Group (CHEM|S# for treated; E3/DMSO for controls)
# -----------------
pick_first_col <- function(cands, df) { x <- intersect(cands, names(df)); if (length(x)) x[1] else NA_character_ }

compound_col_preferred <- c("compound_code","chem4","Chem4","chemical_code","code","compound","ChemCode")
dose_code_col_preferred<- c("dose_code","DoseCode","dose","Dose","dose_tier")
vehicle_col_preferred  <- c("vehicle_type","vehicle","Vehicle","veh")
matched_ctrl_col_pref  <- c("matched_control_id","match_control_id","control_pair_id","ctrl_pair_id",
                            "matched_ctrl","control_sample_id")

compound_col     <- pick_first_col(compound_col_preferred, meta)
dose_code_col    <- pick_first_col(dose_code_col_preferred, meta)
vehicle_col      <- pick_first_col(vehicle_col_preferred, meta)
timepoint_col    <- pick_first_col(timepoint_col_preferred, meta)
matched_ctrl_col <- pick_first_col(matched_ctrl_col_pref, meta)

# Heuristic: locate matched control column if not detected
if (is.na(matched_ctrl_col)) {
  cand <- names(meta)[sapply(meta, function(x) is.character(x) && any(grepl("CTRL", x, ignore.case=TRUE)))]
  cand <- setdiff(cand, c("sample_id", vehicle_col))
  if (length(cand)) matched_ctrl_col <- cand[1]
}

# Normalize case + trim whitespace on key columns
if (!is.na(compound_col)) meta[[compound_col]] <- toupper(trimws(meta[[compound_col]]))
if (!is.na(dose_code_col)) meta[[dose_code_col]] <- toupper(trimws(meta[[dose_code_col]]))
if (!is.na(vehicle_col))   meta[[vehicle_col]]   <- toupper(trimws(meta[[vehicle_col]]))
if (!is.na(matched_ctrl_col)) meta[[matched_ctrl_col]] <- trimws(meta[[matched_ctrl_col]])

# Normalize vehicle values to E3/DMSO
if (!is.na(vehicle_col)) {
  v <- meta[[vehicle_col]]
  meta[[vehicle_col]] <- ifelse(grepl("^E3$", v), "E3",
                                ifelse(grepl("^DMSO$|^DM$", v), "DMSO", v))
}

# Timepoint handling (dynamic T-levels)
if (!is.na(timepoint_col)) {
  meta[[timepoint_col]] <- toupper(trimws(meta[[timepoint_col]]))
  tp_nums <- suppressWarnings(as.integer(sub("^T","", meta[[timepoint_col]])))
  tp_levels <- paste0("T", sort(unique(na.omit(tp_nums))))
  meta[[timepoint_col]] <- factor(meta[[timepoint_col]], levels = tp_levels)
} else {
  tp <- stringr::str_match(meta$sample_id, "-(T[0-9]{1,2})[-_]")[,2]
  tp_nums <- suppressWarnings(as.integer(sub("^T","", tp)))
  tp_levels <- paste0("T", sort(unique(na.omit(tp_nums))))
  meta$timepoint_code <- factor(tp, levels = tp_levels)
  timepoint_col <- "timepoint_code"
}

# Build Group:
is_control_row <- !is.na(dose_code_col) & grepl("^C$", meta[[dose_code_col]], ignore.case=TRUE)
meta$Group <- NA_character_
if (!is.na(vehicle_col)) meta$Group[is_control_row] <- meta[[vehicle_col]][is_control_row]

if (!is.na(compound_col) & !is.na(dose_code_col)) {
  meta$Group[!is_control_row] <- paste0(meta[[compound_col]][!is_control_row], "|", meta[[dose_code_col]][!is_control_row])
} else {
  parsed <- stringr::str_match(meta$sample_id, "-T[0-9]{1,2}-([A-Z0-9.-]{2,15})-(S[1-3]|C)\\b")
  ok <- is.na(meta$Group) & !is.na(parsed[,2])
  meta$Group[ok & parsed[,3] != "C"] <- paste0(parsed[ok & parsed[,3]!="C",2], "|", parsed[ok & parsed[,3]!="C",3])
  if (!is.na(vehicle_col)) meta$Group[is.na(meta$Group)] <- meta[[vehicle_col]][is.na(meta$Group)]
}

meta$Group <- trimws(meta$Group)

# Dose order S1<S2<S3 (C for bookkeeping)
meta$DoseTier <- NA_character_
if (!is.na(dose_code_col)) meta$DoseTier <- meta[[dose_code_col]]
meta$DoseTier <- factor(meta$DoseTier, levels = c("S1","S2","S3","C"))

# Subset groups if requested
if (!is.null(groups_include)) meta <- meta %>% filter(Group %in% groups_include)

# Reorder counts to match meta
counts <- counts[, meta$sample_id, drop=FALSE]
stopifnot(all(colnames(counts) == meta$sample_id))

# -----------------
# Helpers
# -----------------
is_integer_matrix <- function(m) { if (!is.numeric(as.matrix(m))) return(FALSE); all(abs(m - round(m)) < 1e-6, na.rm=TRUE) }

safe_group_dir <- function(g) {
  g <- trimws(g)
  g <- gsub("\\s+", "", g)
  g <- gsub("\\|", "_", g)
  g <- gsub("[^A-Za-z0-9._-]+", "_", g)
  g <- sub("^\\.+", "", g)
  if (nchar(g) > 60) {
    cs <- sprintf("%06d", sum(utf8ToInt(g)) %% 1e6)
    g <- paste0(substr(g, 1, 50), "_", cs)
  }
  if (g == "") g <- "Group"
  g
}

save_heatmap_pdf <- function(mat, file, ann_row=NULL, ann_col=NULL,
                             show_rownames=TRUE, main=NULL,
                             fontsize_row=NULL, fontsize_col=10,
                             cluster_rows=FALSE, cluster_cols=FALSE,
                             angle_col=90, ...) {
  nR <- nrow(mat); nC <- ncol(mat)
  rn_chars <- if (show_rownames && !is.null(rownames(mat))) suppressWarnings(max(nchar(rownames(mat)), na.rm = TRUE)) else 0
  left_extra <- if (is.finite(rn_chars) && rn_chars > 0) pmin(8, 0.085 * rn_chars) else 0
  height <- max(5,  min(80, 2 + 0.12 * nR))
  width  <- max(6 + left_extra, min(80, 3 + 0.22 * nC))
  if (is.null(fontsize_row)) fontsize_row <- max(5, 12 - 0.03 * nR)
  pheatmap::pheatmap(
    mat, filename=file, width=width, height=height, scale="none",
    show_rownames=show_rownames, annotation_row=ann_row, annotation_col=ann_col,
    cluster_rows=cluster_rows, cluster_cols=cluster_cols, na_col="grey90",
    main=main, fontsize_row=fontsize_row, fontsize_col=fontsize_col,
    border_color=NA, angle_col=angle_col, silent=TRUE, ...
  )
}

# ssGSEA wrapper (supports both GSVA APIs)
run_gsva_ssgsea <- function(expr_mat, gene_sets, min_size, max_size) {
  res <- if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    par <- GSVA::ssgseaParam(as.matrix(expr_mat), gene_sets, minSize = min_size, maxSize = max_size)
    GSVA::gsva(par, verbose = FALSE)
  } else {
    GSVA::gsva(expr=as.matrix(expr_mat), gset.idx.list=gene_sets, method="ssgsea",
               min.sz=min_size, max.sz=max_size, ssgsea.norm=TRUE)
  }
  if (inherits(res, "SummarizedExperiment")) {
    res <- SummarizedExperiment::assay(res)
  }
  as.matrix(res)
}

parse_group_time <- function(cols) {
  data.frame(
    Group     = sub("\\|[^|]*$", "", cols),
    Timepoint = sub("^.*\\|", "", cols),
    stringsAsFactors = FALSE
  )
}
parse_group_delta_time <- function(cols) {
  data.frame(
    Group     = sub("\\|[^|]+\\|[^|]+$", "", cols),
    Delta     = sub("^.*\\|([^|]+)\\|[^|]+$", "\\1", cols),
    Timepoint = sub("^.*\\|", "", cols),
    stringsAsFactors = FALSE
  )
}

choose_ctrl_matrix <- function(g_name, ctrl_mats) {
  vtype <- group_vehicle_map[[g_name]]
  if (!is.null(vtype) && !is.na(vtype) && vtype %in% names(ctrl_mats)) return(ctrl_mats[[vtype]])
  if (length(ctrl_mats) > 1) {
    rows <- Reduce(intersect, lapply(ctrl_mats, rownames))
    tps  <- Reduce(intersect, lapply(ctrl_mats, colnames))
    avg  <- Reduce(`+`, lapply(ctrl_mats, function(m) m[rows, tps, drop=FALSE])) / length(ctrl_mats)
    rownames(avg) <- rows; colnames(avg) <- tps
    return(avg)
  }
  ctrl_mats[[1]]
}

label_go_terms <- function(goids) {
  tm <- suppressMessages(AnnotationDbi::mapIds(GO.db::GO.db, keys = goids, column = "TERM",
                                               keytype = "GOID", multiVals = "first"))
  labels <- ifelse(is.na(tm), goids, tm)
  dup <- duplicated(labels) | duplicated(labels, fromLast=TRUE)
  labels[dup] <- paste0(labels[dup], " (", goids[dup], ")")
  labels
}

# -----------------
# GLOBAL: build gene sets + expanded NTC GO selection (no obsolete filtering)
# -----------------
message("Building GO:BP gene sets (org.Dr.eg.db)…")
build_gobp_gene_sets_ensembl <- function() {
  eg2go <- suppressMessages(AnnotationDbi::select(
    x      = org.Dr.eg.db,
    keys   = AnnotationDbi::keys(org.Dr.eg.db, keytype = "ENSEMBL"),
    columns= c("ENSEMBL","GO","ONTOLOGY"),
    keytype= "ENSEMBL"
  ))
  eg2go <- eg2go %>% dplyr::filter(ONTOLOGY == "BP", !is.na(GO))
  split(eg2go$ENSEMBL, eg2go$GO)
}
gobp_sets_all <- build_gobp_gene_sets_ensembl()

collapse_unique <- function(x, sep = " | ") {
  x <- unique(unlist(x, use.names = FALSE))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  paste(x, collapse = sep)
}

expand_ntc_go <- function(seed_ids,
                          keyword_patterns,
                          pathway_patterns = NULL,
                          restrict_ontology = c("BP"),
                          drop_obsolete = TRUE) {  # kept for API compatibility; ignored
  all_ids <- AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID")
  
  meta_basic <- suppressMessages(
    AnnotationDbi::select(GO.db::GO.db, keys = all_ids,
                          keytype = "GOID",
                          columns = c("GOID","TERM","ONTOLOGY"))
  ) %>% dplyr::distinct(GOID, .keep_all = TRUE)
  
  # Robust synonym retrieval across GO.db versions
  avail_cols <- tryCatch(AnnotationDbi::columns(GO.db::GO.db), error = function(e) character(0))
  if ("SYNONYM" %in% avail_cols) {
    syn_tbl <- suppressMessages(
      AnnotationDbi::select(GO.db::GO.db, keys = all_ids,
                            keytype = "GOID", columns = "SYNONYM")
    ) %>%
      dplyr::group_by(GOID) %>%
      dplyr::summarize(SYNONYMS = collapse_unique(SYNONYM), .groups = "drop")
  } else if (exists("GOSYNONYM", envir = asNamespace("GO.db"), inherits = FALSE)) {
    syn_map <- AnnotationDbi::mget(all_ids, GO.db::GOSYNONYM, ifnotfound = NA)
    syn_tbl <- tibble::tibble(
      GOID = names(syn_map),
      SYNONYMS = vapply(syn_map, collapse_unique, character(1))
    )
  } else {
    syn_tbl <- tibble::tibble(GOID = all_ids, SYNONYMS = NA_character_)
  }
  
  meta <- meta_basic %>% dplyr::left_join(syn_tbl, by = "GOID")
  if (!is.null(restrict_ontology)) {
    meta <- meta %>% dplyr::filter(ONTOLOGY %in% restrict_ontology)
  }
  
  # Build regexes
  rx_any  <- function(patts) paste0("(", paste(patts, collapse = "|"), ")")
  kw_rx   <- stringr::regex(rx_any(keyword_patterns), ignore_case = TRUE)
  path_rx <- if (!is.null(pathway_patterns))
    stringr::regex(rx_any(pathway_patterns), ignore_case = TRUE) else NULL
  # Guarded “pathway” hits must also mention a morphogen axis term
  axis_rx <- stringr::regex("(retinoic acid|retinoid|RALDH|CYP26|RDH10|DHRS3|RAR|RXR|sonic hedgehog|shh|hedgehog|bmp|wnt|fgf|notch|tgf\\s*β|tgf\\s*b|smad|erbb|egf|vegf|pdgf|eph|ephrin|semaphorin|netrin|slit|robo)", ignore_case = TRUE)
  
  text_blob <- function(term, syn) paste(term, ifelse(is.na(syn), "", syn), sep = " | ")
  meta$TEXT <- text_blob(meta$TERM, meta$SYNONYMS)
  
  # Keyword hits
  kw_hits <- meta %>% dplyr::filter(stringr::str_detect(TEXT, kw_rx)) %>%
    dplyr::pull(GOID) %>% unique()
  
  # Broad pathway hits (guarded by axis keywords)
  broad_hits <- character(0)
  if (!is.null(path_rx)) {
    broad_hits <- meta %>%
      dplyr::filter(stringr::str_detect(TERM, path_rx) & stringr::str_detect(TERM, axis_rx)) %>%
      dplyr::pull(GOID) %>% unique()
  }
  
  # Seeds + BP offspring
  .safe_mget <- function(keys, env) unique(unlist(mget(keys, envir = env, ifnotfound = NA), use.names = FALSE))
  seed_bp <- intersect(seed_ids, meta$GOID)
  bp_offs <- intersect(.safe_mget(seed_bp, GO.db::GOBPOFFSPRING), meta$GOID)
  
  union_ids <- unique(c(seed_bp, bp_offs, kw_hits, broad_hits))
  
  prov <- dplyr::bind_rows(
    tibble::tibble(GOID = seed_bp,    source = "seed"),
    tibble::tibble(GOID = bp_offs,    source = "offspring"),
    tibble::tibble(GOID = kw_hits,    source = "keyword"),
    tibble::tibble(GOID = broad_hits, source = "broad_pathway")
  ) %>% dplyr::distinct() %>%
    dplyr::group_by(GOID) %>%
    dplyr::summarize(Source = paste(sort(unique(source)), collapse = ","), .groups = "drop")
  
  out_tbl <- meta %>% dplyr::filter(GOID %in% union_ids) %>%
    dplyr::left_join(prov, by = "GOID") %>%
    dplyr::select(GOID, TERM, ONTOLOGY, Source, SYNONYMS)
  
  list(go_ids = union_ids, table = out_tbl)
}

message("Selecting GO terms via expanded NTC builder…")
ntc_expanded <- expand_ntc_go(
  seed_ids          = unique(ntc_goids_seed),
  keyword_patterns  = ntc_keywords,
  pathway_patterns  = broad_pathway_terms,
  restrict_ontology = "BP",
  drop_obsolete     = TRUE  # ignored
)

# This is the full (broader) set we’ll score everywhere:
chosen_goids_global  <- sort(ntc_expanded$go_ids)
chosen_labels_global <- label_go_terms(chosen_goids_global)

# Save provenance table for transparency
write.csv(ntc_expanded$table %>% dplyr::arrange(TERM),
          file.path(out_dir, "NTC_GO_expanded_BP.csv"), row.names = FALSE)

# For NTC-specific testing we reuse the same broadened set:
ntc_goids_only <- chosen_goids_global
message("Expanded NTC GO count: ", length(ntc_goids_only))


# -----------------
# MAIN: iterate groups, compute ssGSEA, collect timepoint means
# -----------------
all_groups <- unique(meta$Group)
tp_means_by_group <- list()
present_rows_union <- character(0)
sample_scores_long <- list()
group_vehicle_map <- list()

for (g in all_groups) {
  message("\n=============================="); message("Group: ", g); message("==============================")
  gdir <- file.path(out_dir, "Group", safe_group_dir(g))
  dir.create(gdir, showWarnings = FALSE, recursive = TRUE)
  
  md <- meta %>% filter(Group == g)
  if (nrow(md) < min_samples_per_group) { message("Skipping '", g, "' (only ", nrow(md), " samples; need >= ", min_samples_per_group, ")."); next }
  
  veh_in_group <- if (!is.na(vehicle_col)) unique(na.omit(toupper(md[[vehicle_col]]))) else character(0)
  group_vehicle_map[[g]] <- if (length(veh_in_group)) veh_in_group[1] else NA_character_
  
  cnt <- counts[, md$sample_id, drop=FALSE]
  keep <- rowSums(cnt > 0) >= max(2, ceiling(0.2*ncol(cnt)))
  cnt  <- cnt[keep, , drop=FALSE]
  if (nrow(cnt) == 0) { message("No genes after filtering. Skipping ", g); next }
  
  if (is_integer_matrix(cnt)) {
    dds <- DESeqDataSetFromMatrix(countData = round(cnt), colData = md, design = ~ 1)
    dds <- suppressMessages(estimateSizeFactors(dds))
    vsd <- suppressMessages(vst(dds, blind = TRUE))
    expr_mat <- assay(vsd)
  } else {
    expr_mat <- log2(cnt + 1)
  }
  rownames(expr_mat) <- sub("\\..*$", "", rownames(expr_mat))
  
  gsets_focus <- gobp_sets_all[chosen_goids_global]
  gsets_focus <- lapply(gsets_focus, function(s) intersect(s, rownames(expr_mat)))
  keep_set <- vapply(gsets_focus, function(x) { lx <- length(x); !is.na(lx) && lx >= gsva_min_size && lx <= gsva_max_size }, logical(1))
  gsets_focus <- gsets_focus[keep_set]
  if (!length(gsets_focus)) { message("No pathway sets passed size filter for ", g); next }
  
  sizes <- vapply(gsets_focus, length, integer(1))
  term_names <- label_go_terms(names(gsets_focus))
  write.csv(data.frame(GOID = names(gsets_focus), TERM = term_names, SIZE = sizes),
            file.path(gdir, "Focused_GO_sets.csv"), row.names = FALSE)
  
  present_rows_union <- union(present_rows_union, names(gsets_focus))
  
  message("Running GSVA (ssGSEA) on ", nrow(expr_mat), " genes x ", ncol(expr_mat), " samples...")
  gsva_scores <- run_gsva_ssgsea(expr_mat, gsets_focus, gsva_min_size, gsva_max_size)
  
  write.csv(gsva_scores, file.path(gdir, "GSVA_scores_per_sample.csv"))
  
  md_map <- md %>% dplyr::select(sample_id, dplyr::all_of(c(timepoint_col, "Group"))) %>% as.data.frame(stringsAsFactors = FALSE)
  names(md_map)[names(md_map) == timepoint_col] <- "Timepoint"
  
  sample_scores_long[[g]] <- gsva_scores %>%
    as.data.frame() %>% tibble::rownames_to_column("GOID") %>%
    dplyr::mutate(TERM = label_go_terms(GOID)) %>%
    tidyr::pivot_longer(-c(GOID, TERM), names_to = "sample_id", values_to = "ssGSEA") %>%
    dplyr::left_join(md_map, by = "sample_id")
  
  tp_levels_g <- levels(meta[[timepoint_col]])
  tp_levels_g <- tp_levels_g[tp_levels_g %in% md[[timepoint_col]]]
  tp_means <- sapply(tp_levels_g, function(tp) {
    cols <- md$sample_id[md[[timepoint_col]] == tp]
    if (length(cols) == 1) gsva_scores[, cols, drop = FALSE] else rowMeans(gsva_scores[, cols, drop = FALSE])
  })
  if (is.null(dim(tp_means))) tp_means <- matrix(tp_means, ncol = length(tp_levels_g))
  rownames(tp_means) <- rownames(gsva_scores); colnames(tp_means) <- tp_levels_g
  
  if (order_pathways_by == "variance") {
    v_path <- matrixStats::rowVars(as.matrix(tp_means), na.rm = TRUE)
    tp_means <- tp_means[order(v_path, decreasing = TRUE), , drop = FALSE]
  } else {
    tp_means <- tp_means[names(gsets_focus), , drop = FALSE]
  }
  rownames(tp_means) <- label_go_terms(rownames(tp_means))
  
  write.csv(tp_means, file.path(gdir, "GSVA_scores_timepoint_means.csv"))
  n_paths <- nrow(tp_means); fs_row <- max(font_row_min, 12 - 0.03 * n_paths)
  ann_col <- data.frame(Timepoint = factor(colnames(tp_means), levels = colnames(tp_means)))
  row.names(ann_col) <- colnames(tp_means)
  
  save_heatmap_pdf(
    mat = tp_means,
    file = file.path(gdir, "GSVA_focus_heatmap.pdf"),
    ann_col = ann_col,
    main = paste0("Signaling pathways - Group: ", g, " (", n_paths, " pathways)"),
    show_rownames = TRUE, fontsize_row = fs_row, fontsize_col = 10,
    cluster_rows = TRUE, cluster_cols = FALSE, angle_col = 0
  )
  
  tp_means_by_group[[g]] <- tp_means
  message("Done group: ", g, " -> ", gdir)
}

# -----------------
# COMBINED heatmap across groups
# -----------------
if (length(tp_means_by_group)) {
  master_rows <- label_go_terms(sort(present_rows_union))
  combo_list <- lapply(names(tp_means_by_group), function(g) {
    x <- tp_means_by_group[[g]]
    x2 <- matrix(NA_real_, nrow=length(master_rows), ncol=ncol(x),
                 dimnames = list(master_rows, colnames(x)))
    common <- intersect(rownames(x), master_rows)
    x2[common, ] <- x[common, , drop=FALSE]
    colnames(x2) <- paste0(g, "|", colnames(x2))  # "Group|Timepoint"
    x2
  })
  combo <- do.call(cbind, combo_list)
  
  v_all <- matrixStats::rowVars(as.matrix(combo), na.rm = TRUE)
  ord <- order(v_all, decreasing = TRUE, na.last = TRUE)
  combo <- combo[ord, , drop = FALSE]
  
  write.csv(combo, file.path(out_dir, "GSVA_signaling_pathways_ALL_groups_timepoints.csv"))
  
  parts <- parse_group_time(colnames(combo))
  tp_order_all <- levels(meta[[timepoint_col]])
  ann_col_all <- data.frame(
    Group     = factor(parts$Group),
    Timepoint = factor(parts$Timepoint, levels = tp_order_all[tp_order_all %in% parts$Timepoint])
  )
  row.names(ann_col_all) <- colnames(combo)
  
  min_non_na <- max(2, ceiling(0.2 * ncol(combo)))
  keep <- rowSums(!is.na(combo)) >= min_non_na
  v2   <- matrixStats::rowVars(as.matrix(combo), na.rm = TRUE)
  keep <- keep & is.finite(v2) & (v2 > 0)
  combo_plot <- combo[keep, , drop = FALSE]
  
  do_cluster_rows <- nrow(combo_plot) > 1
  dist_rows <- NULL
  if (do_cluster_rows) {
    cor_mat <- suppressWarnings(stats::cor(t(combo_plot), use = "pairwise.complete.obs", method = "pearson"))
    if (all(is.na(cor_mat))) {
      do_cluster_rows <- FALSE
    } else {
      cor_mat[is.na(cor_mat)] <- 0
      dist_rows <- stats::as.dist(1 - cor_mat)
    }
  }
  
  n_paths_all <- nrow(combo_plot)
  fs_row_all  <- max(font_row_min, 12 - 0.03 * n_paths_all)
  
  save_heatmap_pdf(
    mat = combo_plot,
    file = file.path(out_dir, "GSVA_signaling_pathways_COMBINED.pdf"),
    ann_col = ann_col_all[colnames(combo_plot), , drop = FALSE],  # exact order
    main = paste0("Signaling pathways – combined (", n_paths_all, " pathways; columns = Group|Timepoint)"),
    show_rownames = TRUE, fontsize_row = fs_row_all, fontsize_col = 8,
    cluster_rows = do_cluster_rows, cluster_cols = FALSE, angle_col = 90,
    clustering_distance_rows = if (!is.null(dist_rows)) dist_rows else "euclidean"
  )
}

# -----------------
# Δ #1: Matched-vehicle Δ (treated − vehicle), per-vehicle type (E3/DMSO)
# -----------------
if (make_delta_heatmap_avg && length(tp_means_by_group)) {
  control_names <- names(tp_means_by_group)[grepl(infer_controls_regex, names(tp_means_by_group))]
  if (length(control_names)) {
    ctrl_mats <- tp_means_by_group[control_names]
    treated_groups <- setdiff(names(tp_means_by_group), control_names)
    if (length(treated_groups)) {
      tp_order_all <- levels(meta[[timepoint_col]])
      for (g in treated_groups) {
        M <- tp_means_by_group[[g]]
        C <- choose_ctrl_matrix(g, ctrl_mats)
        if (is.null(M) || is.null(C)) next
        rows <- intersect(rownames(M), rownames(C))
        tps  <- intersect(colnames(M), colnames(C))
        if (!length(rows) || !length(tps)) next
        
        D <- M[rows, tps, drop=FALSE] - C[rows, tps, drop=FALSE]
        v_delta <- matrixStats::rowVars(as.matrix(D), na.rm = TRUE)
        D <- D[order(v_delta, decreasing = TRUE), , drop = FALSE]
        
        vlab <- ifelse(is.null(group_vehicle_map[[g]]) || is.na(group_vehicle_map[[g]]), "Ctrl", group_vehicle_map[[g]])
        gdir <- file.path(out_dir, "Group", safe_group_dir(g))
        
        write.csv(D, file.path(gdir, sprintf("GSVA_focus_DELTA_vs_%s_timepoint_means.csv", vlab)))
        ann_col <- data.frame(Timepoint = factor(colnames(D), levels = tps))
        row.names(ann_col) <- colnames(D)
        fs_row_delta <- max(font_row_min, 12 - 0.03 * nrow(D))
        save_heatmap_pdf(
          mat = D,
          file = file.path(gdir, sprintf("GSVA_focus_heatmap_DELTA_vs_%s.pdf", vlab)),
          ann_col = ann_col,
          main = sprintf("Δ ssGSEA (Group %s − %s)", g, vlab),
          show_rownames = TRUE, fontsize_row = fs_row_delta, fontsize_col = 8,
          cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 90
        )
      }
    }
  }
}

# -----------------
# Δ #2: Paired Δ (treated − its matched control sample ID)
# -----------------
if (make_delta_heatmap_paired && length(sample_scores_long)) {
  if (!is.na(matched_ctrl_col) && matched_ctrl_col %in% names(meta)) {
    
    long_all <- dplyr::bind_rows(sample_scores_long)
    tp_order_all <- levels(meta[[timepoint_col]])
    
    for (g in names(sample_scores_long)) {
      if (grepl(infer_controls_regex, g)) next
      md_g <- meta %>% dplyr::filter(Group == g)
      treated_rows <- md_g %>%
        dplyr::filter(!is.na(.data[[dose_code_col]]),
                      grepl("^S[123]$", .data[[dose_code_col]], ignore.case = TRUE),
                      !is.na(.data[[matched_ctrl_col]]),
                      .data[[matched_ctrl_col]] %in% meta$sample_id) %>%
        dplyr::transmute(sample_id,
                         MatchCtrl = .data[[matched_ctrl_col]],
                         Timepoint = .data[[timepoint_col]])
      if (!nrow(treated_rows)) next
      
      Lg <- sample_scores_long[[g]] %>% dplyr::select(GOID, TERM, sample_id, ssGSEA)
      DC <- long_all %>% dplyr::select(GOID, sample_id, ssGSEA) %>%
        dplyr::rename(MatchCtrl = sample_id, ssGSEA_ctrl = ssGSEA)
      
      D <- dplyr::inner_join(Lg, treated_rows, by = "sample_id") %>%
        dplyr::inner_join(DC, by = c("GOID","MatchCtrl")) %>%
        dplyr::mutate(delta = ssGSEA - ssGSEA_ctrl)
      if (!nrow(D)) next
      
      D_sum <- D %>%
        dplyr::group_by(TERM, GOID, Timepoint) %>%
        dplyr::summarise(delta = mean(delta, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(Timepoint)
      
      delta_wide <- tidyr::pivot_wider(D_sum, id_cols = c(TERM, GOID),
                                       names_from = Timepoint, values_from = delta)
      
      rw  <- delta_wide %>% dplyr::select(-TERM, -GOID)
      if (!ncol(rw)) next
      mat <- as.matrix(rw)
      rownames(mat) <- label_go_terms(delta_wide$GOID)
      
      v <- matrixStats::rowVars(mat, na.rm = TRUE)
      mat <- mat[order(v, decreasing = TRUE), , drop = FALSE]
      
      gdir <- file.path(out_dir, "Group", safe_group_dir(g))
      write.csv(cbind(TERM = delta_wide$TERM, GOID = delta_wide$GOID, as.data.frame(mat)),
                file.path(gdir, "GSVA_focus_DELTA_PAIRED_vs_MATCHED_CONTROL.csv"),
                row.names = FALSE)
      
      ann_col_delta <- data.frame(Timepoint = factor(colnames(mat), levels = tp_order_all[tp_order_all %in% colnames(mat)]))
      row.names(ann_col_delta) <- colnames(mat)
      fs_row_delta <- max(font_row_min, 12 - 0.03 * nrow(mat))
      
      save_heatmap_pdf(
        mat = mat,
        file = file.path(gdir, "GSVA_focus_heatmap_DELTA_PAIRED_vs_MATCHED_CONTROL.pdf"),
        ann_col = ann_col_delta,
        main = sprintf("Δ ssGSEA (paired) – %s", g),
        show_rownames = TRUE, fontsize_row = fs_row_delta, fontsize_col = 8,
        cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 90
      )
    }
  } else {
    message("Per-group paired Δ requested, but no matched control column found; skipping per-group paired Δ.")
  }
}

# -----------------
# EXPORT: tidy long CSV across all groups (sample-level) + TIMEPOINT heatmaps
# -----------------
if (length(sample_scores_long)) {
  long_all <- dplyr::bind_rows(sample_scores_long)
  write.csv(long_all, file.path(out_dir, "GSVA_signaling_pathways_scores_long.csv"), row.names = FALSE)
  
  tp_dir <- file.path(out_dir, "Time_plots")
  dir.create(tp_dir, showWarnings = FALSE, recursive = TRUE)
  tp_levels_all <- levels(meta[[timepoint_col]])
  tp_present <- tp_levels_all[tp_levels_all %in% unique(long_all$Timepoint)]
  
  build_ann_col <- function(sample_ids) {
    ann <- meta %>% dplyr::filter(sample_id %in% sample_ids) %>%
      dplyr::mutate(
        Dose = if ("DoseTier" %in% names(meta)) as.character(DoseTier) else NA_character_,
        Vehicle = if (!is.na(vehicle_col) && vehicle_col %in% names(meta)) meta[[vehicle_col]][match(sample_id, meta$sample_id)] else NA_character_
      ) %>%
      dplyr::select(sample_id, Group, Dose, Vehicle)
    rownames(ann) <- ann$sample_id
    ann$sample_id <- NULL
    ann$Group   <- factor(ann$Group)
    if (!all(is.na(ann$Dose)))    ann$Dose    <- factor(ann$Dose, levels = c("S1","S2","S3","C"))
    if (!all(is.na(ann$Vehicle))) ann$Vehicle <- factor(ann$Vehicle, levels = c("E3","DMSO"))
    as.data.frame(ann, stringsAsFactors = FALSE)
  }
  
  for (tp in tp_present) {
    df_tp <- long_all %>% dplyr::filter(Timepoint == tp)
    wide <- df_tp %>%
      dplyr::select(TERM, GOID, sample_id, ssGSEA) %>%
      tidyr::pivot_wider(id_cols = c(TERM, GOID),
                         names_from = sample_id, values_from = ssGSEA)
    if (nrow(wide) == 0 || ncol(wide) <= 2) next
    
    mat <- as.matrix(wide[, -(1:2), drop = FALSE])
    min_non_na <- max(2, ceiling(0.3 * ncol(mat)))
    keep_rows <- rowSums(!is.na(mat)) >= min_non_na
    mat <- mat[keep_rows, , drop = FALSE]
    if (nrow(mat) == 0) next
    
    v <- matrixStats::rowVars(mat, na.rm = TRUE)
    ord <- order(v, decreasing = TRUE)
    mat <- mat[ord, , drop = FALSE]
    rownames(mat) <- label_go_terms(wide$GOID[keep_rows][ord])
    
    ann_col <- build_ann_col(colnames(mat))
    
    out_csv <- file.path(tp_dir, sprintf("GSVA_Time_%s_per-sample.csv", tp))
    row_info <- data.frame(
      TERM = wide$TERM[keep_rows][ord],
      GOID = wide$GOID[keep_rows][ord],
      stringsAsFactors = FALSE
    )
    write.csv(cbind(row_info, as.data.frame(mat)),
              out_csv, row.names = FALSE)
    
    fs_row <- max(font_row_min, 12 - 0.03 * nrow(mat))
    save_heatmap_pdf(
      mat = mat,
      file = file.path(tp_dir, sprintf("GSVA_Time_%s_per-sample_heatmap.pdf", tp)),
      ann_col = ann_col,
      main = sprintf("ssGSEA per sample – %s (%d pathways × %d samples)", tp, nrow(mat), ncol(mat)),
      show_rownames = TRUE, fontsize_row = fs_row, fontsize_col = 6,
      cluster_rows = TRUE, cluster_cols = TRUE, angle_col = 90
    )
  }
  message("Time-specific per-sample heatmaps written to: ", tp_dir)
}

# -----------------
# NTC significance (CONSERVATIVE): joint scoring + limma + broader FDR + effect-size floor
# -----------------
if (length(chosen_goids_global)) {
  message("\n[NTC] Conservative significance testing: joint scoring (Chem×TP) + limma + broader FDR…")
  
  get_chem <- function(g) sub("\\|.*$", "", g)
  get_dose <- function(g) sub("^.*\\|", "", g)
  
  control_names <- unique(meta$Group[grepl(infer_controls_regex, meta$Group)])
  treated_names <- unique(meta$Group[!grepl(infer_controls_regex, meta$Group) & grepl("\\|S[123]$", meta$Group)])
  chems <- sort(unique(get_chem(treated_names)))
  
  if (!length(chems)) {
    message("No treated chemicals found; skipping NTC.")
  } else {
    limma_results_all       <- list()  # original S1/S2/S3 contrasts
    limma_results_all_S12   <- list()  # NEW pooled S1+S2 contrasts
    
    for (chem in chems) {
      message("\n--- NTC: ", chem, " ---")
      dose_groups_chem <- intersect(treated_names, paste0(chem, "|", c("S1","S2","S3")))
      if (!length(dose_groups_chem)) next
      
      tp_levels_all <- levels(meta[[timepoint_col]])
      tps_present <- unique(meta[[timepoint_col]][meta$Group %in% c(dose_groups_chem, control_names)])
      tps_present <- tp_levels_all[tp_levels_all %in% tps_present]
      if (!length(tps_present)) next
      
      res_chem      <- list()
      res_chem_S12  <- list()  # NEW: pooled
      
      for (tp in tps_present) {
        md_tp <- meta %>% dplyr::filter(.data[[timepoint_col]] == tp)
        
        veh_for_chem <- if (!is.na(vehicle_col))
          unique(na.omit(toupper(md_tp[[vehicle_col]][md_tp$Group %in% dose_groups_chem])))
        else character(0)
        ctrl_tp  <- unique(md_tp$Group[grepl(infer_controls_regex, md_tp$Group)])
        ctrl_use <- if (length(veh_for_chem) && any(ctrl_tp %in% veh_for_chem)) intersect(ctrl_tp, veh_for_chem) else ctrl_tp
        
        trt_tp <- intersect(md_tp$Group, dose_groups_chem)
        if (!length(trt_tp) || !length(ctrl_use)) next
        
        sample_ctrl    <- md_tp$sample_id[md_tp$Group %in% ctrl_use]
        sample_by_dose <- lapply(trt_tp, function(g) md_tp$sample_id[md_tp$Group == g])
        names(sample_by_dose) <- trt_tp
        
        if (length(sample_ctrl) < ntc_min_samples_per_arm) next
        sample_by_dose <- sample_by_dose[vapply(sample_by_dose, length, integer(1)) >= ntc_min_samples_per_arm]
        if (!length(sample_by_dose)) next
        
        sample_ids <- c(sample_ctrl, unlist(sample_by_dose, use.names = FALSE))
        cnt <- counts[, sample_ids, drop = FALSE]
        if (ncol(cnt) < (ntc_min_samples_per_arm * 2)) next
        
        # VST / log transform
        if (is_integer_matrix(cnt)) {
          dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(cnt),
                                                colData = meta[match(colnames(cnt), meta$sample_id), ],
                                                design = ~ 1)
          dds <- suppressMessages(DESeq2::estimateSizeFactors(dds))
          vsd <- suppressMessages(DESeq2::vst(dds, blind = TRUE))
          expr <- SummarizedExperiment::assay(vsd)
        } else {
          expr <- log2(cnt + 1)
        }
        rownames(expr) <- sub("\\..*$", "", rownames(expr))
        
        # Gene sets present at this TP
        gsets <- gobp_sets_all[chosen_goids_global]
        gsets <- lapply(gsets, function(s) intersect(s, rownames(expr)))
        keep_set <- vapply(gsets, function(x) { lx <- length(x); !is.na(lx) && lx >= gsva_min_size && lx <= gsva_max_size }, logical(1))
        gsets <- gsets[keep_set]
        if (!length(gsets)) next
        
        # Score all selected samples together (joint scoring)
        gsva_scores <- run_gsva_ssgsea(expr, gsets, gsva_min_size, gsva_max_size)
        
        # Drop constant rows for stability
        rv <- matrixStats::rowVars(gsva_scores, na.rm = TRUE)
        gsva_scores <- gsva_scores[rv > 0, , drop = FALSE]
        if (!nrow(gsva_scores)) next
        
        # ------------------------------
        # Original contrasts: S1/S2/S3 vs Ctrl (kept exactly as before)
        # ------------------------------
        cond <- rep("Ctrl", length(sample_ctrl))
        for (g in names(sample_by_dose)) cond <- c(cond, rep(get_dose(g), length(sample_by_dose[[g]])))
        cond <- factor(cond, levels = c("Ctrl","S1","S2","S3")[c("Ctrl","S1","S2","S3") %in% unique(cond)])
        design <- stats::model.matrix(~ 0 + cond)
        colnames(design) <- sub("^cond", "", colnames(design))
        
        fit <- limma::lmFit(gsva_scores, design)
        doses_present <- setdiff(colnames(design), "Ctrl")
        if (length(doses_present)) {
          contr <- limma::makeContrasts(
            contrasts = vapply(doses_present, function(d) paste0(d, "-Ctrl"), character(1)),
            levels = design
          )
          fit2 <- limma::contrasts.fit(fit, contr)
          fit2 <- limma::eBayes(fit2)
          
          for (d in doses_present) {
            tt <- limma::topTable(fit2, coef = paste0(d, "-Ctrl"), number = Inf, sort.by = "none")
            tt$GOID  <- rownames(tt)
            tt$TERM  <- label_go_terms(tt$GOID)
            
            ctrl_cols <- sample_ctrl
            dose_cols <- sample_by_dose[[paste0(chem,"|", d)]]
            mean_ctrl <- rowMeans(gsva_scores[, ctrl_cols, drop = FALSE]); names(mean_ctrl) <- rownames(gsva_scores)
            mean_trt  <- rowMeans(gsva_scores[, dose_cols, drop = FALSE]); names(mean_trt)  <- rownames(gsva_scores)
            
            res <- tt %>%
              dplyr::transmute(
                Group     = paste0(chem, "|", d),
                Timepoint = tp,
                GOID, TERM,
                n_ctrl    = length(ctrl_cols),
                n_treat   = length(dose_cols),
                mean_ctrl = unname(mean_ctrl[GOID]),
                mean_treat= unname(mean_trt[GOID]),
                delta     = logFC,          # (dose - Ctrl)
                t_stat    = t,
                p_value   = P.Value,
                p_adj_local = adj.P.Val
              )
            res_chem[[paste(tp, d, sep="|")]] <- res
          }
        }
        
        # ------------------------------
        # NEW pooled contrast: (S1 + S2) vs Ctrl  → label as S12
        # ------------------------------
        if (ntc_pool_S1_S2) {
          s1_name <- paste0(chem, "|S1")
          s2_name <- paste0(chem, "|S2")
          s12_cols <- c(
            if (s1_name %in% names(sample_by_dose)) sample_by_dose[[s1_name]] else character(0),
            if (s2_name %in% names(sample_by_dose)) sample_by_dose[[s2_name]] else character(0)
          )
          s12_cols <- unique(s12_cols)
          
          if (length(sample_ctrl) >= ntc_min_samples_per_arm && length(s12_cols) >= ntc_min_samples_per_arm) {
            # Subset to Ctrl + S12 for this pooled test
            keep_cols <- c(sample_ctrl, s12_cols)
            gsva_pooled <- gsva_scores[, keep_cols, drop = FALSE]
            rv2 <- matrixStats::rowVars(gsva_pooled, na.rm = TRUE)
            gsva_pooled <- gsva_pooled[rv2 > 0, , drop = FALSE]
            if (nrow(gsva_pooled)) {
              cond2 <- factor(c(rep("Ctrl", length(sample_ctrl)), rep("S12", length(s12_cols))),
                              levels = c("Ctrl","S12"))
              design2 <- stats::model.matrix(~ 0 + cond2)
              colnames(design2) <- sub("^cond2", "", colnames(design2))
              
              fitp <- limma::lmFit(gsva_pooled, design2)
              contrp <- limma::makeContrasts(S12 - Ctrl, levels = design2)
              fitp2 <- limma::contrasts.fit(fitp, contrp)
              fitp2 <- limma::eBayes(fitp2)
              
              ttp <- limma::topTable(fitp2, coef = 1, number = Inf, sort.by = "none")
              ttp$GOID <- rownames(ttp)
              ttp$TERM <- label_go_terms(ttp$GOID)
              
              mean_ctrl_p <- rowMeans(gsva_pooled[, seq_len(length(sample_ctrl)), drop = FALSE]); names(mean_ctrl_p) <- rownames(gsva_pooled)
              mean_s12_p  <- rowMeans(gsva_pooled[, (length(sample_ctrl)+1):ncol(gsva_pooled), drop = FALSE]); names(mean_s12_p) <- rownames(gsva_pooled)
              
              res_p <- ttp %>%
                dplyr::transmute(
                  Group     = paste0(chem, "|S12"),
                  Timepoint = tp,
                  GOID, TERM,
                  n_ctrl    = length(sample_ctrl),
                  n_treat   = length(s12_cols),
                  mean_ctrl = unname(mean_ctrl_p[GOID]),
                  mean_treat= unname(mean_s12_p[GOID]),
                  delta     = logFC,     # (S12 - Ctrl)
                  t_stat    = t,
                  p_value   = P.Value,
                  p_adj_local = adj.P.Val
                )
              res_chem_S12[[paste(tp, "S12", sep="|")]] <- res_p
            }
          }
        }
        
      } # TP loop
      
      # ---------- finish per-chem: adjust p-values & save ----------
      if (length(res_chem)) {
        df_chem <- dplyr::bind_rows(res_chem)
        if (ntc_fdr_family == "by_group_by_time") {
          df_chem$p_adj <- df_chem$p_adj_local
        } else if (ntc_fdr_family == "by_group_across_time") {
          df_chem <- df_chem %>%
            dplyr::group_by(Group) %>%
            dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
            dplyr::ungroup()
        } else if (ntc_fdr_family == "by_chemical") {
          df_chem$p_adj <- p.adjust(df_chem$p_value, method = "BH")
        } else {
          df_chem$p_adj <- df_chem$p_adj_local
        }
        df_chem <- df_chem %>% dplyr::mutate(pass_effect = abs(delta) >= ntc_effect_size_floor)
        out_chem <- df_chem %>% dplyr::arrange(Group, Timepoint, p_adj)
        write.csv(out_chem, file.path(out_dir, paste0("NTC_GO_significance_limma_", chem, ".csv")), row.names = FALSE)
        limma_results_all[[chem]] <- out_chem
      }
      
      if (length(res_chem_S12)) {
        df_chem_S12 <- dplyr::bind_rows(res_chem_S12)
        if (ntc_fdr_family == "by_group_by_time") {
          df_chem_S12$p_adj <- df_chem_S12$p_adj_local
        } else if (ntc_fdr_family == "by_group_across_time") {
          df_chem_S12 <- df_chem_S12 %>%
            dplyr::group_by(Group) %>%
            dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
            dplyr::ungroup()
        } else if (ntc_fdr_family == "by_chemical") {
          df_chem_S12$p_adj <- p.adjust(df_chem_S12$p_value, method = "BH")
        } else {
          df_chem_S12$p_adj <- df_chem_S12$p_adj_local
        }
        df_chem_S12 <- df_chem_S12 %>% dplyr::mutate(pass_effect = abs(delta) >= ntc_effect_size_floor)
        out_chem_S12 <- df_chem_S12 %>% dplyr::arrange(Group, Timepoint, p_adj)
        write.csv(out_chem_S12, file.path(out_dir, paste0("NTC_GO_significance_limma_S12_", chem, ".csv")), row.names = FALSE)
        limma_results_all_S12[[chem]] <- out_chem_S12
      }
      
    } # chem loop
    
    # ---------- combined CSVs ----------
    if (length(limma_results_all)) {
      comb <- dplyr::bind_rows(limma_results_all) %>%
        dplyr::arrange(Group, Timepoint, p_adj)
      write.csv(comb, file.path(out_dir, "NTC_GO_significance_ALL_groups_LIMMA.csv"), row.names = FALSE)
      
      comb_sig <- comb %>% dplyr::filter(!is.na(p_adj), p_adj < 0.05, pass_effect)
      if (nrow(comb_sig) > 0) {
        write.csv(comb_sig %>% dplyr::arrange(Group, Timepoint, p_adj),
                  file.path(out_dir, "NTC_GO_significance_ONLY_significant_LIMMA.csv"),
                  row.names = FALSE)
        
        term_summary_all <- comb_sig %>%
          dplyr::mutate(Direction = ifelse(delta > 0, "up", "down")) %>%
          dplyr::group_by(GOID, TERM) %>%
          dplyr::summarise(
            best_FDR  = min(p_adj, na.rm = TRUE),
            n_sig     = dplyr::n(),
            directions= paste(sort(unique(Direction)), collapse = "/"),
            groups    = paste(sort(unique(Group)), collapse = ";"),
            timepoints= paste(sort(unique(Timepoint)), collapse = ";"),
            .groups   = "drop"
          ) %>%
          dplyr::arrange(best_FDR)
        write.csv(term_summary_all,
                  file.path(out_dir, "GSVA_GO_Treated_vs_Vehicle_TERM_summary_ALL_TP_CHEM_LIMMA.csv"),
                  row.names = FALSE)
      }
    }
    
    # ---------- NEW: combined CSVs for pooled S12 ----------
    if (length(limma_results_all_S12)) {
      comb_p <- dplyr::bind_rows(limma_results_all_S12) %>%
        dplyr::arrange(Group, Timepoint, p_adj)
      write.csv(comb_p, file.path(out_dir, "NTC_GO_significance_ALL_groups_LIMMA_S12.csv"), row.names = FALSE)
      
      comb_p_sig <- comb_p %>% dplyr::filter(!is.na(p_adj), p_adj < 0.05, pass_effect)
      if (nrow(comb_p_sig) > 0) {
        write.csv(comb_p_sig %>% dplyr::arrange(Group, Timepoint, p_adj),
                  file.path(out_dir, "NTC_GO_significance_ONLY_significant_LIMMA_S12.csv"),
                  row.names = FALSE)
        
        term_summary_all_p <- comb_p_sig %>%
          dplyr::mutate(Direction = ifelse(delta > 0, "up", "down")) %>%
          dplyr::group_by(GOID, TERM) %>%
          dplyr::summarise(
            best_FDR  = min(p_adj, na.rm = TRUE),
            n_sig     = dplyr::n(),
            directions= paste(sort(unique(Direction)), collapse = "/"),
            groups    = paste(sort(unique(Group)), collapse = ";"),
            timepoints= paste(sort(unique(Timepoint)), collapse = ";"),
            .groups   = "drop"
          ) %>%
          dplyr::arrange(best_FDR)
        write.csv(term_summary_all_p,
                  file.path(out_dir, "GSVA_GO_Treated_vs_Vehicle_TERM_summary_ALL_TP_CHEM_LIMMA_S12.csv"),
                  row.names = FALSE)
      }
    }
  }
}

message("\nAll done (conservative NTC, with pooled S1+S2). Outputs in: ", out_dir)