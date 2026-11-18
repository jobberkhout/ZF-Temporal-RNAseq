#######################################################
### DESeq2 contrasts: each treated group vs vehicle ###
#######################################################
# Writes the DESeq2 result tables in DEGs_NTC/ that GSEA_NTD.R ranks, and the all-sample PCA
# (PCA_Time.png, supplementary figure).
#
# Based on the R-ODAF DESeq2 template (Verheijen et al. 2022, Regul Toxicol Pharmacol 131:105143,
# doi:10.1016/j.yrtph.2022.105143). Steps:
#   1. Drop libraries with <= 5 million counts (removes 219-1305-T2-2-4D-S2, which has ~26,000).
#   2. For each contrast, fit one DESeq2 model (design ~ Label) on the treated group (2 samples)
#      and the controls of the same vehicle at the same stage (4 samples).
#   3. Keep genes with >= 1 CPM in at least 75% of the samples of at least one of the two groups.
#   4. Wald test, Benjamini-Hochberg FDR; Cook's distance cutoff and independent filtering are off.
#   5. DEG tables (padj < 0.01) additionally apply the R-ODAF quantile and spike rules. These DEG
#      tables are not used downstream; GSEA uses the full result tables.
#
# Paths are relative to the repository root (open ZF-Temporal-RNAseq.Rproj, or setwd() to the repo root).
# INPUT
#   - ZF_RNA_Raw_Count.txt     merged RSEM expected counts (genes x samples); rounded before DESeq2
#   - metadata/Samplekey.txt   sample metadata; the Label column defines the groups (e.g. CYCL_S1_T3, DMSO_C_T3)
# OUTPUT
#   - DEGs_NTC/<group>_vs_<control>_FDR_0.01_all_genes.txt   all tested genes, input for GSEA_NTD.R
#   - Output/DESeq2/   the same tables in R-ODAF format (no header for the ID column), DEG tables,
#                      per-contrast normalised counts, heatmaps of the top 50 DEGs, all-sample PCAs,
#                      Norm_Data_all.txt and report.txt (number of DEGs per contrast)
# The heatmaps look up gene symbols at Ensembl (biomaRt), so the script needs an internet connection.
# Sample 315-2005-T5-TEBU-S1 is kept, as in the published results. It is an outlier in a PCA of the
# 24 hpf samples alone; a later run that excluded it (and the TEBU_S1_T5 contrast) was not used.

library(DESeq2)
require("DESeq2")
require("edgeR")
require("pheatmap")
library(ggplot2)
library(ggrepel)
library(reshape2)
library(biomaRt)     # map IDs → symbols


###################################################################################
###################################################################################
# PARAMETERS

# File locations (relative to the repository root)
	repo_dir  <- getwd()
	outputdir <- file.path(repo_dir, "Output", "DESeq2")
	degdir    <- file.path(repo_dir, "DEGs_NTC")
	dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)
	dir.create(degdir, recursive = TRUE, showWarnings = FALSE)
# Names of files to load
	SampleDataFile <- file.path(repo_dir, "ZF_RNA_Raw_Count.txt")      # tab-delimited, merged RSEM genes.results (expected counts)
	SampleKeyFile  <- file.path(repo_dir, "metadata", "Samplekey.txt") # tab-delimited; sample_id = column names of the count matrix, Label = group

# Groups to compare: every chemical x dose (S3 low, S2 mid, S1 high) x stage (T1-T5),
# each against the controls of its vehicle at the same stage.
	contrasts <- expand.grid(dose = c("S3", "S2", "S1"), stage = paste0("T", 1:5),
	                         chem = c("GLYP", "CAFF", "2-4D", "IMID", "CYCL", "RACI", "TEBU", "THIR"),
	                         stringsAsFactors = FALSE)
	contrasts$vehicle <- ifelse(contrasts$chem %in% c("GLYP", "CAFF", "2-4D", "IMID"), "E3", "DMSO")
	# 2-4D_S2_T2 is skipped: one of its two libraries (219-1305-T2-2-4D-S2) fails the read-depth filter below.
	contrasts <- contrasts[!(contrasts$chem == "2-4D" & contrasts$dose == "S2" & contrasts$stage == "T2"), ]
	Samp4compare <- paste(contrasts$chem, contrasts$dose, contrasts$stage, sep = "_")  # treated group (Label)
	Cont4compare <- paste(contrasts$vehicle, "C", contrasts$stage, sep = "_")         # matching control group (Label)
	DESIGN<- "Label"	# column of the sample key that defines the groups to be compared

# Set analysis ID. This ID will be used as prefix for the output files
	analysisID <-"Job"
# Specify used platform/technology for data generation:
	Platform <- "RNA-Seq" # Specify "RNA-Seq" or "TempO-seq"

###################################################################################
###################################################################################

#Set parameters according to platform
if (Platform=="RNA-Seq"){
	minCoverage <- 5000000
	MinCount<- 1
	pAdjValue<- 0.01 
} else if (Platform=="TempO-seq") {
	minCoverage <- 1000000
	MinCount<- 0.5
	pAdjValue<- 0.01 
} else { print("Platform/technology not recognized") }

# Load input files 
sampleData <- read.table(SampleDataFile, sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=1, check.names=F)
DESeqDesign <- read.table(SampleKeyFile, stringsAsFactors=FALSE, sep="\t", header=TRUE, row.names="sample_id", check.names=F)

# The supplementary all-sample PCA in the paper (Output/Figures/PCA_Time.png) came from a later run
# that excluded the 24 hpf outlier below. Uncomment these two lines to reproduce that PCA; the
# TEBU_S1_T5 contrast then has one sample left and cannot be tested.
# sampleData <- sampleData[,-grep("315-2005-T5-TEBU-S1",colnames(sampleData))]   #remove a clear outlier
# DESeqDesign <- DESeqDesign[-grep("315-2005-T5-TEBU-S1",rownames(DESeqDesign)),] #remove a clear outlier

# All plots and tables below are written to Output/DESeq2/ (DEGs_NTC/ files go to degdir)
setwd(outputdir)

# savePlot() only works on Windows; this saves the current plot as a PNG on any system.
# With Rscript the default plot device does not record plots, so open one that does.
if (!interactive()) { pdf(NULL); dev.control("enable") }
save_png <- function(filename, width = 10, height = 7) {
	dev.copy(png, filename = filename, width = width, height = height, units = "in", res = 150)
	invisible(dev.off())
}

# DEGs_NTC/ gets the full results table with a named ID column (input for GSEA_NTD.R).
write_deg_ntc <- function(res, FileName) {
	out <- data.frame(ENSEMBL = rownames(res), as.data.frame(res), check.names = FALSE)
	write.table(out, file = file.path(degdir, paste0(FileName, "_all_genes.txt")), sep = "\t", quote = FALSE, row.names = FALSE)
}

NORM_TYPE<-paste0(analysisID, "_DESeq2_", Platform)
print(NORM_TYPE)

# First data clean-up: replace NA & remove samples with total readcount < threshold
sampleData[ is.na(sampleData) ] <- 0 # this remove a sample with almost no read.
low_cov <- colnames(sampleData)[colSums(sampleData) <= minCoverage]
if (length(low_cov)) print(paste("Removed for", minCoverage, "or fewer counts:", paste(low_cov, collapse = ", ")))
sampleData <- sampleData[,(colSums(sampleData)> minCoverage)]

DESeqDesign <- DESeqDesign[rownames(DESeqDesign) %in% colnames(sampleData),]

sampleData <- sampleData[,order(colnames(sampleData))]
DESeqDesign <- DESeqDesign[order(rownames(DESeqDesign)),]

samplekey <- DESeqDesign
mRNAtot <- sampleData

raw_count <- colSums(mRNAtot)
barplot(raw_count,las=2,cex.names=0.6, main="Raw_Mapped_Counts")
save_png("Raw_Mapped_Counts.png")

report <- matrix(data=NA, ncol=2, nrow=length(Cont4compare))
rownames(report) <- paste((Samp4compare),(Cont4compare),sep="_vs_")

dds <- DESeqDataSetFromMatrix(countData = round(mRNAtot), colData = samplekey, design = ~ compound_code)
dds <- DESeq(dds, quiet=TRUE)

# perform a PCA plot. You can change the condition assessed by the pca plot by changing the value between [] from 1 to 3 (corresponding to the column number of the metadata file)
vst <- vst(dds)
plotPCA(vst,intgroup=c(colnames(samplekey)[2]))
save_png("PCA_Compound.png")

pcaData <- plotPCA(vst, intgroup = c(colnames(samplekey)[2]), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Create the PCA plot with sample names annotated
p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) +
     geom_point(size = 3) +
     geom_text_repel(aes(label = name), size = 3) +  # Annotate with sample names
     labs(x = paste0("PC1: ", percentVar[1], "% variance"),
          y = paste0("PC2: ", percentVar[2], "% variance"),
          title = "PCA of VST data") +
     theme_minimal()
p
save_png("PCA_annotated_Compound.png")

# perform a PCA plot. You can change the condition assessed by the pca plot by changing the value between [] from 1 to 3 (corresponding to the column number of the metadata file)
vst <- vst(dds)
plotPCA(vst,intgroup=c(colnames(samplekey)[6]))
save_png("PCA_Time.png")

pcaData <- plotPCA(vst, intgroup = c(colnames(samplekey)[6]), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Create the PCA plot with sample names annotated
p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) +
     geom_point(size = 3) +
     geom_text_repel(aes(label = name), size = 3) +  # Annotate with sample names
     labs(x = paste0("PC1: ", percentVar[1], "% variance"),
          y = paste0("PC2: ", percentVar[2], "% variance"),
          title = "PCA of VST data") +
     theme_minimal()
p
save_png("PCA_annotated_Time.png")

# perform a PCA plot. You can change the condition assessed by the pca plot by changing the value between [] from 1 to 3 (corresponding to the column number of the metadata file)
vst <- vst(dds)
plotPCA(vst,intgroup=c(colnames(samplekey)[4]))
save_png("PCA_Dose.png")

pcaData <- plotPCA(vst, intgroup = c(colnames(samplekey)[4]), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Create the PCA plot with sample names annotated
p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) +
     geom_point(size = 3) +
     geom_text_repel(aes(label = name), size = 3) +  # Annotate with sample names
     labs(x = paste0("PC1: ", percentVar[1], "% variance"),
          y = paste0("PC2: ", percentVar[2], "% variance"),
          title = "PCA of VST data") +
     theme_minimal()
p
save_png("PCA_annotated_Dose.png")

norm_data <- counts(dds[rownames(mRNAtot)],normalized=TRUE) 
write.table(norm_data,file=paste0("Norm_Data_all.txt"), sep="\t", quote=FALSE)

##########
# DESeq2 #
##########

for (x in 1:length(Samp4compare)){	## for all comparisons to be done	
	condition1<- Cont4compare[x]	    		
	condition2<- Samp4compare[x]  

	DE_Design <- matrix(data=NA, ncol=2)
	DE_Design <- DESeqDesign [c(grep(condition1,DESeqDesign[,DESIGN[1]]), grep(condition2,DESeqDesign[,DESIGN[1]])),]
	samples <- sampleData[, rownames(DE_Design) ]

	########### 
	print(paste(condition2, " vs ", condition1, ":", NORM_TYPE))		

	colnames(samples)<-NULL
	if (DESIGN[1] == "Label") {
		dds <- DESeqDataSetFromMatrix(countData = round(samples), colData = as.data.frame(DE_Design), design = ~ Label)
	} else {print("Setting of DESIGN failed. Please check code: dds <- DESeqDataSetFromMatrix(countData = round(samples), colData = as.data.frame(DE_Design), design = ~ Compound)}.")}

	print("Wait... (dds step executing)")				
	dds <- DESeq(dds, quiet=TRUE)

#DEmRNA_fdr_FC <- subset(DEmRNA_fdr,abs(log2FoldChange) > 1.5)
#nrow(DEmRNA_fdr_FC)

#Filter low readcounts (genes not meeting the condition : at least one condition with 75% of the samples above 1 CPM)
	
SampPerGroup<-table(DE_Design[,DESIGN])
idx<-FlagSpike<-NameRows<-NULL
Counts<-counts(dds, normalized=TRUE)
CPMdds<-cpm(counts(dds, normalized=FALSE))
	
Filter <- matrix(data=NA, ncol=3, nrow= nrow(Counts))
rownames(Filter) <- rownames(Counts)
colnames(Filter) <- c("Low","quantile","spike")
		
print(paste0("Filtering genes: 75% of at least 1 group need to be above ", MinCount, " CPM"))
		
	for (gene in 1:nrow(dds)) {
		
		# Check the low count condition
		
		CountsPass<-NULL
			for (group in 1:length(SampPerGroup)) { 
				sampleCols<-grep(dimnames(SampPerGroup)[[1]][group],DE_Design[,DESIGN])
				Check<-sum(CPMdds[gene,sampleCols] >= MinCount)>= 0.75*SampPerGroup[group]
				CountsPass<-c(CountsPass, Check)
		}
# bypass low count filtering for testing : all condition to <- 1. Change back the else to 0		
		
	if ( sum(CountsPass) > 0 ) {Filter[gene,1] <- 1 }	else { Filter[gene,1] <- 0 }
			
		}	

		compte <- Counts[Filter[,1] == 1,]
		Filter <- Filter[rownames(Filter) %in% rownames(compte),]

print(paste("Low counts filterting removed ",nrow(dds)- nrow(Filter)," genes from the ",nrow(dds)," assessed. ", nrow(Filter)," genes remaining",sep=""))

print("Obtaining the DESeq2 results")

	res <- results(dds[rownames(compte),], alpha=0.05, cooksCutoff =F, independentFiltering=F, contrast=c(DESIGN[1], condition2, condition1), pAdjustMethod= 'fdr')   

	setwd(outputdir)
	
	#Save output tables		
	norm_data <<- counts(dds[rownames(compte)],normalized=TRUE) 
	DEsamples <<- subset(res,res$padj < pAdjValue)	

	DECounts <- compte[rownames(compte) %in% rownames(DEsamples),]
	Filter <- Filter[rownames(Filter) %in% rownames(DECounts),]
	
	
if (nrow(DEsamples) > 1) {
	
		print("Check median against third quartile" )
		print("AND")
		print("Check the presence of a spike" )

		sampleColsg1 <- grep(dimnames(SampPerGroup)[[1]][1],DE_Design[,DESIGN])
		sampleColsg2 <- grep(dimnames(SampPerGroup)[[1]][2],DE_Design[,DESIGN])
		
		for (gene in 1:nrow(DECounts)) {
		
		# Check the median against third quantile
		quantilePass <-NULL
		
		Check <- median(DECounts[gene,sampleColsg1]) > quantile(DECounts[gene,sampleColsg2], 0.75)[[1]]
		quantilePass <-c(quantilePass, Check)
		Check <- median(DECounts[gene,sampleColsg2]) > quantile(DECounts[gene,sampleColsg1], 0.75)[[1]]
		quantilePass <-c(quantilePass, Check)
			
		if ( sum(quantilePass) > 0 ) {Filter[gene,2] <- 1 }	else { Filter[gene,2] <- 0 }

		# Check for spike 
		spikePass <- NULL		
		
if (sum(DECounts[gene,sampleColsg1]) > sum(DECounts[gene,sampleColsg2])) {spikePass <- (max(DECounts[gene,sampleColsg1])/sum(DECounts[gene,sampleColsg1])) >= 1.4*(SampPerGroup[1])^(-0.66)} 
	else {spikePass <- (max(DECounts[gene,sampleColsg2])/sum(DECounts[gene,sampleColsg2])) >= 1.4*(SampPerGroup[2])^(-0.66)} 
		if ( sum(spikePass) > 0 ) {Filter[gene,3] <- 0 } else { Filter[gene,3] <- 1 }
	}		

	DECounts_real <- subset(DEsamples, rownames(DEsamples) %in% rownames(Filter[rowSums(Filter) == 3,]))
	
	DECounts_relev <- DEsamples[Filter[,1] == 0 ,]
	DECounts_no_quant <- DEsamples[Filter[,2] == 0 ,]
	DECounts_spike <- DEsamples[Filter[,3] == 0,]
	
	relev <- subset(Filter,Filter[,1] == 0)
	no_quant <-  subset(Filter,Filter[,2] == 0)
	spike <-  subset(Filter,Filter[,3] == 0)

	
print(paste("A total of ",nrow(DECounts_real), " DEGs were selected, after ",nrow(DECounts_no_quant)," genes(s) removed by the quantile rule and ", nrow(DECounts_spike)," gene(s) with a spike",sep=""))
	
DECounts_real <- DECounts_real[order(DECounts_real$padj),]

res <- res[order(res$padj),]

FileName<-paste(condition2,"vs",condition1, "FDR", pAdjValue, sep="_")
	
write.table(DECounts_real,file=paste0(FileName,"_DEG_table.txt"), sep="\t", quote=FALSE)

print(paste("A total of ",nrow(DECounts_real), " DEGs were selected, after ",nrow(DECounts_no_quant)," genes(s) removed by the quantile rule and ", nrow(DECounts_spike)," gene(s) with a spike",sep=""))
	
# perform a heat Map (for all DEGs):
DECounts_real <- DECounts_real[order(DECounts_real$padj),]

norm_data_heatmap <- subset(
  norm_data,
  rownames(norm_data) %in% rownames(DECounts_real)[1:min(nrow(DECounts_real), 50)]
)


## ---- Look the IDs up at Ensembl ------------------------------------------
ensembl_ids <- sub("\\..*$", "", rownames(norm_data_heatmap))
mart <- useEnsembl(biomart = "genes", dataset = "drerio_gene_ensembl")

gene_map <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters    = "ensembl_gene_id",
    values     = unique(ensembl_ids),
    mart       = mart
)

gene_map$external_gene_name[gene_map$external_gene_name == "" | is.na(gene_map$external_gene_name)] <- gene_map$ensembl_gene_id[gene_map$external_gene_name == "" | is.na(gene_map$external_gene_name)]

## ---- 4. Build a vector we can index with ------------------------------------
# named vector: names = IDs, values = symbols
id2sym <- setNames(gene_map$external_gene_name, gene_map$ensembl_gene_id)

## ---- 5. Replace rownames, keep the ID if no symbol was found ----------------
new_names <- id2sym[ensembl_ids]             # symbols where they exist, NA otherwise
new_names[is.na(new_names)] <- ensembl_ids[is.na(new_names)]

# Optional: make any duplicates unique so pheatmap doesn’t complain
dup <- duplicated(new_names)
new_names[dup] <- paste0(new_names[dup], "_", ave(dup, new_names, FUN = seq_along)[dup])

rownames(norm_data_heatmap) <- new_names

## ---- 6. Plot exactly as before ----------------------------------------------
norm_data_heatmap <- norm_data_heatmap + 0.01            # your own code
Log_heatmap        <- log2(norm_data_heatmap)

my_annot           <- data.frame(Label = DE_Design$Label,
                                 row.names = rownames(DE_Design))

pheatmap(Log_heatmap,
         show_rownames   = TRUE,
         annotation_col  = my_annot)

norm_data_heatmap = norm_data_heatmap + 0.01 # pseudo_count for allowing log when value is 0

my_annot <- as.data.frame(DE_Design$Label)
Log_heatmap <- log2(norm_data_heatmap)
rownames(my_annot) <- rownames(DE_Design)
colnames(my_annot) <- "Label"
pheatmap(Log_heatmap, show_rownames=T, annotation_col = my_annot)

save_png(paste(FileName,"Heatmap_with_names.png",sep="_"))	
	
write.table(norm_data,file=paste0(FileName, "_Norm_Data.txt"), sep="\t", quote=FALSE)
write.table(norm_data,file=paste0("_Norm_Data.txt"), sep="\t", quote=FALSE)
write.table(res,file=paste0(FileName,"_all_genes.txt"), sep="\t", quote=FALSE)
write_deg_ntc(res, FileName)

report[x,1] <- nrow(DECounts_real)
}

else {

res <- res[order(res$padj),]

FileName<-paste(condition2,"vs",condition1, "FDR", pAdjValue, sep="_")
write.table(norm_data,file=paste0(FileName, "_Norm_Data.txt"), sep="\t", quote=FALSE)
write.table(res,file=paste0(FileName,"_all_genes.txt"), sep="\t", quote=FALSE)
write_deg_ntc(res, FileName)
report[x,1] <- nrow(DEsamples)

} 


	print("DESeq2 Done")
		
}

	write.table(report,file="report.txt", sep="\t", quote=FALSE)
