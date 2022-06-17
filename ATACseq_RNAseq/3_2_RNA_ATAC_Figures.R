# Author: Joern Pezoldt
# Advice: Maria Litovchenko, Vincent Gardeaux
# 12.08.2018
# Function:
# 1) 
# 

#####
#Libraries & PATHS & Data & Global_Variables
#####

require(data.table)
require(ggplot2)
#require(ggfortify)
require(limma)
require(DESeq2)
library(pheatmap)
library(dplyr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPpeakAnno)
library(org.Mm.eg.db)
library(biomaRt)
library(foreign)
library(GO.db)
library(topGO)
library(org.Mm.eg.db)
library(GenomicFeatures)

# Input required: Set directory
# 1) Counts per peak per replicate
path_input <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/homer/Overlap_Group_Merged/Run_3_Outlier_excluded"
PATH_input_TFdb_Riken <- "/home/pezoldt/NAS2/pezoldt/Data/Databases/Riken_TFdB/2019-02-08_Riken_TFdb_curated.txt"
# 2) Common Peak regions
path_common_peak_regions <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/peaks/broad/MACS2_broad/Run_3/ATAC_FSC_all_broad_merged_peaks.bed"
# 3) Output path
path_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/DESeq2/Run_3_Outlier_excluded"
path_output_Motif <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Regions_of_Interest"
path_output_GO_DEG <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Submission/Heatmap/Expression_GO"
path_output_DARs <- "/home/pezoldt/NAS2/pezoldt//Analysis/ATACseq/ATAC_FSC_all/Submission"

# 4) RNAseq DESeq2 analysis
path_RNAseq_DESeq2 <- "/home/pezoldt/NAS2/pezoldt/Data/RNAseq/2017_FSC_LEC_BEC_mLN_pLN_GF_SPF/FSC"
# 5) Results of DAR and DEG overlap
path_DAR_DEG <- "/home/pezoldt/NAS2/pezoldt/Analysis/01_Integrated/DARs_DEGs"
# Input required: 
name = "ATAC_FSC_all"
paste(path_input, "/",name,".txt",sep="")
#read count table
counts <- fread(paste(path_input, "/",name,".txt",sep=""))

#Global variables
#Param
log2FC_RNA = 1.0
padj = 0.05

log2FC_ATAC = 1.0

####################
#Load expression data
####################
#load all DESeq2 files from directory to list
setwd(path_RNAseq_DESeq2)
names = list.files(pattern="*.csv")
myfiles = lapply(names, read.csv)

###############
#Global variables
###############
#minimal replicate number of sample
min_rep_num = 3
log2FC = 1.0
padj = 0.05
min_RPKM = 1

#####
#QC
#####

# peaks never found
which(rowSums(counts)==0) ## none (expected)
# peaks found in < 50 % of samples
# Input required: Number of samples
number_of_samples <- 13

# peaks with zero count per sample 
zero <- colSums(counts[,-1]==0)
zero.df <- data.frame(colnames(counts[,-1]), zero)
colnames(zero.df) <- c("cell", "zero")

ggplot(zero.df) + geom_col(aes(x=(reorder(zero.df$cell, zero.df$zero)), y= zero.df$zero)) +
  theme(axis.text.x = element_text(angle = 50, hjust=1))

#Store rownames
row_counts <- rownames(counts)

# transform into count matrix
counts.mat <-(counts[,-1])
counts.mat <- as.matrix(counts.mat)
rownames(counts.mat) <- row_counts

#####
#DEseq2 ATAC-seq
#####
#Define conditons 
condition <- factor(substr(colnames(counts.mat),1,nchar(colnames(counts.mat))-1))
dds <- DESeqDataSetFromMatrix(counts.mat, DataFrame(condition),~condition)
dds.factors <- estimateSizeFactors(dds)

#Perform DESeq2
dds <- DESeq(dds)
res <- results(dds)

#generate dds-resuls objects for relevant comparisons
mLN_SPF_GF <- results(dds, contrast=c("condition","ATAC_mLNSPF","ATAC_mLNGF"))
pLN_SPF_GF <- results(dds, contrast=c("condition","ATAC_pLNSPF","ATAC_pLNGF"))
SPF_mLN_pLN <- results(dds, contrast=c("condition","ATAC_mLNSPF","ATAC_pLNSPF"))
GF_mLN_pLN <- results(dds, contrast=c("condition","ATAC_mLNGF","ATAC_pLNGF"))


#Perform RPKM calculation
peak_regions <- read.delim(path_common_peak_regions, header = FALSE)
colnames(peak_regions) <- c("chr","start","end")
#make Granges object all samples are computed on the same Granges object thus all have the same 
gr_peak_regions <- toGRanges(peak_regions, seqnames.field=c("chr"), start.field="start", end.field="end")
rowRanges(dds) <- gr_peak_regions

fpkm_dds <- fpkm(dds)
fpm_dds <- fpm(dds)

#QC
#MA
par(mfrow = c(2,2))
plotMA(mLN_SPF_GF, ylim = c(-6, 6), main = "mLN_SPF_GF")
plotMA(pLN_SPF_GF, ylim = c(-6, 6), main = "pLN_SPF_GF" )
plotMA(SPF_mLN_pLN, ylim = c(-6, 6), main = "SPF_mLN_pLN" )
plotMA(GF_mLN_pLN, ylim = c(-6, 6), main = "GF_mLN_pLN" )

#DispEsts
plotDispEsts( dds, ylim = c(1e-6, 1e1), main = "All" )


#Hits pValue
par(mfrow = c(2,2))
hist( mLN_SPF_GF$pvalue, breaks=20, col="grey", main = "mLN_SPF_GF"  )
hist( pLN_SPF_GF$pvalue, breaks=20, col="grey", main = "pLN_SPF_GF"  )
hist( SPF_mLN_pLN$pvalue, breaks=20, col="grey", main = "SPF_mLN_pLN"  )
hist( GF_mLN_pLN$pvalue, breaks=20, col="grey", main = "GF_mLN_pLN"  )

par(mfrow = c(1,1))
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

#Get log2FC and pValue for all peaks across all chosen comparisons
dds <- DESeq(dds)
res <- results(dds)
DAR_mLN_SPF_GF <- as.data.frame(results(dds, contrast=c("condition","ATAC_mLNSPF","ATAC_mLNGF"), format = c("DataFrame")))
DAR_mLN_SPF_GF <- DAR_mLN_SPF_GF[,c("log2FoldChange","padj")]
colnames(DAR_mLN_SPF_GF) <- c(paste("log2FC_","mLN_SPF_GF",sep=""),paste("padj_","mLN_SPF_GF",sep=""))

DAR_pLN_SPF_GF <- as.data.frame(results(dds, contrast=c("condition","ATAC_pLNSPF","ATAC_pLNGF"), format = c("DataFrame")))
DAR_pLN_SPF_GF <- DAR_pLN_SPF_GF[,c("log2FoldChange","padj")]
colnames(DAR_pLN_SPF_GF) <- c(paste("log2FC_","pLN_SPF_GF",sep=""),paste("padj_","pLN_SPF_GF",sep=""))

DAR_SPF_mLN_pLN <- as.data.frame(results(dds, contrast=c("condition","ATAC_mLNSPF","ATAC_pLNSPF"), format = c("DataFrame")))
DAR_SPF_mLN_pLN <- DAR_SPF_mLN_pLN[,c("log2FoldChange","padj")]
colnames(DAR_SPF_mLN_pLN) <- c(paste("log2FC_","SPF_mLN_pLN",sep=""),paste("padj_","SPF_mLN_pLN",sep=""))


DAR_GF_mLN_pLN <- as.data.frame(results(dds, contrast=c("condition","ATAC_mLNGF","ATAC_pLNGF"), format = c("DataFrame")))
DAR_GF_mLN_pLN <- DAR_GF_mLN_pLN[,c("log2FoldChange","padj")]
colnames(DAR_GF_mLN_pLN) <- c(paste("log2FC_","GF_mLN_pLN",sep=""),paste("padj_","GF_mLN_pLN",sep=""))

#Table of all conditions
DARs <- cbind(id = rownames(DAR_mLN_SPF_GF), DAR_mLN_SPF_GF, DAR_pLN_SPF_GF, DAR_SPF_mLN_pLN, DAR_GF_mLN_pLN)
rownames(DARs) <- c()

#####
#Associate genes with DESeq results
#####
#Load the genomic regions associated with ids
# Common peak track for Experiment is annotated using homer
regions <- read.delim(path_common_peak_regions, header = FALSE)
colnames(regions) <- c("chr","start","end")
DARs_regions <- cbind(regions, DARs)
#replace chromosome numbers with 1 -> chr1
DARs_regions$chr <- paste("chr",DARs_regions$chr, sep="")
rownames(DARs_regions) <- c()
#DARs with FPKM
DARs_regions <- cbind(DARs_regions,fpkm_dds)

#####
#Peak/gene location
#####
library("GenomicFeatures")

#All features
mm10_TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#TSS sites BioMart
mm10 = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl") 
TSS.mouse.mm10 = getAnnotation(mart=mm10, featureType="TSS")
Exon.mouse.mm10 = getAnnotation(mart=mm10, featureType="Exon")
#mm10 whole database
#annoData_exon <- toGRanges(mm10_TxDb, feature=c("exon"))
#Grange object from peak file
gr_DARs_regions <- toGRanges(DARs_regions, names = DARs_regions$id)
#Annotate Peaks
gr_DARs_regions_anno_TSS <- annotatePeakInBatch(gr_DARs_regions, AnnotationData=TSS.mouse.mm10)
#add gene name
gr_DARs_regions_anno_TSS <- addGeneIDs(annotatedPeak=gr_DARs_regions_anno_TSS, 
                                       feature_id_type="ensembl_gene_id",
                                       orgAnn="org.Mm.eg.db", 
                                       IDs2Add="symbol")
#Generate dataframe
DARs_features <- as.data.frame(gr_DARs_regions_anno_TSS)
#Number of features
length(DARs_features$insideFeature[DARs_features$insideFeature == "upstream"])

#Commensal specific DARs
DARs_features_mLN_Commensal <- subset(DARs_features, abs(log2FC_mLN_SPF_GF) >= log2FC_ATAC & padj_mLN_SPF_GF <= 0.05)
write.table(DARs_features_mLN_Commensal, paste(path_output_DARs,"/DARs_features_mLN_Commensal.txt", sep = ""), dec = ".", sep = "\t")
unique(DARs_features_mLN_Commensal$symbol)
DARs_features_pLN_Commensal <- subset(DARs_features, abs(log2FC_pLN_SPF_GF) >= log2FC_ATAC & padj_pLN_SPF_GF <= 0.05)
write.table(DARs_features_pLN_Commensal, paste(path_output_DARs,"/DARs_features_pLN_Commensal.txt", sep = ""), dec = ".", sep = "\t")
unique(DARs_features_pLN_Commensal$symbol)

####################
#Heatmaps for DARs
####################
# Note: Figures for Publication
#####
#SPF
#####
DARs_features_SPF <- subset(DARs_features, abs(log2FC_SPF_mLN_pLN) >= log2FC_ATAC & padj_SPF_mLN_pLN <= padj)
DARs_features_SPF_NA <- DARs_features_SPF[!is.na(DARs_features_SPF$symbol),]
#Heatmap
data_heatmap <- DARs_features_SPF_NA[,c("id","ATAC_mLNSPF1","ATAC_mLNSPF2","ATAC_mLNSPF3","ATAC_mLNSPF4",
                                        "ATAC_pLNSPF1","ATAC_pLNSPF2","ATAC_pLNSPF3")]
rownames(data_heatmap) <- data_heatmap$id
data_heatmap <- as.matrix(data_heatmap[,2:ncol(data_heatmap)])
data_heatmap[data_heatmap == 0] <- 0.4
data_heatmap <- log2(data_heatmap)

min(data_heatmap)
pheatmap(data_heatmap, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = TRUE, cellwidth = 10,
         scale = "row", border_color = "black", color = colorRampPalette(c("white","grey", "darkgreen"), space="rgb")(128),
         main = "DARs TSS associated SPF")

#####
#GF
#####
DARs_features_GF <- subset(DARs_features, abs(log2FC_GF_mLN_pLN) >= log2FC_ATAC & padj_GF_mLN_pLN <= padj)
DARs_features_GF_NA <- DARs_features_GF[!is.na(DARs_features_GF$symbol),]
#Heatmap
data_heatmap <- DARs_features_GF_NA[,c("id",
                                       "ATAC_mLNGF1","ATAC_mLNGF2","ATAC_mLNGF3",
                                       "ATAC_pLNGF1","ATAC_pLNGF3")]
rownames(data_heatmap) <- data_heatmap$id
data_heatmap <- as.matrix(data_heatmap[,2:ncol(data_heatmap)])
data_heatmap[data_heatmap == 0] <- 0.4
min(data_heatmap)
data_heatmap <- log2(data_heatmap)

min(data_heatmap)
pheatmap(data_heatmap, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = FALSE, cluster_cols = TRUE, cellwidth = 10,
         scale = "row", border_color = "black", color = colorRampPalette(c("white","grey", "darkgreen"), space="rgb")(128),
         main = "DARs TSS associated SPF")

#####
#GF & SPF
#####
DARs_features_GF_SPF <- subset(DARs_features, (abs(log2FC_GF_mLN_pLN) >= log2FC_ATAC & padj_GF_mLN_pLN <= padj) |
                                 (abs(log2FC_GF_mLN_pLN) >= log2FC_ATAC & padj_GF_mLN_pLN <= padj) )
DARs_features_GF_SPF_NA <- DARs_features_GF_SPF[!is.na(DARs_features_GF_SPF$symbol),]
#Heatmap
data_heatmap <- DARs_features_GF_SPF_NA[,c("id",
                                           "ATAC_mLNSPF1","ATAC_mLNSPF2","ATAC_mLNSPF3","ATAC_mLNSPF4",
                                           "ATAC_pLNSPF1","ATAC_pLNSPF2","ATAC_pLNSPF3",
                                           "ATAC_mLNGF1","ATAC_mLNGF2","ATAC_mLNGF3",
                                           "ATAC_pLNGF1","ATAC_pLNGF3")]
rownames(data_heatmap) <- data_heatmap$id
data_heatmap <- as.matrix(data_heatmap[,2:ncol(data_heatmap)])
data_heatmap[data_heatmap == 0] <- 0.4
min(data_heatmap)
data_heatmap <- log2(data_heatmap)

min(data_heatmap)
pheatmap(data_heatmap, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = FALSE, cluster_cols = TRUE, cellwidth = 10,
         scale = "row", border_color = "black", color = colorRampPalette(c("white","grey", "darkgreen"), space="rgb")(128),
         main = "DARs TSS associated SPF & GF")




#####
#Meaned GF & SPF
#####
DARs_features_GF_SPF <- subset(DARs_features, (abs(log2FC_GF_mLN_pLN) >= log2FC_ATAC & padj_GF_mLN_pLN <= padj) |
                                 (abs(log2FC_SPF_mLN_pLN) >= log2FC_ATAC & padj_SPF_mLN_pLN <= padj) )
DARs_features_GF_SPF_NA <- DARs_features_GF_SPF[!is.na(DARs_features_GF_SPF$symbol),]
#DARs_features_GF_SPF_NA_Prom <- subset(DARs_features_GF_SPF_NA, distancetoFeature >= -400 & distancetoFeature <= 10000)
DARs_features_GF_SPF_NA_unique <- DARs_features_GF_SPF_NA[!duplicated(DARs_features_GF_SPF_NA[,c("symbol")]),]


#Heatmap
data_heatmap <- DARs_features_GF_SPF_NA_unique[,c("symbol",
                                                       "ATAC_mLNSPF1","ATAC_mLNSPF2","ATAC_mLNSPF3","ATAC_mLNSPF4",
                                                       "ATAC_pLNSPF1","ATAC_pLNSPF2","ATAC_pLNSPF3",
                                                       "ATAC_mLNGF1","ATAC_mLNGF2","ATAC_mLNGF3",
                                                       "ATAC_pLNGF1","ATAC_pLNGF3")]

rownames(data_heatmap) <- data_heatmap$symbol
data_heatmap <- as.matrix(data_heatmap[,2:ncol(data_heatmap)])
mLN_SPF <- rowMeans(data_heatmap[,1:4])
pLN_SPF <- rowMeans(data_heatmap[,5:7])
mLN_GF <- rowMeans(data_heatmap[,8:10])
pLN_GF <- rowMeans(data_heatmap[,11:12])
data_heatmap <- cbind(mLN_SPF,pLN_SPF,mLN_GF,pLN_GF)
#data_heatmap[data_heatmap == 0] <- 0.4
min(data_heatmap)
data_heatmap <- log2(data_heatmap)
data_heatmap[data_heatmap == -Inf] <- -2
min(data_heatmap)

pheatmap(data_heatmap, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 20, show_rownames = FALSE, cluster_cols = TRUE, cellwidth = 10,
         scale = "none", border_color = "black", color = colorRampPalette(c("black","dimgrey","gray48","grey","honeydew2","white"), space="rgb")(128),
         main = "DARs TSS associated SPF")
outs <- pheatmap(data_heatmap, cluster_rows = TRUE, legend = TRUE,
                 treeheight_row = 0, treeheight_col = 20, show_rownames = TRUE, cluster_cols = TRUE, cellwidth = 10,
                 scale = "none", border_color = "black", color = colorRampPalette(c("black","dimgrey","gray48","grey","honeydew2","white"), space="rgb")(128),
                 main = "DARs TSS associated SPF")

#####
#Volcano plots 
#####
# Note: Figures for Publication
DARs_volcano <- DARs_features
DARs_volcano$log10_padj_pLN_SPF_GF <- -log10(DARs_volcano$padj_pLN_SPF_GF)
DARs_volcano$log10_padj_mLN_SPF_GF <- -log10(DARs_volcano$padj_mLN_SPF_GF)
DARs_volcano$log10_padj_SPF_mLN_pLN <- -log10(DARs_volcano$padj_SPF_mLN_pLN)
DARs_volcano$log10_padj_GF_mLN_pLN <- -log10(DARs_volcano$padj_GF_mLN_pLN)

# mLN SPF / pLN SPF ------------------------------
#Count and order
# Left
Left <- as.numeric(DARs_volcano$log2FC_SPF_mLN_pLN < -1.0 &
                     DARs_volcano$padj_SPF_mLN_pLN < 0.05)
Left_n_Genes <- unique(subset(DARs_volcano, 
                                DARs_volcano$log2FC_SPF_mLN_pLN < -1.0 &
                                DARs_volcano$padj_SPF_mLN_pLN < 0.05)$symbol)
Left_n_Regions <- sum(Left)

# Left
Right <- as.numeric(DARs_volcano$log2FC_SPF_mLN_pLN > 1.0 &
                     DARs_volcano$padj_SPF_mLN_pLN < 0.05)
Right_n_Genes <- unique(subset(DARs_volcano, 
                              DARs_volcano$log2FC_SPF_mLN_pLN > 1.0 &
                              DARs_volcano$padj_SPF_mLN_pLN < 0.05)$symbol)
Right_n_Regions <- sum(Right)

Right <- replace(Right, Right == 1, 2)

#Categorize for coloring
ordering <- Left + Right
data_DARs_volcano <- cbind(DARs_volcano, Left, Right, ordering)
data_DARs_volcano <- arrange(data_DARs_volcano, ordering)


#Numbers of genes
Right_n <- length(Right_n_Genes)
Left_n <- length(Left_n_Genes)

#Plot
d <- ggplot(data_DARs_volcano, aes(log2FC_SPF_mLN_pLN, log10_padj_SPF_mLN_pLN)) #, colour = factor(ordering))
d + geom_point(aes(colour = factor(ordering)), size = 1.0, show.legend = FALSE) +
  scale_color_manual(values = c("black","red", "blue")) +
  theme_classic(base_size = 12, base_family = "") +
  #Lab
  xlab("log2(FC) mLNSPF/pLNSPF") +
  ylab("-log10(padj) mLNSPF/pLNSPF") +
  
  #Line
  geom_hline(yintercept=1, colour = "black", linetype="dashed") +
  geom_hline(yintercept=-1, colour = "black", linetype="dashed") +
  geom_vline(xintercept=1, colour = "black", linetype="dashed") +
  geom_vline(xintercept=-1, colour = "black", linetype="dashed") +
  
  #Numbers
  geom_text(data = data.frame(), aes(-5, 18, label = Left_n), colour = "darkred", size = 6) + 
  geom_text(data = data.frame(), aes(-5, 20, label = Left_n_Regions), colour = "red", size = 6) + 
  geom_text(data = data.frame(), aes(5, 18, label = Right_n), colour = "navyblue", size = 6) +
  geom_text(data = data.frame(), aes(5, 20, label = Right_n_Regions), colour = "blue", size = 6) +

  #Genenames
  #Range
  scale_x_continuous(limits = c(-6, 6)) +
  scale_y_continuous(limits = c(0, 20))

# mLN SPF / mLN GF ------------------------------
#Count and order
# Left
Left <- as.numeric(DARs_volcano$log2FC_mLN_SPF_GF < -1.0 &
                     DARs_volcano$padj_mLN_SPF_GF < 0.05)
Left[is.na(Left)] <- 0
Left_n_Genes <- na.omit(unique(subset(DARs_volcano, 
                              DARs_volcano$log2FC_mLN_SPF_GF < -1.0 &
                              DARs_volcano$padj_mLN_SPF_GF < 0.05)$symbol))

Left_n_Regions <- sum(Left)

# Left
Right <- as.numeric(DARs_volcano$log2FC_mLN_SPF_GF > 1.0 &
                      DARs_volcano$padj_mLN_SPF_GF < 0.05)
Right[is.na(Right)] <- 0

Right_n_Genes <- na.omit(unique(subset(DARs_volcano, 
                               DARs_volcano$log2FC_mLN_SPF_GF > 1.0 &
                                 DARs_volcano$padj_mLN_SPF_GF < 0.05)$symbol))
Right_n_Regions <- sum(Right)

Right <- replace(Right, Right == 1, 2)

#Categorize for coloring
ordering <- Left + Right
data_DARs_volcano <- cbind(DARs_volcano, Left, Right, ordering)
data_DARs_volcano <- arrange(data_DARs_volcano, ordering)


#Numbers of genes
Right_n <- length(Right_n_Genes)
Left_n <- length(Left_n_Genes)

#Plot
d <- ggplot(data_DARs_volcano, aes(log2FC_mLN_SPF_GF, log10_padj_mLN_SPF_GF)) #, colour = factor(ordering))
d + geom_point(aes(colour = factor(ordering)), size = 1.0, show.legend = FALSE) +
  scale_color_manual(values = c("black","red", "blue")) +
  theme_classic(base_size = 12, base_family = "") +
  #Lab
  xlab("log2(FC) mLNSPF/mLNGF") +
  ylab("-log10(padj) mLNGF/mLNGF") +
  
  #Line
  geom_hline(yintercept=1, colour = "black", linetype="dashed") +
  geom_hline(yintercept=-1, colour = "black", linetype="dashed") +
  geom_vline(xintercept=1, colour = "black", linetype="dashed") +
  geom_vline(xintercept=-1, colour = "black", linetype="dashed") +
  
  #Numbers
  geom_text(data = data.frame(), aes(-5, 18, label = Left_n), colour = "darkred", size = 6) + 
  geom_text(data = data.frame(), aes(-5, 20, label = Left_n_Regions), colour = "red", size = 6) + 
  geom_text(data = data.frame(), aes(5, 18, label = Right_n), colour = "navyblue", size = 6) +
  geom_text(data = data.frame(), aes(5, 20, label = Right_n_Regions), colour = "blue", size = 6) +
  
  #Genenames
  #Range
  scale_x_continuous(limits = c(-6, 6)) +
  scale_y_continuous(limits = c(0, 20))

# pLN SPF / pLN GF ------------------------------
#Count and order
# Left
Left <- as.numeric(DARs_volcano$log2FC_pLN_SPF_GF < -1.0 &
                     DARs_volcano$padj_pLN_SPF_GF < 0.05)
Left[is.na(Left)] <- 0
Left_n_Genes <- na.omit(unique(subset(DARs_volcano, 
                                      DARs_volcano$log2FC_pLN_SPF_GF < -1.0 &
                                        DARs_volcano$padj_pLN_SPF_GF < 0.05)$symbol))

Left_n_Regions <- sum(Left)

# Left
Right <- as.numeric(DARs_volcano$log2FC_pLN_SPF_GF > 1.0 &
                      DARs_volcano$padj_pLN_SPF_GF < 0.05)
Right[is.na(Right)] <- 0

Right_n_Genes <- na.omit(unique(subset(DARs_volcano, 
                                       DARs_volcano$log2FC_pLN_SPF_GF > 1.0 &
                                         DARs_volcano$padj_pLN_SPF_GF < 0.05)$symbol))
Right_n_Regions <- sum(Right)

Right <- replace(Right, Right == 1, 2)

#Categorize for coloring
ordering <- Left + Right
data_DARs_volcano <- cbind(DARs_volcano, Left, Right, ordering)
data_DARs_volcano <- arrange(data_DARs_volcano, ordering)


#Numbers of genes
Right_n <- length(Right_n_Genes)
Left_n <- length(Left_n_Genes)

#Plot
d <- ggplot(data_DARs_volcano, aes(log2FC_pLN_SPF_GF, log10_padj_pLN_SPF_GF)) #, colour = factor(ordering))
d + geom_point(aes(colour = factor(ordering)), size = 1.0, show.legend = FALSE) +
  scale_color_manual(values = c("black","red", "blue")) +
  theme_classic(base_size = 12, base_family = "") +
  #Lab
  xlab("log2(FC) pLNSPF/pLNGF") +
  ylab("-log10(padj) pLNSPF/pLNGF") +
  
  #Line
  geom_hline(yintercept=1, colour = "black", linetype="dashed") +
  geom_hline(yintercept=-1, colour = "black", linetype="dashed") +
  geom_vline(xintercept=1, colour = "black", linetype="dashed") +
  geom_vline(xintercept=-1, colour = "black", linetype="dashed") +
  
  #Numbers
  geom_text(data = data.frame(), aes(-5, 18, label = Left_n), colour = "darkred", size = 6) + 
  geom_text(data = data.frame(), aes(-5, 20, label = Left_n_Regions), colour = "red", size = 6) + 
  geom_text(data = data.frame(), aes(5, 18, label = Right_n), colour = "navyblue", size = 6) +
  geom_text(data = data.frame(), aes(5, 20, label = Right_n_Regions), colour = "blue", size = 6) +
  
  #Genenames
  #Range
  scale_x_continuous(limits = c(-6, 6)) +
  scale_y_continuous(limits = c(0, 20))


####################
#RNAseq analysis
####################
#Load data
setwd(path_RNAseq_DESeq2)
names = list.files(pattern="*.csv")
myfiles = lapply(names, read.csv)

#merrge Pair-wise comparisons
names <- unlist(strsplit(names, split='DESeq2_genes_diffexp_by_fc_', fixed=TRUE))
names <- unlist(strsplit(names, split='.csv', fixed=TRUE))
#vector of length 2 x
names_group <- unlist(strsplit(names, split='_vs_', fixed=TRUE))
#take even and uneven elements
#uneven
names_group_1 <- names_group[seq(1, length(names_group), 2)]
#even
names_group_2 <- names_group[seq(2, length(names_group), 2)]
#new names with propper divivsion order
names <- paste(names_group_2, "_VS_", names_group_1, sep = "")
#get relevant columns from list of tables (RPKM, FC, padj)
#assign list names
names(myfiles) <- names

#get list name and attach to column names
list_all <- lapply(seq_along(myfiles),function(i){
  name_i <- names(myfiles[i])
  #colnames(myfiles[[i]])
  table_i <- cbind(myfiles[[i]][,c("GeneSymbol","log2FoldChange","padj")],myfiles[[i]][,grepl("RPKMcounts", colnames(myfiles[[i]]))])
  colnames(table_i)[2] <- paste("log2FC_", name_i, sep = "")
  colnames(table_i)[3] <- paste("padj_", name_i, sep = "")
  table_i
}
)

#make table from list
DF <- list_all[[1]]
for (.df in list_all) {
  DF <-merge(DF,.df,by = "GeneSymbol", all=T)
  DF <- DF[!duplicated(DF$GeneSymbol),]
}
data_all <- cbind(DF$GeneSymbol, 
                  #DF[,grepl(c(".x"), colnames(DF))],
                  DF[,grepl(c("padj"), colnames(DF))],
                  DF[,grepl(c("log2FC"), colnames(DF))],
                  DF[,c("RPKMcounts.mLN_GF_FRC_1.x","RPKMcounts.mLN_GF_FRC_2.x","RPKMcounts.mLN_GF_FRC_3.x",
                        "RPKMcounts.pLN_GF_FRC_1.x","RPKMcounts.pLN_GF_FRC_2.x","RPKMcounts.pLN_GF_FRC_3.x",
                        "RPKMcounts.mLN_SPF_FRC_1.x","RPKMcounts.mLN_SPF_FRC_2.x","RPKMcounts.mLN_SPF_FRC_3.x",
                        "RPKMcounts.pLN_SFP_FRC_1.x","RPKMcounts.pLN_SFP_FRC_2.x","RPKMcounts.pLN_SFP_FRC_3.x")])
colnames(data_all)[1] <- "GeneSymbol"
colnames(data_all)[(ncol(data_all)-11):ncol(data_all)] <- c("mLN_GF_1","mLN_GF_2","mLN_GF_3",
                                                            "pLN_GF_1","pLN_GF_2","pLN_GF_3",
                                                            "mLN_SPF_1","mLN_SPF_2","mLN_SPF_3",
                                                            "pLN_SFP_1","pLN_SFP_2","pLN_SFP_3")
#Eliminate columns
data_all <- data_all[,c(1,2,4:7,9:ncol(data_all))]
data_all <- data_all[!(is.na(data_all$GeneSymbol)),]
data_all <- data_all[!duplicated(data_all$GeneSymbol),]
#reorder
data_all <- data_all[,c(1,6:9,2:5,10:ncol(data_all))]
colnames(data_all)[1:9] <- c("GeneSymbol",
                        "GF_mLN_pLN_log2FoldChange","mLN_SPF_GF_log2FoldChange",
                        "SPF_mLN_pLN_log2FoldChange","pLN_SPF_GF_log2FoldChange",
                        "GF_mLN_pLN_padj","mLN_SPF_GF_padj",
                        "SPF_mLN_pLN_padj","pLN_SPF_GF_padj")
#Invert FC to align to renaming of columns
data_all_test <- data_all[,2:6] * (-1)

#store
expression <- data_all
expression_NA <- na.omit(expression)

#####
#FC SPF plot DAR and DEG
#####
# Note: Figure for paper
# Check consistency of DARs per TSS region for significant ones
DARs_features_NA <- DARs_features[!is.na(DARs_features$symbol),]
# Average number of DARs per TSS
mean(unlist(lapply(l_DARs_features_NA, nrow)))
# Take all Genes that have a significant DAR
DARs_features_NA_sig <- subset(DARs_features_NA, abs(log2FC_SPF_mLN_pLN) >= 1.0 & padj_SPF_mLN_pLN <= 0.05)
print(paste("Number of significant DARs: ", nrow(DARs_features_NA_sig)))
sig_Genes_DARs <- unique(DARs_features_NA_sig$symbol)
print(paste("Number of significant Genes with a DAR: ", length(sig_Genes_DARs)))
DARs_features_NA_sig_plusX <- subset(DARs_features_NA, symbol %in% sig_Genes_DARs)
print(paste("Number of DARs for Genes that contain at least one significant DAR: ", nrow(DARs_features_NA_sig_plusX)))


#Merge DARs and DEGs
DARs_DEGs <- merge(DARs_features_NA, expression, by.x = "symbol", by.y = "GeneSymbol")
plot(DARs_DEGs$log2FC_SPF_mLN_pLN, DARs_DEGs$SPF_mLN_pLN_log2FoldChange)
#Change orientation of FC for expression
DARs_DEGs$SPF_mLN_pLN_log2FoldChange <- DARs_DEGs$SPF_mLN_pLN_log2FoldChange * (-1)
plot(DARs_DEGs$log2FC_SPF_mLN_pLN, DARs_DEGs$SPF_mLN_pLN_log2FoldChange)

#Set Quadrants
#mLN open and UP - DEG regulated by openess
right_up <- as.numeric(DARs_DEGs$log2FC_SPF_mLN_pLN > 1.0 &
                         DARs_DEGs$SPF_mLN_pLN_log2FoldChange > 1.0 &
                         (DARs_DEGs$padj_SPF_mLN_pLN < 0.05 |
                            DARs_DEGs$SPF_mLN_pLN_padj < 0.05))
right_up_n_Genes <- unique(subset(DARs_DEGs, DARs_DEGs$log2FC_SPF_mLN_pLN > 1.0 &
                                    DARs_DEGs$SPF_mLN_pLN_log2FoldChange > 1.0 &
                                    (DARs_DEGs$padj_SPF_mLN_pLN < 0.05 |
                                       DARs_DEGs$SPF_mLN_pLN_padj < 0.05))$symbol)

#pLN open and UP - - DEG regulated by openess
left_down <- as.numeric(DARs_DEGs$log2FC_SPF_mLN_pLN < -1.0 &
                          DARs_DEGs$SPF_mLN_pLN_log2FoldChange < -1.0 &
                          (DARs_DEGs$padj_SPF_mLN_pLN < 0.05 |
                             DARs_DEGs$SPF_mLN_pLN_padj < 0.05))
left_down_n_Genes <- unique(subset(DARs_DEGs, DARs_DEGs$log2FC_SPF_mLN_pLN < -1.0 &
                                     DARs_DEGs$SPF_mLN_pLN_log2FoldChange < -1.0 &
                                     (DARs_DEGs$padj_SPF_mLN_pLN < 0.05 |
                                        DARs_DEGs$SPF_mLN_pLN_padj < 0.05))$symbol)
left_down <- replace(left_down, left_down == 1, 2)

#mLN open no DEG - activatable
right <- as.numeric(DARs_DEGs$log2FC_SPF_mLN_pLN > 1.0 &
                      DARs_DEGs$SPF_mLN_pLN_log2FoldChange > -1.0 &
                      DARs_DEGs$SPF_mLN_pLN_log2FoldChange < 1.0 &
                      DARs_DEGs$padj_SPF_mLN_pLN < 0.05)
right_n_Genes <- unique(subset(DARs_DEGs, DARs_DEGs$log2FC_SPF_mLN_pLN > 1.0 &
                                 DARs_DEGs$SPF_mLN_pLN_log2FoldChange > -1.0 &
                                 DARs_DEGs$SPF_mLN_pLN_log2FoldChange < 1.0 &
                                 DARs_DEGs$padj_SPF_mLN_pLN < 0.05)$symbol)
right <- replace(right, right == 1, 3)

#pLN open no DEG - activatable
left <- as.numeric(DARs_DEGs$log2FC_SPF_mLN_pLN < -1.0 &
                      DARs_DEGs$SPF_mLN_pLN_log2FoldChange > -1.0 &
                      DARs_DEGs$SPF_mLN_pLN_log2FoldChange < 1.0 &
                      DARs_DEGs$padj_SPF_mLN_pLN < 0.05)
left_n_Genes <- unique(subset(DARs_DEGs, DARs_DEGs$log2FC_SPF_mLN_pLN < -1.0 &
                                 DARs_DEGs$SPF_mLN_pLN_log2FoldChange > -1.0 &
                                 DARs_DEGs$SPF_mLN_pLN_log2FoldChange < 1.0 &
                                 DARs_DEGs$padj_SPF_mLN_pLN < 0.05)$symbol)
left <- replace(left, left == 1, 4)

#mLN noDAR but UP DEG - TF regulated
top <- as.numeric(DARs_DEGs$SPF_mLN_pLN_log2FoldChange > 1.0 &
                      DARs_DEGs$log2FC_SPF_mLN_pLN > -1.0 &
                      DARs_DEGs$log2FC_SPF_mLN_pLN < 1.0 &
                      DARs_DEGs$SPF_mLN_pLN_padj < 0.05)
top_n_Genes <- unique(subset(DARs_DEGs,DARs_DEGs$SPF_mLN_pLN_log2FoldChange > 1.0 &
                                 DARs_DEGs$log2FC_SPF_mLN_pLN > -1.0 &
                                 DARs_DEGs$log2FC_SPF_mLN_pLN < 1.0 &
                                 DARs_DEGs$SPF_mLN_pLN_padj < 0.05)$symbol)
top <- replace(top, top == 1, 5)

#pLN noDAR but UP DEG - TF regulated
bottom <- as.numeric(DARs_DEGs$SPF_mLN_pLN_log2FoldChange < -1.0 &
                    DARs_DEGs$log2FC_SPF_mLN_pLN > -1.0 &
                    DARs_DEGs$log2FC_SPF_mLN_pLN < 1.0 &
                    DARs_DEGs$SPF_mLN_pLN_padj < 0.05)
bottom_n_Genes <- unique(subset(DARs_DEGs,DARs_DEGs$SPF_mLN_pLN_log2FoldChange < -1.0 &
                               DARs_DEGs$log2FC_SPF_mLN_pLN > -1.0 &
                               DARs_DEGs$log2FC_SPF_mLN_pLN < 1.0 &
                               DARs_DEGs$SPF_mLN_pLN_padj < 0.05)$symbol)
bottom <- replace(bottom, bottom == 1, 6)


#Categorize for coloring
ordering <- right_up + left_down + right + left + top + bottom
data_DEGs_DARs <- cbind(DARs_DEGs, right_up, left_down, right, left, top, bottom, ordering)
data_DEGs_DARs <- arrange(data_DEGs_DARs, ordering)


#Numbers of genes
right_up_n <- length(right_up_n_Genes)
left_down_n <- length(left_down_n_Genes)
left_n <- length(left_n_Genes)
right_n <- length(right_n_Genes)
top_n <- length(top_n_Genes)
bottom_n <- length(bottom_n_Genes)

#Plot
d <- ggplot(data_DEGs_DARs, aes(log2FC_SPF_mLN_pLN, SPF_mLN_pLN_log2FoldChange)) #, colour = factor(ordering))
d + geom_point(aes(colour = factor(ordering)), size = 2.0, show_guide = FALSE) +
  scale_color_manual(values = c("black","gray20", "gray40", "#7B337E", "#9C70AA","#206396","#5099BB")) +
  #Line
  geom_hline(yintercept=1, colour = "black", linetype="dashed") +
  geom_hline(yintercept=-1, colour = "black", linetype="dashed") +
  geom_vline(xintercept=1, colour = "black", linetype="dashed") +
  geom_vline(xintercept=-1, colour = "black", linetype="dashed") +
  theme_classic(base_size = 12, base_family = "") +
  
  #Numbers
  geom_text(data = data.frame(), aes(5, 8, label = right_up_n), colour = "gray20", size = 6) + 
  geom_text(data = data.frame(), aes(-3, -8, label = left_down_n), colour = "gray40", size = 6) +
  geom_text(data = data.frame(), aes(-3, 0, label = left_n), colour = "#9C70AA", size = 6) + 
  geom_text(data = data.frame(), aes(5, 0, label = right_n), colour = "#7B337E", size = 6) +
  geom_text(data = data.frame(), aes(0, 8, label = top_n), colour = "#206396", size = 6) + 
  geom_text(data = data.frame(), aes(0, -8, label = bottom_n), colour = "#5099BB", size = 6) +
  
  #Lab
  xlab("ATAC log(FC) mLN/pLN") +
  ylab("RNAseq log(FC) mLN/pLN") +
  
  #Genenames
  #Range
  scale_x_continuous(limits = c(-3, 5)) +
  scale_y_continuous(limits = c(-8, 8))

#####
#FC GF plot DAR and DEG
#####
# Check consistency of DARs per TSS region for significant ones
DARs_features_NA <- DARs_features[!is.na(DARs_features$symbol),]
# Average number of DARs per TSS
#mean(unlist(lapply(l_DARs_features_NA, nrow)))
# Take all Genes that have a significant DAR
DARs_features_NA_sig <- subset(DARs_features_NA, abs(log2FC_GF_mLN_pLN) >= 1.0 & padj_GF_mLN_pLN <= 0.05)
print(paste("Number of significant DARs: ", nrow(DARs_features_NA_sig)))
sig_Genes_DARs <- unique(DARs_features_NA_sig$symbol)
print(paste("Number of significant Genes with a DAR: ", length(sig_Genes_DARs)))
DARs_features_NA_sig_plusX <- subset(DARs_features_NA, symbol %in% sig_Genes_DARs)
print(paste("Number of DARs for Genes that contain at least one significant DAR: ", nrow(DARs_features_NA_sig_plusX)))
# mean FC and SD per Gene
#mean_FC_per_DAR <- lapply(l_DARs_features_NA, function(x){
#  FC_mean <- mean(x$log2FC_GF_mLN_pLN)
#  FC_mean
#})
#hist(unlist(mean_FC_per_DAR))

#Merge DARs and DEGs
DARs_DEGs <- merge(DARs_features_NA, expression, by.x = "symbol", by.y = "GeneSymbol")
plot(DARs_DEGs$log2FC_GF_mLN_pLN, DARs_DEGs$GF_mLN_pLN_log2FoldChange)
#Change orientation of FC for expression
DARs_DEGs$GF_mLN_pLN_log2FoldChange <- DARs_DEGs$GF_mLN_pLN_log2FoldChange * (-1)
plot(DARs_DEGs$log2FC_GF_mLN_pLN, DARs_DEGs$GF_mLN_pLN_log2FoldChange)

#Set Quadrants

#mLN open and UP - DEG regulated by openess
right_up <- as.numeric(DARs_DEGs$log2FC_GF_mLN_pLN > 1.0 &
                         DARs_DEGs$GF_mLN_pLN_log2FoldChange > 1.0 &
                         (DARs_DEGs$padj_GF_mLN_pLN < 0.05 |
                            DARs_DEGs$GF_mLN_pLN_padj < 0.05))
right_up_n_Genes <- unique(subset(DARs_DEGs, DARs_DEGs$log2FC_GF_mLN_pLN > 1.0 &
                                    DARs_DEGs$GF_mLN_pLN_log2FoldChange > 1.0 &
                                    (DARs_DEGs$padj_GF_mLN_pLN < 0.05 |
                                       DARs_DEGs$GF_mLN_pLN_padj < 0.05))$symbol)

#pLN open and UP - - DEG regulated by openess
left_down <- as.numeric(DARs_DEGs$log2FC_GF_mLN_pLN < -1.0 &
                          DARs_DEGs$GF_mLN_pLN_log2FoldChange < -1.0 &
                          (DARs_DEGs$padj_GF_mLN_pLN < 0.05 |
                             DARs_DEGs$GF_mLN_pLN_padj < 0.05))
left_down_n_Genes <- unique(subset(DARs_DEGs, DARs_DEGs$log2FC_GF_mLN_pLN < -1.0 &
                                     DARs_DEGs$GF_mLN_pLN_log2FoldChange < -1.0 &
                                     (DARs_DEGs$padj_GF_mLN_pLN < 0.05 |
                                        DARs_DEGs$GF_mLN_pLN_padj < 0.05))$symbol)
left_down <- replace(left_down, left_down == 1, 2)

#mLN open no DEG - activatable
right <- as.numeric(DARs_DEGs$log2FC_GF_mLN_pLN > 1.0 &
                      DARs_DEGs$GF_mLN_pLN_log2FoldChange > -1.0 &
                      DARs_DEGs$GF_mLN_pLN_log2FoldChange < 1.0 &
                      DARs_DEGs$padj_GF_mLN_pLN < 0.05)
right_n_Genes <- unique(subset(DARs_DEGs, DARs_DEGs$log2FC_GF_mLN_pLN > 1.0 &
                                 DARs_DEGs$GF_mLN_pLN_log2FoldChange > -1.0 &
                                 DARs_DEGs$GF_mLN_pLN_log2FoldChange < 1.0 &
                                 DARs_DEGs$padj_GF_mLN_pLN < 0.05)$symbol)
right <- replace(right, right == 1, 3)

#pLN open no DEG - activatable
left <- as.numeric(DARs_DEGs$log2FC_GF_mLN_pLN < -1.0 &
                     DARs_DEGs$GF_mLN_pLN_log2FoldChange > -1.0 &
                     DARs_DEGs$GF_mLN_pLN_log2FoldChange < 1.0 &
                     DARs_DEGs$padj_GF_mLN_pLN < 0.05)
left_n_Genes <- unique(subset(DARs_DEGs, DARs_DEGs$log2FC_GF_mLN_pLN < -1.0 &
                                DARs_DEGs$GF_mLN_pLN_log2FoldChange > -1.0 &
                                DARs_DEGs$GF_mLN_pLN_log2FoldChange < 1.0 &
                                DARs_DEGs$padj_GF_mLN_pLN < 0.05)$symbol)
left <- replace(left, left == 1, 4)

#mLN noDAR but UP DEG - TF regulated
top <- as.numeric(DARs_DEGs$GF_mLN_pLN_log2FoldChange > 1.0 &
                    DARs_DEGs$log2FC_GF_mLN_pLN > -1.0 &
                    DARs_DEGs$log2FC_GF_mLN_pLN < 1.0 &
                    DARs_DEGs$GF_mLN_pLN_padj < 0.05)
top_n_Genes <- unique(subset(DARs_DEGs,DARs_DEGs$GF_mLN_pLN_log2FoldChange > 1.0 &
                               DARs_DEGs$log2FC_GF_mLN_pLN > -1.0 &
                               DARs_DEGs$log2FC_GF_mLN_pLN < 1.0 &
                               DARs_DEGs$GF_mLN_pLN_padj < 0.05)$symbol)
top <- replace(top, top == 1, 5)

#pLN noDAR but UP DEG - TF regulated
bottom <- as.numeric(DARs_DEGs$GF_mLN_pLN_log2FoldChange < -1.0 &
                       DARs_DEGs$log2FC_GF_mLN_pLN > -1.0 &
                       DARs_DEGs$log2FC_GF_mLN_pLN < 1.0 &
                       DARs_DEGs$GF_mLN_pLN_padj < 0.05)
bottom_n_Genes <- unique(subset(DARs_DEGs,DARs_DEGs$GF_mLN_pLN_log2FoldChange < -1.0 &
                                  DARs_DEGs$log2FC_GF_mLN_pLN > -1.0 &
                                  DARs_DEGs$log2FC_GF_mLN_pLN < 1.0 &
                                  DARs_DEGs$GF_mLN_pLN_padj < 0.05)$symbol)
bottom <- replace(bottom, bottom == 1, 6)


#Categorize for coloring
ordering <- right_up + left_down + right + left + top + bottom
data_DEGs_DARs <- cbind(DARs_DEGs, right_up, left_down, right, left, top, bottom, ordering)
data_DEGs_DARs <- arrange(data_DEGs_DARs, ordering)


#Numbers of genes
right_up_n <- length(right_up_n_Genes)
left_down_n <- length(left_down_n_Genes)
left_n <- length(left_n_Genes)
right_n <- length(right_n_Genes)
top_n <- length(top_n_Genes)
bottom_n <- length(bottom_n_Genes)

#Plot
d <- ggplot(data_DEGs_DARs, aes(log2FC_GF_mLN_pLN, GF_mLN_pLN_log2FoldChange)) #, colour = factor(ordering))
d + geom_point(aes(colour = factor(ordering)), size = 2.0, show_guide = FALSE) +
  scale_color_manual(values = c("black","gray48", "gray32", "deepskyblue", "deeppink3","deepskyblue4","deeppink4")) +
  #Line
  geom_hline(yintercept=1, colour = "black", linetype="dashed") +
  geom_hline(yintercept=-1, colour = "black", linetype="dashed") +
  geom_vline(xintercept=1, colour = "black", linetype="dashed") +
  geom_vline(xintercept=-1, colour = "black", linetype="dashed") +
  theme_classic(base_size = 20, base_family = "") +
  
  #Numbers
  geom_text(data = data.frame(), aes(4, 4, label = right_up_n), colour = "gray48", size = 8) + 
  geom_text(data = data.frame(), aes(-4, -4, label = left_down_n), colour = "gray32", size = 8) +
  geom_text(data = data.frame(), aes(-4, 0, label = left_n), colour = "deeppink3", size = 8) + 
  geom_text(data = data.frame(), aes(4, 0, label = right_n), colour = "deepskyblue", size = 8) +
  geom_text(data = data.frame(), aes(0, 4, label = top_n), colour = "deepskyblue4", size = 8) + 
  geom_text(data = data.frame(), aes(0, -4, label = bottom_n), colour = "deeppink4", size = 8) +
  
  #Lab
  xlab("ATAC log(FC) mLN/pLN") +
  ylab("RNAseq log(FC) mLN/pLN") +
  
  #Genenames
  #Range
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(-5, 5))

#####
#Compare RNA and ATAC
#####
#Comparisons of  all
#Prep RNAseq tables
#for the four relevant comparisons
#generate a list that contains vectors with the UP and DOWN regulated genes

RNAseq_UP <- list()
RNAseq_DOWN <- list()
print("Number of DEGs per comparison")
for(i in 1:4){
  index_log2FC = i + 1
  index_padj = i + 5
  print(colnames(expression[index_log2FC]))
  UP <- as.character(expression[(expression[,index_log2FC] <= -log2FC_RNA & !is.na(expression[,index_log2FC])) & 
                                  (expression[,index_padj] <= padj & !is.na(expression[,index_padj])),]$GeneSymbol)
  print(length(UP))
  DOWN <- as.character(expression[(expression[,index_log2FC] >= log2FC_RNA & !is.na(expression[,index_log2FC])) & 
                                    (expression[,index_padj] <= padj & !is.na(expression[,index_padj])),]$GeneSymbol)
  print(length(DOWN))
  RNAseq_UP[[i]] <- UP
  RNAseq_DOWN[[i]] <- DOWN
}

#Name list elements for later utilization
names(RNAseq_UP) <- c("SPF_mLN_pLN","GF_mLN_pLN",
                      "mLN_SPF_GF","pLN_SPF_GF")
names(RNAseq_DOWN) <- c("SPF_mLN_pLN","GF_mLN_pLN",
                        "mLN_SPF_GF","pLN_SPF_GF")

#Prep ATACseq tables
#for the four relevant comparisons
#generate a list that contains vectors with the UP and DOWN regulated genes
DARs_features_NA <- DARs_features[!is.na(DARs_features$symbol),]
ATACseq_UP <- list()
ATACseq_DOWN <- list()
k = 5
for(i in 1:4){
  k = k + 2 
  index_log2FC = k
  print(colnames(DARs_features_NA[index_log2FC]))
  index_padj = k + 1
  UP <- as.character(DARs_features_NA[(DARs_features_NA[,index_log2FC] >= log2FC_ATAC & !is.na(DARs_features_NA[,index_log2FC])) & 
                                        (DARs_features_NA[,index_padj] <= padj & !is.na(DARs_features_NA[,index_padj])),]$symbol)
  print(length(UP))
  #print(UP)
  DOWN <- as.character(DARs_features_NA[(DARs_features_NA[,index_log2FC] <= -log2FC_ATAC & !is.na(DARs_features_NA[,index_log2FC])) & 
                                          (DARs_features_NA[,index_padj] <= padj & !is.na(DARs_features_NA[,index_padj])),]$symbol)
  print(length(DOWN))
  #print(DOWN)
  ATACseq_UP[[i]] <- UP
  ATACseq_DOWN[[i]] <- DOWN
}
#Name list elements for later utilization
names(ATACseq_UP) <- c("mLN_SPF_GF","pLN_SPF_GF",
                       "SPF_mLN_pLN","GF_mLN_pLN")
ATACseq_UP <- ATACseq_UP[c("SPF_mLN_pLN","GF_mLN_pLN",
                           "mLN_SPF_GF","pLN_SPF_GF")]

names(ATACseq_DOWN) <- c("mLN_SPF_GF","pLN_SPF_GF",
                         "SPF_mLN_pLN","GF_mLN_pLN")
ATACseq_DOWN <- ATACseq_DOWN[c("SPF_mLN_pLN","GF_mLN_pLN",
                               "mLN_SPF_GF","pLN_SPF_GF")]
#reorder lists according to RNAseq
ATACseq_UP <- ATACseq_UP[names(RNAseq_UP)]
ATACseq_DOWN <- ATACseq_DOWN[names(RNAseq_DOWN)]

#Check the overlaps
#empty table for storage
summary_comparison = data.frame(matrix("", ncol = 8, nrow = 0))
#empty lists for storage of GeneNames
Plus_Plus <- list()
Minus_Minus <- list()
Plus_Zero <- list()
Zero_Plus <- list()
Minus_Zero <- list()
Zero_Minus <- list()
Present_Present <- list()

for(i in 1:4){
  #Grab corresponding Gene lists
  RNA_UP_i <- RNAseq_UP[[i]]
  RNA_DOWN_i <- RNAseq_DOWN[[i]]
  ATAC_UP_i <- ATACseq_UP[[i]]
  ATAC_DOWN_i <- ATACseq_DOWN[[i]]
  print("index")
  print(i)
  #Track
  print(names(RNAseq_UP)[i])
  print(names(ATACseq_UP)[i])
  #Overlap in line
  Overlap_in_UP <- sort(RNA_UP_i[RNA_UP_i %in% ATAC_UP_i])
  print(length(Overlap_in_UP))
  Overlap_in_DOWN <- sort(RNA_DOWN_i[RNA_DOWN_i %in% ATAC_DOWN_i])
  print(length(Overlap_in_DOWN))
  
  #Overlap out line
  Overlap_out_UP <- RNA_UP_i[RNA_UP_i %in% ATAC_DOWN_i]
  print(length(Overlap_out_UP))
  Overlap_out_DOWN <- RNA_DOWN_i[RNA_DOWN_i %in% ATAC_UP_i]
  print(length(Overlap_out_DOWN))
  
  #RNA UP but no peak
  NonOverlap_UP_RNA <- RNA_UP_i[!RNA_UP_i %in% c(Overlap_in_UP,Overlap_out_UP)]
  print(length(NonOverlap_UP_RNA))
  #ATAC UP but no expression
  NonOverlap_UP_ATAC <- ATAC_UP_i[!ATAC_UP_i %in% c(Overlap_in_UP,Overlap_out_UP)]
  print(length(NonOverlap_UP_ATAC))
  
  #RNA DOWN but no peak
  NonOverlap_DOWN_RNA <- RNA_DOWN_i[!RNA_DOWN_i %in% c(Overlap_in_DOWN,Overlap_out_DOWN)]
  print(length(NonOverlap_DOWN_RNA))
  #ATAC DOWN but no expression
  NonOverlap_DOWN_ATAC <- ATAC_DOWN_i[!ATAC_DOWN_i %in% c(Overlap_in_DOWN,Overlap_out_DOWN)]
  print(length(NonOverlap_DOWN_ATAC))
  
  #Store counts
  comparison <- c(length(Overlap_in_UP),
                  length(Overlap_in_DOWN),
                  length(Overlap_out_DOWN),
                  length(Overlap_out_UP),
                  length(NonOverlap_UP_ATAC),
                  length(NonOverlap_UP_RNA),
                  length(NonOverlap_DOWN_ATAC),
                  length(NonOverlap_DOWN_RNA))
  #compile into table
  summary_comparison <- rbind(summary_comparison, comparison)
  
  #Store GeneSymbols
  Plus_Plus[[i]] <- Overlap_in_UP
  Minus_Minus[[i]] <- Overlap_in_DOWN
  Plus_Zero[[i]] <- NonOverlap_UP_ATAC
  Zero_Plus[[i]] <- NonOverlap_UP_RNA
  Minus_Zero[[i]] <- NonOverlap_DOWN_ATAC
  Zero_Minus[[i]] <- NonOverlap_DOWN_RNA
  
}

#Annotion of +, - and o combinations
# +/+ open/DEG
# -/- closed/suppressed
# o/- noDAR/suppressed
# +/o DAR/noDEG
colnames(summary_comparison) <- c("+/+","-/-",
                                  "+/-","-/+",
                                  "+/o","o/+",
                                  "-/o","o/-")
rownames(summary_comparison) <- names(RNAseq_UP)
print(summary_comparison)
write.table(summary_comparison, paste(path_output_Motif, "/","summaryOverlaps_ATAC_lowRNAseq.txt", sep = ""),
            row.names = TRUE,quote=FALSE,sep="\t")

#####
#Genes expressed and open
#####
DARs_features_Open_GeneSymbol <- unique(DARs_features$symbol)
length(DARs_features_Open_GeneSymbol)
expression_GeneSymbol <- unique(expression$GeneSymbol)
length(expression_GeneSymbol)
overlap_expression_open <- expression_GeneSymbol[DARs_features_Open_GeneSymbol %in% expression_GeneSymbol]
length(overlap_expression_open)

#####
#Get genes with defined peaks for Motif enrichment
#####
# Rational: 
#1) DEGs that don't overlap with a DAR are potentially influenced by TF networks
#   But Promotor regions needs to contain a peak, aka be open
#   get o/+ and o/-
#SPF
No_DAR_UP_SPF <- unique(c(Zero_Plus[[1]]))
No_DAR_DOWN_SPF <- unique(c(Zero_Minus[[1]]))
#GF
No_DAR_UP_GF <- unique(c(Zero_Plus[[2]]))
No_DAR_DOWN_GF <- unique(c(Zero_Minus[[2]]))

# Rational: 
#2) Regions that are open/closed could regulate expression upon activation/inflammation
#   driven by TF networks
#   But Promotor regions needs to contain a peak, aka be open
#   get +/o and -/o
#SPF
Open_NoDEG_SPF <- unique(c(Plus_Zero[[1]]))
Closed_NoDEG_SPF <- unique(c(Minus_Zero[[1]]))
#GF
Open_NoDEG_GF <- unique(c(Plus_Zero[[2]]))
Closed_NoDEG_GF <- unique(c(Minus_Zero[[2]]))

# Rational:
#3) DEGs that are also DARs are potentially addressed by Regulons (TF->>>> genes)
#   get +/+ and -/-
#SPF
Open_DEG_SPF <- unique(c(Plus_Plus[[1]]))
Closed_DEG_SPF <- unique(c(Minus_Minus[[1]]))
#GF
Open_DEG_GF <- unique(c(Plus_Plus[[2]]))
Closed_DEG_GF <- unique(c(Minus_Minus[[2]]))

# Rational:
#4) General DARs
Open_mLN_DAR_SPF <- subset(DARs_features_NA, log2FC_SPF_mLN_pLN >= log2FC_ATAC & padj_SPF_mLN_pLN <= 0.05)$symbol
Closed_mLN_DAR_SPF <- subset(DARs_features_NA, log2FC_SPF_mLN_pLN <= log2FC_ATAC & padj_SPF_mLN_pLN <= 0.05)$symbol
Open_mLN_DAR_GF <- subset(DARs_features_NA, log2FC_GF_mLN_pLN >= log2FC_ATAC & padj_GF_mLN_pLN <= 0.05)$symbol
Closed_mLN_DAR_GF <- subset(DARs_features_NA, log2FC_GF_mLN_pLN <= log2FC_ATAC & padj_GF_mLN_pLN <= 0.05)$symbol

#Which have a peak in the TSS region
downstream <- 200
upstream <- 10000
# or "overlapStart"

#Extract promotor region
#genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
prom <- promoters(TSS.mouse.mm10, upstream=upstream, downstream=downstream)

#Get a symbol to TSS and define promoter region
gr_TSS.mouse.mm10_ensembl <- annotatePeakInBatch(TSS.mouse.mm10, AnnotationData=TSS.mouse.mm10)
gr_TSS.mouse.mm10_symbol <- addGeneIDs(annotatedPeak=gr_TSS.mouse.mm10_ensembl, 
                                       feature_id_type="ensembl_gene_id",
                                       orgAnn="org.Mm.eg.db", 
                                       IDs2Add="symbol")
gr_prom <- promoters(gr_TSS.mouse.mm10_symbol, upstream = downstream, downstream = upstream)
#Make table
t_promotors <- as.data.frame(gr_prom)
rownames(t_promotors) <- c()


#Find TSS for identified genes
# 1) background
#all genes expressed but not differentially
expression_NA <- na.omit(expression)
only_expressed <- as.character(subset(expression_NA, 
                                      (abs(SPF_mLN_pLN_log2FoldChange) < log2FC_RNA & SPF_mLN_pLN_padj > padj) &
                                        (abs(GF_mLN_pLN_log2FoldChange) < log2FC_RNA & GF_mLN_pLN_padj > padj) &
                                        (abs(mLN_SPF_GF_log2FoldChange) < log2FC_RNA & mLN_SPF_GF_padj > padj) &
                                        (abs(pLN_SPF_GF_log2FoldChange) < log2FC_RNA & pLN_SPF_GF_padj > padj))$GeneSymbol)

#all Peaks detected in TSS region but differentially
only_open <- as.character(subset(DARs_features, 
                                 (abs(log2FC_mLN_SPF_GF) < log2FC_ATAC & padj_mLN_SPF_GF > padj) &
                                   (abs(log2FC_pLN_SPF_GF) < log2FC_ATAC & padj_pLN_SPF_GF > padj) &
                                   (abs(log2FC_SPF_mLN_pLN) < log2FC_ATAC & padj_SPF_mLN_pLN > padj) &
                                   (abs(log2FC_GF_mLN_pLN) < log2FC_ATAC & padj_GF_mLN_pLN > padj))$symbol)

background_genes <- na.omit(unique(only_expressed, only_open))
background_promotors_open <- subset(t_promotors, symbol %in% background_genes)

# 2) Key genes
#   a) list gene modules not regulated by accessibility
TF_DARs <- list(No_DAR_UP_SPF,No_DAR_DOWN_SPF,No_DAR_UP_GF,No_DAR_DOWN_GF)
#   b) list gene modules potentially regulated by activation
TF_Activity <- list(Open_NoDEG_SPF,Closed_NoDEG_SPF,Open_NoDEG_GF,Closed_NoDEG_GF)
#   c) list gene modules potentially regulating differential expression
TF_DEG_DAR <- list(Open_DEG_SPF,Closed_DEG_SPF,Open_DEG_GF,Closed_DEG_GF)
#   d) Genes generally open
TF_DARs_only <- list(Open_mLN_DAR_SPF,Closed_mLN_DAR_SPF,Open_mLN_DAR_GF,Closed_mLN_DAR_GF)

#instigate storage list for:
# a) not accessibility
t_Prom_Feat_TF_reg <- list()
for(i in seq(length(TF_DARs))){
  TF_DARs_i <- TF_DARs[[i]]
  t_Prom_Feat_i <- subset(t_promotors, symbol %in% TF_DARs_i)
  t_Prom_Feat_TF_reg[[i]] <- t_Prom_Feat_i
}
names(t_Prom_Feat_TF_reg) <- c("No_DAR_UP_SPF","No_DAR_DOWN_SPF","No_DAR_UP_GF","No_DAR_DOWN_GF")

#instigate storage list for:
# b) Activation dependent
t_Prom_Feat_TF_act <- list()
for(i in seq(length(TF_Activity))){
  TF_Activity_i <- TF_Activity[[i]]
  t_Prom_Feat_i <- subset(t_promotors, symbol %in% TF_Activity_i)
  t_Prom_Feat_TF_act[[i]] <- t_Prom_Feat_i
}
names(t_Prom_Feat_TF_act) <- c("NoDEG_Open_SPF","NoDEG_Closed_SPF","NoDEG_Open_GF","NoDEG_Closed_GF")

#instigate storage list for:
# c) Generally expressed and open genes
t_Prom_Feat_TF_act_reg <- list()
for(i in seq(length(TF_DEG_DAR))){
  TF_DEG_DAR_i <- TF_DEG_DAR[[i]]
  t_Prom_Feat_i <- subset(t_promotors, symbol %in% TF_DEG_DAR_i)
  t_Prom_Feat_TF_act_reg[[i]] <- t_Prom_Feat_i
}
names(t_Prom_Feat_TF_act_reg) <- c("Open_DEG_SPF","Closed_DEG_SPF","Open_DEG_GF","Closed_DEG_GF")

#instigate storage list for:
# d) Genes generally open
t_TF_DARs_only <- list()
for(i in seq(length(TF_DARs_only))){
  TF_DARs_only_i <- TF_DARs_only[[i]]
  t_Prom_Feat_i <- subset(t_promotors, symbol %in% TF_DARs_only_i)
  t_TF_DARs_only[[i]] <- t_Prom_Feat_i
}
names(t_TF_DARs_only) <- c("Open_mLN_DAR_SPF","Closed_mLN_DAR_SPF","Open_mLN_DAR_GF","Closed_mLN_DAR_GF")


#Regions
# a) not accessibility
for(i in seq(length(t_Prom_Feat_TF_reg))){
  t_Prom_Feat_TF_reg_i <- t_Prom_Feat_TF_reg[[i]]
  t_Prom_Feat_TF_reg_i <- t_Prom_Feat_TF_reg_i[,c("feature","symbol","seqnames","start","end")]
  names(t_Prom_Feat_TF_reg_i)[1] <- c("Acc")
  name_i <- names(t_Prom_Feat_TF_reg)[i]
  print(nrow(t_Prom_Feat_TF_reg_i))
  print(name_i)
  rownames(t_Prom_Feat_TF_act_i) <- c()
  write.table(t_Prom_Feat_TF_reg_i, paste(path_output_Motif, "/", name_i,"_lowRNAseq.txt", sep = ""),
              row.names = FALSE,quote=FALSE,sep="\t")
}

# b) Activation dependent
for(i in seq(length(t_Prom_Feat_TF_act))){
  t_Prom_Feat_TF_act_i <- t_Prom_Feat_TF_act[[i]]
  t_Prom_Feat_TF_act_i <- t_Prom_Feat_TF_act_i[,c("feature","symbol","seqnames","start","end")]
  names(t_Prom_Feat_TF_act_i)[1] <- c("Acc")
  name_i <- names(t_Prom_Feat_TF_act)[i]
  print(nrow(t_Prom_Feat_TF_act_i))
  print(name_i)
  rownames(t_Prom_Feat_TF_act_i) <- c()
  write.table(t_Prom_Feat_TF_act_i, paste(path_output_Motif, "/", name_i,"_lowRNAseq.txt", sep = ""),
              row.names = FALSE,quote=FALSE,sep="\t")
}

# c) Generally expressed and open genes
for(i in seq(length(t_Prom_Feat_TF_act_reg))){
  t_Prom_Feat_TF_act_reg_i <- t_Prom_Feat_TF_act_reg[[i]]
  t_Prom_Feat_TF_act_reg_i <- t_Prom_Feat_TF_act_reg_i[,c("feature","symbol","seqnames","start","end")]
  names(t_Prom_Feat_TF_act_reg_i)[1] <- c("Acc")
  name_i <- names(t_Prom_Feat_TF_act_reg)[i]
  print(nrow(t_Prom_Feat_TF_act_reg_i))
  print(name_i)
  rownames(t_Prom_Feat_TF_act_reg_i) <- c()
  write.table(t_Prom_Feat_TF_act_reg_i, paste(path_output_Motif, "/", name_i,"_lowRNAseq.txt", sep = ""),
              row.names = FALSE,quote=FALSE,sep="\t")
}

# d) Peaks generally open
for(i in seq(length(t_TF_DARs_only))){
  t_TF_DARs_only_i <- t_TF_DARs_only[[i]]
  t_TF_DARs_only_i <- t_TF_DARs_only_i[,c("feature","symbol","seqnames","start","end")]
  names(t_TF_DARs_only_i)[1] <- c("Acc")
  name_i <- names(t_TF_DARs_only)[i]
  print(nrow(t_TF_DARs_only_i))
  print(name_i)
  rownames(t_TF_DARs_only_i) <- c()
  write.table(t_TF_DARs_only_i, paste(path_output_Motif, "/", name_i,"_lowRNAseq.txt", sep = ""),
              row.names = FALSE,quote=FALSE,sep="\t")
}

#####
#Get peaks of defined genes for Motif enrichment
#####
# a) not accessibility
#SPF
No_DAR_UP_SPF <- unique(c(Zero_Plus[[1]]))
No_DAR_DOWN_SPF <- unique(c(Zero_Minus[[1]]))
#GF
No_DAR_UP_GF <- unique(c(Zero_Plus[[2]]))
No_DAR_DOWN_GF <- unique(c(Zero_Minus[[2]]))
#list gene modules not regulated by accessibility
TF_DARs <- list(No_DAR_UP_SPF,No_DAR_DOWN_SPF,No_DAR_UP_GF,No_DAR_DOWN_GF)
names(TF_DARs) <- c("No_DAR_UP_SPF","No_DAR_DOWN_SPF","No_DAR_UP_GF","No_DAR_DOWN_GF")
#storing list for regions
DARs_features_TSS <- list()
for(i in seq(length(TF_DARs))){
  #load genes
  DARs_genes_i <- TF_DARs[[i]]
  #identify peaks for genes
  DARs_features_i <- subset(DARs_features_NA, symbol %in% DARs_genes_i)
  #Use only genes within Promotor regions
  DARs_features_TSS_i <- subset(DARs_features_i, (distancetoFeature <= downstream & distancetoFeature >= -upstream))
  print(names(TF_DARs)[i])
  print(nrow(DARs_features_TSS_i))
  #Generate Granges object
  gr_DARs_features_TSS_i <- toGRanges(DARs_features_TSS_i)
  #adjust width
  start(gr_DARs_features_TSS_i) <- start(gr_DARs_features_TSS_i) #- 50
  end(gr_DARs_features_TSS_i) <- end(gr_DARs_features_TSS_i) #+ 50
  #Save generate .bed compatible for homer
  DARs_features_TSS_i <- as.data.frame(gr_DARs_features_TSS_i)
  DARs_features_TSS_bed_i <- DARs_features_TSS_i[,c("seqnames","start","end","id","feature_strand")] 
  #Make chr1 -> 1
  v_chr <- as.character(DARs_features_TSS_bed_i$seqnames)
  l_chr <- strsplit(v_chr, "r", fixed = FALSE, perl = FALSE, useBytes = FALSE)
  v_chr_integers <- unlist(lapply(l_chr, function(x) x[[2]]))
  DARs_features_TSS_bed_i$seqnames <- v_chr_integers
  #Generate empty column required by homer
  DARs_features_TSS_bed_i$blankVar <- NA
  DARs_features_TSS_bed_i <- DARs_features_TSS_bed_i[,c("seqnames","start","end","id","blankVar","feature_strand")] 
  DARs_features_TSS_bed_i$blankVar <- c("")
  #Store Feature BEDs 
  DARs_features_TSS[[i]] <- DARs_features_TSS_bed_i
  #Write .bed compatible for homer
  print(head(DARs_features_TSS_bed_i))
  write.table(DARs_features_TSS_bed_i, file=paste(path_output_Motif, "/", names(TF_DARs)[i],"_lowRNAseq.bed", sep = ""),
              quote=F, sep="\t", row.names=F, col.names=F)
}

# b) Activation dependent
#SPF
Open_NoDEG_SPF <- unique(c(Plus_Zero[[1]]))
Closed_NoDEG_SPF <- unique(c(Minus_Zero[[1]]))
#GF
Open_NoDEG_GF <- unique(c(Plus_Zero[[2]]))
Closed_NoDEG_GF <- unique(c(Minus_Zero[[2]]))
#   b) list gene modules potentially regulated by activation
TF_Activity <- list(Open_NoDEG_SPF,Closed_NoDEG_SPF,Open_NoDEG_GF,Closed_NoDEG_GF)
names(TF_Activity) <- c("NoDEG_Open_SPF","NoDEG_Closed_SPF","NoDEG_Open_GF","NoDEG_Closed_GF")
#storing list for regions
PotAct_features_TSS <- list()
for(i in seq(length(TF_Activity))){
  #load genes
  PotAct_genes_i <- TF_Activity[[i]]
  #identify peaks for genes
  PotAct_features_i <- subset(DARs_features_NA, symbol %in% PotAct_genes_i)
  #Use only genes within Promotor regions
  PotAct_features_TSS_i <- subset(PotAct_features_i, (distancetoFeature <= downstream & distancetoFeature >= -upstream))
  print(names(TF_Activity)[i])
  print(nrow(PotAct_features_TSS_i))
  #Generate Granges object
  gr_PotAct_features_TSS_i <- toGRanges(PotAct_features_TSS_i)
  #adjust width
  start(gr_PotAct_features_TSS_i) <- start(gr_PotAct_features_TSS_i) #- 50
  end(gr_PotAct_features_TSS_i) <- end(gr_PotAct_features_TSS_i) #+ 50
  #Save generate .bed compatible for homer
  PotAct_features_TSS_i <- as.data.frame(gr_PotAct_features_TSS_i)
  PotAct_features_TSS_bed_i <- PotAct_features_TSS_i[,c("seqnames","start","end","id","feature_strand")] 
  #Make chr1 -> 1
  v_chr <- as.character(PotAct_features_TSS_bed_i$seqnames)
  l_chr <- strsplit(v_chr, "r", fixed = FALSE, perl = FALSE, useBytes = FALSE)
  v_chr_integers <- unlist(lapply(l_chr, function(x) x[[2]]))
  PotAct_features_TSS_bed_i$seqnames <- v_chr_integers
  #Generate empty column required by homer
  PotAct_features_TSS_bed_i$blankVar <- NA
  PotAct_features_TSS_bed_i <- PotAct_features_TSS_bed_i[,c("seqnames","start","end","id","blankVar","feature_strand")] 
  PotAct_features_TSS_bed_i$blankVar <- c("")
  #Store Feature BEDs 
  PotAct_features_TSS[[i]] <- PotAct_features_TSS_bed_i
  #Write .bed compatible for homer
  #print(head(PotAct_features_TSS_bed_i))
  write.table(PotAct_features_TSS_bed_i, file=paste(path_output_Motif, "/", names(TF_Activity)[i],"_lowRNAseq.bed", sep = ""),
              quote=F, sep="\t", row.names=F, col.names=F)
}

# c) Generally expressed and open genes
#SPF
Open_DEG_SPF <- unique(c(Plus_Plus[[1]]))
Closed_DEG_SPF <- unique(c(Minus_Minus[[1]]))
#GF
Open_DEG_GF <- unique(c(Plus_Plus[[2]]))
Closed_DEG_GF <- unique(c(Minus_Minus[[2]]))
#   b) list gene modules potentially regulated by activation
TF_DEG_DAR <- list(Open_DEG_SPF,Closed_DEG_SPF,Open_DEG_GF,Closed_DEG_GF)
names(TF_DEG_DAR) <- c("Open_DEG_SPF","Closed_DEG_SPF","Open_DEG_GF","Closed_DEG_GF")

#storing list for regions
ActReg_features_TSS <- list()
for(i in seq(length(TF_DEG_DAR))){
  #load genes
  ActReg_genes_i <- TF_DEG_DAR[[i]]
  #identify peaks for genes
  ActReg_features_i <- subset(DARs_features_NA, symbol %in% ActReg_genes_i)
  #Use only genes within Promotor regions
  ActReg_features_TSS_i <- subset(ActReg_features_i, (distancetoFeature <= downstream & distancetoFeature >= -upstream))
  print(names(TF_DEG_DAR)[i])
  print(nrow(ActReg_features_TSS_i))
  #Generate Granges object
  gr_ActReg_features_TSS_i <- toGRanges(ActReg_features_TSS_i)
  #adjust width
  start(gr_ActReg_features_TSS_i) <- start(gr_ActReg_features_TSS_i) #- 50
  end(gr_ActReg_features_TSS_i) <- end(gr_ActReg_features_TSS_i) #+ 50
  #Save generate .bed compatible for homer
  ActReg_features_TSS_i <- as.data.frame(gr_ActReg_features_TSS_i)
  ActReg_features_TSS_i_bed_i <- ActReg_features_TSS_i[,c("seqnames","start","end","id","feature_strand")] 
  #Make chr1 -> 1
  v_chr <- as.character(ActReg_features_TSS_i_bed_i$seqnames)
  l_chr <- strsplit(v_chr, "r", fixed = FALSE, perl = FALSE, useBytes = FALSE)
  v_chr_integers <- unlist(lapply(l_chr, function(x) x[[2]]))
  ActReg_features_TSS_i_bed_i$seqnames <- v_chr_integers
  #Generate empty column required by homer
  ActReg_features_TSS_i_bed_i$blankVar <- NA
  ActReg_features_TSS_i_bed_i <- ActReg_features_TSS_i_bed_i[,c("seqnames","start","end","id","blankVar","feature_strand")] 
  ActReg_features_TSS_i_bed_i$blankVar <- c("")
  #Store Feature BEDs 
  ActReg_features_TSS[[i]] <- ActReg_features_TSS_i_bed_i
  #Write .bed compatible for homer
  #print(head(ActReg_features_TSS_i_bed_i))
  write.table(ActReg_features_TSS_i_bed_i, file=paste(path_output_Motif, "/", names(TF_DEG_DAR)[i],"_lowRNAseq.bed", sep = ""),
              quote=F, sep="\t", row.names=F, col.names=F)
}

# d) DARs without any 
#SPF
Open_mLN_DAR_SPF <- subset(DARs_features_NA, log2FC_SPF_mLN_pLN >= log2FC_ATAC & padj_SPF_mLN_pLN <= 0.05)$symbol
Closed_mLN_DAR_SPF <- subset(DARs_features_NA, log2FC_SPF_mLN_pLN <= log2FC_ATAC & padj_SPF_mLN_pLN <= 0.05)$symbol
Open_mLN_DAR_GF <- subset(DARs_features_NA, log2FC_GF_mLN_pLN >= log2FC_ATAC & padj_GF_mLN_pLN <= 0.05)$symbol
Closed_mLN_DAR_GF <- subset(DARs_features_NA, log2FC_GF_mLN_pLN <= log2FC_ATAC & padj_GF_mLN_pLN <= 0.05)$symbol

#   b) list gene modules potentially regulated by activation
TF_DARs_only <- list(Open_mLN_DAR_SPF,Closed_mLN_DAR_SPF,Open_mLN_DAR_GF,Closed_mLN_DAR_GF)
names(TF_DARs_only) <- c("Open_mLN_DAR_SPF","Closed_mLN_DAR_SPF","Open_mLN_DAR_GF","Closed_mLN_DAR_GF")

#storing list for regions
OpenClosed_features_TSS <- list()
for(i in seq(length(TF_DARs_only))){
  #load genes
  ActReg_genes_i <- TF_DARs_only[[i]]
  #identify peaks for genes
  ActReg_features_i <- subset(DARs_features_NA, symbol %in% ActReg_genes_i)
  #Use only genes within Promotor regions
  OpenClosed_features_TSS_i <- subset(ActReg_features_i, (distancetoFeature <= downstream & distancetoFeature >= -upstream))
  print(names(TF_DARs_only)[i])
  print(nrow(OpenClosed_features_TSS_i))
  #Generate Granges object
  gr_OpenClosed_features_TSS_i <- toGRanges(OpenClosed_features_TSS_i)
  #adjust width
  start(gr_OpenClosed_features_TSS_i) <- start(gr_OpenClosed_features_TSS_i) #- 50
  end(gr_OpenClosed_features_TSS_i) <- end(gr_OpenClosed_features_TSS_i) #+ 50
  #Save generate .bed compatible for homer
  OpenClosed_features_TSS_i <- as.data.frame(gr_OpenClosed_features_TSS_i)
  OpenClosed_features_TSS_i_bed_i <- OpenClosed_features_TSS_i[,c("seqnames","start","end","id","feature_strand")] 
  #Make chr1 -> 1
  v_chr <- as.character(OpenClosed_features_TSS_i_bed_i$seqnames)
  l_chr <- strsplit(v_chr, "r", fixed = FALSE, perl = FALSE, useBytes = FALSE)
  v_chr_integers <- unlist(lapply(l_chr, function(x) x[[2]]))
  OpenClosed_features_TSS_i_bed_i$seqnames <- v_chr_integers
  #Generate empty column required by homer
  OpenClosed_features_TSS_i_bed_i$blankVar <- NA
  OpenClosed_features_TSS_i_bed_i <- OpenClosed_features_TSS_i_bed_i[,c("seqnames","start","end","id","blankVar","feature_strand")] 
  OpenClosed_features_TSS_i_bed_i$blankVar <- c("")
  #Store Feature BEDs 
  OpenClosed_features_TSS[[i]] <- OpenClosed_features_TSS_i_bed_i
  #Write .bed compatible for homer
  #print(head(OpenClosed_features_TSS_i_bed_i))
  write.table(OpenClosed_features_TSS_i_bed_i, file=paste(path_output_Motif, "/", names(TF_DARs_only)[i],".bed", sep = ""),
              quote=F, sep="\t", row.names=F, col.names=F)
}

#Average peak region size
hist(OpenClosed_features_TSS_i_bed_i$end - OpenClosed_features_TSS_i_bed_i$start, freq = NULL)

#####
#Get Background Peaks
#####
#Three types of background are tested
# 1) All genes names with open promotors
# 2) All DARs not used for motif analsis
# 3) All TSS regions of open promotors 
#Homer Motif only needs Gene identifier in the context of TSS/promotor dependent Motif identification
# 1) Background
background_gene_symbol <- background_promotors_open[,c("feature","symbol","seqnames","start","end")]
names(background_gene_symbol)[1] <- c("Acc")
rownames(background_gene_symbol) <- NULL
write.table(background_gene_symbol, paste(path_output_Motif, "/", "background_genes_lowRNAseq.txt", sep = ""),
            row.names = FALSE,quote=FALSE,sep="\t")

# 2) ####Obtain BED for all peak files for HOMER as background
#Obtain BED for all peak files for HOMER as background
# All Peaks that are not used in any enrichment analysis
PotAct_id <- unlist(lapply(PotAct_features_TSS, function(x){x$id}))
TFs_id <- unlist(lapply(DARs_features_TSS, function(x){x$id}))
background <- DARs_regions[!(DARs_regions$id %in% c(TFs_id,PotAct_id)),]
#Generate Granges object
gr_background <- toGRanges(background)
#adjust width
start(gr_background) <- start(gr_background) #- 50
end(gr_background) <- end(gr_background) #+ 50
#Save generate .bed compatible for homer
background <- as.data.frame(gr_background)
background <- background[,c("seqnames","start","end","id","strand")]
background$blankVar <- NA
colnames(background) <- c("seqnames","start","end","id","feature_strand","blankVar")
background <- background[,c("seqnames","start","end","id","blankVar","feature_strand")] 
background$blankVar <- c("")
#annotate all strand as positive, information is not provided in the context of macs2 and ATAC-seq
background$feature_strand <- rep("+", nrow(background))
#Make chr1 -> 1
v_chr <- as.character(background$seqnames)
l_chr <- strsplit(v_chr, "r", fixed = FALSE, perl = FALSE, useBytes = FALSE)
v_chr_integers <- unlist(lapply(l_chr, function(x) x[[2]]))
background$seqnames <- v_chr_integers
#Export BED file for background
write.table(background, file=paste(path_output_Motif, "/background_all_peaks_minus_input_lowRNAseq",".bed", sep = ""),
            quote=F, sep="\t", row.names=F, col.names=F)

# 3) All TSS regions of open promotors 
#DARs_features contains all peaks
#Only the ones annotated to TSS "carry" gene name in column Symbol 
#Grab TSS that are open
DARs_features_TSS <- subset(DARs_features, (!is.na(DARs_features$symbol)))
background_genes <- unique(DARs_features_TSS[!(DARs_features_TSS$id %in% c(TFs_id,PotAct_id)),]$symbol)
#add gene name
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
listFilters(mouse)
listAttributes(mouse)
ENSEMBLE_to_GeneSymbol <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol"),
                                filters="mgi_symbol",
                                values=background_genes,
                                mart=mouse)
#Grab only background genes from TSS file
#Get Genomic locations of TSS
d_TSS.mouse.mm10 <- as.data.frame(TSS.mouse.mm10)
background_TSS <- subset(d_TSS.mouse.mm10, rownames(d_TSS.mouse.mm10) %in% ENSEMBLE_to_GeneSymbol$ensembl_gene_id)
#Make export bed file
# id is not relevant for homer in the context of background files
background_TSS <- background_TSS[,c("seqnames","start","end","width","strand")]
background_TSS$blankVar <- NA
colnames(background_TSS) <- c("seqnames","start","end","id","feature_strand","blankVar")
background_TSS <- background_TSS[,c("seqnames","start","end","id","blankVar","feature_strand")] 
background_TSS$blankVar <- c("")
#Eliminate all intel with CHR_..._PATCH
background_TSS <- background_TSS[- grep("CHR_", background_TSS$seqnames),]
#chr to integers in seqname
#background_TSS$seqnames <- paste("chr",background_TSS$seqnames,sep="")
write.table(background_TSS, file=paste(path_output_Motif, "/background_all_TSS_minus_input_genes_lowRNAseq",".bed", sep = ""),
            quote=F, sep="\t", row.names=F, col.names=F)

#############
#GO analysis Quadrants
############
#####
#Make a gene2GO list
#####
x <- org.Mm.egGO
# Get the entrez gene identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(x)
# Build a list GeneID and GOs
GeneID2GO <- list()

xx <- as.list(x[mapped_genes])

for(i in 1:length(xx)){
  #Initiate vector to collect GOIDs for geneID i
  GO_vector <- c()
  #Get geneID
  Gene_ID <- ls(xx[i])
  #grab data for geneID_i
  temp <- xx[[i]]
  
  #check Ontology category
  for(i in 1:length(temp)){
    category <- as.character(temp[[i]]["Ontology"])
    #if ontology category matches collect GOIDs  (flex)
    if(category == "BP"){
      temp_GOID <- as.character(temp[[i]]["GOID"])
      GO_vector <- c(GO_vector, temp_GOID)
    }
    #print(GO_vector)
  }
  #Generate list name geneID, content
  GO_IDs <- GO_vector 
  GeneID2GO[[Gene_ID]] <- GO_IDs 
}
GO2GeneID <- inverseList(GeneID2GO)

###
#GeneSymbol to GeneID
###
idfound <- expression$GeneSymbol %in% mappedRkeys(org.Mm.egSYMBOL)
SYMBOL <- toTable(org.Mm.egSYMBOL)
head(SYMBOL)

#Expression DEGs
m <- match(expression_NA$GeneSymbol, SYMBOL$symbol)
GENE_ID <- SYMBOL$gene_id[m]
expression <- cbind(expression_NA, GENE_ID)

#Accessibility DARs
m <- match(DARs_features_NA$symbol, SYMBOL$symbol)
GENE_ID <- SYMBOL$gene_id[m]
DARs_features_NA <- cbind(DARs_features_NA, GENE_ID)

#Define Gene Universe for which the pValue is set to 1.0
geneNames <- unique(c(levels(as.factor(DARs_features_NA$GENE_ID)),levels(as.factor(expression$GENE_ID))))
#genes_pval_background <- rep(1.0, length(geneNames))
#names(genes_pval_background) <- geneNames


###
#Perform GO analysis single
###
#put all gene vectors in one list
# Note: Input required
# Add the gene sets you are interested in
my_gene_groups <- list(Plus_Plus[[1]],Plus_Plus[[2]],
                       Minus_Minus[[1]],Minus_Minus[[2]],
                       Plus_Zero[[1]],Plus_Zero[[2]],
                       Zero_Plus[[1]],Zero_Plus[[2]],
                       Minus_Zero[[1]],Minus_Zero[[2]],
                       Zero_Minus[[1]],Zero_Minus[[2]])

names(my_gene_groups) <- c("SPF_mp_PP","GF_mp_PP",
                           "SPF_mp_MM","GF_mp_MM",
                           "SPF_mp_PZ","GF_mp_PZ",
                           "SPF_mp_ZP","GF_mp_ZP",
                           "SPF_mp_MZ","GF_mp_MZ",
                           "SPF_mp_ZM","GF_mp_ZM")

#Replace GeneSymbols with GENE_IDs and standard p-value of 0.01
#This is only legit if Fisher-Test is used for GO analysis
#only SPF
my_gene_groups <- my_gene_groups[c("SPF_mp_PP",
                                   "SPF_mp_MM",
                                   "SPF_mp_PZ",
                                   "SPF_mp_ZP",
                                   "SPF_mp_MZ",
                                   "SPF_mp_ZM")]


l_gene_pval <- list()
for(i in seq(length(my_gene_groups))){
  my_gene_group_i <- my_gene_groups[[i]]
  #replace Symbol with GENE_ID
  m <- match(my_gene_group_i, SYMBOL$symbol)
  GENE_ID <- SYMBOL$gene_id[m]
  my_gene_group_i <- GENE_ID
  
  #add standard p-Value of 0.01
  genes_pval_i <- rep(0.01, length(my_gene_group_i))
  
  #make names vector
  names(genes_pval_i) <- my_gene_group_i
  
  l_gene_pval[[i]] <- genes_pval_i
}


#get one gene list
#l_gene_pval <- l_all_gene_pval[[1]]

#gene lists
l_gene_List <- list()
for(i in 1:length(l_gene_pval)){
  geneList_i <- factor(as.integer(geneNames %in% names(l_gene_pval[[i]])))
  names(geneList_i) <- geneNames
  l_gene_List[[i]] <- geneList_i 
}


List_allRes <- list()
#Do  GO statistics for all gene lists
for(i in 1:length(l_gene_List)){
  #Access the gene lists and p-values for the differentially expressed genes
  geneList <- l_gene_List[[i]]
  pvalue_of_interest <- l_gene_pval[[i]]
  
  #build GOdata object
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel = pvalue_of_interest,
                annot = annFUN.gene2GO , gene2GO = GeneID2GO, nodeSize = 5)
  
  #get number of differentially expressed genes in the GOdata object
  sg <- sigGenes(GOdata)
  numSigGenes(GOdata)
  #get the number of GO_IDs that are within the applied GeneUniverse
  #graph(GOdata)
  number_GOIDs <- usedGO(GOdata)
  number_nodes <- length(number_GOIDs)
  
  
  #Run statistics
  
  #Fisher's works with weight
  test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
  resultWeight <- getSigGroups(GOdata, test.stat)
  #KS works with elim but not with weight
  test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
  resultElim <- getSigGroups(GOdata, test.stat)
  #runTest a high level interface for testing Fisher
  resultFis <- runTest(GOdata, statistic = "fisher")
  #Kolmogorov-Smirnov includes p-values
  test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
  resultKS <- getSigGroups(GOdata,  test.stat)
  #runTest a high level interface for testing KS
  elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  
  #make table
  allRes <- GenTable(GOdata, classic = resultFis, KS = resultKS, weight = resultWeight,
                     orderBy = "weight", ranksOf = "KS", topNodes = number_nodes)
  #make list of result tables
  List_allRes[[i]] <- allRes
}


saveRDS(List_allRes,paste(path_DAR_DEG,"/List_allRes_TopGO.Rds",sep=""))
List_allRes <- readRDS(paste(path_DAR_DEG,"/List_allRes_TopGO.Rds",sep=""))

####
#Find the Top GOs from the lists
####
#Build table with Top GOs according to "weight"
List_TopGOs <- list()
for(i in 1:length(List_allRes)){
  table_i <- List_allRes[[i]]
  table_tophits <- subset(table_i, weight < 0.005)
  #print(i)
  #print(nrow(table_tophits))
  List_TopGOs[[i]] <- table_tophits
}

#collect all TopGos
TopGOs_vector <- c()
for(i in 1:length(List_TopGOs)){
  table_i <- List_TopGOs[[i]]
  TopGOs_i <- as.character(table_i$GO.ID)
  TopGOs_vector <- c(TopGOs_vector, TopGOs_i) 
}

#Condense tables and sort by GO.ID
l_topGO <- list()
for(i in 1:length(List_allRes)){
  table_i <- List_allRes[[i]]
  table_i_subset <- subset(table_i, table_i$GO.ID %in% TopGOs_vector)
  table_i_subset_GO_weight <- table_i_subset[,c("GO.ID","weight", "Term")]
  k = i - 1
  colnames(table_i_subset_GO_weight)[2] <- names(my_gene_groups)[i]
  table_i_subset_GO_weight[,2] <- as.numeric(table_i_subset_GO_weight[,2])
  table_i_subset_GO_weight <- table_i_subset_GO_weight[order(table_i_subset_GO_weight$"GO.ID", decreasing=TRUE), ]
  l_topGO[[i]] <- table_i_subset_GO_weight
}
#merger
data_TopGO_weight = Reduce(function(...) merge(..., all=T), l_topGO)
#Delete NAs
data_TopGO_weight <- na.omit(data_TopGO_weight)

####
#Make GO comparison heatmap
####
data_heatmap <- data_TopGO_weight
row.names(data_heatmap) <- paste(data_heatmap$"GO.ID", data_heatmap$Term, sep = "_")
data_heatmap <- data_heatmap[,3:ncol(data_heatmap)]

min(data_heatmap)
data_heatmap <- -log10(data_heatmap)
#Common Minimum at 5
data_heatmap[data_heatmap > 5] <- 5
data_heatmap_matrix <- data.matrix(data_heatmap)
title <- c("GO ATAC/RNA-seq overlap")

#reorder according to Heatmap in main figure
#data_heatmap_matrix <- data_heatmap_matrix[,subset_order]

pheatmap(data_heatmap_matrix, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = FALSE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white", "grey","deepskyblue2"), space="rgb")(128),
         main = title)

#Zoom into GO heatmap
data_heatmap_zoom <- subset(data_heatmap, SPF_mp_PZ > 2 |
                                          SPF_mp_MZ > 2 |
                                          SPF_mp_ZM > 2 |
                                          SPF_mp_ZP > 2)

data_heatmap_zoom <- data_heatmap_zoom[,c("SPF_mp_PZ","SPF_mp_MZ","SPF_mp_ZM","SPF_mp_ZP")]

pheatmap(data_heatmap_zoom, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = FALSE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white", "grey","deepskyblue2"), space="rgb")(128),
         main = title)


#########
#Heatmap for set of genes in GO group
#########
Data_ready <- expression
List_GOs <- unlist(lapply(strsplit(rownames(data_heatmap_zoom), "_"), function(x) {x[[1]]}))
Names_GOs <- unlist(lapply(strsplit(rownames(data_heatmap_zoom), "_"), function(x) {x[[2]]}))
#Generate named vector to access name of GO by ID
names(Names_GOs) <- List

#Selected GOs
List_GOs_interest <- c("GO:0061314","GO:0006955","GO:0030593",
                       "GO:0034097","GO:0042127")


#Named vector GOs of interest
Names_GOs_interest[names(Names_GOs) %in% List_GOs_interest]


Genes_for_GOs_interest <- subset(GO2GeneID, ls(GO2GeneID) %in% List_GOs_interest)
number_input <- length(Genes_for_GOs_interest)
AllGenes_interest <- unlist(l_gene_pval)

#Save Output Tables in List
l_Genes_GO_per_cluster_interest <- list()

for(i in 1:number_input){
  #i = 2
  #Total number of Genes in GO
  number_genes_total_in_GO <- length(unlist(GO2GeneID[names(Genes_for_GOs_interest)[i]]))
  
  #Get genes in GO group
  GO_interest_Genes <- AllGenes_interest[names(AllGenes_interest) %in% Genes_for_GOs_interest[[i]]]
  
  #print(GO_interest_Genes)
  number_genes_in_GO <- length(GO_interest_Genes)
  print(number_genes_in_GO)
  
  #Get gene names
  gene_names <- GO_interest_Genes
  print(gene_names)
  
  #Generate list of genes of interest contained in expression data
  genes_expression <- subset(Data_ready, Data_ready$GENE_ID %in% names(gene_names))
  
  #Identify DEGs among expressed genes of GO
  genes_expression_DEG <- subset(genes_expression, abs(SPF_mLN_pLN_log2FoldChange) >= 1.0 & SPF_mLN_pLN_padj <= 0.05)
  
  number_genes_expressed_in_GO <- nrow(genes_expression_DEG)
  print(number_genes_expressed_in_GO)
  if(number_genes_expressed_in_GO >= 4){
    #format data for pheatmap
    GO_name_i <- names(Genes_for_GOs_interest)[i]
    print(GO_name_i)
    #Get name to corresponding GO ID
    
    Names_GOs_interest_i <- Names_GOs_interest[names(Names_GOs_interest) == GO_name_i]
    print(Names_GOs_interest_i)
    #Setup data for heatmap
    data_heatmap <- genes_expression_DEG[,c(1,10:21)]
    data_heatmap_matrix <- as.matrix(log2(data_heatmap[,2:ncol(data_heatmap)]))
    rownames(data_heatmap_matrix) <- data_heatmap$GeneSymbol
    
    #Find minimum generically
    find_min_data <- data_heatmap_matrix
    find_min_data[find_min_data == -Inf] <- 1000
    gen_minimum <- floor(min(find_min_data))
    data_heatmap_matrix[data_heatmap_matrix == -Inf] <- gen_minimum
    
    #Mean per sample
    mLN_GF <- rowMeans(data_heatmap_matrix[,1:3])
    pLN_GF <- rowMeans(data_heatmap_matrix[,4:6])
    mLN_SPF <- rowMeans(data_heatmap_matrix[,7:9])
    pLN_SPF <- rowMeans(data_heatmap_matrix[,10:12])
    
    #Title
    title <- paste(GO_name_i,"_",number_genes_in_GO,"/",number_genes_total_in_GO,"_",Names_GOs_interest_i,sep="")
    print(title)
    #make meaned table
    data_heatmap_matrix <- as.matrix(cbind(mLN_SPF, mLN_GF, pLN_SPF, pLN_GF))
    
    #Expression rownames
    mean_expression = round(rowMeans(data_heatmap_matrix))
    rownames_expression = paste(rownames(data_heatmap_matrix)," ", "(",mean_expression,")",sep = "")
    rownames(data_heatmap_matrix) <- rownames_expression
    
    #height of PDF
    #height_pdf <- 2 + nrow(data_heatmap_matrix) * 10 / 72
    #file_name <- paste(path_output_GO_DEG, "/",GO_name_i,"_",number_genes_expressed_in_GO,".pdf", sep="")

    #pdf(file_name, height = height_pdf, width = 5)
    
    pheatmap(data_heatmap_matrix, cluster_rows = TRUE, legend = TRUE,
             treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = FALSE,
             scale = "row", border_color = "black", cellwidth = 10,
             cellheigth = 10, color = colorRampPalette(c("blue", "white","red"), space="rgb")(128),
             main = title, fontsize = 8, fontsize_col = 12, fontsize_row = 12, fontsize_number = 12)
    #dev.off()
  }
}


####################
#Perform GO for DARs
####################
# Note: Figure for Paper
#####
#Make a gene2GO list
#####
x <- org.Mm.egGO
# Get the entrez gene identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(x)
# Build a list GeneID and GOs
GeneID2GO <- list()

xx <- as.list(x[mapped_genes])

for(i in 1:length(xx)){
  #Initiate vector to collect GOIDs for geneID i
  GO_vector <- c()
  #Get geneID
  Gene_ID <- ls(xx[i])
  #grab data for geneID_i
  temp <- xx[[i]]
  
  #check Ontology category
  for(i in 1:length(temp)){
    category <- as.character(temp[[i]]["Ontology"])
    #if ontology category matches collect GOIDs  (flex)
    if(category == "BP"){
      temp_GOID <- as.character(temp[[i]]["GOID"])
      GO_vector <- c(GO_vector, temp_GOID)
    }
    #print(GO_vector)
  }
  #Generate list name geneID, content
  GO_IDs <- GO_vector 
  GeneID2GO[[Gene_ID]] <- GO_IDs 
}
GO2GeneID <- inverseList(GeneID2GO)

#####
#GeneSymbol to GeneID
#####
DARs_features_NA <- DARs_features[!is.na(DARs_features$symbol),]
#DARs
idfound <- DARs_features_NA$symbol %in% mappedRkeys(org.Mm.egSYMBOL)
SYMBOL <- toTable(org.Mm.egSYMBOL)
head(SYMBOL)
m <- match(DARs_features_NA$symbol, SYMBOL$symbol)
GENE_ID <- SYMBOL$gene_id[m]
DARs_features_NA <- cbind(GENE_ID, DARs_features_NA)



#####
#Get gene groups
#####
DARs_Open_mLN_SPF <- subset(DARs_features_NA, log2FC_SPF_mLN_pLN >= 1.0 & padj_SPF_mLN_pLN <= 0.05)
DARs_Closed_mLN_SPF <- subset(DARs_features_NA, log2FC_SPF_mLN_pLN <= -1.0 & padj_SPF_mLN_pLN <= 0.05)



#Demethylated in mLN
genes_of_interest_1 <- unique(as.character(DARs_Open_mLN_SPF$GENE_ID))
# Note: pvalue is not used for fisher's test and filled with 0.05 for lack of p-value from Bsmooth
pvalue_of_interest_1 <- rep(0.05, length(genes_of_interest_1))
names(pvalue_of_interest_1) <- genes_of_interest_1

#Demethylated in pLN
genes_of_interest_2 <- unique(as.character(DARs_Closed_mLN_SPF$GENE_ID))
# Note: pvalue is not used for fisher's test and filled with 0.05 for lack of p-value from Bsmooth
pvalue_of_interest_2 <- rep(0.05, length(genes_of_interest_2))
names(pvalue_of_interest_2) <- genes_of_interest_2

#####
#GO analysis
#####
geneNames <- unique(c(DARs_features_NA$GENE_ID))

#gene lists
geneList_1 <- factor(as.integer(geneNames %in% genes_of_interest_1))
geneList_2 <- factor(as.integer(geneNames %in% genes_of_interest_2))

#named gene lists
names(geneList_1) <- geneNames
names(geneList_2) <- geneNames


#list the lists for looping
list_gene_List <- list(geneList_1, geneList_2)
list_pvalue_of_interest <- list(pvalue_of_interest_1, pvalue_of_interest_2)
List_allRes <- list()

#Do  GO statistics for all gene lists
for(i in 1:2){
  #Access the gene lists and p-values for the differentially expressed genes
  geneList <- list_gene_List[[i]]
  pvalue_of_interest <- list_pvalue_of_interest[[i]]
  
  #build GOdata object
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, geneSel = pvalue_of_interest,
                annot = annFUN.gene2GO , gene2GO = GeneID2GO, nodeSize = 5)
  
  #get number of differentially expressed genes in the GOdata object
  sg <- sigGenes(GOdata)
  numSigGenes(GOdata)
  #get the number of GO_IDs that are within the applied GeneUniverse
  #graph(GOdata)
  number_GOIDs <- usedGO(GOdata)
  number_nodes <- length(number_GOIDs)
  
  
  #Run statistics
  #Fishe test
  test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
  resultWeight <- getSigGroups(GOdata, test.stat)
  #KS 
  test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
  resultElim <- getSigGroups(GOdata, test.stat)
  #runTest
  resultFis <- runTest(GOdata, statistic = "fisher")
  #Kolmogorov-Smirnov -> used for further downstream analysis
  test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
  resultKS <- getSigGroups(GOdata,  test.stat)
  #runTest
  elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  
  #make table
  allRes <- GenTable(GOdata, classic = resultFis, KS = resultKS, weight = resultWeight,
                     orderBy = "weight", ranksOf = "KS", topNodes = number_nodes)
  #make list of result tables
  List_allRes[[i]] <- allRes
}


#####
#Find the Top GOs from the lists
#####
#Build table with Top GOs according to "weight"
List_TopGOs <- list()
for(i in 1:2){
  table_i <- List_allRes[[i]]
  table_tophits <- subset(table_i, weight < 0.02)
  #print(i)
  #print(nrow(table_tophits))
  List_TopGOs[[i]] <- table_tophits
}

#collect all TopGos
TopGOs_vector <- c()
for(i in 1:2){
  table_i <- List_TopGOs[[i]]
  TopGOs_i <- as.character(table_i$GO.ID)
  TopGOs_vector <- c(TopGOs_vector, TopGOs_i) 
}

table_1 <- List_allRes[[1]]
table_1_subset <- subset(table_1, table_1$GO.ID %in% TopGOs_vector)
table_1_subset_GO_weight <- table_1_subset[,c("GO.ID","weight", "Term")]
table_1_subset_GO_weight[,2] <- as.numeric(table_1_subset_GO_weight$weight)

table_2 <- List_allRes[[2]]
table_2_subset <- subset(table_2, table_2$GO.ID %in% TopGOs_vector)
table_2_subset_GO_weight <- table_2_subset[,c("GO.ID","weight", "Term")]
table_2_subset_GO_weight[,2] <- as.numeric(table_2_subset_GO_weight$weight)

#merger
a <- merge(table_1_subset_GO_weight, table_2_subset_GO_weight, by = "GO.ID")
data_TopGOs_weight <- a
#sort table
data_TopGOs_weight <- data_TopGOs_weight[,c(1,3,2,4)]
colnames(data_TopGOs_weight) <- c("GOID","GO_Term",
                                  "DAR_Open_mLN_SPF", "DMR_Closed_mLN_SPF")
#setwd("C:/Users/jpe12/PowerFolders/R/Paper/2017_Pezoldt_Pasztoi/GO_Tx_FRCs/Output")
#write.table(data_TopGOs_weight, "TopGOs_weight.txt", sep = "\t", dec = ",")

#####
#Make GO comparison heatmap
#####
GO_DARs <- data_TopGOs_weight
#Exclude biologically unconnected GO terms
GO_DARs <- subset(GO_DARs, !(GO_DARs$GOID %in% c("GO:0003151","GO:0007605",
                                                 "GO:0010470","GO:0048592",
                                                 "GO:0060004","GO:0097306",
                                                 "GO:0010460")))
row_names <- paste(GO_DARs$"GOID", GO_DARs$GO_Term, sep = "_")
row.names(GO_DARs) <- row_names
GO_DARs <- GO_DARs[,3:4]

min(GO_DARs)
GO_DARs <- -log10(GO_DARs)
#Common Minimum at 5
GO_DARs[GO_DARs > 5] <- 5
GOs_maintained_curated_matrix <- data.matrix(GO_DARs)
title <- c("GOs DARs")

pheatmap(GOs_maintained_curated_matrix, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white", "grey","deepskyblue2"), space="rgb")(128),
         main = title)

