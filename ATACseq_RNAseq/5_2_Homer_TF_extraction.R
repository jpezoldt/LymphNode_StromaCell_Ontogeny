# dissect compiled and curated Output from Homer
# Author: Joern Pezoldt
# Date: 12.10.2018
# Type of script: Exploratory

#Libraries
library("UpSetR")
library("pheatmap")

#####
#Global & PATHs
#####
table_ID <- "SPF_known_homer_compilation_curated.txt"
sample_ID <- "SPF"
PATH_input <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/", sample_ID, sep = "")
PATH_input_TFdb_Riken <- "/home/pezoldt/NAS2/pezoldt/Data/Databases/Riken_TFdB/2019-02-08_Riken_TFdb_curated.txt"
table_ID <- "SPF_known_homer_compilation_curated.txt"
# RNAseq DESeq2 analysis
path_RNAseq_DESeq2 <- "/home/pezoldt/NAS2/pezoldt/Data/RNAseq/2017_FSC_LEC_BEC_mLN_pLN_GF_SPF/FSC"
PATH_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Submission/Scatterplot"

#Global variables
log2FC_RNA = 1.0
padj = 0.05

#####
### Check for Unique TFs ---------------------------------------
#####
#Load Riken
TF_table <- read.delim(PATH_input_TFdb_Riken)
Riken_TFs <- as.character(TF_table$Symbol)

# load
t_homer_TF <- read.table(paste(PATH_input,"/",table_ID, sep=""), sep="\t",dec=".",header=TRUE)


# split according to comparisons
l_homer_TF <- split(t_homer_TF, t_homer_TF$group)
#rename
l_names <- strsplit(names(l_homer_TF), "noExt_")
names(l_homer_TF) <- unlist(sapply(l_names, function(x) {x[2]}))
#Delete global DAR sets
l_homer_TF <- l_homer_TF[c("Closed_DEG_SPF_lowRNAseq","Open_DEG_SPF_lowRNAseq","No_DAR_DOWN_SPF_lowRNAseq","No_DAR_UP_SPF_lowRNAseq",
                "NoDEG_Closed_SPF_lowRNAseq","NoDEG_Open_SPF_lowRNAseq")]
names(l_homer_TF) <- c("pLN_Open_UP","mLN_Open_UP","pLN_peak_UP","mLN_peak_UP",
                          "pLN_Open_None","mLN_Open_None")

# grab TFs per comparison
l_motifs_per_group <- lapply(l_homer_TF, function(x){
  #x <- l_homer_TF[[1]]                              
  TFs_i <- as.character(unlist(strsplit(as.character(x$Motifs_TF), ";")))
                              })
#Generate binary Matrix for upset()
t_TFs_binary <- as.data.frame.matrix(table(stack(l_motifs_per_group)))
 #Replace 
t_TFs_binary[t_TFs_binary >= 2] <- 1
t_TFs_binary <- t_TFs_binary[,c("pLN_Open_UP","mLN_Open_UP","pLN_peak_UP","mLN_peak_UP","pLN_Open_None","mLN_Open_None")]
# Plot intersection
# Note: Figure for paper
upset(t_TFs_binary, sets = colnames(t_TFs_binary),
      mainbar.y.label = "TFs ATACseq",
      number.angles = 30, point.size = 3,
      text.scale = 2, keep.order = TRUE)
outs <- upset(t_TFs_binary, sets = colnames(t_TFs_binary),
              mainbar.y.label = "TFs ATACseq",
              number.angles = 30, point.size = 3,
              text.scale = 2, keep.order = TRUE)

# Grab intersection of interest
# Note: Only one intersection interesting
Unique_TFs_NoDEG_Open <- subset(t_TFs_binary, 
                                  pLN_Open_UP == 0 &
                                  mLN_Open_UP == 0 &
                                  pLN_peak_UP == 1 &
                                  mLN_peak_UP == 0 &
                                  pLN_Open_None == 0 &
                                  mLN_Open_None == 0)

# Write binary table
write.table(t_TFs_binary, paste(PATH_input,"/",sample_ID,"_existence_matrix_known_homer_compile.txt",sep = "") ,dec=".", sep="\t")

#####Make TF binding site heatmap----------------------------------------------
#Get TF and pValue
l_homer_TF_pValue <- list()
for(i in 1:length(l_homer_TF)){
  name_i <- names(l_homer_TF)[i]
  print(name_i)
  t_i <- l_homer_TF[[i]][,c("Motifs_TF","q.value")]
  colnames(t_i) <- c("Motifs_TF", paste("pValue_",name_i,sep=""))
  l_homer_TF_pValue[[i]] <- t_i
  names(l_homer_TF_pValue)[i] <- name_i
}

# Merge tables by Motifs
t_homer_TF_pValue <- Reduce(function(...) merge(..., all=T), l_homer_TF_pValue)

# Condense per TF by meaning accross conditions
l_homer_TF_pValue <- split(t_homer_TF_pValue, t_homer_TF_pValue$Motifs_TF)
# Replace NAs with 1 and mean
l_homer_TF_pValue_mean <- lapply(l_homer_TF_pValue, function(x){
  x[is.na(x)] <- 1
  l_homer_TF_pValue_i <- colMeans(x[,2:ncol(x)])
  l_homer_TF_pValue_i
})

#rowbind list and add TF names
t_homer_TF_pValue_mean <- do.call("rbind", l_homer_TF_pValue_mean)
rownames(t_homer_TF_pValue_mean) <- names(l_homer_TF_pValue)
#eliminate first row as it is not a TF
t_homer_TF_pValue_mean <- t_homer_TF_pValue_mean[2:nrow(t_homer_TF_pValue_mean),]
#eliminate NaN which are derived from the orignial input, that contained the two lists of DARs only
t_homer_TF_pValue_mean <- t_homer_TF_pValue_mean[complete.cases(t_homer_TF_pValue_mean), ]

#-log10 of pValue for plotting
t_homer_TF_pValue_mean <- -log10(t_homer_TF_pValue_mean)
# Replace outlier high pValues with 12
t_homer_TF_pValue_mean[t_homer_TF_pValue_mean == Inf] <- 4
t_homer_TF_pValue_mean <- t_homer_TF_pValue_mean[order(rownames(t_homer_TF_pValue_mean)), ]

# Note: Figure for paper
output_pheatmap <- pheatmap(t_homer_TF_pValue_mean, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = FALSE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white","darkolivegreen1","darkolivegreen2","darkolivegreen3","darkolivegreen3","darkolivegreen4","darkgreen","darkgreen"), space="rgb")(128),
         main = "TF by motif -log10(qValue)")


####################
#RNAseq data
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
                  DF[,grepl(c("padj"), colnames(DF))],
                  DF[,grepl(c("log2FC"), colnames(DF))],
                  DF[,c("RPKMcounts.mLN_GF_FRC_1.x","RPKMcounts.mLN_GF_FRC_2.x","RPKMcounts.mLN_GF_FRC_3.x",
                        "RPKMcounts.pLN_GF_FRC_1.x","RPKMcounts.pLN_GF_FRC_2.x","RPKMcounts.pLN_GF_FRC_3.x",
                        "RPKMcounts.mLN_SPF_FRC_1.y","RPKMcounts.mLN_SPF_FRC_2.y","RPKMcounts.mLN_SPF_FRC_3.y",
                        "RPKMcounts.pLN_SFP_FRC_1.x","RPKMcounts.pLN_SFP_FRC_2.x","RPKMcounts.pLN_SFP_FRC_3.x")])
colnames(data_all)[1] <- "GeneSymbol"
#Eliminate columns
data_all <- data_all[,c(1,2,4:7,9:ncol(data_all))]
data_all <- data_all[!(is.na(data_all$GeneSymbol)),]
data_all <- data_all[!duplicated(data_all$GeneSymbol),]
#reorder
data_all <- data_all[,c(1,6:9,2:5,10:ncol(data_all))]
colnames(data_all) <- c("GeneSymbol",
                        "GF_mLN_pLN_log2FoldChange","mLN_SPF_GF_log2FoldChange",
                        "SPF_mLN_pLN_log2FoldChange","pLN_SPF_GF_log2FoldChange",
                        "GF_mLN_pLN_padj","mLN_SPF_GF_padj",
                        "SPF_mLN_pLN_padj","pLN_SPF_GF_padj",
                        "RPKM_mLN_GF_1","RPKM_mLN_GF_2","RPKM_mLN_GF_3",
                        "RPKM_pLN_GF_1","RPKM_pLN_GF_2","RPKM_pLN_GF_3",
                        "RPKM_mLN_SPF_1","RPKM_mLN_SPF_2","RPKM_mLN_SPF_3",
                        "RPKM_pLN_SPF_1","RPKM_pLN_SPF_2","RPKM_pLN_SPF_3")

#Reorder dataframe
data_all <- data_all[,c("GeneSymbol",
                        "SPF_mLN_pLN_log2FoldChange","GF_mLN_pLN_log2FoldChange",
                        "mLN_SPF_GF_log2FoldChange","pLN_SPF_GF_log2FoldChange",
                        "SPF_mLN_pLN_padj","GF_mLN_pLN_padj",
                        "mLN_SPF_GF_padj","pLN_SPF_GF_padj",
                        "RPKM_mLN_GF_1","RPKM_mLN_GF_2","RPKM_mLN_GF_3",
                        "RPKM_pLN_GF_1","RPKM_pLN_GF_2","RPKM_pLN_GF_3",
                        "RPKM_mLN_SPF_1","RPKM_mLN_SPF_2","RPKM_mLN_SPF_3",
                        "RPKM_pLN_SpF_1","RPKM_pLN_SPF_2","RPKM_pLN_SPF_3")]
#store
expression <- data_all

######
#Check expression of TFs
######
#Get rownames in order of pValue heatmap
TFs_by_Motif <- rownames(t_homer_TF_pValue_mean[output_pheatmap$tree_row[["order"]],])
TFs_by_Motif <- unlist(lapply(TFs_by_Motif, function(x){
  unlist(strsplit(x, ";", fixed = TRUE))
}))
TFs_by_Motif <- unique(TFs_by_Motif)
#Get expression values TFs identified
expression_TFs <- subset(expression, GeneSymbol %in% TFs_by_Motif)

#Prep data for heatmap
data_heatmap <- expression_TFs
rownames(data_heatmap) <- as.character(data_heatmap$GeneSymbol)
data_heatmap <- data_heatmap[,10:21]
data_heatmap <- as.matrix(data_heatmap)
data_heatmap <- log2(data_heatmap)
#Find minimum generically
find_min_data <- data_heatmap
find_min_data[find_min_data == -Inf] <- 1000
gen_minimum <- floor(min(find_min_data))
data_heatmap[data_heatmap == -Inf] <- gen_minimum
#Mean replicates
meanRPKM_mLN_GF <- rowMeans(data_heatmap[,1:3])
meanRPKM_pLN_GF <- rowMeans(data_heatmap[,4:6])
meanRPKM_mLN_SPF <- rowMeans(data_heatmap[,7:9])
meanRPKM_pLN_SPF <- rowMeans(data_heatmap[,10:12])
data_heatmap <- cbind(meanRPKM_mLN_GF,meanRPKM_pLN_GF,meanRPKM_mLN_SPF,meanRPKM_pLN_SPF)

pheatmap(data_heatmap, cluster_rows = FALSE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 30, show_rownames = TRUE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("blue", "white","red"), space="rgb")(128))

#####
#Compare TFs
#####
data_heatmap <- expression_TFs
rownames(data_heatmap) <- as.character(data_heatmap$GeneSymbol)
data_heatmap <- data_heatmap[,16:21]

#Calculate pValues
p_Values_expression_TF <- c()
for(i in 1:nrow(data_heatmap)){
  data_heatmap_row_i <- data_heatmap[i,]
  p_Value_i <- t.test(x = data_heatmap_row_i[,1:3],y =  data_heatmap_row_i[,4:6])$p.value
  p_Values_expression_TF[i] <- p_Value_i
}

#FDR adjust pValues
pAdjust_expression_TF <- p.adjust(p_Values_expression_TF, method = "fdr")
names(pAdjust_expression_TF) <- rownames(data_heatmap)

#Select significantly DEGs
sig_expression_TFs <- names(pAdjust_expression_TF[pAdjust_expression_TF < 0.10])

#Select expression data for significant pValues
data_heatmap_sig <- subset(data_heatmap, rownames(data_heatmap) %in% sig_expression_TFs)

#Set groups
group_1 = "mLN"
group_2 = "pLN"
v_classifier <- c(rep(group_1, 3),rep(group_2,3))
#List to store graphs
l_ggplot_images <- list()
for(i in 1:nrow(data_heatmap_sig)){
  TF_expression_i <- data_heatmap_sig[i,]
  print(i)
  name_TF_i <- rownames(TF_expression_i)
  print(name_TF_i)
  #vector of expression values
  Expression_i <- as.numeric(TF_expression_i)

  #Jitterplot
  t_ggplot_i <- data.frame(v_classifier = v_classifier, Expression_i = Expression_i)
  p <- ggplot(t_ggplot_i, aes(x=v_classifier, y=Expression_i, shape=v_classifier)) + geom_jitter(position=position_jitter(0.1), cex=3) + 
    scale_color_manual(values = c("blue", "red")) +
    scale_y_continuous(limits = c(floor(min(Expression_i)), ceiling(max(Expression_i)))) + 
    scale_shape_manual(values=c(1,19)) +
    xlab(name_TF_i) +
    ylab("RPKM") +
    theme(legend.position="none")
  l_ggplot_images[[i]] <- p
}
names(l_ggplot_images) <- rownames(data_heatmap_sig)

# Note Figure for paper

multiplot(plotlist = l_ggplot_images, cols = 5)
for(i in 1:length(l_ggplot_images)){
  i = 24
  setwd(paste(PATH_output,sep=""))
  gene_i <- names(l_ggplot_images)[i]
  print(gene_i)
  filename_i <- paste("Scatter_TF_", gene_i, ".eps", sep="")
  print(filename_i)
  l_ggplot_images[[i]]
  ggsave(filename_i, width = 4, height = 6, units = "cm")
  dev.off()
}

#TFs but not DEGs
TFs_not_DEG <- TFs_by_Motif[!(TFs_by_Motif %in% rownames(data_heatmap_sig))]
