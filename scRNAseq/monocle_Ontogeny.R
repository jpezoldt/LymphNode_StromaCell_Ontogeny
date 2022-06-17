#!/usr/bin/env Rscript
#args <- commandArgs(trailingOnly = T)
#blah = args[1]

# Autor: Joern Pezoldt
# 08.08.2018
# Function:
#1) Compile several datasets from kinetic into on cde object
#2) Classify cells using transcriptional signatures
#3) Check classification using scmap
#4) Perform Trajectory analysis
#5) Extract DEGs (e.g. TFs) per key cluster at branching point

#Libraries
library(monocle)
library(pheatmap)
library(reshape2)
library(cellrangerRkit)
library(biomaRt)
library(RcisTarget)
library(scran)
set.seed(123)

#install.packages('devtools')
# Replace '2.3.0' with your desired version
#library(devtools)
#devtools::install_version(package = 'Seurat', version = package_version())
#####
#Notes - Important
#####
#Script only works if all datasets have the same phenotypcial data (gene id annotation files)

#####
#Global variable
#####
condition = "day0_10_24_56_300"
datasets <- c("D000","D010","D024","D056","D300")
#Number of cells in which a gene should be expressed to be included in analysis
num_cells_exp <- 20
#Import PATH
seurat_dataset <- "_1500_1000_1_12_SC_minus_feeder1_2_3_PvC.Rds"
PATH_input_seurat <- "/home/pezoldt/NAS2/pezoldt/X_Analysis_Before2020june/Analysis/scRNAseq/seurat/Ontogney_Multialign/D0_to_D300"
#Export PATH
PATH_output <- "/home/pezoldt/NAS2/pezoldt/X_Analysis_Before2020june/Analysis/scRNAseq/monocle/day0_10_24_56_300/Submission"
PATH_CDS_objects <- "/home/pezoldt/NAS2/pezoldt/X_Analysis_Before2020june/Analysis/scRNAseq/monocle/day0_10_24_56_300/Submission/CDS_objects"

#Signatures from mLN
NC_Pezoldt_Pasztoi_2018 <- read.csv("/home/pezoldt/NAS2/pezoldt/X_Analysis_Before2020june/Analysis/scRNAseq/2018_NatComm_Pezoldt_Pasztoi/NC_2018_mLN_Clusters.csv",sep = ";")
#Riken TF database
PATH_input_TFdb_Riken <- "/home/pezoldt/NAS2/pezoldt/Data/Databases/Riken_TFdB/2019-02-08_Riken_TFdb_curated.txt"
#TF from DARs
PATH_input_TFBS_DARs <- "/home/pezoldt/NAS2/pezoldt/X_Analysis_Before2020june/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/SPF/SPF_existence_matrix_known_homer_compile.txt"
#Secreted proteins Uniprot
PATH_input_secreted_proteins <- "/home/pezoldt/NAS2/pezoldt/X_Analysis_Before2020june/Data/Databases/Uniprot/2019-03-23_Uniprot_secreted.tab"
#TF from DMRs
PATH_input_DMRs <- "/home/pezoldt/NAS2/pezoldt/X_Analysis_Before2020june/Analysis/WGBS/FSC/Submission/DMR_TFs_TSS_vicinity.csv"

#####
#Load RDS objects
#####
cde_SC <- readRDS(file=paste(PATH_CDS_objects,"/day0_10_24_56_300_Seurat_1500_1000_1_12_SC_minus_feeder1_2_3_PvC.Rds",sep=""))
cde_NonAdv <- readRDS(paste(PATH_CDS_objects,"/day0_10_24_56_300_Seurat_1500_1000_1_12_SC_minus_feeder1_2_3_PvC_only_NonAdventi_.Rds",sep=""))
cde_Adv <- readRDS(paste(PATH_CDS_objects,"/day0_10_24_56_300_Seurat_1500_1000_1_12_SC_minus_feeder1_2_3_PvC_only_Adventi_.Rds",sep=""))


#####
#Load and compile data into cde object
#####
#Timepoint 1: e.g. neonatal/d0
dir_d0 <- "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/C11_neo_mLN"
d0 <- load_cellranger_matrix(dir_d0, genome = "mm10")

#Timepoint 2: e.g day10
dir_d10 <- "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/C12_d10_mLN"
d10 <- load_cellranger_matrix(dir_d10, genome = "mm10")

#Timepoint 3: e.g. day24
dir_d24 <- "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/270_2018-01-03_scRNA-Seq_mLN_ontogeny/Data/D1_d24_mLN"
d24 <- load_cellranger_matrix(dir_d24, genome = "mm10")

#Timepoint 4: e.g. day56 (merge two datasets)
#Load as seurat
# Cell IDs: 1 to nrow(Experiment1)
dir_d56_1 <- "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/244_scRNA-Seq_mLN_pLN_SPF/Data/10X_results/L1700567_mLN_SPF"
d56_1 <- load_cellranger_matrix(dir_d56_1, genome = "mm10")
#sample_4_1.data@Dimnames[[2]] <- paste("D056_", c(1:ncol(sample_4_1.data)), sep = "")
#sample_4_1_seurat <- CreateSeuratObject(raw.data = sample_4_1.data, min.cells = num_cells_exp, min.genes = 1000)
#Experiment 2
# Cell IDs: nrow(Experiment1) to nrow(Experiment2)
#sample_4.data <- Read10X(data.dir = "C:/Users/jpe12/PowerFolders/R/2017_252_scRNAseq/Data/mLN_SPF/B6/filtered_gene_bc_matrices/mm10")
dir_d56_2 <- "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/252_2017-07-19_scRNA-Seq_mLN_pLN_SPF_GF/Data/10x/mLN_SPF/B6"
d56_2 <- load_cellranger_matrix(dir_d56_2, genome = "mm10")
#Merge datasets
#pd <- rbind(pData(d56_1),pData(d56_2))
fd <- fData(d56_1)
names(fd) <- c('id', 'gene_short_name')
expr_merged <- as(cbind(exprs(d56_1), exprs(d56_2)), "dgTMatrix")
colnames(expr_merged) <- paste("d056_", 1:ncol(expr_merged), sep = "")
pd <- data.frame(barcode = colnames(expr_merged))
rownames(pd) <- colnames(expr_merged)
d56 <- newCellDataSet(expr_merged, 
                      phenoData = new("AnnotatedDataFrame", data = pd),
                      featureData = new("AnnotatedDataFrame", data = fd))

#Timepoint 6: e.g. day300
dir_d300 <- "/home/pezoldt/NAS2/pezoldt/Data/scRNAseq/275_2018-04-23_scRNA-Seq_mLN_SPF_45wk/Data/E7_SPF_mLN_42wk"
d300 <- load_cellranger_matrix(dir_d300, genome = "mm10")

#list datasets
l_GBCMs <- list(d0,d10,d24,d56,d300)
#,d24,d56)
names(l_GBCMs) <- datasets

#quick glance at data (first entry number of genes, second entry number of cells)
lapply(l_GBCMs, function(x) dim(exprs(x)))

#Rename first columns
l_cds <- lapply(l_GBCMs, function(x) {
  #x <- l_GBCMs[[1]]
  feat <- fData(x)
  names(feat) <- c('id', 'gene_short_name')
  newCellDataSet(exprs(x),
                 phenoData = new("AnnotatedDataFrame", data = pData(x)),
                 featureData = new("AnnotatedDataFrame", data = feat),
                 lowerDetectionLimit = 1,
                 expressionFamily = negbinomial.size())
})

#####
#Keep only specific cells for trajectory analysis dataset
#####
# Note: Selection performed via Seurat
# At this stage of data processing both d0 Seurat and monocle object have the same number of cells
# Choose cells over index of row
for(i in 1:length(datasets)){
  # i = 4
  dataset_time_point_i <- datasets[i]
  print(dataset_time_point_i)
  cells_include_i <- readRDS(paste(PATH_input_seurat,"/",dataset_time_point_i,seurat_dataset, sep=""))
  print(length(cells_include_i))
  #Row index
  index <- as.numeric(unlist(lapply(strsplit(cells_include_i, "_"), function(x){x[[2]]})))
  print(length(index))
  #make cds object for selection
  print(ncol(exprs(l_cds[[i]])))
  expression <- exprs(l_cds[[i]])[,index]
  phenoData <- data.frame(barcode = pData(l_cds[[i]])[index,])
  rownames(phenoData) <- phenoData$barcode.barcode
  feat <- fData(l_cds[[i]])
  names(feat) <- c('id', 'gene_short_name')
  cds_cleared_i <- newCellDataSet(expression,
                                  phenoData = new("AnnotatedDataFrame", data = phenoData),
                                  featureData = new("AnnotatedDataFrame", data = feat),
                                  lowerDetectionLimit = 1,
                                  expressionFamily = negbinomial.size())
  
  #replace element in list of cds
  l_cds[[i]] <- cds_cleared_i
}


#Add time_point and Cell_id column to pData
for(i in seq(length(datasets))){
  dataset_time_point_i <- datasets[i]
  print(datasets[i])
  pData(l_cds[[i]])$Time_point <- rep(datasets[i], nrow(pData(l_cds[[i]])))
  Cell_id <- readRDS(paste(PATH_input_seurat,"/",dataset_time_point_i,seurat_dataset, sep=""))
  pData(l_cds[[i]])$Cell_id <- Cell_id
}

#####Merge datasets to obtain one cde object
#build lists for pData
l_pData <- list()
#l_fData <- list()
l_expData <- list()
for(i in seq(length(datasets))){
  l_pData[[i]] <- pData(l_cds[[i]])
  l_expData[[i]] <- l_cds[[i]]@assayData$exprs
}
lapply(l_pData, nrow)
lapply(l_expData, nrow)
lapply(l_expData, ncol)

#rbind pData
pData_all_t <- do.call("rbind", l_pData)
rownames(pData_all_t) <- pData_all_t$Cell_id
#sparse Matrices
#Note: Only possible if the colnames (gene identifiers are identical across all conditions)
expData_all_t <- do.call("cbind", l_expData)
colnames(expData_all_t) <- pData_all_t$Cell_id
nrow(pData_all_t)
ncol(expData_all_t)
nrow(expData_all_t)
nrow(fData(l_cds[[1]]))

#Full data cde object
cde_all <- newCellDataSet(expData_all_t,
                          phenoData = new("AnnotatedDataFrame", data = pData_all_t),
                          featureData = new("AnnotatedDataFrame", data = fData(l_cds[[1]])),
                          lowerDetectionLimit = 1,
                          #for UMI
                          expressionFamily = negbinomial.size())
dim(cde_all@assayData$exprs)
dim(pData(cde_all))
dim(fData(cde_all))

#####
#Attach MetaData of Seurat to pData moncole
#####
Meta_Data_Seurat <- readRDS(paste(PATH_input_seurat,"/MetaData",seurat_dataset,sep=""))
head(Meta_Data_Seurat)
print(paste("Number of Cells:",nrow(Meta_Data_Seurat)))
head(pData(cde_all))
print(paste("Number of Cells:",nrow(pData(cde_all))))
#Merge by rownames
pData_new <- merge(pData(cde_all),Meta_Data_Seurat, by=0, all=TRUE)
#Drop columns only containing seurat data
pData_new <- subset(pData_new, !(barcode.barcode == "<NA>"))
#Cluster_names_seurat <- pData_new$Cluster
rownames(pData_new) <- pData_new$Row.names
colnames(pData_new)[ncol(pData_new)] <- "Cluster_Seurat"
#Get order of cells in expression table
names_exprs <- colnames(exprs(cde_all))
#Sort tables according to cell IDs and time points
pData_new <- pData_new[match(names_exprs, pData_new$Row.names),]

#pData_new$Row.names <- as.numeric(unlist(lapply(strsplit(pData_new$Row.names, "_"), function(x){x[2]})))
#l_pData_new <- split(pData_new, pData_new$Time_point)
#for(i in 1:length(l_pData_new)){
# pData_new_i <- l_pData_new[[1]]
#pData_new_i <- pData_new_i[order(pData_new_i$Row.names),]
#}
#Replace pData
pData(cde_all) <- pData_new

#####
#QC
#####
#Estimate factor matrices
cde_all <- estimateSizeFactors(cde_all)
cde_all <- estimateDispersions(cde_all)

#Filter low qualitiy cells
cde_all <- detectGenes(cde_all, min_expr = 0.1)
print(head(fData(cde_all)))
nrow(fData(cde_all))

#Thresh by number of cells gene is expressed in
expressed_genes <- row.names(subset(fData(cde_all), num_cells_expressed >= num_cells_exp))
saveRDS(expressed_genes, paste(PATH_output,"/expressed_genes.Rds",sep = ""))
print(head(pData(cde_all)))

#Check distribution of expression across condition
pData(cde_all)$Total_mRNAs <- Matrix::colSums(exprs(cde_all))
pData(cde_all)$Total_n_genes <- nrow(fData(cde_all)) - Matrix::colSums(exprs(cde_all)==0)
print(head(pData(cde_all)))

#Estimate cutoff
hist(log2(Matrix::colSums(exprs(cde_all))))

#Plot detected gene number range over conditions
upper_bound <- 250
lower_bound <- 5000
qplot(Total_n_genes, data = pData(cde_all), color = Time_point, geom = "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)
levels(as.factor(pData(cde_all)$Time_point))

#Ribosomal protein genes
rpl.genes <- grep(pattern = "^Rpl", x = fData(cde_all)$gene_short_name, value = TRUE)
rps.genes <- grep(pattern = "^Rps", x = fData(cde_all)$gene_short_name, value = TRUE)
ribo.genes <- c(rpl.genes, rps.genes)
#get ensemble Ids
ribo.genes.ensembl <- subset(fData(cde_all), gene_short_name %in% ribo.genes)$id
#subset expression matrix and perform colmeans per cell
ribo.genes_expr <- colSums(as.matrix(exprs(cde_all)[rownames(exprs(cde_all))%in% ribo.genes.ensembl,]))/colSums(as.matrix(exprs(cde_all)))
pData(cde_all)$ribo_expression <- ribo.genes_expr
#Add Gene signature averages to pData
#Mitochondrial genes
mito.genes <- grep(pattern = "^mt-", x = fData(cde_all)$gene_short_name, value = TRUE)
#get ensemble Ids
mito.genes.ensembl <- subset(fData(cde_all), gene_short_name %in% mito.genes)$id
#subset expression matrix and perform colmeans per cell
mito.genes_expr <- colSums(as.matrix(exprs(cde_all)[rownames(exprs(cde_all))%in% mito.genes.ensembl,]))/colSums(as.matrix(exprs(cde_all)))
pData(cde_all)$mito_expression <- mito.genes_expr

#Get an overview of distributions
#Take only cells were number of mRNA and genes correlate linear on main axis

# to do so: calculate regression per data point -> eliminate cells with too steep regression -> see below
#Regression of n_gene and total
pData(cde_all)$exp_regression <- pData(cde_all)$Total_mRNAs / pData(cde_all)$Total_n_genes
dim(cde_all@assayData$exprs)
dim(pData(cde_all))
dim(fData(cde_all))

#####
#Define subsets
#####
#2) Via variable genes
#Identify variable genes
disp_table <- dispersionTable(cde_all)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cde_all <- setOrderingFilter(cde_all, unsup_clustering_genes$gene_id)
#black dots are used for ordering
nrow(unsup_clustering_genes)

#substract/regress factors
#Notes
#ribo and number of expressed genes make no big impact on clustering
#Based on the Seurat Output containing only the non-endothelial stromal cells cde_all is already preselected
cde_SC <- cde_all
cde_SC <- reduceDimension(cde_SC, max_components = 2, num_dim = 20,
                          reduction_method = 'tSNE',
                          residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression  + Time_point",
                          verbose = T)
cde_SC <- clusterCells(cde_SC, num_clusters = 15)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster")
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Il15")) #+ facet_wrap(~Time_point)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Pecam1","Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                "Cd248","Ackr3","Cxcl13","Nes","Cdk1"),cell_size = 0.5)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster") + facet_wrap(~Time_point)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster") + facet_wrap(~Time_point)
#plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Cd34")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                "Cd248","Ackr3","Cdk1","Nes"),cell_size = 0.5)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Tagln"),
                   cell_size = 0.5)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Cxcl13"),
                   cell_size = 0.5)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Tnfsf11"),
                   cell_size = 0.5)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Il21"),
                   cell_size = 0.5)

plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Acta2","Tagln","Lmod1","Tinagl1"),
                   cell_size = 0.5)
plot_cell_clusters(cde_SC, 1, 2, color = "Cluster", markers = c("Vcam1","Icam1","Tnfsf11","Cxcl13"),
                   cell_size = 0.5)

###############
#Analyze SC cell
###############
#########
#Building trajectories SC------------------------------
#########
cde_SC <- reduceDimension(cde_SC, max_components = 2, num_dim = 20,
                          reduction_method = 'tSNE',
                          residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression + Time_point",
                          verbose = T)
cde_SC <- clusterCells(cde_SC, num_clusters = 15)

#Reduce dimensionality
cde_SC <- reduceDimension(cde_SC, max_components = 2, method = 'DDRTree')
#Build trajectory
cde_SC <- orderCells(cde_SC)

#Inferring the genes----------------------------------------------
plot_cell_trajectory(cde_SC, color_by = "State",
                     cell_size = 0.25, show_branch_points = FALSE)
ggsave("FSC_All_State.eps", width = 5.0, height = 7.4, units = "cm")
ggsave(paste(PATH_output,"/Figures/Trajectories_SC/FSC_All_State.png",sep=""), width = 5.0, height = 7.7, units = "cm")
plot_cell_trajectory(cde_SC, markers = c("Cdk1","Aldh1a2","Cxcl13","Cd34","Acta2"), use_color_gradient = TRUE, cell_size = 0.5)

#Note: Figure for paper
plot_cell_trajectory(cde_SC, markers = c("Bst1"), use_color_gradient = TRUE, cell_size = 0.5)
plot_cell_trajectory(cde_SC, markers = c("Vcam1"), use_color_gradient = TRUE, cell_size = 0.5)
plot_cell_trajectory(cde_SC, markers = c("Cd34"), use_color_gradient = TRUE, cell_size = 0.25)


#Ludewig 2019, NC UP
plot_cell_trajectory(cde_SC, markers = c("Isl1","Nkx2-5","Pdgfra","Ptdfrb",
                                         "Ltbr","Acta2","Ccl19","Cd248",
                                         "Itga3","Gli1","Nkx3-2","Nkx2-5",
                                         "Wt1","Col14a1","Tnxb","Ly6a",
                                         "Eng","Vcam1","Itgb1",
                                         "Cdk1","Cxcl13","Cd34"), use_color_gradient = TRUE, cell_size = 0.25, show_branch_points = FALSE)



gene_list_map_on_trajectory <- c("Acta2",
                                 "Bst1",
                                 "Cd34",
                                 "Cdk1",
                                 "Icam1",
                                 "Vcam1",
                                 "Cxcl13",
                                 #Ludewig 2019, NC
                                 "Isl1","Nkx2-5","Pdgfra","Ptdfrb",
                                 "Ltbr","Acta2","Ccl19","Cd248",
                                 "Itga3","Gli1","Nkx3-2","Nkx2-5",
                                 "Wt1","Col14a1","Tnxb","Ly6a",
                                 "Eng","Vcam1","Itgb1")
for(i in 1:length(gene_list_map_on_trajectory)){
  # i = 1
  setwd(paste(PATH_output,"/Figures/Trajectories_SC",sep=""))
  gene_i <- gene_list_map_on_trajectory[i]
  print(gene_i)
  filename_i <- paste("FSC_", gene_i, ".png", sep="")
  print(filename_i)
  plot_cell_trajectory(cde_SC, markers = c(gene_i), use_color_gradient = TRUE,
                       cell_size = 0.25, show_branch_points = FALSE)
  
  ggsave(filename_i, width = 5.0, height = 7.9, units = "cm")
}
for(i in 1:length(gene_list_map_on_trajectory)){
  # i = 1
  setwd(paste(PATH_output,"/Figures/Trajectories_SC",sep=""))
  gene_i <- gene_list_map_on_trajectory[i]
  print(gene_i)
  filename_i <- paste("FSC_", gene_i, ".eps", sep="")
  print(filename_i)
  plot_cell_trajectory(cde_SC, markers = c(gene_i), use_color_gradient = TRUE,
                       cell_size = 0.25, show_branch_points = FALSE)
  
  ggsave(filename_i, width = 5.0, height = 7.9, units = "cm")
}
for(i in 1:length(gene_list_map_on_trajectory)){
  # i = 1
  setwd(paste(PATH_output,"/Figures/Trajectories_SC",sep=""))
  gene_i <- gene_list_map_on_trajectory[i]
  print(gene_i)
  filename_i <- paste("FSC_", gene_i, "_Scale", ".eps", sep="")
  print(filename_i)
  plot_cell_trajectory(cde_SC, markers = c(gene_i), use_color_gradient = TRUE,
                       cell_size = 0.25, show_branch_points = FALSE)
  
  ggsave(filename_i, width = 7.0, height = 7.9, units = "cm")
}







#Set root state
cde_SC <- orderCells(cde_SC, root_state = 5)
nrow(subset(pData(cde_SC), State == 3))
plot_cell_trajectory(cde_SC, color_by = "State")
plot_cell_trajectory(cde_SC,cell_size = 0.2, color_by = "State", show_branch_points = FALSE) + 
  facet_wrap(~Time_point)
#Note: Figure for paper
plot_cell_trajectory(cde_SC, color_by = "State", cell_size = 0.2) + 
  facet_wrap(~Cluster_Seurat) +
  stat_density2d(color='black', h = 6, alpha=I(0.5), size=I(0.5))
plot_cell_trajectory(cde_SC, color_by = "State", markers = c("Vcam1"), use_color_gradient = TRUE)  + facet_wrap(~Time_point)
plot_cell_trajectory(cde_SC, markers = c("Cd34"), use_color_gradient = TRUE) + facet_wrap(~Time_point)

plot_complex_cell_trajectory(cde_SC, color_by = 'State', show_branch_points = T,
                             cell_size = 0.5, cell_link_size = 0.3, root_states = c(1))  + facet_wrap(~Time_point)
plot_complex_cell_trajectory(cde_SC, color_by = 'State', show_branch_points = T,
                             cell_size = 0.5, cell_link_size = 0.3, root_states = c(2))  + facet_wrap(~Cluster_Seurat)

#Save RDS
saveRDS(cde_SC, file=paste(PATH_CDS_objects,"/day0_10_24_56_300_Seurat_1500_1000_1_12_SC_minus_feeder1_2_3_PvC.Rds",sep=""))
cde_SC <- readRDS(file=paste(PATH_CDS_objects,"/day0_10_24_56_300_Seurat_1500_1000_1_12_SC_minus_feeder1_2_3_PvC.Rds",sep=""))


#####
#Extract RCM of CDE per Stage
#####
# Note: Input required
cde_RCM_export <- cde_nonPvC
name <- "RCM_cde_nonPvC"

#Prep Biomart intel
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=genes,mart= mart)

#Stages
States <- levels(pData(cde_RCM_export)$State)
for(i in 1:length(States)){
  #Get State
  State_i <- as.numeric(States[i])
  print(State_i)
  #Get cells in state and obtain RCM
  Cells_state_i <- row.names(subset(pData(cde_RCM_export), Time_point == "d024" & State == State_i))
  RCM_export_i <- as.matrix(exprs(cde_RCM_export)[,Cells_state_i])
  
  #Rename RCM
  genes <- rownames(RCM_export_i)
  RCM_export_i <- cbind(genes, RCM_export_i)
  print(ncol(RCM_export_i))
  rownames(RCM_export_i) <- c()
  RCM_export_i <- merge(G_list,RCM_export_i,by.x="ensembl_gene_id",by.y="genes")
  RCM_export_i <- RCM_export_i[,2:ncol(RCM_export_i)]
  #Export RCM for usage in progenitor profile
  write.table(RCM_export_i, paste(PATH_output,"/",name,"_d24_",State_i,".txt",sep=""), sep = "\t")
}

save_pData <- subset(pData(cde_RCM_export), Time_point == "d024")
write.table(save_pData, paste(PATH_output,"/",name,"_d24_phenotypicalData",".txt",sep=""), sep = "\t")

#Generate ReadCountMatrix with GeneSymbols for export
RCM_export <- as.matrix(exprs(cde_RCM_export))
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- rownames(RCM_export)
RCM_export <- cbind(genes, RCM_export)
rownames(RCM_export) <- c()
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=genes,mart= mart)
RCM_export <- merge(G_list,RCM_export,by.x="ensembl_gene_id",by.y="genes")
RCM_export <- RCM_export[,2:ncol(RCM_export)]
write.table(RCM_export, paste(PATH_output,"/",name,"txt",sep=""), sep = "\t")

#####
#Extract RCM of CDE per Cluster
#####
# Note: Input required
cde_RCM_export <- cde_selected
name <- "RCM_Progenitors"

#Prep Biomart intel
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=genes,mart= mart)

#Stages
Clusters <- levels(as.factor(pData(cde_RCM_export)$Cluster))
for(i in 1:length(Clusters)){
  #Get State
  Clusters_i <- as.numeric(Clusters[i])
  #Get cells in state and obtain RCM
  Cells_cluster_i <- row.names(subset(pData(cde_RCM_export), Cluster == Clusters_i))
  RCM_export_i <- as.matrix(exprs(cde_RCM_export)[,Cells_cluster_i])
  print(Clusters_i)
  print(ncol(RCM_export_i))
  #Rename RCM
  genes <- rownames(RCM_export_i)
  RCM_export_i <- cbind(genes, RCM_export_i)
  rownames(RCM_export_i) <- c()
  RCM_export_i <- merge(G_list,RCM_export_i,by.x="ensembl_gene_id",by.y="genes")
  RCM_export_i <- RCM_export_i[,2:ncol(RCM_export_i)]
  #Export RCM for usage in progenitor profile
  write.table(RCM_export_i, paste(PATH_output,"/",name,"_Cluster_","_d56_d300_State7",".txt",sep=""), sep = "\t")
}

#combine defined clusters
Cluster_3 <- read.delim(paste(PATH_output,"/",name,"_Cluster_","3",".txt",sep=""))
Cluster_4 <- read.delim(paste(PATH_output,"/",name,"_Cluster_","4",".txt",sep=""))
Cluster_3_4 <- cbind(Cluster_3,Cluster_4)
write.table(Cluster_3_4, paste(PATH_output,"/",name,"_Cluster_","3_4",".txt",sep=""), sep = "\t")


#Generate ReadCountMatrix with GeneSymbols for export
RCM_export <- as.matrix(exprs(cde_RCM_export))
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- rownames(RCM_export)
RCM_export <- cbind(genes, RCM_export)
rownames(RCM_export) <- c()
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=genes,mart= mart)
RCM_export <- merge(G_list,RCM_export,by.x="ensembl_gene_id",by.y="genes")
RCM_export <- RCM_export[,2:ncol(RCM_export)]
write.table(RCM_export, paste(PATH_output,"/",name,".txt",sep=""), sep = "\t")

################
#select Stage / Cluster and rerun Trajectory
################
#Select stage 2 (earliest developmental stage)
cells_NonAdv <- row.names(subset(pData(cde_SC),
                                 State %in% c("4","5","6")))
cells_Adv <- row.names(subset(pData(cde_SC),
                              State %in% c("1","2","3","7")))

##### NonAdv
cde_NonAdv<- cde_SC[,cells_NonAdv]
# Clustering
cde_NonAdv <- reduceDimension(cde_NonAdv, max_components = 2, num_dim = 12,
                              reduction_method = 'tSNE',
                              residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression + Time_point",
                              verbose = T)
cde_NonAdv <- clusterCells(cde_NonAdv, num_clusters = 6)
plot_cell_clusters(cde_NonAdv, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10)
plot_cell_clusters(cde_NonAdv, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10) + facet_wrap(~Time_point)
plot_cell_clusters(cde_NonAdv, 1, 2, color = "State", show_cell_names = TRUE, cell_name_size = 10) + facet_wrap(~Time_point)
plot_cell_clusters(cde_NonAdv, 1, 2, color = "Cluster", markers = c("Cxcl13")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_NonAdv, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                    "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                    "Cd248","Gdf10","Cxcl9","Has1","Ackr3",
                                                                    "Cdk1","Il6"), cell_size = 0.5)

#Ludewig 2019, NC UP
plot_cell_trajectory(cde_NonAdv, markers = c("Isl1","Nkx2-5","Pdgfra","Ptdfrb",
                                             "Ltbr","Acta2","Ccl19","Cd248",
                                             "Itga3","Gli1","Nkx3-2","Nkx2-5",
                                             "Wt1","Col14a1","Tnxb","Ly6a",
                                             "Eng","Vcam1","Itgb1"), use_color_gradient = TRUE, cell_size = 0.25, show_branch_points = FALSE)

#Calculate % per cluster across age
meta_data <- pData(cde_Adv)
l_condition <- split(meta_data, meta_data$tech)
l_cells_per_condition <- lapply(l_condition, function(x) {
  #x <- l_condition[[1]]
  x_cells <- unlist(lapply(split(x, x$Cluster_Seurat), nrow))
  x_cells
})
names(l_cells_per_condition) <- names(l_condition)
t_cells_per_condition <- do.call(rbind.fill.matrix, lapply(l_cells_per_condition, function(x) {t(as.data.frame(x))}))
t_cells_per_condition[is.na(t_cells_per_condition)] <- 0
rownames(t_cells_per_condition) <- names(l_condition)
#Thresh for subsets with at least defined cell number
min_num_cell_cluster = 30
cells_per_cluster <- colSums(t_cells_per_condition)
clusters_present <- cells_per_cluster[cells_per_cluster > min_num_cell_cluster]
t_cells_per_condition <- t_cells_per_condition[,names(clusters_present)]
#Normalize to 1000 cells
norm_factor_vector <- rowSums(t_cells_per_condition) / 1000
#Divide cell number per cluster by normalization factor
t_cells_per_condition_norm <- t_cells_per_condition / norm_factor_vector
#Calculate Frequencies
total <- colSums(t_cells_per_condition_norm)
freq_cells_per_condition <- (t(t_cells_per_condition_norm) / total) * 100
freq_cells_per_condition <- as.matrix(freq_cells_per_condition)
pheatmap(freq_cells_per_condition, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = FALSE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white","grey","grey23","black"), space="rgb")(128))


#Trajectories
cde_NonAdv <- reduceDimension(cde_NonAdv, max_components = 2, method = 'DDRTree')
#Build trajectory
cde_NonAdv <- orderCells(cde_NonAdv)
#Note: Figure for Paper
plot_cell_trajectory(cde_NonAdv, color_by = "Pseudotime", cell_size = 0.25, show_tree = TRUE,
                     show_branch_points = FALSE)
#cde_NonAdv root state
#cde_SC <- orderCells(cde_NonAdv, root_state = 5)
plot_cell_trajectory(cde_NonAdv, color_by = "State")
plot_cell_trajectory(cde_NonAdv, color_by = "State") + facet_wrap(~Time_point)
#Note: Figure for Paper
plot_cell_trajectory(cde_NonAdv, color_by = "State", cell_size = 0.25) + 
  facet_wrap(~Cluster_Seurat) +
  stat_density2d(color='black', h = 6, alpha=I(0.5), size=I(0.5))

#Note: Figure for Paper
plot_cell_trajectory(cde_NonAdv, color_by = "State")  + facet_wrap(~Time_point)
plot_cell_trajectory(cde_NonAdv,
                     color_by = "State",
                     cell_size = 0.25,
                     show_branch_points = FALSE) + 
  facet_wrap(~Time_point) +
  scale_color_manual(breaks = c("1", "2", "3","4","5"), 
                     values=c("#AB886C", "#A33C33", "#C4847B","#600B00","#522EDB")) +
  stat_density2d(color='black', h = 6, alpha=I(0.5), size=I(1))
plot_cell_trajectory(cde_NonAdv,
                     color_by = "State",
                     cell_size = 0.25,
                     show_branch_points = FALSE) + 
  facet_wrap(~Cluster_Seurat) +
  scale_color_manual(breaks = c("1", "2", "3","4","5"), 
                     values=c("#522EDB", "#A33C33", "#C4847B","#600B00","#AB886C")) +
  stat_density2d(color='black', h = 6, alpha=I(0.5), size=I(0.5))
plot_cell_trajectory(cde_NonAdv, markers = c("Inmt","Cxcl13","Cxcl1","Ccl19","Il6","Cxcl12","Ccl19","Inmt","Il7","Tnfsf11"), use_color_gradient = TRUE, cell_size = 0.25, show_branch_points = FALSE)
plot_cell_trajectory(cde_NonAdv, markers = c("Stat1"),
                     use_color_gradient = TRUE,
                     cell_size = 0.25,
                     show_branch_points = FALSE) + 
  facet_wrap(~Cluster_Seurat) +
  #scale_color_manual(breaks = c("1", "2", "3","4","5"), 
   #                  values=c("#522EDB", "#A33C33", "#C4847B","#600B00","#AB886C")) +
  stat_density2d(color='black', h = 6, alpha=I(0.5), size=I(0.5))

plot_cell_trajectory(cde_NonAdv, markers = c("Il6"), use_color_gradient = TRUE, cell_size = 0.25, show_branch_points = FALSE)
plot_cell_trajectory(cde_NonAdv, markers = c("Apoe","Cxcl1"), use_color_gradient = TRUE, cell_size = 0.25, show_branch_points = FALSE)



plot_cell_trajectory(cde_NonAdv, color_by = "State", markers = c("Cdk1","Inmt","Cxcl13","Cxcl1","Ccl19","Il6"), use_color_gradient = TRUE)
plot_cell_trajectory(cde_NonAdv, markers = c("Cdk1","Aldh1a2","Cxcl13","Cd34","Acta2"), use_color_gradient = TRUE, cell_size = 0.5)
plot_cell_trajectory(cde_NonAdv, markers = c("Bst1","Pdpn","Pdgfrb","Cd34","Il6","Cxcl1"), use_color_gradient = TRUE, cell_size = 0.25)
nrow(subset(pData(cde_NonAdv)))
plot_complex_cell_trajectory(cde_NonAdv, color_by = 'State', show_branch_points = T,
                             cell_size = 0.5, cell_link_size = 0.3, root_states = c(1))  + facet_wrap(~Time_point)
plot_complex_cell_trajectory(cde_NonAdv, color_by = 'State', show_branch_points = T,
                             cell_size = 0.5, cell_link_size = 0.3, root_states = c(1))  + facet_wrap(~Cluster_Seurat)

#Save RDS
saveRDS(cde_NonAdv, file=paste(PATH_CDS_objects,"/day0_10_24_56_300_Seurat_1500_1000_1_12_SC_minus_feeder1_2_3_PvC_only_NonAdventi_.Rds",sep=""))

##### Adv
cde_Adv<- cde_SC[,cells_Adv]
# Clustering
cde_Adv <- reduceDimension(cde_Adv, max_components = 2, num_dim = 12,
                           reduction_method = 'tSNE',
                           residualModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression + Time_point",
                           verbose = T)
cde_Adv <- clusterCells(cde_Adv, num_clusters = 6)
plot_cell_clusters(cde_Adv, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10)
plot_cell_clusters(cde_Adv, 1, 2, color = "Cluster", show_cell_names = TRUE, cell_name_size = 10) + facet_wrap(~Time_point)
plot_cell_clusters(cde_Adv, 1, 2, color = "State", show_cell_names = TRUE, cell_name_size = 10) + facet_wrap(~Time_point)
plot_cell_clusters(cde_Adv, 1, 2, color = "Cluster", markers = c("Cxcl13")) + facet_wrap(~Time_point)
plot_cell_clusters(cde_Adv, 1, 2, color = "Cluster", markers = c("Cd34","Tnfsf11","Aldh1a2","Cxcl13",
                                                                 "Nfkb1","Ltbr","Icam1","Vcam1","Acta2",
                                                                 "Cd248","Gdf10","Cxcl9","Has1","Ackr3",
                                                                 "Cdk1","Il6"), cell_size = 0.5)

#Calculate % per cluster across age
meta_data <- pData(cde_Adv)
l_condition <- split(meta_data, meta_data$tech)
l_cells_per_condition <- lapply(l_condition, function(x) {
  #x <- l_condition[[1]]
  x_cells <- unlist(lapply(split(x, x$Cluster_Seurat), nrow))
  x_cells
})
names(l_cells_per_condition) <- names(l_condition)
t_cells_per_condition <- do.call(rbind.fill.matrix, lapply(l_cells_per_condition, function(x) {t(as.data.frame(x))}))
t_cells_per_condition[is.na(t_cells_per_condition)] <- 0
rownames(t_cells_per_condition) <- names(l_condition)
#Thresh for subsets with at least defined cell number
min_num_cell_cluster = 30
cells_per_cluster <- colSums(t_cells_per_condition)
clusters_present <- cells_per_cluster[cells_per_cluster > min_num_cell_cluster]
t_cells_per_condition <- t_cells_per_condition[,names(clusters_present)]
#Normalize to 1000 cells
norm_factor_vector <- rowSums(t_cells_per_condition) / 1000
#Divide cell number per cluster by normalization factor
t_cells_per_condition_norm <- t_cells_per_condition / norm_factor_vector
#Calculate Frequencies
total <- colSums(t_cells_per_condition_norm)
freq_cells_per_condition <- (t(t_cells_per_condition_norm) / total) * 100
freq_cells_per_condition <- as.matrix(freq_cells_per_condition)
pheatmap(freq_cells_per_condition, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = TRUE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white","grey","grey23","black"), space="rgb")(128))

#Trajectories
cde_Adv <- reduceDimension(cde_Adv, max_components = 2, method = 'DDRTree')
#Build trajectory
cde_Adv <- orderCells(cde_Adv)
plot_cell_trajectory(cde_Adv, color_by = "Pseudotime", cell_size = 0.25, show_tree = TRUE,
                     show_branch_points = FALSE)
#cde_NonAdv root state
cde_Adv <- orderCells(cde_Adv, root_state = 3)
#Note: Figure for Paper
plot_cell_trajectory(cde_NonAdv, color_by = "State")  + facet_wrap(~Time_point)
plot_cell_trajectory(cde_Adv,
                     color_by = "State",
                     cell_size = 0.25,
                     show_branch_points = FALSE) + 
  facet_wrap(~Time_point) +
  scale_color_manual(breaks = c("1", "2", "3","4","5"), 
                     values=c("#006A8E", "#04C4BB", "#522EDB","grey40","grey")) +
  stat_density2d(color='black', h = 6, alpha=I(0.5), size=I(0.5))
#046EC4


plot_cell_trajectory(cde_Adv,
                     color_by = "State",
                     cell_size = 0.25,
                     show_branch_points = FALSE) + 
  facet_wrap(~Cluster_Seurat) +
  scale_color_manual(breaks = c("1", "2", "3","4","5"), 
                     values=c("#006a8e", "#04C4BB", "#522EDB","grey40","#39DDCA")) +
  stat_density2d(color='black', h = 6, alpha=I(0.5), size=I(0.5))
plot_cell_trajectory(cde_Adv, markers = c("Ptgis","Aldh1a2","Cd248","Flt3l","Col4a1","Ackr3"), use_color_gradient = TRUE, cell_size = 0.25, show_branch_points = FALSE)
plot_cell_trajectory(cde_Adv, markers = c("Aldh1a2"), use_color_gradient = TRUE, cell_size = 0.25, show_branch_points = FALSE)
plot_cell_trajectory(cde_Adv, markers = c("Apoe","Cxcl1","H19","Ma","Cd248","Ackr3"), use_color_gradient = TRUE, cell_size = 0.25, show_branch_points = FALSE)



plot_cell_trajectory(cde_Adv, color_by = "State") + 
  facet_wrap(~Cluster_Seurat) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25))
plot_cell_trajectory(cde_Adv, color_by = "State", markers = c("Cdk1","Ptgis","Aldh1a2","Cd248","Flt3l","Col4a1"), use_color_gradient = TRUE)
plot_cell_trajectory(cde_Adv, markers = c("Cdk1","Aldh1a2","Cxcl13","Cd34","Ackr3","Cd248"), use_color_gradient = TRUE, cell_size = 0.5)
plot_cell_trajectory(cde_Adv, markers = c("Bst1","Pdpn","Pdgfrb","Cd34"), use_color_gradient = TRUE, cell_size = 0.5)
nrow(subset(pData(cde_NonAdv)))
plot_complex_cell_trajectory(cde_Adv, color_by = 'State', show_branch_points = T,
                             cell_size = 0.5, cell_link_size = 0.3, root_states = c(3))  + facet_wrap(~Time_point)
plot_complex_cell_trajectory(cde_Adv, color_by = 'State', show_branch_points = T,
                             cell_size = 0.5, cell_link_size = 0.3, root_states = c(3))  + facet_wrap(~Cluster_Seurat)

#Save RDS
saveRDS(cde_Adv, file=paste(PATH_CDS_objects,"/day0_10_24_56_300_Seurat_1500_1000_1_12_SC_minus_feeder1_2_3_PvC_only_Adventi_.Rds",sep=""))

#Load RDS
cde_NonAdv <- readRDS(paste(PATH_CDS_objects,"/day0_10_24_56_300_Seurat_1500_1000_1_12_SC_minus_feeder1_2_3_PvC_only_NonAdventi_.Rds",sep=""))
cde_Adv <- readRDS(paste(PATH_CDS_objects,"/day0_10_24_56_300_Seurat_1500_1000_1_12_SC_minus_feeder1_2_3_PvC_only_Adventi_.Rds",sep=""))

#####
#Compile Tables for dynGENIE3: per Age
#####

# Select Stages of trajectory
# Adventitial
stages_adventi_cells <- subset(pData(cde_Adv))
# Nonadventitial
stages_nonadventi_cells <- subset(pData(cde_NonAdv))

# Dissect into clusters
time_points <- levels(as.factor(pData(cde_nonPvC)$Time_point))
l_adventi_cells_per_Age <- split(stages_adventi_cells, stages_adventi_cells$Time_point)
names(l_adventi_cells_per_Age) <- time_points
l_nonadventi_cells_per_Age <- split(stages_nonadventi_cells, stages_nonadventi_cells$Time_point)
names(l_nonadventi_cells_per_Age) <- time_points
# Average expression over time point and compile table
#' average_expr
#' average expression across each element of list
#' @param l_cell_substs named list of pData tables containing the cells per subset
#' @param cds_interest CDS object to be sampled from
#' @return t_average_expression table of averaged expression with columns being time-point and rows gene names

f_average_expr <- function(l_cell_subsets, cds_interest){
  #Get gene names
  gene_names <- fData(cde_nonPvC)$gene_short_name
  #initiate empty table
  t_averaged <- data.frame(matrix(ncol = length(names(l_cell_subsets)), nrow = length(gene_names)))
  
  #for each table
  for(i in 1:length(l_cell_subsets)){
    cell_subsets_i <- l_cell_subsets[[i]]
    nrow(cell_subsets_i)
    expr_i <- as.matrix(exprs(cds_interest)[,colnames(exprs(cds_interest)) %in% rownames(cell_subsets_i)])
    expr_mean_i <- rowMeans(expr_i)
    t_averaged[,i] <- expr_mean_i
  }
  colnames(t_averaged) <- names(l_cell_subsets)
  t_averaged <- cbind(gene_names,t_averaged)
  t_averaged <- t_averaged[!duplicated(t_averaged[c("gene_names")]),]
  gene_names_final <- t_averaged[,1]
  t_averaged <- t_averaged[,2:ncol(t_averaged)]
  rownames(t_averaged) <- gene_names_final
  #return
  t_averaged
}

# run function
averaged_adventi_dynGENIE3 <- f_average_expr(l_adventi_cells_per_Age,cde_Adv)
averaged_nonadventi_dynGENIE3 <- f_average_expr(l_nonadventi_cells_per_Age,cde_NonAdv)
subset(averaged_adventi_dynGENIE3, rownames(averaged_adventi_dynGENIE3) %in% c("Cd34"))
subset(averaged_nonadventi_dynGENIE3, rownames(averaged_nonadventi_dynGENIE3) %in% c("Cd34"))
#export
write.table(averaged_adventi_dynGENIE3, paste(PATH_output,"/","averaged_adventi_dynGENIE3.txt",sep=""), dec = ".", sep = "\t")
write.table(averaged_adventi_dynGENIE3, paste(PATH_output,"/","averaged_nonadventi_dynGENIE3.txt",sep=""), dec = ".", sep = "\t")

#####
#Compile Tables for dynGENIE3: per Age and State
#####

# Select Stages of trajectory
# Adventitial
stages_adventi_cells <- subset(pData(cde_Adv))
# Nonadventitial
stages_nonadventi_cells <- subset(pData(cde_NonAdv))

# Dissect into clusters
order_adv_a <- c("1","5","2","4")
order_adv_b <- c("1","5","2","3")
#Note: Set dataset
stages_adventi_cells <- stages_nonadventi_cells
cds_of_interest <- cde_NonAdv
#make list by state
l_stages_adventi_cells <- split(stages_adventi_cells, stages_adventi_cells$State)
#make list of state list by time point
l_l_stages_ages_adventi_cells <- lapply(l_stages_adventi_cells, function(x){split(x, x$Time_point)})
#Order by state obtained by trajectory analysis
l_l_stages_ages_adventi_cells_a <- l_l_stages_ages_adventi_cells[order_adv_a]
#Store tables in list to calculate averaged RCMs
l_stages_ages_adventi_cells_a <- list()
#Note: Input required
# If cell number is very low it makes no sense to average the expression, Some states are also not present at certain ages
min_cells_per_State_AND_TimePoint = 10
#index
k = 1
for(i in 1:length(l_l_stages_ages_adventi_cells_a)){
  #take state
  names_i <- names(l_l_stages_ages_adventi_cells_a)[i]
  l_stages_ages_i <- l_l_stages_ages_adventi_cells_a[[i]]
  for(j in 1:length(l_stages_ages_i)){
    #per state check number of cells per Timepoint
    names_j <- names(l_stages_ages_i)[j]
    stages_ages_j <- l_stages_ages_i[[j]]
    if(nrow(stages_ages_j) >= min_cells_per_State_AND_TimePoint){
      name_i_j <- paste(names_i,names_j, sep="_")
      print(name_i_j)
      l_stages_ages_adventi_cells_a[[k]] <- stages_ages_j
      names(l_stages_ages_adventi_cells_a)[k] <- name_i_j
      k = k + 1
    }
  }
}
l_stages_ages_adventi_cells_a <- l_stages_ages_adventi_cells_a[1:(k-1)]


#l_nonadventi_cells_per_Age <- split(stages_nonadventi_cells, stages_nonadventi_cells$Time_point)
#names(l_nonadventi_cells_per_Age) <- datasets
# Average expression over time point and compile table
#' average_expr
#' average expression across each element of list
#' @param l_cell_substs named list of pData tables containing the cells per subset
#' @param cds_interest CDS object to be sampled from
#' @return t_average_expression table of averaged expression with columns being time-point and rows gene names

f_average_expr <- function(l_cell_subsets, cds_interest){
  #Get gene names
  gene_names <- fData(cds_interest)$gene_short_name
  #initiate empty table
  t_averaged <- data.frame(matrix(ncol = length(names(l_cell_subsets)), nrow = length(gene_names)))
  
  #for each table
  for(i in 1:length(l_cell_subsets)){
    cell_subsets_i <- l_cell_subsets[[i]]
    nrow(cell_subsets_i)
    expr_i <- as.matrix(exprs(cds_interest)[,colnames(exprs(cds_interest)) %in% rownames(cell_subsets_i)])
    expr_mean_i <- rowMeans(expr_i)
    t_averaged[,i] <- expr_mean_i
  }
  colnames(t_averaged) <- names(l_cell_subsets)
  t_averaged <- cbind(gene_names,t_averaged)
  t_averaged <- t_averaged[!duplicated(t_averaged[c("gene_names")]),]
  gene_names_final <- t_averaged[,1]
  t_averaged <- t_averaged[,2:ncol(t_averaged)]
  rownames(t_averaged) <- gene_names_final
  #return
  t_averaged
}

# run function
averaged_adventi_dynGENIE3 <- f_average_expr(l_stages_ages_adventi_cells_a,cds_of_interest)
subset(averaged_adventi_dynGENIE3, rownames(averaged_adventi_dynGENIE3) %in% c("Cd34"))

#export
write.table(averaged_adventi_dynGENIE3, paste(PATH_output,"/","averaged_adventi_State321_AND_TIME_dynGENIE3.txt",sep=""), dec = ".", sep = "\t")

write.table(averaged_adventi_dynGENIE3, paste(PATH_output,"/","averaged_nonadventi_State1524_AND_TIME_dynGENIE3.txt",sep=""), dec = ".", sep = "\t")

write.table(averaged_adventi_dynGENIE3, paste(PATH_output,"/","averaged_nonadventi_dynGENIE3.txt",sep=""), dec = ".", sep = "\t")

#####
#Signature identification across subsets
#####
cde_to_check <- cde_NonAdv
GeneSymbol <- fData(cde_to_check)$gene_short_name
#subset expression matrix and perform colmeans per cell
Expression <- as.matrix(exprs(cde_to_check))
#Add GeneSymbol as rownames
rownames(Expression) <- GeneSymbol
#For each cluster calculate the mean expression across all genes
l_meaned_expression <- list()
for(i in seq(length(levels(pData(cde_to_check)$Cluster)))){
  cluster_ID_i <- levels(pData(cde_to_check)$Cluster)[i]
  print(cluster_ID_i) 
  cells_cluster_i <- rownames(subset(pData(cde_to_check), Cluster == cluster_ID_i))
  Expression_mean_i <- rowMeans(Expression[,colnames(Expression) %in% cells_cluster_i])
  l_meaned_expression[[i]] <- Expression_mean_i
}
#Prep table for meaned expression per cluster
t_meaned_expression <- as.data.frame(do.call("cbind", l_meaned_expression))
t_meaned_expression <- cbind(GeneSymbol, t_meaned_expression)
colnames(t_meaned_expression) <- c("GeneSymbol", levels(pData(cde_to_check)$Cluster))

#Load DEG tables and make list of clusters
DEGs_core_40 <- NC_Pezoldt_Pasztoi_2018
DEGs_list <- split(DEGs_core_40, DEGs_core_40$cluster)


i = 1
out = NULL
for(i in i:length(DEGs_list)){
  cluster_DEG_i <- as.character(DEGs_list[[i]]$gene)
  cluster_DEG_n_i <- length(cluster_DEG_i)
  
  cluster_ID_i <- names(DEGs_list)[i]
  y_label_i <- paste("cZ: ", cluster_ID_i)
  
  print(cluster_ID_i)
  print(i)
  print(cluster_DEG_n_i)
  
  aExp_cluster_i <- subset(t_meaned_expression, rownames(t_meaned_expression) %in% cluster_DEG_i)[,2:ncol(t_meaned_expression)]
  
  Scale_aExp_cluster_i <- apply(aExp_cluster_i, 1, scale)
  Zscore_cluster_i <- rowSums(Scale_aExp_cluster_i)
  Zscore_cluster_i_norm <- Zscore_cluster_i / (nrow(aExp_cluster_i)/10)
  
  #mLN maintained Scores
  Zscore_table_cluster_i <- data.frame(module = c(rep(cluster_ID_i,length(Zscore_cluster_i_norm))),
                                       cluster = colnames(t_meaned_expression)[2:ncol(t_meaned_expression)],
                                       cZscore = Zscore_cluster_i_norm)
  
  
  #Generate cZscore bargraph matrix
  p <- ggplot(data = Zscore_table_cluster_i, aes(x = cluster, y = cZscore, fill = module)) +
    geom_bar(stat = "identity", colour = "black") +
    theme(axis.text.x=element_text(angle=270,hjust=1), legend.position="none") +
    scale_fill_manual(values = c("deepskyblue1")) +
    ylab(y_label_i)
  #print(p)
  out[[i]] <- p
}

title <- paste(condition,"DEGs_core_40",sep = "_")
pdf(paste(PATH_output,"/day0_10_24_56_300_nonPvC_NC_clusters.pdf", sep=""))
marrangeGrob(out, nrow = round(length(DEGs_list) / 3 - 1), ncol = 3, top = title)
dev.off()

#####
#Percentage of cell types across Stages
#####
# Cell distribution over tree structure
cde_x <- cde_Adv
state_cluster_stat <- table(pData(cde_x)[, c('State', 'Cluster_Seurat')])
#Thresh for subsets with at least defined cell number
min_num_cell_cluster = 100
cells_per_cluster <- colSums(state_cluster_stat)
clusters_present <- cells_per_cluster[cells_per_cluster > 50]
state_cluster_stat <- state_cluster_stat[,names(clusters_present)]

rowSums(state_cluster_stat)
#Normalize
state_cluster_stat <- apply(state_cluster_stat, 2, function(x) x / sum(x))
state_cluster_stat_ordered <- t(state_cluster_stat)

options(repr.plot.width=3, repr.plot.height=3)
pheatmap(state_cluster_stat_ordered, cluster_rows = TRUE, legend = TRUE,
         treeheight_row = 0, treeheight_col = 0, show_rownames = TRUE, cluster_cols = FALSE,
         scale = "none", border_color = "black", cellwidth = 10,
         cellheigth = 10, color = colorRampPalette(c("white","grey","grey23","black"), space="rgb")(128))

#####
#Percentage of cells per cluster
#####
#Adv n
state_cluster_stat_Adv <- table(pData(cde_Adv)[, c('State', 'Cluster_Seurat')])
min_num_cell_cluster = 100
cells_per_cluster_Adv <- colSums(state_cluster_stat_Adv)
#NonAdv n 
state_cluster_stat_NonAdv <- table(pData(cde_NonAdv)[, c('State', 'Cluster_Seurat')])
min_num_cell_cluster = 100
cells_per_cluster_NonAdv <- colSums(state_cluster_stat_NonAdv)
#Table numbers
total_Numbers <- cells_per_cluster_Adv + cells_per_cluster_NonAdv
freq_NonAdv <-
freq_Adv <- round(cells_per_cluster_Adv / total_Numbers * 100, 2)
freq_NonAdv <- round(cells_per_cluster_NonAdv / total_Numbers * 100, 2)
toPlot <- data.frame(subsets = names(freq_Adv),
                     percentage = c(freq_Adv,freq_NonAdv),
                     lineage = c(rep("CD34FSC",length(freq_Adv)),
                                 rep("FRC",length(freq_NonAdv))))


#order
toPlot$subsets <- factor(toPlot$subsets,levels = c("CD34+Aldh1a2+", "CD34+Ackr3+", "CD34+CD248+", "SCx",
                                         "Cxcl9+","Il6+Cxcl1+",
                                         "Cdk1+","LTolike",
                                         "Ccl19+Il7+","Inmt+","Inmt+Cxcl12+","pSC"))


library(ggplot2)

# Small multiple
ggplot(toPlot, aes(fill=lineage, y=percentage, x=subsets)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Percentage of Cell per Trajectory") +
  scale_fill_manual(values = c("#006a8e", "#600B00")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  xlab("Subset")
  #AB886C", "#A33C33", "#C4847B","#600B00","#522EDB"
  #006a8e", "#04C4BB", "#522EDB","grey40","#39DDCA"
#################
#BEAM
#################
### Load Genes of interest
#Riken TF
TF_table <- read.delim(PATH_input_TFdb_Riken)
Riken_TFs <- as.character(TF_table$Symbol)
# ATAC-seq TFs
TF_DARs <- read.table(PATH_input_TFBS_DARs)
DAR_TFs <- rownames(TF_DARs)
# Secreted extracellular UniProt
Secreted_uniprot <- read.delim(PATH_input_secreted_proteins)
Secreted_uniprot_GeneNames <- unique(Secreted_uniprot$Gene.names...primary..)
# DMR TFs
DMR_TFs <- row.names(read.csv(PATH_input_DMRs))


#Perform BEAM for all branches
#####
#BEAM adventi
#####
#List to store BEAMs
cde_x <- cde_Adv
l_BEAM_res <- list()
for(i in 1:2){
  BEAM_res_i <- BEAM(cde_x, branch_point = i, cores =12)
  l_BEAM_res[[i]] <- BEAM_res_i
}
l_BEAM_res_Adv <- readRDS(paste(PATH_output,"/","BEAM_day0_10_24_56_300_Seurat_1500_1000_1_12_SC_minus_feeder1_2_3_PvC_only_Adventi.rds",sep=""))

#Eliminate ribosomal genes
rpl.genes <- grep(pattern = "^Rpl", x = l_BEAM_res_Adv[[2]]$gene_short_name, value = TRUE)
rps.genes <- grep(pattern = "^Rps", x = l_BEAM_res_Adv[[2]]$gene_short_name, value = TRUE)
ribo.genes <- c(rpl.genes, rps.genes)
#get ensemble Ids
BEAM_res_minus_ribo_2 <- subset(l_BEAM_res_Adv[[2]], !(gene_short_name %in% ribo.genes))
BEAM_res_minus_ribo_TFs_2 <- subset(BEAM_res_minus_ribo_2, gene_short_name %in% Riken_TFs)
BEAM_res_minus_ribo_DARsTFs_Adv_2 <- subset(BEAM_res_minus_ribo_2, gene_short_name %in% DAR_TFs)
BEAM_res_minus_ribo_DMRsTFs_2 <- subset(BEAM_res_minus_ribo_2, gene_short_name %in% DMR_TFs)
BEAM_res_minus_ribo_Secreted_2 <- subset(BEAM_res_minus_ribo_2, gene_short_name %in% Secreted_uniprot_GeneNames)

#Continue here:
# 1) Table of DAR significant in BEAM branchpoint 2
DAR_TFs_sig_in_BEAM_Adv <- subset(BEAM_res_minus_ribo_DARsTFs_Adv_2,qval < 5*1e-2)
write.table(DAR_TFs_sig_in_BEAM_Adv, paste(PATH_output,"/Figures/BEAM/Adv_DAR_TFs_sig_in_BEAM.txt",sep=""), sep = "\t")
# 2) Table of DMR significant in BEAM branchpoint 2
DMRs_TFs_sig_in_BEAM <- subset(BEAM_res_minus_ribo_DMRsTFs_2,qval < 5*1e-2)
write.table(DMRs_TFs_sig_in_BEAM, paste(PATH_output,"/Figures/BEAM/Adv_DMRs_TFs_sig_in_BEAM.txt",sep=""), sep = "\t")


plot_genes_branched_heatmap(cde_x[row.names(subset(BEAM_res_minus_ribo_TFs_2,qval < 5*1e-2)),],
                            branch_point = 2,
                            #branch_states = c(1,5),
                            num_clusters = 5,
                            #length(levels(as.factor(pData(cde_x)$Cluster))),
                            cores = 12,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_colors = c("#522EDB", "#04C4BB", "#046EC4"))#04C4BB#049EC4

plot_cell_trajectory(cde_Adv, markers = c("Egr4","Jund","Nr4a1","Hsf2",
                                          "Hoxa10","Irf8",
                                          "Sox11","Notch3",
                                          "Nkx2-3"), use_color_gradient = TRUE, cell_size = 0.25, show_branch_points = FALSE)


plot_cell_trajectory(cde_Adv, markers = c("Il34","Bmp4","Cxcl10","Bmp2",
                                          "Wnt10b","Col4a3","Cxcl2","Ngf",
                                          "Mmp3","Mmp2","Fgf7","Vtn",
                                          "Col4a1","Col15a1","Wnt5a","Enpp2",
                                          "Gdf10","Isg10","Ackr3","Aldh1a3",
                                          "Tgfbr2",
                                          "Islr",
                                          "Ptgis",
                                          "Ptgs1",
                                          "Tnxb",
                                          "Ogn"), use_color_gradient = TRUE, cell_size = 0.25, show_branch_points = FALSE)

gene_list_map_on_trajectory <- c(#Up
  "Ackr3",
  "Aldh1a3",
  "Tgfbr2",
  "Ptgis",
  "Ptgs1",
  #Down
  "Col15a1",
  "Gdf10",
  "Vtn",
  "Wnt5a",
  "Fgf7",
  "Cxcl9",
  "Stat1",
  "Atf3","Jun","Egr2","Jund","Irf1","Klf9","Egr1")
for(i in 1:length(gene_list_map_on_trajectory)){
  setwd(paste(PATH_output,"/Figures/Trajectories_advFSC",sep=""))
  gene_i <- gene_list_map_on_trajectory[i]
  print(gene_i)
  filename_i <- paste("advFSC_", gene_i, ".png", sep="")
  print(filename_i)
  plot_cell_trajectory(cde_Adv, markers = c(gene_i), use_color_gradient = TRUE,
                       cell_size = 0.25, show_branch_points = FALSE)
  
  ggsave(filename_i, width = 5.0, height = 7.9, units = "cm")
}
for(i in 1:length(gene_list_map_on_trajectory)){
  # i = 1
  setwd(paste(PATH_output,"/Figures/Trajectories_advFSC",sep=""))
  gene_i <- gene_list_map_on_trajectory[i]
  print(gene_i)
  filename_i <- paste("advFSC_", gene_i, ".eps", sep="")
  print(filename_i)
  plot_cell_trajectory(cde_Adv, markers = c(gene_i), use_color_gradient = TRUE,
                       cell_size = 0.25, show_branch_points = FALSE)
  
  ggsave(filename_i, width = 5.0, height = 7.9, units = "cm")
}
for(i in 1:length(gene_list_map_on_trajectory)){
  # i = 1
  setwd(paste(PATH_output,"/Figures/Trajectories_advFSC",sep=""))
  gene_i <- gene_list_map_on_trajectory[i]
  print(gene_i)
  filename_i <- paste("advFSC_", gene_i, "_Scale", ".eps", sep="")
  print(filename_i)
  plot_cell_trajectory(cde_Adv, markers = c(gene_i), use_color_gradient = TRUE,
                       cell_size = 0.25, show_branch_points = FALSE)
  
  ggsave(filename_i, width = 7.0, height = 7.9, units = "cm")
}


#####
#BEAM Nonadventi
#####
#List to store BEAMs
cde_x <- cde_NonAdv
l_BEAM_res <- list()
for(i in 1:2){
  BEAM_res_i <- BEAM(cde_x, branch_point = i, cores =12)
  l_BEAM_res[[i]] <- BEAM_res_i
}
l_BEAM_res_NonAdv <- readRDS(paste(PATH_output,"/","BEAM_day0_10_24_56_300_Seurat_1500_1000_1_12_SC_minus_feeder1_2_3_PvC_only_NonAdventi.rds",sep=""))

#Eliminate ribosomal genes
rpl.genes <- grep(pattern = "^Rpl", x = l_BEAM_res_NonAdv[[2]]$gene_short_name, value = TRUE)
rps.genes <- grep(pattern = "^Rps", x = l_BEAM_res_NonAdv[[2]]$gene_short_name, value = TRUE)
ribo.genes <- c(rpl.genes, rps.genes)
#get ensemble Ids
BEAM_res_minus_ribo_2 <- subset(l_BEAM_res_NonAdv[[2]], !(gene_short_name %in% ribo.genes))
BEAM_res_minus_ribo_TFs_2 <- subset(BEAM_res_minus_ribo_2, gene_short_name %in% Riken_TFs)
BEAM_res_minus_ribo_DARsTFs_NonAdv_2 <- subset(BEAM_res_minus_ribo_2, gene_short_name %in% DAR_TFs)
BEAM_res_minus_ribo_DMRsTFs_2 <- subset(BEAM_res_minus_ribo_2, gene_short_name %in% DMR_TFs)
BEAM_res_minus_ribo_Secreted_2 <- subset(BEAM_res_minus_ribo_2, gene_short_name %in% Secreted_uniprot_GeneNames)

#Continue here:
# 1) Table of DAR significant in BEAM branchpoint 2
DAR_TFs_sig_in_BEAM_NonAdv <- subset(BEAM_res_minus_ribo_DARsTFs_NonAdv_2,qval < 5*1e-2)
write.table(DAR_TFs_sig_in_BEAM_NonAdv, paste(PATH_output,"/Figures/BEAM/NonAdv_DAR_TFs_sig_in_BEAM.txt",sep=""), sep = "\t")
# 2) Table of DMR significant in BEAM branchpoint 2
DMRs_TFs_sig_in_BEAM <- subset(BEAM_res_minus_ribo_DARsTFs_NonAdv_2,qval < 5*1e-2)
write.table(DMRs_TFs_sig_in_BEAM, paste(PATH_output,"/Figures/BEAM/NonAdv_DMRs_TFs_sig_in_BEAM.txt",sep=""), sep = "\t")


plot_genes_branched_heatmap(cde_x[row.names(subset(BEAM_res_minus_ribo_Secreted_2,qval < 5*1e-10)),],
                            branch_point = 2,
                            #branch_states = c(1,5),
                            num_clusters = 7,
                            #length(levels(as.factor(pData(cde_x)$Cluster))),
                            cores = 12,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_colors = c("#522EDB", "#C4847B","#A33C33"))#04C4BB#049EC4


plot_cell_trajectory(cde_NonAdv, markers = c("Egr4","Jund","Nr4a1","Hsf2",
                                             "Hoxa10","Irf8",
                                             "Sox11","Notch3",
                                             "Nkx2-3"), use_color_gradient = TRUE, cell_size = 0.25, show_branch_points = FALSE)


plot_cell_trajectory(cde_NonAdv, markers = c(#Il6Cxcl1
  "Cxcl14","Cxcl5","Cxcl11","Cxcl1","Cxcl9",
  "Il6","Ccl7",
  #Ccl19Il7
  "Il33","Ccl19",
  "Col4a2","Il7","Stat1",
  "Tgfbi",
  "Bmp3",
  "Tnfsf13b",
  "Ccl21b",
  "Fmod", "Il34",
  "Bmp4"), use_color_gradient = TRUE, cell_size = 0.25, show_branch_points = FALSE)

gene_list_map_on_trajectory <- c("Cxcl1",
                                 "Ccl19",
                                 "Il33",
                                 "Tgfbi",
                                 "Tnfsf13b",
                                 "Il6",
                                 "Cxcl5",
                                 "Cxcl2",
                                 "Cxcl10",
                                 "Ccl2",
                                 "Ccl7",
                                 "Cxcl9",
                                 "Stat1",
                                 "Atf3","Jun","Egr2","Jund","Irf1","Klf9","Egr1")
for(i in 1:length(gene_list_map_on_trajectory)){
  # i = 1
  setwd(paste(PATH_output,"/Figures/Trajectories_retFSC",sep=""))
  gene_i <- gene_list_map_on_trajectory[i]
  print(gene_i)
  filename_i <- paste("retFSC_", gene_i, ".png", sep="")
  print(filename_i)
  plot_cell_trajectory(cde_NonAdv, markers = c(gene_i), use_color_gradient = TRUE,
                       cell_size = 0.25, show_branch_points = FALSE)
  
  ggsave(filename_i, width = 5.0, height = 7.9, units = "cm")
}
for(i in 1:length(gene_list_map_on_trajectory)){
  # i = 1
  setwd(paste(PATH_output,"/Figures/Trajectories_retFSC",sep=""))
  gene_i <- gene_list_map_on_trajectory[i]
  print(gene_i)
  filename_i <- paste("retFSC_", gene_i, ".eps", sep="")
  print(filename_i)
  plot_cell_trajectory(cde_NonAdv, markers = c(gene_i), use_color_gradient = TRUE,
                       cell_size = 0.25, show_branch_points = FALSE)
  
  ggsave(filename_i, width = 5.0, height = 7.9, units = "cm")
}
for(i in 1:length(gene_list_map_on_trajectory)){
  # i = 1
  setwd(paste(PATH_output,"/Figures/Trajectories_retFSC",sep=""))
  gene_i <- gene_list_map_on_trajectory[i]
  print(gene_i)
  filename_i <- paste("retFSC_", gene_i, "_Scale", ".eps", sep="")
  print(filename_i)
  plot_cell_trajectory(cde_NonAdv, markers = c(gene_i), use_color_gradient = TRUE,
                       cell_size = 0.25, show_branch_points = FALSE)
  
  ggsave(filename_i, width = 7.0, height = 7.9, units = "cm")
}




#Overlap TFs Adv and NonAdv
# DAR
BEAM_res_minus_ribo_DARsTFs_Adv_2

overlap <- subset(DAR_TFs_sig_in_BEAM_Adv, gene_short_name %in% DAR_TFs_sig_in_BEAM_NonAdv$gene_short_name)
unique_TF_DAR_Adv <- subset(DAR_TFs_sig_in_BEAM_Adv, !(gene_short_name %in% DAR_TFs_sig_in_BEAM_NonAdv$gene_short_name))
unique_TF_DAR_NonAdv <- subset(DAR_TFs_sig_in_BEAM_NonAdv, !(gene_short_name %in% DAR_TFs_sig_in_BEAM_Adv$gene_short_name))
print(paste("n Total Adv:", nrow(DAR_TFs_sig_in_BEAM_Adv)))
print(paste("n Total NonAdv:", nrow(DAR_TFs_sig_in_BEAM_NonAdv)))

print(paste("n Unique Adv:", nrow(unique_TF_DAR_Adv)))
print(paste("n Unique NonAdv:", nrow(unique_TF_DAR_NonAdv)))















# Plotting stuff
options(repr.plot.width=6, repr.plot.height=6)
plot_multiple_branches_pseudotime(cde_x[c(TFs_Branching_from_DARs),],
                                  branches=c(3,4,5), 
                                  color_by = 'Branch',
                                  branches_name=c(3,4,5), nrow = 8, ncol = 5)

options(repr.plot.width=8, repr.plot.height=12)
plot_multiple_branches_heatmap(cde_x[c(TFs_Branching_from_DARs),],
                               branches=c(1,4,5),
                               branches_name=c(1,4,5),
                               show_rownames=T,
                               num_clusters=length(levels(as.factor(pData(cde_x)$Cluster))))

saveRDS(l_BEAM_res, paste(PATH_output,"/","BEAM_day0_10_24_56_300_Seurat_1500_1000_1_12_SC_minus_feeder1_2_3_PvC_only_Adventi.rds",sep=""))
#l_BEAM_res_NonAdv <- readRDS(paste(PATH_output,"/","BEAM_day0_10_24_56_300_Seurat_1500_1000_1_12_SC_minus_feeder1_3_PvC_only_NonAdventi.rds",sep=""))
#l_BEAM_res_Adv <- readRDS(paste(PATH_output,"/","BEAM_day0_10_24_56_300_Seurat_1500_1000_1_12_SC_minus_feeder1_3_PvC_only_Adventi.rds",sep=""))


#Identify TFs at branching points
# Use TF list from Genomatix
Genomatix_murineTFs <- read.table("/home/pezoldt/NAS2/pezoldt/Analysis/WGBS/FSC/Genomatix_TF/ANNa/Mouse/Genomatix_murineTFs_OnlyNames.txt",
                                  sep = "\t", skip = 1)
# DMR TFs
mLN_hypo <- read.table("/home/pezoldt/NAS2/pezoldt/Analysis/WGBS/FSC/Genomatix_TF/ANNa/Mouse/To_Use/Output/Compiled_TFs_Genomatix_hypo_all.txt")
mLN_hypo_TF <- as.character(mLN_hypo$V1)
mLN_hyper <- read.table("/home/pezoldt/NAS2/pezoldt/Analysis/WGBS/FSC/Genomatix_TF/ANNa/Mouse/To_Use/Output/Compiled_TFs_Genomatix_hyper_all.txt")
mLN_hyper_TF <- as.character(mLN_hyper$V1)

# DAR TFs existence matrix
DAR_existence_matrix <- read.table("/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/SPF/SPF_existence_matrix_known_homer_compile.txt")
Unique_pLN_peak_UP <- rownames(subset(DAR_existence_matrix, 
                                      pLN_Open_UP == 0 &
                                        mLN_Open_UP == 0 &
                                        pLN_peak_UP == 1 &
                                        mLN_peak_UP == 0 &
                                        pLN_Open_None == 0 &
                                        mLN_Open_None == 0))
Unique_mLN_Open_None <- rownames(subset(DAR_existence_matrix, 
                                        pLN_Open_UP == 0 &
                                          mLN_Open_UP == 0 &
                                          pLN_peak_UP == 0 &
                                          mLN_peak_UP == 0 &
                                          pLN_Open_None == 0 &
                                          mLN_Open_None == 1))

Common_used <- rownames(subset(DAR_existence_matrix, 
                               !(pLN_Open_UP == 0 &
                                   mLN_Open_UP == 0 &
                                   pLN_peak_UP == 0 &
                                   mLN_peak_UP == 0 &
                                   pLN_Open_None == 0 &
                                   mLN_Open_None == 1) |
                                 !(pLN_Open_UP == 0 &
                                     mLN_Open_UP == 0 &
                                     pLN_peak_UP == 1 &
                                     mLN_peak_UP == 0 &
                                     pLN_Open_None == 0 &
                                     mLN_Open_None == 0)))

#Eliminate ribosomal genes from branch defining list
TFs_genomatix <- Genomatix_murineTFs$V2
#get ensemble Ids
l_BEAM_res <- l_BEAM_res_Adv
branch_choice_index <- 2
BEAM_TFs_genomatix_2 <- subset(l_BEAM_res[[2]], gene_short_name %in% TFs_genomatix)
BEAM_mLN_hypo_TF <- subset(l_BEAM_res[[2]], gene_short_name %in% mLN_hypo_TF)
BEAM_mLN_hyper_TF <- subset(l_BEAM_res[[2]], gene_short_name %in% mLN_hyper_TF)
BEAM_mLN_Unique_pLN_peak_UP <- subset(l_BEAM_res[[2]], gene_short_name %in% Unique_pLN_peak_UP)
BEAM_mLN_Unique_mLN_Open_None <- subset(l_BEAM_res[[2]], gene_short_name %in% Unique_mLN_Open_None)
BEAM_Common_pLN_peak_UP_mLN_Open_None <- subset(l_BEAM_res[[2]], gene_short_name %in% Common_pLN_peak_UP_mLN_Open_None)
BEAM_Common_used <- subset(l_BEAM_res[[2]], gene_short_name %in% Common_used)
plot_genes_branched_heatmap(cde_x[row.names(subset(BEAM_mLN_hypo_TF,
                                                   qval < 5*1e-2)),],
                            branch_point = 2,
                            #num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_x[row.names(subset(BEAM_mLN_hyper_TF,
                                                   qval < 5*1e-2)),],
                            branch_point = 2,
                            #num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_x[row.names(subset(BEAM_mLN_Unique_pLN_peak_UP,
                                                   qval < 5*1e-2)),],
                            branch_point = 2,
                            #num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_x[row.names(subset(BEAM_mLN_Unique_mLN_Open_None,
                                                   qval < 5*1e-2)),],
                            branch_point = 2,
                            #num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_x[row.names(subset(Common_used,
                                                   qval < 5*1e-2)),],
                            branch_point = 2,
                            #num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_x[row.names(subset(BEAM_Common_used,
                                                   qval < 5*1e-2)),],
                            branch_point = 2,
                            #num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))




test <- rbind(BEAM_mLN_Unique_mLN_Open_None,BEAM_mLN_Unique_pLN_peak_UP)

plot_genes_branched_heatmap(cde_x[row.names(subset(test,
                                                   qval < 5*1e-2)),],
                            branch_point = 2,
                            #num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_x[row.names(BEAM_mLN_hyper_TF),],
                            branch_point = 2,
                            num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_x[row.names(BEAM_Common_pLN_peak_UP_mLN_Open_None),],
                            branch_point = 2,
                            num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T,
                            branch_labels = c("NonAdv", "Adv"))

plot_genes_branched_heatmap(cde_x[row.names(subset(l_BEAM_res[[3]],
                                                   qval < 1e-10)),],
                            branch_point = 1,
                            num_clusters = length(levels(as.factor(pData(cde_nonPvC)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = F,
                            branch_labels = c("NonAdv", "Adv"))

#####
#Plot Genes of TF screen
#####
TFs_screen <- read.delim("/home/pezoldt/NAS2/pezoldt/0_Experiments/088_TF_overexpression_C3H10/Analysis/RCM/ID_to_Vector.txt")
TFs_screen_genes <- as.character(TFs_screen$Vector)
#Adv:::Plot Trajectories to check states---------------
plot_cell_trajectory(cde_Adv, color_by = "State") +
  facet_wrap(~Cluster_Seurat, nrow = 2)
clusters_adv <- c("CD34+Ackr3+","CD34+Aldh1a2+","CD34+CD248+",
                  "Cdk1+")
state_CD34Aldh1a2 <- c("2","3","5")
state_CD34CD238 <- c("2","3","1")
#Adv:::Plot Trajectories to check states---------------
plot_cell_trajectory(cde_NonAdv, color_by = "State") +
  facet_wrap(~Cluster_Seurat, nrow = 2)
#Trajectory to Cd34+Alch1a2+
# 3, 2, 5
clusters_nonadv <- c("Ccl19+Il7+","Il6+Cxcl1+",
                  "Inmt+","Inmt+Cxcl12+","LTolike",
                  "Cdk1+","Cxcl9+")
state_InmtCcl19 <- c("1","2","3")
state_FRCIl6Cxcl1 <- c("1","2","4")


levels(as.factor(pData(cde_Adv)$Cluster_Seurat))
#Check out states for given sub-trajectory-------------
cds_x <- cde_NonAdv
genes_x <- TFs_screen_genes
clusters_x <- clusters_nonadv
states_x <- state_InmtCcl19


#Get genes of interest
cds_x_expressed_genes <-  row.names(subset(fData(cds_x),
                                           num_cells_expressed >= 5))
cds_x_filtered <- cds_x[cds_x_expressed_genes,]
my_genes <- row.names(subset(fData(cds_x_filtered),
                             gene_short_name %in% genes_x))
cds_subset <- cds_x_filtered[my_genes,]
#Get cells of given trajectory
cds_subset_cells <- row.names(subset(pData(cds_subset),
                                     State %in% states_x &
                                     Cluster_Seurat %in% clusters_x))
cds_subset<- cds_subset[,cds_subset_cells]
#Plot
plot_genes_in_pseudotime(cds_subset, color_by = "Cluster_Seurat",
                         ncol = 3)

#####
#Plot Genes expression
#####
#set cde of interest
cde_x <- readRDS(file=paste(PATH_CDS_objects,condition,"/CDS_objects","/day0_10_24_56_300_nonLTO.Rds",sep=""))
cde_x <- cde_NonAdv
#Jitter-----------------------------
blast_genes <- row.names(subset(fData(cde_x),
                                gene_short_name %in% c("Nkx2-3", "Ccl9", "Cd34","Gdf10","Ackr3", "Cxcl9")))
blast_genes <- row.names(subset(fData(cde_x),
                                gene_short_name %in% c("Ccnb2", "Myod1", "Myog","Cdk1")))

plot_genes_jitter(cde_x[blast_genes,],
                  grouping = "State",
                  color_by = "Cluster",
                  min_expr = 0.05)

#Plot gene expression in pseudotime
cde_x_expressed_genes <-  row.names(subset(fData(cde_x),
                                           num_cells_expressed >= 10))
cde_x_filtered <- cde_x[cde_x_expressed_genes,]
my_genes <- row.names(subset(fData(cde_x_filtered),
                             gene_short_name %in% c("Aldh1a2","Cd34","Ccl19","Ackr3")))
cde_x_subset <- cde_x_filtered[my_genes,]
plot_genes_in_pseudotime(cde_x_subset, color_by = "State")


#####
#DEGs across clusters
#####
cde_x_cluster_DEG <- differentialGeneTest(cde_x[cde_x_expressed_genes,],
                                          fullModelFormulaStr = "~Total_mRNAs + mito_expression + ribo_expression + Time_point",
                                          cores = 8)

#order
cde_x <- orderCells(cde_x, root_state = 4)

#
plot_cell_trajectory(cde_x)
plot_cell_trajectory(cde_x, color_by = "State") + facet_wrap(~Time_point)
plot_cell_trajectory(cde_x, color_by = "State") + facet_wrap(~Cluster)
plot_cell_trajectory(cde_x, markers = "Cdk1")  + facet_wrap(~Cluster)
plot_cell_trajectory(cde_x, markers = "Aldh1a2", use_color_gradient = TRUE)
plot_complex_cell_trajectory(cde_x, color_by = 'State', show_branch_points = T,
                             cell_size = 0.5, cell_link_size = 0.3, root_states = c(2)) + facet_wrap(~Cluster)

#####
#Multifactorial DEG analyis
#####
#####
#Branch Analysis - BEAM
#####
#Analyse branch
#Input: Ordered according to pseudotime
cde_x <- orderCells(cde_x)

#Perform BEAM for all branches
#List to store BEAMs
l_BEAM_res <- list()
for(i in 1:2){
  BEAM_res_i <- BEAM(cde_x, branch_point = i, cores =12)
  l_BEAM_res[[i]] <- BEAM_res_i
}
names(l_BEAM_res) <- paste("Branch", c(1:2), sep = "_")



saveRDS(l_BEAM_res, paste(PATH_output,"/","BEAM_D0_res_all.rds",sep=""))
l_BEAM_res <- readRDS(paste(PATH_output,"/","BEMA_day0_10_24_56_300_nonPvC.rds",sep=""))
#Eliminate ribosomal genes from branch defining list
rpl.genes <- grep(pattern = "^Rpl", x = l_BEAM_res[[2]]$gene_short_name, value = TRUE)
rps.genes <- grep(pattern = "^Rps", x = l_BEAM_res[[2]]$gene_short_name, value = TRUE)
ribo.genes <- c(rpl.genes, rps.genes)
#get ensemble Ids
BEAM_res_minus_ribo_2 <- subset(l_BEAM_res[[2]], !(gene_short_name %in% ribo.genes))

#Branch 1
plot_genes_branched_heatmap(cde_x[row.names(subset(l_BEAM_res[[1]],
                                                   qval < 1e-10)),],
                            branch_point = 1,
                            num_clusters = 5,
                            #length(levels(as.factor(pData(cde_x)$Cluster_Seurat))),
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = T)
#Branch 2
plot_genes_branched_heatmap(cde_x[row.names(subset(BEAM_res_minus_ribo_2,
                                                   qval < 1e-60)),],
                            branch_point = 2,
                            #num_clusters = length(levels(as.factor(pData(cde_x)$Cluster))),
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T)

#Branch 3
plot_genes_branched_heatmap(cde_x[row.names(subset(l_BEAM_res[[2]],
                                                   qval < 1e-40)),],
                            branch_point = 1,
                            num_clusters = length(levels(as.factor(pData(cde_x)$Cluster))),
                            cores = 16,
                            use_gene_short_name = T,
                            show_rownames = T)

#Branch 3
plot_genes_branched_heatmap(cde_x[row.names(subset(l_BEAM_res[[2]],
                                                   qval < 1e-7)),],
                            branch_point = 4,
                            num_clusters = length(levels(as.factor(pData(cde_x)$Cluster))),
                            cores = 16,
                            use_gene_short_name = T,
                            show_rownames = T)
#####
#Check clusters
#####
expressed_genes <- readRDS(paste(PATH_output,"/expressed_genes.Rds",sep = ""))
diff_test_res <- differentialGeneTest(cde_x[expressed_genes,],
                                      #fullModelFormulaStr = "~Clusters",
                                      cores = 4)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
# Note: time trajectories are usefull to identify key genes
#       but not required
#Set trajectory genes in cde
cde_x <- setOrderingFilter(cde_x, ordering_genes)
plot_ordering_genes(cde_x)

#plot by condition/timepoint
plot_cell_trajectory(cde_x, color_by = "Time_point")
#plot by state
plot_cell_trajectory(cde_x, color_by = "State")
#plot by cluster
plot_cell_trajectory(cde_x, color_by = "State") + facet_wrap(~Cluster)



test_Atf <- subset(l_BEAM_res[[1]], gene_short_name %in% c("Atf1","Atf2","Atf3","Atf4","Atf5",
                                                           "Atf6","Nfkb1","Nkx2-3","Batf","Egr1",
                                                           "Egr2","Bach2",
                                                           "Klf1","Klf2","Klf3","Klf4","Klf5",
                                                           "Stat1","Stat2","Stat3","Stat4",
                                                           "Ctcf","Junb") &
                     qval < 1e-1)

#####
#Ordering based on genes that differ between clusters
#####
# Pick genes that are expressed in >5% of the cell
cde_x <- detectGenes(cde_x, min_expr = 0.1)
fData(cde_x)$use_for_ordering <-
  fData(cde_x)$num_cells_expressed > 0.05 * ncol(cde_x)
#Plot variance explainedby picked genes
plot_pc_variance_explained(cde_x, return_all = F)
#Reduce dimensions using genes
cde_x <- reduceDimension(cde_x,
                         max_components = 2,
                         norm_method = 'log',
                         num_dim = 5,
                         reduction_method = 'tSNE',
                         verbose = T)

#cluster cells
cde_x <- clusterCells(cde_x, num_clusters = 13)
plot_cell_clusters(cde_x, color_by = 'as.factor(Cluster)')

#Decision plot to define Rho and P
plot_rho_delta(cde_x, rho_threshold = 2, delta_threshold = 4 )

#Re-cluster cells implementing cut-off
cde_x <- clusterCells(cde_x,
                      rho_threshold = 70,
                      delta_threshold = 4,
                      skip_rho_sigma = T,
                      verbose = F)

#####
#Distinguish Cluster or State
#####
#Calculation DEGs
# Determine how good a gene explains the state of a cell by checking/building a model
# 1) that knows about the cluster annotation
# 2) that does not know about the cluster annotation
#Compare output of models and choose highest scoring models
#DEG according to Cluster
cde_x_Cluster_DEG <- differentialGeneTest(cde_x,
                                          fullModelFormulaStr = "~Cluster",
                                          cores = 8)
#DEG according to State
cde_x_State_DEG <- differentialGeneTest(cde_x, 
                                        fullModelFormulaStr = "~State",
                                        cores = 8)
#DEG according to pseudotime
cde_x_PseuTim_DEG <- differentialGeneTest(cde_x,
                                          fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                          cores = 8)


# Select genes that are significant at an FDR < 10%
sig_genes_State <- subset(cde_x_State_DEG, qval < 0.01)
sig_genes_Cluster <- subset(cde_x_Cluster_DEG, qval < 0.01)
sig_genes_PseuTim <- subset(cde_x_PseuTim_DEG, qval < 0.01)
saveRDS(sig_genes_State, paste(PATH_output,"/","sig_genes_State.rds",sep=""))
saveRDS(sig_genes_Cluster, paste(PATH_output,"/","sig_genes_Cluster.rds",sep=""))
saveRDS(sig_genes_PseuTim, paste(PATH_output,"/","sig_genes_PseuTim.rds",sep=""))
sig_genes_State <- readRDS(paste(PATH_output,"/","sig_genes_State.rds",sep=""))
sig_genes_Cluster <- readRDS(paste(PATH_output,"/","sig_genes_Cluster.rds",sep=""))
sig_genes_PseuTim <- readRDS(paste(PATH_output,"/","sig_genes_PseuTim.rds",sep=""))

sig_genes_PseuTim_TRUE <- subset(sig_genes_PseuTim, use_for_ordering == "TRUE")
sig_genes_Cluster_TRUE <- subset(sig_genes_Cluster, use_for_ordering == "TRUE")
sig_genes_Cluster_TRUE <- subset(sig_genes_Cluster, use_for_ordering == "TRUE")

#Check intesections
sig_genes_SCP <- Reduce(intersect, list(sig_genes_State$gene_short_name,sig_genes_Cluster$gene_short_name,sig_genes_PseuTim$gene_short_name))
sig_genes_SC <- Reduce(intersect, list(sig_genes_State$gene_short_name,sig_genes_Cluster$gene_short_name))
sig_genes_CP <- Reduce(intersect, list(sig_genes_Cluster$gene_short_name,sig_genes_PseuTim$gene_short_name))
sig_genes_SP <- Reduce(intersect, list(sig_genes_State$gene_short_name,sig_genes_PseuTim$gene_short_name))
print(paste("Number of common sig. SCP:", length(sig_genes_SCP), sep=" "))
print(paste("Number of common sig. SC:", length(sig_genes_SC), sep=" "))
print(paste("Number of common sig. CP:", length(sig_genes_CP), sep=" "))
print(paste("Number of common sig. SP:", length(sig_genes_SP), sep=" "))

# Note: Use genes commonly identified for state and pseudotime
sig_genes_common_SP <- subset(sig_genes_PseuTim, gene_short_name %in% sig_genes_SP)
head(sig_genes_common_SP[,c("gene_short_name", "pval", "qval")])

#Plot genes that describe Pseudotime
#Get Ensemble IDs for genes of interest
sig_genes_SP_ensem <- rownames(subset(sig_genes_PseuTim, gene_short_name %in% sig_genes_SP))

cde_x <- setOrderingFilter(cde_x, sig_genes_SP_ensem)
plot_ordering_genes(cde_x)
plot_pseudotime_heatmap(cde_x[sig_genes_SP_ensem,],
                        num_clusters = length(levels(as.factor(pData(cde_x)$Cluster))),
                        cores = 4,
                        show_rownames = T)
#trend_formula = "~Cluster")


