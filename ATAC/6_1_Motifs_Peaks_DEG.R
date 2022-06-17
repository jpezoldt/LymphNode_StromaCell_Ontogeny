# Annotate motifs to cohorts of DARs and DEGs and explore expression

# Author: Joern Pezoldt
# Date: 26.10.2018

# Functions:
# 1) Compile known motifs from homer output into one motif list
# 2) Use "homer" to annotate motifs to peaks
# 3) Obtain Peaks that contain motif
# 4) Check to which gene Peaks belong 
# 5) Classify TFs into Core and Unique sets per condition
# 6) Plot expression of genes for Core and unique TF sets


#Library
library(biomaRt)
library(pheatmap)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GenomicFeatures)
library(ChIPpeakAnno)

#####
#PATHs and global variables
#####
# Note: Input required
PATH_general <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/SPF/"
PATH_motifs <- "/knownResults/"
homer_find_conditions <- "si500_v46_bgTSS_noExt_"
conditions <- c("Closed_DEG_SPF_lowRNAseq",
                "No_DAR_DOWN_SPF_lowRNAseq",
                "No_DAR_UP_SPF_lowRNAseq",
                "NoDEG_Closed_SPF_lowRNAseq",
                "NoDEG_Open_SPF_lowRNAseq",
                "Open_DEG_SPF_lowRNAseq")
PATH_TF_motifs_in_gene <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Regions_by_Motif"

#Databases
#All features
mm10_TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#TSS sites BioMart
mm10 = useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl") 
TSS.mouse.mm10 = getAnnotation(mart=mm10, featureType="TSS")
Exon.mouse.mm10 = getAnnotation(mart=mm10, featureType="Exon")

#Global variables
# score for motifs obtained for mapping motifs to peaks
# Note: distribution of hist(t_TF_binding_sites_i$MotifScore) yields 5 as a robust cutoff
motif_score = 5

#####
#Functions
#####


#Concatenate motifs -----------------------------------
sampleNames_complete <- paste(homer_find_conditions, conditions, sep = "")
PATH_sampleNames <- paste(PATH_general,sampleNames_complete,PATH_motifs,sep = "")
#' concatenateKnownMotifs
#' loads single motif files from homer (known motifs) and compiles them into one txt file
#' @param sampleNames_complete cahractter vector with path to motif files
#' @param PATH_general string to store all files
#' @param conditions character vector of sample names
#' @return concatenated motif tables within "knownMotifs" and higher folder



concatenateKnownMotifs <- function(PATH_sampleNames,conditions,PATH_general){
  for(i in 1:length(PATH_sampleNames)){
    PATH_sampleNames_i <- PATH_sampleNames[i]
    setwd(PATH_sampleNames_i)
    names_i <- list.files(pattern="*.motif")
    names_i <- paste(PATH_sampleNames_i,names_i,sep = "")
    myfiles_i <- lapply(names_i, read.delim, header = FALSE)
    motif_file_i <- do.call("rbind", myfiles_i)
    motif_file_i[is.na(motif_file_i)] <- ""
    colnames(motif_file_i) <- NULL
    write.table(motif_file_i,paste(PATH_sampleNames_i,conditions[i],"_knownMotif_compiled.txt", sep = ""), sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)
    write.table(motif_file_i,paste(PATH_general,conditions[i],"_knownMotif_compiled.txt", sep = ""), sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)
  }
}

#Run Function
concatenateKnownMotifs(PATH_sampleNames,conditions,PATH_general)


#Identify genes targeted by TFs
# Explanation Input:
# Input 1: TXT with the peak identifiers, the start and the stop regions (and TF bindings per peak)
# Input 2: TXT with altered peak identifier (e.g. preceding BG), offset from peak center and Motif name

#' pinpointTFMotifstoGenes
#' generates RDS:
#' 1) containing a table that counts the incidence number (motif occurences) for each TF within each gene 
#' @param PATH_sampleNames cahractter vector with path to motif files
#' @param TSS.mouse.mm10 bimoart annotation for TSS
#' @return l_TF_binding_freq returns table with incidence rates of TF binding per gene
PATH_samples_TF_Motif_Gene <- paste(PATH_TF_motifs_in_gene,"/",conditions,sep = "")



pinpointTFMotifstoGenes <- function(PATH_samples_TF_Motif_Gene, conditions){
  l_TF_binding_freq <- list()
  for(i in 1:length(PATH_samples_TF_Motif_Gene)){
    PATH_samples_TF_Motif_Gene_i <- PATH_samples_TF_Motif_Gene[i]
    conditions_i <- conditions[i]
    print(i)
    #Load TF binding sits
    #t_TF_binding_sites_i <- read.delim(paste(PATH_samples_TF_Motif_Gene_i,"/",conditions_i,"_MotifInstances.txt",sep = ""),stringsAsFactors=FALSE)
    #apply cutoff
    #motif_score
    #t_TF_binding_sites_i <- subset(t_TF_binding_sites_i, MotifScore >= 5)
    
    #Load Peak location
    t_Peak_location_i <- read.delim(paste(PATH_samples_TF_Motif_Gene_i,"/",conditions_i,"_MotifInstances_Location.txt",sep = ""),stringsAsFactors=FALSE,
                                  header = FALSE)
    #Change Column for TFs
    colnames_i <- as.character(t(as.vector(t_Peak_location_i[1,])))
    colnames_i <- unlist(lapply(strsplit(colnames_i, "/"), function(x){
      x[[1]]
    }))
    #Eliminate SeqBias: columns
    colnames_i <- unlist(lapply(strsplit(colnames_i, "Bias:"), function(x){
      x[[1]]
    }))
    
    colnames(t_Peak_location_i) <- colnames_i
    
    t_Peak_location_i <- t_Peak_location_i[2:nrow(t_Peak_location_i),]
    #t_Peak_location <- t_Peak_location[,1:4]
    colnames(t_Peak_location_i)[1] <- "ID"
    #Annotate Peak to gene
    gr_Peaks_i <- toGRanges(t_Peak_location_i, names = t_Peak_location[,1])
    #Annotate Peaks
    gr_Peaks_TSS_i <- annotatePeakInBatch(gr_Peaks_i, AnnotationData=TSS.mouse.mm10)
    #add gene name
    gr_Peaks_TSS_GeneName_i <- addGeneIDs(annotatedPeak=gr_Peaks_TSS_i, 
                                           feature_id_type="ensembl_gene_id",
                                           orgAnn="org.Mm.eg.db", 
                                           IDs2Add="symbol")
    #Introduces ten new columns
    t_Peaks_TSS_GeneName_i <- as.data.frame(gr_Peaks_TSS_GeneName_i)
    #Rearrange table for counting the TF binding events per TF
    t_Peaks_TSS_GeneName_i <- cbind(t_Peaks_TSS_GeneName_i[,1:22],
                  t_Peaks_TSS_GeneName_i[,(ncol(t_Peaks_TSS_GeneName_i)-11):ncol(t_Peaks_TSS_GeneName_i)],
                  t_Peaks_TSS_GeneName_i[,23:(ncol(t_Peaks_TSS_GeneName_i)-10)])
    #Determine number of TF binding sites per peak/Gene per TF
    t_Peaks_TSS_GeneName_TFs_i <- t_Peaks_TSS_GeneName_i[,35:ncol(t_Peaks_TSS_GeneName_i)]
    #string split for each TF for each Peak
    numbers_i <- matrix(ncol = ncol(t_Peaks_TSS_GeneName_TFs_i), nrow = 0)
    for(j in 1:nrow(t_Peaks_TSS_GeneName_TFs_i)){
      row_j <- as.character(t_Peaks_TSS_GeneName_TFs_i[j,] )
      n_j <- unlist(lapply(strsplit(row_j, '[(]'),function(x){
        length(x) - 1
      }))
      n_j[n_j == -1] <- 0
      numbers_i <- rbind(numbers_i, n_j)
    }
    numbers_i <- as.data.frame(numbers_i)
    t_Peaks_TSS_GeneName_i[,35:ncol(t_Peaks_TSS_GeneName_i)] <- numbers_i
    print(nrow(t_Peaks_TSS_GeneName_i))
    l_TF_binding_freq[[i]] <- t_Peaks_TSS_GeneName_i
  }
  #returned
  l_TF_binding_freq
}

#Run function
l_TF_binding_freq <- pinpointTFMotifstoGenes(PATH_samples_TF_Motif_Gene, conditions)



#####
#Exploratory
#####
#Generate Heatmap of log2 number of TF binding sites per gene
test_TF <- l_TF_binding_freq[[5]]

data_heatmap <- test_TF
drops <- c("Seq","GAGA.repeat")
data_heatmap <- data_heatmap[ , !(names(data_heatmap) %in% drops)]
rownames(data_heatmap) <- paste(data_heatmap$symbol, 1:length(data_heatmap$symbol))
data_heatmap <- data_heatmap[,35:ncol(data_heatmap)]
data_heatmap <- as.matrix(data_heatmap)
data_heatmap <- log2(data_heatmap)
data_heatmap[data_heatmap == -Inf] <- 0
pheatmap(data_heatmap, scale = "none", color = colorRampPalette(c("white", "lightgoldenrod1","chartreuse3","chartreuse4"), space="rgb")(128))


