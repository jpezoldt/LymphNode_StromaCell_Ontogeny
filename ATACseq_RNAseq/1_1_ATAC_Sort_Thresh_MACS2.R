# Author: 
# 12.11.2018
# Function:
# 1) Filter MACS Output
# 2) Filter Peaks not associated to blacklisted regions
# 2) Generate BED for Peaks present in any biological sample

#Libraries
library("GenomicRanges")
library("data.table")
library("stringr")

#####
#Global variables
#####

# Input required: Set path to directory with the MACS2 XLS files
path_input <- setwd("/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_denovo_Treg/peaks/broad/MACS2_XLS")
path_input_blacklist <- "/home/pezoldt/NAS2/pezoldt/Data/ATACseq/mm10.blacklist_kundaje_20180718_final.bed"
path_output <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_denovo_Treg/peaks/broad/MACS2_broad/Run_1_in_all"
path_output_at_least_2 <- "/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_denovo_Treg/peaks/broad/MACS2_broad/Run_1_min2tracks"
# Input required: add group names that allow to distinguish experimental conditions across replicates
experimental_groups <- c("ATAC_TconGF","ATAC_TconSPF","ATAC_Tnaive","ATAC_TregGF","ATAC_TregSPF")

#####
#1) Filter MACS Output
#####
#Per file
#-log10(pValue) <= 1.3
#FC over background 2

#Get file names
MACS2_thresh_peaks_files <- list.files(path_input, pattern=".xls")

#Set storage vectors and lists
count_before <- c()
count_after <- c()
MACS2_thresh_peaks.list <- list()

for(i in seq(length(MACS2_thresh_peaks_files))){
  #Load table
  table_before <- fread(file=MACS2_thresh_peaks_files[i], skip = 24, header = T)
  count_before[i] <- nrow(table_before)
  #Threstable
  table_after <- subset(table_before,  table_before$`-log10(qvalue)` >= 1.3 &  fold_enrichment >= 2)
  count_after[i] <- nrow(table_after)
  
  #Store table in list
  MACS2_thresh_peaks.list[[i]] <- table_after
}
names(MACS2_thresh_peaks.list) <- paste("Thresh_",list.files(pattern="*.xls"),sep="")
#plots number of peaks before vs. after thresholding
plot(count_before, count_after)

#####
#2) Eliminate blacklisted regions
#####
#make Grange object from Blacklist-bed
blacklist <- read.delim(path_input_blacklist, header = FALSE)
gr_blacklist <- GRanges(seqnames = blacklist$V1, IRanges(start = blacklist$V2, end = blacklist$V3))

#Generate empty list to store blacklist depleted peaks
MACS2_thresh_peaks_minus_bl.list <- list()
#loop over each replicate
for(i in seq(length(MACS2_thresh_peaks.list))){
  #fore every replicate
  #take object
  replicate_i <- MACS2_thresh_peaks.list[[i]]
  #generate Granges
  grange_obj_i <- makeGRangesFromDataFrame(replicate_i,
                                           seqnames.field=c("seqnames", "seqname",
                                                            "chromosome", "chrom",
                                                            "chr", "chromosome_name",
                                                            "seqid"),
                                           start.field="start",
                                           end.field=c("end"),
                                           strand.field="strand",
                                           starts.in.df.are.0based=FALSE)
  #output replicate indes
  print(i)
  #and number of peaks
  print(length(grange_obj_i))
  #find overlap
  overlap_index <- suppressWarnings(findOverlaps(grange_obj_i, gr_blacklist, minoverlap = 1))
  #get index numbers for overlapping peaks
  query_overlap <- queryHits(overlap_index)
  #get Granges object for overlapping peaks
  grange_obj_i <- grange_obj_i[-query_overlap]
  #output number of peaks after depletion
  print(length(grange_obj_i))
  MACS2_thresh_peaks_minus_bl.list[[i]] <- grange_obj_i
}

#add names to blacklist depleted replicates
names(MACS2_thresh_peaks_minus_bl.list) <- paste("Thresh_bl_",list.files(pattern="*.xls"),sep="")

#####
#3) Identify Peaks present in all files in all tracks
#####
#for every experimental group
#grab replicates
#make grange objects
#find inner overlap looping over replicates
#compile in one file

#Initiate empty list to store grange object per group
grange_common_group <- list()
for(i in seq(length(experimental_groups))){
  
  #for each experimental group
  group_i <- experimental_groups[i]
  #obtain the threshed MACS2_peak list names
  group_names_i <- names(MACS2_thresh_peaks_minus_bl.list)[str_detect(names(MACS2_thresh_peaks_minus_bl.list),group_i)]
  #grab the respective table
  group_table_i <- MACS2_thresh_peaks_minus_bl.list[group_names_i]
  
  #Intiate empty list to store grange objects per replicate within experimental group
  grange_obj_list_j <- list()
  
  
  for(j in seq(length(group_table_i))){
    #make grange object for each replicate
    grange_obj_j <- group_table_i[[j]]
    
    #store grange in experimental group list
    grange_obj_list_j[[j]] <- grange_obj_j
    
  }
  #name replicates
  names(grange_obj_list_j) <- paste("Grange_",group_names_i,sep="")
  #check if replicates correspond to input
  print(length(grange_obj_list_j))
  print(names(grange_obj_list_j))
  
  #Find common peak regions across replicates per experimental group
  #grange_common_group[[i]] <- Reduce(subsetByOverlaps, grange_obj_list_j)
  
  #loop over Granges List object and iteratively identify the overlap
  #Generate Granges object to start loop
  gr_final <- with(grange_obj_list_j[[1]], GRanges(seqnames = seqnames, IRanges(start = start, end = end)))
  print(length(gr_final))
  for(k in seq(length(grange_obj_list_j)-1)){
    #add plus one to start with the second object
    index = k + 1
    print(index)
    #make Granges object of 2nd in list
    gr_index <- with(grange_obj_list_j[[index]], GRanges(seqnames = seqnames, IRanges(start = start, end = end)))
    #find overlap
    overlap_index <- suppressWarnings(findOverlaps(gr_final, gr_index, minoverlap = 1))
    #get index numbers for overlapping peaks
    query_overlap <- queryHits(overlap_index)
    #get Granges object for overlapping peaks
    gr_final <- gr_final[query_overlap]
    print(length(gr_final))
  }
  grange_common_group[[i]] <- gr_final
}

#Make .bed from grange object
for(i in seq(length(experimental_groups))){
  print(experimental_groups[i])
  print(i)
  merged_table_i <- as.data.frame(grange_common_group[[i]])
  merged_table_i <- cbind(merged_table_i[,1:3],rep(paste("Over_Merged_",experimental_groups[i],sep=""),nrow(merged_table_i)))
  print(head(merged_table_i))
  names(merged_table_i) <- NULL
  #rownames(merged_table_i) <- NULL
  write.table(merged_table_i, paste(path_output,"/", paste("Over_Merged_",experimental_groups[i],".bed",sep=""),sep=""), quote = F, sep="\t", row.names = F)
}

#Take output and pipe it into "20_..."

#####
#3) Identify Overlapping and non-overlapping Peaks
#####
#for every experimental group
#grab replicates
#make grange objects
#find inner overlap looping over replicates
#compile in one file

#Initiate empty list to store grange object per group
grange_common_group <- list()
for(i in seq(length(experimental_groups))){
  
  #for each experimental group
  group_i <- experimental_groups[i]
  #obtain the threshed MACS2_peak list names
  group_names_i <- names(MACS2_thresh_peaks_minus_bl.list)[str_detect(names(MACS2_thresh_peaks_minus_bl.list),group_i)]
  #grab the respective table
  group_table_i <- MACS2_thresh_peaks_minus_bl.list[group_names_i]
  
  #Intiate empty list to store grange objects per replicate within experimental group
  grange_obj_list_j <- list()
  
  
  for(j in seq(length(group_table_i))){
    #make grange object for each replicate
    grange_obj_j <- group_table_i[[j]]
    
    #store grange in experimental group list
    grange_obj_list_j[[j]] <- grange_obj_j
    
  }
  #name replicates
  names(grange_obj_list_j) <- paste("Grange_",group_names_i,sep="")
  #check if replicates correspond to input
  print(length(grange_obj_list_j))
  print(names(grange_obj_list_j))
  
  #Find common peak regions across replicates per experimental group
  #grange_common_group[[i]] <- Reduce(subsetByOverlaps, grange_obj_list_j)
  
  #loop over Granges List object and iteratively identify the overlap
  #Generate Granges object to start loop
  gr_final <- with(grange_obj_list_j[[1]], GRanges(seqnames = seqnames, IRanges(start = start, end = end)))
  print(length(gr_final))
  for(k in seq(length(grange_obj_list_j)-1)){
    #add plus one to start with the second object
    index = k + 1
    print(index)
    #make Granges object of 2nd in list
    gr_index <- with(grange_obj_list_j[[index]], GRanges(seqnames = seqnames, IRanges(start = start, end = end)))
    #find overlap
    overlap_index <- suppressWarnings(findOverlaps(gr_final, gr_index, minoverlap = 1))
    #get index numbers for overlapping peaks
    query_overlap <- queryHits(overlap_index)
    #get Granges object for overlapping peaks
    gr_final <- gr_final[query_overlap]
    #non-overlapping peaks
    gr_nonoverlap <- suppressWarnings(setdiff(gr_final, gr_index))
    #appending overlapping and nonoverlapping regions
    df_final <- as.data.frame(gr_final)
    df_nonoverlap <- as.data.frame(gr_nonoverlap)
    df_final <- rbind(df_final, df_nonoverlap)
    #Make granges object
    gr_final <- makeGRangesFromDataFrame(df_final,
                                             seqnames.field=c("seqnames", "seqname",
                                                              "chromosome", "chrom",
                                                              "chr", "chromosome_name",
                                                              "seqid"),
                                             start.field="start",
                                             end.field=c("end"),
                                             strand.field="strand",
                                             starts.in.df.are.0based=FALSE)
    print(length(gr_final))
  }
  grange_common_group[[i]] <- gr_final
}

#Make .bed from grange object
for(i in seq(length(experimental_groups))){
  print(experimental_groups[i])
  print(i)
  merged_table_i <- as.data.frame(grange_common_group[[i]])
  merged_table_i <- cbind(merged_table_i[,1:3],rep(paste("Over_Merged_",experimental_groups[i],sep=""),nrow(merged_table_i)))
  print(head(merged_table_i))
  names(merged_table_i) <- NULL
  #rownames(merged_table_i) <- NULL
  write.table(merged_table_i, paste(path_output_at_least_2,"/", paste("Over_Merged_",experimental_groups[i],".bed",sep=""),sep=""), quote = F, sep="\t", row.names = F)
}
