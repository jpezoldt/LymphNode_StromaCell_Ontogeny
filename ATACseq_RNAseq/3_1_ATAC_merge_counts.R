# Author: Vicnent Gardeux
# Adapted by: Joern Pezoldt
# 12.07.2018
# Function:
# 1) Merges tables of count per peak from homer ATAC-seq pipeline, by id


#Libraries
require(data.table)

# Input required: Set path to directory with the homer .txt files (output of annotatedpeak.pl)
setwd("/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all")
# Input required: Set name of experiment
name = "ATAC_FSC_all"

#Set directory
# Input required: Set path to peak count tables
path <- paste(getwd(),"/homer/Overlap_Group_Merged/Run_4_in_all",sep="")

#Merge tables over id column
merge_counts <- function(path) {
  count.list <- list()
  countfiles <- list.files(path, pattern=".txt$", full.names = T)
  for (i in seq(length(countfiles))){
    count.list[[i]] <- fread(file=countfiles[i], skip = 1, header = F, select = c(1,20), 
                             col.names = c("id", (strsplit(basename(countfiles[i]), ".txt")[[1]])))
    setkey(count.list[[i]],id)
  }
  count <- Reduce(merge, count.list)
  return(count)
}

#Run function
count <- merge_counts(path)
# Note: For ATAC-seq data multiply by 2x as only one strand is counted
count <- count * 2
#Export table
write.table(count, paste(path,"/", name,  ".txt",sep=""), quote = F, sep="\t", row.names = F)




