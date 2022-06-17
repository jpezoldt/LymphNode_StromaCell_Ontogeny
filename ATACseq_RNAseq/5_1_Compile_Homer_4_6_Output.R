# FILE: Homer_for_Joern.R -----------------------------------------------------
#
# DESCRIPTION : Post-processing of the motif enrichment done with Homer
#
# USAGE: 
#
# OPTIONS:  none
# REQUIREMENTS:  data.table, XML
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
# COMPANY:  EPFL, Lausanne, Switzerland
# VERSION:  1
# CREATED:  21.08.2017
# REVISION: 21.08.2017

library(data.table)
library(XML)

# FUNCTIONS -------------------------------------------------------------------
#' readHomerKnownIn
#' Reads homer results in
#' @param filePath path to the file
#' @param qvalueCut cut off on q value
#' @param toRemove what to remove from file name before parcing
#' @return data frame with results from homer, q-value < 0.05, 
#' motif name parsed
readHomerKnownIn <- function(filePath, qvalueCut = 0.05) {
  homer <- fread(filePath)
  toParse <- strsplit(filePath, '/')[[1]]
  
  # group = folder name
  homer[, group := toParse[length(toParse) - 1]]
  
  # motif len, doesn't work for v2
  homer[, motifLen := NA]
  
  # select needed columns
  homer <- homer[, c(10:11, 1, 2, 3, 5, 7, 9), with = F]
  colnames(homer) <- c('group', 'motifLen', 'MotifName', 'Consensus', 
                       'P-value', 'q-value', 'PercTarg', 'PercBg')
  homer$PercTarg <- as.numeric(gsub('%', '', homer$PercTarg))
  homer$PercBg <- as.numeric(gsub('%', '', homer$PercBg))
  homer$MotifName <- gsub('/.*', '', homer$MotifName)
  homer$MotifName <- gsub('\\(.*', '', homer$MotifName)
  # add info about motif file location
  homer$File <- paste0(gsub('.txt', '/known', filePath), 1:nrow(homer), 
                       '.motif')
  homer$File <- sapply(homer$File, function(x) if(file.exists(x)) {x} 
                       else {NA})
  
  # select motifs passing q-value cut-off
  homer <- homer[`q-value` < qvalueCut & `P-value` < qvalueCut]
  if (nrow(homer) != 0) {
    homer <- homer[order(`q-value`), ]
  }
  
  # sometimes, then the same motifs come from different sourses, we see 
  # duplication of the entries
  homer <- homer[!duplicated(homer[, -9, with = F])]
  # remove all seq biases
  homer <- homer[!grepl('SeqBias', MotifName)]
  # add rank
  homer[, Rank := 1:nrow(homer)]
  # add fold change
  homer[, FC := PercTarg / PercBg]
  # rearrange
  homer[, .(group, motifLen, MotifName, Consensus, `P-value`, `q-value`, 
            PercTarg, PercBg, Rank, FC, File)]
}

#' readHomerKnownAllIn_v2
#' Reads in all knownMotifs.txt from all folders and subfolders of homerResDir
#' Version 2, not compatible with version 1
#' @param homerResDir directory full of homer results
#' @param qValCut cutoff on q-value
#' @return data table with homer results
readHomerKnownAllIn_v2 <- function(homerResDir, qValCut) {
  # read all the results files from all the folders
  resultDirs <- list.dirs(homerResDir, recursive = F)
  knownResultsFiles <- sapply(resultDirs, 
                              function(x) list.files(x, 'knownResults.txt',
                                                     recursive = T, 
                                                     full.names = T))
  knownResultsFiles <- unlist(knownResultsFiles)
  knownResultsFiles <- knownResultsFiles[lapply(knownResultsFiles, 
                                                length) != 0]
  knownResults <- lapply(knownResultsFiles, readHomerKnownIn, qValCut)
  
  # merge
  knownResults <- knownResults[sapply(knownResults, 
                                      function(x) is.data.frame(x))]
  knownResults <- do.call(rbind, knownResults)
}

#' readHomerDeNovoIn
#' Reads-in one .motif file (result of de-novo motif discovery from Homer)
#' @param filePath path to the file
#' @return data table with columns motifLen, MotifName, Consensus, P-value,
#' q-value, PercTarg, PercBg, Rank
readHomerDeNovoIn <- function(filePath) {
  # read-in and leave only 1st line starting with ">" because it contains
  # info about motif
  motifInfo <- read.table(filePath, stringsAsFactors = F, nrows = 1)
  motifInfo <- motifInfo[c(1, 2, length(motifInfo))]
  motifInfo <- c(motifInfo, gsub('', '', gsub('%.*', '', motifInfo[3])))
  Rank <- gsub('.motif', '', gsub('.*/motif', '', filePath))
  Rank <- as.integer(Rank)
  group <- strsplit(filePath, '/')[[1]]
  group <- group[length(group) - 2]
  motifInfoDT <- data.table(group = group,
                            motifLen = nchar(gsub('>', '', motifInfo[1])),
                            MotifName = gsub('>', '', motifInfo[1]), 
                            Consensus = gsub('>', '', motifInfo[1]),
                            Pvalue = as.numeric(gsub('.*:', '', 
                                                     motifInfo[3])),
                            qvalue = NA,
                            PercTarg = as.numeric(gsub('%.*', '', 
                                                       gsub('.*\\(', '', 
                                                            motifInfo[4]))),
                            PercBg = as.numeric(gsub('%.*', '', 
                                                     gsub('.*\\(', '', 
                                                          motifInfo[3]))),
                            Rank = Rank,
                            BestMatch = gsub('\\(.*', '', gsub('.*:', '',
                                                               motifInfo[2])))
  motifInfoDT[, MotifName := BestMatch]
  motifInfoDT[, MotifName := gsub('/.*', '', MotifName)]
  motifInfoDT[, FC := PercTarg / PercBg]
  colnames(motifInfoDT)[5:6] <- c('P-value', 'q-value')
  motifInfoDT$filePath <- filePath
  
  motifInfoDT
}

#' parseHomerDeNovoHTML_v2
#' Parses HTML output of de-novo motifs from homer. Not compatible with v1.
#' @param filePath path to html 
#' @return data frame of if the de novo motif is false positive
parseHomerDeNovoHTML_v2 <- function(filePath) {
  # parse html 
  html <- xmlParse(filePath)
  # get to the table
  html <- xmlChildren(html)$HTML[[2]][['TABLE']]
  isFP <- c() # vector which will contain if it's false positive
  motifPath <- c() # path to motif file, need for merging with rest of data
  i = 2
  while(!is.null(html[[i]])) {
    # get motif paths
    oneMotifPath <- xmlAttrs(html[[i]][2]$TD[['IMG']])
    oneMotifPath <- gsub('homerResults.html', oneMotifPath, filePath)
    oneMotifPath <- gsub('logo.png', 'motif', oneMotifPath)
    motifPath <- c(motifPath, oneMotifPath)
    
    # getting to FP
    if (is.null(html[[i]][1]$TD[['FONT']])) {
      isFP <- c(isFP, F)
    } else {
      if (xmlAttrs(html[[i]][1]$TD[['FONT']])[['color']] == 'red') {
        isFP <- c(isFP, T)
      } else {
        isFP <- c(isFP, NA)
      }
    }
    i <- i + 1
  }
  result <- data.table(filePath = motifPath, FP = isFP)
  result
}

#' readHomerDeNovoAllIn_v2
#' Reads in all .motif from all folders and subfolders of homerResDir. Not 
#' compatible with v1
#' @param homerResDir directory full of homer results
#' @param enrichCode code for the enrichment performed with homer
#' @return data table with homer results
readHomerDeNovoAllIn_v2 <- function(homerResDir) {
  # read all the results files from all the folders
  resultDirs <- list.dirs(homerResDir, recursive = T, full.names = T)
  resultDirs <- resultDirs[grepl('homerResults', resultDirs)]
  deNovoResultsFiles <- sapply(resultDirs, 
                               function(x) list.files(x, 'f\\d+.motif$',
                                                      recursive = T, 
                                                      full.names = T))
  deNovoResultsFiles <- unlist(deNovoResultsFiles) 
  deNovoResultsFiles <- deNovoResultsFiles[lapply(deNovoResultsFiles, 
                                                  length) != 0]
  deNovoResults <- lapply(deNovoResultsFiles, readHomerDeNovoIn)
  
  # merge
  deNovoResults <- deNovoResults[sapply(deNovoResults, 
                                        function(x) is.data.frame(x))]
  deNovoResults <- do.call(rbind, deNovoResults)
  
  # if Best match is empty, then it's seq bias
  deNovoResults <- deNovoResults[BestMatch != '']
  
  # add info about if motif is false-positive
  allHtml <- sapply(homerResDir, 
                    function(x) list.files(x, 'homerResults.html',
                                           recursive = T, 
                                           full.names = T))
  deNovoHtml <- unlist(allHtml)
  names(deNovoHtml) <- deNovoHtml
  deNovoFP <- lapply(deNovoHtml, parseHomerDeNovoHTML_v2)
  deNovoFP <- do.call(rbind, deNovoFP)
  
  # merge
  setkey(deNovoFP, filePath)
  setkey(deNovoResults, filePath)
  deNovoResults <- merge(deNovoResults, deNovoFP)
  setnames(deNovoResults, "filePath", 'File')
  
  deNovoResults 
}

#' readHomerMotif
#' Reads-in result of homer motif map 
#' @param filePath path to the file
#' @return data table with columns FASTA, ID, Offset, Sequence, Motif, Name, 
#' Strand, MotifScore, chr, regStart, regEnd, start, end
#' where start and end is an absolute start and end of motif
readHomerMotif <- function(filePath) {
  motifDF <- fread(filePath, sep = '\t', header = T)
  motifDF[, chr := sapply(motifDF$`FASTA ID`, 
                          function(x) strsplit(x, ':')[[1]][1])]
  motifDF[, regStart := sapply(`FASTA ID`,
                               function(x) as.integer(strsplit(gsub('.*:', '',
                                                                    x), '-')[[1]][1]))]
  motifDF[, regEnd := sapply(`FASTA ID`,
                             function(x) as.integer(strsplit(gsub('.*:', '',
                                                                  x), '-')[[1]][2]))]
  motifDF[, start := ifelse(Strand == '+', 
                            regStart + ((regEnd - regStart) / 2) + Offset + 1,
                            regStart + ((regEnd - regStart) / 2) + Offset + 2 - 
                              nchar(Sequence))]
  motifDF[, end := start + nchar(Sequence) - 1]
  motifDF
}

#' motifToMeme
#' Converts a Homer .motif PWM/PFM into MEME format
#' @param inFile path to input file
#' @return list ready to be printed in MEME format
motifToMeme <- function(filePath) {
  # Reading the input file
  motif <- scan(file = filePath, character(0), sep = "\n", quote = NULL)
  motifName <- motif[1]
  
  # parcing of the matrix
  motifMatr <- motif[2:length(motif)] 
  motifMatr <- strsplit(motifMatr, split = "\t")
  motifMatr <- do.call(rbind, lapply(motifMatr, as.numeric))
  
  # geting name and p-value
  pval <- gsub('.*:', '', motifName)
  concensus <- gsub('>', '', strsplit(motifName, split = "\t")[[1]][1])
  motifName <- strsplit(motifName, split = "\t")[[1]][2]
  if (!grepl('BestGuess', motifName)) {
    motifName <- gsub('/.*', '', motifName)
    motifName <- gsub('\\(.*', '', motifName)
  } else {
    motifName <- gsub('\\(.*', '', gsub('.*:', '', motifName))
    motifName <- gsub('/.*', '', motifName)
  }
  
  result <- list(header = paste0('MOTIF ', concensus, ' ', motifName, '\n',
                                 'letter-probability matrix: alength= 4 w=',
                                 nrow(motifMatr), ' nsites= 20 E= ', pval),
                 matrix = motifMatr)
  result
}

#' readInTomTom
#' Reads in and parces tomtom output
#' @param filePath path to tomtom output file
#' @param eValCut cut off on e value
#' @return data table with MotifName, Consensus, DrosoTF, E-value, 
#' TomTom_consensus TomTom_orientation
readInTomTom <- function(filePath, eValCut) {
  tomtom <- fread(filePath, header = T)
  tomtom <- tomtom[`E-value` < eValCut]
  colnames(tomtom)[1] <- 'Consensus'
  tomtom <- tomtom[order(Consensus, `E-value`)]
  tomtom <- tomtom[,.SD[1], Consensus]
  tomtom[, `Target ID` := gsub('_.*', '', `Target ID`)]
  
  tomtom <- tomtom[,.(Consensus, `Target ID`, `E-value`, `Target consensus`,
                      Orientation)]
  colnames(tomtom)[4:5] <- c('TomTom_consensus', 'TomTom_orientation')
  tomtom
}

#####
# INPUTS ----------------------------------------------------------------------
#####
# folder with all homer results
sample_ID <- "SPF"
homerInDir <- paste("/home/pezoldt/NAS2/pezoldt/Analysis/ATACseq/ATAC_FSC_all/Motif/Homer_Output/allHomerResults/size_500bp/", sample_ID, sep = "")
setwd(homerInDir)
# cutoff on q-value for the motif enrichment of known motifs 
qValueCutoff <- 0.1

# HOMER MOTIF DISCOVERY read-in -----------------------------------------------
# KNOWN MOTIFS
known <- readHomerKnownAllIn_v2(homerInDir, qValueCutoff)
# put here length which you used. Unfortunately, I don't fish it out of the 
# files. But you don't have a lot of samples, so it's not a big deal
known$motifLen <- as.numeric(known$motifLen)
known$motifLen <- 10
# save RDS, because sometimes it takes long time to read-in
saveRDS(known, 'known_SPF.Rds')
write.table(known, paste(sample_ID, "known_homer_compile.txt") ,dec=".", sep="\t")
#known_10 <- readRDS('known.Rds')

# DE-NOVO MOTIFS
deNovo <- readHomerDeNovoAllIn_v2(homerInDir) 
deNovo$motifLen <- as.numeric(deNovo$motifLen)
deNovo$motifLen <- 10
# save RDS, because it takes long time to read-in
saveRDS(deNovo, 'deNovo.Rds')
write.table(deNovo, paste(sample_ID, "deNovo_homer_compile.txt") ,dec=".", sep="\t")
#deNovo <- readRDS('deNovo.Rds')

# COMPARISON OF MOTIFS WITH TOMTOM --------------------------------------------
# 1) Install MEME suit: http://meme-suite.org/doc/install.html?man_type=web
# 2) Download DBs of all known motifs: 
# meme-suite.org/meme-software/Databases/motifs/motif_databases.12.17.tgz
# Select the ones for mouse and put then in the folder tomtom_db_mouse

# comparison of motifs for the motifs enriched as "known"
# known - PWM should be the same for the same MotifName
printInMemeFormat(known[!duplicated(MotifName)]$File,  "known.meme")
system(paste('tomtom known.meme tomtom_db_mouse/*.meme -oc', 
             'known_tomtom'))
system('mv known_tomtom/tomtom.txt known_tomtom.txt')
system('rm -r known_tomtom')
# read in results of TomTom, cut off on p-value 0.05
known_tomtom <- readInTomTom('known_tomtom.txt', 0.05)
# merge with homer results
setkey(known_tomtom, Consensus)
setkey(known, Consensus)
known <- merge(known, known_tomtom, all = T)
known_10[, MouseTF := sapply(1:nrow(known),
                             function(x) if(is.na(known$`Target ID`[x])) {
                               paste0(known$MotifName[x], '(not in Mouse DB)')
                             } else {known$`Target ID`[x]})]
saveRDS(known_10, 'known_withTomTomTF.Rds')

# repeat the code above for de-novo motifs