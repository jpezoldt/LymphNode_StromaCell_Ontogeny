#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_FSC_all/ATAC_FSC_all.err
#BSUB -o /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_FSC_all/ATAC_FSC_all.out
#BSUB -J ATAC_pilot[1-12]
#BSUB -M 40000000
#BSUB -R rusage[mem=40000]
#BSUB -n 12
#BSUB -u joern.pezoldt@epfl.ch

# Author: Vincent Gardeaux
# Adapted by: Joern Pezoldt
# 17.07.2018
# Function:
# 1) MAPs paired-end reads using STAR
#		Note: Forward and Reverse read need to pe stored in two different folders
# 2) Generates index files
# 3) Marks and removes duplicates via picard.jar
# 4) Fastqc on aligned and de-duplicated reads
# 5) Removes blacklistes regions
# 6) Performs MACS2
#		Note: 	- Check quality of Peak-calling using the bigwig files in the igv-browser
#				- Choose one set of identifier condition (e.g. broad peaks)
#				- Pipe into "2_b_ATAC_Sort_Thresh_MACS2.R"
#				- Use output of "2_b_ATAC_Sort_Thresh_MACS2.R" (.bed-file with threshed peaks)
# 7) Counts reads per peak using the .bed-track from "2_b_ATAC_Sort_Thresh_MACS2.R"
# Note: Run on VitalIT


module use /software/module/;
module add UHTS/Quality_control/fastqc/0.11.2;
module add Development/java;
module add UHTS/Aligner/STAR/2.5.3a;
module add UHTS/Analysis/samtools/1.3;
module add UHTS/Analysis/HTSeq/0.6.1;
module add UHTS/Analysis/macs/2.1.1.20160309;
module add UHTS/Analysis/deepTools;
module add UHTS/Analysis/BEDTools/2.26.0;
module add UHTS/Analysis/homer/4.9;

rootdir=/scratch/el/monthly/bdeplanc/pezoldt;
fastqdir=$rootdir/Data/ATAC_FSC_all/Trimmed;
resultdir=$rootdir/Analysis/ATAC_FSC_all;
genomedir=$rootdir/Data/Genome_mouse/GRCm38.87;
indexdir=$rootdir/Data/Genome_mouse/GRCm38.87/STARIndex;
blacklist=$genomedir/mm10.blacklist_kundaje_20180718_final.bed;

#get file names
#name=`ls $fastqdir/R1 | grep '.fastq.gz$' | sed "${LSB_JOBINDEX}q;d" | cut -d '_' -f 1-4`;
name=`ls $fastqdir/R1 | grep '.fastq.gz$' | sed "${LSB_JOBINDEX}q;d" | cut -d '_' -f 1-5`;

#sample names
#sample=`echo $name |cut -d '_' -f 1,2`;
sample=`echo $name |cut -d '_' -f 1-2`;

echo $name
echo $sample

#Count peaks
#makeTagDirectory $resultdir/aligned/$sample/tags $resultdir/aligned/${sample}/${sample}.nodup.bam
annotatePeaks.pl $resultdir/peaks/broad/Overlap_Merged/Run_4_in_all/ATAC_FSC_all_broad_merged_peaks.bed mm10 -size given -noadj -organism mouse -d $resultdir/aligned/$sample/tags > $resultdir/homer/Overlap_Group_Merged/Run_4_in_all/${sample}.txt
