#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_FSC_all/ATAC_FSC_all_enrich.err
#BSUB -o /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_FSC_all/ATAC_FSC_all_enrich.out
#BSUB -J ATAC_FSC_all
#BSUB -M 40000000
#BSUB -R rusage[mem=40000]
#BSUB -n 1
#BSUB -u joern.pezoldt@epfl.ch


module use /software/module/;
module add UHTS/Quality_control/fastqc/0.11.2;
module add Development/java/1.8.0_121;
module add UHTS/Aligner/STAR/2.5.3a;
module add UHTS/Analysis/samtools/1.3;
module add UHTS/Analysis/HTSeq/0.6.1;
module add UHTS/Analysis/BEDTools/2.26.0;
module add UHTS/Analysis/deepTools;

rootdir=/scratch/el/monthly/bdeplanc/pezoldt;
resultdir=$rootdir/Analysis/ATAC_FSC_all/peaks;
genomedir=$rootdir/Genome_mouse/GRCm38.87;

method='broad'
experiment='ATAC_FSC_all'
methodwet='ATAC'


#provide names
names=`ls /scratch/el/monthly/bdeplanc/pezoldt/Analysis/$experiment/aligned/*/*.nodup.bam | cut -d '.' -f1 | cut -d '/' -f10`
echo $names

#Provide MACS2 peak annotation files (.bed) per sample
#Provide MACS2 peaks present in all replicates
cat $resultdir/$method/Overlap_Merged/Run_4_min2tracks/*.bed | sort -k 1,1 -k2,2n > $resultdir/$method/Overlap_Merged/Run_4_min2tracks/ATAC_FSC_all_${method}_concat_peaks.bed

#merge all genomic regions if overlap at lest 1bp
bedtools merge -i $resultdir/$method/Overlap_Merged/Run_4_min2tracks/ATAC_FSC_all_${method}_concat_peaks.bed > $resultdir/$method/Overlap_Merged/Run_4_min2tracks/ATAC_FSC_all_${method}_merged_peaks.bed


