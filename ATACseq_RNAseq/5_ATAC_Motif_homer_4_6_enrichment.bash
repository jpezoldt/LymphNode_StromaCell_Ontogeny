#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_FSC_all/std_err_out/ATAC_FSC_all_v46_si500_bgTSS_noExt_lowRNA.err
#BSUB -o /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_FSC_all/std_err_out/ATAC_FSC_all_v46_si500_bgTSS_noExt_lowRNA.out
#BSUB -J Motif[1-16]
#BSUB -M 40000000
#BSUB -R rusage[mem=40000]
#BSUB -n 16
#BSUB -u joern.pezoldt@epfl.ch

# Author by: Joern Pezoldt
# 25.07.2018
# Input:
# a) background bed (e.g. regions not DEG but expressed and not DAR but open)
# Function:
# 1) Find Motifs 
# Note: Run on VitalIT

module use /software/module/;
module add UHTS/Analysis/homer/4.6;

parameter=si500_v46_bgTSS_noExt_;
rootdir=/scratch/el/monthly/bdeplanc/pezoldt;
motifdir=$rootdir/Analysis/ATAC_FSC_all/homer/Motif_Input_R/Run_8_low_all/BED_Peaks;
backgrounddir=$rootdir/Analysis/ATAC_FSC_all/homer/Motif_Input_R/Run_8_low_all/background_all_TSS_minus_input_genes_lowRNAseq.bed;
resultdir=$rootdir/Analysis/ATAC_FSC_all/homer/Motif_Input_R/Run_8_low_all/Output_Regions;
genomedir=$rootdir/Data/Genome_mouse/GRCm38.87/Mus_musculus.GRCm38.87.dna.primary_assembly.fa;
preparseddir=$rootdir/tmp;
name=`ls $motifdir | grep '.bed' | sed "${LSB_JOBINDEX}q;d"`;
sample=`ls $motifdir | grep '.bed' | sed "${LSB_JOBINDEX}q;d" | cut -d '.' -f 1`;

echo $sample
echo $name

#mkdir -p $resultdir/${sample};

#Find Motifs
findMotifsGenome.pl $motifdir/${name} $genomedir $resultdir/$parameter${sample} -bg $backgrounddir -size 500 -len 10


