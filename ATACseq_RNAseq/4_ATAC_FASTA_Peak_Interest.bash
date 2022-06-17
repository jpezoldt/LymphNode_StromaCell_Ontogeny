#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_FSC_all/std_err_out/ATAC_FSC_make_Fasta.err
#BSUB -o /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_FSC_all/std_err_out/ATAC_FSC_make_Fasta.out
#BSUB -J Motif[1-12]
#BSUB -M 40000000
#BSUB -R rusage[mem=40000]
#BSUB -n 12
#BSUB -u jorn.pezoldt@epfl.ch

# Author by: Joern Pezoldt
# 12.10.2018
# Input:
# a) peak bed (e.g. regions not DEG but expressed and not DAR but open)
# Function:
# 1) obtain sequence in FASTA format
# Note: Run on VitalIT

# Modules to use
module use /software/module/;
module add UHTS/Analysis/BEDTools/2.26.0;


rootdir=/scratch/el/monthly/bdeplanc/pezoldt;
motifdir=$rootdir/Analysis/ATAC_FSC_all/homer/Motif_Input_R/Run_7_lowRNAseq_noExt/BED_peaks;
# backgrounddir=$rootdir/Analysis/ATAC_FSC_all/homer/Motif_Input_R/Run_7_lowRNAseq_noExt/background_all_TSS_minus_input_genes_lowRNAseq.bed;
resultdir=$rootdir/Analysis/ATAC_FSC_all/homer/Motif_Input_R/Run_7_lowRNAseq_noExt/Output_Regions;
genomedir=$rootdir/Data/Genome_mouse/GRCm38.87/Mus_musculus.GRCm38.87.dna.primary_assembly.fa
preparseddir=$rootdir/tmp;
name=`ls $motifdir | grep '.bed' | sed "${LSB_JOBINDEX}q;d"`;
sample=`ls $motifdir | grep '.bed' | sed "${LSB_JOBINDEX}q;d" | cut -d '.' -f 1`;

echo $sample
echo $name

#mkdir -p $resultdir/${sample};

# Generate Fasta files using bedtools
bedtools getfasta -fo $resultdir/${sample}.fa -fi $genomedir -bed $motifdir/${name}