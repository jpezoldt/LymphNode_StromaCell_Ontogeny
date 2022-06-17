#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_FSC_all/std_err_out/ATAC_motifs_in_regions.err
#BSUB -o /scratch/el/monthly/bdeplanc/pezoldt/Analysis/ATAC_FSC_all/std_err_out/ATAC_motifs_in_regions.out
#BSUB -J Motif[1-8]
#BSUB -M 40000000
#BSUB -R rusage[mem=40000]
#BSUB -n 8
#BSUB -u joern.pezoldt@epfl.ch

# Author by: Joern Pezoldt
# 17.10.2018
# Input:
# a) background bed (e.g. regions not DEG but expressed and not DAR but open)
# Function:
# 1) Find Motifs 
# Note: Run on VitalIT

module use /software/module/;
#homer 4.2 is not on vital-it anymore
module add UHTS/Analysis/homer/4.6;


#parameter=si500_v46_bgTSS_noExt_;
rootdir=/scratch/el/monthly/bdeplanc/pezoldt;
beddir=$rootdir/Analysis/ATAC_FSC_all/homer/Find_Region_by_Motif/Input/BED_Peak;
motifdir=$rootdir/Analysis/ATAC_FSC_all/homer/Find_Region_by_Motif/Input/Motif_Files;
backgrounddir=$rootdir/Analysis/ATAC_FSC_all/homer/Motif_Input_R/Run_8_low_all/background_all_TSS_minus_input_genes_lowRNAseq.bed;
resultdir=$rootdir/Analysis/ATAC_FSC_all/homer/Find_Region_by_Motif/Output;
genomedir=$rootdir/Data/Genome_mouse/GRCm38.87/Mus_musculus.GRCm38.87.dna.primary_assembly.fa
preparseddir=$rootdir/tmp;
name=`ls $beddir | grep '.bed' | sed "${LSB_JOBINDEX}q;d"`;
sample=`ls $beddir | grep '.bed' | sed "${LSB_JOBINDEX}q;d" | cut -d '.' -f 1`;
motif_name=`ls $motifdir | grep '_compiled.txt' | sed "${LSB_JOBINDEX}q;d"`;

echo $sample
echo $name
echo $motif_name

mkdir -p $resultdir/${sample};

#Find Motifs in regions
findMotifsGenome.pl $beddir/${name} $genomedir  $resultdir/${sample} -find $motifdir/${motif_name} -bg $backgrounddir -size 500 -len 10 > $resultdir/${sample}/${sample}_MotifInstances.txt

#Annotate to peak
annotatePeaks.pl $beddir/${name} $genomedir -m $motifdir/${motif_name} > $resultdir/${sample}/${sample}_MotifInstances_Location.txt

