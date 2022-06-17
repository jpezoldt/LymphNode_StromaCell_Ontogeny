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

#Test for bitbucket sync
#get file names
#name=`ls $fastqdir/R1 | grep '.fastq.gz$' | sed "${LSB_JOBINDEX}q;d" | cut -d '_' -f 1-4`;
name=`ls $fastqdir/R1 | grep '.fastq.gz$' | sed "${LSB_JOBINDEX}q;d" | cut -d '_' -f 1-5`;

#sample names
#sample=`echo $name |cut -d '_' -f 1,2`;
sample=`echo $name |cut -d '_' -f 1-2`;

echo $name
echo $sample

mkdir -p $resultdir/aligned/${sample};
mkdir -p $resultdir/fastqc/${sample};
mkdir -p $resultdir/TSSplots;
mkdir -p $resultdir/bigwig;
mkdir -p $resultdir/homer;
mkdir -p $resultdir/aligned/$sample/tags;

# Align with STAR and sort reads
#STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --runThreadN 4 --genomeDir $indexdir --readFilesIn $fastqdir/R1/${sample}_1.fastq.gz $fastqdir/R2/${sample}_2.fastq.gz  --readFilesCommand zcat --outFilterMultimapNmax 1 --outFileNamePrefix $resultdir/aligned/${sample}/${sample}.

#mv $resultdir/aligned/${sample}/${sample}.Aligned.sortedByCoord.out.bam $resultdir/aligned/${sample}/${sample}.bam

# Generate index files
#samtools index $resultdir/aligned/${sample}/${sample}.bam 

#4. Mark and remove duplicates
#java -jar $rootdir/Scripts/picard.jar MarkDuplicates I=$resultdir/aligned/${sample}/${sample}.bam O=$resultdir/aligned/${sample}/${sample}.nodup.bam CREATE_INDEX=true  M= $resultdir/aligned/${sample}/${sample}.output.metrics REMOVE_DUPLICATES=true

#5.Fastqc on aligned and de-duplicated reads
#fastqc $resultdir/aligned/${sample}/${sample}.bam -o $resultdir/fastqc/${sample}
#fastqc $resultdir/aligned/${sample}/${sample}.nodup.bam -o $resultdir/fastqc/${sample}

#6. Heatmap of fragment distribution around TSS
#bamCoverage -b $resultdir/aligned/$sample/${sample}.nodup.bam -o $resultdir/bigwig/${sample}.bw -of "bigwig"
#computeMatrix reference-point -S $resultdir/bigwig/${sample}.bw -R $genomedir/Mus_musculus.GRCm38.87.gtf -a 3000 -b 3000 -o $resultdir/TSSplots/${sample}.matrixTSS.gz
#plotHeatmap -m $resultdir/TSSplots/${sample}.matrixTSS.gz -out $resultdir/TSSplots/${sample}.heatmap.png

#7.Peak calling 
mkdir -p $resultdir/peaks/default/${sample};
mkdir -p $resultdir/peaks/paired-end/${sample};
mkdir -p $resultdir/peaks/broad/${sample};
mkdir -p $resultdir/peaks/nolambda/${sample};
mkdir -p $resultdir/peaks/ext73/${sample};
mkdir -p $resultdir/peaks/ext150/${sample};

macs2 callpeak -t $resultdir/aligned/${sample}/${sample}.nodup.bam -f BAM -n ${sample}_default --outdir $resultdir/peaks/default/${sample} -g mm -q 0.05
macs2 callpeak -t $resultdir/aligned/${sample}/${sample}.nodup.bam -f BAMPE -n ${sample}_PE --outdir $resultdir/peaks/paired-end/${sample} -g mm -q 0.05
macs2 callpeak -t $resultdir/aligned/${sample}/${sample}.nodup.bam -f BAM --broad -n ${sample}_broad --outdir $resultdir/peaks/broad/${sample} -g mm -q 0.05
macs2 callpeak -t $resultdir/aligned/${sample}/${sample}.nodup.bam -f BAM --nolambda -n ${sample}_nolambda --outdir $resultdir/peaks/nolambda/${sample} -g mm -q 0.05
macs2 callpeak -t $resultdir/aligned/${sample}/${sample}.nodup.bam -f BAM -n ${sample}_ext73 --outdir $resultdir/peaks/ext73/${sample} --extsize 73 -g mm -q 0.05 
macs2 callpeak -t $resultdir/aligned/${sample}/${sample}.nodup.bam -f BAM -n ${sample}_ext150 --outdir $resultdir/peaks/ext150/${sample} --extsize 150 -g mm -q 0.05 

#Take .xls macs2 output and process with 15_...

