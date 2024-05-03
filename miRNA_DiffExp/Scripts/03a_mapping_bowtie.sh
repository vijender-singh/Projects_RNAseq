#!/bin/bash
#SBATCH --job-name=Map2Count
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --partition=xeon
#SBATCH --qos=general
#SBATCH --array=[0-7]%4
#SBATCH -o ./logs/%x_%A_%a.out

echo "HOSTNAME: `hostname`"
echo "Start Time: `date`"

mkdir -p ./.tmp

export TMPDIR=./.tmp


PROJECTDIR="/core/projects/GAP/CBC/ADelany_UCHC_miRNAseq_Apr2024/"

#SAMPLES=(sample-1 sample-2 sample-3)
SAMPLES=(oldfrac_1 oldfrac_2 oldfrac_3 oldfrac_4 oldfrac_5 youngfrac_6 youngfrac_8 youngfrac_10)

sampleID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

BOWTIE_DIR=/core/projects/GAP/CBC/ADelany_UCHC_miRNAseq_Apr2024/tools/bowtie-1.2.3

mkdir -p ../03-mapping/tmp_${sampleID}

cd ../03-mapping/tmp_${sampleID}

#${BOWTIE_DIR}/bowtie-build mm10.fa ./BowtieIndex/mm10
#${BOWTIE_DIR}/bowtie 

${BOWTIE_DIR}/bowtie-align-s -p 6 \
	../../resources/BowtieIndex/mm10 \
	../../02_qualityQC/TrimReads_fastp/trim_fp_${sampleID}_R1.fastq.gz \
	-S ${sampleID}.sam

module load samtools/1.9

samtools view -@ 6 -bhS ${sampleID}.sam -o ${sampleID}_mapped.bam

samtools sort -@ 6 ${sampleID}_mapped.bam -o ${sampleID}_mapped_sort.bam


mv ${sampleID}_mapped_sort.bam ../

