#!/bin/bash
#SBATCH -J ResourceProcess
#SBATCH -p mcbstudent
#SBATCH -q mcbstudent
#SBATCH --mem=150G
#SBATCH -c 6
#SBATCH -o ./logs/%x_%A.out

echo "Hostname : " `hostname`

mkdir -p ../TMPDIR

export TMPDIR=../TMPDIR

cd ../resources

BOWTIE_DIR=/core/projects/GAP/CBC/ADelany_UCHC_miRNAseq_Apr2024/tools/bowtie-1.2.3

mkdir -p BowtieIndex
# Creating bowtie Index
${BOWTIE_DIR}/bowtie-build mm10.fa ./BowtieIndex/mm10

module load STAR/2.7.11a

mkdir -p StarIndex

# Creating STAR Index
STAR \
        --runThreadN 6 \
        --runMode genomeGenerate \
        --genomeDir ./StarIndex \
        --sjdbGTFfile mm10.gtf \
        --sjdbOverhang 1 \
        --genomeFastaFiles mm10.fa

