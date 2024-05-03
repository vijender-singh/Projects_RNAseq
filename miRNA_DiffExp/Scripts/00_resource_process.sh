#!/bin/bash
#SBATCH -J ResourceProcess
#SBATCH -p general
#SBATCH -q general
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH -o ./logs/%x_%A.out

echo "Hostname : " `hostname`

mkdir -p ../TMPDIR

export TMPDIR=../TMPDIR

###############################
# module load 
###############################
module load samtools/1.16.1
module load seqtk/1.3

mkdir -p ../00_resources

cd ../00_resources

###############################
#     GENOME
###############################
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
# Uncompress Genome fasta file
gunzip mm10.fa.gz

# Select complete chromosmes (chr1...chr19,chrX and chrY)
# Selected chromosomes are listed in chromosomes_selected.txt with one chromosome ID per line.
# It is worth removing MT sequences here if study doesnot involve MT features.
seqtk mm10.fa chromosomes_selected.txt > mm10_selected.fa

# Deleting the downloaded fasta file.
rm mm10.fa

# Renaming selected sequence fasta file
mv mm10_selected.fa mm10.fa

# Index fasta file
samtools faidx mm10.fa

#Get chromosome sizes
cut -f1,2 mm10.fa.fai > chromosome.sizes.genome

# Download miRNA gff3 file from miRBAse
wget https://www.mirbase.org/download/mmu.gff3

