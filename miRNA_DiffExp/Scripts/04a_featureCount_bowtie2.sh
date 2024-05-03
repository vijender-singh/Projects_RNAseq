#!/bin/bash
#SBATCH --job-name=count_bt
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH -o ./logs/%x_%A.out

echo "HOSTNAME: `hostname`"
echo "Start Time: `date`"

mkdir -p ./.tmp

export TMPDIR=./.tmp


PROJECTDIR="/core/projects/GAP/CBC/ADelany_UCHC_miRNAseq_Apr2024/"

#SAMPLES=(sample-1 sample-2 sample-3)
SAMPLES=(oldfrac_1 oldfrac_2 oldfrac_3 oldfrac_4 oldfrac_5 youngfrac_6 youngfrac_8 youngfrac_10)

miRNAgff3=../resources/mmu.gff3

bam=""
for i in ${SAMPLES[@]}; do  bam="${bam} ${i}_mapped_sort.bam"; done

OUTDIR=../05_featureCount
mkdir -p ${OUTDIR}

cd ../03-mapping

module load subread/2.0.3

featureCounts -F GTF \
        -t miRNA \
        -g 'Name' \
        -a ${miRNAgff3} \
        -o ${OUTDIR}/bt-counts.txt \
        ${bam}



