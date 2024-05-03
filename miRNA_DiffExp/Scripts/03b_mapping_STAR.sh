#!/bin/bash
#SBATCH --job-name=Star_align
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

mkdir -p ../03-mapping/star_${sampleID}

cd ../03-mapping/star_${sampleID}

module load STAR/2.7.11a

params=' --runThreadN 6
--sjdbGTFfile ../../resources/mm10.gtf
--alignEndsType EndToEnd
--outFilterMismatchNmax 1
--outFilterMultimapScoreRange 0
--quantMode TranscriptomeSAM GeneCounts
--outReadsUnmapped Fastx
--outFilterMultimapNmax 10
--outSAMunmapped Within
--outFilterScoreMinOverLread 0
--outFilterMatchNminOverLread 0
--outFilterMatchNmin 16
--alignSJDBoverhangMin 1000
--alignIntronMax 1
--outWigType wiggle
--outWigStrand Stranded
--outWigNorm RPM
'
STAR --genomeDir ../../resources/StarIndex --readFilesIn ../../02_qualityQC/TrimReads_fastp/trim_fp_${sampleID}_R1.fastq $params --outFileNamePrefix ${sampleID} --outSAMtype BAM SortedByCoordinate

mv *.bam ../
