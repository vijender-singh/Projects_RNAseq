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

ADAPTERFILE="/isg/shared/apps/Trimmomatic/0.39/adapters/NexteraPE-PE.fa"


PROJECTDIR="/core/projects/GAP/CBC/ADelany_UCHC_miRNAseq_Apr2024/"

#SAMPLES=(sample-1 sample-2 sample-3)
SAMPLES=(oldfrac_1 oldfrac_2 oldfrac_3 oldfrac_4 oldfrac_5 youngfrac_6 youngfrac_8 youngfrac_10)

sampleID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

mkdir -p ${PROJECTDIR}/rawdata


cd ${PROJECTDIR}

mkdir ${PROJECTDIR}/02_qualityQC

##########################################################

cd ${PROJECTDIR}/02_qualityQC

module load fastqc/0.11.7

mkdir -p ./{RAWfastqc_OUT,TrimReads_trimmo,TrimReads_fastp,fastp_qc,TRIMfastqc_tm_OUT,TRIMfastqc_fp_OUT}

fastqc -t 6 -o ./RAWfastqc_OUT ../rawdata/${sampleID}_R1.fastq.gz #../rawdata/${sampleID}_R2.fastq.gz


mkdir -p ./{TrimReads_trimmo,TrimReads_fastp,fastp_qc}

module load Trimmomatic/0.39

java -jar $Trimmomatic SE -threads 6 \
        ../rawdata/${sampleID}_R1.fastq.gz \
        ./TrimReads_trimmo/trim_tm_${sampleID}_R1.fastq.gz \
        ILLUMINACLIP:${ADAPTERFILE}:2:30:10:5 \
        SLIDINGWINDOW:4:25 MINLEN:10

module rm Trimmomatic/0.39

module load fastp/0.23.2

fastp --thread 4 \
        --in1 ../rawdata/${sampleID}_R1.fastq.gz \
        --out1 ./TrimReads_fastp/trim_fp_${sampleID}_R1.fastq.gz \
        --json ./fastp_qc/${sampleID}_fastp.json \
        --html ./fastp_qc/${sampleID}_fastp.html


module load fastqc

mkdir -p ./TRIMfastqc_tm_OUT
mkdir -p ./TRIMfastqc_fp_OUT

fastqc -t 6 -o ./TRIMfastqc_tm_OUT ./TrimReads_trimmo/trim_tm_${sampleID}_R1.fastq.gz

fastqc -t 6 -o ./TRIMfastqc_fp_OUT ./TrimReads_fastp/trim_fp_${sampleID}_R1.fastq.gz

module purge
