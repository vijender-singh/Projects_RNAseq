## Repository Details

Welcome to this GitHub repository dedicated to RNAseq data analysis using R! 

Repository Structure:

    Each project within the repository is neatly organized into its own directory, containing scripts in Markdown format for easy readability and instruction. Accompanying each Markdown file, you will find its corresponding HTML and PDF formats, allowing for accessibility and convenience in reviewing the analysis steps and outcomes.

Workflow Overview:

    Initial Stages (Pre-processing and Quality Control): The journey from raw fastq files to the generation of counts is orchestrated using Nextflow, a workflow tool that enables scalable and reproducible scientific workflows. We start with quality control checks using FastQC, followed by trimming and cleaning the reads using a choice of Trim Galore/Cutadapt/Fastp based on the dataset's specific needs. A second round of FastQC ensures the high quality of our data post-processing.
    Alignment and Quantification: The clean reads are then aligned to the reference genome using Hisat2, creating SAM files that are further processed by Samtools. For quantification, we employ either HTSeq-count or StringTie, depending on the downstream analysis requirements, to generate counts data that serve as the foundation for the subsequent analysis stages.

Downstream Analysis in R:

    The core of our repository lies in the downstream analysis performed in R. Here, we delve into the counts data to uncover the biological significance behind the numbers. Our scripts guide users through normalization, differential expression analysis, and various methods of visualization and interpretation of the results, leveraging the rich ecosystem of R packages tailored for genomics data analysis.
