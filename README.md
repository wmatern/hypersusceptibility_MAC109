# Tn-seq analysis of Mycobacterium Avium antibiotic hypersusceptibility
This workflow, implemented in a Jupyter notebook, seeks to reproduce the results presented in our upcoming paper: "Genetic determinants of antibiotic hypersusceptibility in Mycobacterium avium". Following the steps below should allow you to get identical results to those reported in our paper.

## Authors

* William Matern (@wmatern)

## Dependencies

You will need to install the following packages (and their dependencies) in order to run this workflow:
* conda (tested with version: 4.8.3)

## Raw Data
You can find the Tn-seq data in NCBI's SRA under BioProject: PRJNA559896. You will need to download a number of data files and place them into input folder. A list of commands using the SRA Toolkit is provided in each input folder as a guide for how to download .fastq files from SRA. Genbank file (.gb, .gbff) along with corresponding FASTA files (.fa, .fna) can be downloaded at ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/408/535/GCA_003408535.1_ASM340853v1. You will need to rename the files to use the prefix MAC109. The expected file names in order to run this workflow are:

    input/Genome/MAC109.fa
    input/Genome/MAC109.gb
    input/.fastq
    input/.fastq
    input/.fastq
    input/.fastq
    input/.fastq
    input/.fastq
    input/.fastq
    input/.fastq
    input/.fastq
    input/.fastq

## Usage

### Step 1: Download workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/).

### Step 2: Execute workflow

## Report Issue

If you have any questions or issues reproducing our results please send me an email at maternwill@gmail.com.
