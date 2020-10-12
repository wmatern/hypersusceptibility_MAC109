# Tn-seq analysis of Mycobacterium Avium antibiotic hypersusceptibility
This workflow, implemented in a Jupyter notebook and two scripts (make_barplots.py and make_venn.R), seeks to reproduce the results presented in our upcoming paper: "Genetic determinants of antibiotic hypersusceptibility in Mycobacterium avium". Following the steps below should allow you to get identical results to those reported in our paper.

## Authors

* William Matern (@wmatern)

## Dependencies

You will need to install the following packages (and their dependencies) in order to run this workflow:
* conda (tested with version: 4.8.3)
* R and eulerr

## Raw Data
The necessary Tn-seq data can be found in NCBI SRA under BioProject: PRJNA559896. You will need to add all of the associated fastq files (SRR9953606-SRR953664, 118 files, with files ending in \*\_1.fastq and \*\_2.fastq) into the input/ folder before running the Jupyter notebook. Additionally, the genome information contained in the Genbank file (.gb, .gbff) along with corresponding FASTA files (.fa, .fna) can be downloaded at ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/408/535/GCA_003408535.1_ASM340853v1. You will need to rename the genome information files to use the prefix MAC109. 

The expected file names in order to run this workflow are:

    input/Genome/MAC109.fa
    input/Genome/MAC109.gb

You can view input/Samples\_Labels.csv to see the list of sample information.

## Usage

### Step 1: Download workflow
If you simply want to run this workflow, download and extract the [latest release](https://github.com/).

### Step 2: Execute workflow
Run each cell of Jupyter notebook in order.

### Step 3: Run make_barplots.py to create barplots
`python3 make_barplots.py`

### Step 4: Change working directory and run make_venn.R
`vim make_venn.R`
Set the first line of make_venn.R to the appropriate location.
`Rscript make_venn.R`

## Report Issue
If you have any questions or issues reproducing our results please send me an email at maternwill@gmail.com.
