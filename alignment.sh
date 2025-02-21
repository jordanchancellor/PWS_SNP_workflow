#!/bin/bash

# a conda-dependent bash script to run trimmoimatic on paired fastq files 
# usage: bash scripts/trimming.sh
# assumes you have conda already installed and environments created

# define resources and batch sizes (us batch size of 1 to run one file at a time in loop)
THREADS="2"
your_email="jordan.chancellor.15@gmail.com"
batch_size=24 
count=0
jobs=()

# define global variables
wd="/home/jordan/project"
log_out="${wd}/trimmomatic.out"
log_err="${wd}/trimmomatic.err"
samples_file="${wd}/samples_file.txt"
readsdir="${wd}/reads/trimmed"
outdir="${wd}/"

# Load conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate trimming

# create output directory (if needed)
if [[ ! -d $outdir ]]
then
	printf "Creating output directory %s\n" "${outdir}"
	mkdir $outdir 
fi

# BWA checkpoint
trimmomatic --version 2>> "$log_err" | tee -a "$log_out"
if [[ $? -ne 0 ]]; then
    echo "Error: trimmomatic not installed correctly." | tee -a "$log_err"
    exit 1
fi
