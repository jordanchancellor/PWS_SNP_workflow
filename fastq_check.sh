#!/bin/bash

# a conda-dependent bash script to run fastq_utils paired fastq files 
# usage: bash scripts/fastq_check.sh
# assumes you have conda already installed and environments created for running fastqc and multiqc

# define resources and log files
# THREADS="12"
# MEM="10000"
your_email="jordan.chancellor.15@gmail.com"

# define global variables
wd="/home/jordan/project"
log_out="${wd}/fastq_utils.out"
log_err="${wd}/fastq_utils.err"
samples_file="${wd}/samples_file.txt"
readsdir="${wd}/reads/raw/F7SAN"
outdir="${wd}/reads_qc/"

# create output directory (if needed)
if [[ ! -d $outdir ]]
then
	printf "Creating output directory %s\n" "${outdir}"
	mkdir $outdir 
fi

# Load conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate fastqc

# fastq_utils checkpoint
fastq_info -h 2>> "$log_err" | tee -a "$log_out"
if [[ $? -ne 0 ]]; then
    echo "Error: fastq_utils not installed correctly." | tee -a "$log_err \n"
    exit 1
fi

# run fastq_info on fastq file pairs to check file integrity
while read -r sample_id; do
    read1="${readsdir}/${sample_id}_1.fastq.gz"
    read2="${readsdir}/${sample_id}_2.fastq.gz"

    # Check if both read files exist
    if [[ ! -f "$read1" || ! -f "$read2" ]]; then
        echo "Error: Missing files for ${sample_id}, skipping." | tee -a "$log_err"
        continue
    fi

    echo "Running fastq_info on files associated with ${sample_id}" | tee -a "$log_out"
    
    # Run fastq_info
    fastq_info "$read1" "$read2" >> "$log_out" 2>> "$log_err"
    if [[ $? -ne 0 ]]; then
        echo "Error: fastq_info failed for ${sample_id}" | tee -a "$log_err"
        exit 1
    fi
done < "$samples_file"

echo "fastq_utils finished running" | ssmtp "$your_email"
