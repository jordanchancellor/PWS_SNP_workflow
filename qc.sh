#!/bin/bash

# a conda-dependent bash script to run fastqc and multiqc on paired fastq files 
# usage: bash scripts/qc.sh
# assumes you have conda already installed and environments created for running fastqc and multiqc

# define resources and log files
THREADS="48"
MEM="10000"
your_email="jordan.chancellor.15@gmail.com"

# define global variables
wd="/home/jordan/project"
log_out="${wd}/rawreads_fastqc.out"
log_err="${wd}/rawreads_fastqc.err"
# samples_file="${wd}/samples_file.txt"
readsdir="${wd}/reads/raw/F7SAN"
outdir="${wd}/reads_qc/fastqc_raw"

# create output directory (if needed)
if [[ ! -d $outdir ]]
then
	printf "Creating output directory %s\n" "${outdir}"
	mkdir $outdir 
fi

# Load conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate fastqc

# fastqc checkpoint
fastqc --version 2>> "$log_err" | tee -a "$log_out"
if [[ $? -ne 0 ]]; then
    echo "Error: FastQC not installed correctly." | tee -a "$log_err \n"
    exit 1
fi

##### run fastqc in loop for fastq files one-by-one
# while read -r sample_id; do
#     read1="${readsdir}/${sample_id}_1.fastq.gz"
#     read2="${readsdir}/${sample_id}_2.fastq.gz"

#     # Check if both read files exist
#     if [[ ! -f "$read1" || ! -f "$read2" ]]; then
#         echo "Error: Missing files for ${sample_id}, skipping." | tee -a "$log_err"
#         continue
#     fi

#     echo "Running FASTQC on files associated with ${sample_id}" | tee -a "$log_out"
    
#     # Run FastQC
#     fastqc -t "$THREADS" -o "$outdir" "$read1" "$read2" >> "$log_out" 2>> "$log_err" --memory="$FASTQC_MEMORY"
#     if [[ $? -ne 0 ]]; then
#         echo "Error: FASTQC failed for ${sample_id}" | tee -a "$log_err"
#         exit 1
#     fi
# done < "$samples_file"


##### run in parallel using ${readsdir} as input 

# check if output directory is empty, if so, process all fastq files in ${readsdir}
if [[ -z "$(ls -A "$outdir")" ]]; then
    echo "Output directory is empty. Running FastQC on all files in ${readsdir}." | tee -a "$log_out"
    fastqc -t "$THREADS" -o "$outdir" ${readsdir}/*.fastq.gz --memory ${MEM} >> "$log_out" 2>> "$log_err"
else
    echo "Output directory is not empty. Checking for unprocessed files..." | tee -a "$log_out"
    files_to_process=()

    for fastq in "${readsdir}"/*.fastq.gz; do
        base_name=$(basename "$fastq" .fastq.gz)
        if [[ ! -f "${outdir}/${base_name}_fastqc.html" ]]; then
            files_to_process+=("$fastq")
        fi
    done

    # ff no files need processing, exit
    if [[ ${#files_to_process[@]} -eq 0 ]]; then
        echo "All FASTQ files have already been processed. Exiting." | tee -a "$log_out"
        exit 0
    fi

    # run FastQC on unprocessed files
    echo "Running FastQC on the following files:" | tee -a "$log_out"
    printf "%s\n" "${files_to_process[@]}" | tee -a "$log_out"

    fastqc -t "$THREADS" -o "$outdir" "${files_to_process[@]}" --memory ${MEM} >> "$log_out" 2>> "$log_err"
fi

echo "FASTQC completed." | tee -a "$log_out"
echo "FASTQC finished running" | ssmtp "$your_email"

# run multiqc on fastqc results in output directory

# multiqc checkpoint
multiqc --version 2>> "$log_err" | tee -a "$log_out"
if [[ $? -ne 0 ]]; then
    echo "Error: MultiQC not installed correctly." | tee -a "$log_err"
    exit 1
fi

# If input directory exists, run MultiQC
echo "Running MultiQC on ${outdir}" | tee -a "$log_out"
multiqc -o "$outdir" "$outdir" >> "$log_out" 2>> "$log_err"
if [[ $? -ne 0 ]]; then
    echo "Error: MultiQC failed." | tee -a "$log_err"
    exit 1
fi

echo "MultiQC completed." | tee -a "$log_out"
echo "MultiQC finished running" | ssmtp "$your_email"
