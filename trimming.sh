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
readsdir="${wd}/reads/raw/F7SAN"
outdir="${wd}/reads/trimmed"


# define desired trimming options
phred="phred33"
ADAPTERS="/home/jordan/anaconda3/envs/trimming/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10"
LEADING="3"
TRAILING="3"
SLIDINGWINDOW="4:20"
MINLEN="50"

# Load conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate trimming

# create output directory (if needed)
if [[ ! -d $outdir ]]
then
	printf "Creating output directory %s\n" "${outdir}"
	mkdir $outdir 
fi

# trimmomatic checkpoint
trimmomatic --version 2>> "$log_err" | tee -a "$log_out"
if [[ $? -ne 0 ]]; then
    echo "Error: trimmomatic not installed correctly." | tee -a "$log_err"
    exit 1
fi

# run trimmomatic on paired reads in loop one-by-one 
while read -r sample_id; do
    read1="${readsdir}/${sample_id}_1.fastq.gz"
    read2="${readsdir}/${sample_id}_2.fastq.gz"

    # Check if both read files exist
    if [[ ! -f "$read1" || ! -f "$read2" ]]; then
        echo "Error: Missing files for ${sample_id}, skipping." | tee -a "$log_err"
        continue
    fi
    # Set up log files for each individual job
    job_out="${sample_id}_trimmomatic.out"
    job_err="${sample_id}_trimmomatic.err"

    # Run trimmomatic on paired-end reads
    echo "Running trimmomatic on files associated with ${sample_id}: ${read1} and ${read2}" | tee -a "$log_out"

    trimmomatic PE -threads "$THREADS" -"${phred}" \
    "${read1}" "${read2}" \
    ${outdir}/"${sample_id}"_F_trimmed_paired.fq.gz ${outdir}/"${sample_id}"_F_trimmed_unpaired.fq.gz \
    ${outdir}/"${sample_id}"_R_trimmed_paired.fq.gz ${outdir}/"${sample_id}"_R_trimmed_unpaired.fq.gz \
    ILLUMINACLIP:"${ADAPTERS}" \
    LEADING:"${LEADING}" \
    TRAILING:"${TRAILING}" \
    SLIDINGWINDOW:"${SLIDINGWINDOW}" \
    MINLEN:"${MINLEN}" > "$job_out" 2> "$job_err" &

    jobs+=($!)

    count=$((count + 1))

    if [[ $count -ge $batch_size ]]; then
        wait "${jobs[@]}"
        jobs=()
        count=0
    fi

done < "$samples_file"

wait
echo "PE Trimming complete. Trimmed files can be found in ${outdir}." | tee -a "$log_out"
echo "PE Trimming complete. Trimmed files can be found in ${outdir}." | ssmtp "$your_email"
