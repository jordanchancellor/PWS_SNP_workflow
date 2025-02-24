#!/bin/bash

# a conda-dependent bash script to map fastq files to reference genome using bwa
# usage: bash PWS_SNP_workflow/alignment.sh
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
outdir="${wd}/alignment"
index="${wd}/indices/pvannamei/GCF_042767895.1_ASM4276789v1_genomic.fna.gz"

# Load conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate mapping

# create output directory (if needed)
if [[ ! -d $outdir ]]
then
	printf "Creating output directory %s\n" "${outdir}"
	mkdir $outdir 
fi

# BWA checkpoint
bwa mem -version 2>> "$log_err" | tee -a "$log_out"
if [[ $? -ne 0 ]]; then
    echo "Error: trimmomatic not installed correctly." | tee -a "$log_err"
    exit 1
fi

# align read pairs using bwa in batches 
while read -r sample_id; do
    read1="${readsdir}/${sample_id}_1.fastq.gz"
    read2="${readsdir}/${sample_id}_2.fastq.gz"

    # Check if both read files exist
    if [[ ! -f "$read1" || ! -f "$read2" ]]; then
        echo "Error: Missing files for ${sample_id}, skipping." | tee -a "$log_err"
        continue
    fi
    # Set up log files for each individual job
    job_out="${sample_id}_mapping.out"
    job_err="${sample_id}_mapping.err"

    # Run bwa on paired-end reads
    echo "Running bwa on files associated with ${sample_id}: ${read1} and ${read2}" | tee -a "$log_out"

    bwa mem -t "$THREADS" \
    "${index}" "${read1}" "${read2}" \
    > "${outdir}/${sample_id}.sam" 2> "$job_out" 3> "$job_err" &

    jobs+=($!)

    count=$((count + 1))

    if [[ $count -ge $batch_size ]]; then
        wait "${jobs[@]}"
        jobs=()
        count=0
    fi

done < "$samples_file"

wait
echo "BWA mapping completed. Mapping files can be found in ${outdir}." | tee -a "$log_out"
echo "BWA mapping completed. Mapping files can be found in ${outdir}." | ssmtp "$your_email"

# convert sam files to sorted bam files, remove sam files (for the sake of saving space)


############## below is scratch from an old script, I am working on updating the above

bwacmd="bwa mem -t $SLURM_CPUS_PER_TASK -R $(echo "@RG\tID:$id\tSM:$sm\tLB:$lb\tPL:ILLUMINA") $index $reads1 $reads2"

echo "$bwacmd" 

$bwacmd > "${outdir}/${sample_id}.sam"

conda deactivate
conda activate samtools

echo "Converting sam files to sorted bam files with samtools."

samfile="${outdir}/${sample_id}.sam"

cat ${samfile} | samtools view -bSh | samtools sort -o ${outdir}/${sample_id}.sorted.bam

echo "Creating samtools index file from ${sortedbamfile}"

samtools index -b ${sortedbamfile} 

echo "Running idxstats on sorted, indexed bam file from ${sortedbamfile}"

samtools idxstats ${sortedbamfile} | tee ${outdir}/${sample_id}.idxstats.txt

echo "Creating mapping rate file with samtools flagstat from ${bamfile}"

samtools flagstat ${sortedbamfile} > ${outdir}/${sample_id}.flagstat.txt

conda deactivate
conda activate gatk

# Run GATK4's ValidateSamFile in summary mode
gatk ValidateSamFile --INPUT ${sortedbamfile} \
--REFERENCE_SEQUENCE ${index} \
--MODE SUMMARY \
--OUTPUT ${outdir}/${sample_id}.validate.summary
