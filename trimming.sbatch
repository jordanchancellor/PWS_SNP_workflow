#!/usr/bin/env bash
#SBATCH -J trim
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=05:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=jordan.chancellor.15@gmail.com
#SBATCH --array=1-45

# a module-dependent bash script array to run trimmoimatic on paired fastq files from input samples_file.txt sample list
# usage: sbatch PWS_SNP_workflow/trimming.sbatch
# set --array vairable to number of samples you wish to process

# set global variables
THREADS=$SLURM_CPUS_PER_TASK
wd="/mnt/scratch/users/agr_3216_ifapa/alorenzo/vannamei_snp_workflow"
samples_file="${wd}/samples_file.txt"
sample_id=$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)
reads="${wd}/reads/raw"
outdir="${wd}/reads/trimmed"
read1="${reads}/${sample_id}_1.fastq.gz"
read2="${reads}/${sample_id}_2.fastq.gz"

# define desired trimming options
phred="phred33"
ADAPTERS="/mnt/home/soft/trimmomatic/programs/x86_64/0.39/adapters/TruSeq3-PE.fa:2:30:10"
LEADING="3"
TRAILING="3"
SLIDINGWINDOW="4:20"
MINLEN="50"

# load enrivonment
module load trimmomatic/0.36

# Create output directories (if needed)
if [[ ! -d $outdir ]]
then
	printf "Creating output directory ${outdir}."
	mkdir $outdir 
fi

# trimmomatic checkpoint
trimmomatic -version
if [[ $? -ne 0 ]]
then
	printf "Error - check trimmomatic installation.\n"
	exit 1
fi

# check if both read files exist
    if [[ ! -f "$read1" || ! -f "$read2" ]]; then
        echo "Error: Missing files for ${sample_id}, skipping." | tee -a "$log_err"
        continue
    fi

# run trimmomatic on paired reads

echo "Running Trimmomatic on paired reads associated with ${sample_id}: ${read1} and ${read2}"

trimmomatic PE -threads "$THREADS" "-${phred}" \
    "${read1}" "${read2}" \
    "${outdir}/${sample_id}_F_trimmed_paired.fq.gz" "${outdir}/${sample_id}_F_trimmed_unpaired.fq.gz" \
    "${outdir}/${sample_id}_R_trimmed_paired.fq.gz" "${outdir}/${sample_id}_R_trimmed_unpaired.fq.gz" \
    ILLUMINACLIP:"${ADAPTERS}" \
    LEADING:"${LEADING}" \
    TRAILING:"${TRAILING}" \
    SLIDINGWINDOW:"${SLIDINGWINDOW}" \
    MINLEN:"${MINLEN}"

echo "PE Trimming complete. Trimmed files can be found in ${outdir}."