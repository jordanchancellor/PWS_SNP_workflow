#!/usr/bin/env bash
#SBATCH -J qc_trimmed
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH --array=1-45

# a module-dependent bash script array to run fastqc and multiqc on paired fastq files from input samples_file.txt sample list
# usage: sbatch PWS_SNP_workflow/qc.sbatch
# set --array vairable to number of samples you wish to process

# set global variables
THREADS=$SLURM_CPUS_PER_TASK
wd="/mnt/scratch/users/agr_3216_ifapa/alorenzo/vannamei_snp_workflow"
samples_file="${wd}/samples_file.txt"
sample_id=$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)
reads="${wd}/reads/trimmed"
outdir="${wd}/reads_qc/fastqc_trimmed"
read1="${reads}/${sample_id}_F_trimmed_paired.fq.gz"
read2="${reads}/${sample_id}_R_trimmed_paired.fq.gz"

# Load  environment
module purge
module load fastqc/0.11.9
module load multiqc/1.13a

# create output directory (if needed)
if [[ ! -d "$outdir" ]]
then
    printf "Creating output directory %s\n" "${outdir}"
    mkdir "$outdir"
fi

# fastqc checkpoint
fastqc --version
if [[ $? -ne 0 ]]; then
    echo "Error: FastQC not installed correctly."
    exit 1
fi

# check if both read files exist
if [[ ! -f "$read1" || ! -f "$read2" ]]; then
    echo "Error: Missing files for ${sample_id}."
    exit 1
fi

# run fastqc
echo "Running FASTQC on files associated with ${sample_id}"

fastqc -t "$THREADS" -o "$outdir" "$read1" "$read2"

# check exit status
if [[ $? -ne 0 ]]; then
    echo "Error: FASTQC failed for ${sample_id}"
    exit 1
else
    # log the existence of the FASTQC output files for this sample and track file existence
    expected_files=(
        "${outdir}/${sample_id}_fastqc.html"
        "${outdir}/${sample_id}_fastqc.zip"
    )
    for file in "${expected_files[@]}"; do
        if [[ -f "$file" ]]; then
            echo "${sample_id} $file exists" >> "$log_fastq_file_tracking"
        else
            echo "${sample_id} $file missing" >> "$log_fastq_file_tracking"
        fi
    done
fi

# run multiqc if final array job
if [[ "$SLURM_ARRAY_TASK_ID" -eq "$SLURM_ARRAY_TASK_COUNT" ]]; then
    echo "Final job: Checking if all FASTQC files exist before running MultiQC..."

    # Check if all expected files are present
    missing_files=$(grep "missing" "$log_fastq_file_tracking" | wc -l)

    if [[ "$missing_files" -eq 0 ]]; then
        # multiqc checkpoint
        multiqc --version
        if [[ $? -ne 0 ]]; then
            echo "Error: MultiQC not installed correctly."
            exit 1
        fi

        # check that outdir exists
        if [[ -d "$outdir" ]]; then
            echo "Running MultiQC on ${outdir}"
            multiqc -o "$outdir" "$outdir"
            echo "MultiQC completed. Results can be found in ${outdir}."
        else
            echo "Error: Output directory ${outdir} not found. MultiQC failed."
            exit 1
        fi
    else
        echo "Some expected FASTQC output files are missing. MultiQC will not run."
        exit 1
    fi
fi