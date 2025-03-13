#!/usr/bin/env bash
#SBATCH -J cleanbams
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=2-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH --array=1-45

# a module-dependent bash script to clean bam files based on flags in output GATK4's ValidateSamFile
# requires input list of bam files to process in array
# usage: sbatch PWS_SNP_workflow/cleanbamfiles.sbatch <flag>

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
flag=$1
wd="/mnt/scratch/users/agr_3216_ifapa/alorenzo/vannamei_snp_workflow"
samples_file="${wd}/samples_file.txt"
sample_id=$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)
outdir="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/alignment"
reference="${wd}/indices/pvannamei/GCF_042767895.1_ASM4276789v1_genomic.fna" # Genome assembly ASM4276789v1
sortedbamfile="${outdir}/${sample_id}_sorted.bam"
cleanbam="${outdir}/${sample_id}_sorted_cleaned.bam"

if [ -z "$1" ]; then
    echo "Error: No flag specified."
    echo "Usage: PWS_SNP_workflow/cleanbamfiles.sbatch <flag>"
    echo "Options: INVALID_TAG_NM, "
    exit 1
fi

# create output directory (if needed)
if [[ ! -d $outdir ]]; then
	printf "Creating output directory %s\n" "${outdir}"
	mkdir $outdir 
fi

# clean bam files according to specified flags

case "$flag" in
    INVALID_TAG_NM) 
        echo "Running calmd on bam files using samtools"
        module purge
        module load samtools/1.21
        # clean "INVALID_TAG_NM" flag with samtools calmd
        samtools calmd -bAr ${sortedbamfile} ${reference} > ${cleanbam}
        echo "samtools calmd completed. Cleaned bam files can be found in ${outdir}."
        ;;
    *)
    echo "Error: No flag specified."
    echo "Usage: PWS_SNP_workflow/cleanbamfiles.sbatch <flag>"
    echo "Options: INVALID_TAG_NM, "
    exit 1
        ;;
esac

# re-validate bam files using GATK4's ValidateSamFile

if [[ -f ${cleanbam} ]]
then
	echo "Validating ${cleanbam} with GATK's ValidateSamFile."
    module purge
    module load gatk/4.4.0.0
    gatk ValidateSamFile --INPUT ${cleanbam} \
    --REFERENCE_SEQUENCE ${reference} \
    --MODE SUMMARY \
    --OUTPUT ${outdir}/${sample_id}_cleaned_validate.summary
else
    echo "Error in cleaning bam files. GATK ValidateSamFile not re-run"
fi

