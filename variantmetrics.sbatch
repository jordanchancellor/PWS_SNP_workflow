#!/usr/bin/env bash
#SBATCH -J variantmetrics
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# a module-dependent bash script to collect variant calling metrics using Picard within GATK4
# requires that reference genome dictionary has been created ahead of time
# usage: sbatch PWS_SNP_workflow/variantmetrics.sbatch <vcf_file> <reference_fasta> <output_file>

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
wd="/mnt/scratch/users/agr_3216_ifapa/alorenzo/vannamei_snp_workflow"
variant_file=$1
prefix=$(basename "${variant_file}")
reference=$2
# reference="${wd}/indices/pvannamei/GCF_042767895.1_ASM4276789v1_genomic.fna" # Genome assembly ASM4276789v1
ref_path=$(dirname "$reference")
index=$(basename "${reference}" .fna)
output=$3
# outdir="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/variants"

if [ "$#" -ne 3 ]; then
    echo "Error: input files not specified."
    echo "Usage: sbatch PWS_SNP_workflow/variantmetrics.sbatch <vcf_file> <reference_fasta> <output_file>"
    exit 1
fi

# Load  environment
module purge
module load gatk/4.4.0.0

# collect variant calling 
gatk --java-options "-Xmx80G" CollectVariantCallingMetrics \
    -I ${variant_file} \
    -SD ${ref_path}/${index}.dict \
    -O ${output}