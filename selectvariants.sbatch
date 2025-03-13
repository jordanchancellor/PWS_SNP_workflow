#!/usr/bin/env bash
#SBATCH -J selectvariants
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# a module-dependent bash script to select variants using GATK4
# usage: sbatch PWS_SNP_workflow/variantmetrics.sbatch <vcf_file> <variant_type> <output_file>
# output is a .vcf.gz file containing only selected variants. Include file ending in <output_file>

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
variant_file=$1
variant_type=$2
output_file=$3

if [ "$#" -ne 3 ]; then
    echo "Error: input files not specified."
    echo "Usage: sbatch PWS_SNP_workflow/variantmetrics.sbatch <vcf_file> <variant_type> <output_file>"
    echo "variant_type options: INDEL, SNP, MIXED, MNP, SYMBOLIC, NO_VARIATION"
    exit 1
fi

# Load  environment
module purge
module load gatk/4.4.0.0

# Select variants
gatk --java-options "-Xmx80G" SelectVariants \
    -V ${variant_file} \
    -select-type ${variant_type} \
    -O ${output_file}