#!/usr/bin/env bash
#SBATCH -J filterSNPs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# a module-dependent bash script to filter monomorphic, biallelic SNPs as well as input MAF
# usage: sbatch PWS_SNP_workflow/filterSNPs.sbatch <vcf_file> <MAF>
# output is a filtered .vcf.gz file

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
variant_file=$1
outdir=$(dirname "$variant_file")
MAF=$2

if [ "$#" -ne 2 ]; then
    echo "Error: input files not specified."
    echo "Usage: sbatch PWS_SNP_workflow/filterSNPs.sbatch <vcf_file> <MAF>"
    exit 1
fi

# Load  environment
module purge
module load bcftools/1.16

# Remove monomorphic SNPs (retain only polymorphic sites): reference-only variants (AC==0), fixed alternative variants (AC == AN)
# Retain only bielleic snps
bcftools filter \
-e 'AC==0 || AC==AN' ${variant_file} -O z -o out1.vcf.gz | \
bcftools view -m2 -M2 -v snps -O z -o ${outdir}/biallelic_qualfiltered_snps.vcf.gz

# Get stats file
bcftools stats -s - ${outdir}/biallelic_qualfiltered_snps.vcf.gz > ${outdir}/bcftools_stats_biallelic_qualfiltered_snps.txt

# filter for MAF (set in command)
bcftools filter \
-e "MAF <= ${MAF}" ${outdir}/biallelic_qualfiltered_snps.vcf.gz \
-O z -o ${outdir}/MAF${MAF}_biallelic_qualfiltered_snps.vcf.gz

# Get stats file
bcftools stats -s - ${outdir}/MAF${MAF}_biallelic_qualfiltered_snps.vcf.gz > ${outdir}/bcftools_stats_MAF${MAF}_biallelic_qualfiltered_snps.txt
