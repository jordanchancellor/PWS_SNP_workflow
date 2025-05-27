#!/usr/bin/env bash
#SBATCH -J simprocessing
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# usage: sbatch PWS_SNP_workflow/simulationprocessing.sbatch

wd="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/"
vcf="${wd}/variants/qualfiltered_snps.vcf.gz"
outdir="${wd}/simulations"
chroms="${outdir}/GCF_042767895.1_ASM4276789v1_chroms.bed"

module purge
module load bcftools/1.16

# select biallelic snps
bcftools view -m2 -M2 ${vcf} -O z -o ${outdir}/biallelic.vcf.gz

# select only chromosomes

bcftools index ${outdir}/biallelic.vcf.gz

bcftools view -R ${chroms} ${outdir}/biallelic.vcf.gz -O z -o ${outdir}/biallelic_chromsonly.vcf.gz

bcftools index ${outdir}/biallelic_chromsonly.vcf.gz