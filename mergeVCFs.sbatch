#!/usr/bin/env bash
#SBATCH -J mergeVCFs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# a module-dependent bash script to merge VCF files using GATK4
# requires that gvcf files were imported by interval using genomicsDBimport.sbatch and genotyped using genotypeGVCFs.sbatch
# usage: sbatch PWS_SNP_workflow/mergeVCFs.sbatch

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
wd="/mnt/scratch/users/agr_3216_ifapa/alorenzo/vannamei_snp_workflow"
reference="${wd}/indices/pvannamei/GCF_042767895.1_ASM4276789v1_genomic.fna" # Genome assembly ASM4276789v1
ref_path=$(dirname "$reference")
index=$(basename "${reference}" .fna)
outdir="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/variants"
tempdir="${outdir}/tmp"
input_list="${outdir}/vcf_files.list"

# Load  environment
module purge
module load samtools/1.21
module load gatk/4.4.0.0

# check GATK4 install
gatk --version
if [[ $? -ne 0 ]]
then
	printf "Error - check GATK4 installation.\n"
	exit 1
fi

# create output directories (if needed)
if [[ ! -d $outdir ]]
then
	mkdir $outdir
fi
if [[ ! -d $tempdir ]]
then
        mkdir $tempdir
fi

# Run GATK4's mergeVcfs
gatk --java-options "-Xmx80G" MergeVcfs \
    -R ${reference} \
	--SEQUENCE_DICTIONARY ${ref_path}/${index}.dict \
    -I ${input_list} \
    -O ${outdir}/merged_variants.vcf.gz \
    --TMP_DIR ${tempdir}

echo "VCF files merged. Final VCF file can be found at ${outdir}/merged_variants.vcf.gz."

# Validate Merged VCF file
echo "Validating merged VCF file"
 gatk ValidateVariants \
   -R ${reference}\
   -V ${outdir}/merged_variants.vcf.gz \
   --validation-type-to-exclude ALLELES > \
    ${outdir}/merged_variants.validate.summary 2>&1