#!/usr/bin/env bash
#SBATCH -J genotypeGVCFs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH --array=1-43

# a module-dependent bash script to perform joint genotyping on a multi-sample GVCF from GenomicsDBImportjoint genotype intervals using GATK4's genotypeGVCFs
# requires that gvcf files were imported by interval using genomicsDBimport.sbatch
# set array size to number of intervals: $(ls ${split_intervals_dir}/*list | wc -l)
# usage: sbatch PWS_SNP_workflow/genotypeGVCFs.sbatch

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
wd="/mnt/scratch/users/agr_3216_ifapa/alorenzo/vannamei_snp_workflow"
intervals_dir="${wd}/indices/pvannamei/splitintervals"
interval=$(ls ${intervals_dir}/*.interval_list | sed -n ${SLURM_ARRAY_TASK_ID}p)
reference="${wd}/indices/pvannamei/GCF_042767895.1_ASM4276789v1_genomic.fna" # Genome assembly ASM4276789v1
ref_path=$(dirname "$reference")
index=$(basename "${reference}" .fna)
outdir="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/variants"
tempdir="${outdir}/tmp"
# gvcf_list="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/alignment/F7SAN_gvcf.sample.map"

# Load  environment
module purge
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

# Run GATK4's genotypeGVCFs
gatk --java-options "-Xmx80G" GenotypeGVCFs \
    -R ${reference} \
	-L "${interval}" \
	--sequence-dictionary ${ref_path}/${index}.dict \
    -V gendb://${outdir}/interval_${SLURM_ARRAY_TASK_ID} \
    -O ${outdir}/interval_variants_${SLURM_ARRAY_TASK_ID}.vcf.gz \
	--tmp-dir ${tempdir}
