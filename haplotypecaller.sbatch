#!/usr/bin/env bash
#SBATCH -J haplotypecaller
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=7-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH --array=1-45

# a module-dependent bash script to genotype snps using GATK4's HaplotypeCaller
# usage: sbatch PWS_SNP_workflow/haplotypecaller.sbatch

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
wd="/mnt/scratch/users/agr_3216_ifapa/alorenzo/vannamei_snp_workflow"
samples_file="${wd}/samples_file.txt"
sample_id=$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)
outdir="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/alignment"
reference="${wd}/indices/pvannamei/GCF_042767895.1_ASM4276789v1_genomic.fna" # Genome assembly ASM4276789v1
ref_path=$(dirname "$reference")
index=$(basename "${reference}" .fna)
markedsortedbamfile="${outdir}/${sample_id}_sorted_marked.bam"
tempdir="${outdir}/tmp"

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

# run GATK4's HaplotypeCaller
gatk HaplotypeCaller -ERC GVCF \
-R ${reference} --sequence-dictionary ${ref_path}/${index}.dict \
-I ${markedsortedbamfile} \
-O ${outdir}/${sample_id}.g.vcf.gz \
--tmp-dir $tempdir
