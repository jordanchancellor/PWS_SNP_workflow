#!/usr/bin/env bash
#SBATCH -J splitintervals
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# a module-dependent bash script to split reference genome into intervals for joint genotyping using GATK4's SplitIntervals
# usage: sbatch PWS_SNP_workflow/splitintervals.sbatch

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
scatter=50
wd="/mnt/scratch/users/agr_3216_ifapa/alorenzo/vannamei_snp_workflow"
samples_file="${wd}/samples_file.txt"
sample_id=$(cat $samples_file | sed -n ${SLURM_ARRAY_TASK_ID}p)
reference="${wd}/indices/pvannamei/GCF_042767895.1_ASM4276789v1_genomic.fna" # Genome assembly ASM4276789v1
ref_path=$(dirname "$reference")
index=$(basename "${reference}" .fna)
outdir="${ref_path}/splitintervals"
# interval_list="${wd}/indices/pvannamei/GCF_042767895.1_ASM4276789v1_chromlist.list"

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

# run GATK4's SplitIntervals
gatk SplitIntervals \
-R ${reference} --sequence-dictionary ${ref_path}/${index}.dict \
-O ${outdir} --scatter-count $scatter \
--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION
