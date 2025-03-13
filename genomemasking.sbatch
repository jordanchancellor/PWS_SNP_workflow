#!/usr/bin/env bash
#SBATCH -J masking
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=4-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# a module-dependent bash script to mask repetitive regions in input genome
# usage: sbatch PWS_SNP_workflow/genomemasking.sbatch <reference genome>

# check if aligner is specified
if [ -z "$1" ]; then
    echo "Error: No genome specified."
    echo "Usage: PWS_SNP_workflow/genomemasking.sbatch <reference genome>"
    exit 1
fi

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
reference=$1
ref_path=$(dirname "$reference")
ref_base=$(basename "${reference}")
outdir="${ref_path}/${ref_base}_masked" 

# Load  environment
module purge
module load samtools/1.21
# module load repeatmasker/4.1.5
module load bbmap/38.92


# create output directory (if needed)
if [[ ! -d $outdir ]]; then
	printf "Creating output directory %s\n" "${outdir}"
	mkdir $outdir 
fi

# run repeat masker on input genome
# RepeatMasker -pa 4 -species "your_species" -dir ${outdir}/${ref_base}_masked.fasta $reference

# run bbmask on input genome
echo "Masking input genome $reference with bbmask.sh"
bbmask.sh in=$reference out=${outdir}/${ref_base}_masked.fasta kr=10 threads=$THREADS
