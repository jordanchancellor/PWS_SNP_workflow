#!/usr/bin/env bash
#SBATCH -J SNeP
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# Set input prefix
NeEstimator="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/simulations/"
wd="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/simulations/"
INPUT="${wd}/data_filtered"
THREADS=$SLURM_CPUS_PER_TASK

# Chromosome list
CHROMS=$(cut -f1 ${INPUT}.bim | sort | uniq)

# Load modules
module purge
module load bcftools/1.16
module load vcftools/20160205
module load plink/1.90_beta