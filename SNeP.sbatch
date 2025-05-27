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
SNEP="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/simulations/SNeP_111"
wd="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/simulations/SNeP_wholegenome"
INPUT="${wd}/data_filtered"
THREADS=$SLURM_CPUS_PER_TASK

# Chromosome list
CHROMS=$(cut -f1 ${INPUT}.bim | sort | uniq)

# Load modules
module purge
module load bcftools/1.16
module load vcftools/20160205
module load plink/1.90_beta

# Run SNeP per-chromosome
for CHR in ${CHROMS}; do
    echo "Processing chromosome $CHR..."

    # Extract chromosome
    plink --bfile "${INPUT}" --double-id --allow-extra-chr --chr "${CHR}" --make-bed --out "${wd}/chr${CHR}"

    # Convert to PED/MAP
    plink --bfile "${wd}/chr${CHR}" --allow-extra-chr --recode --out "${wd}/chr${CHR}"

    # Run SNeP with custom parameters
    ${SNEP} \
        -ped ${wd}/chr${CHR}.ped \
        -map ${wd}/chr${CHR}.map \
        -maxdist 1000 \
        -mindist 0 \
        -binwidth 25 \
        -maf 0.01 \
        -threads ${THREADS} \
        -out ./chr${CHR}_Ne_output

    Run SNeP with default parameters
    ${SNEP} \
    -ped ${wd}/chr${CHR}.ped \
    -map ${wd}/chr${CHR}.map \
    -maf 0 \
    -threads ${THREADS} \
    -out ./chr${CHR}_Ne_output
done

wait
echo "All SNeP jobs complete."
