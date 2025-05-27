#!/usr/bin/env bash
#SBATCH -J plink_LD
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --time=08:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH --array=1-44

# usuage: sbatch PWS_SNP_workflow/LD.sbatch

wd="/mnt2/fscratch/users/agr_3216_ifapa/alorenzo/simulations"
chrom_list="${wd}/chromosomes.txt"
# vcf="${wd}/biallelic_chromsonly.vcf.gz"
chrom=$(cat $chrom_list | sed -n ${SLURM_ARRAY_TASK_ID}p)

module load plink/1.90_beta 

plink --bfile ${wd}/data_filtered \
--chr-set 44 \
--allow-extra-chr \
--chr ${chrom} \
--double-id \
--r2 \
--ld-window-kb 100 \
--ld-window 99999 \
--ld-window-r2 0 \
--out ${wd}/plink_LD_${chrom}

