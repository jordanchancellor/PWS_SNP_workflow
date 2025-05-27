#!/usr/bin/env bash
#SBATCH -J gzip
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# usage: sbatch PWS_SNP_workflow/gzip.sbatch <file>

# check that command input is correct
if [ "$#" -ne 1 ]; then
    echo "Error: input file not specified."
    echo "sbatch PWS_SNP_workflow/gzip.sbatch <file>"
    exit 1
fi

# Define variables
file=$1

# zip input file 
gzip ${file}