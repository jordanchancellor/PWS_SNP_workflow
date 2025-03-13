#!/usr/bin/env bash
#SBATCH -J selfblast
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x_%a.err
#SBATCH --output=%x_%a.out
#SBATCH --array=1-44

# a module-dependent bash script to self blast each chromosome
# change array to number of chromosomes to run self-blast on 
# usage: sbatch PWS_SNP_workflow/selfblast.sbatch <list>

# check if aligner is specified
if [ "$#" -ne 1 ]; then
    echo "Error: files not specified."
    echo "Usage: sbatch PWS_SNP_workflow/selfblast.sbatch <list>"
    exit 1
fi

# load environment
module purge
module load blast_plus/2.15.0+

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
list=$1
chrom=$(cat $list | sed -n ${SLURM_ARRAY_TASK_ID}p)
path=$(dirname "$list")
fasta="${path}/${chrom}.fasta"
outdir="${path}/selfblast"

# create output directory (if needed)
if [[ ! -d $outdir ]]; then
	printf "Creating output directory %s\n" "${outdir}"
	mkdir $outdir 
fi

# self-blast to determine repetitive sequences
makeblastdb -in $fasta -dbtype nucl
blastn -query $fasta -db $fasta -outfmt 7 -evalue 1e-10 -out ${outdir}/${chrom}_self_blast.txt
