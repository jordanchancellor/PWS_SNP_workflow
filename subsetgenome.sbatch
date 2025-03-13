#!/usr/bin/env bash
#SBATCH -J splitgenome
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# a bash script to split a genome fasta file based on input list
# usage: sbatch PWS_SNP_workflow/subsetgenome.sbatch <genomefasta> <splitlist>

# check if aligner is specified
if [ "$#" -ne 2 ]; then
    echo "Error: files not specified."
    echo "Usage: sbatch PWS_SNP_workflow/buildindex.sbatch <genomefasta> <splitlist>"
    exit 1
fi

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
genome=$1
path=$(dirname "$genome")
list=$2

##### use this section to output a fasta file for each individual chromosome

# while read -r CHROM; do
#     output="${path}/${CHROM}.fasta"
#     awk -v chrom=">$CHROM" '
#         $0 ~ chrom { print_seq=1; print $0; next }
#         /^>/ { print_seq=0 }
#         print_seq
#     ' "$genome" > "$output"

#     echo "Extracted $CHROM into $output"
# done < "$list"

##### use this section to output a single fasta file containing only those in <splitlist>

genome_base=$(basename "${genome}" .fna)
output="${path}/${genome_base}_chroms.fasta"

module load seqkit/2.2.0 
seqkit grep -f $list $genome > $output

echo "Extraction complete. Output saved to '$output'."
