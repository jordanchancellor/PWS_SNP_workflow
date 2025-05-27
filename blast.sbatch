#!/usr/bin/env bash
#SBATCH -J blast
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# usage: sbatch PWS_SNP_workflow/blast.sbatch <query> <db> <type>

# Define variables
query=$1
db=$2
type=$3
prefix="${db%.*}"
THREADS=$SLURM_CPUS_PER_TASK

# Load required modules
module purge
module load blast_plus/2.15.0+

# check that command input is correct
if [ "$#" -ne 3 ]; then
    echo "Error: input files not specified."
    echo "Usage: sbatch PWS_SNP_workflow/blast.sbatch <query> <db> <type>"
    echo "type options: nucl, prot"
    exit 1
fi

# check if db exists
if [[ "$type" == "nucl" ]]; then
    db_extensions=("nhr" "nin" "nsq")
elif [[ "$type" == "prot" ]]; then
    db_extensions=("phr" "pin" "psq")
else
    echo "Error: type must be 'nucl' or 'prot'"
    exit 1
fi

db_missing=false
for ext in "${db_extensions[@]}"; do
    if [[ ! -f "${prefix}.${ext}" ]]; then
        db_missing=true
        break
    fi
done

# Create blastdb if necessary
if $db_missing; then
    echo "Creating BLAST database from $db..."
    makeblastdb -in ${db} -dbtype ${type} -out ${prefix}
else
    echo "BLAST database already exists: ${prefix}"
fi

# Run BLAST
if [[ "$type" == "nucl" ]]; then
    echo "Running blastn..."
    blastn -query ${query} -db ${prefix} \
        -outfmt '6 qseqid sseqid sstart send qstart qend sstrand evalue bitscore' \
        -evalue 1e-5 \
        -num_threads ${THREADS} \
        -out blast_results.tsv
else
    echo "Running blastp..."
    blastp -query ${query} -db ${prefix} \
        -outfmt '6 qseqid sseqid sstart send qstart qend evalue bitscore' \
        -evalue 1e-5 \
        -num_threads ${THREADS} \
        -out blast_results.tsv
fi

# Convert results to .bed format
awk '{
  chrom=$2;
  start=($3 < $4 ? $3 : $4) - 1;
  end=($3 > $4 ? $3 : $4);
  print chrom"\t"start"\t"end
}' blast_results.tsv > blast_results.bed

# Extract best matches based on e-score
sort -k1,1 -k6,6g blast_results.tsv | awk '!seen[$1]++' > best_hits.tsv

# Convert top hit results to .bed format
awk '{
  chrom=$2;
  start=($3 < $4 ? $3 : $4) - 1;
  end=($3 > $4 ? $3 : $4);
  print chrom"\t"start"\t"end
}' best_hits.tsv > top_blast_matches.bed

