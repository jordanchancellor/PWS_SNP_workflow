#!/usr/bin/env bash
#SBATCH -J duplicate_stats
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# usuage: sbatch PWS_SNP_workflow/duplicate_stats.sbatch

wd="/mnt/scratch/users/agr_3216_ifapa/alorenzo/vannamei_snp_workflow/"
dir="${wd}/outfiles/"

# Define output CSV file and write the header
output_file="markduplicates_summary.csv"
echo "sample_id,READ_PAIR_DUPLICATES,PERCENT_DUPLICATION" > "$output_file"

# Loop through each input file as a command-line argument
for file in ${dir}/*_markdupes.metrics; do
    # Extract the sample ID by removing the suffix ".fq.gz_trimming_report.txt"
    sample_id=$(basename "$file" | sed 's/_markdupes.metrics//g')

    READ_PAIR_DUPLICATES=$(grep -A 1 "READ_PAIR_DUPLICATES" "$file" | awk '{print $7}' | awk 'NR==2')
    PERCENT_DUPLICATION=$(grep -A 1 "PERCENT_DUPLICATION" "$file" | awk '{print $9}' | awk 'NR==2')

    # Write the results to the output file
    echo "$sample_id,$READ_PAIR_DUPLICATES,$PERCENT_DUPLICATION" >> "$output_file"
done
