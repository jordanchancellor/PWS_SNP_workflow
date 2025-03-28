#!/usr/bin/env bash
#SBATCH -J trimming_stats
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# usuage: sbatch PWS_SNP_workflow/trimmingstats.sbatch

wd="/mnt/scratch/users/agr_3216_ifapa/alorenzo/vannamei_snp_workflow/"
dir="${wd}/outfiles/"

# Define output CSV file and write the header
output_file="trimming_summary.csv"
echo "sample_id,total_reads,both_surviving,both_surviving_percent,forward_only,forward_only_percent,reverse_only,reverse_only_percent" > "$output_file"

# Loop through each input file as a command-line argument
for file in ${dir}/trim_*.err; do
    # Extract the sample ID by removing the suffix ".fq.gz_trimming_report.txt"
    sample_id=$(basename "$file" | sed 's/trim_//g' | sed 's/\.err//g')

    total_reads=$(grep "Input Read Pairs:" "$file" | awk '{print $4}' | tr -d '')

    both_surviving=$(grep "Both Surviving:" "$file" | awk '{print $7}' | tr -d '')
    both_surviving_percent=$(grep "Both Surviving:" "$file" | awk '{print $8}' | tr -d ''| sed 's/(*)//g' | sed 's/(//g')
    
    forward_only=$(grep "Forward Only Surviving:" "$file" | awk '{print $12}' | tr -d '')
    forward_only_percent=$(grep "Forward Only Surviving:" "$file" | awk '{print $13}' | tr -d '' | sed 's/(*)//g' | sed 's/(//g')
    
    reverse_only=$(grep "Reverse Only Surviving:" "$file" | awk '{print $17}' | tr -d '')
    reverse_only_percent=$(grep "Reverse Only Surviving:" "$file" | awk '{print $18}' | tr -d '' | sed 's/(*)//g' | sed 's/(//g')

    # Write the results to the output file
    echo "$sample_id,$total_reads,$both_surviving,$both_surviving_percent,$forward_only,$forward_only_percent,$reverse_only,$reverse_only_percent" >> "$output_file"
done
