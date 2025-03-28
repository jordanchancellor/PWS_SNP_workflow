#!/usr/bin/env bash
#SBATCH -J mapping_stats
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# usuage: sbatch PWS_SNP_workflow/mappingstats.sbatch

wd="/mnt/scratch/users/agr_3216_ifapa/alorenzo/vannamei_snp_workflow/"
dir="${wd}/outfiles/"

# Define output CSV file and write the header
output_file="mapping_summary.csv"
echo "sample_id,concordant_zero,concordant_zero_percent,concordant_one,concordant_one_percent,concordant_greater,concordant_greater_percent,total_mapped" > "$output_file"

# Loop through each input file as a command-line argument
for file in ${dir}/bowtie_*.err; do
    # Extract the sample ID by removing the suffix ".fq.gz_trimming_report.txt"
    sample_id=$(basename "$file" | sed 's/bowtie_//g' | sed 's/\.err//g')

    concordant_zero=$(grep "aligned concordantly 0 times" "$file" | awk '{print $1}' | tr -d '')
    concordant_zero_percent=$(grep "aligned concordantly 0 times" "$file" | awk '{print $2}' | tr -d '' | sed 's/(*)//g' | sed 's/(//g')

    concordant_one=$(grep "aligned concordantly exactly" "$file" | awk '{print $1}' | tr -d '')
    concordant_one_percent=$(grep "aligned concordantly exactly" "$file" | awk '{print $2}' | tr -d '' | sed 's/(*)//g' | sed 's/(//g')

    concordant_greater=$(grep "aligned concordantly >1" "$file" | awk '{print $1}' | tr -d '')
    concordant_greater_percent=$(grep "aligned concordantly >1" "$file" | awk '{print $2}' | tr -d '' | sed 's/(*)//g' | sed 's/(//g')

    total_mapped=$(grep "overall alignment rate" "$file" | awk '{print $1}' | tr -d '')

    # Write the results to the output file
    echo "$sample_id,$concordant_zero,$concordant_zero_percent,$concordant_one,$concordant_one_percent,$concordant_greater,$concordant_greater_percent,$total_mapped" >> "$output_file"
done
