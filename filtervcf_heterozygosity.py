import sys
import subprocess
import os
import gzip

# Usage: python3 filtervcf_heterozygosity.py <vcf_file> (can be gzipped) <heterozygosity_threshold> (numeric)
# Note: must have bcftools and python version > 3.7 loaded into path prior to running

# Heterozygosity threshold (adjust as needed)
heterozygosity_threshold = float(sys.argv[2])

# Define input and output VCF files
vcf_file = sys.argv[1]

if vcf_file.endswith(".vcf.gz"):
    output_vcf_file = vcf_file.rsplit(".vcf.gz", 1)[0] + "_filtered_heterozygosity_" + str(heterozygosity_threshold) + ".vcf"
else:
    output_vcf_file = vcf_file.rsplit(".vcf", 1)[0] + "_filtered_heterozygosity_" + str(heterozygosity_threshold) + ".vcf"

# Run bcftools query to extract genotypes
cmd = f"bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%FORMAT\t[ %GT]\n' {vcf_file}"
process = subprocess.run(cmd, shell=True, capture_output=True, text=True)

# Prepare to store filtered VCF lines
filtered_lines = []

# Process bcftools output
for line in process.stdout.strip().split("\n"):
    cols = line.split("\t")
    chrom, pos, id_field, ref, alt, qual, filter_field, info, format_field = cols[:9]
    genotypes = " ".join(cols[9:]).split()
    
    # Filter out missing genotypes (./.)
    valid_genotypes = [gt for gt in genotypes if gt != "./."]
    
    # Count heterozygous genotypes (0/1, 1/0, 0|1, 1|0)
    heterozygous_count = sum(gt in {"0/1", "1/0", "0|1", "1|0"} for gt in valid_genotypes)

    # Calculate heterozygosity per site
    total_samples = len(valid_genotypes)
    heterozygosity = heterozygous_count / total_samples

    # Apply heterozygosity threshold filtering
    if heterozygosity <= heterozygosity_threshold:
        # Keep the line
        genotypes_str = "\t".join(genotypes)
        filtered_lines.append(f"{chrom}\t{pos}\t{id_field}\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info}\t{format_field}\t{genotypes_str}")

print(f"Number of sites passing the threshold: {len(filtered_lines)}")

# Output the filtered VCF with header (this assumes the original VCF file has a header)
# First, get the header from the original VCF file
cmd_header = f"bcftools view -h {vcf_file}"
header_process = subprocess.run(cmd_header, shell=True, capture_output=True, text=True)

# Create the new VCF file
with open(output_vcf_file, "w") as output_file:
    # Write header to output file
    output_file.write(header_process.stdout)
    
    # Write the filtered VCF lines to the output file
    for line in filtered_lines:
        output_file.write(line + "\n")

print(f"Filtered VCF file saved to: {output_vcf_file}")
