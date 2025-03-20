import sys
import subprocess

# Usage: python3 computeheterozygosity.py <vcf_file> (can be gzipped)
# Note: must have bcftools and python version > 3.7 loaded into path prior to running

# Define input and output files
vcf_file = sys.argv[1]

if vcf_file.endswith(".vcf.gz"):
    output_file = vcf_file.rsplit(".vcf.gz", 1)[0] + "_heterozygosity.txt"
else:
    output_file = vcf_file.rsplit(".vcf", 1)[0] + "_heterozygosity.txt"

# Run bcftools query to extract genotypes
cmd = f"bcftools query -f '%CHROM\t%POS\t[ %GT]\n' {vcf_file}"
process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
results = []

# Process bcftools output
for line in process.stdout.strip().split("\n"):
    cols = line.split("\t")
    chrom, pos = cols[:2]
    genotypes = " ".join(cols[2:]).split()

    # Filter out missing genotypes
    valid_genotypes = [gt for gt in genotypes if gt != "./."]

    # Count heterozygous genotypes (0/1, 1/0, 0|1, 1|0)
    heterozygous_count = sum(gt in {"0/1", "1/0", "0|1", "1|0"} for gt in valid_genotypes)
    
    print(f"Site {chrom}:{pos} -> Valid Genotypes: {valid_genotypes}, Heterozygous Count: {heterozygous_count}")

    total_samples = len(valid_genotypes)
    heterozygosity = heterozygous_count / total_samples

    # Store result
    results.append(f"{chrom}\t{pos}\t{heterozygosity:.4f}")

# Write all results to file
with open(output_file, "w") as file:
    file.write("CHROM\tPOS\tHETEROZYGOSITY\n")  # Header
    file.write("\n".join(results) + "\n")

print(f"Heterozygosity results saved to: {output_file}")