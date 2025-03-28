import sys
import gzip
import re

# Usage: python3 snp_proximity_filtering.py <input.vcf.gz> <snp_distance (bp)> <output.vcf>
# Note: must have python version > 3.7 loaded into path prior to running
# a thank you to chatgpt for help with this one

def parse_vcf(input_vcf, bp_distance, output_vcf):
    def get_qd(info_field):
        match = re.search(r'QD=([0-9\\.]+)', info_field)
        return float(match.group(1)) if match else -1  # Default to -1 if QD is missing
    
    snps = []  # Store SNPs within the range
    prev_pos = -bp_distance  # Initialize with a large gap
    prev_chr = None

    input_open = gzip.open if input_vcf.endswith('.gz') else open
    with input_open(input_vcf, 'rt') as vcf_in, open(output_vcf, 'w') as vcf_out:
        for line in vcf_in:
            if line.startswith('#'):
                vcf_out.write(line)
                continue

            fields = line.strip().split('\t')
            chrom, pos = fields[0], int(fields[1])
            info = fields[7]
            full_line = fields

            # If new chromosome or position is outside the range, process SNPs
            if chrom != prev_chr or pos - prev_pos > bp_distance:
                if len(snps) > 1:
                    best_snp = max(snps, key=lambda x: get_qd(x[7]))  # Select best SNP
                    vcf_out.write('\t'.join(best_snp) + '\n')
                elif snps:
                    vcf_out.write('\t'.join(snps[0]) + '\n')
                snps = []

            snps.append(full_line)
            prev_chr, prev_pos = chrom, pos

        # Process remaining SNPs
        if len(snps) > 1:
            best_snp = max(snps, key=lambda x: get_qd(x[7]))
            vcf_out.write('\t'.join(best_snp) + '\n')
        elif snps:
            vcf_out.write('\t'.join(snps[0]) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 snp_proximity_filtering.py <input.vcf.gz> <snp_distance (bp)> <output.vcf>")
        sys.exit(1)

    input_vcf = sys.argv[1]
    bp_distance = int(sys.argv[2])
    output_vcf = sys.argv[3]

    parse_vcf(input_vcf, bp_distance, output_vcf)
