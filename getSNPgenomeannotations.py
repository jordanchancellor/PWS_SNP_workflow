#!/usr/bin/env python

import io
import subprocess
import sys
import pandas as pd
import pybedtools

# This script will generate a file of SNP annotations from a genome GFF file using an invput VCF file
# note: you must load python version  > 3.7 in order for script to execute properly
# note: you must also have bcftools and bedtools installed and loaded in path for script to execute properly
# Usage: python3 getSNPgenomeannotations.py <input_vcf> <input_gff>

# define functions
def get_snp_positions(vcf_file):
    """Extracts SNP positions from a VCF file and creates a BED file."""

    snp_bed_file = vcf_file.rsplit(".", 1)[0] + "_snp_positions.bed"

    command = f"bcftools view -H {vcf_file} | awk '{{print $1 \"\t\" $2 \"\t\" $2}}' > {snp_bed_file}"
    subprocess.run(command, shell=True, check=True)

    return snp_bed_file

def get_snp_annotations(snp_bed_file, input_gff):
    """Finds SNP overlaps with genome annotations."""

    # Convert GFF to BED format
    gff_bed_file = input_gff.rsplit(".", 1)[0] + ".bed"

    gff_df = pd.read_csv(input_gff, comment='#', sep="\t", header=None,
                         names=["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"])

    gff_df[["chrom", "start", "end", "feature", "attributes"]].to_csv(gff_bed_file, sep="\t", index=False, header=False)

    # Perform BEDTools intersection
    gff_bed = pybedtools.BedTool(gff_bed_file)
    snp_bed = pybedtools.BedTool(snp_bed_file)

    annotated_snps = snp_bed.intersect(gff_bed, wa=True, wb=True)

    output_annotations = snp_bed_file.rsplit(".", 1)[0] + "_annotations.txt"
    with open(output_annotations, "w") as out_file:
        out_file.write(str(annotated_snps))

    print(f"SNP position annotations saved to {output_annotations}")
    return output_annotations

def main():
    """Main function to handle input arguments and workflow."""
    if len(sys.argv) != 3:
        print("Usage: python3 getSNPgenomeannotations.py <input_vcf> <input_gff>")
        sys.exit(1)

    if len(sys.argv) == 3:
        
        input_vcf = sys.argv[1]
        input_gff = sys.argv[2]
        
        if not input_vcf.endswith((".vcf")):
                print("Invalid input. Input variant file must have .vcf extension. Trying unzipping your vcf file.")
        
        if not input_gff.endswith((".gff")):
                print("Invalid input. Input genome annotation file must have .gff extension")
        
        # Generate bed file of SNP positions
        snp_bed_file = get_snp_positions(input_vcf)

        # Annotate SNP positions using generated bed file
        get_snp_annotations(snp_bed_file, input_gff)

if __name__ == "__main__":
    main()