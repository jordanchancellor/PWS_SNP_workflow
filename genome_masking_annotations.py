#!/usr/bin/env python

import gzip
import io
import subprocess
import sys
import pandas as pd
import pybedtools
import subprocess

# This script will identify coding and non-coding regions using an input gff file and output in bed format
# if run with option --overlap, will output intersections of coding and non-coding regions with input masked genome .bed file
# note: you must load python version  > 3.7 in order for script to execute properly
# note: you must also have bedtools installed in path for script to execute properly
# Usage: python3 script.py <gff file> [--overlap <masked bed file>]

def parse_gff(gff_file):
    """Parses the GFF file and returns a DataFrame with genomic features."""
    
    if gff_file.endswith(".gz"):
        gff_file = io.TextIOWrapper(io.BufferedReader(gzip.open(gff_file, "rb")))
    elif gff_file.endswith(".gff"):
        gff_file = open(gff_file, "r")
    else:
        print("Invalid GFF file format. Expected <.gff or .gff.gz file>")
        sys.exit(1)
    
    gff_df = pd.read_csv(gff_file, comment='#', sep="\t", header=None, 
                         names=["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"])
    return gff_df


def generate_genome_bed(gff_file):
    """Generates a genome-wide BED file from the GFF file."""
    gff_df = parse_gff(gff_file)
    
    # Get the unique chromosomes and their genomic ranges
    genome_bed_file = gff_file.rsplit(".", 1)[0] + "_genome.bed"
    genome_bed = gff_df.groupby("chrom").agg({"start": "min", "end": "max"}).reset_index()

    # Save as a BED file
    genome_bed[["chrom", "start", "end"]].to_csv(genome_bed_file, sep="\t", index=False, header=False)
    print(f"Genome-wide BED file saved to {genome_bed_file}")
    
    return genome_bed_file


def generate_bed_files(gff_file, genome_bed_file):
    """Generates three BED files: one for all features, one for coding regions, and one for noncoding regions."""
    gff_df = parse_gff(gff_file)
    
    # 1. BED file for all genomic features
    all_bed_file = gff_file.rsplit(".", 1)[0] + "_all.bed"
    gff_df[["chrom", "start", "end"]].to_csv(all_bed_file, sep="\t", index=False, header=False)
    print(f"All regions saved to {all_bed_file}")

    # 2. BED file for coding regions (genes and exons)
    coding_bed_file = gff_file.rsplit(".", 1)[0] + "_coding.bed"
    coding_df = gff_df[gff_df["feature"].isin(["gene", "exon"])]
    coding_df[["chrom", "start", "end"]].to_csv(coding_bed_file, sep="\t", index=False, header=False)
    print(f"Coding regions saved to {coding_bed_file}")
    
    # 3. Generate noncoding regions by subtracting coding regions from the genome
    genome_bed = pybedtools.BedTool(genome_bed_file)
    coding_bed = pybedtools.BedTool(coding_bed_file)
    
    noncoding_bed = genome_bed.subtract(coding_bed)
    noncoding_bed_file = gff_file.rsplit(".", 1)[0] + "_noncoding.bed"
    noncoding_bed.saveas(noncoding_bed_file)
    print(f"Noncoding regions saved to {noncoding_bed_file}")

    return all_bed_file, coding_bed_file, noncoding_bed_file


def find_overlap(masked_bed_file, bed_file):
    """Finds the overlap between a masked BED file and another BED file (coding or noncoding)."""
    masked_bed = pybedtools.BedTool(masked_bed_file)
    bed_bed = pybedtools.BedTool(bed_file)
    
    overlap = masked_bed.intersect(bed_bed, wa=True, wb=True)
    return overlap


def save_overlap_results(overlap, overlap_file):
    """Saves the overlap results to a file."""
    overlap.saveas(overlap_file)
    print(f"Overlap results saved to {overlap_file}")


def main():
    """Main function to handle input arguments and workflow."""
    if len(sys.argv) < 2:
        print("Usage: python3 script.py <gff file> [--overlap <masked bed file>]")
        sys.exit(1)

    gff_file = sys.argv[1]
    
    # Generate a genome-wide BED file from the GFF file
    genome_bed_file = generate_genome_bed(gff_file)

    # Generate the BED files (all, coding, noncoding)
    all_bed_file, coding_bed_file, noncoding_bed_file = generate_bed_files(gff_file, genome_bed_file)

    # Check if the --overlap flag is provided
    if "--overlap" in sys.argv:
        masked_bed_file = sys.argv[sys.argv.index("--overlap") + 1]

        # Find overlap for coding and noncoding regions
        print("Finding overlap with masked regions...")

        # Overlap with coding regions
        coding_overlap = find_overlap(masked_bed_file, coding_bed_file)
        save_overlap_results(coding_overlap, "masked_in_coding_overlap.bed")

        # Overlap with noncoding regions
        noncoding_overlap = find_overlap(masked_bed_file, noncoding_bed_file)
        save_overlap_results(noncoding_overlap, "masked_in_noncoding_overlap.bed")


if __name__ == "__main__":
    main()
