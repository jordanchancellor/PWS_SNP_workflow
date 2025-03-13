#!/usr/bin/env python

import gzip
import io
import subprocess
import sys
import pandas as pd
import pybedtools

# This file will find masked intersections when provided with a .bed file of masked regions and a .gff file
# if a masked .fasta file and .gff file are provided, the script will first generate a .bed file of masked ranges using generate_masked_ranges.py, then find intersections
# if only a masked .fasta file is provided, the script will only generate a .bed file of masked ranges using generate_masked_ranges.py
# note: you must load python version  > 3.7 in order for script to execute properly
# note: you must also have bedtools installed in path for script to execute properly
# Usage: python3 find_masked_intersections.py <fasta file OR bed file> [gff file]


# define functions
def run_generate_masked_ranges(fasta_file):
    """Runs generate_masked_ranges.py and returns the output BED file name."""
    output_bed = fasta_file.rsplit(".", 1)[0] + "_masked.bed"
    print(f"Running generate_masked_ranges.py on {fasta_file}...")

    result = subprocess.run(["python3", "PWS_SNP_workflow/generate_masked_ranges.py", fasta_file], capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Error running generate_masked_ranges.py: {result.stderr}")
        sys.exit(1)

    # Save output to a BED file
    with open(output_bed, "w") as f:
        f.write(result.stdout)

    print(f"Masked regions saved to {output_bed}")
    return output_bed

def find_intersections(bed_file, gff_file):
    """Finds intersections between masked regions in BED and annotated regions in GFF."""
    print(f"Finding intersections between {bed_file} and {gff_file}...")

    if gff_file.endswith(".gz"):
        gff_file = io.TextIOWrapper(io.BufferedReader(gzip.open(gff_file, "rb")))
    elif gff_file.endswith(".gff"):
        gff_file = open(gff_file, "r")
    else:
        print("Invalid GFF file format. Expected <.gff or .gff.gz file>")
        sys.exit(1)

    # Load GFF and convert to BED
    gff_df = pd.read_csv(gff_file, comment='#', sep="\t", header=None, 
                         names=["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"])
    
    gff_df = gff_df[gff_df["feature"] == "gene"]  # Keep only gene annotations
    gff_bed_file = gff_file.name.rsplit(".", 1)[0] + ".bed"
    
    gff_df[["chrom", "start", "end", "attributes"]].to_csv(gff_bed_file, sep="\t", index=False, header=False)
    
    # Run intersection using pybedtools
    masked_bed = pybedtools.BedTool(bed_file)
    annotations_bed = pybedtools.BedTool(gff_bed_file)
    
    overlap = masked_bed.intersect(annotations_bed, wa=True, wb=True)
    
    # Save results
    overlap_file = "intersections.txt"
    overlap.saveas(overlap_file)
    
    print(f"Intersections saved to {overlap_file}")
    return overlap_file

def main():
    """Main function to handle input arguments and workflow."""
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python find_masked_intersections.py <fasta file OR bed file> [gff file]")
        sys.exit(1)

    input_file = sys.argv[1]
    
    # If only one file is given, assume it's a FASTA and generate a BED file
    if len(sys.argv) == 2:
        if input_file.endswith((".fa", ".fasta", ".fa.gz", ".txt")):
            run_generate_masked_ranges(input_file)
        else:
            print("Invalid input. Please provide a FASTA file (.fa, .fasta, .fa.gz) or a BED file.")
            sys.exit(1)

    # If two files are given, assume first is BED or FASTA, second is GFF
    elif len(sys.argv) == 3:
        gff_file = sys.argv[2]

        if input_file.endswith((".fa", ".fasta", ".fa.gz", ".txt")):
            bed_file = run_generate_masked_ranges(input_file)
            find_intersections(bed_file, gff_file)
        elif input_file.endswith(".bed"):
            find_intersections(input_file, gff_file)
        else:
            print("Invalid input. Expected <fasta OR bed file> <gff file>")
            sys.exit(1)

if __name__ == "__main__":
    main()
