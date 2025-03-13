#!/usr/bin/env python

import gzip
import io
import sys

# This file will generate a bedfile of the masked regions a fasta file.
# Adapted from Daniel E. Cook at https://www.danielecook.com/generate-a-bedfile-of-masked-ranges-a-fasta-file/
# Usage: python3 generate_masked_ranges.py <fasta file | .fa or .fa.gz> > masked_ranges.bed

# STDIN or arguments
if len(sys.argv) > 1:

    # Check file type
    if sys.argv[1].endswith(".fa.gz"):
        input_fasta = io.TextIOWrapper(io.BufferedReader(gzip.open(sys.argv[1])))
    elif sys.argv[1].endswith(".fa") or sys.argv[1].endswith(".txt") or sys.argv[1].endswith(".fasta"):
        input_fasta = open(sys.argv[1],'r')
    else:
        raise Exception("Unsupported File Type")
else:
    print("Usage: generate_masked_ranges.py <fasta file | .fa or .fa.gz>")

# initialize varables
n, state = 0, 0 # line, character, state (0=Out of gap; 1=In Gap)
chrom, start, end = None, None, None

# process input FASTA files
with input_fasta as f:
    for line in f:
        line = line.replace("\n","")
        if line.startswith(">"):
            # Print end range
            if state == 1:
                print('\t'.join([chrom ,str(start), str(n)]))
                start, end, state  = 0, 0, 0
            n = 0 # Reset character
            chrom = line.split(" ")[0].replace(">","")

        else:
            for char in line:
                if state == 0 and char == "N":
                    state = 1
                    start = n
                elif state == 1 and char != "N":
                    state = 0
                    end = n
                    print('\t'.join([chrom ,str(start), str(end)]))
                else:
                    pass

                n += 1 # First base is 0 in bed format.

# Print mask close if on the last chromosome.
if state == 1:
            print('\t'.join([chrom ,str(start), str(n)]))
            start, end, state  = 0, 0, 0