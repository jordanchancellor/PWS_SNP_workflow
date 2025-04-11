#!/usr/bin/env python

import sys
import subprocess
import os
from pybedtools import BedTool
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


# This script will generate a file of flanking sequences from a genome fasta file using an invput VCF file
# note: you must load python version  > 3.7 in order for script to execute properly
# Usage: python3 getSNPflankingregions.py <input_vcf> <input_fasta> <flankingregionupstreamlength> <flankingregiondownstreamlength> OR only supply one length that will be used for both up and downstream

# define functions

def parse_vcf(vcf_file):
    snps = []
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            if "," in alt:
                alt = alt.split(",")[0]  # Take only first ALT allele
            snps.append((chrom, pos, ref, alt))
    return snps

def make_bed(snps, flank_up, flank_down):
    bed_lines = []
    for chrom, pos, ref, alt in snps:
        start = pos - int(flank_up) - 1  # 0-based BED format
        end = pos + int(flank_down)
        name = f"{chrom}:{start}-{end}"
        bed_lines.append(f"{chrom}\t{start}\t{end}\t{name}")
    
    bed_str = "\n".join(bed_lines)
    print("BED FILE CONTENTS:")
    print(bed_str)  # Debug: Print the BED content to check
    return BedTool(bed_str, from_string=True)

def extract_sequences(bed, fasta, snps, flank_up, flank_down):
    bedfile = bed.fn
    tmp_fasta = bedfile + ".tmp.fasta"
    output_fasta = bedfile + "_snp_flank.fasta"

    print(f"Extracting sequences from FASTA file: {fasta}")
    try:
        # Extract sequences from the FASTA file using the BED regions
        bed.sequence(fi=fasta, s=False, name=False, fo=tmp_fasta)
    except Exception as e:
        print(f"[ERROR] Failed to extract sequences: {e}")
        return None

    # Debug: Check if tmp_fasta file is created and readable
    print(f"Temporary FASTA file created at: {tmp_fasta}")

    # Make sure we can read the sequences from the temporary FASTA file
    seqs = SeqIO.to_dict(SeqIO.parse(tmp_fasta, "fasta"))
    print(f"[INFO] {len(seqs)} sequences were extracted.")
    
    records = []

    for i, (chrom, pos, ref, alt) in enumerate(snps):
        key = f"{chrom}:{pos - int(flank_up) - 1}-{pos + int(flank_down)}"

        if key not in seqs:
            print(f"[WARNING] Sequence for {key} not found in extracted sequences.")
            continue

        seq = str(seqs[key].seq)
        snp_index = int(flank_up)
        if snp_index >= len(seq):
            print(f"[WARNING] SNP index {snp_index} out of range for {key}")
            continue

        # Modify sequence to include SNP information
        modified_seq = seq[:snp_index] + f"[{ref}/{alt}]" + seq[snp_index + 1:]
        record_id = f"{chrom}_{pos}_{ref}_{alt}"
        record = SeqRecord(Seq(modified_seq), id=record_id, description="")
        records.append(record)

    if records:
        SeqIO.write(records, output_fasta, "fasta")
        print(f"[INFO] Flanking sequences saved to: {output_fasta}")
        os.remove(tmp_fasta)
        return output_fasta
    else:
        print("[ERROR] No sequences were extracted.")
        return None
    

def main():
    if len(sys.argv) not in [4, 5]:
        print("Usage: python3 getSNPgenomeannotations.py <input_vcf> <input_fasta> <flankingregionupstreamlength> <flankingregiondownstreamlength>")
        print("OR")
        print("Usage: python3 getSNPgenomeannotations.py <input_vcf> <input_fasta> <flankingregionlength>")
        sys.exit(1)

    input_vcf = sys.argv[1]
    input_fasta = sys.argv[2]

    if len(sys.argv) == 4:
        flank_up = flank_down = sys.argv[3]
    else:
        flank_up = sys.argv[3]
        flank_down = sys.argv[4]

    if not input_vcf.endswith(".vcf"):
        print("Error: Input VCF file must have a .vcf extension.")
        sys.exit(1)

    if not input_fasta.endswith((".fa", ".fna")):
        print("Error: Input FASTA file must have .fa or .fna extension.")
        sys.exit(1)

    print("[INFO] Parsing VCF...")
    snps = parse_vcf(input_vcf)

    print("[INFO] Creating BED entries for flanking regions...")
    bed = make_bed(snps, flank_up, flank_down)

    print("[INFO] Extracting sequences using pybedtools...")
    output_fasta = extract_sequences(bed, input_fasta, snps, flank_up, flank_down)

    if output_fasta:
        print(f"[DONE] Flanking sequences saved to: {output_fasta}")
    else:
        print("[ERROR] No sequences were extracted.")

if __name__ == "__main__":
    main()