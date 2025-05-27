#!/usr/bin/env python

import sys
import subprocess
import os
from pybedtools import BedTool
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


# This script will extract snp effect annotations (i.e. snpeff) from input annotated vcf file and output a tab-delimited table
# note: you must load python version  > 3.7 in order for script to execute properly
# Usage: python3 getSNPeffectannotations.py <input_vcf>

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
            info = fields[7]
            snps.append((chrom, pos, ref, alt, info))
    return snps

def create_annotation_file(snps):
    output = []

    header = [
        "CHROM", "POS", "REF", "ALT",
        "Effect", "Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Transcript_ID",
        "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA_pos_len",
        "CDS_pos_len", "AA_pos_len", "Distance", "Errors"
    ]
    output.append("\t".join(header))

    for chrom, pos, ref, alt, info in snps:
        ann_field = ""
        for entry in info.split(";"):
            if entry.startswith("ANN="):
                ann_field = entry.replace("ANN=", "")
                break

        if ann_field:
            annotations = ann_field.split(",")
            for ann in annotations:
                ann_parts = ann.split("|")
                
                def safe_get(index):
                    return ann_parts[index] if len(ann_parts) > index else "NA"

                effect = safe_get(1)
                impact = safe_get(2)
                gene_name = safe_get(3)
                gene_id = safe_get(4)
                feature_type = safe_get(5)
                transcript_id = safe_get(6)
                biotype = safe_get(7)
                rank = safe_get(8)
                hgvs_c = safe_get(9)
                hgvs_p = safe_get(10)
                cdna_pos = safe_get(11)
                cds_pos = safe_get(12)
                aa_pos = safe_get(13)
                distance = safe_get(14)
                error = safe_get(15)

                line = [
                    chrom, str(pos), ref, alt,
                    effect, impact, gene_name, gene_id, feature_type, transcript_id,
                    biotype, rank, hgvs_c, hgvs_p, cdna_pos,
                    cds_pos, aa_pos, distance, error
                ]
                output.append("\t".join(line))
        else:
            # If no ANN field, still add one row with "NA" annotation info
            line = [chrom, str(pos), ref, alt] + ["NA"] * (len(header) - 4)
            output.append("\t".join(line))

    output_str = "\n".join(output)
    print("ANNOTATION FILE CONTENTS:")
    print(output_str)
    return output_str

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 getSNPeffectannotations.py <input_vcf>")
        sys.exit(1)

    input_vcf = sys.argv[1]

    if not input_vcf.endswith(".vcf"):
        print("Error: Input VCF file must have a .vcf extension.")
        sys.exit(1)

    print("[INFO] Parsing VCF...")
    snps = parse_vcf(input_vcf)

    print("[INFO] Extracting SNP effect annotation from vcf file...")
    annotations = create_annotation_file(snps)

    if annotations:
        with open("snp_annotations.tsv", "w") as f:
            f.write(annotations)
        print(f"[DONE] SNP annotations saved to: snp_annotations.tsv")
    else:
        print("[ERROR] No position effects were extracted.")

if __name__ == "__main__":
    main()
