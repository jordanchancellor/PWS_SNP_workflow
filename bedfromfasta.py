import sys
from Bio import SeqIO

def create_bed(fasta, output_bed):
    with open(output_bed, 'w') as bedfile:
        for record in SeqIO.parse(fasta, "fasta"):
            chrom = record.id  # Extract chromosome/contig name
            length = len(record.seq)  # Get sequence length
            bedfile.write(f"{chrom}\t0\t{length}\n")  # Write BED format

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py genome.fasta output.bed")
        sys.exit(1)

    fasta_file = sys.argv[1]
    bed_file = sys.argv[2]

    create_bed(fasta_file, bed_file)
    print(f"BED file saved as {bed_file}")
