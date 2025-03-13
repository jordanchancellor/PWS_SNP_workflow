#!/usr/bin/env bash
#SBATCH -J selfblast
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=7-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# a module-dependent bash script to search for repeat regions in a specified reference genome
# usage: sbatch PWS_SNP_workflow/repeats.sbatch <reference> [optional] <bamfile>

# check if aligner is specified
if [ -z "$1" ]; then
    echo "Error: No genome specified."
    echo "Usage: PWS_SNP_workflow/repeats.sbatch <genome> [optional] <bamfile>"
    exit 1
fi

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
reference=$1
ref_path=$(dirname "$reference")
ref_base=$(basename "${reference}")
outdir="$ref_path/repeatregionsearch"  # Define the new directory path
bamfile=$2 # (optional)
bam_base=$(basename "${bamfile}" .bam)

# load environment
module purge
module load jellyfish/2.3.0
module load bbmap/38.92
module load blast_plus/2.15.0+
module load samtools/1.21

# create output directory (if needed)
if [[ ! -d $outdir ]]; then
	printf "Creating output directory %s\n" "${outdir}"
	mkdir $outdir 
fi

# Count k-mer Frequencies with Jellyfish
jellyfish count -m 21 -s 4G -t $THREADS -C $reference -o ${outdir}/${ref_base}_mer_counts.jf
jellyfish histo ${outdir}/${ref_base}_mer_counts.jf > ${outdir}/${ref_base}_kmer_histo.txt

# Compute depth at each position with samtools (if BAM file is provided)
if [ -f "$bamfile" ]; then
    # samtools depth -a "$bamfile" > "${outdir}/${bam_base}_depth.txt"

    # Get coverage distribution
    awk '{print $3}' "${outdir}/${bam_base}_depth.txt" | sort -n | uniq -c > "${outdir}/coverage_distribution.txt"

    # Compute average depth
    awk '{sum+=$3; count++} END {print "Avg Depth:", sum/count}' "${outdir}/${bam_base}_depth.txt"

    # Compute median depth
    awk '{print $3}' "${outdir}/${bam_base}_depth.txt" | sort -n | awk '{a[i++]=$1} END{print "Median Depth:", (i%2==1)? a[int(i/2)] : (a[int(i/2)-1] + a[int(i/2)])/2}'

    # Identify positions with coverage >2 SD above mean
    awk '{sum+=$3; sumsq+=$3*$3; count++} 
         END { mean=sum/count; stddev=sqrt(sumsq/count - mean*mean); 
               print "Mean Depth:", mean; 
               print "Standard Deviation:", stddev; 
               print "Threshold for repeats:", mean + 2*stddev }' "${outdir}/${bam_base}_depth.txt"

    # Extract potential repetitive regions
    threshold=$(awk '{sum+=$3; sumsq+=$3*$3; count++} 
                     END { mean=sum/count; stddev=sqrt(sumsq/count - mean*mean); 
                           print mean + 2*stddev }' "${outdir}/${bam_base}_depth.txt")

    awk -v threshold="$threshold" '$3 > threshold' "${outdir}/${bam_base}_depth.txt" > "${outdir}/${bam_base}_high_coverage.txt"
fi

# self-blast to determine repetitive sequences
makeblastdb -in $reference -dbtype nucl
blastn -query $reference -db $reference -outfmt 7 -perc_identity 90 -qcov_hsp_perc 80 -evalue 1e-10 -out ${outdir}/${ref_base}_self_blast.txt
