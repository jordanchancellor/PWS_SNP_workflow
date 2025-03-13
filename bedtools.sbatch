#!/usr/bin/env bash
#SBATCH -J convertmasked
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=1-00:00:00
#SBATCH --constraint=cal
#SBATCH --error=%x.err
#SBATCH --output=%x.out

# a module-dependent bash script extract masked regions in a fasta file and high coverage regions and convert to .bed format
# compares two resulting .bed file intersection regions and outputs overlap files
# usage: sbatch PWS_SNP_workflow/bedtools.sbatch <fasta> <highcoveragefile>

# define global variables
THREADS=${SLURM_CPUS_PER_TASK:-12}
maskedfasta=$1
ref_path=$(dirname "$maskedfasta")
ref_base=$(basename "${maskedfasta}")
highcoveragefile=$2 
highcoveragefile_path=$(dirname "$highcoveragefile")
highcoveragefile_base=$(basename "${highcoveragefile}" .txt)

# load environment
module purge
module load bedtools/2.31.0_gcc12

# convert masked fasta to .bed format
echo "Converting masked fasta ${maskedfasta} to .bed format"
bedtools maskfasta -fi ${makedfasta} -bed ${ref_path}/${ref_base}_masked_regions.bed -fo ${ref_path}/tmp.fasta

# convert High-Coverage Regions in .txt format to .bed format
echo "Converting High-Coverage Regions in ${highcoveragefile} to .bed format"
awk '
  BEGIN { OFS="\t" }
  NR==1 { start=$2; prev=$2; chr=$1; coverage=$3; next }
  {
    if ($1 == chr && $2 == prev + 1) {
      prev = $2  # Extend the interval
    } else {
      print chr, start, prev+1  # Output previous interval
      start = $2  # Start new interval
    }
    prev = $2
  }
  END { print chr, start, prev+1 }  # Output last interval
' "${highcoveragefile}" > "${highcoveragefile_path}/${highcoveragefile_base}_high_coverage.bed"

# find masked regions that overlap high-coverage regions
echo "Finding masked regions that overlap high-coverage regions."
bedtools intersect -a ${ref_path}/${ref_base}_masked_regions.bed -b "${highcoveragefile_path}/${highcoveragefile_base}_high_coverage.bed" > overlap_masked_highcov.bed

# find high-coverage regions which were NOT masked
echo "Finding high-coverage regions that were not masked."
bedtools intersect -v -a "${highcoveragefile_path}/${highcoveragefile_base}_high_coverage.bed" -b ${ref_path}/${ref_base}_masked_regions.bed > highcov_not_masked.bed
