#!/usr/bin/env python

import sys
import pybedtools
import pandas as pd

# This script will calculate the percentage of genome masked in both coding and noncoding regions
# requires that .bed files have been created for coding and noncoding regions using genome_masking_annotations.py
# note: you must load python version  > 3.7 in order for script to execute properly
# note: you must also have bedtools installed in path for script to execute properly
# Usage: python3 script.py <masked regions .bed file> <coding regions .bed file> <noncoding regions .bed file>


def calculate_percentage_coverage(region_bed_file, overlap_bed_file):
    """Calculates the percentage of the region covered by masked regions."""
    region_bed = pybedtools.BedTool(region_bed_file)
    overlap_bed = pybedtools.BedTool(overlap_bed_file)
    
    # Get the total length of the regions (coding or noncoding)
    total_length = sum([interval.end - interval.start for interval in region_bed])

    # Get the total length of the overlap
    overlap_length = sum([interval.end - interval.start for interval in overlap_bed])
    
    # Calculate the percentage of the region that is masked
    coverage_percentage = (overlap_length / total_length) * 100
    return coverage_percentage, total_length


def write_to_file(filename, content):
    """Writes the calculated values to a text file."""
    with open(filename, "w") as file:
        file.write(content)


def calculate_genome_coverage(masked_bed_file, coding_bed_file, noncoding_bed_file, output_file):
    """Calculate total genome coverage in coding and noncoding regions."""
    
    # Calculate overlap for coding and noncoding regions
    coding_overlap = pybedtools.BedTool(masked_bed_file).intersect(pybedtools.BedTool(coding_bed_file), wa=True, wb=True)
    noncoding_overlap = pybedtools.BedTool(masked_bed_file).intersect(pybedtools.BedTool(noncoding_bed_file), wa=True, wb=True)
    
    # Calculate percentage coverage in coding regions
    coding_coverage, coding_total = calculate_percentage_coverage(coding_bed_file, coding_overlap.fn)
    coding_result = f"Percentage of coding genome covered by masked regions: {coding_coverage:.2f}%\n"
    coding_result += f"Total coding genome length: {coding_total / 1e6:.2f} Mb\n"
    
    # Calculate percentage coverage in noncoding regions
    noncoding_coverage, noncoding_total = calculate_percentage_coverage(noncoding_bed_file, noncoding_overlap.fn)
    noncoding_result = f"Percentage of noncoding genome covered by masked regions: {noncoding_coverage:.2f}%\n"
    noncoding_result += f"Total noncoding genome length: {noncoding_total / 1e6:.2f} Mb\n"
    
    # Calculate total genome length (coding + noncoding)
    total_genome_length = coding_total + noncoding_total
    total_genome_result = f"Total genome length: {total_genome_length / 1e6:.2f} Mb\n"
    
    # Calculate the total overlap length (coding + noncoding)
    total_overlap_length = sum([interval.end - interval.start for interval in coding_overlap]) + \
                           sum([interval.end - interval.start for interval in noncoding_overlap])
    
    # Calculate the total percentage of genome masked in coding and noncoding
    total_genome_coverage = (total_overlap_length / total_genome_length) * 100
    total_genome_result += f"Total genome coverage by masked regions: {total_genome_coverage:.2f}%\n"
    
    # Calculate the percentage of the genome that is both masked and coding
    coding_masked_genome_percentage = (sum([interval.end - interval.start for interval in coding_overlap]) / total_genome_length) * 100
    coding_masked_result = f"Percentage of total genome covered by masked coding regions: {coding_masked_genome_percentage:.2f}%\n"
    
    # Calculate the percentage of the genome that is both masked and noncoding
    noncoding_masked_genome_percentage = (sum([interval.end - interval.start for interval in noncoding_overlap]) / total_genome_length) * 100
    noncoding_masked_result = f"Percentage of total genome covered by masked noncoding regions: {noncoding_masked_genome_percentage:.2f}%\n"
    
    # Calculate the total masked regions (coding + noncoding) as a percentage of the genome
    total_masked_region_percentage = (total_overlap_length / total_genome_length) * 100
    total_masked_result = f"Total percentage of genome masked (coding + noncoding): {total_masked_region_percentage:.2f}%\n"
    
    # Combine all the results
    results = coding_result + noncoding_result + total_genome_result + \
             coding_masked_result + noncoding_masked_result + total_masked_result
    
    # Write results to file
    write_to_file(output_file, results)


def main():
    """Main function to handle input arguments and calculate genome coverage."""
    if len(sys.argv) < 5:
        print("Usage: python3 calculate_genome_coverage.py <masked bed file> <coding bed file> <noncoding bed file> <output file>")
        sys.exit(1)

    masked_bed_file = sys.argv[1]
    coding_bed_file = sys.argv[2]
    noncoding_bed_file = sys.argv[3]
    output_file = sys.argv[4]

    # Calculate coverage for coding and noncoding regions and write to output file
    calculate_genome_coverage(masked_bed_file, coding_bed_file, noncoding_bed_file, output_file)


if __name__ == "__main__":
    main()
