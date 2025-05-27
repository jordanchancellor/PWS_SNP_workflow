#!/usr/bin/env python

import sys
import subprocess
import os
import pandas as pd
import numpy as np

# This script will take provided input files, identify outliers and output sites
# If an annotation file is included (gtf), script will output annotations located within 100bp of outlier sites
# note: you must load python version  > 3.7 in order for script to execute properly
# note: input files must contain four columns: CHROM, START, END, VALUE
# Usage: python3 selectionoutliers.py <input_files> [gtf_file]

# define functions

def identify_outlier_positions(df):
    """Function to identify outlier positions from input files based on percentiles."""
    # calculate upper and lower limits
    Q1 = df['VALUE'].quantile(0.25)
    Q3 = df['VALUE'].quantile(0.75)
    IQR = Q3 - Q1
    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR

    # Select outlier rows
    outliers = df[(df['VALUE'] < lower) | (df['VALUE'] > upper)]
    return outliers

def identify_nearby_annotations(outliers, annotation_file):
    """Function to identify annotations located within 100bp of outlier sites based on input gtf genome annotation file."""
    # read  annotation file
    gtf_cols = ['CHROM', 'source', 'feature', 'START', 'END', 'score', 'strand', 'frame', 'attribute']
    gtf = pd.read_csv(annotation_file, sep='\t', comment='#', header=None, names=gtf_cols)

    nearby_annotations = []

    for _, outlier in outliers.iterrows():
        chrom = outlier['CHROM']
        start = outlier['START']
        end = outlier['END']

        # Filter annotations on the same chromosome and within ±100bp of the outlier region
        nearby = gtf[
            (gtf['CHROM'] == chrom) &
            (gtf['END'] >= start - 100) &
            (gtf['START'] <= end + 100)
        ]
        if not nearby.empty:
            nearby_annotations.append(nearby)

    # output annotations
    if nearby_annotations:
        return pd.concat(nearby_annotations).drop_duplicates()
    else:
        return pd.DataFrame()

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 selectionoutliers.py <input_file1.txt> [input_file2.txt ...] [annotations.gtf]")
        sys.exit(1)

    input_files = [f for f in sys.argv[1:] if f.endswith('.txt')]
    annotation_file = next((f for f in sys.argv[1:] if f.endswith('.gtf')), None)

    all_outliers = []

    for file in input_files:
        print(f"Processing {file}...")
        df = pd.read_csv(file, sep='\t', header=None, names=["CHROM", "START", "END", "VALUE"])
        outliers = identify_outlier_positions(df)
        all_outliers.append(outliers)

    if not all_outliers:
        print("No outliers found in any input file.")
        sys.exit(0)

    combined_outliers = pd.concat(all_outliers).drop_duplicates(subset=["CHROM", "START", "END"])
    combined_outliers = combined_outliers.sort_values(by=["CHROM", "START", "END"])
    
    combined_outliers_file = "all_outliers_combined.txt"
    combined_outliers.to_csv(combined_outliers_file, sep='\t', index=False)
    print(f"Combined outliers saved to {combined_outliers_file}")

    if annotation_file:
        annotations = identify_nearby_annotations(combined_outliers, annotation_file)
        if not annotations.empty:
            annotation_output_file = "annotations_near_combined_outliers.txt"
            annotations.to_csv(annotation_output_file, sep='\t', index=False)
            print(f"Nearby annotations saved to {annotation_output_file}")
        else:
            print("No nearby annotations found.")

######### if you have combined results with ROH values using bedtools intersect, comment out the above section and use below to extract annotations

    # annotation_file = next((f for f in sys.argv[1:] if f.endswith('.gtf')), None)
    # file = sys.argv[1]
    # combined_outliers = pd.read_csv(file, sep='\t', header=None, names=["CHROM", "START", "END"])
    
    # if annotation_file:
    #     annotations = identify_nearby_annotations(combined_outliers, annotation_file)
    #     if not annotations.empty:
    #         annotation_output_file = "roh_combined_results_annotations.txt"
    #         annotations.to_csv(annotation_output_file, sep='\t', index=False)
    #         print(f"Nearby annotations saved to {annotation_output_file}")
    #     else:
    #         print("No nearby annotations found.")

if __name__ == "__main__":
    main()