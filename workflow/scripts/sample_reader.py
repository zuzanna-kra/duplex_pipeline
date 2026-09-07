#!/usr/bin/env python3

"""
sample_utils.py

This script provides three utility functions for sequencing data preprocessing:
1. Extract the sample name from a CSV given a FASTQ file.
2. Extract sample name with UDI (barcode/index).
3. Extract library name from a FASTQ filename.

Usage:
    python sample_utils.py --mode sample --samples_csv_file samples.csv --fastq read1.fq
    python sample_utils.py --mode sample_UDI --samples_csv_file samples.csv --fastq read1.fq
    python sample_utils.py --mode library --fastq read1.fq
"""

import csv
import argparse
import sys
import os

# Function to extract the UDI portion from the FASTQ file name 
def extract_udi(fastq_file):
    # Assumes filename format like: SLX12345.UDI-A01.sample.fq
    return os.path.basename(fastq_file).split('.')[1]

# Function to get sample name from CSV using the extracted UDI 
def extract_sample_name(csv_file, fastq_file):
    udi = extract_udi(fastq_file)

    # Open the csv file
    with open(csv_file) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # Assumes CSV format: sample_name, UDI
            if row[1] == udi:
                return row[0]
        
    # If UDI was not found, exit with an error
    sys.exit(f"Error: UDI '{udi}' not found in {csv_file}")

# Function to get sample name and append the UDI 
def extract_sample_name_with_UDI(csv_file, fastq_file):
    sample = extract_sample_name(csv_file, fastq_file)
    udi = extract_udi(fastq_file)
    return f"{sample}_{udi}"

# Extract the library name from the file name 
def extract_library_name(fastq_file):
    return os.path.basename(fastq_file).split('.')[0]

# 
def main():
    # Set up the command-line argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples_csv_file", help="CSV file mapping sample name to UDI")
    parser.add_argument("--fastq", help="FASTQ file")
    parser.add_argument("--mode", choices=["sample", "sample_UDI", "library"], required=True, help="What to extract")

    # Parse arguments from the command line 
    args = parser.parse_args()

    # Dispatch based on mode
    if args.mode == "sample":
        print(extract_sample_name(args.samples_csv_file, args.fastq))
    elif args.mode == "sample_UDI":
        print(extract_sample_name_with_UDI(args.samples_csv_file, args.fastq))
    elif args.mode == "library":
        print(extract_library_name(args.fastq))

# Run main if this script is executed directly 
if __name__ == "__main__":
    main()
