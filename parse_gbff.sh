#!/bin/bash

# check if file to be parsed is provided
if [ "$#" -ne 1 ]; then
    printf "No input file provides %s\n Usage: sbatch $0 <file>"
    exit 1
fi

gbff_file=$1

# check if the file exists
if [ ! -f "$gbff_file" ]; then
    echo "Error: The file $gbff_file does not exist."
    exit 1
fi

# check if the file is gzipped
if [[ "$gbff_file" == *.gz ]]; then
    is_gzipped="True"
else
    is_gzipped="False"
fi

# run biopython to parse through gbff file provided
python3 - <<END
import sys
import importlib.util

# Check if Bio package is installed
package_name = 'Bio'
spec = importlib.util.find_spec(package_name)
if spec is None:
    print(f"Error: {package_name} is not installed.")
    sys.exit(1)
else:
    print(f"{package_name} is installed. Version is: ")
    import Bio
    print(Bio.__version__)

from Bio import SeqIO
import gzip

gbff_file = "$gbff_file"
is_gzipped = "$is_gzipped" == "True"  # Convert the passed string to boolean

# parse file
try:
    if is_gzipped:
        # Open the gzipped file with gzip.open() in text mode ('rt')
        with gzip.open(gbff_file, "rt") as file_handle:
            for seq_record in SeqIO.parse(file_handle, "genbank"):
                print(f"ID: {seq_record.id}")
                print(f"Sequence: {repr(seq_record.seq)}")
                print(f"Length: {len(seq_record)}")
    else:
        # Open the regular (uncompressed) file
        with open(gbff_file, "r") as file_handle:
            for seq_record in SeqIO.parse(file_handle, "genbank"):
                print(f"ID: {seq_record.id}")
                print(f"Sequence: {repr(seq_record.seq)}")
                print(f"Length: {len(seq_record)}")
except Exception as e:
    print(f"Error parsing the file: {e}")
END
