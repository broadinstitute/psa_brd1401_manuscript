#!/bin/bash

# This script is to be run for a single plat file (i.e. ./PO2Q01.fastq)

### Variable setup ###

input_path="$1" #full path to the file, including any directories
input_base=$(basename "$input_path" .fastq) #filename without extension
input_dir=$(dirname "$input_path") #the directory which contains the file

#get optional results_dir from $2, with default
step01_results_dir="$2" #i.e. ./results/step01
if [ -z "$step01_results_dir" ]; then
  step01_results_dir="$input_dir/results/step01"
fi
plate_results_dir="$step01_results_dir/$input_base" #i.e. ./results/step01/PO2Q01

echo "---- $0 $1 $2 ----"

### Input Validations ###

#make sure input_path was provided
if [ -z "$input_path" ]; then
  >&2 echo "Usage: $0 <name.fastq>"
  exit 1
fi
#make sure input_path exits
if [ ! -f "$input_path" ]; then
  >&2 echo "Error: Input file '$input_path' not found"
  exit 1
fi
#make sure input_path is fastq
if [[ ! "$input_path" == *.fastq ]]; then
  >&2 echo "Error: Input file '$input_path' must end with .fastq"
  exit 1
fi
if [ -d "$plate_results_dir" ]; then
  >&2 echo "Error: Directory $plate_results_dir already exists"
  exit 1
fi

echo "mkdir -p $plate_results_dir"
mkdir -p "$plate_results_dir" || exit 1 #exit on fail (already exists)

#Step 1. Take each fastq file (plate separated) and separate based on strain barcodes. Put in new folder which I precreated as the filename prefix.
cutadapt -a file:strain_barcodes.fasta -e 0 --overlap 13 -m 34 -M 37 -o "$plate_results_dir/{name}.fastq" "$input_path" || exit 1
