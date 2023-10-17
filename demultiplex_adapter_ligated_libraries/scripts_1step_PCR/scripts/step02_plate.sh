#!/bin/bash

#Step 2. Take each fastq file (plate and strain separated) and separate based on well barcodes
# This script runs step 2 on a single plate folder (i.e. output/step01/P02Q1)

### Variable setup ###

input_dir="$1" #i.e. ./results/step01/P02Q1/
base=$(basename "$input_dir") #i.e. P02Q1
#TODO read result_dir from $2 with default
step01_dir=$(dirname "$input_dir") #i.e. ./results/step01/
results_dir=$(dirname "$step01_dir")

echo "---- $0 $1"

output_dir="$results_dir/step02/$base"

### Input Validations ###

if [ -z "$input_dir" ]; then
  >&2 echo "Usage: $0 <path/to/input/folder>"
  exit 1
fi
if [ -e "$output_dir" ]; then
  >&2 echo "Error: Directory $output_dir already exists"
  exit 1
fi

echo "mkdir -p $output_dir"
mkdir -p "$output_dir" || exit 1 #exit on fail (already exists)

for file in $( find "$input_dir" -type f -name '*.fastq'); do
  # echo "Step2: $file"
  base=$(basename "$file" .fastq)
  cutadapt -g file:well_barcodes.fasta -e 0 --overlap 8 -m 26 -M 29 -o "$output_dir/${base}_{name}.fastq" "$file" || exit 1
done
