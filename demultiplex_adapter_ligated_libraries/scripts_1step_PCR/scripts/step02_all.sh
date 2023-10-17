#!/bin/bash

#Step 2. Take each fastq file (plate and strain separated) and separate based on well barcodes

# This script runs step 2 on all subfolders of the provided folder (i.e. ./results/* )

script_dir=$(dirname "$0")
step02_plate_script="$script_dir/step02_plate.sh"

results_dir="$1"

# if results_dir is not set, default to ./results
if [ -z "$results_dir" ]; then
  results_dir="./results"
#else if the results_dir was specified, make sure it exists
elif [[ ! -d "$results_dir" ]]; then
  >&2 echo "$results_dir does not exist or is not a directory"
  exit 1
fi

#determine step01 and step02 dirs after defaulting results_dir
step01_dir="$results_dir/step01"
step02_dir="$results_dir/step02"

if [ -d "$step02_dir" ]; then
  >&2 echo "Error: Directory $step02_dir already exists"
  exit 1
fi

#loop through the folders within step01
for plate_dir in $( find "$step01_dir" -type d -mindepth 1); do
  sh "$step02_plate_script" "$plate_dir" || exit 1 #exit on single error
done
