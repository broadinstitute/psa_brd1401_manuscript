#!/bin/bash

### Variable setup ###

script_dir=$(dirname "$0")
step01_plate_sh="$script_dir/step01_plate.sh"

input_dir="$1" #directory to process. defaults to current

# if input_dir is not set, default to the current directory
if [ -z "$input_dir" ]; then
  input_dir="."
#else if the input_dir was specified, make sure it exists
elif [[ ! -d "$input_dir" ]]; then
  >&2 echo "Error:  $input_dir does not exist or is not a directory"
  exit 1
fi

#get optional results_dir from $2, with default
results_dir="$2" #i.e. ./results/step01
if [ -z "$results_dir" ]; then
  results_dir="$input_dir/results/step01"
fi

if [ -d "$results_dir" ]; then
  >&2 echo "Error: Directory $results_dir already exists"
  exit 1
fi

for file in $( find "$input_dir" -name '*.fastq' -maxdepth 1); do
  sh "$step01_plate_sh" "$file" "$results_dir" || exit 1 #exit on single error
done
