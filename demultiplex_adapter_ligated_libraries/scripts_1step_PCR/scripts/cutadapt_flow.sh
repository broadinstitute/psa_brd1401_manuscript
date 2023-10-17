#!/bin/bash

### Variable setup ###

input_dir="$1" #full path to the file, including any directories
script_dir=$(dirname "$0")

step_from=$2
step_to=$3

# if input_dir is not set, default to the current directory
if [ -z "$input_dir" ]; then
  input_dir="."
#else if the input_dir was specified, make sure it exists
elif [[ ! -d "$input_dir" ]]; then
  >&2 echo "$input_dir does not exist or is not a directory"
  exit 1
fi

#initialize results_dir after defaulting input_dir
results_dir="$input_dir/results"
results_file="$results_dir/Reads.tsv"

if [ -z "$step_from" ]; then
 step_from=1
fi

if [ -z "$step_to" ]; then
 step_to=9999
fi

function checkResultDir() {
  if [ -d "$1" ]; then
    >&2 echo "Error: Directory $1 already exists"
    exit 1
  fi
}

#see if some results are there before we start. this check also depends on the desired step from and to.

step=1
if [ $step -ge $step_from ] && [ $step -le $step_to ]; then
  checkResultDir "$results_dir/step01"
fi

step=2
if [ $step -ge $step_from ] && [ $step -le $step_to ]; then
  checkResultDir "$results_dir/step02"
fi

step=3
if [ $step -ge $step_from ] && [ $step -le $step_to ]; then
  if [ -e "$results_file" ]; then
    >&2 echo "Error: File $results_file already exists"
    exit 1
  fi
fi

# now we start the actual steps

step=1
if [ $step -ge $step_from ] && [ $step -le $step_to ]; then
  sh "$script_dir/step01_all.sh" "$input_dir" || exit 1
fi

step=2
if [ $step -ge $step_from ] && [ $step -le $step_to ]; then
  sh "$script_dir/step02_all.sh" "$results_dir" || exit 1
fi

step=3
if [ $step -ge $step_from ] && [ $step -le $step_to ]; then
  sh "$script_dir/step03_all.sh" "$results_dir" || exit 1
fi
