#!/bin/bash

#Step 3. Determine number of lines of each file (wc -l), divide by 4. This is your count. There are 96 per strain, per initial .fastq file.

results_dir="$1"

# script_dir=$(dirname "$0")

# if results_dir is not set, default to ./results
if [ -z "$results_dir" ]; then
  results_dir="./results"
#else if the results_dir was specified, make sure it exists
elif [[ ! -d "$results_dir" ]]; then
  >&2 echo "$results_dir does not exist or is not a directory"
  exit 1
fi

#determine step02 and step03 dirs after defaulting results_dir
step02_dir="$results_dir/step02"

results_file="$results_dir/Reads.tsv"

#Command Explanation:
# This command performs way better than looping through the files and running wc individually and echoing out the results.
# 1. wc: Do a word count of all *fastq files in the step02 folder
# 2. grep: Get rid of the total line
# 3. sed: Strip the leading spaces from the wc result, so we can easily split it into linecount and filename with awk
# 4. awk: Split the lines by dot (extension separator), space (wc -l separator), and slash (directory separator).
#         Then print the 3rd to last, 2nd to last, 1st (linecount), and 1st/4 (reads)
#         For example, the wc result "4 results/step02/P02Q1/ostA_H08.fastq"
#          will be broken down into: 4, results, step02, P02Q1, ostA_H08, fastq
#          and then will be printed as: "P02Q1  ostA_H08  4  1"
wc -l "$step02_dir"/*/*.fastq | grep -v total | sed "s/^[ \t]*//" | awk 'BEGIN{FS="/|[ ]"; OFS="\t"}{print $(NF-1),$(NF),$1,$1/4}' | sed -e 's/.fastq//g' > "$results_file"

echo ""
echo "$results_file generated"

# for file in $( find "$step02_dir" -type f -name '*.fastq' ); do
#   #Example file is ./results/step02/P02Q1/gcp_G05.fastq
#   basename=$(basename "$file" .fastq) #i.e. gcp_G05
#   dirname=$(dirname "$file") #i.e. ./results/step02/P02Q1
#   strain=$(basename "$dirname") #i.e. P02Q1
#
#   linecount=$(wc -l < "$file")
#   reads=$(( linecount / 4 ))
#   echo "$strain,$basename,$reads"
#   echo "$strain,$basename,$reads" >> Results.txt
# done
