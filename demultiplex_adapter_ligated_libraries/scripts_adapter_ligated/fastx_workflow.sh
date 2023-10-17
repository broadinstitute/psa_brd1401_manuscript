#!/bin/bash

INPUT=$1
STEP_FROM=$2
STEP_TO=$3
DIR=$(dirname "$INPUT")
FILENAME=$(basename "$INPUT")

echo "DIR: $DIR"
echo "FILE: '$FILENAME'"

PATH=.:$PATH

DEBUG_CMD_LOG=debug.log
DEBUG_RESULT_LOG=debug.log
DEBUG_TO_STDOUT=true
RESULTS_TO_STDOUT=true

#Note. The script will cd to $DIR so all paths should be relative to there.
STRAIN_BCFILE="strain_barcodes.txt"
WELL_BCFILE="well_barcodes.txt"
PLATE_BCFILE="plate_barcodes.txt"

#STEP0_FOLDER="$DIR" #currently use the same as input folder
STEP1_FOLDER="Step_1"
STEP2_FOLDER="Step_2"
STEP3_FOLDER="Step_3"
STEP4_FOLDER="Step_4"
STEP5_FOLDER="Step_5"
STEP6_FOLDER="Step_6"
STEP7_FOLDER="Step_7"
STEP8_FOLDER="Step_8"

STEP0_INPUT="$FILENAME"
STEP1_TRIMMED="$STEP1_FOLDER/trimmed92"
STEP2_BCFILE=$PLATE_BCFILE
STEP2_PREFIX="$STEP2_FOLDER/Split"
STEP2_UNMATCHED="${STEP2_PREFIX}unmatched"
STEP3_REVERSED="$STEP3_FOLDER/reversecomp"
STEP3_PREFIX="$STEP3_FOLDER/Split"
STEPS_23_COMBINED_PREFIX="$STEP2_FOLDER/Combined"

#Note: for Step 4 prefix, only include the filename prefix, not the directory. The runtime prefix will be $FOLDER/$PLATE/$PREFIX
STEP4_BCFILE=$WELL_BCFILE
STEP4_PREFIX="Well#"
STEP6_BCFILE=$STRAIN_BCFILE
#STEP6_PREFIX=""
#STEP6_SUMMARY_TXT="Summary.txt" #this variable is not used.
#STEP6_SUMMARY_CSV="Summary.csv" #this variable is not used.

rm $DEBUG_CMD_LOG 2>> /dev/null
rm $DEBUG_RESULT_LOG 2>> /dev/null

function usage() {
  echo "Usage: $0 <input_file> [step_from] [step_to]"
  exit -1
}
function debugResult() {
  if [[ "$1" = "" ]]; then
    return
  fi
  #comment out the echo line to disable debug
  if [[ "$RESULTS_TO_STDOUT" = "true" ]]; then
    echo "$1";
  fi
  echo "$1" >> $DEBUG_RESULT_LOG 2>> /dev/null
}
function debug() {
  if [[ "$DEBUG_TO_STDOUT" = "true" ]]; then
    echo "$1";
  fi
  echo "$1" >> $DEBUG_CMD_LOG 2>> /dev/null
}
#runs a command string and sends results to the filename specified in the 2nd param
#if no 2nd param is passed, results are sent to DEBUG_RESULT_LOG
function runCmd() {
  if [ -n "$2" ]; then
    FWD="$2"
  else
    FWD="$DEBUG_RESULT_LOG"
  fi
  debug "$1 > $FWD"
  stderr=$(eval "$1" 2>&1 1>>"$FWD")  # Run command. stderr is captured.
  if [ ! -z "$stderr" ]; then
    debugResult "ERROR: $stderr"
    echo "An error occurred. Exiting script."
    exit 1
  fi
}
#this is the old version of runCmd which sets stdout variable.
#it does not have any error handling
function runCmdOld() {
  debug "$1"
  stdout=$(eval "$1") #Run command. stdout is captured. stderr is ignore
  debugResult "$stdout"
}
function zmkdir() {
  if [ ! -d "$1" ]; then
    runCmd "mkdir -p '$1'"
  fi
}

#temporarily change the IFS to allow looping over filenames with spaces
OIFS="$IFS"
IFS=$'\n'

if [[ ! -f "$INPUT" ]]; then
  echo "ERROR: input file $INPUT does not exist"
  usage
fi

if [[ -z "$STEP_FROM" ]]; then
  STEP_FROM=0
else
  debug "Starting at Step $STEP_FROM"
fi

if [[ -z "$STEP_TO" ]]; then
  STEP_TO=99
else
  debug "Ending at Step $STEP_TO"
fi

cd "$DIR" || exit

debug "Changed working directory to $(pwd)"

#Step 0
STEP=0
##figure out STEP0_OUTPUT regardless of whether step 0 will be run
extension0="${FILENAME##*.}" #get the extension, if any
if [ "$extension0" == "fastq" ]; then
  #the input file already has the fastq extension. don't need to rename
  STEP0_INPUT=$FILENAME
else
  #the input file does not have the fastq extension. need to rename
  STEP0_INPUT="$FILENAME.fastq" #if there is no extension on input, use this
fi
#find/replace to add _cutadapt to the output filename right before the extension
STEP0_OUTPUT=${STEP0_INPUT/.fastq/_cutadapt.fastq}

#run step 0 unless STEP_FROM is higher
if [[ "$STEP_FROM" -le "$STEP" ]]; then
  debug ""
  debug "Step 0:"
  if [[ "$FILENAME" != "$STEP0_INPUT" ]]; then
    runCmd "mv '$FILENAME' '$STEP0_INPUT'"
  fi
  runCmd "cutadapt -a GCAGATCGGA -o '$STEP0_OUTPUT' '$STEP0_INPUT'"
fi

#Step 1
STEP=1
if [[ ("$STEP_FROM" -le "$STEP") && ("$STEP_TO" -ge "$STEP") ]]; then
  debug ""
  debug "Step 1:"
  zmkdir $STEP1_FOLDER
  runCmd "fastx_trimmer -l 92 -i '$STEP0_OUTPUT' -o '$STEP1_TRIMMED' -Q33"
fi

#Step 2
STEP=2
if [[ ("$STEP_FROM" -le "$STEP") && ("$STEP_TO" -ge "$STEP") ]]; then
  debug ""
  debug "Step 2:"
  zmkdir $STEP2_FOLDER
  #The number of output files from Step 2 should match the number of Plate barcodes (12) plus an optional unmatched file
  runCmd "cat '$STEP1_TRIMMED' | fastx_barcode_splitter.pl --bcfile '$STEP2_BCFILE' --prefix '$STEP2_PREFIX' --eol --exact"
fi

#Step 3
STEP=3
if [[ (-f "$STEP2_UNMATCHED") ]]; then

  if [[ ("$STEP_FROM" -le "$STEP") && ("$STEP_TO" -ge "$STEP") ]]; then
    zmkdir $STEP3_FOLDER

    debug ""
    debug "Unmatched Check: '$STEP2_UNMATCHED' file exists. Will reverse and combine."

    debug ""
    debug "Step 3a:"
    runCmd "fastx_reverse_complement -i '$STEP2_UNMATCHED' -o '$STEP3_REVERSED' -Q33"

    #The number of output files from Step 3 should match the number of Plate barcodes (12)
    debug ""
    debug "Step 3b:"
    runCmd "cat '$STEP3_REVERSED' | fastx_barcode_splitter.pl --bcfile '$STEP2_BCFILE' --prefix '$STEP3_PREFIX' --eol --exact"

    debug ""
    debug "Step 3c: (Using File Search '*$STEP3_PREFIX*', excluding unmatched)"
    #combine each corresponding Split Plate file from steps 2 & 3 (i.e. Step2/Plate1 + Step3/Plate1)
    for step3plate in $(find . -wholename "*$STEP3_PREFIX*" | grep -v unmatched); do
      step2plate=${step3plate/${STEP3_PREFIX}/${STEP2_PREFIX}} #find/replace to get the corresponding step 2 filename
      combined=${step3plate/${STEP3_PREFIX}/${STEPS_23_COMBINED_PREFIX}} #find/replace to determine the combined filename
      runCmd "cat '${step3plate}' '${step2plate}' > '${combined}'"
    done
  fi
  STEP4_SEARCH=$STEPS_23_COMBINED_PREFIX
else
  debug ""
  debug "Check: $STEP2_UNMATCHED does not exist. No reversing and combining required."
  STEP4_SEARCH=$STEP2_PREFIX
fi

#Step 4
STEP=4 #Split wells for each plate file
if [[ ("$STEP_FROM" -le "$STEP") && ("$STEP_TO" -ge "$STEP") ]]; then
  SEARCH="*${STEP4_SEARCH}*"
  debug ""
  debug "Step 4: (Using File Search '$SEARCH', excluding '$STEP2_BCFILE' and unmatched)"
  for step4plate in $(find . -wholename "$SEARCH" | grep -v "'$STEP2_BCFILE'" | grep -v unmatched); do
if [[ ! -s "$step4plate"  ]]; then
continue
fi
    PLATENAME=${step4plate/${STEP4_SEARCH}/} #remove the Step2/3 prefixes
    PLATENAME=$(basename $PLATENAME) #remove any excess directory info
    #debug "Plate=$PLATENAME"
    debug ""
    zmkdir $STEP4_FOLDER/${PLATENAME}
    runCmd "cat '${step4plate}' | fastx_barcode_splitter.pl --bcfile '$STEP4_BCFILE' --prefix '$STEP4_FOLDER/${PLATENAME}/$STEP4_PREFIX' --bol --exact"
  done
fi

#Step 5
STEP=5 #Trim to keep strain barcode
if [[ ("$STEP_FROM" -le "$STEP") && ("$STEP_TO" -ge "$STEP") ]]; then
  SEARCH="*${STEP4_FOLDER}/*/${STEP4_PREFIX}*"
  debug ""
  debug "Step 5: (Using File Search '$SEARCH', excluding unmatched)"
  for STEP5_INPUT in $(find . -wholename "$SEARCH" | grep -v unmatched); do
  if [[ -s "$STEP5_INPUT" ]] ; then
  STEP5_OUTFILE=${STEP5_INPUT/$STEP4_FOLDER/$STEP5_FOLDER} #find/replace
  STEP5_OUTDIR=$(dirname "$STEP5_OUTFILE")
  debug ""
  zmkdir "$STEP5_OUTDIR"
  runCmd "fastx_trimmer -f 35 -l 58 -i '$STEP5_INPUT' -o '$STEP5_OUTFILE' -Q33"
  fi
  done
fi

#Step 6
STEP=6
if [[ ("$STEP_FROM" -le "$STEP") && ("$STEP_TO" -ge "$STEP") ]]; then

  #delete the existing summary files (txt and and csv files).
  find ${STEP6_FOLDER} -wholename '*Summary*' -exec rm {} \;

  SEARCH="*${STEP5_FOLDER}/*/${STEP4_PREFIX}*"
  debug ""
  debug "Step 6: (Using File Search '$SEARCH')"
  for STEP6_INPUT in $(find . -wholename "$SEARCH"); do
if [[ ! -s "$STEP6_INPUT"  ]]; then
continue
fi
    STEP6_OUTFILE=${STEP6_INPUT/$STEP5_FOLDER/$STEP6_FOLDER} #find/replace
    #STEP6_OUTFILE="$STEP6_OUTFILE/" #add a slash to the end to create subfolders
    PLATEDIR=$(dirname $STEP6_OUTFILE) #removes the well portion of the directory
    PLATENAME=$(basename $PLATEDIR)

    debug ""
    zmkdir $PLATEDIR
    #use runCmdOld to assign stdout
    runCmdOld "cat '$STEP6_INPUT' | fastx_barcode_splitter.pl --bcfile '$STEP6_BCFILE' --prefix '${STEP6_OUTFILE}-' --bol --exact"

    #Append to the summary file, removing the header and total rows
    echo "$stdout" | grep -v Barcode | grep -v total >> "${PLATEDIR}Summary.txt"

  done

  #clear stdout variable when done the loop
  stdout=

  SEARCH="*${STEP6_FOLDER}/*Summary.txt"
  debug ""
  debug "Step 6b: (Using File Search '$SEARCH' | grep -v unmatched'"
  for TXT in $(find . -wholename "$SEARCH"); do
    #Convert the default Summary file to CSV. Use tabs(\t), slash(/), hash(#), and dash(-) as the awk delimeters
    #Parses the Well from the 2nd to last part of filename
    #    Input Example: ACINE   123     /path/Plate1/Well#A1#ACINE
    #    Awk Variables: $1=ACINE,$2=123,...,$(NF-1)=A1,$NF=ACINE, where NF means the last item
    #    Output Example: ACINE,Plate1,A1,123
    CSV=${TXT/txt/csv}
    debug "converting '$TXT' to '$CSV'"
    cat "$TXT" | awk -F'\t|/|#|-' '{print $1","$(NF-3)","$(NF-1)","$2}' > "$CSV"

  done

fi

#Step 7
STEP=7
if [[ ("$STEP_FROM" -le "$STEP") && ("$STEP_TO" -ge "$STEP") ]]; then
  debug ""
  debug "Step 7. Merging csv files into $STEP7_FOLDER folder"
  zmkdir "$STEP7_FOLDER"

  sort -k 1,1 -k 3,3 -f -o "$STEP7_FOLDER"/all_raw_data_combined.csv "$STEP6_FOLDER"/*.csv
  echo "Wrote $STEP7_FOLDER/all_raw_data_combined.csv"
  CWD=$(pwd) #save the CWD to cd back to. Note, don't use DIR
  cd "$STEP7_FOLDER" #cd to Step 7 folder for the next couple commands
  cat all_raw_data_combined.csv | awk -F, '{ print $0","log($4+1)/log(2)}' > all_data.csv
  echo "Wrote $STEP7_FOLDER/all_data.csv"
  cat all_data.csv | grep -E 'CIP_128|DMSO' > data_for_Z_factor.csv
  echo "Wrote data_for_Z_factor.csv"
  cd "$CWD" || exit #cd back to the main directory. Note. don't use $DIR for this since it may not be absolute
fi

STEP=8
if [[ ("$STEP_FROM" -le "$STEP") && ("$STEP_TO" -ge "$STEP") ]]; then
  debug ""
  debug "Step 8. ."
  zmkdir "$STEP8_FOLDER"
  barcodes=$(cat "$STRAIN_BCFILE" | awk '{print $1}')
  for barcode in $barcodes; do
    runCmd "cat '$STEP7_FOLDER/all_data.csv' | grep -E '$barcode' > '$STEP8_FOLDER/$barcode.csv'"
    #echo "Wrote $barcode.csv"
  done
fi

IFS="$OIFS" #revert IFS to the original value
