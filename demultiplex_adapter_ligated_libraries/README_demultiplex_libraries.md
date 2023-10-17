# Demultiplexing Adapter-Ligated Libraries

Libraries are constructed using short PCR primers containing annealing regions with ‘well’ or ‘plate’ barcodes built in. These short primers result in a small PCR fragment (one end contains the well barcode and the other contains the plate barcode) that is prepared for Illumina sequencing with the ligation of a Y-adapter. This Y-adapter is made of 2 primers that contain either the Illumina P7 or P5 adapters (each at the top of the Y) and a short annealing region that also serves as the sequencing primer (forming the bottom of the Y). 

Since the Y adapters become ligated at both ends of the PCR fragment, a read can begin with either the well barcode or plate barcode during the Illumina sequencing, making the demultiplexing steps slightly complex (see Step 3 below for more details). However, from a single 100 base read, the single output fastq data file contains the plate barcode, well barcode, and strain barcode and no demultiplexing is done by the Illumina machine. 

## Software required:
cutadapt
fastx-toolkit-0.0.13

## Files required (all in the same folder):
fastq file from Illumina
./fastx_workflow.sh
./strain_barcodes.txt
./well_barcodes.txt
./plate_barcodes.txt

## Command-line usage:
```./fastx_workflow.sh inputfilename.fastq```

## Explanation of processing:
Executing the above command will run a series of processing steps beginning with the fastq input file, outlined below.

### Step 0
This pre-processing step uses cutadapt to remove tailing adapter sequence from the read and writes a new fastq file:
```cutadapt -a GCAGATCGGA -o '$STEP0_OUTPUT' '$STEP0_INPUT```

### Step 1
This step uses a trimmer function from the fastx-toolkit to ensure all reads have been properly processed by Step 0 and the resulting fragment is 92 bases long:
```fastx_trimmer -l 92 -i '$STEP0_OUTPUT' -o '$STEP1_TRIMMED' -Q33```

### Step 2
This step begins the demultiplexing using the fastx-toolkit barcode splitter and the plate barcodes file:
```cat '$STEP1_TRIMMED' | fastx_barcode_splitter.pl --bcfile '$STEP2_BCFILE' --prefix '$STEP2_PREFIX' --eol –exact```

--eol -exact indicates that exact matching of barcodes at the end of line is required

### Step 3
This step has three substeps (A-C)

#### Step 3A finds the “unmatched” file from Step 2. This file contains half the usable reads that were not found in Step 2 because these reads begin with the plate barcode in reverse complement instead of ending with the plate barcodes. This is due to the Y-adapter, which is ligated at both ends of the PCR fragment, thus sequencing occurs from both ends. 3A uses the unmatched file from Step 2 and uses fastx reverse complement tool:
```fastx_reverse_complement -i '$STEP2_UNMATCHED' -o '$STEP3_REVERSED' -Q33```

#### Step 3B is the same as Step 2, but on this newly created file:
```cat '$STEP3_REVERSED' | fastx_barcode_splitter.pl --bcfile '$STEP2_BCFILE' --prefix '$STEP3_PREFIX' --eol --exact```

#### Step 3C combines both the files for each plate using UNIX commands:
```cat '${step3plate}' '${step2plate}' > '${combined}'```

### Step 4
This step uses the combined file from Step 3C and splits the wells based on the well barcode file using the fastx barcode splitter:
```cat '${step4plate}' | fastx_barcode_splitter.pl --bcfile '$STEP4_BCFILE' --prefix '$STEP4_FOLDER/${PLATENAME}/$STEP4_PREFIX' --bol --exact```

--bol –exact tells the barcode_splitter that the barcodes are at the beginning of the line and exact matching is required.

### Step 5
This step prepares the fastq file for the strain demultiplexing. The strain barcodes are in the middle of the read; the fastx trimmer function trims both ends so that only the strain barcode remains:
```fastx_trimmer -f 35 -l 58 -i '$STEP5_INPUT' -o '$STEP5_OUTFILE' -Q33```

### Step 6
This step now uses the input file and splits it based on the strain barcode file using barcode splitter, as used in the steps above:
```cat '$STEP6_INPUT' | fastx_barcode_splitter.pl --bcfile '$STEP6_BCFILE' --prefix '${STEP6_OUTFILE}-' --bol --exact```

### Step 7
This file counts the fastq files and combines the data into a single file “all_data.csv.” This file is structured in columns: strain, plate, well, read counts, log reads. A file “data_for_Z_factor.csv” is also created by selecting all DMSO and CIP_128 wells. This can be modified by changing the script labels and well barcode naming:
```cat all_data.csv | grep -E 'CIP_128|DMSO' > data_for_Z_factor.csv```

### Step 8
This step writes a .csv file for each strain, using the “all_data.csv” file:
```cat '$STEP7_FOLDER/all_data.csv' | grep -E '$barcode' > '$STEP8_FOLDER/$barcode.csv```


# Demultiplexing 1-Step PCR Libraries

Libraries are constructed using long primers containing standard Illumina P7 and P5 adapters and built-in index read primers. Since these primers are long, they need to be highly pure for sequencing. This is cost-prohibitive for high-throughput sequencing (96 plate primers and 96 well primers needed in large quantities). 

During the Illumina sequencing, a single read (50 or 100 bases) using the standard sequencing primer and the index reads are performed and the “plate” demultiplexing is performed by the Illumina machine. The result is 96 plate files well and strain information. These files are demultiplexed using the custom script “cutadaptflow.sh” and the cutadapt software.

## Software required:
cutadapt

## Files required (scripts in subfolder named scripts):
96 fastq plate files from Illumina
./scripts/cutadapt_flow.sh
./scripts/step01_all.sh
./scripts/step01_plate.sh
./scripts/step02_all.sh
./scripts/step02_plate.sh
./scripts/step03_all.sh
./strain_barcodes.fasta
./well_barcodes.fasta

## Command-line usage:
```./scripts/cutadapt_flow.sh ./input_directory_containing_fastq_files```

## Explanation of processing:
Executing the above command will run the remaining .sh files sequentially, beginning with Step 1.

### Step 1
The first script run is ```./scripts/step01_all.sh``` that runs step 1 for all plates by running ```./scripts/step01_plate.sh``` for each fastq file in the directory. This was designed so the option of running only the “…plate.sh” command on a single plate is possible, see the readme.md for further details.

The sequence read processing step that is run on each fastq file is:
```cutadapt -a file:strain_barcodes.fasta -e 0 --overlap 13 -m 34 -M 37 -o "$plate_results_dir/{name}.fastq" "$input_path" || exit 1```

This command uses cutadapt to find each strain barcode with the following parameters:
-e 0 ensures exact matching; this can be altered if mismatches in the strain barcodes is desired
--overlap 13 is chosen for 50bp reads, since not all the barcode is sequenced. For longer reads (e.g., 100bp reads) the full barcode length can be matched by selecting --overlap 24 in the case of 24bp barcodes.
-m 34 and -M37 ensures that only full reads are kept. The products in this 50bp read example after trimming the 13bp barcode should have 34-37bp remaining.

### Step 2
Run (```./scripts/step02_all.sh```), which again, runs ```./scripts/step02_plate.sh``` on each file.

The sequence read processing step that is run on each fastq file is:
```cutadapt -g file:well_barcodes.fasta -e 0 --overlap 8 -m 26 -M 29 -o "$output_dir/${base}_{name}.fastq" "$file" || exit 1```

This command uses cutadapt to find each well barcode with the following parameters:
-e 0 ensures exact matching; this can be altered if mismatches in the strain barcodes is desired
--overlap 8 is chosen because the well barcodes are 8bp
-m 26 and -M29 ensures that only full reads are kept. The products in this 50bp read example after trimming the 13bp strain barcode and 8bp well barcode should have 26-29bp remaining.

### Step 3 
Count the barcodes for each plate/strain/well combination and writes the results to a file. A detailed explanation is found in ```step03_all.sh```
