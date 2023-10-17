# cutadapt-demultiplexing

Some scripts for demultiplexing with cutadapt.

## Usage

	./scripts/cutadapt_flow.sh [dir_to_search_for_fastq='.'] [from=1] [to=*]
	
* or alternatively run with sh: `sh scripts/cutadapt_flow.sh ...`

### Examples

To run the whole flow by search the current directory for *.fastq files, use either of these:

	./scripts/cutadapt_flow.sh
	
Or to explicitly set the current directory:

	./scripts/cutadapt_flow.sh .
	
Run with a different input directory:

	./scripts/cutadapt_flow.sh ./input

Run from step 2 (the directory is required here for parameter ordering):

	./scripts/cutadapt_flow.sh . 2
	
Run from from 2 to 3 (the directory is required here for parameter ordering):

	./scripts/cutadapt_flow.sh . 2 3
	
## Alternative / Underlying Scripts

There are some underlying scripts that can also be called directly.

**step01_all.sh**

This script runs step 1. The search directory will default to the current.

	./scripts/step01_all.sh [dir_to_search_for_fastq='.']

**step01_plate.sh**

This script runs step 1 for a single plate (.fastq) file. The file path is required.

	./scripts/step01_plate.sh <path/to/plate.fastq>
	
**step02_all.sh**

This script runs step 2. The path to the results folder can be specified, and defaults to the local ./results folder. 
It expects a step01 folder inside the results directory.

	./scripts/step02_all.sh [/path/to/results="./results"]

**step02_plate.sh**

This script runs step 2 for a single plate. The path to the plate's step01 results is required.

	./scripts/step02_plate.sh <path/to/results/step01/plate/>

**step03_all.sh**

Runs step 3. The path to the results folder can be specified, and defaults to the local ./results folder. 
It expects a step02 folder inside the results directory.


