#!/bin/bash
#$ -N SLF_Pipeline
#$ -l h_vmem=16G
#$ -l h_rt=72:00:00

# Use this script to run SLFPipeline.R to calculate SLF and ZZ-scores
# 10/28/20, jbagnall@broadinstitute.org

source /broad/software/scripts/useuse
use UGER
use R-3.5

##############################################################################################
#        Please type your inputs below: (careful not to put extra spaces)                   #
##############################################################################################
#Path to input summary count table
countdatapath1='/idi/cgtb/jbagnall/Brad_pseudomonas/screen1/screen1_rawcounts_for_SLF_pipeline_230117.csv'

#Path to output where you want to save your final SLF and ZZ-score file. No extension necessary. Pipeline will save an .rds and .csv file.
savefilepath="/idi/cgtb/jbagnall/Brad_pseudomonas/screen1/screen1_processed_test_230117"


###########################################################################################
#               Other parameters (shouldn't need to be changed for each run):             #
###########################################################################################
#Do you want to use exact matches for counts (count_exact)? If so, type "T" otherwise type "F". 
#"F" means that it will use count_all, which includes mismatches (depends on count script)
count_exact1="T"

#What is the name of the negative control compound? Usually is "DMSO"
untreated_name="DMSO"

#What is the name of the internal spike-in controls. Should be a common word shared by them all, but not by hypomorphs, e.g. "control" or "PCR"
intcon_name="PCR"

#What column names to keep in the final .csv output file. Default is to keep as written here. No spaces, and commas between column names
keep_colnames="strain,compound,concentration,plate_name,row,column,count,rep,wellcount,wellcountfrac,std_lf,zscore_stdlf,zscore_stdlf2,correlation"

#Remove wells with counts less than or equal to this value (lowcountfilter), consider level of spike-ins. 
#Otherwise set to -1 to not filter out any wells based on total well count.
lowcountfilter=-1

#Remove DMSO control wells with counts less than or equal to this value (lowcountfilter_untreated), consider level of spike-ins.
lowcountfilter_untreated=100

#########################################################################################
#                      Run SLF/ZZ pipeline                                              #
#########################################################################################
Rscript SLFPipeline.R "$countdatapath1" "$savefilepath" "$count_exact1" "$untreated_name" "$intcon_name" $lowcountfilter $lowcountfilter_untreated "$keep_colnames" 

