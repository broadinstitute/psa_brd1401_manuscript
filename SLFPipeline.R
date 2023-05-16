#10/28/20, jbagnall@broadinstitute.org
#Script for running SLF and ZZ-score pipeline
#Assumes input data is counts

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(source('/Functions_CompoundScreenPipelineSLF_210927.R'))

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  countdatapath1 <- args[1]
  savefilepath <- args[2]
  count_exact1 <- args[3]
  untreated_name <- args[4]
  intcon_name = args[5]
  lowcountfilter = args[6]
  lowcountfilter_untreated = args[7]
  keep_colnames = unlist(strsplit(args[8], split = ','))

  #Clean count data from Concensus2
  countdata = cleanfromConcensus2(rawcountpath = countdatapath1, count_exact = count_exact1)

  #Calculate SLF and ZZ-scores
  savefilepath_rds = paste0(gsub("\\..*","",savefilepath),".rds")
  compScreenPipeline(countdata, untreatedname = untreated_name, intconname = intcon_name, comp_conc_separator = ":", lowwellcount = lowcountfilter, low_untreated_count = lowcountfilter_untreated, medianLFC = FALSE, newSchema = F, savefilename = savefilepath_rds)
  
  #Save CSV file with chosen columns
  savefilepath_csv = paste0(gsub("\\..*","",savefilepath),".csv")
  subsetSLF(screendata_path = savefilepath_rds, keep_columns = keep_colnames, savefilename = savefilepath_csv)
} 

main()
