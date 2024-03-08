#Last edited 9/27/21
#jbagnall@broadinstitute.org

library(dplyr)
library(tidyr)
library(stringr)

cleanfromConcensus2 <- function(rawcountpath, count_exact = T){
  #Takes output from count pipeline with the following header:
  #pool,i5_id,i5_plate,plate_barcode,i7_id,plate_quadrant,well_pos,Broad_Sample,mg_per_ml,mmoles_per_liter,strain_id,strain_gene,count_exact,count_nonexact,count_all
  #if count_exact == T, then uses count_exact. Otherwise, uses count_all (includes barcode mismatches)
  #Adjusts column names
  #outputs data frame to be read by compScreenPipeline for calculation of ZZ and SLF
  rawcounts = read.csv(rawcountpath, stringsAsFactors = F)
  
  if(is.na(count_exact) | count_exact %in% c("f", "F", "false", "FALSE")){
    count_exact = F
  }else if(count_exact %in% c("t", "true", "TRUE", "T")){
    count_exact = T
  }
  
  if(count_exact){
    colnames(rawcounts)[colnames(rawcounts) == "count_exact"] = "count"
  }else{
    colnames(rawcounts)[colnames(rawcounts) == "count_all"] = "count"
  }
  # rawcounts = select(rawcounts, Broad_Sample, mmoles_per_liter, plate_barcode, well_pos, strain_gene, count)
  rawcounts = select(rawcounts, Broad_Sample, mmoles_per_liter, plate_barcode, well_pos, strain_gene, count, pool) 
  
  #Concatenate pool to plate_name to make sure it's unique
  rawcounts = mutate(rawcounts, plate_barcode = paste(plate_barcode, pool, sep = "_"))
  rawcounts = select(rawcounts, -pool)
  
  rawcounts = dplyr::filter(rawcounts, !(strain_gene %in% c("", "#N/A", "NA"))) #remove cases with no strain_gene name (empty string)
  rawcounts = dplyr::filter(rawcounts, !(Broad_Sample == "")) #remove cases with no compound (empty string)
  rawcounts$row = gsub("\\d+", '', rawcounts$well_pos)
  rawcounts$column = gsub("\\D+", '', rawcounts$well_pos)
  rawcounts = select(rawcounts, -well_pos)
  colnames(rawcounts)[colnames(rawcounts) == "Broad_Sample"] = "compound"
  colnames(rawcounts)[colnames(rawcounts) == "mmoles_per_liter"] = "concentration"
  colnames(rawcounts)[colnames(rawcounts) == "plate_barcode"] = "plate_name"
  colnames(rawcounts)[colnames(rawcounts) == "strain_gene"] = "strain"
  return(rawcounts)
}

compScreenPipeline <- function(rawcounts_df, untreatedname = "untreated", intconname = "control", comp_conc_separator = ":", lowwellcount = 10, low_untreated_count = 100, medianLFC = FALSE, newSchema = F, savefilename = NA, saveCSV = F){ 
  #Input raw count data frame, 
  #Returns a data frame with compound (BRD-xxxxxxx), strain, plate_name, row, column, replicate number, 
  # strain counts per well (effect size?), total counts in wells, log fraction, 
  #"corrected" log fraction (normalized by plate and untreated), standardized residuals
  #p-value, replicate correlation, rounded concentration
  #If newSchema == T, then it will save the file with the new schema column names
  #low_untreated_count is the count threshold below which the untreated wells will be removed from the SLF calculation. Assumed to be a technical error, and not representative of the variation within the plate
  
  library(dplyr)
  library(tidyr)
  
  if(is.na(savefilename)){
    stop('No output file path specified')
  }

  convertSchema(rawcounts_df, to = "old")
  
  #Subset raw data frame
  if("well" %in% colnames(rawcounts_df)){
    rawcounts_df$row = gsub("\\d+", '', rawcounts_df$well)
    rawcounts_df$column = gsub("\\D+", '', rawcounts_df$well)
    rawcounts_df = select(rawcounts_df, -well)
  }
  
  rawcounts_subset <- dplyr::select(rawcounts_df, compound, concentration, strain, plate_name,  row, column, count)

  #Convert factors to characters
  factorcols <- sapply(rawcounts_subset, is.factor)
  rawcounts_subset[factorcols] <- lapply(rawcounts_subset[factorcols], as.character)
  remove(factorcols)
  
  rawcounts_subset <- as.data.frame(rawcounts_subset)
  
  #Remove cases where compound is NA
  print("Remove cases where compound is NA")
  rawcounts_subset <- dplyr::filter(rawcounts_subset, !is.na(compound))
  
  #Remove cases where strain is empty
  print("Remove cases where strain name is empty or NA")
  rawcounts_subset = dplyr::filter(rawcounts_subset, !(strain %in% c("", "#N/A"))) #remove cases with no strain_gene name
  
  #Check if control spike-in name is included
  if(!any(grepl(intconname, unique(rawcounts_subset$strain)))){
    warning('No spike-in controls found in strain names')
  }
  
  #Sum counts from multiple sequencing attempts on same plate & well
  print("Sum counts from multiple sequencing attempts")
  rawcounts_subset <- rawcounts_subset %>%
    group_by(strain, compound, concentration, plate_name, row, column) %>%
    summarise(count = sum(count)) %>%
    ungroup()
  
  #Make sure negative controls exist
  if(!(untreatedname %in% rawcounts_subset$compound)){
    stop("Error: No DMSO wells found. Check input name of negative control compound")
  }
  
  #BRD-ID only available if registered with compound management, should maintain original labels
  #Added relevant BRD_IDs for annotation purposes
  if(any(rawcounts_subset$compound == "rifampin" | rawcounts_subset$compound == "rifampicin" | rawcounts_subset$compound == "Rifampicin" |rawcounts_subset$compound == "Rifampin" )){
    print("Renaming rifampin to BRD-K01507359-rifampin")
    #Check if K01507359 should also be changed, But don't want to affect batch number
    rawcounts_subset$compound <- replace(rawcounts_subset$compound, rawcounts_subset$compound == "rifampin", "BRD-K01507359-rifampin") #Replace rifampin with Broad ID
    rawcounts_subset$compound <- replace(rawcounts_subset$compound, rawcounts_subset$compound == "rifampicin", "BRD-K01507359-rifampin") #Replace rifampin with Broad ID
    rawcounts_subset$compound <- replace(rawcounts_subset$compound, rawcounts_subset$compound == "Rifampicin", "BRD-K01507359-rifampin") #Replace rifampin with Broad ID
    rawcounts_subset$compound <- replace(rawcounts_subset$compound, rawcounts_subset$compound == "Rifampin", "BRD-K01507359-rifampin") #Replace rifampin with Broad ID
    
  }
  if(any(rawcounts_subset$compound == "trimethoprim")){
    print("Renaming trimethoprim to BRD-K07208025-trimethoprim")
    rawcounts_subset$compound <- replace(rawcounts_subset$compound, rawcounts_subset$compound == "trimethoprim", "BRD-K07208025-trimethoprim") #Replace trimethoprim with Broad ID
  }
  if(any(rawcounts_subset$compound == "methotrexate")){
    print("Renaming methotrexate to BRD-K59456551-methotrexate")
    rawcounts_subset$compound <- replace(rawcounts_subset$compound, rawcounts_subset$compound == "methotrexate", "BRD-K59456551-methotrexate") #Replace trimethoprim with Broad ID
  }
  if(any(rawcounts_subset$compound == "drug_4592")){
    print("Renaming drug_4592 to drug-4592") #to remove the underscore
    rawcounts_subset$compound <- replace(rawcounts_subset$compound, rawcounts_subset$compound == "drug_4592", "drug-4592") #Replace trimethoprim with Broad ID
  }
  
  #Removing wells with "none" treatment
  print("Removing wells with 'none' or '' as compound")
  rawcounts_subset <- dplyr::filter(rawcounts_subset, compound != "none")
  rawcounts_subset <- dplyr::filter(rawcounts_subset, compound != "") #make sure '' is not DMSO
  
  #Removing NA plates
  if(any(is.na(rawcounts_subset$plate_name))){
    print("Removing plates with name NA")
    rawcounts_subset <- dplyr::filter(rawcounts_subset, !is.na(plate_name))
  }
  
  #Removing plates with no name
  if(any(rawcounts_subset$plate_name == "")){
    print("Removing plates with no name (empty string)")
    rawcounts_subset <- dplyr::filter(rawcounts_subset, plate_name != "")
  }
  
  #Calcualate rounded concentrations (may add ranks of concentrations)
  rawcounts_subset <- mutate(rawcounts_subset, conc_round = signif(concentration, digits=1))
  rawcounts_subset <- mutate(rawcounts_subset, comp_conc = paste(compound, conc_round, sep = comp_conc_separator)) #changed from _ to : on 9/25/19
  
  #Remove wells with technical error
  #Remove all wells where intcon counts are 0
  print("Removing wells where either spike-in control has 0 counts")
  if((intconname != "") & (intconname != " ") & !is.na(intconname)){ #only run if intconname is not empty string and not NA
    intcons <- dplyr::filter(rawcounts_subset, grepl(intconname, strain))
    intcons_remove <- intcons %>%
      group_by(compound, concentration, comp_conc, plate_name, row, column) %>%
      mutate(eithermissing = any(count == 0)) %>%
      ungroup() %>%
      filter(eithermissing == TRUE) %>%
      select(compound, concentration, comp_conc, plate_name, row, column) %>%
      distinct()
    print(paste("Num wells removed due to no spike-in counts:", dim(intcons_remove)[1]))
    print(intcons_remove)
    filename1 = paste(gsub(".rds", "", savefilename), "_WellsNoSpikeIns_", format(Sys.time(), "%m_%d_%Y_%s"), ".rds", sep = "")
    saveRDS(intcons_remove, file = filename1)
    print(paste("intcons_remove saved here:", filename1, sep = " "))
    intcons_remove <- mutate(intcons_remove, platename_rowcol = paste(plate_name, row, column, sep = ":"))
    rawcounts_subset <- mutate(rawcounts_subset, platename_rowcol = paste(plate_name, row, column, sep = ":"))
    rawcounts_subset <- dplyr::filter(rawcounts_subset, !(platename_rowcol %in% intcons_remove$platename_rowcol))
    rawcounts_subset <- select(rawcounts_subset, -platename_rowcol)
  }else{
    print("No spike-in control name specified")
  }
  
  #Remove strains where all untreated counts are < num of untreated wells in entire screen
  untreated_raw <- dplyr::filter(rawcounts_subset, compound == untreatedname & !(strain %in% c("PCR-control", "Lysis-control", "intcon1", "intcon2")))
  untreated_raw <- dplyr::filter(untreated_raw, !grepl(intconname, strain))
  untreated_sumstrain <- untreated_raw %>%
    group_by(strain) %>%
    summarise(sumcount = sum(count, na.rm = T), numwells = n()) %>%
    ungroup() %>%
    filter(sumcount < numwells) #some arbitrarily small number, e.g. # of untreated wells
  print(paste("Removing strains with low counts (<", unique(untreated_sumstrain$numwells), "total) in all untreated: "))
  if(dim(untreated_sumstrain)[1] > 0){
    print(paste(unique(untreated_sumstrain$strain), collapse = "; "))
  }
  rawcounts_subset <- dplyr::filter(rawcounts_subset, !(strain %in% untreated_sumstrain$strain))
  untreated_raw <- dplyr::filter(untreated_raw, !(strain %in% untreated_sumstrain$strain))
  
  #Remove wells where all the counts are sum < lowwellcount (after removing strains with low count)--added 9/17/19
  print(paste0("Removing wells with low well count (<=", lowwellcount,")"))
  rawcounts_subset <- rawcounts_subset %>%
    group_by(compound, concentration, comp_conc, plate_name, row, column) %>%
    mutate(wellcount = sum(count, na.rm = T)) %>% #includes spike-in controls
    ungroup()
  
  #add wellcount removing spike-in controls
  wellcount_nointcon = rawcounts_subset %>%
    filter(!grepl(intconname, strain)) %>%
    group_by(compound, concentration, comp_conc, plate_name, row, column) %>%
    summarise(wellcount_remove_intcon = sum(count, na.rm = T)) %>%
    ungroup()
  rawcounts_subset = left_join(rawcounts_subset, wellcount_nointcon, by = c("compound", "concentration", "comp_conc", "plate_name", "row", "column"))
  
  lowwellcount_subset = filter(rawcounts_subset, wellcount <= lowwellcount)
  rawcounts_subset = filter(rawcounts_subset, wellcount > lowwellcount) #some what arbitrary small number, remove all wells that have low well counts (untreated or has compound) #check if throws error no replicates
  
  #Remove additional untreated wells where counts are < low_untreated_count (after removing intcon and removed strains)
  lowwellcount_untreated = filter(rawcounts_subset, compound == untreatedname) %>%
    filter(wellcount_remove_intcon <= low_untreated_count)
  rawcounts_subset = filter(rawcounts_subset, !((compound == untreatedname) & (wellcount_remove_intcon <= low_untreated_count)))
  
  lowwellcount_subset = rbind(lowwellcount_subset, lowwellcount_untreated)
  
  #Save which wells were removed due to low total counts
  #Add DMSO wells with low count (different total well count threshold) to the "removed" list
  print(paste0("Number of unique comps removed due to low total well count: ", length(unique(lowwellcount_subset$compound))))
  filename_lowcount =  paste(gsub(".rds", "", savefilename), "_LowWellCount_", format(Sys.time(), "%m_%d_%Y_%s"), ".rds", sep = "")
  saveRDS(lowwellcount_subset, file = filename_lowcount)
  remove(lowwellcount_subset)
  rawcounts_subset = select(rawcounts_subset, -wellcount, -wellcount_remove_intcon)

  #Remove plates with no untreated information (could have been removed due to low counts or no intcon)
  #Since we removed untreated wells with low counts, now make the criteria to have 3 or more "good' wells per plate, otherwise remove.
  noUntreatedonPlate = rawcounts_subset %>%
    select(plate_name, compound, concentration, row, column) %>%
    distinct() %>%
    group_by(plate_name) %>%
    summarise(countuntreated = sum(compound == untreatedname)) %>%
    ungroup() %>%
    filter(countuntreated < 3)
  
  print("Remove plates with < 3 untreated wells:")
  print(paste(noUntreatedonPlate$plate_name, collapse = "; "))
  
  rawcounts_subset <- filter(rawcounts_subset, !(plate_name %in% noUntreatedonPlate$plate_name))
  remove(noUntreatedonPlate)
  
  
  #Throw a warning if a strain always has zero counts in untreated wells on a given plate
  temp <- rawcounts_subset %>%
    dplyr::filter(compound == untreatedname) %>%
    group_by(strain, plate_name) %>%
    summarise(sumcount = sum(count)) %>%
    ungroup()
  if(any(temp$sumcount == 0)){
    print("WARNING: The following strain(s) has zero counts in all untreated wells in same plate")
    nountreated <- unique(select(dplyr::filter(temp, sumcount==0), plate_name, strain))
    print(paste(unique(nountreated$strain), collapse = "; "))
    saveRDS(nountreated, file = paste(gsub(".rds", "", savefilename), "_StrainPlate_nountreated", format(Sys.time(), "%m_%d_%Y_%s"), ".rds", sep = ""))
    print(nountreated)
  }
  
  
  #Moved replicate_id assignment to after removing all anomolous wells
  #Also Changed replicate_id to be unique for each well (rather than the same for those on the same plate)
  #Identify replicates--based on plate, row, column, compound, and concentration
  #Replicate wells on same plate will be averaged when calculating replicate correlation between plates
  repid <- rawcounts_subset %>%
    dplyr::select(compound, concentration, plate_name, row, column) %>%
    dplyr::distinct() %>%
    group_by(compound, concentration) %>%
    mutate(rep = paste("rep",row_number(), sep = "_")) %>%
    ungroup()
  
  rawcounts_subset <- left_join(rawcounts_subset, repid, by = c("compound", "concentration", "plate_name", "row", "column"))
  
  
  #Adding a pseudocount of 1 
  print ("Adding pseudocount of 1")
  rawcounts_subset$count <- rawcounts_subset$count + 1
  #rawcounts_subset$count <- replace(rawcounts_subset$count, rawcounts_subset$count == 0, 1)
  
  #Remove DNA spike-ins ("intcon1", "intcon2", "PCR-control", "Lysis-control") just for summing of the wells. Then return them. 
  print("Summing counts for each well (without internal controls)")
  sumwells <- rawcounts_subset %>%
    dplyr::filter(!grepl("intcon", strain) & !grepl("control", strain) & !grepl(intconname, strain)) %>%
    group_by(plate_name, row, column) %>%
    summarise(wellcount = sum(count)) %>%
    ungroup()
  
  #Add back the internal controls, and fill them with well counts for the same wells
  rawcounts_subset <- left_join(rawcounts_subset, sumwells, by = c("plate_name", "row", "column"))
  remove(sumwells)
  
  #Calculate log fraction for each strain
  print ("Calculating log fraction for each strain")
  rawcounts_subset <- rawcounts_subset %>%
    group_by(strain, plate_name, row, column) %>% 
    mutate(lf = log10(count/wellcount)) %>%
    ungroup()
  
  
  #Extract untreated wells for normalization
  rawcounts_untreated <- dplyr::filter(rawcounts_subset, compound == untreatedname)
  
  #Normalization and batch correction (per plate)
  #Either find averages for each strain for each plate, or use a linear regression 
  #untreated_lm <- lm(formula = lf ~ plate_name + 0, data = rawcounts_untreated)
  #use predict.lm on all the other treated data
  print("Batch correction and standardization")
  avg_untreateddf <- rawcounts_untreated %>%
    group_by(strain, plate_name) %>%
    summarise(avg_untreated = mean(lf), sd_untreated = sd(lf), avg_untreatedcount = mean(count), med_untreatedcount = median(count), sd_untreatedcount = sd(count)) %>%
    ungroup()
  
  rawcounts_subset <- left_join(rawcounts_subset, avg_untreateddf, by = c("strain", "plate_name"))
  rawcounts_subset <- rawcounts_subset %>%
    dplyr::mutate(corrected_lf = lf - avg_untreated) %>%
    dplyr::mutate(std_lf = corrected_lf/sd_untreated)
  
  #Replace infinites with NAs
  rawcounts_subset$std_lf <- replace(rawcounts_subset$std_lf, which(is.nan(rawcounts_subset$std_lf) | is.infinite(rawcounts_subset$std_lf)), NA)
  
  #Replace std log fractions with NAs if avg_untreatedcount is 1, and sd_untreatedcount = 0 (if all the untreated wells are 0 raw counts on a plate)
  rawcounts_subset$std_lf <- replace(rawcounts_subset$std_lf, which(rawcounts_subset$avg_untreatedcount == 1 & rawcounts_subset$sd_untreatedcount == 0), NA)
  
  #log2FC
  if(medianLFC){
    rawcounts_subset <- mutate(rawcounts_subset, log2FC = log2(count/med_untreatedcount))
  }else{
    rawcounts_subset <- mutate(rawcounts_subset, log2FC = log2(count/avg_untreatedcount))
  }
  
  #Calculate p-values from left tail of normal distribution
  #Use a standard normal, since values were standardized by stdev of untreated
  print("Calculating p-values")
  rawcounts_subset <- rawcounts_subset %>%
    dplyr::mutate(pval = pnorm(std_lf, 0, 1, lower.tail = TRUE))
  
  #Remove PCR-control and Lysis-control from replicate correlations
  rawcounts_subset_nointcon <- dplyr::filter(rawcounts_subset, !grepl("intcon", strain)) #remove intcons for downstream calculations
  rawcounts_subset_nointcon <- dplyr::filter(rawcounts_subset_nointcon, !grepl(intconname, strain)) #remove intcons for downstream calculations
  rawcounts_subset_nointcon <- dplyr::filter(rawcounts_subset_nointcon, strain != "Lysis-control" & strain != "PCR-control")
  
  #Calculate correlations 
  print("Calculating replicate correlations (without internal controls)")
  tempcorr <- replicatecorrelation(screendata = rawcounts_subset_nointcon, metric = "std_lf", comp_conc_separator = comp_conc_separator)
  rawcounts_subset <- left_join(rawcounts_subset, tempcorr, by = c("comp_conc", "rep")) #rep accounts for plate_name
  
  #Calculate z-score of averaged z-scores (std_lf)
  print("Calculating zscore of averaged standardized log fractions (without internal controls)")
  tempzscore <- zscore_avgreplicate(screendata = rawcounts_subset_nointcon, metric = "std_lf", comp_conc_separator = comp_conc_separator)
  rawcounts_subset <- left_join(rawcounts_subset, tempzscore, by = c("comp_conc", "rep", "strain")) #intcons should be NA
  
  #Calculate individual zscore of zscores for each replicate
  print("Calculating zscore of individual standardized log fractions (without internal controls)")
  tempzscore2 <- zscore_individual(screendata = rawcounts_subset_nointcon, metric = "std_lf", comp_conc_separator = comp_conc_separator)
  rawcounts_subset <- left_join(rawcounts_subset, tempzscore2, by = c("comp_conc", "rep", "strain"))
  
  #1/15/20 calculate well count fractions
  untreated = filter(rawcounts_subset, compound == untreatedname)
  untreated_byplate = untreated %>% 
    select(compound, concentration, wellcount, plate_name) %>%
    group_by(plate_name) %>% 
    summarise(medianwellcount = median(wellcount, na.rm = T)) %>%
    ungroup()
  rawcounts_subset = left_join(rawcounts_subset, untreated_byplate, by = "plate_name")
  rawcounts_subset = mutate(rawcounts_subset, wellcountfrac = wellcount/medianwellcount)
  
  if(newSchema){
    convertSchema(rawcounts_subset, to = "new")
  }
  
  if(!is.na(savefilename)){
    saveRDS(rawcounts_subset, file = savefilename)
    if(saveCSV){
      write.csv(rawcounts_subset, file = paste0(gsub(".rds", "", savefilename), ".csv"), row.names = F)
    }
  }else{
    return(rawcounts_subset)
  }
}

replicatecorrelation <- function(screendata, metric = "std_lf", comp_conc_separator = ":", randomize = FALSE, removemin = FALSE, savefilename = NA){
  #Calculate Pearson correlation between replicates.
  #Input
  #screendata: dataframe similar to screen3_adj, having compound, concentration, strain, rep
  #metric: the metric that correlations will be calculated from (std_lf is the standardized log fraction)
  #randomize: if true, will make a randomized column of rep_1 and find replicate correlations to assess null distribution
  #removemin: if true, will remove the minimum normalized residual and recalculate Pearson correlations,
  #will save the name of the strain that was removed
  #savefilename: if a file name is supplied, the function will save the dataframe into a file, else returns it
  #Finds replicate correlations between rep_n and rep_n+1
  #Output
  #Returns dataframe with compound_concentrations, correlation values or saves it as an rds. Will use this to then left_join to original data frame
  # Returns NA for correlation if there is no replicate
  
  #Add rounded concentration column "conc_round" if it doesn't exist
  if(!("conc_round" %in% colnames(screendata)))
  {
    screendata <- mutate(screendata, conc_round = signif(concentration, digits=1))
  }
  
  screendata <- dplyr::mutate(screendata, comp_conc = paste(compound, conc_round, sep= comp_conc_separator))
  
  #Rename metric column to target_metric
  colnames(screendata)[colnames(screendata)== metric] <- "target_metric"
  
  temp <- separate(screendata, col = rep, into = c("repname", "repnum"), sep = "_")
  temp$repnum <- as.numeric(temp$repnum)
  
  corlist = list()
  idx = 1
  startidx = min(temp$repnum)
  endidx = max(temp$repnum)
  if(startidx == endidx){
    combinerep <- select(screendata, comp_conc, rep)
    combinerep <- distinct(combinerep)
    combinerep$correlation <- NA
  }else{
    
    if(((endidx - startidx) %% 2) == 0){ 
      print("Odd number of replicates")
      if(endidx > startidx){
        endidx <- endidx - 1
      }
    }
    for(numrep in seq(startidx, endidx, 2)){ 
      name_rep1 = paste("rep", numrep, sep = "_")
      name_rep2 = paste("rep", numrep+1, sep = "_")
      subset <- screendata %>%
        dplyr::filter(rep == name_rep1  | rep == name_rep2 ) %>%
        group_by(comp_conc, strain, rep) %>% 
        summarise(avg_stdlf = mean(target_metric)) %>% #average replicate wells over plate if needed
        ungroup() %>%
        group_by(comp_conc) %>%
        spread(key = rep, value = avg_stdlf, drop = FALSE) %>%
        ungroup()
      subset <- as.data.frame(subset)
      
      #Remove rows with NAs
      if(sum(!complete.cases(subset)) > 0){
        print("replicatecorrelation- Removed strains with NA values:")
        print(unique(select(subset[!complete.cases(subset),], comp_conc, strain)))
        subset <- subset[complete.cases(subset), ]
      }
      
      colnames(subset)[colnames(subset) == name_rep1] <- "namerep1"
      colnames(subset)[colnames(subset) == name_rep2] <- "namerep2"
      
      corlist[[idx]] <- subset %>%
        group_by(comp_conc) %>%
        summarise(correlation = cor(namerep1, namerep2)) %>%
        mutate(repnames = paste(name_rep1, name_rep2, sep = "&"))
      
      idx = idx + 1
    }
    
    combinerep = as.data.frame(bind_rows(corlist))
    combinerep <- separate(combinerep, repnames, into= c("rep1", "rep2"), sep = "&" )
    combinerep <- tidyr::gather(combinerep, key = "repnum", value = "rep", rep1, rep2)
    combinerep <- select(combinerep, comp_conc, correlation, rep)
  }
  
  
  if(!is.na(savefilename)){
    saveRDS(combinerep, file = savefilename)
  }else{
    return(combinerep)
  }
  
}

zscore_avgreplicate <- function(screendata, metric = "std_lf", comp_conc_separator = ":"){
  #Averages of "metric" are taken between 2 replicates
  #The z-scores of these averaged metrics are calculated
  #and returns data frame
  
  #Add rounded concentration column "conc_round" if it doesn't exist
  if(!("conc_round" %in% colnames(screendata)))
  {
    screendata <- mutate(screendata, conc_round = signif(concentration, digits=1))
  }
  
  screendata <- dplyr::mutate(screendata, comp_conc = paste(compound, conc_round, sep=comp_conc_separator))
  
  #Rename metric column to target_metric
  colnames(screendata)[colnames(screendata)== metric] <- "target_metric"
  
  temp <- separate(screendata, col = rep, into = c("repname", "repnum"), sep = "_")
  temp$repnum <- as.numeric(temp$repnum)
  
  zscorelist = list()
  idx = 1
  for(numrep in seq(min(temp$repnum), max(temp$repnum), 2)){ 
    name_rep1 = paste("rep", numrep, sep = "_")
    name_rep2 = paste("rep", numrep+1, sep = "_")
    subset <- screendata %>%
      dplyr::filter(rep == name_rep1  | rep == name_rep2 ) %>%
      group_by(comp_conc, strain, rep) %>% 
      summarise(stdlf2 = mean(target_metric)) %>% #average replicate wells over plate
      ungroup()
    
    #Remove rows with NAs
    if(sum(!complete.cases(subset)) > 0){
      print("zscore_avgreplicate- Some strains may not have replicate information:")
      print("zscore_avgreplicate- Removed strains with NA values:")
      print(unique(select(subset[!complete.cases(subset),], comp_conc, strain)))
      subset <- subset[complete.cases(subset), ]
    }
    
    subset <- subset %>%
      group_by(comp_conc, strain) %>%
      mutate(avg_stdlf = mean(stdlf2)) %>% #average stdlf over 2 replicates
      ungroup() %>% 
      group_by(comp_conc) %>%
      mutate(zscore_stdlf = (avg_stdlf - mean(avg_stdlf))/sd(avg_stdlf)) %>%
      ungroup()
    subset <- as.data.frame(subset)
    
    zscorelist[[idx]] <- subset 
    
    idx = idx + 1
  }
  
  combinezscore = as.data.frame(bind_rows(zscorelist))
  combinezscore <- select(combinezscore, comp_conc, strain, zscore_stdlf, rep)
  return(combinezscore)
  
}

zscore_individual <- function(screendata, metric = "std_lf", comp_conc_separator = ":"){
  #Individual units of "metric" are taken for each replicate
  #The z-scores between strains for a given compound are calculated
  #and returns data frame
  
  #Add rounded concentration column "conc_round" if it doesn't exist
  if(!("conc_round" %in% colnames(screendata)))
  {
    screendata <- mutate(screendata, conc_round = signif(concentration, digits=1))
  }
  
  screendata <- dplyr::mutate(screendata, comp_conc = paste(compound, conc_round, sep=comp_conc_separator))
  
  #Rename metric column to target_metric
  colnames(screendata)[colnames(screendata)== metric] <- "target_metric"
  
  temp <- separate(screendata, col = rep, into = c("repname", "repnum"), sep = "_")
  temp$repnum <- as.numeric(temp$repnum)
  
  zscorelist = list()
  idx = 1
  for(numrep in seq(min(temp$repnum), max(temp$repnum), 1)){ 
    name_rep1 = paste("rep", numrep, sep = "_")
    subset <- screendata %>%
      dplyr::filter(rep == name_rep1) %>%
      group_by(comp_conc, strain, rep) %>% 
      summarise(avgmetric = mean(target_metric)) %>% #average replicate wells within plate
      ungroup()
    
    #Remove rows with NAs
    if(sum(!complete.cases(subset)) > 0){
      print("zscore_individual_replicate- Some strains may not have replicate information:")
      print("zscore_individual_replicate- Removed strains with NA values:")
      print(unique(select(subset[!complete.cases(subset),], comp_conc, strain)))
      subset <- subset[complete.cases(subset), ]
    }
    
    subset <- subset %>%
      group_by(comp_conc) %>%
      mutate(zscore_stdlf2 = (avgmetric - mean(avgmetric))/sd(avgmetric)) %>%
      ungroup()
    subset <- as.data.frame(subset)
    
    zscorelist[[idx]] <- subset 
    
    idx = idx + 1
  }
  
  combinezscore = as.data.frame(bind_rows(zscorelist))
  combinezscore <- select(combinezscore, comp_conc, strain, zscore_stdlf2, rep)
  return(combinezscore)
  
}

convertSchema <- function(screendata, to = "new"){
  #Converts column labels between old and new schema
  #Input data frame
  #Output same data frame but with column names converted to new schema (to = "new") (or back to original to = "old")
  if(to == "new"){
    if("compound" %in% colnames(screendata)){
      print("Converting to new schema")
      colnames(screendata)[colnames(screendata) == "compound"] = "broad_id"
      colnames(screendata)[colnames(screendata) == "concentration"] = "pert_dose_act"
      colnames(screendata)[colnames(screendata) == "row"] = "well_row_id"
      colnames(screendata)[colnames(screendata) == "column"] = "well_col_id"
      colnames(screendata)[colnames(screendata) == "plate_name"] = "plate_id"
      colnames(screendata)[colnames(screendata) == "strain"] = "strain_name"
      colnames(screendata)[colnames(screendata) == "rep"] = "replicate_id"
      colnames(screendata)[colnames(screendata) == "conc_round"] = "x_conc_round"
    }
  }
  
  if(to == "old"){
    if("broad_id" %in% colnames(screendata)){
      print("Converting to old naming convention")
      colnames(screendata)[colnames(screendata) == "broad_id"] = "compound"
      colnames(screendata)[colnames(screendata) == "pert_dose_act"] = "concentration"
      colnames(screendata)[colnames(screendata) == "well_id"] = "well"
      colnames(screendata)[colnames(screendata) == "well_row_id"] = "row"
      colnames(screendata)[colnames(screendata) == "well_col_id"] = "column"
      colnames(screendata)[colnames(screendata) == "plate_id"] = "plate_name"
      colnames(screendata)[colnames(screendata) == "strain_name"] = "strain" #was strain_id?
      colnames(screendata)[colnames(screendata) == "replicate_id"] = "rep"
      colnames(screendata)[colnames(screendata) == "x_conc_round"] = "conc_round"
    }
  }
}

subsetSLF <- function(screendata_path, keep_columns = c("strain", "compound", "concentration", "plate_name", "row", "column", "count", "rep", "wellcount", "wellcountfrac", "std_lf", "zscore_stdlf", "zscore_stdlf2", "correlation"), savefilename = NA){
  #Takes .rds file as input, which was output by compScreenPipeline
  #Returns a subset of the columns defined by keep_columns
  #Saves .csv file in savefilename, otherwise uses the screendatapath name
  
  screendata = readRDS(screendata_path)
  convertSchema(screendata, to = "old")
  screendata_sub = select(screendata, one_of(keep_columns))
  if(!is.na(savefilename)){
    write.csv(screendata_sub, savefilename, row.names = F)
  }else{
    savefilename = paste0(gsub(".rds", "", screendata_path), "_subset.csv")
    write.csv(screendata_sub, savefilename, row.names = F)
  }
  
}

