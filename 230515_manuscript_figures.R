######Manuscript figures######
#January 17, 2023

library(dplyr)
library(ggplot2)
library(tidyr)

######Function definition######
NPI_vs_binned_SLF = function(screen_df, collapse_reps = "max", SLF_colname = "std_lf", NPI_colname = "NPI_robust_median", exclude_controls = T, PSA_only = T, binwidth1 = 0.02, NPI_min = -999, SLF_bin_min = -8){
  #input: screen_df is a data frame with SLF info and NPI info
  #collapse_reps allows you to choose how to collapse SLF replicates, by "max" (most conservative), "min", or "mean"
  #Assumes compounds that need to be excluded were already excluded, and any strains that should be excluded were also already excluded
  #Print table summary
  #Print both boxplot version and histogram version
  #NPI_min defines the lowest NPI (whatever column defined by NPI_colname) to plot. All NPIs < NPI_min get converted to NPI_min. Set to -999 for it to have no effect
  #SLF_bin_min defines the lowest SLF bin. It always goes up to -3.
  
  
  #Find appropriate columns
  colnames(screen_df)[colnames(screen_df) == NPI_colname] = "NPI_forplot"
  colnames(screen_df)[colnames(screen_df) == SLF_colname] = "SLF_forplot"
  
  #Exclude control compounds if desired
  if(exclude_controls){
    neg_control = "DMSO"
    pos_control = "BRD-K04804440" #CIP
    trimethoprim = "BRD-K07208025" #Trimethoprim
    intconname = "PCR"
    
    compexclude = c(neg_control, pos_control, trimethoprim)
    screen_df = filter(screen_df, !(grepl(paste(compexclude, collapse = "|"), compound)))
  }
  
  if(PSA_only){
    nonpsa = c("Acine", "Ecoli", "Acine", "Klebs", "Staph")
    screen_df = filter(screen_df, !(strain %in% nonpsa))
  }

  #Collapse SLF to max SLF per treatment
  if(collapse_reps == "max"){
    screen_df2 <- screen_df %>%
      group_by(strain, compound, concentration) %>%
      dplyr::slice(which.max(SLF_forplot)) %>%
      ungroup() 
  }else if(collapse_reps == "min"){
    screen_df2 <- screen_df %>%
      group_by(strain, compound, concentration) %>%
      dplyr::slice(which.min(SLF_forplot)) %>%
      ungroup() 
  }else if(collapse_reps == "mean"){
    screen_df2 <- screen_df %>%
      group_by(strain, compound, concentration) %>%
      summarise(SLF_forplot = mean(SLF_forplot, na.rm = T), NPI_forplot = unique(NPI_forplot)
      ) %>%
      ungroup() 
  }

  ######Convert lowest NPI to NPI_min###################
  screen_df2 = mutate(screen_df2, NPI_forplot = ifelse(NPI_forplot < NPI_min, NPI_min, NPI_forplot))
  
  ######PLOT: BINNED SLF vs NPI###########################
  complist2 = filter(screen_df2, SLF_forplot <= -3) %>% arrange(SLF_forplot)
  complist2 = mutate(complist2, comp_strain = paste(compound, strain, sep = ":"))
  length(unique(complist2$compound))
  
  ####Try negative SLF labels 
  SLF_bins = seq(from = SLF_bin_min, to = -3, by = 1)
  screen_df2_filt = screen_df2 %>%
    filter(!is.na(NPI_forplot))
  screen_df2_filt$threshold = '0'
  
  #Lowest bin
  screen_df2_filt = mutate(screen_df2_filt, threshold = ifelse(SLF_forplot < SLF_bins[1], paste0("< ", SLF_bins[1]), threshold))
  for(binnum in 2:length(SLF_bins)){
    screen_df2_filt = mutate(screen_df2_filt, threshold = ifelse(threshold == '0' & (SLF_forplot < SLF_bins[binnum]), as.character(SLF_bins[binnum]), threshold))
  }
  screen_df2_filt = mutate(screen_df2_filt, threshold = ifelse(threshold == '0' & (SLF_forplot >= SLF_bins[length(SLF_bins)]), paste0(">= ", as.character(SLF_bins[length(SLF_bins)])), threshold))
  screen_df2_filt$threshold = factor(screen_df2_filt$threshold, levels = c(paste0('< ', SLF_bins[1]), as.character(SLF_bins[2:length(SLF_bins)]), paste0(">= ", SLF_bins[length(SLF_bins)])))
  
  ####Summary of counts in all the bins
  screen_df2_filt_summary = screen_df2_filt %>%
    group_by(threshold, strain) %>%
    mutate(compsperstrain = length(unique(compound)), strain_numcomp = paste0(unique(strain), "(", unique(compsperstrain), ")")) %>%
    ungroup() %>%
    distinct() %>%
    group_by(threshold) %>%
    arrange(desc(compsperstrain)) %>%
    summarise(numcomps_strain = n(), numcomps = length(unique(compound)), numstrains = length(unique(strain)), strains = paste(unique(strain_numcomp), collapse = "; ")) %>%
    ungroup()
  
  View(screen_df2_filt_summary)
  
  p2 = ggplot(screen_df2_filt, aes(x = threshold, y = NPI_forplot * 100))+
    geom_boxplot()+
    geom_hline(yintercept = 0, color = "gray", linetype = "dashed")+
    geom_hline(yintercept = 100, color = "gray", linetype = "dashed")+
    # labs(x = "SLF Bin", y = NPI_colname)+
    scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100), limits = c(-20, 120))+
    labs(x = "SLF Bin", y = "% Growth Inhibition")+
    theme_bw(base_size = 16)
  print(p2)
  
  

  #histograms, subsample >= -3 bin
  temp = filter(screen_df2_filt, threshold == ">= -3")
  temp_sub = sample_n(temp, 100)
  # temp_sub$threshold = ">= -3"
  screen_df2_filt_subsample = rbind(temp_sub, filter(screen_df2_filt, !(threshold %in% c(">= -3"))))
  
  p4 = ggplot(screen_df2_filt_subsample, aes(x = NPI_forplot*100))+
    geom_histogram(binwidth = binwidth1)+
    geom_vline(xintercept = 0, color = "gray", linetype = "dashed")+
    geom_vline(xintercept = 100, color = "gray", linetype = "dashed")+
    labs(x = "% Growth Inhibition", y = "# Compound-Strain")+
    scale_y_continuous(breaks = c(0, 5, 10))+
    scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100), limits = c(-20, 120))+
    # xlim(-1.15, 1.15)+
    coord_flip()+
    facet_wrap(~threshold, nrow = 1)+
    theme_bw(base_size = 16)
  print(p4)
  
  plotlist = list()
  plotlist[[1]] = p2
  plotlist[[2]] = p4
  return(plotlist)
  
}

######Load Data############
# screen_hits_demultiplex = readRDS('/idi/cgtb/jbagnall/Brad_pseudomonas/screen1/screen1_processed_NPI_filt_201027.rds')
screen_hits_demultiplex = readRDS('/idi/cgtb/jbagnall/Brad_pseudomonas/manuscript_figures/screen1_slf_npi_230515.rds')


######Plot NPI vs SLF########
plots = NPI_vs_binned_SLF(screen_hits_demultiplex, collapse_reps = "max", NPI_colname = "NPI_robust_median", exclude_controls = T, PSA_only = T, binwidth = 2, NPI_min = -0.2, SLF_bin_min = -7)
# plots2 = NPI_vs_binned_SLF(screen_hits_demultiplex2, collapse_reps = "max", NPI_colname = "NPI_robust_median", binwidth = 2, NPI_min = -0.2, SLF_bin_min = -7)
# print(plots[[1]])
# print(plots[[2]])


#########BRD1401 SLF vs demux########
brd1401= filter(screen_hits_demultiplex, grepl("K47991401", compound))
brd1401_collapse = brd1401 %>%
  group_by(strain, compound) %>%
  summarise(std_lf = max(std_lf, na.rm = T), NPI_robust_median = unique(NPI_robust_median)) %>%
  ungroup()

p0 = ggplot(brd1401_collapse, aes(x = std_lf, y = NPI_robust_median*100))+
  geom_point() +
  ggrepel::geom_text_repel(data = filter(brd1401_collapse, strain == "oprL"), aes(x = std_lf, y = NPI_robust_median *100, label = strain))+
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed")+
  geom_hline(yintercept = 100, color = "gray", linetype = "dashed")+
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed")+
  scale_x_continuous(breaks = c(-10, -5, 0, 5), limits=c(-10,5))+
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits=c(-5, 100))+
  labs(x = "SLF", y = "% Growth Inhibition") +
  theme_bw(base_size = 16)+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"))




