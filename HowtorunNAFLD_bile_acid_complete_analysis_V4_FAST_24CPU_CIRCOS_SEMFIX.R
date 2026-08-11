if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
setwd("C:/Users/Rainjie/Desktop/NAFLD/NAFLD_Data/NAFLD_BASummary20260805/BA_CompleteAnalysis")
rm(list = ls())
gc()
source("NAFLD_bile_acid_complete_analysis_V4_FAST_24CPU_CIRCOS_SEM_MANUAL.R")

