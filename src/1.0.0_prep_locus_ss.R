#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
setwd("/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_cluster")
source("src/finemapping_functions.R")

library(readr)
library(data.table)
library(dplyr)
library(Matrix)

#read in arguments
path_leadSNP <- as.character(args[1]) #locus needs to be saved as "chr.pos"
sumstats_name <- as.character(args[2])
path_df_sumstats <- as.character(args[3])
ld_pop <- as.character(args[4])
window_mb <- as.numeric(args[5])

df.leadSNP <- fread(path_leadSNP, colClasses = c("numeric","numeric", "character"))

window_bp <- window_mb*1000000
df.leadSNP <- df.leadSNP %>%
    mutate(LOWER = ifelse((BP-window_bp) < 0, 1, BP-window_bp), #if lead snp is too close to chr beginning
           UPPER = BP+window_bp)
  
dir.create(paste0("output/", sumstats_name))
dir.create(paste0("output/", sumstats_name, "/", ld_pop, "_", window_mb, "Mb"))
dir.create(paste0("output/", sumstats_name, "/", ld_pop, "_", window_mb, "Mb/ss"))
dir.create(paste0("output/", sumstats_name, "/", ld_pop, "_", window_mb, "Mb/ld"))

ss <- fread(path_df_sumstats, na.strings=c(""," ","NA"))

for (i in 1:nrow(df.leadSNP)) {
  LOCUS <- as.character(df.leadSNP$locus[i])
  CHR <- df.leadSNP$CHR[i]
  START <- df.leadSNP$LOWER[i]
  END <- df.leadSNP$UPPER[i]
  save_ss_per_locus(sumstats_name, ss, LOCUS, CHR, START, END, window_mb)
}

