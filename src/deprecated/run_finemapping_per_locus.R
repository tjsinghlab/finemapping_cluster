#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
setwd("/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_cluster")
source("src/finemapping_functions.R")

sumstats_name <- as.character(args[1])
path_df_sumstats <- as.character(args[2])
ld_pop <- as.character(args[3])
LOCUS <- as.character(args[4])
CHR <- as.numeric(args[5])
START <- as.numeric(args[6])
END <- as.numeric(args[7])
N_tot <- as.numeric(args[8])
N_cases <- as.numeric(args[9])
window_mb <- as.numeric(args[10])

library(susieR)
library(readr)
library(CARMA)
library(data.table)
library(dplyr)
library(Matrix)

ss_per_locus <- get_ss_per_locus(path_df_sumstats, 
                                 CHR,
                                 LOCUS,
                                 START, 
                                 END)

ld_per_locus_res <- get_ld_per_locus(ss_per_locus, LOCUS, CHR, START, END)
ss_filtered <- ld_per_locus_res[[1]]
ld_filtered <- ld_per_locus_res[[2]]

write_tsv(ss_filtered, paste0("output/", sumstats_name, "/", ld_pop, "_",
                              window_mb, "Mb_window/ss/", sumstats_name, "_ss_",
                              window_mb, "Mb_locus_", LOCUS, "_matched.z"))
write_tsv(ld_filtered, paste0("output/", sumstats_name, "/", ld_pop, "_",
                              window_mb, "Mb_window/ld/", sumstats_name, "_ss_",
                              window_mb, "Mb_locus_", LOCUS, "_matched.ld"))

run_susie(ss_filtered, ld_filtered, N_tot, N_cases,  sumstats_name, ld_pop, window_mb, LOCUS, "UKBB")
run_CARMA(ss_filtered, ld_filtered, sumstats_name, ld_pop, window_mb, LOCUS,  "UKBB") 






