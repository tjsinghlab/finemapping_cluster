#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
setwd("/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_cluster")
source("src/finemapping_functions.R")

library(susieR)
library(readr)
library(CARMA)
library(data.table)
library(dplyr)
library(Matrix)
library(tidyr)
library(stringr)

sumstats_name <- as.character(args[1])
ld_pop <- as.character(args[2])
locus <- as.character(args[3])
N_tot <- as.numeric(args[4])
N_cases <- as.numeric(args[5])
window_mb <- as.numeric(args[6])
  

dat <- bring_in_UKBB_data(sumstats_name, ld_pop, window_mb, locus)
ss_filtered <- dat[[1]]
ld_filtered <- dat[[2]]

run_susie(ss_filtered, ld_filtered, N_tot, N_cases,  sumstats_name, ld_pop, window_mb, locus, "UKBB")
run_CARMA(ss_filtered, ld_filtered, sumstats_name, ld_pop, window_mb, locus,  "UKBB") 






