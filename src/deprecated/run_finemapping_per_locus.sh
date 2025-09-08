#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --output=log/%x.%j.out

cd /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_cluster
module purge
module load R/4.2.2 

# sumstats_name <- as.character(args[1])
# path_df_sumstats <- as.character(args[2])
# ld_pop <- as.character(args[3])
# LOCUS <- as.character(args[4])
# CHR <- as.numeric(args[5])
# START <- as.numeric(args[6])
# END <- as.numeric(args[7])
# N_tot <- as.numeric(args[8])
# N_cases <- as.numeric(args[9])
# window_mb <- as.numeric(args[10])


Rscript src/run_finemapping_per_locus.R $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}
