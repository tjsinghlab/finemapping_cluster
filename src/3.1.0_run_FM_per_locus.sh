#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --output=log/%x.%j.out

cd /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune

module purge
module load R/4.2.2 

sumstats_name=$1
ld_pop=$2
locus=$3
N_tot=$4
N_cases=$5
window_mb=$6

Rscript src/3.1.1_run_FM_per_locus.R $1 $2 $3 $4 $5 ${6}

