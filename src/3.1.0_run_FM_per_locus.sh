#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --output=log/%x.%j.out

cd /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune

module purge
module load R/4.2.2 

sumstats_name=$1
path_df_sumstats=$2
ld_pop=$3
locus=$4
N_tot=$5
N_cases=$6
window_mb=$7

Rscript src/3.1.0_run_FM_per_locus.sh $1 $2 $3 $4 $5 $6 ${7} 

