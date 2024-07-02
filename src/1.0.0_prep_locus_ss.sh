#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --output=log/%x.%j.out

cd /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune
module purge
module load R/4.2.2 

#path_leadSNP <- as.character(args[1]) #locus needs to be saved as "chr.pos"
#sumstats_name <- as.character(args[2])
#path_df_sumstats <- as.character(args[3])
#ld_pop <- as.character(args[4])
#window_mb <- as.numeric(args[5])
                 
Rscript src/1.0.0_get_locus_ss.R $1 $2 $3 $4 ${5} 

proj_dir="/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune"
ss_path="${proj_dir}/output/${sumstats_name}/${ld_pop}_${window_mb}Mb/ss/*.txt"

#copy locus-specific SS files to cloud
gsutil -m cp ${ss_path} gs://nygc-comp-d-95c4-tf/finemapping_analysis/autoimmune/data/${sumstats_name}/ss/

#copy lead SNP file to cloud
gsutil -m cp ${1} gs://nygc-comp-d-95c4-tf/finemapping_analysis/autoimmune/data/${sumstats_name}/