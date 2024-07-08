#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --output=log/%x.%j.out
module purge
module load R/4.2.2 
module load gcloud

cd /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune

path_leadSNP=$1
sumstats_name=$2
path_df_sumstats=$3
ld_pop=$4
window_mb=$5

ld_files="gs://nygc-comp-d-95c4-tf/finemapping_analysis/autoimmune/data/${sumstats_name}/LD/*${window_mb}*.bgz"
ss_files="gs://nygc-comp-d-95c4-tf/finemapping_analysis/autoimmune/data/${sumstats_name}/LD/*${window_mb}*.tsv"
nygc_dir="/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/output/${sumstats_name}/${ld_pop}_${window_mb}Mb/ld"

#gsutil -m cp ${ld_files} ${nygc_dir}
gsutil -m cp ${ss_files} ${nygc_dir}

loci=($(awk '{print $3}' ${path_leadSNP}))

for locus in "${loci[@]}"; do
  
  echo "Third column value: $locus"
done