#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --output=log/%x.%j.out

module purge
module load R/4.2.2 
module load gcloud

cd /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune

sumstats_name=$1
ld_pop=$2
N_tot=$3
N_cases=$4
window_mb=$5

path_leadSNP="/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/data/${sumstats_name}_leadSNPs.tsv"

ld_files="gs://nygc-comp-d-95c4-tf/finemapping_analysis/autoimmune/data/${sumstats_name}/LD/*${window_mb}*.bgz"
ss_files="gs://nygc-comp-d-95c4-tf/finemapping_analysis/autoimmune/data/${sumstats_name}/LD/*${window_mb}*.tsv"
nygc_dir="/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune/output/${sumstats_name}/${ld_pop}_${window_mb}Mb/ld"

#gsutil -m cp ${ld_files} ${nygc_dir}
#gsutil -m cp ${ss_files} ${nygc_dir}

loci=($(awk '{print $3}' ${path_leadSNP}))

for locus in "${loci[@]}"; do
  echo $locus
  sbatch -J "${sumstats_name}_${locus}" src/3.1.0_run_FM_per_locus.sh $sumstats_name $ld_pop $locus $N_tot $N_cases ${window_mb}
done