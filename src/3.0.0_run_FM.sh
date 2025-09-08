#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --output=log/%x.%j.tmp.out
#SBATCH --error=log/%x.%j.tmp.err


# Function to initialize Conda and activate an environment
initialize_conda() {
  # Locate the conda executable using `which conda`
  local conda_exe=$(which conda)

  # Check if conda was found
  if [ -z "$conda_exe" ]; then
    # If not found, print an error message and exit with status 1
    echo "Conda not found in PATH."
    exit 1
  fi

  # Determine the base directory of the Conda installation using `conda info --base`
  local conda_base=$($conda_exe info --base)

  # Check if the conda.sh initialization script exists in the determined base directory
  if [ -f "$conda_base/etc/profile.d/conda.sh" ]; then
    # If it exists, source the script to initialize Conda
    source "$conda_base/etc/profile.d/conda.sh"
  else
    # If not found, print an error message and exit with status 1
    echo "Conda initialization script not found."
    exit 1
  fi

  # Activate the specified Conda environment
  conda activate hail_env
}

# Call the function to initialize Conda and activate the environment
initialize_conda

#module purge
#module load R/4.2.2 
#module load gcloud

cd /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_cluster

sumstats_name=$1
ld_pop=$2
N_tot=$3
N_cases=$4
window_mb=$5

path_leadSNP="/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_cluster/data/${sumstats_name}_leadSNPs.tsv"

ld_files="gs://nygc-comp-d-95c4-tf/finemapping_analysis/autoimmune/data/${sumstats_name}/LD/*${window_mb}*.bgz"
ss_files="gs://nygc-comp-d-95c4-tf/finemapping_analysis/autoimmune/data/${sumstats_name}/LD/*${window_mb}*.tsv"
nygc_dir="/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_cluster/output/${sumstats_name}/${ld_pop}_${window_mb}Mb/ld"

gsutil -m cp ${ld_files} ${nygc_dir}
gsutil -m cp ${ss_files} ${nygc_dir}

loci=($(awk '{print $3}' ${path_leadSNP}))

for locus in "${loci[@]}"; do
  echo $locus
  sbatch -J "${sumstats_name}_${locus}" src/3.1.0_run_FM_per_locus.sh $sumstats_name $ld_pop $locus $N_tot $N_cases ${window_mb}
done

mv log/%x.%j.tmp.out log/${sumstats_name}_%j.out
mv log/%x.%j.tmp.err log/${sumstats_name}_%j.err
