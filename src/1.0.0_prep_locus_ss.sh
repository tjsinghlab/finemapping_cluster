#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
##SBATCH --output=log/%x.%j.tmp.out


#to submit:
##sbatch 1.0.0_prep_locus_ss.sh $path_to_lead_SNP $sumstats_name $path_to_sumstats $ld_pop $window_mb
cd /gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune

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

initialize_conda

#module purge
#module load R/4.2.2 
#module load gcloud

#path_leadSNP <- as.character(args[1]) #locus needs to be saved as "chr.pos"
#sumstats_name <- as.character(args[2])
#path_df_sumstats <- as.character(args[3])
#ld_pop <- as.character(args[4])
#window_mb <- as.numeric(args[5])
                 
Rscript src/1.0.0_prep_locus_ss.R $1 $2 $3 $4 ${5} 

sumstats_name=$2
ld_pop=$4
window_mb=$5

proj_dir="/gpfs/commons/groups/nygcfaculty/lappalainen_singh/finemapping_autoimmune"
ss_path="${proj_dir}/output/${sumstats_name}/${ld_pop}_${window_mb}Mb/ss/*.txt"

#copy locus-specific SS files to cloud
gsutil -m cp ${ss_path} gs://nygc-comp-d-95c4-tf/finemapping_analysis/autoimmune/data/${sumstats_name}/ss/

#copy lead SNP file to cloud
gsutil -m cp ${1} gs://nygc-comp-d-95c4-tf/finemapping_analysis/autoimmune/data/${sumstats_name}/

mv log/%x.%j.tmp.out log/${sumstats_name}_%j.out
