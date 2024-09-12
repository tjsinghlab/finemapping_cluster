#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --output=log/3.1.0_run.%j.tmp.out
#SBATCH --error=log/3.1.0_run.%j.tmp.err
#
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

# Call the function to initialize Conda and activate the environment
initialize_conda



#module purge
#module load R/4.2.2 

sumstats_name=$1
ld_pop=$2
locus=$3
N_tot=$4
N_cases=$5
window_mb=$6

echo "Running sumstats_name: $sumstats_name, locus: $locus, ld_pop: $ld_pop, window_mb: $window_mb"

Rscript src/3.1.1_run_FM_per_locus.R $1 $2 $3 $4 $5 ${6}


# Rename the log files to include sumstats_name and locus, preserving the job ID
mv log/3.1.0_run.${SLURM_JOB_ID}.tmp.out log/3.1.0_run_${sumstats_name}_${locus}_${SLURM_JOB_ID}.out
mv log/3.1.0_run.${SLURM_JOB_ID}.tmp.err log/3.1.0_run_${sumstats_name}_${locus}_${SLURM_JOB_ID}.err
