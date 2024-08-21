#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --output=./logs/preprocess_%A_%a.out
#SBATCH --array=1-68
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G

# details in README.md
export OPENAI_ORGANIZATION=""
export OPENAI_PROJECT=""
export OPENAI_API_KEY=""
export GEMINI_API_KEY=""

export PREPROC_INPUT_DIR="/gpfs/commons/groups/sanjana_lab/mdrabkin/gwas_data/raw/European"
export PREPROC_OUTPUT_DIR="/gpfs/commons/groups/sanjana_lab/mdrabkin/gwas_data/preprocessed/European/v5"
export PREPROC_FILEPATHS="./filepaths.txt"

module load python
source .venv/bin/activate
python3 batch.py ${PREPROC_INPUT_DIR} ${PREPROC_FILEPATHS} ${PREPROC_OUTPUT_DIR}
