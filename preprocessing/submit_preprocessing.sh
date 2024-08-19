#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --output=./logs/preprocess_%A_%a.out
#SBATCH --array=1-68
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G

# details in README.md
export OPENAI_ORGANIZATION="org-7Tz4ejgn3Rw1nC3LMSjTcdVN"
export OPENAI_PROJECT="proj_JwOi33NfPxndlxwTCQL9UIiK"
export OPENAI_API_KEY="sk-None-DvEi4YfrZBudVnrZK1nTT3BlbkFJQsggOLv29hV8bq6Sm8Mf"

base_dir="/gpfs/commons/groups/sanjana_lab/mdrabkin/gwas_data"
export PREPROC_INPUT_DIR="${base_data}/raw/European"
export PREPROC_OUTPUT_DIR="${base_data}/preprocessed/European/v5"
export PREPROC_FILEPATHS="./filepaths.txt"

module load python
source .venv/bin/activate
python3 batch.py ${PREPROC_INPUT_DIR} ${PREPROC_FILEPATHS} ${PREPROC_OUTPUT_DIR}
