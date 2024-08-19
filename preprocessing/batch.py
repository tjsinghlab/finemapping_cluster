import pandas as pd
import os
import glob
from openai import OpenAI
import multiprocessing
import argparse
import importlib

import preprocessing
import cols
importlib.reload(cols)
importlib.reload(preprocessing)
from cols import Cols
from preprocessing import Preprocess

def main(output_directory, file_path):
    client = OpenAI(
        organization=os.getenv('OPENAI_ORGANIZATION'),  # see README.md
        project=os.getenv('OPENAI_PROJECT'),
        api_key=os.getenv('OPENAI_API_KEY')
    )

    ft = Preprocess(client=client, out_dir=output_directory)
    res = ft.loadmap_sumstats_table(file_path, verbose=False)
    if res != 0:
        return
    ft.create_leadsnp_table(verbose=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process sumstats files in batch')
    parser.add_argument('input_directory', type=str, help='Path to the input directory')
    parser.add_argument('file_paths', type=str, help='Path to the file paths')
    parser.add_argument('output_directory', type=str, help='Path to the output directory')
    args = parser.parse_args()

    # cols_instance = Cols('header.yaml')
    # print(cols_instance)

    # Get the SLURM_ARRAY_TASK_ID from environment variables
    slurm_array_task_id = os.environ.get("SLURM_ARRAY_TASK_ID")
    if not slurm_array_task_id:
        raise ValueError("SLURM_ARRAY_TASK_ID environment variable is not set")
    with open(args.file_paths, "r") as f:
        file_path = f.readlines()[int(slurm_array_task_id) - 1].strip()

    main(args.output_directory, os.path.join(args.input_directory, file_path))