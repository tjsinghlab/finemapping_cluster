import pandas as pd
import os
import glob
from openai import OpenAI

from cols import Cols
from preprocessing import Preprocess

def main():

    directory_of_sumstats = '/gpfs/commons/groups/sanjana_lab/mdrabkin/gwas_data/raw/European' # YOUR INPUT PATH
    output_directory = '/gpfs/commons/groups/sanjana_lab/mdrabkin/gwas_data/preprocessed/European/v2' # YOUR OUTPUT PATH

    ft = Preprocess(
        client = OpenAI(
            organization = os.getenv('OPENAI_ORGANIZATION'), # see README.md
            project = os.getenv('OPENAI_PROJECT'),
            api_key = os.getenv('OPENAI_API_KEY')
        ), 
        out_dir=output_directory # sets output path
    )

    # optional parameters
    ft.significance_threshold = 5e-8 # default value is already 5e-8
    ft.ancestry = 'EUR' # default value is 'EUR' but can be statically

    for path in glob.glob(directory_of_sumstats + '/*'):

        # loop through all files in input path
        print(path)

        # run loadmap_sumstats_table
        res = ft.loadmap_sumstats_table(
            path,
            verbose=True # set to True if you want it to print what it's doing
        ) # returns 0 upon success

        if res != 0:
            continue

        # create leadsnp table if GPT was able to detect, map, and rearrange columns
        ft.create_leadsnp_table(verbose=True)


if __name__ == '__main__':
    cols_instance = Cols('header.yaml')
    print(cols_instance)

    main()
