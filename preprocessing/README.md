# Fine-Mapping Preprocessing

## Major caveat

tl;dr it doesn't always work and when it does, always check that the columns mapped correctly. **loadmap_sumstats_table** creates a log file with stated mappings, which should be verified by human inspection.

We cannot be certain from the columns alone which allele is the effect allele if it is not explicitly stated. For these cases, consult the README from the study or the paper itself. 

This preprocessing fails when chr & position are omitted. In such cases, rsIDs need to be mapped manually. 

This preprocessing fails if the PI didn't include enough of the necessary columns.

## Functionality

1. **ask_sumstats_col_order**: Compares the given file column names to the expected fields. This method is used for standardizing summary stat files by mapping column name and ordering to an expected format.
    - Determines which columns from sumstats file are needed for the finemapping pipeline using a 2-shot prompt.
    - Returns those column names in the order expected by the pipeline.
    - See caveat above.

2. **loadmap_sumstats_table**: Loads and maps the summary statistics file to the expected column order for a finemapping pipeline.
    - Loads in sumstats dataframe
        - expects one of the following formats: **.csv, .tsv, .txt, .csv.gz, .tsv.gz, .txt.gz**
    - Calls ask_sumstats_col_order to get column mappings.
    - Reorders and renames dataframe columns. Possibly excludes some columns.
    - Save the results to a TSV file, based on input file name and user set output directory.
    - Returns 1 on success and 0 on fail, most likely because the column mappings didn't fully work. (See caveat above.)
    
3. **create_leadsnp_table**: Creates a table of lead SNPs from the summary statistics file.
    - Filter to only significant SNPs.
    - Convert chromosome "X" to 23, "Y" to 24 and ensure numeric types.
    - Sort by chromosome and position.
    - Assign previous chromosome and position for comparison.
    - Assign locus numbers.
        - First assign truth values based on conditions:
            ```python
                (sumstats_df['chromosome'] != sumstats_df['prev_Chr']) | (abs(sumstats_df['position'] - sumstats_df['prev_Position']) > 5e5)
                # chromosome has changed or the prior SNP is > 5e5 BP away --> new locus reached
            ```
        - Then create running sum of truths -> +1 for every true indicates new locus.
    - Remove temporary columns.
    - Group by locus and filter for minimum p-value.
    - Remove duplicated loci.
    - Transform columns and format locus identifier.
    - Save the results to a TSV file based on input file name and user set output directory.


## ChatGPT Integration: Setup OpenAI API Endpoint

Set environment variables in your terminal, as follows:
```bash
export OPENAI_ORGANIZATION='<KEY>'        # https://platform.openai.com/settings/organization/general
export OPENAI_PROJECT='<KEY>'             # https://platform.openai.com/settings/ -> project
export OPENAI_API_KEY='<KEY>'             # https://platform.openai.com/settings/profile?tab=api-keys
```

Then you will be able to load the keys into your project, as follows:
```python
from openai import OpenAI
import os

organization = os.getenv('OPENAI_ORGANIZATION')
project = os.getenv('OPENAI_PROJECT')
api_key = os.getenv('OPENAI_API_KEY')

openai_client = OpenAI(
    organization=organization,
    project=project,
    api_key=api_key
)
```

## Load in sumstats

Use the path to folder containing your summary statistics files.

Also, set up your output directory for
1) Restructured sumstats files
2) Lead SNP files

```python
import glob

# set input directory
directory_of_sumstats_path = '/path/to/sumstats'
my_input_files = glob.glob(directory_of_sumstats_path + '/*')

for path in my_input_files:
    print(path)

# set output directory
output_directory = 'path/to/my/output/dir'
```

## Instantiate Preprocess class

```python
from preprocessing import Preprocess

ft = Preprocess(
    client=openai_client, 
    out_dir=output_directory
)

ft.significance_threshold = 5e-8 # default value is 5e-8
ft.ancestry = 'EUR' # default value is 'EUR' but can be statically
# changed here or dynamically changed in the loop below.
# Ancestry may play a role in GPT's selection of a column, e.g., say there
# are two pval columns (1) pval_afr (2) pval_eur, then GPT will select
# the column matching the ancestry variable set by user.
```

## Iterate over each file and run preprocessing

```python
for path in my_input_files:
    print(f'==> {path}')

    # 1) restructure sumstats file
    ft.loadmap_sumstats_table(
        path, 
        verbose=True
    )

    # 2) create lead SNP table
    ft.create_leadsnp_table(verbose=True)
```

## Verify files were outputted to desired location

```python
fps = glob.glob(f'{output_directory}/*')
print(fps)
```

## Run batch job in cluster using Slurm

Running `submit_preprocessing.sh` involves a few steps:

1. Get all filenames of interest into one folder and create a file that lists these filenames, as follows:

```{bash}
# first move all sumstats files of interest into one input folder
# cd into that folder
(base) [sfriedman@pe2-login01 ~]$ cd /gpfs/commons/groups/sanjana_lab/mdrabkin/gwas_data/raw/European

# use ls <pattern>* to get the filenames of interest into a file:
(base) [sfriedman@pe2-login01 European]$ ls GCST* > ~/filepaths.txt # all files with GCST pattern should be listed here

# if you need to manually edit this file afterwards, here's one method:
(base) [sfriedman@pe2-login01 European]$ vim ~/filepaths.txt # manual filtering - only include desired files

```

2. Edit `submit_preprocessing.sh`:

```{bash}
#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --output=./logs/preprocess_%A_%a.out    # set this to your desired location
#SBATCH --array=1-68                            # use wc -l ~/filepaths.txt to get the count
#SBATCH --cpus-per-task=2                       #   & use as many #'s as files you have
#SBATCH --mem=64G

# fill in the following fields:

export OPENAI_ORGANIZATION=""
export OPENAI_PROJECT=""
export OPENAI_API_KEY=""

export PREPROC_INPUT_DIR=""
export PREPROC_OUTPUT_DIR=""
export PREPROC_FILEPATHS="~/filepaths.txt" # however you named your filepaths file

...

```

3. Now run with `sbatch submit_preprocessing.sh`
- Use `squeue -u <your-username>` for updates.
- Use `scancel <job-id>` to cancel.