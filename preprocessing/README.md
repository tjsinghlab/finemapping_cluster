# Fine-Mapping Preprocessing

**major caveat**: We cannot be certain from the columns alone which allele is the effect allele if it is not explicitly stated. For these cases, consult the README from the study or the paper itself. Additionally, this fails when chr & position are omitted. In such cases, rsIDs need to be mapped manually.

This package has three main methods:
1. **ask_sumstats_col_order**: Compares the given file column names to the expected fields. This method is used for standardizing summary stat files by mapping column name and ordering to an expected format.
    - Determines which columns from sumstats file are needed for the finemapping pipeline using a 2-shot prompt.
    - Returns those column names in the order expected by the pipeline.

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

## Instantiate FiMTools class

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

## Iterate over each file and run FiMTools

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
