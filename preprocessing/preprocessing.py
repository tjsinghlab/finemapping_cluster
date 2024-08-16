import os, pandas as pd
import shutil
import numpy as np
import subprocess
from openai import OpenAI

class Cols:
    chromosome = 'chromosome'
    position = 'position'
    allele1 = 'allele1'
    allele2 = 'allele2'
    beta = 'beta'
    se = 'se'
    pval = 'pval'
    odds_ratio = 'or'

class Preprocess:

    DEFAULT_SUMSTATS_COLS = {
        Cols.chromosome : 'str', 
        Cols.position   : 'int32', 
        Cols.allele1    : 'str', 
        Cols.allele2    : 'str', 
        Cols.beta       : 'float64', 
        Cols.se         : 'float64',
        Cols.pval       : 'float64',
        Cols.odds_ratio : 'float64'
    }

    DELIM_MAPPINGS = {
        '.tsv': '\t',
        '.csv': ',',
        '.txt': ' ',
        '.tbl': '\t'
    }

    DEFAULT_ANCESTRY = 'EUR'
    DEFAULT_SIG_THRESHOLD = 5e-8
    DEFAULT_WINDOW = 5e5
    CHUNK_SIZE = 2e6

    def __init__(self, client, out_dir, **kwargs):
        """
        Initializes the Preprocess class.

        Parameters:
            client : object
                Client object used for communication.
            out_dir : str
                Output directory where results will be saved.
            kwargs : dict
                Additional arguments:
                    desired_columns_dict : dict, optional
                        These are the column names your downstream pipeline expects
                        of your sumstats file as keys and their corresponding dtypes
                        as values (defaults to DEFAULT_SUMSTATS_COLS).
                    significance_threshold : float, optional
                        Significance threshold for filtering SNPs.
                        - loadmap_sumstats_table defaults to None.
                        - create_leadsnp_table defaults to DEFAULT_SIG_THRESHOLD.
                    ancestry : str, optional
                        Set a default ancestry in case column names include this meta information
                        and more information is needed to discern column placement (default is 'EUR').
        
        """
        self.client = client
        self.out_dir = out_dir
        self.desired_columns_dict = kwargs['desired_columns_dict'] \
                                    if 'desired_columns_dict' in kwargs \
                                    else Preprocess.DEFAULT_SUMSTATS_COLS
        self.desired_columns = list(self.desired_columns_dict.keys())
        self.significance_threshold = kwargs['significance_threshold'] \
                                        if 'significance_threshold' in kwargs \
                                        else Preprocess.DEFAULT_SIG_THRESHOLD
        self.ancestry : str = kwargs['ancestry'] \
                                        if 'ancestry' in kwargs \
                                        else Preprocess.DEFAULT_ANCESTRY


    def ask_sumstats_col_order(self, actual_columns : list, **kwargs) -> list:
        """
        Compares the given file column names to the expected fields and order.
        Determines which columns from the sumstats file are needed for the finemapping pipeline
        and returns those column names in the order expected by the pipeline.

        Parameters:
            actual_columns : list
                These are the column names as found in your sumstats file.
                The goal is to preserve these names but pick out the ones of interest
                and order them as needed to match the finemapping pipeline.
            kwargs : dict
                Additional arguments:
                    print : bool, optional
                        If True, print the generated query (default is False).

        Returns:
            filtered_ordered_cols : list
                List of columns from sumstats file ordered and mapping directly to
                the finemapping pipeline.

        """
        actual_columns_str = ','.join(actual_columns)
        desired_columns_str = ','.join(self.desired_columns)

        query = f"""
            I have a GWAS summary statistics file with columns: {actual_columns_str}. 
            
            Map the file to my sequence of column names: {desired_columns_str}.

            Use NA to represent columns that have no mapping in the GWAS summary statistics file.

            If there are columns that are distinguishable by references to ancestry (e.g.,EUR,EAS,AFR,etc.),
            select {self.ancestry} and disregard other ancestries.
            
            How to discern the effect allele? 
            Typically when you see a "ref" (reference) / "alt" (alternate) allele combination, the "alt" allele = allele1 and the "ref" allele = allele2.
            Typically when you see an "effect allele" / "other allele" combination, the "effect" allele = allele1 and the "other" allele = allele2.
            The effect allele can also be described as the coding allele (and therefore the other allele could be called the non-coding allele).

            Example 1:
            My sequence:    "{desired_columns_str}"
            GWAS file:      "chr,RSID,pos,ref,alt,MAFreq,p,error,b,odds_ratio"

            Mapping:        "chr,pos,ref,alt,b,error,p,or"
            Notes:
                chromosome ->  chr: Chromosome number.
                position  ->  pos: BP stands for base pair,which corresponds to the position of the SNP.
                allele1  ->  alt: Alternate allele, a.k.a. minor allele (alt), often written as Allele 1 (A1).
                allele2  ->  ref: Major allele, a.k.a. reference allele (ref), often written as Allele 2 (A2).
                beta    ->  b: Beta value (B).
                se  ->  error: Standard error (SE).
                pval   ->  p: P-value (P).
                or -> odds_ratio: Odds Ratio (OR).

            Example 2:
            My sequence:    "{desired_columns_str}"
            GWAS file:      "snp,chr,bp,id,maf,noncoding_snp,effect_snp,pval_afr,pval_eur,pval_cas,beta,standard_error"

            Mapping:        "chr,bp,NA,NA,beta,standard_error,pval_{self.ancestry.lower()},NA"
            Notes:
                chromosome ->  chr: Chromosome number.
                position  ->  bp: BP stands for base pair,which corresponds to the position of the SNP.
                allele1  ->  effect_snp: Effect allele corresponds to allele1.
                allele2  ->  noncoding_snp: The non-coding SNP, i.e., NOT the coding SNP, is NOT the effect allele and therefore corresponds to allele2.
                beta    ->  beta: Beta value (B).
                se  ->  standard_error: Standard error (SE).
                pval   ->  pval_{self.ancestry.lower()}: P-value (P). pval_{self.ancestry.lower()} selected among pvalues
                            over other similar options because we are focusing on {self.ancestry}.
                or  -> NA: Missing value.

            Other examples of p-value columns include: neg_log_10_p_value.
            Other examples of beta columns include logOR since beta = log(OR).
            Other examples of position columns include: base_pair_location.

            REMEMBER: the effect (i.e., coding, alternate, alt, ea) allele ALWAYS preceeds the non-effect (i.e. non-coding, other, reference, ref, oa) allele.
            
            Do not respond with anything other than the column name mapping as a comma separated string. 
            Make sure this string output is separated by commas and nothing else.

        """

        if 'print' in kwargs and kwargs['print']:
            print(query)
        
        completion = self.client.chat.completions.create(
            model="gpt-4o",
            messages=[
                {"role": "system", "content": "You are a helpful assistant."},
                {"role": "user", "content": query}
            ]
        )

        filtered_ordered_cols = completion.choices[0].message.content.split(',')
        return filtered_ordered_cols

    def gpt_qc(self, gpt_output):
        ## QC
        supposed_ea = gpt_output[2]
        for m in ['ref', 'other', 'oa']:
            if m in supposed_ea.lower():
                gpt_output[2] = gpt_output[3]
                gpt_output[3] = supposed_ea

        if supposed_ea == 'allele1' or supposed_ea == 'allele2' or supposed_ea == 'allele0':
            gpt_output[2] = 'NA'

        return gpt_output
    
    def remap_dataframe(self, df : pd.DataFrame, name_mappings : dict, cols_in_order : list) -> pd.DataFrame:
        """
            Renames and rearranges the columns of a DataFrame according to specified mappings and desired order.

            Parameters:
            ----------
            df : pd.DataFrame
                The DataFrame to be modified.
            name_mappings : dict
                A dictionary where keys are the current column names and values are the new column names.
            cols_in_order : list
                A list specifying the desired order of the columns in the resulting DataFrame.
                
            Returns:
            -------
            pd.DataFrame
                A DataFrame with columns renamed and rearranged according to the specified mappings and desired order.

            Example:
            -------
            > df = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6], 'c': [7, 8, 9]})
            > name_mappings = {'a': 'x', 'b': 'y', 'c': 'z'}
            > cols_in_order = ['z', 'y', 'x']
            > remap_dataframe(df, name_mappings, cols_in_order)
               z  y  x    
            0  7  4  1
            1  8  5  2
            2  9  6  3
        """
        df.rename(columns=name_mappings, inplace=True) # rename
        return df[cols_in_order] # rearrange

    def sort_input_file(self, sumstats_mapped_columns):
        """
        Sorts the input file by chromosome and position columns.
        
        """
        new_old_column_mappings = {
            s: e 
            for s, e in zip(self.desired_columns, sumstats_mapped_columns) if e != 'NA'
        }

        # Sort input file by chromosome and position column names
        chrom_col = new_old_column_mappings[Cols.chromosome]
        chrom_col_idx = self.input_cols.index(chrom_col) + 1
        pos_col = new_old_column_mappings[Cols.position]
        position_col_idx = self.input_cols.index(pos_col) + 1

        print(f"Sorting table by chromosome and position columns:")
        print(f"{chrom_col}:{chrom_col_idx}, {pos_col}:{position_col_idx}")

        # Create a temporary copy of the compressed input file in the output directory
        temp_copy_path = os.path.join(self.out_dir, self.filename) + self.ext_full
        shutil.copy(self.dataset_path, temp_copy_path)

        # Decompress the temporary file
        if self.ext2 == '.gz':
            print(f"Decompressing temporary copy {temp_copy_path}")
            subprocess.run(f"gzip -d {temp_copy_path}", shell=True, check=True)
            temp_file = os.path.join(self.out_dir, self.filename) + self.ext1
        else:
            temp_file = temp_copy_path

        # Sort the decompressed file while keeping the header
        sorted_input_file = os.path.join(self.out_dir, self.filename) + '_sorted' + self.ext1

        # Extract the header
        subprocess.run(f"head -n 1 {temp_file} > {sorted_input_file}", shell=True, check=True)

        # Sort the rest of the file
        subprocess.run(f"tail -n +2 {temp_file} | sort -k{chrom_col_idx},{chrom_col_idx}n -k{position_col_idx},{position_col_idx}n >> {sorted_input_file}", shell=True, check=True)

        # Clean up the decompressed temporary file if it was created
        if self.ext2 == '.gz':
            os.remove(temp_file)

        print(f"Sorted file saved to {sorted_input_file}")
        return sorted_input_file

    def loadmap_sumstats_table(self, dataset_path, **kwargs) -> int:
        """
        Loads and maps the summary statistics file to the expected column order for the finemapping pipeline.

        Parameters:
            dataset_path : str
                Path to the summary statistics file.
            kwargs : dict
                Additional arguments:
                    verbose : bool, optional
                        If True, print verbose output (default is False).

        Returns:
            None, however sumstats_df is saved as self.sumstats_df

        """
        self.dataset_path = dataset_path
        v = True if ('verbose' in kwargs and kwargs['verbose']) else False

        ## Get the base filename
        filename_with_ext = os.path.basename(dataset_path)

        ## Split the filename and extension
        filename, ext1 = os.path.splitext(filename_with_ext)

        ## Read in column names
        ## Handle the case for double extension (e.g., .txt.gz)
        if ext1 == '.gz':
            filename, ext2 = os.path.splitext(filename)
            ext_full = ext2 + ext1
            tmp = ext2
            ext2 = ext1
            ext1 = tmp
        else:
            ext_full = ext1
            ext2 = None


        ## Handle .sumstats ext
        prev_f = filename
        filename, ext0 = os.path.splitext(filename)
        if ext0 == '.sumstats':
            ext_full = ext0 + ext_full
        else:
            filename = prev_f
            # self.sumstats_df = None
            # print('[loadmap_sumstats_table] Skipping: Not a sumstats file (.sumstats.<ext>[.gz]).')
            # return
        

        if v:
            print(f'found {filename_with_ext}')
            # print(ext0, ext1, ext2, ext_full)
            print()

        self.ext1 = ext1
        self.ext2 = ext2
        self.ext_full = ext_full
        self.filename = filename

        dm_keys = Preprocess.DELIM_MAPPINGS.keys()
        if ext1 not in dm_keys and ext1 != '.gz':
            print(f"Only supports {dm_keys} (+ .gz)")
            return 0

        # ChatGPT integration
        ## Determine column names in sumstats file
        cols_df = pd.read_table(dataset_path, nrows=0)
        input_cols = cols_df.columns.tolist()
        self.input_cols = input_cols

        ## Map column names using GPT
        gpt_output = self.ask_sumstats_col_order(input_cols, print=False)
        sumstats_mapped_columns = self.gpt_qc(gpt_output)[:-1] # TODO: check if last column needed to fill in info

        ## Count NA's
        ct = sum(1 for n in sumstats_mapped_columns if n == 'NA')

        if v:
            print(f"Sumstat cols:\t{input_cols} ->\nGPT mapping:\t{sumstats_mapped_columns}")
            
        os.makedirs(self.out_dir, exist_ok=True)
        with open(self.out_dir + '/column_mappings.log', 'a') as out:
            fnout = filename if ct != 0 else f'{filename} *'
            out.write(f"{fnout}\nUSR\t{input_cols}\nGPT\t{sumstats_mapped_columns}\n\n")
        
        ## Don't process file if there are NAs
        if ct != 0:
            print("GPT couldn't make out all the columns.")
            return 0
    
        if v:
            print()
        
        # Get mappings
        try:
            sorted_input_file = self.sort_input_file(sumstats_mapped_columns)
        except:
            print('Failed to sort input file.')
            return 0

        old_new_column_mappings = {
            s: e 
            for s, e in zip(sumstats_mapped_columns, self.desired_columns) if s != 'NA'
        }
        
        cols_in_order = [c for c in self.desired_columns if c in old_new_column_mappings.values()]

        ## Load in dataframe
        save_sumstats_as = os.path.join(self.out_dir, filename) + f'_preprocessed.tsv' # {ext_full}
        mode = 'w'
        header = True

        print(input_cols)
        print(sumstats_mapped_columns)

        for chunk in pd.read_table(
                sorted_input_file,
                sep=self.DELIM_MAPPINGS[ext1],
                usecols=sumstats_mapped_columns,
                dtype=self.desired_columns_dict,
                chunksize=Preprocess.CHUNK_SIZE
            ):
            # print(f'loaded in {df.shape[0]} SNPs')
            print('='*15)
            print(f'Processing chunk of {chunk.shape[0]} SNPs')
            print('before mapping:', list(chunk.columns))

            # rename columns to format desired by pipeline
            chunk = self.remap_dataframe(chunk, old_new_column_mappings, cols_in_order)
            print('after mapping:', list(chunk.columns))

            # filter
            chunk[Cols.allele1] = chunk[Cols.allele1].apply(lambda s : s.upper() if len(s) == 1 else np.nan)
            chunk[Cols.allele2] = chunk[Cols.allele2].apply(lambda s : s.upper() if len(s) == 1 else np.nan)
            chunk.dropna(subset=cols_in_order, inplace=True)
            print(f'Filtered chunk to {chunk.shape[0]} SNPs (drop SNPs with NA alleles -> alleles of len > 1 set to NA)')

            chunk.sort_values(by=Cols.pval, ascending=True, inplace=True)
            chunk.drop_duplicates(subset=[Cols.chromosome, Cols.position], keep='first', inplace=True)
            chunk.sort_values(by=[Cols.chromosome, Cols.position], inplace=True)
            print(f'Filtered chunk to {chunk.shape[0]} SNPs (drop duplicate SNPs -> keep SNP with lowest P-value)')

            chunk_filt = chunk[(abs(chunk[Cols.beta]) < np.inf) & (abs(chunk[Cols.beta]) > 0)]
            print(f'Filtered chunk to {chunk_filt.shape[0]} SNPs (drop irregular betas -> keep SNPs with abs(beta) > 0 and abs(beta)  inf)')
            # thrown_away = chunk[~ ((abs(chunk[Cols.beta]) < np.inf) & (abs(chunk[Cols.beta]) > 0))]
            # print(set(list(thrown_away.beta)))

            chunk_filt.to_csv(save_sumstats_as, mode=mode, header=header, index=False, sep='\t') # sep=',' if ext1 == '.csv' else '\t', compression='gzip' if ext2 else None
            
            mode = 'a'
            header = False

            if v:
                print(f'Saved chunk to \n{save_sumstats_as}')
            
        self.save_sumstats_as = save_sumstats_as
        
        if v:
            print(f'Saved reformatted sumstats file to: {filename}_preprocessed.tsv')
            print()

        return 1


    def create_leadsnp_table(self, **kwargs) -> pd.DataFrame:
        """
        Creates a table of lead SNPs from the summary statistics file.

        Parameters:
            kwargs : dict
                Additional arguments:
                    verbose : bool, optional
                        If True, print verbose output (default is False).

        Returns:
            leadsnp_df : DataFrame or None
                A DataFrame containing lead SNPs with columns 'CHR', 'BP', and 'locus'. 
                Returns None in case of an error.
                The DataFrame is also stored in self.leadsnp_df.
        """
        v = True if ('verbose' in kwargs and kwargs['verbose']) else False

        chunks = []
        for chunk in pd.read_table(
                self.save_sumstats_as,
                sep='\t',
                chunksize=Preprocess.CHUNK_SIZE
            ):

            chunk = chunk[chunk['pval'] <= self.significance_threshold].copy() # 5e-8
            chunk = chunk.dropna()
            chunks.append(chunk)

        sumstats_df = pd.concat(chunks)
        sumstats_df = sumstats_df[(abs(sumstats_df['beta']) < np.inf) & (abs(sumstats_df['beta']) > 0)]
        sumstats_df.reset_index(drop=True, inplace=True)

        if v:
            print(f'filtered down to {sumstats_df.shape[0]} SNPs with P < {self.significance_threshold}')

        ## Convert chromosome "X" to 23, "Y" to 24 and ensure numeric types
        ## Sort by chromosome and position
        sumstats_df[Cols.chromosome] = sumstats_df[Cols.chromosome].replace('X', 23)
        sumstats_df[Cols.chromosome] = sumstats_df[Cols.chromosome].infer_objects(copy=False)
        sumstats_df[Cols.chromosome] = sumstats_df[Cols.chromosome].replace('Y', 24).astype(int)
        sumstats_df[Cols.chromosome] = sumstats_df[Cols.chromosome].infer_objects(copy=False)
        sumstats_df.sort_values(by=[Cols.chromosome, Cols.position], inplace=True)
        sumstats_df.reset_index(drop=True, inplace=True)

        ## Assign previous chromosome and position for comparison
        sumstats_df['prev_Chr'] = sumstats_df[Cols.chromosome].shift(1, fill_value=sumstats_df.chromosome.iloc[0])
        sumstats_df['prev_Position'] = sumstats_df[Cols.position].shift(1, fill_value=sumstats_df.position.iloc[0])

        ## Assign locus numbers
        ### First assign truth values based on conditions below
        ### Then create running sum of truths -> +1 for every true indicates new locus
        sumstats_df['locus'] = (sumstats_df[Cols.chromosome] != sumstats_df['prev_Chr']) | \
                        (abs(sumstats_df[Cols.position] - sumstats_df['prev_Position']) > Preprocess.DEFAULT_WINDOW)
        sumstats_df['locus'] = sumstats_df['locus'].cumsum() + 1

        ## Remove temporary columns
        sumstats_df.drop(columns=['prev_Chr', 'prev_Position'], inplace=True)

        ## Group by locus and filter for minimum p-value
        leadsnp_df = sumstats_df.loc[sumstats_df.groupby('locus')['pval'].idxmin()]
        leadsnp_df.reset_index(drop=True, inplace=True)

        ## Remove duplicated loci
        leadsnp_df.drop_duplicates(subset='locus', inplace=True)

        ## Transform columns and format locus identifier
        leadsnp_df.rename(columns={Cols.chromosome:'CHR', Cols.position:'BP'}, inplace=True)
        leadsnp_df['CHR'] = leadsnp_df['CHR'].replace(23, 'X')
        leadsnp_df['CHR'] = leadsnp_df['CHR'].replace(24, 'Y').astype(str)
        leadsnp_df['CHR'] = leadsnp_df['CHR'].infer_objects(copy=False)
        leadsnp_df['locus'] = leadsnp_df['CHR'].astype(str) + '.' + leadsnp_df['BP'].astype(str)
        leadsnp_df = leadsnp_df[['CHR', 'BP', 'locus']]

        ## Save the results to a TSV file
        save_leadsnp_as = os.path.join(self.out_dir, self.filename) + f'_preprocessed_leadSNPs.tsv'
        leadsnp_df.to_csv(save_leadsnp_as, index=False, sep='\t')

        if v:
            print(leadsnp_df.head())
            print()

            print(f'Saved lead SNP file to \n{save_leadsnp_as}')
            print()

        return leadsnp_df
