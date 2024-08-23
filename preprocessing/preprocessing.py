import os, sys
import shutil
import subprocess
import random
import time
import gzip
import bz2
import csv
import chardet
import yaml
import pandas as pd
import numpy as np

import openai
import google.generativeai as genai
from google.api_core import exceptions as google_exceptions

from cols import Cols

class Prints:
    """
    Prints the arguments if verbose is True.
    Appends the arguments to a log file 'preprocess.log'.
    """
    def __init__(self, verbose, **kwargs):
        self.verbose = verbose
        self.out_dir = kwargs['out_dir'] if 'out_dir' in kwargs else '.'
        self.filename = kwargs['filename'] if 'filename' in kwargs else 'preprocess'

    def __call__(self, *args):

        if self.verbose:
            print(*args)
            
        with open(f'{self.out_dir}/{self.filename}.log', 'a') as f:
            print(*args, file=f)


class Preprocess:

    DEFAULT_SUMSTATS_COLS = {
        Cols.chromosome : 'str', 
        Cols.position   : 'int32', 
        Cols.allele1    : 'str', 
        Cols.allele2    : 'str', 
        Cols.beta       : 'float64', 
        Cols.se         : 'float64',
        Cols.pval       : 'float64'
    }
    VALID_DELIMITERS = [',', ' ', '\t', ';', '|', ':']
    DEFAULT_ANCESTRY = 'EUR'
    DEFAULT_SIG_THRESHOLD = 5e-8
    DEFAULT_WINDOW = 5e5
    CHUNK_SIZE = 2e6
    MAX_RETRIES = 10

    def __init__(self, client, out_dir, **kwargs):
        """
        Initializes the Preprocess class, which preprocesses summary statistics files for fine-mapping.

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
                    desired_columns : list, optional
                        These are the column names your downstream pipeline expects
                        of your sumstats file (defaults to list of keys in desired_columns_dict).
                    significance_threshold : float, optional
                        Significance threshold for filtering SNPs.
                        - loadmap_sumstats_table defaults to None.
                        - create_leadsnp_table defaults to DEFAULT_SIG_THRESHOLD.
                    ancestry : str, optional
                        Set a default ancestry in case column names include this meta information
                        and more information is needed to discern column placement (default is 'EUR').
        
        """
        self.client = client
        self.out_dir = os.path.expanduser(out_dir)
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

    @staticmethod
    def check_allele(s : str) -> str | None:
        """
        Checks if the input string is a valid allele.
        If the input is not a string or has length different from 1, returns np.nan.
        Otherwise, returns the input string converted to uppercase.

        Parameters:
            s : str
                The input string to be checked.

        Returns:
            str or np.nan
                The input string converted to uppercase if it is a valid allele, or np.nan otherwise.
        """
        if type(s) != str or len(s) != 1:
            return np.nan
        return s.upper()

    def __retry(self):
        if self.retries >= Preprocess.MAX_RETRIES:
            return
        # Exponential backoff with some jitter
        delay = (2 ** self.retries) + random.uniform(60, 70)
        print(f"Retrying in {delay:.2f} seconds...")
        time.sleep(delay)
        self.retries += 1
        
    def __ask_sumstats_col_order(self, actual_columns : list, **kwargs) -> list:
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
        actual_columns_str = ','.join(list(actual_columns))
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
            GWAS file:      "chr,RSID,pos,ref,alt,MAFreq,p,error,b"

            Mapping:        "chr,pos,ref,alt,b,error,p"
            Notes:
                chromosome ->  chr: Chromosome number.
                position  ->  pos: BP stands for base pair,which corresponds to the position of the SNP.
                allele1  ->  alt: Alternate allele, a.k.a. minor allele (alt), often written as Allele 1 (A1).
                allele2  ->  ref: Major allele, a.k.a. reference allele (ref), often written as Allele 2 (A2).
                beta    ->  b: Beta value (B).
                se  ->  error: Standard error (SE).
                pval   ->  p: P-value (P).

            Example 2:
            My sequence:    "{desired_columns_str}"
            GWAS file:      "snp,chr,bp,id,maf,noncoding_snp,effect_snp,pval_afr,pval_eur,pval_cas,beta,standard_error"

            Mapping:        "chr,bp,effect_snp,noncoding_snp,beta,standard_error,pval_{self.ancestry.lower()}"
            Notes:
                chromosome ->  chr: Chromosome number.
                position  ->  bp: BP stands for base pair,which corresponds to the position of the SNP.
                allele1  ->  effect_snp: Effect allele corresponds to allele1.
                allele2  ->  noncoding_snp: The non-coding SNP, i.e., NOT the coding SNP, is NOT the effect allele and therefore corresponds to allele2.
                beta    ->  beta: Beta value (B).
                se  ->  standard_error: Standard error (SE).
                pval   ->  pval_{self.ancestry.lower()}: P-value (P). pval_{self.ancestry.lower()} selected among pvalues
                            over other similar options because we are focusing on {self.ancestry}.

            Other examples of p-value columns include: neg_log_10_p_value, p_value.
            Other examples of beta columns include logOR since beta = log(OR).
            Other examples of chromosome columns include: chr, chromosome_number, chromosome. Note that variant_id is NOT a chromosome column.
            Other examples of position columns include: base_pair_location.

            REMEMBER: the effect (i.e., coding, alternate, alt, ea) allele ALWAYS preceeds the non-effect (i.e. non-coding, other, reference, ref, oa) allele.
            
            Do not respond with anything other than the column name mapping as a comma separated string. 
            Make sure this string output is separated by commas and nothing else.

        """

        if 'print' in kwargs and kwargs['print']:
            print(query)
        
        try:
            model_gemini = self.client.model_name == "models/gemini-pro"
        except:
            model_gemini = False

        self.retries = 0
        while self.retries < self.MAX_RETRIES:
            try:

                if model_gemini:
                    # GEMINI
                    completion = self.client.generate_content(query).candidates[0].content.parts[0].text.split(',')
                    # print(completion)
                    return completion

                else:
                    # GPT
                    completion = self.client.chat.completions.create(
                        model="gpt-4o",
                        messages=[
                            {"role": "system", "content": "You are a helpful assistant."},
                            {"role": "user", "content": query}
                        ]
                    )

                    # Process the response
                    filtered_ordered_cols = completion.choices[0].message.content.split(',')
                    return filtered_ordered_cols

            except (ValueError, google_exceptions.GoogleAPICallError) as e:
                if "No API_KEY or ADC found" in str(e):
                    print("API Key Error: No API_KEY or ADC found. Please set up your API key.")
                    print("You can set it by:")
                    print("1. Setting the GOOGLE_API_KEY environment variable, or")
                    print("2. Using genai.configure(api_key=your_api_key) in your code.")
                    self.retries = self.MAX_RETRIES + 1
                    break
                
            except openai.AuthenticationError as e:
                print("Authentication Error: Invalid or missing API key.")
                print("Please ensure you have set your OpenAI API key correctly.")
                print("You can set it by:")
                print("1. Setting the OPENAI_API_KEY environment variable, or")
                print("2. Passing it directly to the OpenAI client: OpenAI(api_key='your-api-key')")
                self.retries = self.MAX_RETRIES + 1
                break
                
            except google_exceptions.RetryError as e:
                # Handle rate limit error
                print(f"Rate limit exceeded: {e}")
                self.__retry()
            
            except openai.RateLimitError as e:
                print(f"Rate limit exceeded: {e}")
                self.__retry()
                
            except Exception as e:
                # Handle other unexpected errors
                print(f"An unexpected error occurred: {e}")
                self.__retry()
                
        if self.retries >= self.MAX_RETRIES:
            print("Max retries exceeded.")
            
        return []

    def __determine_file_delimiter(self, file_path : str):
        """
        Determine the delimiter of a given file (compressed or not).
        
        Args:
        file_path (str): Path to the input file.
        
        Returns:
        str: Detected delimiter or None if not determined.
        """
        def open_file(file_path):
            if file_path.endswith('.gz'):
                return gzip.open(file_path, 'rt')
            elif file_path.endswith('.bz2'):
                return bz2.open(file_path, 'rt')
            else:
                return open(file_path, 'r')
        
        try:
            with open_file(file_path) as file:
                # Read a sample of the file (first 1000 bytes)
                sample = file.read(1000)
                
                # Detect the file encoding
                encoding = chardet.detect(sample.encode())['encoding']
                
                # Reset file pointer
                file.seek(0)
                
                # Read the first few lines
                lines = [file.readline() for _ in range(5)]
                
                # Use csv.Sniffer to detect the dialect
                dialect = csv.Sniffer().sniff(''.join(lines), delimiters = self.VALID_DELIMITERS)
                
                return dialect.delimiter
        
        except Exception as e:
            print(f"Error determining delimiter: {str(e)}")
            sys.exit(1)

    def __copy_file(self, source_path : str, destination_path : str) -> bool:
        """
        Copy a file from source to destination with error checking.

        :param source_path: Path to the source file
        :param destination_path: Path where the file should be copied
        :return: True if successful, False otherwise
        """
        prints = self.prints
        try:
            # Check if source file exists
            if not os.path.exists(source_path):
                prints(f"Error: Source file does not exist: {source_path}")
                return False

            # Check if source is a file (not a directory)
            if not os.path.isfile(source_path):
                prints(f"Error: Source is not a file: {source_path}")
                return False

            # Check if we have read permissions on the source file
            if not os.access(source_path, os.R_OK):
                prints(f"Error: No read permission for source file: {source_path}")
                return False

            # Check if destination directory exists, if not, try to create it
            dest_dir = os.path.dirname(destination_path)
            if not os.path.exists(dest_dir):
                os.makedirs(dest_dir)

            # Check if destination file already exists, if so, try to delete it
            if os.path.exists(destination_path):
                prints(f"Destination file already exists: {destination_path}")
                prints(f"Warning: Make sure your input and output directories are different to avoid overwriting files.")
                prints(f"Delete the existing file to continue or quit.")
                cont = input("Continue? (Y/N): ")
                if cont.lower() != 'y':
                    sys.exit(1)
                os.remove(destination_path)
                prints(f"Deleted existing file.")


            # Perform the copy operation
            prints(f"Copying {source_path} to {destination_path}.")
            shutil.copy(source_path, destination_path)

            # Verify that the copy was successful
            if os.path.exists(destination_path) and os.path.getsize(source_path) == os.path.getsize(destination_path):
                prints(f"Successfully copied {source_path} to {destination_path}")
                return True
            else:
                prints(f"Error: Copy verification failed for {destination_path}")
                return False

        except PermissionError:
            prints(f"Error: Permission denied when copying {source_path} to {destination_path}")
        except shutil.SameFileError:
            prints(f"Error: Source and destination are the same file: {source_path}")
        except IOError as e:
            prints(f"I/O error occurred: {e}")
        except Exception as e:
            prints(f"Unexpected error occurred: {e}")

        return False
    
    def __run_command(self, command) -> bool:
        """
        Run a shell command and handle potential errors.
        
        :param command: The shell command to run
        :return: True if successful, False otherwise
        """
        prints = self.prints
        try:
            subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
            return True
        except subprocess.CalledProcessError as e:
            prints(f"Command failed with error: {e}")
            prints(f"Command output: {e.output}")
            return False
        except Exception as e:
            prints(f"An unexpected error occurred: {e}")
            return False

    def __is_gzip(self, filename):
        """
        Check if a file is a valid gzip file.

        :param filename: Path to the file
        :return: True if the file is a valid gzip file, False otherwise
        """
        # Check extension
        if not filename.lower().endswith('.gz'):
            return False
        
        # Check file size (gzip files are at least 18 bytes)
        if os.path.getsize(filename) < 18:
            return False
        
        # Check magic numbers
        with open(filename, 'rb') as f:
            if f.read(2) != b'\x1f\x8b':
                return False
        
        # Try to open with gzip
        try:
            with gzip.open(filename, 'rb') as f:
                f.read(1)
            return True
        except gzip.BadGzipFile:
            return False

    def __get_file_path_elements(self, dataset_path : str) -> bool:

        # Get absolute dataset path
        self.dataset_path = os.path.expanduser(dataset_path)

        # Get the base filename
        filename_with_ext = os.path.basename(dataset_path)

        # Split the filename and extension
        filename, ext1 = os.path.splitext(filename_with_ext)

        # Read in column names
        # Handle the case for double extension (e.g., .txt.gz)
        if self.__is_gzip(dataset_path):
            filename, ext2 = os.path.splitext(filename)
            ext_full = ext2 + ext1
            tmp = ext2
            ext2 = ext1
            ext1 = tmp
        else:
            ext_full = ext1
            ext2 = None

        # Handle .sumstats ext
        prev_f = filename
        filename, ext0 = os.path.splitext(filename)
        if ext0 == '.sumstats':
            ext_full = ext0 + ext_full
        else:
            filename = prev_f

        self.ext1 = ext1
        self.ext2 = ext2
        self.ext_full = ext_full
        self.filename = filename
        self.filename_with_ext = filename_with_ext

        return True


    def sort_input_file(self, sumstats_mapped_columns : list) -> bool:
        """""
        Sorts the input file by chromosome and position columns.
        
        :param sumstats_mapped_columns : list, List of column names in the input file after mapping.
        :returns: str, Path to the sorted input file.
        """""
        prints = self.prints
        new_old_column_mappings = {
            s: e 
            for s, e in zip(self.desired_columns, sumstats_mapped_columns) if e != 'NA'
        }

        # Sort input file by chromosome and position column names
        chrom_col = new_old_column_mappings[Cols.chromosome]
        chrom_col_idx = self.input_cols.index(chrom_col) + 1
        pos_col = new_old_column_mappings[Cols.position]
        position_col_idx = self.input_cols.index(pos_col) + 1

        prints(f"Sorting table by chromosome and position columns:")
        prints(f"{chrom_col}:{chrom_col_idx}, {pos_col}:{position_col_idx}")

        # Create a temporary copy of the compressed input file in the output directory
        temp_copy_path = os.path.join(self.out_dir, self.filename) + self.ext_full

        if not self.__copy_file(self.dataset_path, temp_copy_path):
            prints("Failed to copy file")
            return False

        # Decompress the temporary file
        if self.__is_gzip(self.dataset_path):
            temp_file = os.path.join(self.out_dir, self.filename) + self.ext1
            # Check if decompressed file already exists, if so, try to delete it
            if os.path.exists(temp_file):
                try:
                    print("Removing existing decompressed file.")
                    print(f"Warning: Make sure your input and output directories are different to avoid overwriting files.")
                    cont = input("Continue? (Y/N): ")
                    if cont.lower() != 'y':
                        sys.exit(1)
                    os.remove(temp_file)
                except PermissionError:
                    prints(f"Permission denied: Unable to delete {temp_file}")
                    return False
                except Exception as e:
                    prints(f"An error occurred: {e}")
                    return False

            prints(f"Decompressing temporary copy {temp_copy_path}")
            if not self.__run_command(f"gzip -d {temp_copy_path}"):
                prints("Failed to decompress file")
                # Handle the error (e.g., clean up, exit)
                return False
            
        else:
            temp_file = temp_copy_path

        # Sort the decompressed file while keeping the header
        sorted_input_file = os.path.join(self.out_dir, self.filename) + '_sorted' + self.ext1
    
        # Extract the header
        if not self.__run_command(f"head -n 1 {temp_file} > {sorted_input_file}"):
            prints("Failed to extract header")
            # Handle the error (e.g., clean up, exit)
            return False

        unsorted_head = subprocess.run(f"head -n 1 {temp_file}", shell=True, check=True, capture_output=True, text=True).stdout
        unsorted_head_list = unsorted_head.strip().split(self.delimiter)
        prints("File head:", unsorted_head_list)
        if len(unsorted_head_list) == 1:
            prints(f"Warning: Your file ends in {self.ext1} but are you sure your input file is separated correctly?")

        # Sort the rest of the file
        sort_command = f"tail -n +2 {temp_file} | sort -k{chrom_col_idx},{chrom_col_idx}n -k{position_col_idx},{position_col_idx}n >> {sorted_input_file}"
        if not self.__run_command(sort_command):
            prints("Failed to sort the file")
            # Handle the error
            return False

        prints(f"Sorted file saved to {sorted_input_file}")
        self.sorted_input_file = sorted_input_file
        
        # Clean up temporary file
        try:
            os.remove(temp_file)
            prints(f"Successfully removed temporary file: {temp_file}")
        except FileNotFoundError:
            prints(f"Warning: Temporary file not found: {temp_file}")
        except PermissionError:
            prints(f"Error: Permission denied when trying to remove {temp_file}")
        except Exception as e:
            prints(f"An unexpected error occurred while removing {temp_file}: {e}")

        return True

    def gpt_qc(self, gpt_output : list) -> list:
        """
        Checks the GPT output for common errors and corrects them.
        If the effect allele is not mentioned in the output, the function assumes that the effect allele is allele1.
        If the effect allele is mentioned as 'allele1', 'allele2', or 'allele0', the function assumes that the effect allele is 'NA'.

        Parameters:
        ----------
        gpt_output : list
            The output from the GPT model.

        Returns:
        -------
        list
            The corrected output from the GPT model.
        """
        supposed_ea = gpt_output[2]
        for m in ['ref', 'other', 'oa']:
            if m in supposed_ea.lower():
                gpt_output[2] = gpt_output[3]
                gpt_output[3] = supposed_ea

        # if supposed_ea.lower() in [
        #         'allele1',
        #         'allele2',
        #         'allele0',
        #         'a1',
        #         'a2',
        #         'a0'
        #     ]:
        #     gpt_output[2] = 'NA'

        return [x.strip() for x in gpt_output]
    
    def remap_dataframe(self, df : pd.DataFrame, name_mappings : dict, cols_in_order : list) -> pd.DataFrame:
        """
        Renames and rearranges the columns of a DataFrame according to specified mappings and desired order.
        Drops rows with NA in any column.

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
        > df = pd.DataFrame({'REF': [A, T, C], 'ALT': [G, G, A], 'PVAL': [0.7, 0.08, 0.009]})

           REF ALT PVAL    
        0  A  G  0.7
        1  T  G  0.08
        2  C  A  0.009
        3  NA A  NA

        > name_mappings = {'REF': 'allele2', 'ALT': 'allele1', 'PVAL': 'pval'}
        > df.rename(columns=name_mappings, inplace=True)              # rename columns

           allele2 allele1 pval    
        0  A  G  0.7
        1  T  G  0.08
        2  C  A  0.009
        3  NA A  NA

        > cols_in_order = ['allele1', 'allele2', 'pval']
        > df_rearranged = df[cols_in_order].copy()                    # rearrange columns

           allele1 allele2 pval    
        0  G  A  0.7
        1  G  T  0.08
        2  A  C  0.009
        3  A  NA NA
        
        > df_rearranged.dropna(subset=cols_in_order, inplace=True)
        
           allele1 allele2 pval    
        0  G  A  0.7
        1  G  T  0.08
        2  A  C  0.009
        """
        df.rename(columns=name_mappings, inplace=True)              # rename columns
        df_rearranged = df[cols_in_order].copy()                    # rearrange columns
        df_rearranged.dropna(subset=cols_in_order, inplace=True)    # drop rows with NA in any column
        return df_rearranged
    
    def beta_filter(self, df) -> pd.DataFrame:
        """
        Filters the DataFrame by removing rows with irregular beta values.

        Parameters:
            df : pd.DataFrame
                The DataFrame to be filtered.

        Returns:
            pd.DataFrame
                The filtered DataFrame.
        """
        sumstats_df = df[(abs(df['beta']) < np.inf) & (abs(df['beta']) > 0)]
        sumstats_df.reset_index(drop=True, inplace=True)
        return sumstats_df
    
    def get_columns(self, dataset_path : str) -> list:
        """
        Return column names in sumstats file.
        """
        if not self.__get_file_path_elements(dataset_path):
            return []
        
        self.delimiter = self.__determine_file_delimiter(dataset_path)
        cols_df = pd.read_table(dataset_path, sep=self.delimiter, nrows=0)
        input_cols = cols_df.columns.tolist()
        return input_cols

    def suggest_columns(self, dataset_path : str) -> list | None:
        """ 
        Suggest column names in sumstats file using GPT.
        """
        if self.client is None:
            print('Client must be set for GPT integration.')
            return None
            
        ## Map column names using GPT
        gpt_output = self.__ask_sumstats_col_order(self.get_columns(dataset_path), print=False)
        if gpt_output:
            sumstats_mapped_columns = self.gpt_qc(gpt_output)
            return sumstats_mapped_columns
        return None
    
    def loadmap_sumstats_table(self, dataset_path, **kwargs) -> None:
        """
        Loads and maps the summary statistics file to the expected column order for the finemapping pipeline.

        Parameters:
            dataset_path : str
                Path to the summary statistics file.
            kwargs : dict
                Additional arguments:
                    manual_columns : list, optional
                        If list provided, GPT integration will be avoided and the columns list will
                        follow the selection of the list provided.
                    verbose : bool, optional
                        If True, print verbose output (default is False).

        Returns:
            int -> zero on success, != 0 on failure

        """
        if not self.__get_file_path_elements(dataset_path):
            return
        
        self.delimiter = self.__determine_file_delimiter(dataset_path)

        v = True if ('verbose' in kwargs and kwargs['verbose']) else False
        prints = Prints(v, out_dir = self.out_dir, filename = self.filename)
        self.prints = prints

        if self.out_dir == os.path.dirname(dataset_path):
            prints('Output directory is the same as input directory. Please specify a different output directory.')
            return

        prints(f'found {self.filename_with_ext}')
        prints()

        # ChatGPT integration
        # Determine column names in sumstats file
        input_cols = self.get_columns(dataset_path)
        self.input_cols = input_cols
        
        mc = True if ('manual_columns' in kwargs and kwargs['manual_columns']) else False
        if mc:
            prints('Using manual columns:')
            sumstats_mapped_columns = kwargs['manual_columns']
            prints(sumstats_mapped_columns)
        else:
            if self.client is None:
                prints('Client not set. Client needed to access GPT. Input manual_columns to proceed manually.')
                return
                
            # Map column names using GPT
            gpt_output = self.__ask_sumstats_col_order(input_cols, print=False)
            sumstats_mapped_columns = self.gpt_qc(gpt_output)

            # Count NA's
            ct = sum(1 for n in sumstats_mapped_columns if n == 'NA')

            prints(f"Sumstat cols:\t{input_cols} ->\nGPT mapping:\t{sumstats_mapped_columns}")
                
            os.makedirs(self.out_dir, exist_ok=True)
            with open(self.out_dir + '/column_mappings.log', 'a') as out:
                fnout = self.filename if ct != 0 else f'{self.filename} *'
                out.write(f"{fnout}\nUSR\t{input_cols}\nGPT\t{sumstats_mapped_columns}\n\n")
            
            # Don't process file if there are NAs
            if ct != 0:
                prints("GPT couldn't make out all the columns.")
                return
    
        prints()

        # Get mappings
        rv = self.sort_input_file(sumstats_mapped_columns)
        if not rv:
            return

        old_new_column_mappings = {
            s: e 
            for s, e in zip(sumstats_mapped_columns, self.desired_columns) if s != 'NA'
        }
        
        cols_in_order = [c for c in self.desired_columns if c in old_new_column_mappings.values()]

        # Load in dataframe
        save_sumstats_as = os.path.join(self.out_dir, self.filename) + f'_preprocessed.tsv' # {ext_full}
        mode = 'w'
        header = True

        prints(self.input_cols)
        prints(sumstats_mapped_columns)

        for chunk in pd.read_table(
                self.sorted_input_file,
                sep=self.delimiter,
                usecols=sumstats_mapped_columns,
                dtype=self.desired_columns_dict,
                chunksize=self.CHUNK_SIZE
            ):
            # prints(f'loaded in {df.shape[0]} SNPs')
            prints('='*15)
            prints(f'Processing chunk of {chunk.shape[0]} SNPs')
            prints('before mapping:', list(chunk.columns))
            chunk_size = chunk.shape[0]

            # rename columns to format desired by pipeline
            chunk = self.remap_dataframe(chunk, old_new_column_mappings, cols_in_order)
            prints('after mapping:', list(chunk.columns))

            # filter
            chunk[Cols.allele1] = chunk[Cols.allele1].astype(str)
            chunk[Cols.allele1] = chunk[Cols.allele1].apply(self.check_allele)

            chunk[Cols.allele2] = chunk[Cols.allele2].astype(str)
            chunk[Cols.allele2] = chunk[Cols.allele2].apply(self.check_allele)

            chunk.dropna(subset=cols_in_order, inplace=True)
            if chunk.shape[0] < chunk_size / 2:
                prints('[Quitting] Dropped more than half of the SNPs due to NA alleles and indels removed.')
                return
            chunk_size = chunk.shape[0]

            prints(f'Filtered chunk to {chunk.shape[0]} SNPs (drop SNPs with NA alleles -> alleles of len > 1 set to NA)')

            chunk.sort_values(by=Cols.pval, ascending=True, inplace=True)
            chunk.drop_duplicates(subset=[Cols.chromosome, Cols.position], keep='first', inplace=True)
            chunk.sort_values(by=[Cols.chromosome, Cols.position], inplace=True)
            
            if chunk.shape[0] < chunk_size * 0.75:
                prints('[Quitting] Dropped more than a quarter of the SNPs due to duplicates.')
                return
            chunk_size = chunk.shape[0]
            prints(f'Filtered chunk to {chunk.shape[0]} SNPs (drop duplicate SNPs -> keep SNP with lowest P-value)')

            chunk_filt = self.beta_filter(chunk) # chunk[(abs(chunk[Cols.beta]) < np.inf) & (abs(chunk[Cols.beta]) > 0)]
            if chunk_filt.shape[0] < chunk_size * 0.75:
                prints('[Quitting] Dropped more than a quarter of the SNPs due to irregular betas.')
                return
            chunk_size = chunk_filt.shape[0]
            prints(f'Filtered chunk to {chunk_filt.shape[0]} SNPs (drop irregular betas -> keep SNPs with abs(beta) > 0 and abs(beta) < inf)')
            # thrown_away = chunk[~ ((abs(chunk[Cols.beta]) < np.inf) & (abs(chunk[Cols.beta]) > 0))]
            # prints(set(list(thrown_away.beta)))
            
            # Ensure chromosome is string and replace 'X' with 23, 'Y' with 24
            chunk_filt.loc[:, Cols.chromosome] = chunk_filt.loc[:, Cols.chromosome].astype(str) \
                    .replace('23', 'X') \
                    .replace('24', 'Y') \
                    .infer_objects(copy=False)
            # A value is trying to be set on a copy of a slice from a DataFrame.  Try using .loc[row_indexer,col_indexer] = value instead

            chunk_filt.to_csv(save_sumstats_as, mode = mode, header = header, index = False, sep = '\t') 
            # sep=',' if ext1 == '.csv' else '\t', compression='gzip' if ext2 else None
            
            mode = 'a'
            header = False

            prints(f'Saved chunk to \n{save_sumstats_as}')
            
        prints(f'Saved reformatted sumstats file to: {save_sumstats_as}')
        self.save_sumstats_as = save_sumstats_as


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
        prints = Prints(v, out_dir = self.out_dir, filename = self.filename)
        self.prints = prints

        chunks = []
        for chunk in pd.read_table(
                self.save_sumstats_as,
                sep='\t',
                chunksize=Preprocess.CHUNK_SIZE
            ):

            chunk = chunk[chunk[Cols.pval] <= self.significance_threshold].copy() # 5e-8
            chunk = chunk.dropna()
            chunks.append(chunk)

        sumstats_df = pd.concat(chunks)
        sumstats_df = self.beta_filter(sumstats_df)

        if sumstats_df.shape[0] == 0:
            prints('No SNPs found.')
            return None

        prints(f'filtered down to {sumstats_df.shape[0]} SNPs with P < {self.significance_threshold}')

        # Convert chromosome "X" to 23, "Y" to 24 and ensure numeric types
        # Sort by chromosome and position
        sumstats_df[Cols.chromosome] = sumstats_df[Cols.chromosome].replace('X', 23)
        sumstats_df[Cols.chromosome] = sumstats_df[Cols.chromosome].infer_objects(copy=False)
        sumstats_df[Cols.chromosome] = sumstats_df[Cols.chromosome].replace('Y', 24).astype(int)
        # sumstats_df[Cols.chromosome] = sumstats_df[Cols.chromosome].infer_objects(copy=False)
        sumstats_df.sort_values(by=[Cols.chromosome, Cols.position], inplace=True)
        sumstats_df.reset_index(drop=True, inplace=True)

        # Assign previous chromosome and position for comparison
        sumstats_df['prev_Chr'] = sumstats_df[Cols.chromosome].shift(1, fill_value=sumstats_df.chromosome.iloc[0])
        sumstats_df['prev_Position'] = sumstats_df[Cols.position].shift(1, fill_value=sumstats_df.position.iloc[0])

        # Assign locus numbers
        # First assign truth values based on conditions below
        # Then create running sum of truths -> +1 for every true indicates new locus
        sumstats_df['locus'] = (sumstats_df[Cols.chromosome] != sumstats_df['prev_Chr']) | \
                        (abs(sumstats_df[Cols.position] - sumstats_df['prev_Position']) > Preprocess.DEFAULT_WINDOW)
        sumstats_df['locus'] = sumstats_df['locus'].cumsum() + 1

        # Remove temporary columns
        sumstats_df.drop(columns=['prev_Chr', 'prev_Position'], inplace=True)

        # Group by locus and filter for minimum p-value
        leadsnp_df = sumstats_df.loc[sumstats_df.groupby('locus')['pval'].idxmin()]
        leadsnp_df.reset_index(drop=True, inplace=True)

        # Remove duplicated loci
        leadsnp_df.drop_duplicates(subset='locus', inplace=True)

        # Transform columns and format locus identifier
        leadsnp_df.rename(columns={Cols.chromosome:'CHR', Cols.position:'BP'}, inplace=True)
        leadsnp_df['CHR'] = leadsnp_df['CHR'].replace(23, 'X')
        leadsnp_df['CHR'] = leadsnp_df['CHR'].replace(24, 'Y').astype(str)
        leadsnp_df['CHR'] = leadsnp_df['CHR'].infer_objects(copy=False)
        leadsnp_df['locus'] = leadsnp_df['CHR'].astype(str) + '.' + leadsnp_df['BP'].astype(str)
        leadsnp_df = leadsnp_df[['CHR', 'BP', 'locus']]

        # Save the results to a TSV file
        save_leadsnp_as = os.path.join(self.out_dir, self.filename) + f'_preprocessed_leadSNPs.tsv'
        leadsnp_df.to_csv(save_leadsnp_as, index=False, sep='\t')

        prints(leadsnp_df.head())
        prints()
        prints(f'Saved lead SNP file to \n{save_leadsnp_as}')
        prints()

        return leadsnp_df


if __name__ == "__main__":
    cols = Cols('header.yaml')

    for column, dtype in cols.columns.items():
        print(f"{column}: {dtype}")
