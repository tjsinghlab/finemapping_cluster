"""
PowerLift: A tool for lifting over RSIDs to genomic coordinates in input files.

Samuel M. Friedman
sfriedman@nygenome.org
"""

# outside libraries
import functools
import logging
from multiprocessing import Pool, cpu_count, Manager
import numpy as np
import os
import pandas as pd
from queue import Empty
import tempfile
from tqdm import tqdm
import unittest
from unittest.mock import patch, MagicMock
import yaml

# project libraries
import pathutils
from rsidb import RSIdb

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def log(func):
    """Decorator to log the execution of a function."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        logging.debug(f"Starting execution of {func.__name__}...")
        result = func(*args, **kwargs)
        logging.debug(f"Finished execution of {func.__name__}.")
        return result
    return wrapper

class PowerLift:
    
    DEFAULT_KEYS = ['pl_rsid', 'pl_chrom', 'pl_position', 'pl_ref', 'pl_alt']

    def __init__(self, config_path, output_dir = None, **kwargs):
        """Initialize the PowerLift class with a configuration file."""
        self.config = self.load_config(config_path)['powerlift']
        self.output_dir = os.path.expanduser(output_dir) if output_dir else os.path.dirname('~')
        self.db_path = kwargs['db_path'] if 'db_path' in kwargs else self.config.get('db_path')
        self.keys = list(self.config.get('columns').values()) if self.config.get('columns') else self.DEFAULT_KEYS
        self.chunk_size = kwargs['chunk_size'] if 'chunk_size' in kwargs else self.config.get('chunk_size', 25000)
        self.build = kwargs['build'] if 'build' in kwargs else self.config.get('genome_build', 'hg19')
        self.null_count = 0

    def load_config(self, config_path):
        """Load the configuration from a YAML file."""
        with open(config_path, "r") as file:
            return yaml.safe_load(file)

    @log
    def fetch_rsids(self, input_file, rsid_column):
        self.input_file = input_file
        self.rsid_column = rsid_column
        self.cmprsn = 'gzip' if pathutils.is_gzip(input_file) else None
        self.delim = pathutils.get_file_delimiter(input_file)
        output_filepath, ext = self.__get_output_file_name(input_file)
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)

        # Create a temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            self.temp_dir = temp_dir
            logging.info(f"Utilizing temporary directory {temp_dir}")

            with Manager() as manager:
                self.file_queue = manager.Queue()
                self.progress_queue = manager.Queue()

                with Pool(processes=cpu_count()) as pool:
                    futures = []
                    chunk_iterator = pd.read_csv(input_file, sep=self.delim, compression=self.cmprsn, chunksize=self.chunk_size)
                    self.total_chunks = sum(1 for _ in chunk_iterator)
                    chunk_iterator = pd.read_csv(input_file, sep=self.delim, compression=self.cmprsn, chunksize=self.chunk_size)
                    
                    for i, chunk in enumerate(chunk_iterator):
                        logging.debug(f"Processing chunk {i + 1}/{self.total_chunks}")
                        output_chunk_file = os.path.join(self.temp_dir, f"chunk_{i + 1:04d}{ext}")
                        futures.append(pool.apply_async(self._process_chunk, (chunk, output_chunk_file, i + 1)))

                    # Monitor progress
                    self.__monitor_progress()

                    # Close the pool and wait for the work to finish
                    pool.close()
                    pool.join()

                    # Check for any errors in the worker processes
                    for future in futures:
                        try:
                            future.get()
                        except Exception as e:
                            logging.error(f"Error in worker process: {e}")

                # Merge output files after all chunks have been processed
                self.__merge_output_files(output_filepath)
                logging.info(f"Processed data has been saved to {output_filepath}")

    def apply_powerlift(self, db, df_value):
        """Apply powerlifting to a single RSID value."""
        try:
            values = db.lift_rsid(rsid=df_value, build=self.build, verbose=False)
            if not values:
                self.null_count += 1
                return pd.Series([np.nan] * len(self.keys), index=self.keys)
            return pd.Series(values[0], index=self.keys)
        except Exception as e:
            self.null_count += 1
            logging.error(f"Error processing RSID {df_value}: {e}")
            return pd.Series([np.nan] * len(self.keys), index=self.keys)

    @log
    def _process_chunk(self, chunk, output_chunk_file, chunk_num):
        """Process a single chunk of data and write to the output file."""
        try:
            db = RSIdb(self.db_path)
            new_columns = chunk[self.rsid_column].apply(lambda x: self.apply_powerlift(db, x))
            chunk = pd.concat([chunk, new_columns], axis=1)
            chunk.to_csv(output_chunk_file, sep=self.delim, header=False, index=False)
            db.close()
            self.file_queue.put(output_chunk_file)
            self.progress_queue.put(1)  # Indicate progress
        except Exception as e:
            logging.error(f"Error processing chunk {chunk_num}: {e}")
            self.progress_queue.put(0)  # Indicate failure

    def __monitor_progress(self):
        """Monitor and display progress of chunk processing."""
        with tqdm(total=self.total_chunks, desc="Processing Chunks") as pbar:
            processed_chunks = 0
            while processed_chunks < self.total_chunks:
                try:
                    progress = self.progress_queue.get(timeout=1)
                    if progress:
                        pbar.update(1)
                        processed_chunks += 1
                except Empty:
                    continue

    @log
    def __merge_output_files(self, output_filepath):
        """Merge all processed chunk files into a single output file."""
        with open(output_filepath, 'w') as outfile:
            
            header_written = False
            processed_files = 0

            with tqdm(total=self.total_chunks, desc="Merging Files") as pbar:
                while processed_files < self.total_chunks:
                    try:
                        chunk_file = self.file_queue.get(timeout=5)
                        with open(chunk_file, 'r') as infile:
                            if not header_written:
                                self.__write_header(outfile)
                                header_written = True
                            outfile.write(infile.read())
                        processed_files += 1
                        pbar.update(1)
                        # No need to remove files here, they'll be deleted with the temp directory
                    except Empty:
                        continue
                    except Exception as e:
                        logging.error(f"Error merging files: {e}")
        
        logging.info(f"Number of NULL rsid mappings: {self.null_count}")

    def __get_output_file_name(self, input_file):
        """Generate the output file name based on the input file.
        :param input_file: Path to the input file.
        :return: Tuple containing the extension, and output file name.
        """
        base, ext = os.path.splitext(os.path.basename(input_file))
        self.cmprsn = 'gzip' if pathutils.is_gzip(input_file) else None
        if self.cmprsn == 'gzip':
            base, ext = os.path.splitext(base)
        out = f"{base}_powerlift{ext}"
        out_path = os.path.join(self.output_dir, out)
        return os.path.expanduser(out_path), ext
    
    def __write_header(self, file_handle):
        """Write the header to the output file."""
        columns = list(pd.read_csv(self.input_file, delimiter=self.delim, nrows=0).columns) + self.keys
        file_handle.write(f"{self.delim.join(columns)}\n")



class TestPowerLift(unittest.TestCase):
    def setUp(self):
        self.pl = PowerLift('config.yaml')

    def test_initialization(self):
        self.assertIsNotNone(self.pl.db_path)
        self.assertIsInstance(self.pl.keys, list)
        self.assertIsInstance(self.pl.chunk_size, int)
        self.assertIsInstance(self.pl.chunk_size, range(10, 1000000))
        self.assertIn(self.pl.build, ['hg19', 'hg38'])

    def load_db_path(self, config_path):
        """Load the configuration from a YAML file."""
        with open(config_path, "r") as file:
            config = yaml.safe_load(file)
        return config['database']['path']

    def test_initialization(self):
        self.assertIsNotNone(self.pl.db_path)
        self.assertIsInstance(self.pl.keys, list)

    def test_apply_powerlift(self):
        db = RSIdb(self.pl.db_path)
        result = self.pl.apply_powerlift(db, 'rs12345')
        self.assertIsInstance(result, pd.Series)
        self.assertEqual(len(result), len(self.pl.keys))
        db.close()

    def test_get_output_file_name(self):
        self.pl.output_dir = '/gpfs/commons/home/sfriedman/outest'
        input_file = '/gpfs/commons/home/sfriedman/GCST90093040_buildGRCh37.tsv'
        out, ext = self.pl._PowerLift__get_output_file_name(input_file)
        self.assertEqual(ext, '.tsv')
        expected_output = os.path.join(self.pl.output_dir, 'GCST90093040_buildGRCh37_powerlift.tsv')
        self.assertEqual(out, expected_output)

    @patch('pandas.read_csv')
    def test_write_header(self, mock_read_csv):
        # Mock the pandas read_csv function
        mock_df = pd.DataFrame(columns=['header1', 'header2'])
        mock_read_csv.return_value = mock_df

        self.pl.keys = ['pl_key1', 'pl_key2']
        self.pl.delim = ','
        self.pl.input_file = 'dummy_input_file.csv'

        mock_file_handle = MagicMock()

        self.pl._PowerLift__write_header(mock_file_handle)
        mock_read_csv.assert_called_once_with('dummy_input_file.csv', delimiter=',', nrows=0)

        expected_header = 'header1,header2,pl_key1,pl_key2\n'
        mock_file_handle.write.assert_called_once_with(expected_header)



if __name__ == "__main__":
    # Run tests
    unittest.main(exit=False)

    # Main execution
    pl = PowerLift(config_path = 'config.yaml', output_dir = '/gpfs/commons/home/sfriedman/outest')
    pl.fetch_rsids(input_file = '/gpfs/commons/home/sfriedman/GCST90092944_buildGRCh37.tsv.gz', rsid_column = 'variant_id')