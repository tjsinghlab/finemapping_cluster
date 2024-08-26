import os, sys
import gzip
import bz2
import csv
import chardet

VALID_DELIMITERS = [',', ' ', '\t', ';', '|', ':']

def is_gzip(file_name : str) -> bool:
    """
    Check if a file is a valid gzip file.

    :param file_name (str): Input file name.
    :return: True if the file is a valid gzip file, False otherwise
    """

    file_name = os.path.expanduser(file_name)

    if not os.path.exists(file_name):
        raise FileNotFoundError(f"The file path '{file_name}' does not exist.")

    # Check extension
    if not file_name.lower().endswith('.gz'):
        return False
    
    # Check file size (gzip files are at least 18 bytes)
    if os.path.getsize(file_name) < 18:
        return False
    
    # Check magic numbers
    with open(file_name, 'rb') as f:
        if f.read(2) != b'\x1f\x8b':
            return False
    
    # Try to open with gzip
    try:
        with gzip.open(file_name, 'rb') as f:
            f.read(1)
        return True
    except gzip.BadGzipFile:
        return False
    

def get_file_delimiter(file_path : str) -> str:
    """
    Determine the delimiter of a given file (compressed or not).
    
    :param file_path (str): Path to the input file.
    :return: Detected delimiter or '\t' by default if not determined.
    """

    file_path = os.path.expanduser(file_path)

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
            dialect = csv.Sniffer().sniff(''.join(lines), delimiters = VALID_DELIMITERS)
            
            return dialect.delimiter
    
    except Exception as e:
        print(f"Error determining delimiter: {str(e)}")
        return '\t'  # Default to tab delimiter