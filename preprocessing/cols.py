import yaml
import numpy as np

class Cols:
    chromosome = 'chromosome'
    position = 'position'
    allele1 = 'allele1'
    allele2 = 'allele2'
    beta = 'beta'
    se = 'se'
    pval = 'pval'
    odds_ratio = 'odds_ratio'

    def __init__(self, yaml_file):
        try:
            with open(yaml_file, 'r') as file:
                self.config = yaml.safe_load(file)
        except FileNotFoundError:
            print(f"Error: The file {yaml_file} was not found.")
            self.config = {}
        except yaml.YAMLError as exc:
            print(f"Error parsing YAML file: {exc}")
            self.config = {}

        # Set class attributes based on the YAML configuration
        self.columns = {}
        for column, dtype in self.config.get('columns', {}).items():
            try:
                # Handle string type separately
                if dtype == 'str':
                    self.columns[column] = str
                else:
                    # Use numpy data types for consistency with data processing libraries
                    self.columns[column] = getattr(np, dtype)
                setattr(self, column, column)
            except AttributeError:
                print(f"Error: {dtype} is not a valid numpy data type.")
                self.columns[column] = None

    def __repr__(self):
        return f"Cols(columns={self.columns})"

    def get_dtype(self, column_name):
        """Get the data type for a given column name."""
        return self.columns.get(column_name, None)

if __name__ == "__main__":
    cols = Cols('header.yaml')

    for column, dtype in cols.columns.items():
        print(f"{column}: {dtype}")

    # print(cols.get_dtype('chromosome'))