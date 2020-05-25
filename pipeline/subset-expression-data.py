## Subset a directory of CSV files to a pre-defined cell list

import os
import pandas as pd
import numpy as np
import loompy
import argparse

parser = argparse.ArgumentParser(description='Convert a directory of CSV files to a loom file and add a batch variable')

## argparse treats the options we give it as strings by default

parser.add_argument("--input_dir", help="Directory containing expression data in CSV format", type=str)
parser.add_argument("--cell_labels", help="CSV file with column 'cell' ", type=str)
parser.add_argument("--output_dir", help="Output directory to write reduced expression data to ", type=str)

args = parser.parse_args()

# Step 1: get cell IDs we want to keep
cells = pd.read_csv(args.cell_labels)
to_select = list(cells.cell)

# Step 2: find csv files

files = [f for f in os.listdir(args.input_dir) if f.endswith(".csv")]

# Step 3: iterate through expression files and save if non-zero
for f in files:
    # read in 
    df = pd.read_csv(os.path.join(args.input_dir, f), index_col=0)
    if any([x in df.index for x in to_select]):
        ts = [i for i in to_select if i in list(df.index)]
        df = df.loc[ts]
        df.to_csv(os.path.join(args.output_dir, f))

print("Done")