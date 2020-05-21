## Convert a directory of CSV files to a loom file and add a batch variable

import os
import pandas as pd
import numpy as np
import loompy
import argparse

parser = argparse.ArgumentParser(description='Convert a directory of CSV files to a loom file and add a batch variable')

## argparse treats the options we give it as strings by default

parser.add_argument("input_dir", help="Directory containing CSV files", type=str)
parser.add_argument("output_loom", help="Path to output loom file")

# parser.add_argument("--random_seed",
#                     help="Random seed",
#                     type=int,
#                     default=1234)

args = parser.parse_args()

print("Parsing CSV files to single pd.DataFrame")
csv_files = [f for f in os.listdir(args.input_dir) if f.endswith(".csv")]

dfs = [pd.read_csv(os.path.join(args.input_dir,f), index_col=0) for f in csv_files]
df_gex = pd.concat(dfs, axis=0)

batch = list(np.concatenate([np.repeat(csv_files[i].replace(".csv", ""), dfs[i].shape[0]) for i in range(len(dfs))], axis=0))

row_attrs = {"protein": list(df_gex.columns)}
col_attrs = {"cell_name": list(df_gex.index),
            "batch": batch}

print("Writing loom")

loompy.create(args.output_loom, df_gex.to_numpy().T, row_attrs, col_attrs)

print("Done")