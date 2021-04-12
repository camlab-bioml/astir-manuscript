## Convert a directory of CSV files to a loom file and add a batch variable

import os
import pandas as pd
import numpy as np
import argparse
from anndata import AnnData
import scanpy as sc

parser = argparse.ArgumentParser(description='Convert a directory of CSV files to a loom file and add a batch variable')

## argparse treats the options we give it as strings by default

parser.add_argument("input_dir", help="Directory containing CSV files", type=str)
parser.add_argument("output_h5ad", help="Path to output loom file", type=str)

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


obs = pd.DataFrame({'cell_id': list(df_gex.index),
'batch': batch})

adata = AnnData(df_gex.values, obs)
adata.var_names = list(df_gex.columns)
adata.obs_names = list(df_gex.index)


print("Writing h5ad")

adata.write(args.output_h5ad)

print("Done")