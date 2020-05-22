## Subset basel dataset to the 15k patient samples

import os
import pandas as pd
import numpy as np
import loompy
import argparse

parser = argparse.ArgumentParser(description='Convert a directory of CSV files to a loom file and add a batch variable')

## argparse treats the options we give it as strings by default

parser.add_argument("cell_assignments", help="CSV file containing astir assignments", type=str)
parser.add_argument("cell_labels", help="CSV file with column 'cell' ", type=str)
parser.add_argument("output_csv", help="Output_csv file ", type=str)

args = parser.parse_args()

g = pd.read_csv(args.cell_assignments, index_col=0)
cells = pd.read_csv(args.cell_labels)
to_select = list(cells.cell)

g2 = g.loc[to_select]
g2.to_csv(args.output_csv)

print("Done")