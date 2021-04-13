

import pandas as pd
import numpy as np
from collections import Counter

import matplotlib.pyplot as plt
import argparse
import scanpy as sc 

import yaml

from sklearn.metrics import accuracy_score, confusion_matrix
import phenograph
from sklearn.model_selection import StratifiedKFold


import os

print(os.getcwd())

# assert os.getcwd() == "hi"

import sys
sys.path.insert(len(sys.path), 'utils')



from acdc.ACDC.random_walk_classifier import * 
from acdc.ACDC.cell_type_annotation import *


n_neighbor = 10
thres = 0.5

parser = argparse.ArgumentParser()

parser.add_argument("--input_h5ad", type=str)
parser.add_argument("--input_yaml", type=str)
parser.add_argument("--output_assignments", type=str)
parser.add_argument("--method", type=str, choices=['absent', 'no-consider'], default='absent')
parser.add_argument("--cohort",type=str)


args = parser.parse_args()

adata = sc.read_h5ad(args.input_h5ad)


# csv_files = csv_files[0:10]

# dfs = [pd.read_csv(os.path.join(args.input_dir,f), index_col=0) for f in csv_files]
# df = pd.concat(dfs, axis=0)


# df = df.drop(['X', 'Y'], axis=1)


with open(args.input_yaml, "r") as stream:
    marker_dict = yaml.safe_load(stream)
    
marker_dict = marker_dict['cell_types']

features = sorted(
    list(set([l for s in marker_dict.values() for l in s]))
)
classes = list(marker_dict.keys())

adata = adata[:,features]

G,C = [len(features), len(classes)]

m = -1.
if args.method == "no-consider":
    m = 0.

marker_mat = m * np.ones(shape=(G, C))

for g, feature in enumerate(features):
    for c, cell_class in enumerate(classes):
        if feature in marker_dict[cell_class]:
            marker_mat[g, c] = 1.0





table2 = pd.DataFrame(marker_mat.T)
table2.index = classes # + ['Other']
table2.columns = features


idx2ct = [key for idx, key in enumerate(table2.index)]
idx2ct.append('unknown')

ct2idx = {key:idx for idx, key in enumerate(table2.index)}
ct2idx['unknown'] = len(table2.index)
        
ct_score = np.abs(table2.to_numpy()).sum(axis = 1)

y0 = np.zeros(adata.shape[0])




X = adata.X
df = pd.DataFrame(X)
df.columns = adata.var.index
df.index = adata.obs.index

mk_model =  compute_marker_model(df, table2, 0.0)


score = get_score_mat(X, [], table2, [], mk_model)
score = np.concatenate([score, 1.0 - score.max(axis = 1)[:, np.newaxis]], axis = 1)   


ct_index = get_unique_index(X, score, table2, thres)


score = get_score_mat(X, [], table2, [], mk_model)
score = np.concatenate([score, 1.0 - score.max(axis = 1)[:, np.newaxis]], axis = 1)    


ct_index = get_unique_index(X, score, table2, thres)
    

y_pred_index = np.argmax(score, axis = 1)
    
res_c = get_landmarks(X, score, ct_index, idx2ct, phenograph, thres=0.9)

landmark_mat, landmark_label = output_feature_matrix(res_c, [idx2ct[i] for i in range(len(idx2ct))]) 

landmark_label = np.array(landmark_label)

lp, y_pred = rm_classify(X, landmark_mat, landmark_label, n_neighbor)



df_output = pd.DataFrame({'cell_id': adata.obs.index, 'cell_type': y_pred,
'annotator': 'None', 'cohort': args.cohort, 'method': 'acdc-' + args.method}).set_index('cell_id')

df_output.to_csv(args.output_assignments, sep="\t")