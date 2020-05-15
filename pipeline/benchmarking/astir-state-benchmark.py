# astir: Automated aSsignmenT sIngle-cell pRoteomics

# VSCode tips for python:
# ctrl+` to open terminal
# cmd+P to go to file


import torch
from torch.autograd import Variable
from torch.distributions import Normal
import torch.nn.functional as F
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

import pandas as pd
import numpy as np
import random
import yaml

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

import argparse

parser = argparse.ArgumentParser(description='Run astir')

## argparse treats the options we give it as strings by default

parser.add_argument("expr_csv", help="CSV expression matrix with cells as rows and proteins as columns. First column is cell ID")
parser.add_argument("marker_yaml", help="YAML file of cell markers")
parser.add_argument("output_csv", help="Output CSV of cell assignment probabilities")
parser.add_argument("--epochs", 
                    help="Number of training epochs",
                    type=int,
                    default=100)
parser.add_argument("--learning_rate",
                    help="Learning rate",
                    type=float,
                    default=1e-2)
parser.add_argument("--batch_size",
                    help="Batch size",
                    type=int,
                    default=1024)
parser.add_argument("--random_seed",
                    help="Random seed",
                    type=int,
                    default=1234)
parser.add_argument("--print_loss_every",
                    help="Print the loss this many iterations",
                    type=int,
                    default=1234)

args = parser.parse_args()

## Set seeds
torch.manual_seed(args.random_seed)
torch.cuda.manual_seed(args.random_seed)
np.random.seed(args.random_seed)
random.seed(args.random_seed)
torch.backends.cudnn.benchmark = False
torch.backends.cudnn.deterministic = True

# https://github.com/pytorch/pytorch/issues/7068
def _init_fn(worker_id):
    return args.random_seed

# print((args.expr_csv, args.marker_yaml, args.output_csv)) 
# -2619.1226096122696
# -2645.941610884511
# -2664.729080970788

## Load the gene expression data
df_gex = pd.read_csv(args.expr_csv, index_col = 0)

core_names = list(df_gex.index)
gene_names = list(df_gex.columns)

## Load the marker info data
with open(args.marker_yaml, 'r') as stream:
    markers_states = yaml.safe_load(stream)
    
states = markers_states['cell_states']

state_names = list(states.keys())
genes = [l for s in states.values() for l in s]
genes = list(set(genes)) # Make unique

## Infer dimensions
N = df_gex.shape[0]
G = len(genes)
C = len(state_names)

## Construct marker matrix
state_mat = np.zeros(shape = (G,C))

for g in range(len(genes)):
    for ct in range(len(state_names)):
        gene = genes[g]
        state = state_names[ct]
        
        if gene in states[state]:
            state_mat[g,ct] = 1

Y_np = df_gex[genes].to_numpy()
Y_np = Y_np / (Y_np.std(0))

# ## Dataset class: for loading IMC datasets
# class IMCDataSet(Dataset):
    
#     def __init__(self, Y_np):
             
#         self.Y = torch.from_numpy(Y_np)
#         X = StandardScaler().fit_transform(Y_np)
#         self.X = torch.from_numpy(X)
    
#     def __len__(self):
#         return self.Y.shape[0]
    
#     def __getitem__(self, idx):
#         return self.Y[idx,:], self.X[idx,:]

# ## Recognition network
# class RecognitionNet(nn.Module):
#     def __init__(self, C, G, h=6):
#         super(RecognitionNet, self).__init__()
#         self.hidden_1 = nn.Linear(G, h).double()
#         self.hidden_2 = nn.Linear(h, C+1).double()

#     def forward(self, x):
#         x = self.hidden_1(x)
#         x = F.relu(x)
#         x = self.hidden_2(x)
#         x = F.softmax(x, dim=1)
#         return x

## Parameters and initialization

# log_sigma_init = np.log(Y_np.std(0))
log_sigma_init = np.log(Y_np.std()).reshape(1,-1)
mu_init = Y_np.mean(0).reshape(1, -1)

# def get_pc_init(pathway, df_gex, states):
#     gex_pt = df_gex[ states[ pathway ] ].to_numpy()
#     # gex_pt = np.log(gex_pt)
#     pca = PCA(n_components=1, whiten=True)
#     pca.fit(gex_pt)
#     a = pca.transform(gex_pt)
#     if pca.components_.mean() < 0:
#         a = -a
#     return a
    
# log_alpha_init = np.log(np.stack([get_pc_init(p, df_gex, states) for p in state_names], axis=0).reshape(C, N).T)
alpha_init = np.zeros((N,C))

## Initialize log_alpha


# log_sigma = Variable(torch.from_numpy(log_sigma_init.copy()), requires_grad = True)
log_sigma = Variable(torch.from_numpy(log_sigma_init.copy()), requires_grad = True)
mu = Variable(torch.from_numpy(mu_init.copy()), requires_grad = True)
alpha = Variable(torch.from_numpy(alpha_init.copy()), requires_grad = True)
log_beta = Variable(torch.zeros((C,G)), requires_grad = True)

rho = torch.from_numpy(state_mat.T).double()
Y = torch.from_numpy(Y_np)

## Construct optimizer
optimizer = torch.optim.Adam([log_sigma, mu, alpha, log_beta], lr=args.learning_rate)

## Declare pytorch forward fn
def forward(log_sigma, mu, alpha, log_beta, rho, Y):
    mean = mu + torch.matmul(alpha, rho * torch.exp(log_beta))
    dist = Normal(mean, torch.exp(log_sigma).reshape(1, -1))
    
    log_p_y = dist.log_prob(Y)
    prior_alpha = Normal(torch.zeros(1),  1 * torch.ones(1)).log_prob(alpha)
    prior_beta = Normal(torch.zeros(1), 1 * torch.ones(1)).log_prob(log_beta)
    prior_sigma = Normal(torch.zeros(5), 0.5 * torch.ones(1)).log_prob(log_sigma)
    
    loss = log_p_y.sum() + prior_alpha.sum() + prior_beta.sum() + prior_sigma.sum() 
    
    return -loss


## Make dataloader
# dset = IMCDataSet(Y_np)
# dataloader = DataLoader(dset, 
#                         batch_size=min(args.batch_size,N), shuffle=True,
#                         num_workers=0,
#                         worker_init_fn=_init_fn,
#                         pin_memory=True)


## Run training loop
epochs = args.epochs
losses = np.empty(epochs)

for it in range(epochs):
    optimizer.zero_grad()
    loss = forward(log_sigma, mu, alpha, log_beta, rho, Y)
    loss.backward()
    optimizer.step()

    l = loss.detach().numpy()
    losses[it] = loss
    
    if it % args.print_loss_every == 0:       
        print(l)

## Save output
b = (log_beta).detach().numpy()
print(b)

g = (alpha).detach().numpy()

assignments = pd.DataFrame(g)
assignments.columns = state_names
assignments.index = core_names

assignments.to_csv(args.output_csv)

losses = pd.DataFrame({'iteration': np.arange(epochs), 'loss': losses})
losses.to_csv('loss.csv')


print("Done!")