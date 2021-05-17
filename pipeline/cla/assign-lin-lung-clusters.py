import scanpy as sc
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("--input_h5ad", type=str)
parser.add_argument("--output_h5ad", type=str)


args = parser.parse_args()

adata = sc.read_h5ad(args.input_h5ad)

exprs = adata.X.copy()

sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)

sc.tl.leiden(adata, resolution=1.2)

adata.X = exprs

adata.write_h5ad(args.output_h5ad)

