import scanpy as sc
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("--input_h5ad", type=str)
parser.add_argument("--output_h5ad", type=str)


args = parser.parse_args()

adata = sc.read_h5ad(args.input_h5ad)

adata2 = adata[adata.obs.batch.isin(['TMA-27', 'TMA-28', 'TMA-29']),:]

sc.pp.scale(adata2, max_value=10)

sc.tl.pca(adata2, svd_solver='arpack')

sc.pp.neighbors(adata2, n_neighbors=10, n_pcs=20)

sc.tl.leiden(adata2, resolution=1.2)

adata2.write_h5ad(args.output_h5ad)

