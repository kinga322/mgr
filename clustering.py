import scanpy as sc
from helpers import load_adatas

adatas = load_adatas()

for adata in adatas:
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)
    sc.pl.embedding(adata, "X_phate", color="leiden")
