import scanpy as sc
import anndata

pr1 = sc.read_10x_mtx("/home/kszyman/Downloads/GSE161529_RAW/GSM4909321_mER-MH0068-LN-matrix.mtx/")

sc.pp.pca(pr1)
sc.pp.neighbors(pr1)
sc.tl.leiden(pr1, flavor="igraph", n_iterations=2, directed=False)
sc.external.tl.phate(pr1)
sc.pl.embedding(pr1, "X_phate", color="leiden")
