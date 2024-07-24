import scanpy as sc
import os
import time
import matplotlib.pyplot as plt

pca_times={}
phate_times={}

subfolders = [f.path for f in os.scandir("/home/kszyman/Downloads/GSE161529_RAW") if f.is_dir()]
for subdir in subfolders:
    adata = sc.read_10x_mtx(subdir)
    sc.pp.filter_cells(adata, min_genes=10)
    start = time.time()
    sc.external.tl.phate(adata)
    end = time.time()
    phate_time = end - start
    start = time.time()
    sc.pp.pca(adata)
    end = time.time()
    pca_time = end - start
    pca_times[adata.n_obs] = pca_time
    phate_times[adata.n_obs] = phate_time

print(pca_times)
# sc.pp.neighbors(pr1)
# sc.tl.leiden(pr1, flavor="igraph", n_iterations=2, directed=False)

# sc.pl.embedding(pr1, "X_phate", color="leiden")
plt.scatter(*zip(*sorted(pca_times.items())))
plt.show()
plt.scatter(*zip(*sorted(phate_times.items())))
plt.show()