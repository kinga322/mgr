import scanpy as sc
import os
import time
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

pca_times={}
phate_times={}

main_folder_path = "/home/kszyman/Downloads/GSE161529_RAW"

subfolders = [f.path for f in os.scandir(main_folder_path) if f.is_dir()]
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



keys = np.fromiter(pca_times.keys(), dtype=float).reshape(-1, 1)
vals = np.fromiter(pca_times.values(), dtype=float)
reg = LinearRegression().fit(keys, vals)
pred = reg.predict(keys)
plt.plot(keys, pred)
plt.scatter(*zip(*sorted(pca_times.items())))
plt.show()

phate_times = dict(sorted(phate_times.items()))
keys_phate = np.fromiter(phate_times.keys(), dtype=float).reshape(-1, 1)
vals_phate = np.fromiter(phate_times.values(), dtype=float)
poly = PolynomialFeatures(degree=2, include_bias=False)
poly_feat = poly.fit_transform(keys_phate, vals_phate)
poly_reg_model = LinearRegression()
poly_reg_model.fit(poly_feat, vals)
pred2 = poly_reg_model.predict(poly_feat)
plt.scatter(keys_phate, vals_phate)
plt.plot(sorted(keys_phate), pred2)
plt.show()
