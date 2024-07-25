import scanpy as sc
import os
import time
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from SIMLR import SIMLR

methods = ["SIMLR", "PCA", "PHATE", "tSNE"]
times = {method: {} for method in methods}
adatas = []

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
    times['PCA'][adata.n_obs] = pca_time
    times["PHATE"][adata.n_obs] = phate_time
    start = time.time()
    sc.tl.tsne(adata)
    end = time.time()
    tsne_time = end - start
    times['tSNE'][adata.n_obs] = tsne_time
    X = adata.X

    X_simlr = SIMLR(X, 5)
    adatas.append(adata)

print(times['PCA'])
# sc.pp.neighbors(pr1)
# sc.tl.leiden(pr1, flavor="igraph", n_iterations=2, directed=False)

# sc.pl.embedding(pr1, "X_phate", color="leiden")


def scatter_and_line(times_dict, label):
    times_dict = dict(sorted(times_dict.items()))
    keys_phate = np.fromiter(times.keys(), dtype=float).reshape(-1, 1)
    vals_phate = np.fromiter(times.values(), dtype=float)
    poly = PolynomialFeatures(degree=2, include_bias=False)
    poly_feat = poly.fit_transform(keys_phate, vals_phate)
    poly_reg_model = LinearRegression()
    poly_reg_model.fit(poly_feat, vals_phate)
    pred2 = poly_reg_model.predict(poly_feat)
    plt.scatter(keys_phate, vals_phate, label=label)
    plt.plot(sorted(keys_phate), pred2)
    plt.legend()


scatter_and_line(times['PCA'], "PCA")
scatter_and_line(times["PHATE"], "PHATE")
scatter_and_line(times['tSNE'], "t-SNE")

for method in methods:
    scatter_and_line(times[method], method)
plt.show()
