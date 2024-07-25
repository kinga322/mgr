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


def run_method(adata, method):
    start_time = time.time()
    if method == "SIMLR":
        x = adata.X
        simlr = SIMLR(x, 5)
        adata.obsm['X_simlr'] = simlr
    elif method == "PCA":
        sc.pp.pca(adata)
    elif method == "PHATE":
        sc.external.tl.phate(adata)
    elif method == "tSNE":
        sc.tl.tsne(adata)
    end_time = time.time()
    return adata, end_time - start_time


def read_adata(folder):
    adata = sc.read_10x_mtx(subdir)
    return adata


for subdir in subfolders:
    adata = read_adata(subdir)
    sc.pp.filter_cells(adata, min_genes=10)
    phate_time = run_method(adata, "PHATE")
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
    start = time.time()
    X_simlr = SIMLR(X, 5)
    end = time.time()
    simlr_time = end - start
    times['SIMLR'][adata.n_obs] = simlr_time
    adatas.append(adata)

print(times['PCA'])
# sc.pp.neighbors(pr1)
# sc.tl.leiden(pr1, flavor="igraph", n_iterations=2, directed=False)

# sc.pl.embedding(pr1, "X_phate", color="leiden")


def scatter_and_line(times_dict, label):
    times_dict = dict(sorted(times_dict.items()))
    keys = np.fromiter(times_dict.keys(), dtype=float).reshape(-1, 1)
    vals = np.fromiter(times_dict.values(), dtype=float)
    poly = PolynomialFeatures(degree=2, include_bias=False)
    poly_feat = poly.fit_transform(keys, vals)
    poly_reg_model = LinearRegression()
    poly_reg_model.fit(poly_feat, vals)
    pred2 = poly_reg_model.predict(poly_feat)
    plt.scatter(keys, vals, label=label)
    plt.plot(sorted(keys), pred2)
    plt.legend()


scatter_and_line(times['PCA'], "PCA")
scatter_and_line(times["PHATE"], "PHATE")
scatter_and_line(times['tSNE'], "t-SNE")

for method in methods:
    scatter_and_line(times[method], method)
plt.show()
