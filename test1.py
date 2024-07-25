import scanpy as sc
import scanpy.external as sce
import os
import time
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
# from SIMLR import SIMLR


def run_method(adata, method):
    start_time = time.time()

    """    
    if method == "SIMLR":
    x = adata.X
    simlr = SIMLR(x, 5)
    adata.obsm['X_simlr'] = simlr
    """

    if method == "PCA":
        sc.pp.pca(adata)
    elif method == "PHATE":
        sce.tl.phate(adata)
    elif method == "tSNE":
        sc.tl.tsne(adata)
    elif method == "UMAP":
        sc.tl.umap(adata)
    elif method == "Diffmap":
        sc.pp.neighbors(adata)
        sc.tl.diffmap(adata)
    end_time = time.time()
    elapsed = end_time - start_time
    print(elapsed)
    return adata, elapsed


def read_adata(folder):
    adata = sc.read_10x_mtx(folder)
    return adata


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


def run_all(main_folder, methods_list):
    adatas_list = []
    times_dict = {method: {} for method in methods_list}
    subfolders = [f.path for f in os.scandir(main_folder) if f.is_dir()]
    for subdir in subfolders:
        adata = read_adata(subdir)
        sc.pp.filter_cells(adata, min_genes=10)
        for method in methods_list:
            adata, time_mes = run_method(adata, method)
            times_dict[method][adata.n_obs] = time_mes

        adatas_list.append(adata)

    for method in methods_list:
        scatter_and_line(times_dict[method], method)
    plt.show()

    return times_dict, adatas_list


methods = ["PCA", "PHATE", "tSNE", "Diffmap"]
main_folder_path = "/home/kszyman/Downloads/GSE161529_RAW"
times, adatas = run_all(main_folder_path, methods)
