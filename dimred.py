import scanpy as sc
import scanpy.external as sce
import time
import matplotlib.pyplot as plt
from helpers import load_adatas, scatter_and_line
# from SIMLR import SIMLR

# assumes there's "h5ad" directory in current workdir containing h5ad files


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


def run_all(methods_list):
    adatas_list = load_adatas()
    times_dict = {method: {} for method in methods_list}

    for adata in adatas_list:
        print(adata.filename)
        for method in methods_list:
            adata, time_mes = run_method(adata, method)
            times_dict[method][adata.n_obs] = time_mes

    for method in methods_list:
        scatter_and_line(times_dict[method], method)
    plt.show()

    return times_dict, adatas_list


methods = ["PCA", "PHATE", "tSNE", "Diffmap"]
times, adatas = run_all(methods)
