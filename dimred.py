import scanpy as sc
import scanpy.external as sce
import time
from scvi.model import SCVI
from simlr_py3 import SIMLR_LARGE
from sklearn.decomposition import FastICA
import tracemalloc


def run_method(adata, method, n_neighbors=5):
    start_time = time.time()
    tracemalloc.start()
    if method == "SIMLR":
        x = adata.X
        simlr = SIMLR_LARGE(x, 5)
        adata.obsm['X_simlr'] = simlr
    if method == "PCA":
        sc.pp.pca(adata)
    elif method == "PHATE":
        sce.tl.phate(adata)
    elif method == "tSNE":
        sc.tl.tsne(adata)
    elif method == "UMAP":
        sc.pp.neighbors(adata, n_neighbors=n_neighbors)
        sc.tl.umap(adata)
    elif method == "scVI":
        SCVI.setup_anndata(adata)
        model = SCVI(adata)
        model.train()
        adata.obsm["X_scvi"] = model.get_latent_representation()
    elif method == "ICA":
        x = adata.X.toarray()
        res = FastICA().fit_transform(x)
        adata.obsm["X_ica"] = res
    else:
        raise ValueError("Method must be SIMLR, PCA, t-SNE, PHATE, ICA, UMAP or scVI")
    end_time = time.time()
    elapsed = end_time - start_time
    print(method)
    memory = max(tracemalloc.get_traced_memory())
    print("Memory:" + str(memory))
    tracemalloc.stop()
    print("Time:" + str(elapsed))
    return adata, elapsed, memory
