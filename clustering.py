import scanpy as sc
import matplotlib.pyplot as plt
import time
from sklearn.cluster import KMeans, HDBSCAN
from helpers import load_adatas
import numpy as np


def cluster(adata, method="leiden"):
    start = time.time()
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata)
    if method == "leiden":
        sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)
    elif method == "louvain":
        sc.tl.louvain(adata, flavor="igraph", n_iterations=2)
    elif method == "kmeans":
        X_pca = adata.obsm['X_pca']
        kmeans = KMeans(n_clusters=2, random_state=42).fit(X_pca)
        adata.obs['kmeans'] = kmeans.labels_.astype(str)
    elif method == "HDBSCAN":
        dense = np.asarray(adata.X.todense())
        adata.X = dense
        x_pca = adata.obsm['X_pca']
        hdb = HDBSCAN(min_cluster_size=10)
        hdb.fit(x_pca)
        labels = [str(x) for x in hdb.labels_]
        adata.obs['hdbscan'] = labels
    end = time.time()
    elapsed = end - start
    return adata, elapsed


def plot_embedding(adata, basis, color):
    sc.pl.embedding(adata, basis, color=color)
    plt.show()

