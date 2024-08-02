import scanpy as sc
import matplotlib.pyplot as plt
import time
import numpy as np
import c_index
from sklearn.cluster import KMeans, HDBSCAN
from validclust.validclust import ValidClust
from gapstatistics.gapstatistics import GapStatistics


def cluster(adata, embedding, method="leiden"):
    start = time.time()
    sc.pp.neighbors(adata, use_rep=embedding, key_added=embedding)
    if method == "leiden":
        sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)
    elif method == "louvain":
        sc.tl.louvain(adata, flavor="igraph", n_iterations=2)
    elif method == "kmeans":
        x_pca = adata.obsm['X_pca']
        kmeans = KMeans(n_clusters=2, random_state=42).fit(x_pca)
        adata.obs['kmeans'] = kmeans.labels_.astype(str)
    elif method == "HDBSCAN":
        # dense = np.asarray(adata.X.todense())
        # adata.X = dense
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

"""

metrics.fowlkes_mallows_score()
GapStatistics()
ValidClust()
"""