from dimred import run_method
from clustering import cluster
import sklearn.metrics as metrics
from helpers import load_adatas
import matplotlib.pyplot as plt
import numpy as np


emb_dict = {"PCA": "X_pca", "PHATE": "X_phate", "tSNE": "x_tsne"}


def compare_dimred(adata, methods_list):
    for method in methods_list:
        run_method(adata, method)
        cluster(adata, emb_dict[method], "leiden")


adatas = load_adatas()
calinski_harabasz_dict = {"X_tsne": [], "X_phate": [], "X_pca": [], "X_umap": []}
for adata in adatas:
    for rep in adata.obsm:
        calinski_harabasz = metrics.calinski_harabasz_score(adata.obsm[rep], adata.obs["leiden"])
        calinski_harabasz_dict[rep].append(calinski_harabasz)
metrics.davies_bouldin_score()
metrics.silhouette_score()
plt.boxplot(np.array(list(calinski_harabasz_dict.values())))
plt.show()
