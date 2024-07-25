from os import listdir
from anndata import read_h5ad


def load_adatas():
    adatas = []
    names = listdir('./h5ad')
    for name in names:
        adata = read_h5ad(f"./h5ad/{name}")
        adatas.append(adata)
    return adatas
