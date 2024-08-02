import matplotlib.pyplot as plt
from os import listdir
from anndata import read_h5ad
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from numpy import fromiter


def load_adatas(folder="h5ad"):
    adatas = []
    names = listdir(f"{folder}")
    for name in names:
        adata = read_h5ad(f"{folder}/{name}")
        adatas.append(adata)
    return adatas


def polynomial_regression(x, y, degree):
    poly = PolynomialFeatures(degree=degree, include_bias=False)
    poly_feat = poly.fit_transform(x, y)
    poly_reg_model = LinearRegression()
    poly_reg_model.fit(poly_feat, y)
    pred = poly_reg_model.predict(poly_feat)
    return pred


def scatter_and_line(times_dict, label, degree):
    times_dict = dict(sorted(times_dict.items()))
    keys = fromiter(times_dict.keys(), dtype=float).reshape(-1, 1)
    vals = fromiter(times_dict.values(), dtype=float)
    pred = polynomial_regression(keys, vals, degree)
    plt.scatter(keys, vals, label=label)
    plt.plot(sorted(keys), pred)
    plt.legend()
