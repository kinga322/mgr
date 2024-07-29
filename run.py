from helpers import load_adatas, scatter_and_line
from dimred import run_method
import matplotlib.pyplot as plt


def run_all(methods_list, regression_degree):
    adatas_list = load_adatas("./h5ad")
    times_dict = {method: {} for method in methods_list}

    for i, adata in enumerate(adatas_list):
        print(i)
        for method in methods_list:
            adata, time_mes = run_method(adata, method)
            times_dict[method][adata.n_obs] = time_mes

    for method in methods_list:
        scatter_and_line(times_dict[method], method, regression_degree)
    plt.show()

    return times_dict, adatas_list


methods = ["PCA", "Diffmap"]
times, adatas = run_all(methods, 1)
