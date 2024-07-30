from helpers import load_adatas, scatter_and_line
from dimred import run_method
import matplotlib.pyplot as plt


def run_all(methods_list):
    adatas_list = load_adatas("./h5ad")
    times_dict = {method: {} for method in methods_list}
    memory_dict = {method: {} for method in methods_list}

    for i, adata in enumerate(adatas_list):
        print(i)
        for method in methods_list:
            adata, time_mes, memory = run_method(adata, method)
            times_dict[method][adata.n_obs] = time_mes
            memory_dict[method][adata.n_obs] = memory

    plt.show()

    return times_dict, memory_dict, adatas_list


methods = ["PCA", "Diffmap"]
times, memory, adatas = run_all(methods)
regression_degree = 1
for method in methods:
    scatter_and_line(times[method], method, regression_degree)
    scatter_and_line(memory[method], method, regression_degree)
plt.show()
