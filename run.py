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


methods = ["PHATE", "UMAP"]
times, memory, adatas = run_all(methods)
regression_degree = 2
for method in methods:
    scatter_and_line(times[method], method, regression_degree)
plt.show()
for method in methods:
    scatter_and_line(memory[method], method, regression_degree)
plt.show()


# times_new.sort_index(inplace=True)
"""for col in times_new.columns:
    plt.plot(times_new[col])
plt.show()    """
