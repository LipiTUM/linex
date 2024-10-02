from LipidNetworkInference import (
    LipidNetwork,
    SumSpecies
)
import os
import re
import pandas as pd
import numpy as np
import time


def timer(func):
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        returns = func(*args, **kwargs)
        end = time.perf_counter()
        print(f"{func.__name__} ran {end - start:1.2f} seconds")
        return returns
    return wrapper


@timer
def init_network(**kwargs):
    return LipidNetwork(**kwargs)


# Simple test
# ls = pd.Series(["PC 34:2", "PE 32:2", "PC 32:2", "LPC 16:0"])
# ld = pd.DataFrame(np.random.randn(ls.size, 20),
#                   index=ls)
# print("testing class from species list")
# ln_ = LipidNetwork(lipid_species=ls)
# print("testing class from lipid data")
# ln = LipidNetwork(lipid_data=ld, sample_axis="columns")
# print("computing network")
# ln.compute_network()
# print("partial correlations")
# ln.compute_correlations()
# print("correlations")
# ln.compute_partial_correlations()
# print("plotting network")
# ln.plot_static_network(edge_colour_attr="correlations",
#                        node_colour_attr="lipid_class")

# More complex test using the data from Purdy et al. 2019
path = "/home/nik/PhD/LipidNetworkInference/ValidationStudies/Ansari2016NatureSD"
ansari = pd.read_csv(os.path.join(path, "lynx_format_data.tsv"),
                     sep="\t")

print("testing class from lipid data")
ansari = ansari.loc[~ansari.index.duplicated(), :]
ln = init_network(lipid_data=ansari,
                  sample_axis="columns",
                  lipid_level=SumSpecies,
                  names_in_convention=True)
print("computing network")
ln.compute_network()
print("partial correlations")
ln.compute_correlations()
print("correlations")
ln.compute_partial_correlations()
print("plotting network")
ln.plot_static_network(edge_colour_attr="correlations",
                       node_colour_attr="lipid_class")
ln.plot_interactive_network(edge_colour_attr="correlations",
                            node_colour_attr="lipid_class")

print("Sucessfully finished test!")
