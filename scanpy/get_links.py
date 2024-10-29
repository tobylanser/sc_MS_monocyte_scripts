# 0. Import

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import celloracle as co
co.__version__
oracle = co.load_hdf5("../../h5ad/Classical_prog_RR.celloracle.oracle")

# Calculate GRN for each population in "louvain_annot" clustering unit.
# This step may take some time.(~30 minutes)
links = oracle.get_links(cluster_name_for_GRN_unit="prog", alpha=10,
                         verbose_level=10)

links.links_dict.keys()
# Save Links object.
links.to_hdf5(file_path="links.celloracle._Classical_prog_RR.celloracle.links")
links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)
plt.rcParams["figure.figsize"] = [9, 4.5]
# Visualize top n-th genes with high scores.
links.plot_scores_as_rank(cluster="progressors", n_gene=30)

