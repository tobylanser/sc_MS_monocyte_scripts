# 0. Import

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import celloracle as co


plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

adata_big = sc.read("../../h5ad/v3_allGenes_myeloid_unintegrated_annotated.h5ad")
adata_big.obs['prog'] = adata_big.obs['prog'].astype('category')

celltypes = ["Classical", "Nonclassical", "Intermediate", "mo-DC"] 


for celltype in celltypes:
    adata = adata_big[adata_big.obs['celltypes'].isin([f"{celltype}"])]
    print("Processing: ", f"{celltype}", "prog vs non-prog")
    adata = adata[adata.obs['prog'].isin(['progressor', 'non-progressor'])]
    
    # Random downsampling into 30K cells if the anndata object include more than 30 K cells.
    n_cells_downsample = 30000
    if adata.shape[0] > n_cells_downsample:
        # Let's dowmsample into 30K cells
        sc.pp.subsample(adata, n_obs=n_cells_downsample, random_state=123)
    adata.obs['prog'] = adata.obs['prog'].astype('category')
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="prog", subset=True)
    base_GRN = co.data.load_human_promoter_base_GRN()

    # Instantiate Oracle object
    oracle = co.Oracle()

    # Instantiate Oracle object.
    oracle.import_anndata_as_raw_count(adata=adata,
                                    cluster_column_name="prog",
                                    embedding_name="X_umap.unintegrated")

    # You can load TF info dataframe with the following code.
    oracle.import_TF_data(TF_info_matrix=base_GRN)

    # Alternatively, if you saved the informmation as a dictionary, you can use the code below.
    # oracle.import_TF_data(TFdict=TFinfo_dictionary)

    # Perform PCA
    oracle.perform_PCA()

    # Select important PCs
    plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
    n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
    plt.axvline(n_comps, c="k")
    plt.show()
    print(n_comps)
    n_comps = min(n_comps, 50)

    n_cell = oracle.adata.shape[0]
    print(f"cell number is :{n_cell}")


    k = int(0.025*n_cell)
    print(f"Auto-selected k is :{k}")

    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                        b_maxl=k*4, n_jobs=28)

    # Save oracle object.
    oracle.to_hdf5("../../h5ad/"f"{celltype}""_BP2_unintegrated_annotated2.celloracle.oracle")

    links = oracle.get_links(cluster_name_for_GRN_unit="prog", alpha=10,
                            verbose_level=10)

    links.links_dict.keys()

    # Save Links object.
    links.to_hdf5(file_path="../../h5ad/links.celloracle._"f"{celltype}""_prog_nonprog.celloracle.links")

    links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)

    plt.rcParams["figure.figsize"] = [9, 4.5]

    # Visualize top n-th genes with high scores.
    links.plot_scores_as_rank(cluster="progressors", n_gene=30)

    # Compare GRN score between two clusters
    links.plot_score_comparison_2D(value="betweenness_centrality",
                                cluster1="progressor", cluster2="non-progressor",
                                percentile=98,
                                save = "../../TF_plots/"f"{celltype}""_progressor_vs_nonprog",
                                plt_show = False)

    # Compare GRN score between two clusters
    links.plot_score_comparison_2D(value="eigenvector_centrality",
                                cluster1="progressor", cluster2="non-progressor",
                                percentile=98,
                                save = "../../TF_plots/"f"{celltype}""_progressor_vs_nonprog",
                                plt_show = False)



    # Compare GRN score between two clusters
    links.plot_score_comparison_2D(value="degree_centrality_all",
                                cluster1="progressor", cluster2="non-progressor",
                                save = "../../TF_plots/"f"{celltype}""_progressor_vs_nonprog",
                                plt_show = False)

