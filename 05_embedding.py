#!/usr/bin/env python
"""
Find a low-dimensional embedding for a reference data set and project samples onto embedding.
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import umap

from sklearn.decomposition import PCA

if __name__ == '__main__':

    ddir = '/path/to/data/'

    # load data
    reference_gene_counts = pd.read_csv(ddir + 'reference--normalized.tsv', sep='\t')
    reference_hvgs = pd.HDFStore(ddir + 'reference--highly_variable_genes.tsv')

    data = pd.read_csv(ddir + 'mydata--normalized.tsv', sep='\t')
    data = data.set_index(['plate', 'well'])

    # Restrict analysis to genes that
    # 1) are highly variable in the reference, and
    # 2) also present in the data set of interest.
    highly_variable_genes = list(set(reference_hvgs[reference_hvgs['q>0.9']].index.values.tolist()) & set(data.columns.values))
    print(f"Number of HVGs: {len(highly_variable_genes)}")
    reference_subset = reference_gene_counts[highly_variable_genes]
    data_subset = data[highly_variable_genes]

    # log+1 transform
    reference_subset = np.log1p(reference_subset)
    data_subset = np.log1p(data_subset)

    # Reduce dimensionality using PCA.
    # 1) Plot the variance explained by the first 100 components or so.
    pca = PCA(n_components=100).fit(reference_subset)
    fig, (ax1, ax2) = plt.subplots(1,2)
    ax1.plot(pca.explained_variance_ratio_)
    ax2.plot(np.cumsum(pca.explained_variance_ratio_))
    ax1.set_ylabel('Fraction of total variance explained')
    ax2.set_ylabel('Cumulative fraction')
    for ax in (ax1, ax2):
        ax.set_xlabel('Component')
        _, ymax = ax.get_ylim()
        ax.set_ylim(0, ymax)
        ax.grid(True)
    fig.savefig(fdir + 'pca_variance_explained.pdf')

    # 2) Choose a sensible number of components and transform the data.
    pca = PCA(n_components=50, whiten=True).fit(reference_subset) # 10, 20, 50
    reference_pca_transformed = pca.transform(reference_subset)
    data_pca_transformed = pca.transform(data_subset)

    # Reduce dimensionality further via UMAP.
    print("Using UMAP to reduce data down to two dimensions...")
    umap_instance = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.75)
    reference_umap_transformed = umap_instance.fit_transform(reference_pca_transformed)
    data_umap_transformed = umap_instance.transform(data_pca_transformed)

    is_data = np.concatenate([np.zeros((len(reference_umap_transformed)), dtype=np.bool),
                                   np.ones((len(data_umap_transformed)), dtype=np.bool)])
    umap_transformed = np.concatenate([reference_umap_transformed, data_umap_transformed], axis=0)

    fig, ax = plt.subplots(1,1)
    scatter = ax.scatter(umap_transformed[:,0], umap_transformed[:,1],
                         s=is_data.astype(np.float) + 0.1,
                         c=is_data.astype(np.int),
                         rasterized=True)
    handles, labels = scatter.legend_elements(num=None)
    ax.legend(handles, ['reference samples', 'query samples'], bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., ncol=1)
    ax.set_axis_off()
    fig.subplots_adjust(right=0.66)
    fig.savefig(fdir+'umap_reference_vs_data.pdf', dpi=400)
    fig.savefig(fdir+'umap_reference_vs_data.png', dpi=400)

    plt.show()
