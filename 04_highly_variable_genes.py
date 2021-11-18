#!/usr/bin/env python
"""
Determine highly variable genes.

Based on:
http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.formula.api as smf


def highly_variable_genes(gene_counts, quantile=0.95, show=True, figure_directory=None):
    means = gene_counts.mean(axis=0)
    variances = gene_counts.var(axis=0)
    normalized_variances = variances / means**2

    output = pd.DataFrame(

        {
            'LM': np.log(means),
            'LV' : np.log(normalized_variances),
        },
        index=gene_counts.columns.values
    )

    model = smf.quantreg('LV ~ LM', output)

    if isinstance(quantile, float):
        quantiles = [quantile]
    else:
        quantiles = quantile

    for quantile in quantiles:
        fit = model.fit(q=quantile)
        a = fit.params['LM']
        b = fit.params['Intercept']
        mask = is_above_line(output[['LM', 'LV']], a, b)
        output = pd.concat([output, pd.Series(mask, index=gene_counts.columns.values, name=f'q>{quantile}')], axis=1)

        if show:
            fig, ax = plt.subplots(1, 1)
            ax.scatter(output['LM'], output['LV'], c='k', s=0.1, alpha=0.5)
            ax.set_xlim(1.2*np.percentile(output['LM'], [0.1, 99.9]))
            ax.set_ylim(1.2*np.percentile(output['LV'], [0.1, 99.9]))
            ax.set_xlabel('log(mean)')
            ax.set_ylabel('log(normalized variance)')
            x = np.linspace(output['LM'].min(), output['LM'].max(), 100)
            y = a*x + b

            ax.plot(x, y, '--', lw=3, label=f'{model.formula}, q = {quantile:.2f}: n = {np.sum(mask)}')
            ax.scatter(output['LM'][mask], output['LV'][mask], s=0.2, rasterized=True)
            ax.legend(loc='lower left')
            fig.tight_layout()

            if figure_directory:
                fig.savefig(figure_directory + f'quantile_regression_{quantile:.2f}.pdf')

    return output


def is_above_line(points, gradient, intercept):
    # https://stackoverflow.com/a/45769740/2912349
    a = np.array([0, intercept])
    b = np.array([1, intercept + gradient])
    return np.cross(points-a, b-a) < 0


def get_smooth_estimate_of_quantile(x, y, q, fraction=0.1):
    # order points by independent measure
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    # determine number of points to smooth over
    total_points = len(x)
    delta = int(total_points * fraction / 2)

    xq = [np.mean(x[ii-delta:ii+delta]) for ii in range(delta, total_points-delta)]
    yq = [np.quantile(y[ii-delta:ii+delta], q) for ii in range(delta, total_points-delta)]

    return xq, yq


if __name__ == '__main__':

    ddir = '/home/paul/wdir/transcriptomics/data/gene_counts/'
    file_base_name = 'mydata'
    file_path = ddir + file_base_name + '--normalized.tsv'
    gene_counts = pd.read_csv(file_path, sep='\t', index_col=['plate', 'well'])
    df = highly_variable_genes(gene_counts, show=True)
    df.to_csv(ddir + file_base_name + '--highly_variable_genes.tsv', sep='\t')
    plt.show()
