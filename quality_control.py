#!/usr/bin/env python
"""
Quality control of gene count matrices.
"""

import functools
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


from utils import _handle_io


import matplotlib as mpl
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'


@_handle_io
def quality_control(gene_abundance,
                    minimum_total_reads =  100000, # absolute number
                    maximum_total_reads = 4000000, # absolute number
                    minimum_genes_expressed = 0.05, # fraction of total unique genes
                    maximum_mitochondrial_reads = 0.1, # fraction of total reads
                    minimum_samples = 0.02, # fraction of total samples
                    show = False,
                    fdir = None,
):
    # TODO: suppress plotting if show is False and fdir is None

    figure_names = [
        'total_reads_per_sample.pdf',
        'fraction_of_genes_expressed_per_sample.pdf',
        'fraction_of_mitochondrial_genes.pdf',
        'fraction_of_genes_expressed_per_sample.pdf',
    ]
    if fdir:
        figure_paths = [fdir + name for name in figure_names]
    else:
        figure_paths = [None for _ in figure_names]

    sufficient_reads       = qc_minimum_total_reads(gene_abundance,           threshold=minimum_total_reads,         figure_path=figure_paths[0])
    reasonable_reads       = qc_maximum_total_reads(sufficient_reads,         threshold=maximum_total_reads,         figure_path=figure_paths[0])
    sufficient_genes       = qc_fraction_of_genes_expressed(reasonable_reads, threshold=minimum_genes_expressed,     figure_path=figure_paths[1])
    sufficiently_healthy   = qc_mitochondrial_read_content(sufficient_genes,  threshold=maximum_mitochondrial_reads, figure_path=figure_paths[2])
    sufficiently_expressed = qc_gene_expression(sufficiently_healthy,         threshold=minimum_samples,             figure_path=figure_paths[-1])

    report_loss(gene_abundance, sufficiently_expressed, "\nAfter QC:")
    report_lost_samples(gene_abundance, sufficiently_expressed)

    if show:
        plt.ion()
        plt.show()
        input("Press any key to close all figures.")
        plt.close('all')

    return sufficiently_expressed


def before_after_distplot(before, after,
                          threshold=None,
                          xlabel = None,
                          ylabel = None,
                          figure_path=None):

    fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)
    sns.distplot(before, kde_kws=dict(clip=(0, np.max(before))), ax=ax1)
    sns.distplot(after,  kde_kws=dict(clip=(0, np.max(after))),  ax=ax2)

    if threshold:
        ax1.axvline(threshold, c='gray', ls='--')

    ax1.set_ylabel(ylabel)
    ax1.set_xlabel(xlabel)
    ax2.set_xlabel(xlabel)
    ax1.set_title('Before')
    ax2.set_title('After')

    fig.tight_layout()

    if figure_path:
        fig.savefig(figure_path)


def _report_loss(preamble=None):
    def outer(func):
        @functools.wraps(func)
        def inner(data, *args, **kwargs):
            remaining_data = func(data, *args, **kwargs)
            report_loss(data, remaining_data, preamble)
            return remaining_data
        return inner
    return outer


def report_loss(data, remaining_data, preamble=None):
    if preamble:
        print(preamble)
    total_samples, total_genes = data.shape
    remaining_samples, remaining_genes = remaining_data.shape
    percent_samples = 100 * remaining_samples / total_samples
    percent_genes = 100 * remaining_genes / total_genes
    print(f'{remaining_samples} / {total_samples} ({percent_samples:.1f}%) of samples remain.')
    print(f'{remaining_genes} / {total_genes} ({percent_genes:.1f}%) of genes remain.')


def report_lost_samples(data, remaining_data):
    data_samples = set(data.index.tolist())
    remaining_samples = set(remaining_data.index.tolist())
    lost_samples = data_samples - remaining_samples

    print('\nThe following samples did not pass QC:')
    for sample in sorted(lost_samples):
        print(f'{sample}')


@_report_loss(preamble="\nEliminating samples with low total reads:")
def qc_minimum_total_reads(data, threshold=100000, figure_path=None):
    """Eliminate samples with too low total reads."""

    total_reads_before = data.sum(axis=1)
    remaining_data = data[total_reads_before > threshold]
    total_reads_after = remaining_data.sum(axis=1)

    before_after_distplot(
        before      = total_reads_before,
        after       = total_reads_after,
        threshold   = threshold,
        xlabel      = 'Total reads',
        ylabel      = 'Frequency of samples',
        figure_path = figure_path,
    )

    return remaining_data


@_report_loss(preamble="\nEliminating samples with high total reads:")
def qc_maximum_total_reads(data, threshold=1000000, figure_path=None):
    """Eliminate samples with too low total reads."""

    total_reads_before = data.sum(axis=1)
    remaining_data = data[total_reads_before < threshold]
    total_reads_after = remaining_data.sum(axis=1)

    before_after_distplot(
        before      = total_reads_before,
        after       = total_reads_after,
        threshold   = threshold,
        xlabel      = 'Total reads',
        ylabel      = 'Frequency of samples',
        figure_path = figure_path,
    )

    return remaining_data


@_report_loss(preamble="\nEliminating samples with low number of expressed genes:")
def qc_fraction_of_genes_expressed(data, threshold=0.05, figure_path=None):
    """Eliminate samples with low number of expressed genes."""

    fraction_genes_expressed_before = (data > 0).mean(axis=1)
    remaining_data = data[fraction_genes_expressed_before > threshold]
    fraction_genes_expressed_after = (remaining_data > 0).mean(axis=1)

    before_after_distplot(
        before      = fraction_genes_expressed_before,
        after       = fraction_genes_expressed_after,
        threshold   = threshold,
        xlabel      = 'Fraction of genes expressed',
        ylabel      = 'Frequency of samples',
        figure_path = figure_path,
    )

    return remaining_data


@_report_loss(preamble="\nEliminating samples with high mitochondrial gene expression:")
def qc_mitochondrial_read_content(data, threshold=0.05, figure_path=None):
    """Eliminate samples with a high proportion of reads mapping to mitochondrial transcripts.
    These cells are thought to undergo apoptosis."""

    is_mitochondrial_gene = np.array([True if 'mt-'==name[:3] else False \
                                      for name in data.columns.values], dtype=bool)

    total_mitochondrial_reads = data.loc[:, is_mitochondrial_gene].sum(axis=1)
    total_reads = data.sum(axis=1)
    fraction_mitochondrial_before = total_mitochondrial_reads / total_reads

    remaining_data = data[fraction_mitochondrial_before < threshold]

    total_mitochondrial_reads = remaining_data.loc[:, is_mitochondrial_gene].sum(axis=1)
    total_reads = remaining_data.sum(axis=1)
    fraction_mitochondrial_after = total_mitochondrial_reads / total_reads

    before_after_distplot(
        before      = fraction_mitochondrial_before,
        after       = fraction_mitochondrial_after,
        threshold   = threshold,
        xlabel      = 'Fraction of mitochondrial reads',
        ylabel      = 'Frequency of samples',
        figure_path = figure_path,
    )

    return remaining_data


@_report_loss(preamble="\nEliminating genes that are only present in few samples:")
def qc_gene_expression(data, threshold=0.02, figure_path=None):
    """Eliminate genes that are only present in very few samples."""

    fraction_samples_expressed_before = (data > 0).mean(axis=0)
    remaining_data = data.loc[:, fraction_samples_expressed_before > threshold]
    fraction_samples_expressed_after = (remaining_data > 0).mean(axis=0)

    before_after_distplot(
        before      = fraction_samples_expressed_before,
        after       = fraction_samples_expressed_after,
        threshold   = threshold,
        xlabel      = 'Fraction of samples with GOI',
        ylabel      = 'Frequency of genes',
        figure_path = figure_path,
    )

    return remaining_data


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Quality control for gene abundance data from SmartSeq2 sequencing aligned with kallisto.")
    parser.add_argument(
        'input',
        help="/path/to/gene_abundance.tsv"
    )
    parser.add_argument(
        '-o', '--output',
        help="/path/to/qced_gene_abundance.tsv",
        default="./qced_gene_abundance.tsv")
    parser.add_argument("-s", "--show", action="store_true", help="Plot the output figures of the script.")
    args = parser.parse_args()

    quality_control(args.input, args.output, show=args.show)
