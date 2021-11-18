#!/usr/bin/env python
"""
Merge samples within a data set that have been run on multiple lanes.
Counts are aggregated to their mean.
"""

import matplotlib.pyplot as plt

from utils import _handle_io

import matplotlib as mpl
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'


# def aggregate_lanes(input_file_path, output_file_path, *args, **kwargs):
#     data = pd.read_csv(input_file_path, sep='\t', index_col=['plate', 'well', 'lane'])
#     merged = _aggregate_lanes(data)
#     merged.to_csv(output_file_path, sep='\t')


@_handle_io
def aggregate_lanes(data, show=False):
    if show:
        correlate_read_counts_from_different_lanes(data)
    return data.groupby(['plate', 'well']).aggregate('mean')


def correlate_read_counts_from_different_lanes(data, figure_path=None):
    samples = data[['plate', 'well']].drop_duplicates()
    lanes = set(data['lane'].values)

    maximum = 10000 # maximum read count
    fig, axes = plt.subplots(len(lanes), len(lanes), sharex=True, sharey=True)
    for ii, lane1 in enumerate(lanes):
        for jj, lane2 in enumerate(lanes):
            ax = axes[ii,jj]
            for plate, well in samples:
                x = data.loc[(plate, well, lane1)]
                y = data.loc[(plate, well, lane2)]
                ax.plot(x, y, 'o', ms=0.1)
                ax.plot([0, maximum], [0, maximum], 'k--')
            ax.set_xlim(0, maximum)
            ax.set_ylim(0, maximum)
            if jj == 0:
                ax.set_ylabel(lane1)
            if ii == len(lanes) -1:
                ax.set_xlabel(lane2)

    fig.tight_layout()

    if figure_path:
        fig.savefig(figure_path)


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Aggregate samples that have been sequenced on multiple lanes.")
    parser.add_argument(
        'input',
        help="/path/to/gene_abundances.tsv"
    )
    parser.add_argument(
        '-o', '--output',
        help="/path/to/aggregated_gene_abundances.tsv",
        default="./aggregated_gene_abundances.tsv")
    args = parser.parse_args()

    aggregate_lanes(args.input, args.output)
