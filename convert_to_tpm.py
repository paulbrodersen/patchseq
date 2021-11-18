#!/usr/bin/env python
"""
Convert gene abundance in absolute counts to tpm.
"""

import numpy as np

from utils import _handle_io

# import pandas as pd

# def convert_to_tpm(input_file_path, output_file_path):
#     gene_abundance = pd.read_csv(input_file_path, sep='\t', index_col=['plate', 'well', 'lane'])
#     converted = _convert_to_tpm(gene_abundance)
#     converted.to_csv(output_file_path, sep='\t')

@_handle_io
def convert_to_tpm(absolute_counts):
    return absolute_counts / absolute_counts.sum(axis=1)[:, np.newaxis] * 1e+6


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Convert gene counts to tpm.")
    parser.add_argument(
        'input',
        help="/path/to/gene_abundance_in_absolute_counts.tsv"
    )
    parser.add_argument(
        '-o', '--output',
        help="/path/to/gene_abundance_in_tpm.tsv",
        default="./gene_abundance_in_tpm.tsv")
    args = parser.parse_args()

    convert_to_tpm(args.input, args.output)
