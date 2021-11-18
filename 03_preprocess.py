#!/usr/bin/env python
"""
Preprocess transcript abundance matrices.
1) Merge data from different plates.
2) Quality control.
2) Merge data from different lanes.
3) Normalize with scran.
"""

import os

from utils import _handle_io

from quality_control import quality_control
from convert_to_tpm import convert_to_tpm
from aggregate_lanes import aggregate_lanes


@_handle_io
def remove_control_samples(data):
    blank_wells = ['B7', 'D12', 'F4']
    buffer_wells = ['C11', 'E7', 'F2']
    well_blacklist = blank_wells + buffer_wells
    return data.drop(index=well_blacklist, level=1)


if __name__ == '__main__':

    # code, figure, and data folders
    cdir = '/path/to/code/'
    fdir = '/path/to/figures/'
    ddir = '/path/to/data/'

    file_base_name = 'mydata'

    input_file_path = ddir + file_base_name + '.tsv'

    # --------------------------------------------------------------------------------
    # remove control samples (blanks, etc)

    print('Removing control samples')
    without_controls_file_path = ddir + file_base_name +'--without_controls.tsv'
    remove_control_samples(input_file_path, without_controls_file_path)

    # --------------------------------------------------------------------------------
    # quality control

    print("Quality control...")
    qced_file_path = ddir + file_base_name +'--quality_control.tsv'
    quality_control(without_controls_file_path, qced_file_path, show=True)

    # --------------------------------------------------------------------------------
    # convert to tpm

    print("Convert to tpm...")
    tpm_file_path = ddir + file_base_name + '--tpm.tsv'
    convert_to_tpm(qced_file_path, tpm_file_path)

    # --------------------------------------------------------------------------------
    # merge lanes if any

    print("Merge lanes (if any)...")
    aggregated_file_path = ddir + file_base_name + '--aggregated.tsv'
    aggregate_lanes(tpm_file_path, aggregated_file_path)

    # --------------------------------------------------------------------------------
    # normalize counts using scran

    print("Normalize counts...")
    normalized_file_path = ddir + file_base_name + '--normalized.tsv'
    conda_env = 'rbase'
    script = 'scran_normalization.R'
    # script = 'scran_normalization_multi_index.R'
    cmd = f'conda run -n {conda_env} Rscript {cdir}{script} {aggregated_file_path} {normalized_file_path}'
    os.system(cmd)
