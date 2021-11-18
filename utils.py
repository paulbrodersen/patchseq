import pathlib
import pandas as pd
import h5py as h5
import functools
import numpy as np
import time

from scipy.sparse import csc_matrix


def timeit(func):
    def wrapper(*args, **kwargs):
        tic = time.time()
        result = func(*args, **kwargs)
        toc = time.time()
        print(f"{func.__name__} completed in {toc-tic:.1f} seconds.")
        return result
    return wrapper


def _handle_io(func):
    # Wrap function such that it has two modus operandi:
    # 1) If we pass it data directly, proceed normally.
    # 2) If we pass it two file paths instead, load the data, then evaluate,
    #    then save out the data again.
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if isinstance(args[0], (str, pathlib.Path)):
            input_filepath, output_filepath = args[:2]
            data = read_data(input_filepath)
            output = func(data, *args[2:], **kwargs)
            write_data(output_filepath, output)
            return # nothing
        else:
            return func(*args, **kwargs)
    return wrapper


def read_data(filepath):
    pathlib_path = pathlib.Path(filepath)
    assert pathlib_path.exists(), "Input file path does not exist!"

    suffix = pathlib_path.suffix
    if suffix == '.tsv':
        data = read_tsv(filepath)
    elif suffix == '.tome':
        data = read_tome(filepath)
    else:
        raise NotImplementedError("Read functions not implemented for file with suffix {suffix} yet!")
    return data


def read_tsv(filepath):
    return pd.read_csv(filepath, sep='\t', index_col=['plate', 'well', 'lane'])


@timeit
def read_tome(filepath, mode='r'):
    with h5.File(filepath, mode=mode) as f:
        gene_names   = f['gene_names'][:]
        sample_names = f['sample_names'][:]
        exons        = _read_csc(f['data']['exon'])
        introns      = _read_csc(f['data']['intron'])
    both = np.array((exons + introns).todense())
    return pd.DataFrame(both, index=sample_names, columns=gene_names)


def _read_csc(parent):
    x = parent['x']
    i = parent['i']
    p = parent['p']

    with x.astype(np.int):
        with i.astype(np.int):
            return csc_matrix((x, i, p))


@timeit
def read_loom(filepath):
    import loompy
    loom = loompy.connect(filepath)
    total_rows, total_columns = loom.shape
    data = loom[:total_rows, :total_columns]
    sample_names = loom.col_attrs['CellID']
    gene_names = loom.row_attrs['Gene']
    return pd.DataFrame(data.T, index=sample_names, columns=gene_names)


def write_data(filepath, data):
    pathlib_path = pathlib.Path(filepath)
    if pathlib_path.exists():
        answer = input('Output file exists. Overwrite? y/N\n')
        if answer.lower() != 'y':
            return

    suffix = pathlib_path.suffix
    if suffix == '.tsv':
        write_tsv(filepath, data)
    elif suffix == '.tome':
        write_tome(filepath, data)
    else:
        raise NotImplementedError(f"Write functions not implemented for file with suffix {suffix} yet!")


def write_tsv(filepath, df):
    df.to_csv(filepath, sep='\t')


def write_tome(filepath, data):
    raise NotImplementedError
