#!/usr/bin/env python
"""
Automate the setup of kallisto analysis.
Creates the necessary directory structure,
and links in the appropriate files into the data subdirectory.
"""

import re
import shutil
import pathlib
import textwrap
import pandas as pd


def main(plate_name, path_to_data, path_to_email, batch_name):
    project_paths = setup_folders(
        project='/ifs/projects/proj093/',
        pipeline='kallisto',
        plate=plate_name,
    )

    # gather and link data sets
    link_in_data(src_data_dir=pathlib.Path(path_to_data),
                 dst_data_dir=project_paths['data'],
                 path_to_email=path_to_email,
                 restrict_to=plate_name.split('-')[-1],
                 batch_name=batch_name,
    )

    # create pipeline.yml file
    create_yml(project_paths['plate'])



def setup_folders(project, pipeline, plate):
    """Creates the following folder structure for RNASeq data analysis:

    project
        /analysis
            /pipeline
                /plate
                    /data

    and returns a dictionary that maps each element to a pathlib Path object.

    Arguments:
    ----------
    project : str
         path/to/project
    pipeline : str
         pipeline name
    plate : str
        plate name

    Returns:
    --------
    paths : dict str : pathlib.Path
        Mapping of folder names to paths.
    """
    # construct all paths
    project_dir = pathlib.Path(project)
    pipeline_dir = project_dir.joinpath(f'analysis/{pipeline}')
    plate_dir = pipeline_dir.joinpath(plate)
    data_dir = plate_dir.joinpath('data')

    # make all directories
    data_dir.mkdir(parents=True, exist_ok=True)

    return dict(project=project_dir,
                pipeline=pipeline_dir,
                plate=plate_dir,
                data=data_dir
    )


def link_in_data(src_data_dir, dst_data_dir, path_to_email,
                 restrict_to=None, batch_name='lane1'):
    """

    Arguments:
    ----------
    src_data_dir
        Backed up data folder.
    dst_data_dir
        Folder where to put the symlinks.
    path_to_email
        Path to the email from the sequencing facility.
    restrict_to : str (optional, default None)
        A plate name such as 'Plate15' or None.
        A sequencing batch may contain samples from more than one plate.
        If this argument is specified (i.e. not None), only samples from that
        plate are processed.
    batch_name: str (default 'lane1')
        The sequencing batch. This information is not contained in the email,
        but should be consistent for all files within one email as each email
        maps to a single batch.

    """

    dataframe = parse_email(path_to_email)

    if restrict_to:
        dataframe = dataframe[dataframe['plate']==restrict_to]

    src_name_to_dst_name = remap_file_names(dataframe, batch_name)

    create_symlinks(src_data_dir, dst_data_dir, src_name_to_dst_name)


def remap_file_names(dataframe, batch_name='lane1'):
    """Map sample IDs to human readable file names.

    The file names assigned by the sequencing facility are not human
    readable. Based on the data provided in the email from the facility,
    we create a mapping from the sample ID to a human readable name with
    the following format:

    <sample ID>_<strand>.fastq.gz : <plate>-<well>-<lane>.fastq.<strand>.gz

    Arguments:
    ----------
    dataframe: pandas.DataFrame
        The data from the parsed email from the sequencing facility.
        Contains the following columns:
        - sample ID
        - plate
        - well
    batch_name: str (default 'lane1')
        The sequencing batch. This information is not contained in the email,
        but should be consistent for all files within one email as each email
        maps to a single batch.

    Returns:
    --------
    mapping : dict
        Mapping from sample IDs to human readable file names.

    TODO:
    -----
    - Handle non-stranded sequencing data.
    """

    mapping = dict()
    for ii, row in dataframe.iterrows():
        for strand in [1, 2]:
            src = f"{row['sample_id']}_{strand}.fastq.gz"
            dst = f"{row['plate']}-{row['well']}-{batch_name}.fastq.{strand}.gz"
            mapping[src] = dst

    return mapping


def parse_email(file_path):
    """
    Given an email with the following format:

    WTCHG_693615_70015361 Sample:Plate15_A10 Count:12577357
    WTCHG_693615_70025362 Sample:Plate15_B10 Count:4835494
    WTCHG_693615_70105370 Sample:Plate15_B11 Count:7296271
    WTCHG_693615_70115371 Sample:Plate15_C11_lysis_ctl Count:732017
    WTCHG_693615_70125372 Sample:Plate15_D11 Count:1453627
    WTCHG_693615_70195379 Sample:Plate15_C12 Count:1842041
    WTCHG_693615_70205380 Sample:Plate15_D12_ceph10pgRNA_ctl Count:756343
    WTCHG_693615_70215381 Sample:Plate15_E12 Count:65603
    WTCHG_693615_71345014 Sample:single_GM12878_ctl_14 Count:776734

    return a pandas dataframe with the following format:

    sample_id,             sample_name,                 plate,   well, count
    WTCHG_693615_70015361, Plate15_A10,                 Plate15, A10,  12577357
    WTCHG_693615_70025362, Plate15_B10,                 Plate15, B10,  4835494
    WTCHG_693615_70105370, Plate15_B11,                 Plate15, B11,  7296271
    WTCHG_693615_70115371, Plate15_C11_lysis_ctl,       Plate15, C11,  732017
    WTCHG_693615_70125372, Plate15_D11,                 Plate15, D11,  1453627
    WTCHG_693615_70195379, Plate15_C12,                 Plate15, C12,  1842041
    WTCHG_693615_70205380, Plate15_D12_ceph10pgRNA_ctl, Plate15, D12,  756343
    WTCHG_693615_70215381, Plate15_E12,                 Plate15, E12,  65603
    WTCHG_693615_71345014, single_GM12878_ctl_14,              ,    ,  776734
    """

    with open(file_path) as f:
        lines = f.readlines()

    # remove all lines that do not appear to contain a sample ID;
    # this should remove
    # a) empty lines
    # b) the header
    lines = [line for line in lines if 'WTCHG' in line]

    # parse remaining lines using regex
    data = []
    problems = []
    for line in lines:
        sample_id    = re.search('^WTCHG_[0-9]+_[0-9]+', line).group(0)
        sample_name  = re.search('Sample:([A-z0-9_]+)', line).group(1)
        count_as_str = re.search('Count:([0-9]+)', line).group(1)

        # Not all sample names follow the <plate>_<well><rest> format,
        # as, for example, shown in the last line of the sample text.
        # re.search may thus return None, which we need to handle.
        match = re.search('(Plate[0-9]{1,2})_([A-Z][0-9]{1,2})', line)
        if match:
            plate, well = match.groups()
        else:
            plate, well = '',''
            problems.append(line)

        data.append((sample_id, sample_name, plate, well, int(count_as_str)))

    if problems:
        msg = "Could not determine plate and/or well assignement in the following lines:\n"
        msg += ''.join(problems)
        import warnings
        warnings.warn(msg)

    return pd.DataFrame(data,
                        columns=['sample_id',
                                 'sample_name',
                                 'plate',
                                 'well',
                                 'count'])


def create_symlinks(src_data_dir, dst_data_dir, src_name_to_dst_name):
    """Create symlinks in the project working directory to the raw sequencing data in the backup folder."""
    assert src_data_dir.exists(), "Source directory does not exist."
    assert dst_data_dir.exists(), "Target directory does not exist."
    successes = 0
    failures = 0
    for src_name, dst_name in src_name_to_dst_name.items():
        src = src_data_dir.joinpath(src_name)
        try:
            assert src.exists(), f"File {src} does not exist."
            dst = dst_data_dir.joinpath(dst_name)
            # shutil.copy(src, dst, follow_symlinks=False)
            dst.symlink_to(src)
            successes += 1
        except AssertionError as e:
            failures += 1
            print(e)
    print(f'{successes} / {successes + failures} files linked successfully.')
    print(f'{failures} / {successes + failures} files linked unsuccessfully.')


# # mouse, all genes
# def create_yml(plate_dir):
#     pipeline_yml_text = f"""
#     input_fastq_glob: {plate_dir}/data/*.gz
#     queue: all.q
#     threads: 12
#     index:
#         fasta: /ifs/projects/proj093/backup/mus_musculus/Mus_musculus.GRCm38.cdna.all.fa.gz
#         gtf: /ifs/projects/proj093/backup/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
#     quant:
#         index: /ifs/projects/proj093/backup/mus_musculus/Mus_musculus.GRCm38.index
#         options: --fr-stranded --genomebam --gtf /ifs/projects/proj093/backup/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
#     """
#     pipeline_yml_text = textwrap.dedent(pipeline_yml_text)
#     pipeline_yml = plate_dir.joinpath('pipeline.yml')
#     pipeline_yml.write_text(pipeline_yml_text)

# mouse, protein coding only
def create_yml(plate_dir):
    pipeline_yml_text = f"""
    input_fastq_glob: {plate_dir}/data/*.gz
    queue: all.q
    threads: 12
    index:
        fasta: /ifs/projects/proj093/backup/mus_musculus/gencode.vM25.pc_transcripts.fa.gz
        gtf: /ifs/projects/proj093/backup/mus_musculus/gencode.vM25.annotation.gtf.gz
    quant:
        index: /ifs/projects/proj093/backup/mus_musculus/gencode.vM25.index
        options: --fr-stranded --genomebam --gtf /ifs/projects/proj093/backup/mus_musculus/gencode.vM25.annotation.gtf.gz
    """
    pipeline_yml_text = textwrap.dedent(pipeline_yml_text)
    pipeline_yml = plate_dir.joinpath('pipeline.yml')
    pipeline_yml.write_text(pipeline_yml_text)


if __name__ == '__main__':

    plate_name = 'mouse-Plate1'
    path_to_email = '/path/to/email/from/sequencing/facility/sample_info.txt'
    path_to_data = '/path/to/folder/with/fastq.gz/files/'
    batch_name = 'lane1'
    main(plate_name, path_to_data, path_to_email, batch_name)
