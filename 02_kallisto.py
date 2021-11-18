#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Simple pipeline using k-mer matching to compute approximate
transcript expression counts (probabilities).
"""

import sys
import pandas as pd
import pathlib

from ruffus import *
from cgatcore import pipeline as P


P.get_parameters('pipeline.yml')


@follows(mkdir('fastqc'))
@transform(P.PARAMS["input_fastq_glob"], regex("^.*/(.+?)\.fastq.([12]).gz$"), r'fastqc/\1.fastq.\2_fastqc.html')
def fastqc(infile, outfile):
    """Perform quality control on the raw reads."""
    statement = 'fastqc -q -t %(threads)s --nogroup %(infile)s --outdir fastqc'
    P.run(statement, job_queue=P.PARAMS['queue'], job_threads=P.PARAMS['threads'])

# TODO implement QC based on read quality

@merge([P.PARAMS['index_fasta'], P.PARAMS['index_gtf']], 'transcripts.fasta')
def create_fasta(infiles, outfile):
    """Create a FASTA file from a (mouse) genome file and a GTF file denoting the individual gene boundaries."""
    fasta, gtf = infiles
    statement = 'zcat %(gtf)s | gffread -w %(outfile)s -g %(fasta)s - '
    P.run(statement, job_queue=P.PARAMS['queue'], job_memory ='8G')

@transform(create_fasta, suffix(r'.fasta'), r'.index')
def create_index(infile, outfile):
    """Convert the annotated FASTA file to an index that kallisto can use."""
    statement = 'kallisto index -i %(outfile)s -g %(infile)s'
    P.run(statement, job_queue=P.PARAMS['queue'], job_memory ='8G')

# # all genes
# @transform(P.PARAMS['index_fasta'], suffix(r'.cdna.all.fa.gz'), r'.index')
# def create_index(infile, outfile):
#     """Convert the annotated FASTA file to an index that kallisto can use."""
#     statement = 'kallisto index -i %(outfile)s -g %(infile)s'
#     P.run(statement, job_queue=P.PARAMS['queue'], job_memory ='8G')

# protein coding genes
@transform(P.PARAMS['index_fasta'], suffix(r'.pc_transcripts.fa.gz'), r'.index')
def create_index(infile, outfile):
    """Convert the annotated FASTA file to an index that kallisto can use."""
    statement = 'kallisto index -i %(outfile)s -g %(infile)s'
    P.run(statement, job_queue=P.PARAMS['queue'], job_memory ='8G')


@follows(mkdir('quant'), create_index)
@collate(P.PARAMS["input_fastq_glob"], regex("^.*/(.+?)\.fastq.([12]).gz$"), r'quant/\1/abundance.tsv')
def get_transcript_counts(infiles, outfile):
    read1, read2 = infiles
    outdir = outfile.replace('abundance.tsv', '')
    sample_name = outdir.split('/')[1]
    stdout_file = outdir + sample_name + '.log'
    # statement = 'kallisto quant -t %(threads)s %(quant_options)s -i transcripts.index -o %(outdir)s %(read1)s %(read2)s &> %(stdout_file)s'
    statement = 'kallisto quant -t %(threads)s %(quant_options)s -i %(quant_index)s -o %(outdir)s %(read1)s %(read2)s &> %(stdout_file)s'
    P.run(statement, job_queue=P.PARAMS['queue'], job_threads=P.PARAMS['threads'], job_memory ='8G')


@merge(get_transcript_counts, 'multiqc_report.html')
def multiqc(infiles, outfile):
    statement = '''multiqc .'''
    P.run(statement)


@merge(get_transcript_counts, 'transcript_abundance.tsv')
def join_counts_tables(infiles, outfile):
    data_frames = [pd.read_csv(infile, sep='\t', index_col='target_id') for infile in infiles]
    sample_names = [pathlib.Path(infile).parents[0].stem for infile in infiles]
    # counts_series = [df['tpm'] for df in data_frames]
    counts_series = [df['est_counts'] for df in data_frames]
    # rename counts variable to sample label so we can keep track of the sample origin in the combined tablwe
    counts_series = [series.rename(sample_name) for sample_name, series in zip(sample_names, counts_series)]
    combined_table = pd.concat(counts_series, axis='columns')
    combined_table.to_csv(outfile, sep='\t')


# @transform(join_counts_tables, suffix(r'transcript_abundance.tsv'), r'gene_abundance.tsv')
# def get_gene_counts(infile, outfile):
#     transcript_abundance = pd.read_csv(infile, sep='\t')
#     transcript_to_gene = pd.read_csv('/ifs/projects/proj093/backup/mus_musculus/transcripts_to_genes.txt',
#                                      sep='\t',
#                                      header=None,
#                                      names=['transcript_id', 'gene_id', 'gene_name'])
#     merged = pd.merge(transcript_to_gene, transcript_abundance, left_on='transcript_id', right_on='target_id', how='inner')
#     funcs = {column_name : 'sum' if 'Plate' in column_name else 'first' for column_name in merged.columns.values}
#     gene_abundance = merged.groupby('gene_id').aggregate(funcs)
#     gene_abundance.to_csv('gene_abundance.tsv', sep='\t')


@transform(join_counts_tables, suffix(r'transcript_abundance.tsv'), r'gene_abundance.tsv')
def get_gene_counts(infile, outfile):
    transcript_abundance = pd.read_csv(infile, sep='\t')
    # split target_id string and expand into columns
    transcript_info = transcript_abundance['target_id'].str.split('|').apply(pd.Series)[[0, 1, 5]].rename(columns={0 : 'gene_id', 1 : 'transcript_id', 5 : 'gene_name'})
    transcript_info['gene_name'] = transcript_info['gene_name'].apply(lambda x : x.split('.')[0]) # remove gene splice (?) variant information
    merged = pd.concat([transcript_info, transcript_abundance], axis=1)
    funcs = {column_name : 'sum' if 'Plate' in column_name else 'first' for column_name in merged.columns.values}
    gene_abundance = merged.groupby('gene_name').aggregate(funcs)
    gene_abundance.to_csv('gene_abundance.tsv', sep='\t')


# TODO deal with batch effects
# - combat
# - MNN

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
