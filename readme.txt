Basic workflow

1. 01_setup_kallisto.py        : setup folder structure expected by kallisto.py
2. 02_kallisto.py              : estimate gene counts using kallisto
3. 03_preprocess.py            : preprocess reference data set; adjust parameters in quality_control.py as necessary; preprocess own samples using the same parameters
4. 04_highly_variable_genes.py : determine highly variable genes in reference data sets
5. 05_embedding.py             : find a low-dimensional embedding of the HVG gene counts of the reference; project own samples onto embedding
