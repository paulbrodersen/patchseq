#!/usr/bin/env Rscript
## suppress messages such that we can pipe the output into a file
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scran))
suppressMessages(library(scater))

# read command line arguments, i.e. the input and output file paths
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("At least one argument must be supplied (input file path).n", call.=FALSE)
} else if (length(args)==1) {
    args[2] = paste("normalized", args[1], sep='_')
}

input_file_path = args[1]
output_file_path = args[2]

# load data and split into sample info (plate, well), gene names, and counts
data <- read.table(file=input_file_path, sep='\t', header=TRUE)
sample_names <- data[,1]
counts <- data[,2:length(data)]
gene_names <- colnames(data[,2:length(data)])

# convert to SCE and normalize
sce <- SingleCellExperiment(assays=list(counts=as.matrix(t(counts))))
sce <- computeSumFactors(sce)
summary(sizeFactors(sce))

## ## pre-clustering results in weird plot, in which for some samples,
## ## the correlation between library size and size factor is perfect,
## ## and for others there is no apparent correlation;
## ## this is potentially due to the small number of cells.
## preclusters <- quickCluster(sce, min.size=20)
## sce2 <- computeSumFactors(sce, clusters=preclusters)
## summary(sizeFactors(sce2))

## plot(librarySizeFactors(sce), sizeFactors(sce), log="xy")
## plot(librarySizeFactors(sce), sizeFactors(sce))

sce <- logNormCounts(sce)

# re-assemble the table and write to std::out
normed_counts <- t(logcounts(sce))
colnames(normed_counts) = gene_names
normed_table <- cbind(sample_names, normed_counts)
write.table(normed_table, file=output_file_path, sep="\t", quote=FALSE)
