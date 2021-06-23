library(zellkonverter)
library(SingleCellExperiment)
library(scRNAseq)


sce <- readRDS(snakemake@input[['sce']])

# Remove raw_imc and unwinsorized counts so that these are not copied to h5ad
assays(sce)$raw_imc <- NULL
assays(sce)$logcounts_unwinsorized <- NULL

writeH5AD(sce, snakemake@output[['h5ad']])