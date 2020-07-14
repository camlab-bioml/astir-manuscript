#!/usr/bin/Rscript
args = commandArgs(trailingOnly = TRUE)
rds <- args[1]

library(SingleCellExperiment)
library(scater)

createSCE <- function(files){
  listSCE <- lapply(files, readRDS)
  sce <- do.call('cbind', listSCE)
  
  sce
}

sce <- createSCE(rds)

sce
