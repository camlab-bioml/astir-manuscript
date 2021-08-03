#!/usr/bin/Rscript
library(SingleCellExperiment)
library(scater)
library(stringr)
library(dplyr)
library(argparse)

parser <- ArgumentParser(description = "test")

parser$add_argument('--rds_files', type = 'character', nargs = '+')
parser$add_argument('--output_file', type = 'character')

args <- parser$parse_args()

rds <- strsplit(args$rds_files, split = " ")[[1]]

createSCE <- function(files){
  listSCE <- lapply(files, readRDS)

  for(i in 1:length(listSCE)){
    rowData(listSCE[[i]]) <- NULL
  }
  sce <- do.call('cbind', listSCE)
  
  sce
}

sce <- createSCE(rds)

saveRDS(sce, file = args$output_file)