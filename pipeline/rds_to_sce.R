#!/usr/bin/Rscript
library(SingleCellExperiment)
library(scater)
library(stringr)
library(dplyr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument('--rds', type = 'character', nargs = '+')
parser$add_argument('--output', type = 'character')
args <- parser$parse_args()

#args = commandArgs(trailingOnly = TRUE)
#rds <- unlist(strsplit(args[1], split = " "))
#print(args[1])

createSCE <- function(files){
  listSCE <- lapply(files, readRDS)

  for(i in 1:length(listSCE)){
    rowData(listSCE[[i]]) <- NULL
  }
  sce <- do.call('cbind', listSCE)
  
  sce
}

sce <- createSCE(args$rds)

saveRDS(sce, file = args$output)
