#!/usr/bin/Rscript
library(SingleCellExperiment)
library(scater)
library(stringr)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
#print(args)
#rds <- args[1] %>%
#  str_replace_all("\\[|\\]|'", "")

rds <- unlist(strsplit(args[1], split = ","))

createSCE <- function(files){
  listSCE <- lapply(files, readRDS)

  for(i in 1:length(listSCE)){
    rowData(listSCE[[i]]) <- NULL
  }
  sce <- do.call('cbind', listSCE)
  
  sce
}

sce <- createSCE(rds)

saveRDS(sce, file = args[2])
