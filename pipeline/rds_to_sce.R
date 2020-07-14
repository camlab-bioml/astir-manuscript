#!/usr/bin/Rscript
library(SingleCellExperiment)
library(scater)
library(stringr)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
rds <- args[1] %>%
  str_replace_all("\\[|\\]|'", "")

rds <- unlist(strsplit(rds, split = ","))

createSCE <- function(files){
  listSCE <- lapply(files, readRDS)
  sce <- do.call('cbind', listSCE)
  
  sce
}

sce <- createSCE(rds)

sce
