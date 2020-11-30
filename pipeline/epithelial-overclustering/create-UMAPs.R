#!/usr/local/bin/Rscript
library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(devtools)
devtools::load_all("~/taproom/")

args <- commandArgs(trailingOnly = TRUE)
cells <- read_csv(args[1])
expression <- readRDS(args[2])[,cells$cell_id]
out_dir <- args[3]
markers <- args[4]
markers_list <- read_markers(args[5])
perc <- args[6]

if(markers == "specified_markers"){
  expression <- expression[unique(unlist(markers_list$cell_types)),]
}

expression <- runUMAP(expression)

umap <- reducedDim(expression) %>% as.data.frame() %>% rownames_to_column("cell_id")
colnames(umap) <- c("cell_id", "UMAP 1", "UMAP 2")
umap$percent_epithelial <- perc

write_csv(umap, paste0(out_dir, "UMAPs/UMAP-luminal-", perc, "-", markers, ".csv"))
