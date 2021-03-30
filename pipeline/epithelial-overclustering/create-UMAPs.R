#!/usr/local/bin/Rscript
library(tidyverse)
library(SingleCellExperiment)
library(scater)

args <- commandArgs(trailingOnly = TRUE)
cells <- read_csv(args[1])
expression <- readRDS(args[2])[,cells$cell_id]
out_dir <- args[3]
markers <- args[4]

cells <- read_csv("output/chipmunk/results/epithelial_overclustering/luminal_80.csv")
out_dir <- "output/chipmunk/results/epithelial_overclustering/"
expression <- readRDS("output/chipmunk/sces/wagner_sce.rds")[,cells$cell_id]
markers <- "all_markers"

expression <- runUMAP(expression)

umap <- reducedDim(expression) %>% as.data.frame() %>% rownames_to_column("cell_id")
colnames(umap) <- c("cell_id", "UMAP 1", "UMAP 2")

write_csv(umap, paste0(out_dir, "UMAPs/UMAP-luminal-0-", markers, ".csv"))
