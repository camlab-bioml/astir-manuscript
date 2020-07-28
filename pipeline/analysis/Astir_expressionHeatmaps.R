#!/usr/local/bin/Rscript

library(devtools)
library(SingleCellExperiment)
library(tidyverse)
library(scater)
source("scripts/functions.R")
devtools::load_all("~/taproom/")

args <- commandArgs(trailingOnly = TRUE)
cells <- args[1]
type <- args[2]
state <- args[3]
markers <- read_markers(args[4])
cohort <- args[5]
output_dir <- args[6]

sce <- assignIdentity(cells, type, state)$sce

pdf(paste0(output_dir, "Astir_expressionHeatmap_allMarkers_", cohort, ".pdf"), height = 14, width = 18)
createHeatmap(sce)
dev.off()

pdf(paste0(output_dir, "Astir_expressionHeatmap_specifiedMarkers_", cohort, ".pdf"), height = 10, width = 18)
createHeatmap(sce[unique(unlist(markers$cell_types)),])
dev.off()