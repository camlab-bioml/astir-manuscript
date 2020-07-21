#!/usr/local/bin/Rscript

library(SingleCellExperiment)
library(tidyverse)
library(scater)
library(ComplexHeatmap)
library(devtools)
library(viridis)
devtools::load_all("~/taproom")
source("scripts/functions.R")

### [READ IN DATA] #####
args <- commandArgs(trailingOnly = TRUE)

type <- args[1] %>% column_to_rownames("X1")
cells <- readRDS(args[2])
markers <- read_markers(args[3])
cohort <- args[4]
output_dir <- args[5]

# assign cell types
type$cell_type <- taproom::get_celltypes(type)

type.mat <- type %>% 
  select(-"cell_type") %>% 
  as.matrix()

ha <- HeatmapAnnotation(`Cell Type` = type$cell_type,
                        which = "column",
                        col = list(`Cell Type` = jackson_basel_colours()))

# Do the plotting
pdf(paste0(output_dir, "celltype_probability_", cohort, ".pdf"))
Heatmap(t(type.mat), 
        col = viridis(100), 
        name = "Expression", 
        top_annotation = ha,
        show_column_names = F,
        cluster_columns = F,
        column_order = order(type$cell_type))
dev.off()

### [EXPRESSION OF MARKERS IN UNKNOWN CELLS]
# add cell types to expression data
cells$cell_type <- type[colnames(cells), ]$cell_type

# subset unknown cells
unknown <- cells[,cells$cell_type == "Unknown"]
thresh <- 2

lc <- t(as.matrix(logcounts(unknown[unique(unlist(markers$cell_types))])))
lc <- scale(lc)

lc[lc > thresh] <- thresh
lc[lc < -thresh] <- -thresh

pdf(paste0(output_dir, "Unknown_celltype_expression_", cohort, ".pdf"))
Heatmap(t(lc), 
        name = "Expression",
        column_title = "Cell",
        col=viridis(100),
        show_column_names = FALSE)
dev.off()
