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

cells <- readRDS(args[1])
type <- args[2] %>% read_csv %>% column_to_rownames("X1")
markers <- read_markers(args[3])
cohort <- args[4]
output_dir <- args[5]

# assign cell types
type$cell_type <- taproom::get_celltypes(type)
type$cell_type[type$cell_type == "B cell"] <- "B cells"
type$cell_type[type$cell_type == "T cell"] <- "T cells"

type.mat <- type %>% 
  select(-"cell_type") %>% 
  as.matrix()

ha <- HeatmapAnnotation(`Cell Type` = type$cell_type,
                        which = "column",
                        col = list(`Cell Type` = jackson_basel_colours()))

# Do the plotting
pdf(paste0(output_dir, "celltype_probability_", cohort, ".pdf"), width = 16)
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
if(cohort == "wagner"){
 cells$cell_type[cells$cell_type == "Fibroblasts"] <- "Stromal" 
}


unknown <- cells[,cells$cell_type == "Unknown" | cells$cell_type == "Stromal"]
thresh <- 2

lc <- t(as.matrix(logcounts(unknown[unique(unlist(markers$cell_types))])))
lc <- scale(lc)

lc[lc > thresh] <- thresh
lc[lc < -thresh] <- -thresh

ha <- HeatmapAnnotation(`Cell Type` = unknown$cell_type,
                        which = "column",
                        col = list(`Cell Type` = jackson_basel_colours()))

pdf(paste0(output_dir, "Unknown_celltype_expression_", cohort, ".pdf"), height = 5, width = 10)
Heatmap(t(lc), 
        name = "Expression",
        column_title = "Cell",
        col=viridis(100),
        top_annotation = ha,
        cluster_columns = F,
        show_column_names = FALSE,
        column_order = order(unknown$cell_type))
dev.off()
