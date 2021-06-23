#!/usr/local/bin/Rscript

library(devtools)
library(SingleCellExperiment)
library(tidyverse)
library(scater)
library(circlize)
library(ComplexHeatmap)
source("scripts/functions.R")
#library(taproom)
library(devtools)
devtools::load_all("../taproom/")

args <- commandArgs(trailingOnly = TRUE)
cells <- args[1]
type <- read_csv(args[2])
#state <- args[3]
markers <- read_markers(args[3])
cohort <- args[4]
output_dir <- args[5]

thresh <- 0.5

sce <- readRDS(cells)
type$cell_type <- get_celltypes(select(type, -X1), thresh) 
type <- select(type, X1, cell_type) %>%
  column_to_rownames("X1")

colData(sce)["cell_type"] <- type[colnames(sce),]

if(cohort == "zurich1" | cohort == "wagner"){
    w1 = 14
}else{
    w1 = 18
}

pdf(paste0(output_dir, "Astir_expressionHeatmap_allMarkers_", cohort, ".pdf"), height = 11, width = w1)
createHeatmap(sce)
dev.off()

if(cohort == "zurich1"){
    w = 15
}else{
    w = 16
}

# if(cohort == "lin_cycif"){
#     cell_type_probs <- read_csv(args[2]) %>% column_to_rownames("X1")
#     types <- colnames(cell_type_probs) 

#     colData(sce)[types] <- cell_type_probs[colnames(sce),]

#     col_fun = colorRamp2(c(0, 0.5, 1), c("#440154FF", "#218F8DFF", "#FDE725FF"))

#     createHeatmap <- function(sce,
#                             cell_type_column = "cell_type",
#                             assay = "logcounts",
#                             thresh = 2) {
#     assay <- "logcounts"
#     thresh <- 0.7
#     lc <- t(as.matrix(assay(sce, assay)))
#     lc <- scale(lc)
    
#     lc[lc > thresh] <- thresh
#     lc[lc < -thresh] <- -thresh
    
#     cell_types = colData(sce)[[ cell_type_column ]]

#     stromal = log(1 - colData(sce)[["Stromal"]] + 0.001)
#     other = log(1 - colData(sce)[["Other"]] + 0.001)
#     epithelial = log(1 - colData(sce)[["Epithelial"]] + 0.001)

#     min.v <- min(c(stromal, other, epithelial))
#     max.v <- max(c(stromal, other, epithelial))

#     print(min.v)
#     print(max.v)

#     col_fun = colorRamp2(c(min.v, max.v), c("#440154FF", "#FDE725FF"))
    
#     celltype_annot <- HeatmapAnnotation(`Cell type` = cell_types, 
#                                         stromal_prob = stromal,
#                                         other_prob = other,
#                                         epithelial_prob = epithelial,
#                                         which="column",
#                                         col = list(`Cell type` = jackson_basel_colours(),
#                                                     stromal_prob = col_fun,
#                                                     other_prob = col_fun,
#                                                     epithelial_prob = col_fun))  
    
#     type_exprs <- Heatmap(t(lc), 
#                             name = "Expression",
#                             column_title = "Cell",
#                             col=viridis(100),
#                             top_annotation = celltype_annot,
#                             show_column_names = FALSE,
#                             column_order = order(cell_types))
#     type_exprs
# }
# }

pdf(paste0(output_dir, "Astir_expressionHeatmap_specifiedMarkers_", cohort, ".pdf"), width = w, height = 6)
    createHeatmap(sce[unique(unlist(markers$cell_types)),])
dev.off()