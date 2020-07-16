library(SingleCellExperiment)
library(tidyverse)
library(devtools)
library(scater)
library(ggalluvial)
devtools::load_all('~/taproom/')
source("~/imc-2020/scripts/functions.R")

args <- commandArgs(trailingOnly = TRUE)

cells <- args[1]
celltypes <- args[2]
cellstates <- args[3]
clusters <- args[4]
markers <- args[5]
clusters <- args[6]
output_dir <- args[7]

### [READ IN DATA] #####
sce <- assignIdentity(cells, celltypes, cellstates)$sce
clusters <- read_csv(clusters) %>% column_to_rownames(var = "id")
markers <- read_markers(params$markers)
clusters <- read_csv(clusters) %>% column_to_rownames(var = "id")

sce$clusters <- clusters[sce$id, 1]

# Just for the time being until we switch everything over to remove s' at the end
sce$cell_type[which(sce$cell_type == "B cell")] <- "B cells"
sce$cell_type[which(sce$cell_type == "T cell")] <- "T cells"

imc.data <- t(logcounts(sce[unique(unlist(markers$cell_types)),]))

### [PLOTTING] #####
lc <- scale(imc.data)
thresh <- 2
lc[lc > thresh] <- thresh
lc[lc < -thresh] <- -thresh

ha <- HeatmapAnnotation(`Cell type` = sce$cell_type,
                        Cluster = as.factor(sce$clusters),
                        which="column",
                        col = list(`Cell type` = jackson_basel_colours()),
                        `Cluster #` = sce$clusters)

exprs <- Heatmap(t(lc), name = "Expression",
                 column_title = "Cell",
                 col=viridis(100),
                 top_annotation = ha,
                 show_column_names = FALSE,
                 column_order = order(sce$clusters))

filename <- paste("ExpressionHeatmap_cohort", cohort, "method", method, 
                  "markers", clustering.params, sep = "_")

pdf(file = paste0(output_dir, filename, ".pdf"))
draw(exprs)
dev.off()