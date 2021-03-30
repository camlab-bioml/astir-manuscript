#!/usr/local/bin/Rscript

library(SingleCellExperiment)
library(tidyverse)
library(devtools)
library(scater)
library(ggalluvial)
devtools::load_all('~/taproom/')
source("~/imc-2020/scripts/functions.R")

args <- commandArgs(trailingOnly = TRUE)

### [GET PARAMETERS] #####
cells <- args[1]
celltypes <- args[2]
#cellstates <- args[3]
clusters <- args[3]
cohort <- args[4]
method <- args[5]
clustering_params <- args[6]
output_dir <- args[7]

### [READ IN DATA] #####
thresh <- 0.5

# read in data 
sce <- readRDS(cells)
if(cohort == "lin_cycif" | cohort == "wagner" | cohort == "schapiro" | cohort == "keren"){
  # add a cell ID
  sce$id <- colnames(sce)
}

type <- read_csv(celltypes)

type$cell_type <- get_celltypes(select(type, -X1), thresh) 
type <- select(type, X1, cell_type) %>%
  column_to_rownames("X1")
colData(sce)["cell_type"] <- type[colnames(sce),]

#sce <- assignIdentity(cells, celltypes, cellstates)$sce
clusters <- read_csv(clusters) %>% column_to_rownames(var = "id")


# Just for the time being until we switch everything over to remove s' at the end
sce$cell_type[which(sce$cell_type == "B cell")] <- "B cells"
sce$cell_type[which(sce$cell_type == "T cell")] <- "T cells"

sce$clusters <- clusters[sce$id, 1]
alluv <- sce %>% colData() %>% as.data.frame() %>% 
  select(cell_type, "clusters") %>% group_by(cell_type, clusters) %>% 
  count()

alluv <- alluv %>% as.data.frame()

### [PLOT RESULTS] ##### 
filename <- paste("Alluv_cohort", cohort, "method", method, "markers", clustering_params, sep = "_")

pdf(file = paste0(output_dir, filename, ".pdf"), height = 3)
create.alluvial(alluv, method)
dev.off()
