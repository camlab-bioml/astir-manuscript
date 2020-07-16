#!/usr/local/bin/Rscript

library(SingleCellExperiment)
library(tidyverse)
library(devtools)
library(scater)
devtools::load_all('~/taproom/')
source("~/imc-2020/scripts/functions.R")

args <- commandArgs(trailingOnly = TRUE)

### [GET ARGUMENTS] #####


### [READ IN DATA] #####
sce <- assignIdentity(params$cells, params$celltypes, params$cellstates)
pathways <- sce$pathways
sce <- sce$sce

clusters <- read_csv(params$clusters) %>% column_to_rownames(var = "id")
sce$clusters <- clusters[sce$id, 1]

# Just for the time being until we switch everything over to remove s' at the end
sce$cell_type[which(sce$cell_type == "B cell")] <- "B cells"
sce$cell_type[which(sce$cell_type == "T cell")] <- "T cells"


### [DO THE PLOTTING] #####
phenograph.bar <- sce %>% 
  colData() %>% 
  as.data.frame() %>% 
  select(pathways, "clusters") %>% 
  pivot_longer(cols = pathways, 
               names_to = "pathway", values_to = "activation") %>% 
  mutate(clusters = as.factor(clusters))


ggplot(phenograph.bar, aes(x = clusters, y = activation)) +
  geom_boxplot() +
  facet_wrap("pathway") +
  astir_paper_theme()
