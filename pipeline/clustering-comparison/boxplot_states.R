#!/usr/local/bin/Rscript 

library(SingleCellExperiment)
library(tidyverse)
library(devtools)
library(scater)
devtools::load_all('~/taproom/')
source("~/imc-2020/scripts/functions.R")

args <- commandArgs(trailingOnly = TRUE)

### [GET ARGUMENTS] #####
cells <- args[1]
types <- args[2]
states <- args[3]
clusters <- args[4]
cohort <- args[5]
method <- args[6]
clustering.params <- args[7]
output_dir <- args[8]

### [READ IN DATA] #####
sce <- assignIdentity(cells, types, states)
pathways <- sce$pathways
sce <- sce$sce

clusters <- read_csv(clusters) %>% column_to_rownames(var = "id")
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


png(paste0(output_dir, "StatesBoxplot_cohort_", cohort, "_method_", method, "_clusters_", clustering.params, ".png"), height = 13, width = 18, units = "in", res = 100)
ggplot(phenograph.bar, aes(x = clusters, y = activation)) +
  geom_boxplot() +
  facet_wrap("pathway", ncol = 3) +
  astir_paper_theme()
dev.off()
