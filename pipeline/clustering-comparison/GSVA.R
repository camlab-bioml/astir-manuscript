#!/usr/local/bin/Rscript

### [LOAD LIBRARIES] #####
library(devtools)
library(SingleCellExperiment)
library(scater)
library(tidyverse)
library(GSVA)
library(ggalluvial)
source("scripts/functions.R")
devtools::load_all("~/taproom/")

### [READ IN DATA] #####
args <- commandArgs(trailingOnly = TRUE)
rem_none <- args[1] %>% read_csv()
rem_stromal <- args[2] %>% read_csv()
rem_macrophages <- args[3] %>% read_csv()
rem_endothelial <- args[4] %>% read_csv()
cohort <- args[5]
method <- args[6]
clusters <- args[7]
markers <- args[8]
output_dir <- args[9]

rem <- bind_rows(rem_none, rem_stromal, rem_macrophages, rem_endothelial)

methodOutputs <- lapply(files, function(f) {
  cluster_df <- read_csv(f)
  names(cluster_df)[2] <- "cluster"
  cluster_df$cluster <- as.character(cluster_df$cluster)
  cluster_df$method <- paste0(cluster_df$method, "-", cluster_df$params)
  cluster_df$params <- NULL
  cluster_df
}) %>% 
  bind_rows()


### [ANALYSIS] #####


# Define markers 
final_markers <- markers[!names(markers) %in% markers_remove]
# markers_stromal <- markers[names(markers) != "Stromal"]
# markers_macrophage <- markers_stromal[names(markers_stromal) != "Macrophage"]
# markers_endothelial <- markers_macrophage[names(markers_macrophage) != "Endothelial"]

# Run GSVA
gsva <- get_enrichment(cells, final_markers, paste(markers_remove, collapse = "_"))
# stromal <- get_enrichment(cells, markers_stromal, paste(markers_rem, collapse = "_"))
# stromal_macrophage <- get_enrichment(cells, markers_macrophage, paste(markers_rem, collapse = "_"))
# s_m_endothelial <- get_enrichment(cells, markers_endothelial, paste(markers_rem, collapse = "_"))

#all <- rbind(all_markers, stromal, stromal_macrophage, s_m_endothelial)
 
write_csv(gsva, paste0(output_dir, "Other_approaches_GSVA_", cohort, "_", paste(markers_remove, collapse = "_"), ".csv"))
