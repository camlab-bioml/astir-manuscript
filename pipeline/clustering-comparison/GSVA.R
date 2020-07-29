#!/usr/local/bin/Rscript

### [LOAD LIBRARIES] #####
library(devtools)
library(SingleCellExperiment)
library(tidyverse)
library(GSVA)
library(ggalluvial)
source("scripts/functions.R")
devtools::load_all("../taproom/")

### [READ IN DATA] #####
args <- commandArgs(trailingOnly = TRUE)
cells <- args[1]
markers <- args[2] %>% read_markers()$cell_types
csvs <- args[3]
cohort <- args[4]
output_dir <- args[5]

files <- unlist(strsplit(csvs, split = " "))

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
get_enrichment <- function(cells, markers_list, condition){
  gsva <- gsva(logcounts(cells),
               markers_list,
               method="gsva")
  
  df_gsva <- t(gsva) %>% 
    as.data.frame() %>% 
    rownames_to_column("id") %>% 
    gather(pathway, score, -id) %>% 
    as_tibble()
  
  methods <- inner_join(methodOutputs, df_gsva, by = "id")
  
  methods.summary <- group_by(methods, pathway, method, cluster) %>% 
    dplyr::summarize(mean_score = mean(score)) %>% 
    ungroup() %>% 
    group_by(pathway, method) %>% 
    mutate(mean_score = (mean_score - mean(mean_score)) / sd(mean_score)) %>% 
    ungroup() %>% 
    group_by(cluster, method) %>% 
    mutate(thresh_score = mean_score == max(mean_score)) %>% 
    filter(thresh_score == TRUE)
  
  cell_assigned <- left_join(select(methods, -score), 
                             select(methods.summary, -mean_score), 
                             by = c("method", "pathway", "cluster")) %>% 
    filter(thresh_score == TRUE) %>% 
    select(-thresh_score)

  colnames(cell_assigned) <- c("id", "cluster", "method", "cell_type")
  cell_assigned$condition <- condition
  
  cell_assigned
}

# Define markers 
markers_stromal <- markers[names(markers) != "Stromal"]
markers_macrophage <- markers_stromal[names(markers_stromal) != "Macrophage"]
markers_endothelial <- markers_macrophage[names(markers_macrophage) != "Endothelial"]

# Run GSVA
all_markers <- get_enrichment(cells, markers, "None")
stromal <- get_enrichment(cells, markers_stromal, "Stromal")
stromal_macrophage <- get_enrichment(cells, markers_macrophage, "Stromal & macrophage")
s_m_endothelial <- get_enrichment(cells, markers_endothelial, "Stromal, Macrophage & Endothelial")

all <- rbind(all_markers, stromal, stromal_macrophage, s_m_endothelial)

write_csv(all, paste0(output_dir, "Other_approaches_GSVA_", cohort, ".csv"))
