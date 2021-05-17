library(tidyverse)
library(SingleCellExperiment)
library(GSVA)
library(devtools)
devtools::load_all("../taproom/")
source("scripts/functions.R")

sce <- readRDS(snakemake@input[['sce']])

FlowSOM <- read_csv(snakemake@input[['FlowSOM']])
colnames(FlowSOM)[2] <- "cluster"

markers <- read_markers(snakemake@input[['markers']])



expression_mat <- logcounts(sce) %>% 
  t() %>% 
  as_tibble() %>% 
  mutate(id = colnames(sce))

expression_mat <- expression_mat %>% 
  left_join(select(FlowSOM, id, cluster)) %>% 
  select(-id)


aggregate_expression <- aggregate(expression_mat[, 1:(ncol(expression_mat) - 1)],
                                  expression_mat[, ncol(expression_mat)],
                                  mean) %>% 
  column_to_rownames("cluster")


### Hierarchical clustering
cl <- hclust(dist(aggregate_expression))

cluster_cut <- cutree(cl, 8) %>% 
  as.data.frame() %>% 
  rownames_to_column("FlowSOM_cluster") %>% 
  mutate(FlowSOM_cluster = as.character(FlowSOM_cluster))

colnames(cluster_cut)[2] <- "hclust_cluster"
cluster_cut$hclust_cluster = as.character(cluster_cut$hclust_cluster)


hclust_expression <- expression_mat %>% 
  mutate(cluster = as.character(cluster)) %>% 
  left_join(cluster_cut, by = c("cluster" = "FlowSOM_cluster")) %>% 
  select(-cluster)

### Assign cell types
hclust_assignment <- hclust_expression %>% 
  dplyr::rename("cluster" = "hclust_cluster") %>% 
  assign_clusters(markers)


# GSVA cell type
assigned <- hclust_assignment$assignment %>% 
  dplyr::rename("hclust_cluster" = "cluster") %>% 
  left_join(cluster_cut) %>% 
  mutate(FlowSOM_cluster = as.numeric(FlowSOM_cluster)) %>% 
  ungroup()

assigned_clustering <- left_join(FlowSOM, select(assigned, FlowSOM_cluster, GSVA_cell_type),
                                 by = c("cluster" = "FlowSOM_cluster"))


# z-score cell type
manual_annotation <- manual_cluster_assignment(hclust_assignment$expression, markers) %>% 
  dplyr::rename("hclust_cluster" = "cluster") %>% 
  left_join(cluster_cut) %>% 
  ungroup() %>% 
  select(-hclust_cluster) %>% 
  mutate(FlowSOM_cluster = as.numeric(FlowSOM_cluster))


# Combine everything
assigned_clustering <- left_join(assigned_clustering, manual_annotation, 
                                 by = c("cluster" = "FlowSOM_cluster"))

write_csv(assigned_clustering, snakemake@output[['csv']])


sce$`FlowSOM cell type` <- assigned_clustering$Manual_cell_type


# Create heatmap

lineage_markers <- markers$cell_types %>% 
  unlist() %>% unique()

pdf(snakemake@output[['heatmap']])
createHeatmap(sce[lineage_markers, ], cell_type_column = 'FlowSOM cell type')
dev.off()
