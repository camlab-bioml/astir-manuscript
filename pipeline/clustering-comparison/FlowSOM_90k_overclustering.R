library(tidyverse)
library(SingleCellExperiment)
library(GSVA)
library(devtools)
library(RColorBrewer)
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
  #select(-hclust_cluster) %>% 
  mutate(FlowSOM_cluster = as.numeric(FlowSOM_cluster))


# Combine everything
assigned_clustering <- left_join(assigned_clustering, manual_annotation, 
                                 by = c("cluster" = "FlowSOM_cluster"))


if(snakemake@wildcards[['markers']] == 'all_markers'){
  hand_annotation <- tibble(hclust_cluster = as.character(seq(1:8)),
                            Hand_annotation = c('Epithelial (luminal)', 'Epithelial (luminal)', 'Epithelial (luminal)',
                              'Epithelial (basal)', 'T cells', 'Macrophage', 'Epithelial (luminal)', 'Stromal'))
}else if(snakemake@wildcards[['markers']] == 'specified_markers'){
  hand_annotation <- tibble(hclust_cluster = as.character(seq(1:8)),
                            Hand_annotation = c('Epithelial (luminal)', 'Macrophage', 'Stromal', 'Epithelial (luminal)',
                              'Epithelial (basal)', 'Epithelial (luminal)', 'Stromal', 'Epithelial (luminal)'))
}

head(assigned_clustering)

hand_annotation

left_join(assigned_clustering, hand_annotation)

save_assignments <- left_join(assigned_clustering, hand_annotation) %>%
  dplyr::select(id, Hand_annotation)
#save_assignments <- dplyr::select(assigned_clustering, id, hand_annotation)

save_assignments
write_csv(save_assignments, snakemake@output[['csv']])


sce$`FlowSOM cluster` <- assigned_clustering$hclust_cluster


# Create heatmap

lineage_markers <- markers$cell_types %>% 
  unlist() %>% unique()

rainbow_cols <- function(n){
  cluster_col <- brewer.pal(n = n, name = "Set1")
  names(cluster_col) <- 1:n
  cluster_col
}

pdf(snakemake@output[['heatmap']], width = 11, height = 7)
  lc <- t(as.matrix(assay(sce[lineage_markers, ], 'logcounts')))
  lc <- scale(lc)

  lc[lc > 2] <- 2
  lc[lc < -2] <- -2
  
  cell_types = colData(sce)[[ 'FlowSOM cluster' ]] %>% as.character()

  celltype_annot <- HeatmapAnnotation(`Cell type` = cell_types, 
                                      which="column",
                                      annotation_legend_param = list(title = "Metacluster"),
                                      col = list(`Cell type` = rainbow_cols(8)))  
  
  type_exprs <- Heatmap(t(lc), 
                        name = "Expression",
                        column_title = "Cell",
                        col=viridis(100),
                        top_annotation = celltype_annot,
                        show_column_names = FALSE,
                        column_order = order(cell_types))
  type_exprs
  #createHeatmap(sce[lineage_markers, ], cell_type_column = 'FlowSOM cluster')
dev.off()
