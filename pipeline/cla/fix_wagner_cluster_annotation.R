
library(tidyverse)
library(SingleCellExperiment)

get_clusters <- function(f) {
    sce <- readRDS(f)
    tibble(
        cell_id=colnames(sce),
        cluster=sce$cluster
    )
}

files <- snakemake@input[['rds']]

df_cluster <- map_dfr(files, get_clusters)

df_annotation <- read_csv(snakemake@input[['cluster_map']])


df <- inner_join(df_cluster, df_annotation)

write_csv(df, snakemake@output[[1]])




