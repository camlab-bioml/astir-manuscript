library(SingleCellExperiment)
library(ComplexHeatmap)
library(devtools)
library(tidyverse)
library(RColorBrewer)
devtools::load_all("../taproom/")

sce <- readRDS(snakemake@input[['sce']])#"output/phoenix/sces/zurich1_sce.rds")
clusters <- read_csv(snakemake@input[['clusters']])#"output/phoenix/results/GSVA-assignment-Phenograph-zurich1-all_markers-k40.csv")



sce$cell_type <- clusters$Manual_cell_type
sce$cluster <- clusters$cluster



interesting_markers <- c("E-Cadherin", "c-Myc", 
                         "phospho S6", "phospho mTOR", "Cytokeratin 8/18", 
                         "Cytokeratin 19", "pan Cytokeratin",
                         "Cytokeratin 7", "Twist", "Vimentin", "SMA")

exclude_clusters <- c(61, 8, 17, 52, 25, 6, 10, 60, 50)

epi_sce <- sce[interesting_markers, (sce$cell_type == "Epithelial (luminal)" &
                                       !(sce$cluster %in% exclude_clusters))]

lc <- t(as.matrix(assay(epi_sce, 'logcounts')))
lc <- scale(lc)

lc[lc > 2] <- 2
lc[lc < -2] <- -2

cluster = colData(epi_sce)[[ 'cluster' ]]

rainbow_cols <- function(){
  cluster_no <- unique(cluster) %>% length()
  cluster_col <- brewer.pal(n = cluster_no, name = "Set1")
  names(cluster_col) <- unique(cluster)
  cluster_col
}

celltype_annot <- HeatmapAnnotation(`Phenograph Cluster` = as.character(cluster),
                                    which="column",
                                    annotation_legend_param = list(title = "Phenograph\nCluster"),
                                    col = list(`Phenograph Cluster` = rainbow_cols()))  

type_exprs <- Heatmap(t(lc), 
                      name = "Z-Scaled\nExpression",
                      column_title = "Cell",
                      col=viridis(100),
                      top_annotation = celltype_annot,
                      show_column_names = FALSE,
                      column_order = order(cluster))

pdf(snakemake@output[['heatmap']], height = 3.5, width = 8)
type_exprs
dev.off()



