library(SingleCellExperiment)
library(ComplexHeatmap)
library(devtools)
library(tidyverse)
library(RColorBrewer)
devtools::load_all("../taproom/")


### Read in data

sce <- readRDS("output/phoenix/sces/zurich1_sce.rds")
clusters <- read_csv("output/phoenix/results/GSVA-assignment-Phenograph-zurich1-all_markers-k40.csv")


sce$cell_type <- clusters$Manual_cell_type
sce$cluster <- clusters$cluster

# Scale data
lc <- t(as.matrix(assay(sce, 'logcounts')))
lc <- scale(lc)
lc[lc > 2] <- 2
lc[lc < -2] <- -2


# Select markers and clusters
# Clusters to remove
exclude_clusters <- c(61, 8, 17, 52, 25, 6, 10, 60, 50)

# Get all cells I want to plot
interesting_cells <- clusters %>% 
  filter(Manual_cell_type == "Epithelial (luminal)" &
           !(cluster %in% exclude_clusters))

interesting_cells_id <- interesting_cells$id

# Select interesting clusters
interesting_markers <- c("E-Cadherin", "c-Myc", 
                         "phospho S6", "phospho mTOR", "Cytokeratin 8/18", 
                         "Cytokeratin 19", "pan Cytokeratin",
                         "Cytokeratin 7", "Twist", "Vimentin", "SMA")

interesting_expression_mat <- lc[interesting_cells_id, interesting_markers]

# Get clusters for each cell
clusters <- interesting_cells$cluster

# Define heatmap color function 
cols <- function(){
  cluster_no <- unique(clusters) %>% length()
  cluster_col <- brewer.pal(n = cluster_no, name = "Set1")
  names(cluster_col) <- unique(clusters) %>% as.character()
  cluster_col
}


# Create annotation
celltype_annot <- HeatmapAnnotation(`Phenograph Cluster` = as.character(clusters),
                                    which="column",
                                    col = list(cluster = cols()),
                                    legend_direction = "horizontal")  

# Create heatmap
type_exprs <- Heatmap(t(interesting_expression_mat), 
                      name = "Z-Scaled Expression",
                      column_title = "Cell",
                      col=viridis(100),
                      top_annotation = celltype_annot,
                      show_column_names = FALSE,
                      column_order = order(clusters))



type_exprs



# expression change
# 1. cluster = phenograph cluster
# 2. expression = zscaled expression
# 3. annotations above each other
# 4. cluster colors