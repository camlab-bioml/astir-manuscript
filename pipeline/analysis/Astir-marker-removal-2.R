library(tidyverse)
library(devtools)
library(DelayedArray)
library(ggalluvial)
library(SingleCellExperiment)
devtools::load_all("../taproom/")

create_assignments <- function(path, markers){
  print(path)
  df <- read_csv(path)
  print(df)
  cell_names <- df$X1
  df$X1 <- NULL
  cell_types <- get_celltypes(df, snakemake@wildcards[['thresh']])
  
  data.frame(id = cell_names,
             cell_type = cell_types,
             markers_removed = markers)
}

basel_normal <- create_assignments(snakemake@input[['normal']], "None")

basel_rem <- create_assignments(snakemake@input[['removed']], snakemake@wildcards[['marker']])

marker_removal <- bind_rows(basel_normal, basel_rem) %>% 
  mutate(markers_removed = factor(markers_removed, levels = c("None", snakemake@wildcards[['marker']])))

pdf(snakemake@output[['pdf']])
marker_removal %>% 
  ggplot(aes(x = markers_removed, stratum = cell_type, alluvium = id, 
             fill = cell_type)) +
  geom_stratum(color = "black", alpha = 0.5) +
  stat_flow() +
  labs(y = "Cells", x = "Removed Markers", fill = "Assigned\ncell type") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_manual(values = jackson_basel_colours()) +
  astir_paper_theme() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

sce <- readRDS(snakemake@input[['sce']])

# basel_normal <- as.data.frame(basel_normal)
# basel_rem <- as.data.frame(basel_rem)

rownames(basel_normal) <- basel_normal$id
rownames(basel_rem) <- basel_rem$id

sce$cell_types_all <- basel_normal[colnames(sce), 2]
sce$cell_types_rem <- basel_rem[colnames(sce), 2]

markers <- read_markers(snakemake@input[['markers']])
cell_type_markers <- unique(unlist(markers$cell_types))

pdf(snakemake@output[['heatmap']], width = 10)
lc <- t(as.matrix(assay(sce[cell_type_markers,], 'logcounts')))
lc <- apply(lc, 2, function(x) x / max(x))
# lc <- scale(lc)

# thresh <- 2

# lc[lc > thresh] <- thresh
# lc[lc < -thresh] <- -thresh

cell_types_all = colData(sce)[["cell_types_all"]]
cell_types_rem = colData(sce)[['cell_types_rem']]

celltype_annot <- HeatmapAnnotation(`Normal cell type` = cell_types_all, 
                                    `Removed cell type` = cell_types_rem,
                                    which="column",
                                    col = list(`Normal cell type` = jackson_basel_colours(),
                                               `Removed cell type` = jackson_basel_colours()))  

type_exprs <- Heatmap(t(lc), 
                      name = "Expression",
                      column_title = "Cell",
                      col=viridis(100),
                      top_annotation = celltype_annot,
                      show_column_names = FALSE,
                      column_order = order(paste(cell_types_all, cell_types_rem)))
type_exprs

dev.off()