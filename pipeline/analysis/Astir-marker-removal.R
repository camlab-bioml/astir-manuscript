library(tidyverse)
library(devtools)
library(DelayedArray)
library(ggalluvial)
devtools::load_all("../taproom/")

create_assignments <- function(path, markers){
  print(path)
  df <- read_csv(path)
  print(df)
  cell_names <- df$X1
  df$X1 <- NULL
  cell_types <- get_celltypes(df)
  
  tibble(id = cell_names,
         cell_type = cell_types,
         markers_removed = markers)
}

basel_all <- create_assignments(snakemake@input[['none']], 
                                "None")
basel_Cyto7 <- create_assignments(snakemake@input[['cyto7']],
                                  "Cytokeratin 7")
basel_Cyto7_19 <- create_assignments(snakemake@input[['cyto7_19']],
                                     "Cytokeratin 7\nCytokeratin 19")
basel_Cyto7_19_18 <- create_assignments(snakemake@input[['all_luminal']],
                                        "All luminal\nmarkers")

marker_removal <- bind_rows(basel_all, basel_Cyto7, basel_Cyto7_19, basel_Cyto7_19_18) %>% 
  mutate(markers_removed = factor(markers_removed, levels = c("None", "Cytokeratin 7",
                                             "Cytokeratin 7\nCytokeratin 19",
                                             "All luminal\nmarkers")))


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

#sce <- readRDS(snakemake@input[['sce']])
#pdf(snakemake@output[['heatmap']])