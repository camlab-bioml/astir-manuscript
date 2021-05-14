library(tidyverse)
library(devtools)
library(DelayedArray)
library(ggalluvial)
devtools::load_all("../taproom/")

create_assignments <- function(path, markers){
  df <- read_csv(path)
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
basel_Cyto7_19_18 <- create_assignments(snakemake@input[['cyto7_19_8-18']],
                                        "Cytokeratin 7\nCytokeratin 19\nCytokeratin 8/18")
basel_Cyto7_19_18_pan <- create_assignments(snakemake@input[['cyto7_19_8-18_pan']],
                                            "Cytokeratin 7\nCytokeratin 19\nCytokeratin 8/18\nPan Cytokeratin")


marker_removal <- bind_rows(basel_all, basel_Cyto7, basel_Cyto7_19, 
                            basel_Cyto7_19_18, basel_Cyto7_19_18_pan) %>% 
  mutate(markers_removed = factor(markers_removed, levels = c("None", "Cytokeratin 7",
                                             "Cytokeratin 7\nCytokeratin 19",
                                             "Cytokeratin 7\nCytokeratin 19\nCytokeratin 8/18",
                                             "Cytokeratin 7\nCytokeratin 19\nCytokeratin 8/18\nPan Cytokeratin")))


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
