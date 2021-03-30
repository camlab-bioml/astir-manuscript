library(tidyverse)
library(devtools)
devtools::load_all("~/taproom/")

args <- commandArgs(trailingOnly = TRUE)
markers <- args[1]
wagner_tsne <- args[2]
output_dir <- args[3]

markers <- read_markers(markers)

tsne <- read_csv(wagner_tsne) %>% 
  as.data.frame()

### Cell type colouring
png(paste0(output_dir, "wagner-tSNE-cellType.png"), width = 2709, height = 2000, res = 300)
  ggplot(tsne, aes(x = `tSNE Dim1`, y = `tSNE Dim2`)) +
    geom_point(aes(color = cell_type, alpha = 0.00001)) +
    scale_color_manual(values = jackson_basel_colours(), name = "Cell type") +
    scale_alpha(guide = 'none') +
    astir_paper_theme()
dev.off()


tsne_longer <- tsne %>% 
  pivot_longer(-c(`tSNE Dim1`, `tSNE Dim2`, id, cell_type), values_to = "expression", names_to = "marker") %>% 
  filter(marker %in% unique(unlist(markers$cell_types)))

png(paste0(output_dir, "wagner-tSNE-expression.png"), width = 15000, height = 2000, res = 300)
  ### Marker colouring
  ggplot(tsne_longer, aes(x = `tSNE Dim1`, y = `tSNE Dim2`)) +
    geom_point(aes(color = expression, alpha = 0.00001)) +
    scale_color_viridis_c(name = "Expression") +
    scale_alpha(guide = 'none') +
    astir_paper_theme() +
    coord_fixed() +
    facet_wrap(~marker, nrow = 2)
dev.off()
