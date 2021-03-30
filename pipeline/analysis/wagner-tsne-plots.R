library(tidyverse)
library(devtools)
devtools::load_all("../taproom/")

markers <- read_markers("markers/wagner-2019-markers-v4.yml")

tsne <- read_csv("output/squirrel/results/wagner_tsne_expression-subset.csv") %>% 
  as.data.frame()

### Cell type colouring
ggplot(tsne, aes(x = `tSNE Dim1`, y = `tSNE Dim2`)) +
  geom_point(aes(color = cell_type, alpha = 0.001)) +
  scale_color_manual(values = jackson_basel_colours(), name = "Cell type") +
  scale_alpha(guide = 'none') +
  astir_paper_theme()

tsne_longer <- tsne %>% 
  pivot_longer(-c(`tSNE Dim1`, `tSNE Dim2`, id, cell_type), values_to = "expression", names_to = "marker") %>% 
  filter(marker %in% unique(unlist(markers$cell_types)))


### Marker colouring
ggplot(tsne_longer, aes(x = `tSNE Dim1`, y = `tSNE Dim2`)) +
  geom_point(aes(color = expression, alpha = 0.001)) +
  scale_color_viridis_c(name = "Expression") +
  #scale_color_manual(values = jackson_basel_colours(), name = "Cell type") +
  scale_alpha(guide = 'none') +
  astir_paper_theme() +
  coord_fixed() +
  facet_wrap(~marker, nrow = 2)
