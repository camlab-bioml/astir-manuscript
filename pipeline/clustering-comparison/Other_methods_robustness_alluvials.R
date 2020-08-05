library(ggalluvial)
library(tidyverse)
library(devtools)
devtools::load_all("~/taproom/")

#read in all
args <- commandArgs(trailingOnly = TRUE)
all <- args[1] %>% read_csv()
method <- args[2]
cohort <- args[3]
output_dir <- args[4]

# Order factors
all$condition <- factor(all$condition, levels = c("None", "Stromal", 
                                                  "Stromal & Macrophage", 
                                                  "Stromal, Macrophage & Endothelial"))

# Filter to only include one method
plot_all <- all %>% filter(method == method)

# Only select relevant cells
cells_stromal_none <- filter(plot_all, 
                             condition == "none", 
                             cell_type %in% c("Stromal", "Macrophage", "Endothelial")) %>% 
  .$id

# Plot
pdf(paste0(output_dir, "Other_methods_robustness_", method, "_", cohort, ".pdf"), height = 4)
filter(plot_all, id %in% cells_stromal_none) %>% 
  ggplot(aes(x = condition, stratum = cell_type, alluvium = id, 
             fill = cell_type)) +
  geom_stratum(alpha = .5) +
  stat_flow() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_manual(values = jackson_basel_colours()) +
  astir_paper_theme() +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
  labs(y = "Cells", x = "Cell type removed", fill = "Cell type assigned")
dev.off()
