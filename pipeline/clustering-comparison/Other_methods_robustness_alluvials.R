library(ggalluvial)
library(tidyverse)
library(devtools)
library(taproom)
source("scripts/functions.R")
#devtools::load_all("../taproom/")

#read in all 
args <- commandArgs(trailingOnly = TRUE)
rem_none <- args[1] %>% read_csv()
#rem_none <- read_csv("output/v6/results/GSVA-Other-methods-robustness_basel_FlowSOM_k12_all-markers_None.csv")
rem_stromal <- args[2] %>% read_csv()
#rem_stromal <- read_csv("output/v6/results/GSVA-Other-methods-robustness_basel_FlowSOM_k12_all-markers_Stromal.csv")
rem_macrophages <- args[3] %>% read_csv()
#rem_macrophages <- read_csv("output/v6/results/GSVA-Other-methods-robustness_basel_FlowSOM_k12_all-markers_Stromal-Macrophage.csv")
rem_endothelial <- args[4] %>% read_csv()
#rem_endothelial <- read_csv("output/v6/results/GSVA-Other-methods-robustness_basel_FlowSOM_k12_all-markers_Stromal-Macrophage-Endothelial.csv")
cohort <- args[5]
method <- args[6]
clusters <- args[7]
markers <- args[8]
output_dir <- args[9]

rem <- bind_rows(rem_none, rem_stromal, rem_macrophages, rem_endothelial) %>% 
  mutate(condition = factor(condition, levels = c("None", "Stromal", 
                                                 "Stromal-Macrophage",
                                                 "Stromal-Macrophage-Endothelial"))) %>%
  mutate(Manual_cell_type = ifelse(Manual_cell_type == "T cell", "T cells", Manual_cell_type)) %>%
  mutate(Manual_cell_type = ifelse(Manual_cell_type == "B cell", "B cells", Manual_cell_type))

keep_cells <- rem %>% filter(condition == "None") %>% 
  filter(Manual_cell_type == "Stromal" | 
           Manual_cell_type == "Endothelial" |
           Manual_cell_type == "Macrophage") %>% pull(id)

if(length(keep_cells) == 0){
  keep_cells <- unique(rem$id)
}


# Plot
pdf(paste0(output_dir, "GSVA-robustness-alluv_", cohort, "_", method, "_", clusters, "_", markers, ".pdf"), height = 4.3, width = 8.5)
rem %>% filter(id %in% keep_cells) %>% 
  mutate(condition = str_replace_all(condition, "-", "\n")) %>%
  ggplot(aes(x = condition, stratum = Manual_cell_type, alluvium = id, 
             fill = Manual_cell_type)) +
  geom_stratum(color = "black", alpha = 0.5) +
  stat_flow() + 
  labs(y = "Cells", x = "Cell type(s) for which markers were removed",
       fill = "Assigned\ncell type") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_manual(values = jackson_basel_colours()) +
  astir_paper_theme() +
  theme(axis.ticks.x = element_blank(),
        #axis.ticks.y = element_blank(),
        #axis.text.y = element_blank(),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        panel.background = element_blank())
dev.off()
