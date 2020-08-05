#!/usr/local/bin/Rscript

library(ggplot2)
library(devtools)
library(viridis)
library(tidyverse)
devtools::load_all("~/taproom")

args <- commandArgs(trailingOnly = TRUE)

tsne <- readRDS(args[1])
basel_states <- args[2] %>% read_csv()
zurich1_states <- args[3] %>% read_csv()
wagner_states <- args[4] %>% read_csv()
a <- as.numeric(args[5])

output_dir <- args[6]


basel_tsne <- left_join(filter(tsne, cohort == "Basel"), basel_states, by = c("cell_id" = "X1")) %>% 
  select(`tSNE Dim1`, `tSNE Dim2`, cell_type, cohort, mTOR_signalling) 

zurich1_tsne <- left_join(filter(tsne, cohort == "Zurich1"), zurich1_states, by = c("cell_id" = "X1")) %>% 
  select(`tSNE Dim1`, `tSNE Dim2`, cell_type, cohort, mTOR_signalling)

wagner_tsne <- left_join(filter(tsne, cohort == "Wagner"), wagner_states, by = c("cell_id" = "X1")) %>% 
  select(`tSNE Dim1`, `tSNE Dim2`, cell_type, cohort)
  
tsne_plot <- rbind(basel_tsne, zurich1_tsne)


range01 <- function(x){(x - min(x)) / (max(x) - min(x))}

tsne_plot <- tsne_plot %>%
  group_by(cohort) %>%
  group_by(cell_type) %>% 
  mutate(rescaled_mTOR = range01(mTOR_signalling))

tsne_plot <- tsne_plot %>% 
  group_by(cell_type) %>% 
  arrange(mTOR_signalling, .by_group = TRUE) %>% 
  ungroup()

tsne_plot$cell_type <- factor(tsne_plot$cell_type, 
                              levels = c("Epithelial (luminal)", "Epithelial (basal)",
                                         "Epithelial (other)", "Stromal", "Fibroblasts",
                                         "T cells", "B cells", "Macrophage", 
                                         "Endothelial", "Other", "Unknown"))
  

png(paste0(output_dir, "tSNE_mTOR-", a, ".png"), width = 1700, height = 283)
ggplot(NULL, aes(x = `tSNE Dim1`, y = `tSNE Dim2`)) +
  geom_point(data = wagner_tsne[,1:2], color = "grey90", alpha = a) +
  geom_point(data = select(tsne_plot, -cell_type), colour = "grey90", alpha = a) +
  geom_point(data = filter(tsne_plot, cell_type != "Other", cell_type != "Unknown"), 
             aes(color = rescaled_mTOR), alpha = a) +
  scale_color_viridis(name = "mTOR \nSignalling") +
  astir_paper_theme() +
  facet_wrap(~cell_type, ncol = 8) +
  theme(panel.spacing.x = unit(7, "mm")) +
  coord_fixed()
dev.off()




