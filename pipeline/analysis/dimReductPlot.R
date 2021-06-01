#!/usr/local/bin/Rscript 
library(ggplot2)
library(dplyr)
library(devtools)
library(wesanderson)
library(cowplot)
library(taproom)
source("scripts/functions.R")

args <- commandArgs(trailingOnly = TRUE)

rds <- args[1]
a <- as.numeric(args[2])
output_dir <- args[3]

tsne <- readRDS(rds) %>% sample_frac(1) %>% mutate(cohort = as.character(cohort))

tsne$cohort[tsne$cohort == "Zurich1"] <- "Zurich"

### [SAVE PLOT AS PDF] #####
# save as png
png(paste0(output_dir, "tSNE_cellType-", a, ".png"), width = 2709, height = 2000, res = 300)
tsne_p <- ggplot(tsne, aes(x = `tSNE Dim1`, y = `tSNE Dim2`)) +
  geom_point(aes(color = cell_type, alpha = a)) +
  scale_color_manual(values = jackson_basel_colours(), name = "Cell type") +
  scale_alpha(guide = 'none') +
  astir_paper_theme() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18)) +
  coord_fixed()

tsne_legend <- get_legend(tsne_p)
tsne_p <- tsne_p + theme(legend.position = "none")
plot_grid(tsne_p, tsne_legend, ncol = 2, rel_widths = c(1, 0.42))
dev.off()

png(paste0(output_dir, "tSNE_cohort-", a, ".png"), width = 2709, height = 2000, res = 300)
tsne_p <- ggplot(tsne, aes(x = `tSNE Dim1`, y = `tSNE Dim2`)) +
  geom_point(aes(color = cohort, alpha = a)) +
  scale_color_manual(values = cohort_colours(), name = "Cohort") +
  scale_alpha(guide = 'none') +
  astir_paper_theme() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18)) +
  coord_fixed()

tsne_legend <- get_legend(tsne_p)
tsne_p <- tsne_p + theme(legend.position = "none")
plot_grid(tsne_p, tsne_legend, ncol = 2, rel_widths = c(1, 0.42))
dev.off()