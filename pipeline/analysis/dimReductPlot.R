#!/usr/local/bin/Rscript
library(ggplot2)
library(devtools)
library(wesanderson)
devtools::load_all("~/taproom")

args <- commandArgs(trailingOnly = TRUE)

rds <- args[1]
a <- as.numeric(args[2])
output_dir <- args[3]

tsne <- readRDS(rds)


### [SAVE PLOT AS PDF] #####
# save as png
png(paste0(output_dir, "tSNE_cellType-", a, ".png"))
ggplot(tsne, aes(x = `tSNE Dim1`, y = `tSNE Dim2`)) +
  geom_point(aes(color = cell_type, alpha = a), shape = 1) +
  scale_color_manual(values = jackson_basel_colours()) +
  scale_alpha(guide = 'none') +
  astir_paper_theme()
dev.off()

png(paste0(output_dir, "tSNE_cohort-", a, ".png"))
ggplot(tsne, aes(x = `tSNE Dim1`, y = `tSNE Dim2`), shape = 1) +
  geom_point(aes(color = cohort, alpha = a)) +
  scale_color_manual(values = cohort_colours()) +
  scale_alpha(guide = 'none') +
  astir_paper_theme()
dev.off()
