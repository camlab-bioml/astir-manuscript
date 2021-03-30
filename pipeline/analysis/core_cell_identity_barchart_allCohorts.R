#!/usr/local/bin/Rscript
library(ggplot2)
library(tidyverse)
library(wesanderson)
library(cowplot)
library(devtools)
library(DelayedArray)

devtools::load_all("../taproom/")

### [READ IN ARGUMENTS] #####
args <- commandArgs(trailingOnly = TRUE)

basel_type <- args[1] %>% read_csv()
schapiro_type <- args[2] %>% read_csv()
wagner_type <- args[3] %>% read_csv()
zurich1_type <- args[4] %>% read_csv()

output_dir <- args[5]


### [PROCESS DATA] #####
process_celltypes <- function(t){
  cells <- t %>% 
    mutate(cell_type = get_celltypes(select(., -X1))) %>% 
    select(X1, cell_type) %>% 
    dplyr::rename(core = X1) %>% 
    mutate(core = sub("_[^_]+$", "", core))
  
  cells$cell_type[cells$cell_type == "B cell"] <- "B cells"
  cells$cell_type[cells$cell_type == "T cell"] <- "T cells"
  cells$cell_type[cells$cell_type == "Fibroblasts"] <- "Stromal"
  
  coreFreqs <- cells %>% 
    group_by(core) %>% 
    dplyr::count(cell_type)
  
  # create plotting order
  coreOrder <- coreFreqs %>% 
    filter(grepl("Epithelial", cell_type)) %>% 
    group_by(core) %>% 
    summarise(total = sum(n)) %>% 
    arrange(desc(total)) %>% 
    pull(core)
  
  nonEpithelialCores <- unique(cells$core)[!unique(cells$core) %in% coreOrder]
  
  coreOrder <- c(coreOrder, nonEpithelialCores)
  
  # change plotting order
  coreFreqs$core <- factor(coreFreqs$core, levels = coreOrder)
  
  coreFreqs
}

basel <- process_celltypes(basel_type) %>% 
  mutate(Cohort = "Basel")
schapiro <- process_celltypes(schapiro_type) %>% 
  mutate(Cohort = "Schapiro")
wagner <- process_celltypes(wagner_type) %>% 
  mutate(Cohort = "Wagner")
zurich1 <- process_celltypes(zurich1_type) %>% 
  mutate(Cohort = "Zurich")


cores <- rbind(basel, schapiro, wagner, zurich1)

### [PLOT] #####
pdf(paste0(output_dir, "allCohort_core_cell_identity_barchart.pdf"), width = 14, height = 5)
  h1 <- ggplot(cores, aes(x = core, y = n, fill = cell_type)) + 
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = jackson_basel_colours()) +
    scale_y_continuous(expand=c(0,0)) +
    labs(fill = "Cell type") +
    ylab("Frequency") + xlab("Core") +
    astir_paper_theme() +
    theme(axis.title.x = element_text(size = 20),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 20),
          axis.ticks = element_blank(),
          legend.position=c(0.1, 0.8))
  
  h2 <- ggplot(cores, aes(x = core, y = n, fill = Cohort)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = cohort_colours()) +
    astir_paper_theme() +
    xlab("") + ylab("") + 
    theme(plot.margin = unit(c(1,0,-1.5,0), "lines"),
          legend.position=c(0.1, 0.7),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank())
  
  legend <- cowplot::plot_grid(get_legend(h2), get_legend(h1), nrow = 2, axis = "l")
  
  h1 <- h1 + theme(legend.position = "none")
  h2 <- h2 + theme(legend.position = "none")
  plot <- plot_grid(h2, h1, align = "v", ncol = 1, axis = "tb", rel_heights = c(1.8, 20))
  
  plot_grid(plot, legend, nrow = 1, rel_widths = c(10, 1.8))
dev.off()
