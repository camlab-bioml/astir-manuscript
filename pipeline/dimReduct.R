library(SingleCellExperiment)
library(tidyverse)
library(scater)
library(rsvd)
library(devtools)
devtools::load_all("~/taproom")
source("~/imc-2020/scripts/functions.R")
source('~/FIt-SNE/fast_tsne.R', chdir=T)

args <- commandArgs(trailingOnly = TRUE)

basel_cells <- args[1]
basel_types <- args[2]
basel_states <- args[3]

zurich1_cells <- args[4]
zurich1_types <- args[5]
zurich1_states <- args[6]

wagner_cells <- args[7]
wagner_types <- args[8]
wagner_states <- args[9]
output_dir <- args[10]


### [GET SCE OBJECTS] #####
basel <- assignIdentity(basel_cells, basel_types, basel_states)$sce

basel.cellType <- basel$cell_type
colData(basel) <- NULL
basel$cell_type <- basel.cellType
basel$cohort <- "Basel"


zurich1 <- assignIdentity(zurich1_cells, zurich1_types, zurich1_states)$sce

zurich1.cellType <- zurich1$cell_type
colData(zurich1) <- NULL
zurich1$cell_type <- zurich1.cellType
zurich1$cohort <- "Zurich1"

wagner <- assignIdentity(wagner_cells, wagner_types, wagner_states)$sce

wagner.cellType <- wagner$cell_type
colData(wagner) <- NULL
wagner$cell_type <- wagner.cellType
wagner$cohort <- "Wagner"

### [RENAME MARKERS] #####
wagner.names <- wagner %>% rownames()
wagner.names[6] <- "Cytokeratin 8/18"
wagner.names[13] <- "Her2"
wagner.names[28] <- "E-Cadherin"
wagner.names[29] <- "Ki-67"
rownames(wagner) <- wagner.names

### [COMBINE & START CLUSTERING] #####
clustering.markers <- Reduce(intersect, list(basel %>% rownames(),
                                             #schapiro %>% rownames(),
                                             zurich1 %>% rownames(),
                                             wagner %>% rownames()))

sce <- cbind(basel[clustering.markers, ],
             #schapiro[clustering.markers, ],
             wagner[clustering.markers, ],
             zurich1[clustering.markers, ])

expression <- logcounts(sce) %>% 
  t() %>% 
  as.data.frame()

expression$cell_type <- sce[, rownames(expression)]$cell_type
expression$cohort <- sce[, rownames(expression)]$cohort

expression <- expression %>% as.matrix()

tsne <- fftRtsne(expression[, 1:length(clustering.markers)])

### [SAVE PLOT AS PDF] #####
# save as png
png(paste0(output_dir, "tSNE_cellType.png"))
plot(tsne, col = expression$cell_type)
dev.off()

png(paste0(output_dir, "tSNE_cohort.png"))
plot(tsne, col = expression$cohort)
dev.off()
