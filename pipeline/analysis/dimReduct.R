#!/usr/local/bin/Rscript

library(SingleCellExperiment)
library(tidyverse)
library(scater)
library(rsvd)
library(devtools)
devtools::load_all("~/taproom/")
source("~/imc-2020/scripts/functions.R")
source('~/FIt-SNE/fast_tsne.R', chdir=T)

args <- commandArgs(trailingOnly = TRUE)

basel_cells <- args[1]
basel_types <- args[2]
#basel_states <- args[3]

zurich1_cells <- args[3]
zurich1_types <- args[4]
#zurich1_states <- args[6]

wagner_cells <- args[5]
wagner_types <- args[6]
#wagner_states <- args[9]

output_dir_res <- args[7]


### [GET SCE OBJECTS] #####
#basel <- assignIdentity(basel_cells, basel_types, basel_states)$sce
basel <- readRDS(basel_cells)
assays(basel)$raw_imc <- NULL
assays(basel)$logcounts_unwinsorized <- NULL

basel_type <- read_csv(basel_types)

basel_type$cell_type <- get_celltypes(select(basel_type, -X1)) 
basel_type <- select(basel_type, X1, cell_type) %>%
  column_to_rownames("X1")

# basel.cellType <- basel$cell_type
colData(basel) <- NULL
# basel$cell_type <- basel.cellType
basel$cohort <- "Basel"
colData(basel)["cell_type"] <- basel_type[colnames(basel),]


#zurich1 <- assignIdentity(zurich1_cells, zurich1_types, zurich1_states)$sce
zurich1 <- readRDS(zurich1_cells)
assays(zurich1)$raw_imc <- NULL
assays(zurich1)$logcounts_unwinsorized <- NULL

zurich1_type <- read_csv(zurich1_types)
zurich1_type$cell_type <- get_celltypes(select(zurich1_type, -X1)) 
zurich1_type <- select(zurich1_type, X1, cell_type) %>%
  column_to_rownames("X1")

# zurich1.cellType <- zurich1$cell_type
colData(zurich1) <- NULL
# zurich1$cell_type <- zurich1.cellType
zurich1$cohort <- "Zurich1"
colData(zurich1)["cell_type"] <- zurich1_type[colnames(zurich1),]

#wagner <- assignIdentity(wagner_cells, wagner_types, wagner_states)$sce
wagner <- readRDS(wagner_cells)
assays(wagner)$raw_imc <- NULL
assays(wagner)$logcounts_unwinsorized <- NULL

wagner_type <- read_csv(wagner_types)
wagner_type$cell_type <- get_celltypes(select(wagner_type, -X1)) 
wagner_type <- select(wagner_type, X1, cell_type) %>%
  column_to_rownames("X1")

# wagner.cellType <- wagner$cell_type
colData(wagner) <- NULL
# wagner$cell_type <- wagner.cellType
wagner$cohort <- "Wagner"
colData(wagner)["cell_type"] <- wagner_type[colnames(wagner),]

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
expression$cohort <- as.factor(sce[, rownames(expression)]$cohort)

#expression <- expression %>% as.matrix()

tsne <- fftRtsne(expression[, 1:length(clustering.markers)] %>% as.matrix()) %>% 
  as.data.frame()
colnames(tsne) <- c("tSNE Dim1", "tSNE Dim2")

tsne$cell_type <- expression$cell_type
tsne$cohort <- expression$cohort
tsne$cell_id <- rownames(expression)

# add cohort & cell type
saveRDS(tsne, file = paste0(output_dir_res, "tsne_results.rds"))