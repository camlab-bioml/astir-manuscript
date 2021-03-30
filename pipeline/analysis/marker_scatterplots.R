#!/usr/local/bin/Rscript

library(tidyverse)
library(devtools)
library(DelayedArray)
library(ggplot2)
library(ggpubr)
devtools::load_all("../taproom/")

cells <- readRDS("../imc-prototype-analysis/raw_data/zurich1/zurich1_subset_sce.rds")
types <- read_csv("../imc-prototype-analysis/raw_data/zurich1/zurich1_subset_assignments_type.csv") %>% 
  column_to_rownames("X1")
cohort <- "zurich1"

types$cell_type <- get_celltypes(types)

cells$cell_type <- types[colnames(cells), ]$cell_type

unknown.cells <- cells[, cells$cell_type == "Unknown"]
unknown.expr <- t(logcounts(unknown.cells)) %>% as.data.frame()

## Assigned cells
assigned.cells <- cells[, cells$cell_type != "Unknown"]
assigned.expr <- t(logcounts(assigned.cells)) %>% as.data.frame()




max.vals <- list(
  ECadherin = max(max(unknown.expr$`E-Cadherin`), max(assigned.expr$`E-Cadherin`)),
  CD20 = max(max(unknown.expr$CD20), max(assigned.expr$CD20)),
  CD45 = max(max(unknown.expr$CD45), max(assigned.expr$CD45)),
  CD3 = max(max(unknown.expr$CD3), max(assigned.expr$CD3)),
  Vimentin = max(max(unknown.expr$Vimentin), max(assigned.expr$Vimentin)),
  Fibronectin = max(max(unknown.expr$Fibronectin), max(assigned.expr$Fibronectin))
)


CD45_ECadherin.unknown <- ggplot(unknown.expr, aes(x = CD45, y = `E-Cadherin`)) +
  geom_point(alpha = 0.1) +
  xlim(0, max.vals$CD45) + ylim(0, max.vals$ECadherin) +
  geom_density2d() +
  astir_paper_theme() +
  ggtitle(paste("Unknown cells", cohort))

CD3_CD20.unknown <- ggplot(unknown.expr, aes(x = CD3, y = CD20)) +
  geom_point(alpha = 0.1) +
  xlim(0, max.vals$CD3) + ylim(0, max.vals$CD20) +
  geom_density2d() +
  astir_paper_theme() +
  ggtitle(paste("Unknown cells", cohort))

Vimentin_ECadherin.unknown <- ggplot(unknown.expr, aes(x = Vimentin, y = `E-Cadherin`)) +
  geom_point(alpha = 0.1) +
  xlim(0, max.vals$Vimentin) + ylim(0, max.vals$ECadherin) +
  geom_density2d() +
  astir_paper_theme() +
  ggtitle(paste("Unknown cells", cohort))

Fibronectin_ECadherin.unknown <- ggplot(unknown.expr, aes(x = Fibronectin, y = `E-Cadherin`)) +
  geom_point(alpha = 0.1) +
  xlim(0, max.vals$Fibronectin) + ylim(0, max.vals$ECadherin) +
  geom_density2d() +
  astir_paper_theme() +
  ggtitle(paste("Unknown cells", cohort))

Fibronectin_Vimentin.unknown <- ggplot(unknown.expr, aes(x = Fibronectin, y = Vimentin)) +
  geom_point(alpha = 0.1) +
  xlim(0, max.vals$Fibronectin) + ylim(0, max.vals$Vimentin) +
  geom_density2d() +
  astir_paper_theme() +
  ggtitle(paste("Unknown cells", cohort))




CD45_ECadherin.assigned <- ggplot(assigned.expr, aes(x = CD45, y = `E-Cadherin`)) +
  geom_point(alpha = 0.1) +
  xlim(0, max.vals$CD45) + ylim(0, max.vals$ECadherin) +
  geom_density2d() +
  astir_paper_theme() +
  ggtitle(paste("Assigned cells", cohort))

CD3_CD20.assigned <- ggplot(assigned.expr, aes(x = CD3, y = CD20)) +
  geom_point(alpha = 0.1) +
  xlim(0, max.vals$CD3) + ylim(0, max.vals$CD20) +
  geom_density2d() +
  astir_paper_theme() +
  ggtitle(paste("Assigned cells", cohort))

Vimentin_ECadherin.assigned <- ggplot(assigned.expr, aes(x = Vimentin, y = `E-Cadherin`)) +
  geom_point(alpha = 0.1) +
  xlim(0, max.vals$Vimentin) + ylim(0, max.vals$ECadherin) +
  geom_density2d() +
  astir_paper_theme() +
  ggtitle(paste("Assigned cells", cohort))

Fibronectin_ECadherin.assigned <- ggplot(assigned.expr, aes(x = Fibronectin, y = `E-Cadherin`)) +
  geom_point(alpha = 0.1) +
  xlim(0, max.vals$Fibronectin) + ylim(0, max.vals$ECadherin) +
  geom_density2d() +
  astir_paper_theme() +
  ggtitle(paste("Assigned cells", cohort))

Fibronectin_Vimentin.assigned <- ggplot(assigned.expr, aes(x = Fibronectin, y = Vimentin)) +
  geom_point(alpha = 0.1) +
  xlim(0, max.vals$Fibronectin) + ylim(0, max.vals$Vimentin) +
  geom_density2d() +
  astir_paper_theme() +
  ggtitle(paste("Assigned cells", cohort))


ggarrange(CD45_ECadherin.unknown, CD3_CD20.unknown, Vimentin_ECadherin.unknown,
          Fibronectin_ECadherin.unknown, Fibronectin_Vimentin.unknown,
          CD45_ECadherin.assigned, CD3_CD20.assigned,
          Vimentin_ECadherin.assigned, Fibronectin_ECadherin.assigned,
          Fibronectin_Vimentin.assigned,
          nrow = 2, ncol = 5)



#### NEXT TRY

unk <- unknown.expr[unknown.expr$CD45 > 0.25,]
unk$type <- "unknown"
assig <- assigned.expr[assigned.expr$CD45 > 0.25,]
assig$type <- "assigned"

total <- rbind(unk, assig)
total$type <- as.factor(total$type)

ggplot(total, aes(x = `E-Cadherin`, color = type, fill = type)) +
  geom_histogram(position = "identity", bins = 100)

