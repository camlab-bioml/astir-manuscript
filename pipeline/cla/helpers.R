library(yardstick)
library(tidyverse)

acc_wrap <- function(tt) {
  # annotated <- unique(as.character(tt$cell_type_annotated))
  cell_types <- unique(union(tt$cell_type_predicted, tt$cell_type_annotated))
  
  tt$cell_type_annotated <- factor(tt$cell_type_annotated, levels = cell_types)
  tt$cell_type_predicted <- factor(tt$cell_type_predicted, levels = cell_types)
  
  bind_rows(
    tryCatch({kap(tt, cell_type_annotated, cell_type_predicted)}, error=function(e) NULL),
    tryCatch({bal_accuracy(tt, cell_type_annotated, cell_type_predicted)}, error=function(e) NULL),
    tryCatch({mcc(tt, cell_type_annotated, cell_type_predicted)}, error=function(e) NULL),
    tryCatch({sensitivity(tt, cell_type_annotated, cell_type_predicted)}, error=function(e) NULL),
    tryCatch({specificity(tt, cell_type_annotated, cell_type_predicted)}, error=function(e) NULL)
  )
}