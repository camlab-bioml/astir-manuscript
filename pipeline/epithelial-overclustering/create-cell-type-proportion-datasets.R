#!/usr/local/bin/Rscript

library(tidyverse)
library(devtools)
library(DelayedArray)
devtools::load_all("../taproom/")

args <- commandArgs(trailingOnly = TRUE)
astir_assignments <- args[1]
percentage_luminal <- as.numeric(args[2])
output_dir <- args[3]

# Read in data & remove other, Unknown & basal cells
assignments <- read_csv(astir_assignments) %>% 
  column_to_rownames("X1")

assignments$cell_type <- get_celltypes(assignments, 0.95)
assignments <- filter(assignments, 
                      !cell_type %in% c("Unknown", "Other", "Epithelial (basal)",
                             "B cells", "Epithelial (other)", "Endothelial"))

# which cell type has the lowest number of cells? (this will be the max # of
# cells for any group if I want to retain equal proportions)
min_cells <- table(assignments$cell_type) %>% as.numeric() %>% min()

# How many cell types do we have represented?
cell_type_number <- length(unique(assignments$cell_type))

# total number of cells per experiment
cell_number_per_sample <- min_cells * cell_type_number


get_cell_names <- function(type, number){
  assignments %>% 
    filter(cell_type == type) %>% 
    dplyr::slice(1:number) %>% 
    rownames_to_column("cell_ID") %>%
    pull("cell_ID")
}


### CREATE ANY SAMPLE FUNCTION #########
cell_types_wout_luminal <- c("Macrophage", "Stromal", "T cells")

create_sample <- function(percentage_luminal){
  # Calculate number of cells for each type
  luminal_no <- cell_number_per_sample * percentage_luminal
  other_no <- round((cell_number_per_sample - luminal_no) / length(cell_types_wout_luminal))
  
  # Create subset of non luminal cells
  other_cells <- lapply(cell_types_wout_luminal, get_cell_names, number = other_no) %>% 
    bind_cols()
  colnames(other_cells) <- cell_types_wout_luminal
  
  other_cells <- other_cells %>% 
    pivot_longer(everything(), names_to = "type", values_to = "cell_id")
  
  total <- rbind(other_cells, 
                 data.frame(type = rep("Epithelial (luminal)", luminal_no),
                            cell_id = get_cell_names("Epithelial (luminal)", luminal_no)))
  
  total
}


luminal_99 <- create_sample(percentage_luminal)
write_csv(luminal_99, paste0(output_dir, "luminal-", percentage_luminal, ".csv"))