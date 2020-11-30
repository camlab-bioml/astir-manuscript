#!/usr/local/bin/Rscript

library(tidyverse)
library(devtools)
library(DelayedArray)
devtools::load_all("../taproom/")

args <- commandArgs(trailingOnly = TRUE)
astir_assignments <- args[1]
percentage_luminal <- args[2]
output_dir <- args[3]

astir_assignments <- "output/chipmunk/astir_assignments/basel_astir_assignments.csv"

# Read in data & remove other, Unknown & basal cells
assignments <- read_csv(astir_assignments) %>% 
  column_to_rownames("X1")

assignments$cell_type <- get_celltypes(assignments)
assignments <- filter(assignments, 
                      !cell_type %in% c("Unknown", "Other"))

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
    slice(1:number) %>% 
    rownames_to_column("cell_ID") %>%
    pull("cell_ID")
}


### CREATE ANY SAMPLE FUNCTION #########
cell_types_wout_luminal <- c("Epithelial (other)", "Macrophage", "Stromal", 
                             "Endothelial", "T cells", "B cells", "Epithelial (basal)")

create_sample <- function(percent_luminal){
  # Calculate number of cells for each type
  luminal_no <- cell_number_per_sample * percent_luminal
  other_no <- round((cell_number_per_sample - luminal_no) / 4)
  
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



luminal_99 <- create_sample(percent_luminal)
write_csv(luminal_99, file = paste0(output_dir, "luminal-", percentage_luminal, ".csv"))

# luminal_80 <- create_sample(0.8)
# write_csv(luminal_80, file = "output/chipmunk/results/epithelial_overclustering/luminal_80.csv")
# 
# luminal_70 <- create_sample(0.7)
# write_csv(luminal_70, file = "output/chipmunk/results/epithelial_overclustering/luminal_70.csv")
# 
# luminal_60 <- create_sample(0.6)
# write_csv(luminal_60, file = "output/chipmunk/results/epithelial_overclustering/luminal_60.csv")
# 
# luminal_50 <- create_sample(0.5)
# write_csv(luminal_50, file = "output/chipmunk/results/epithelial_overclustering/luminal_50.csv")
# 
# luminal_40 <- create_sample(0.4)
# write_csv(luminal_40, file = "output/chipmunk/results/epithelial_overclustering/luminal_40.csv")
# 
# luminal_30 <- create_sample(0.3)
# write_csv(luminal_30, file = "output/chipmunk/results/epithelial_overclustering/luminal_30.csv")
# 
# luminal_20 <- create_sample(0.2)
# write_csv(luminal_20, file = "output/chipmunk/results/epithelial_overclustering/luminal_20.csv")





