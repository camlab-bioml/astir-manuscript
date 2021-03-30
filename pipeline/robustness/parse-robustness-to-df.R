


suppressPackageStartupMessages({
  library(tidyverse)
  library(devtools)
  library(glue)
  library(ggalluvial)
  library(cowplot)
  library(matrixStats) 
  library(pheatmap)
  library(glue)
  library(taproom)
  library(yaml)
  library(stringdist)
  library(argparse)
})


parser <- ArgumentParser(description = "")

parser$add_argument('--marker_yml', type = 'character',
                    help='Marker yaml')
parser$add_argument('--input_dir', type = 'character',
                    help='Dir of assigned inputs')
parser$add_argument('--cohort', type='character')
parser$add_argument('--output_file', type = 'character',
                    help='Output tsv')


args <- parser$parse_args()

markers <- read_yaml(args$marker_yml)


# Removing cell types

removed_csvs <- dir(args$input_dir, pattern = 'assignments*.*csv', 
                    full.names = TRUE)


read_to_mat <- function(f) {
  df <- read_csv(f)
  df <- rename(df, cell_id = X1)

  condition <- strsplit(f, "-")[[1]]
  condition <- condition[length(condition)]
  condition <- gsub(".csv", "", condition, fixed=TRUE)
 
  df <- mutate(df, condition = condition) %>% 
    select(cell_id, condition, everything())
  
  return(df)
}
 
df <- map_dfr(removed_csvs, read_to_mat)

df$condition <- gsub("_", " ", df$condition)

conditions <- unique(df$condition)

get_celltype_from_df <- function(cond, threshold, df) {
  df2 <- filter(df, condition == cond) 
  
  
  mat <- df2 %>% 
    select(-cell_id, -condition) %>% 
    as.matrix()
  
  mat <- mat[, !is.na(colSums(mat))]
  
  tibble(
    cell_id = df2$cell_id,
    cell_type = taproom::get_celltypes(mat, threshold)
  )
}

th <- 0.95

get_celltype_from_condition <- function(cond) {
  cell_types <- names(markers$cell_types)
  cell_types[which.min(stringdist(cond, cell_types))]
}

get_mean_for_th <- function(th, df, condition) {
  initial_cells <- get_celltype_from_df("None", 0.95, df)
  ct <- get_celltype_from_condition(condition)
  
  initial <- filter(initial_cells, cell_type == ct)$cell_id

  cells_at_th <- get_celltype_from_df(condition, th, df) %>% 
    filter(cell_id %in% initial) 
  cells_at_th <- mean(cells_at_th$cell_type %in% c("Other", "Unknown"))
  
  
  tribble(
    ~ threshold, ~ condition, ~ mean_unknown,
    th, condition, cells_at_th
    )
}


thresholds <- c(0.5, 0.7, 0.9, 0.99, 0.999)


mean_accs <- lapply(setdiff(conditions, "None"), function(cond) {
  map_dfr(thresholds, get_mean_for_th, df, cond)
}) %>% 
  bind_rows()


cmap <- sapply(conditions, get_celltype_from_condition)

mean_accs$condition <- plyr::mapvalues(mean_accs$condition,
                             from = names(cmap), to = cmap)

mean_accs$cohort <- args$cohort

write_tsv(mean_accs, args$output_file)


