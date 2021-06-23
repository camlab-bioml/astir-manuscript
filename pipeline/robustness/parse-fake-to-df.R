


suppressPackageStartupMessages({
  library(tidyverse)
  library(devtools)
  library(glue)
  library(ggalluvial)
  library(cowplot)
  library(matrixStats) 
  library(pheatmap)
  library(glue)
  library(yaml)
  library(stringdist)
  library(argparse)
})

devtools::load_all("../taproom/")


parser <- ArgumentParser(description = "Spatial locations for core")

# parser$add_argument('--marker_yml', type = 'character',
#                     help='Marker yaml')
parser$add_argument('--input_dir', type = 'character',
                    help='Dir of assigned inputs')
parser$add_argument('--cohort', type='character')
parser$add_argument('--output_file', type = 'character',
                    help='Output tsv')


args <- parser$parse_args()

# markers <- read_yaml(args$marker_yml)


# Removing cell types

added_csvs <- dir(args$input_dir, pattern = 'assignments*.*csv', 
                    full.names = TRUE)

stopifnot(length(added_csvs) > 0)

read_to_mat <- function(f) {
  df <- read_csv(f)
  df <- rename(df, cell_id = X1)

  condition <- strsplit(f, "assignments-")[[1]]
  condition <- condition[length(condition)]
  condition <- gsub(".csv", "", condition, fixed=TRUE)
 
  df <- mutate(df, condition = condition) %>% 
    select(cell_id, condition, everything())
  
  return(df)
}
 
dfs <- lapply(added_csvs, read_to_mat)

conditions <- sapply(dfs, function(d) d$condition[1])
names(dfs) <- conditions

get_prop_assigned <- function(cond, threshold, dfs) {
  df <- dfs[[ cond ]]
  df2 <- filter(df, condition == cond) 
  
  mat <- df2 %>% 
    select(-cell_id, -condition) %>% 
    as.matrix()
  
  mat <- mat[, !is.na(colSums(mat))]
  
  cell_assignments <- taproom::get_celltypes(mat, threshold)
  tibble(cohort = args$cohort, condition = cond, threshold = threshold,
         prop_assigned = mean(cell_assignments == cond))
}


thresholds <- c(0.5, 0.7, 0.9, 0.99, 0.999)

dft <- lapply(thresholds, function(th) {
  map_dfr(conditions, get_prop_assigned, th, dfs)
}) %>% 
  bind_rows()


write_tsv(dft, args$output_file)


