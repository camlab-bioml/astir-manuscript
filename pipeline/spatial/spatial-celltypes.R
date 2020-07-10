
suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(matrixStats)
  library(argparse)
})

devtools::load_all("../taproom/")


parser <- ArgumentParser(description = "Spatial locations for core")

parser$add_argument('--input_locations', type = 'character',
                    help='Input location csv')
parser$add_argument('--input_assignments', type = 'character',
                    help='Input assignment csv')
parser$add_argument('--marker_yaml', type='character')
parser$add_argument('--core', type = 'character',
                    help='Core')
parser$add_argument('--output_csv', type = 'character',
                    help='Output csv file')


args <- parser$parse_args()

locations <- read_csv(args$input_locations)
assignments <- read_csv(args$input_assignments)
cell_names <- assignments$X1
assignments$X1 <- NULL
df_cell_types <- tibble(
  id = cell_names,
  cell_type = get_celltypes(assignments)
)

markers <- read_yaml(args$marker_yaml)

df <- inner_join(locations, df_cell_types, by = "id")

df <- filter(df, core == args$core)

cell_types <- names(markers$cell_types)


get_celltype_dist <- function(df, ct1, ct2) {

  
  df1 <- filter(df, cell_type == ct1)
  df2 <- filter(df, cell_type == ct2)
  x1 <- select(df1, Location_Center_X, Location_Center_Y) %>% as.matrix()
  x2 <- select(df2, Location_Center_X, Location_Center_Y) %>% as.matrix()
  if(nrow(x1) == 0 || nrow(x2) == 0) {
    return(NA)
  }

  
  results <- matrix(NA, nrow = nrow(x1), ncol = nrow(x2))
  
  for(i in 1:nrow(x1)) {
    for(j in 1:nrow(x2)) {
      results[i,j] <- sqrt(sum( (x1[i,] - x2[j,])^2 ))
    }
  }
  mean(results)
}

ct_combs <- expand.grid(cell_types, cell_types)

is_dup <- apply(ct_combs, 1, function(x) {
  paste0(sort(x), collapse = "_")
}) %>% duplicated()
ct_combs <- ct_combs[!is_dup,]


dists <- sapply(seq_len(nrow(ct_combs)), function(i) {
  get_celltype_dist(df, ct_combs[i,1], ct_combs[i,2])
})

dist_df <- mutate(ct_combs, avg_distance = dists) %>% 
  rename(Celltype_1 = Var1, Celltype_2 = Var2) %>% 
  mutate(core = args$core,
         n_cells = nrow(df)) %>% 
  select(core, everything())

write_csv(dist_df, args$output_csv)




