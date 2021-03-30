
suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(matrixStats)
  library(argparse)
  library(scater)
  library(SingleCellExperiment)
})

devtools::load_all("../taproom/")


parser <- ArgumentParser(description = "Spatial locations for core")

parser$add_argument('--input_dir', type = 'character',
                    help='Input expression dir')
parser$add_argument('--input_assignments', type = 'character',
                    help='Input assignment csv')
parser$add_argument('--output_tsv', type = 'character',
                    help='Output tsv file')
parser$add_argument('--output_cells_counted', type = 'character',
                    help='Output cells counted')

# args <- list(
#   input_dir = "output/chipmunk/basel_processed",
#   input_assignments = "output/chipmunk/astir_assignments/basel_astir_assignments.csv"
# )


args <- parser$parse_args()


# Get cell types ----------------------------------------------------------
assignments <- read_csv(args$input_assignments)
cell_names <- assignments$X1
assignments$X1 <- NULL
df_cell_types <- tibble(
  id = cell_names,
  cell_type = get_celltypes(assignments, thresh=0.5)
)


# Read in expression data -------------------------------------------------

rds_files <- dir(args$input_dir, pattern="*.rds", full.names=TRUE)

sces <- lapply(rds_files, readRDS)

sce <- do.call('cbind', sces)

df_expr <- as.data.frame(t(logcounts(sce))) %>% 
  rownames_to_column('id') %>% 
  mutate(core = sce$core)

df_expr <- inner_join(df_expr, df_cell_types, by = "id")

df_cell_pct <- select(df_expr, id, cell_type, core)

df_expr <- gather(df_expr, marker, expression, -cell_type, -id, -core)
df_expr <- as_tibble(df_expr)

df_expr2 <- group_by(df_expr, core, cell_type, marker) %>% 
  summarize(mean_expression = mean(expression))

write_tsv(df_expr2, args$output_tsv)

cells_counted <- dplyr::count(df_cell_pct, core, cell_type)

cells_counted <- group_by(cells_counted, core) %>% 
  mutate(frac_cells = n / sum(n))

cells_counted %>% 
  write_tsv(args$output_cells_counted)





