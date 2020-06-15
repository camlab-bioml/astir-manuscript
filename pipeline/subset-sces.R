## Takes a directory of RDS files containing SingleCellExperiments
## and a csv file whose first name is cell IDs, and joins and subsets
## the SingleCellExperiments to return a single one as output

suppressPackageStartupMessages({
  library(scater)
  library(tidyverse)
  library(argparse)
  library(SingleCellExperiment)
})


subset_sce <- function(rds_file, cell_ids) {
  sce <- readRDS(rds_file)
  if(any(cell_ids %in% colnames(sce))) {
    cell_ids <- intersect(cell_ids, colnames(sce))
    rowData(sce) <- NULL
    return( sce[, cell_ids] )
  } else {
    return(NULL)
  }
}

parser <- ArgumentParser(description = "Convert dir of rds to SCE")

parser$add_argument('--input_dir', type = 'character',
                    help='Directory of RDS containing SCE')
parser$add_argument('--input_csv', type = 'character',
                    help='Path to CSV whose first column is cell ID')
parser$add_argument('--output_rds', type = 'character',
                    help='Output unified SCE')

args <- parser$parse_args()

df <- read_csv(args$input_csv)
cell_ids <- df[[1]]

rds_files <- dir(args$input_dir, pattern="rds", full.names = TRUE)

sces <- lapply(rds_files, subset_sce, cell_ids)
non_null_names <- which(!sapply(sces, is.null))
sces <- sces[non_null_names]

sce <- do.call("cbind", sces)

saveRDS(sce, args$output_rds)


