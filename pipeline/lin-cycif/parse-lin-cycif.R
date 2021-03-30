.libPaths(c("/home/ltri/campbell/kcampbel/R/x86_64-redhat-linux-gnu-library/3.6", .libPaths()))

suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(argparser)
  library(devtools)
  library(matrixStats)
})

select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
rename <- dplyr::rename
filter <- dplyr::filter

devtools::load_all("../taproom")



p <- arg_parser("Read lin cycif data")

p <- add_argument(p, "--input_sc", "Input CSV")
p <- add_argument(p, "--id", "Sample id")
p <- add_argument(p, "--output_csv", "Output CSV location")
p <- add_argument(p, "--output_rds", "Output RDS location")

argv <- parse_args(p)

df <- read_csv(argv$input_sc)


cd <- select(df, Area:Y)
df <- select(df, -(Area:Y))
expr_mat_raw <- t(as.matrix(df))

expr_mat_raw <- t(apply(expr_mat_raw, 1, winsorize_one, c(0.01, 0.99)))

rm <- rowMins(expr_mat_raw)
expr_mat_raw <- expr_mat_raw - rm # make minimum as zero

expr_mat <- asinh(expr_mat_raw / 100)

# expr_mat <- t( scale( t (expr_mat ), center = FALSE))

rownames(expr_mat) <- rownames(expr_mat_raw)

sce <- SingleCellExperiment(
  assays = list(raw = expr_mat_raw, logcounts = expr_mat),
  colData = cd
)

colnames(sce) <- paste0(argv$id, "_", seq_len(ncol(sce)))

to_csv(sce, argv$output_csv, include_xy=FALSE)

saveRDS(sce, argv$output_rds)







