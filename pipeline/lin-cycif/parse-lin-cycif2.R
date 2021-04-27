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

df <- read_tsv(argv$input_sc)


cd <- select(df, cell_id:sample)
df <- select(df, -(cell_id:sample))
expr_mat <- t(as.matrix(df))

feature_names <- sapply(rownames(expr_mat), function(s) strsplit(s, "-")[[1]][1])

feature_names <- unlist(feature_names)
names(feature_names) <- NULL

print(feature_names)

rownames(expr_mat) <- feature_names

# expr_mat_raw <- t(apply(expr_mat_raw, 1, function(x) {
#   x / mean(x)
# }))

# rm <- rowMins(expr_mat_raw)
# expr_mat_raw <- expr_mat_raw - rm # make minimum as zero

# expr_mat <- asinh(expr_mat_raw / 5)

# expr_mat <- t( scale( t (expr_mat ), center = FALSE))

# rownames(expr_mat) <- rownames(expr_mat_raw)

colnames(expr_mat) <- paste0(cd$sample, "_", cd$cell_id)

sce <- SingleCellExperiment(
  assays = list(raw = expr_mat, logcounts = expr_mat),
  colData = cd
)

sce <- sce[,logcounts(sce)['Hoechst/DNA4',] > 3]

# colnames(sce) <- paste0(argv$id, "_", seq_len(ncol(sce)))

to_csv(sce, argv$output_csv, include_xy=FALSE)

saveRDS(sce, argv$output_rds)







