.libPaths(c("/home/ltri/campbell/kcampbel/R/x86_64-redhat-linux-gnu-library/3.6", .libPaths()))

suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(argparser)
  library(devtools)
})

select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
rename <- dplyr::rename
filter <- dplyr::filter

devtools::load_all("../taproom")

p <- arg_parser("Parse MIBI")

p <- add_argument(p, "--input_sc", "Input SC_dat.csv")
p <- add_argument(p, "--output_rds", "Output SCE")
p <- add_argument(p, "--output_csv", "Output CSV")

argv <- parse_args(p)


df_sc <- read_csv(argv$input_sc)


df_spread <- select(df_sc, -X1) %>% 
  spread(channel, expression)

cd <- select(df_spread, cell, sample, x, y, size) %>% 
  as.data.frame()
rownames(cd) <- cd$cell

expression <- select(df_spread, - (sample:size) ) %>% 
  as.matrix()
rownames(expression) <- cd$cell

sce <- SingleCellExperiment(
  assays = list(raw_mibi = t(expression)),
  colData = cd
)

logcounts(sce) <- asinh( assay(sce, 'raw_mibi') / 5 )
sce <- winsorize(sce, w_limits = c(0, 0.999))

saveRDS(sce, argv$output_rds)

## Write expression as csv
to_csv(sce, argv$output_csv, include_xy = FALSE)
