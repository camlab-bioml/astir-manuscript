
suppressPackageStartupMessages({
  library(tidyverse)
  library(argparse)
  library(SingleCellExperiment)
})

devtools::load_all("../taproom")

parser <- ArgumentParser(description = "Read schapiro 2017")

parser$add_argument('--input_csv', type = 'character')
parser$add_argument('--sample', type = 'character')
parser$add_argument('--output_sce', type = 'character')
parser$add_argument('--output_csv', type = 'character')


args <- parser$parse_args()

df <- read_csv(args$input_csv)
df_spread <- select(df, -sample) %>% 
  spread(channel, expression)

cell_names <- df_spread$cell

expr_mat <- select(df_spread, -cell, -x, -y, -size) %>% 
  as.matrix() %>% 
  t()

colnames(expr_mat) <- cell_names

sce <- SingleCellExperiment(
  assays = list('raw_imc' = expr_mat)
)


colData(sce)$x <- df_spread$x
colData(sce)$y <- df_spread$y
colData(sce)$size <- df_spread$size

logcounts(sce) <- asinh( assay(sce, 'raw_imc') / 5 )
sce <- winsorize(sce, w_limits = c(0, 0.999))

## Get rid ofthe isotopes
rn <- sapply(strsplit(rownames(sce), "(", fixed=T), `[`, 1)
rownames(sce) <- rn

saveRDS(sce, args$output_sce)

to_csv(sce, args$output_csv, include_xy = FALSE)