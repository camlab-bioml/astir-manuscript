
library(tidyverse)

df_cluster <- read_csv(snakemake@input[['clusters']])
df_annotation <- read_delim(snakemake@input[['annotation']], ";")

names(df_annotation) <- c('cluster', 'cell_type', 'class')
names(df_cluster) <- c('cell_id', 'cluster')

df <- inner_join(df_cluster, df_annotation)

write_csv(df, snakemake@output[[1]])




