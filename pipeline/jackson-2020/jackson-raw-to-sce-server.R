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


write_one_sce <- function(core_test, df_sc, df_loc, output_dir) {
  message(paste("Writing core", core_test))
  df_loc <- filter(df_loc, core == core_test)
  
  df_count <- df_sc %>% 
    filter(core == core_test) %>% 
    select(id, channel, mc_counts) %>% 
    spread(channel, mc_counts)
  
  df_count_loc <- inner_join(df_loc, df_count, by = "id")
  
  
  ## Parse to an SCE
  ids <- df_count_loc$id
  
  # bad heuristic -- antibodies have space in name
  df_count_only <- select_at(df_count_loc, vars(contains(" "))) %>% 
    as.matrix()
  df_pd <- select_at(df_count_loc, vars(-contains(" "))) %>% 
    as.data.frame()
  
  ## Cell names
  rownames(df_pd) <- rownames(df_count_only) <- ids
  
  sce <- SingleCellExperiment(
    assays = list(raw_imc = t(df_count_only)),
    colData = df_pd
  )

  sce <- tidy_rownames_jackson(sce)

  #logcounts(sce) <- log( assay(sce, 'raw_imc') + 1)
  #sce <- winsorize(sce, w_limits = c(0.01, 0.99))
  logcounts(sce) <- asinh( assay(sce, 'raw_imc') / 5 )
  sce <- winsorize(sce, w_limits = c(0, 0.999))

  ## Write SingleCellExperiment as rds
  saveRDS(sce, file.path(output_dir, paste0(core_test, ".rds")))

  ## Write expression as csv
  to_csv(sce, file.path(output_dir, paste0(core_test, ".csv")))
  
}

p <- arg_parser("Read Jackson 2020 data")

p <- add_argument(p, "--input_sc", "Input SC_dat.csv")
p <- add_argument(p, "--input_loc", "Input location csv")
p <- add_argument(p, "--output_dir", "Output directory")

argv <- parse_args(p)

df_sc <- read_csv(argv$input_sc)
df_loc <- read_csv(argv$input_loc)

cores <- sort(unique(df_sc$core))

for(core in cores) {
  write_one_sce(core, df_sc, df_loc, argv$output_dir)
}





