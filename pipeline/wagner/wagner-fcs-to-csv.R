
suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(scater)
  library(argparser)
  library(devtools)
  library(flowCore)
  library(stringr)
})

select <- dplyr::select
mutate <- dplyr::mutate
arrange <- dplyr::arrange
rename <- dplyr::rename
filter <- dplyr::filter

devtools::load_all("../taproom")

# fields <- c("Center", "Event_length", "Offset", "Residual", "Time", "Width")
# fields <- c("DNA1", "DNA2")


p <- arg_parser("Read Wagner 2019 data")

p <- add_argument(p, "--input_fcs", "Input fcs file")

# p <- add_argument(p, "--cluster_identities", "Cluster identity TSV")
p <- add_argument(p, "--output_rds", "Output rds file")
p <- add_argument(p, "--output_csv", "Output csv file")

argv <- parse_args(p)

fcs <- read.FCS(argv$input_fcs)


clusters <- fcs@exprs[,'cluster']

# unique_id <- apply(fcs@exprs[,fields], 1, paste, collapse="-")


rn <- colnames(exprs(fcs))
is_protein <- grepl("_", rn, fixed = TRUE)
fcs <- fcs[, is_protein]
rn_new <- sapply(strsplit(colnames(fcs), "_", fixed = TRUE), `[`, 2)
colnames(fcs) <- rn_new


fcs <- fcs[, c(13,14:53)]

sce <- SingleCellExperiment(
  assays = list('raw_imc' = t(exprs(fcs)))
)

# logcounts(sce) <- log( assay(sce, 'raw_mc') + 1)
# sce <- winsorize(sce, w_limits = c(0.01, 0.99))
logcounts(sce) <- asinh( assay(sce, 'raw_imc') / 5 )
# sce <- winsorize(sce, w_limits = c(0, 0.999))

## Add in information

guid <- fcs@description$GUID
colData(sce)$guid <- gsub(".fcs", "", guid, fixed = TRUE)

colnames(sce) <- paste0(sce$guid, "_cell_", seq_along(sce$guid))

# colData(sce)$unique_id <- unique_id
colData(sce)$cluster <- clusters

bb_pos <- str_locate(guid, "BB")[1, 'start']
patient_id <- str_sub(guid, bb_pos, bb_pos + 4) # patient ID is 5 chars long
colData(sce)$patient_id <- patient_id

plate_pos <- str_locate(guid, "Plate")[1,'start']
plate <- str_sub(guid, plate_pos, plate_pos + 5) # plate pos is 6 chars long
colData(sce)$plate <- plate



# sce <- calculateQCMetrics(sce, exprs_values = "logcounts")
# sce <- sce[, sce$total_logcounts > 55]


## Write SingleCellExperiment as rds
saveRDS(sce, argv$output_rds)

## Write expression as csv
to_csv(sce, argv$output_csv, include_xy = FALSE) # No spatial data
  



