
library(scater)
library(SingleCellExperiment)
library(argparse)
library(tidyverse)

parser <- ArgumentParser(description = "Run cytofLDA")

parser$add_argument('--input_h5ad', type = 'character',
                    help='Input h5ad')
parser$add_argument('--input_traintest', type='character',
                    help='Train test splits')
parser$add_argument('--cytofLDA_path', type='character',
                    help='Path to cytofLDA implementation')
parser$add_argument('--input_labels', type='character',
                    help='Path to human labeled tsv')
parser$add_argument('--output_assignments', type = 'character',
                    help='Path to output assignments')
parser$add_argument('--annotator', type='character',
                    help='ID of annotator')
parser$add_argument('--cohort', type='character',
                    help='ID of annotator')
parser$add_argument('--method', type='character',
                    help='method (cytofLDA)')

args <- parser$parse_args()



df_annot <- readr::read_tsv(args$input_labels) 


df_annot <- dplyr::rename(df_annot, annotated_cell_type = cell_type)


df_train_test <- readr::read_tsv(args$input_traintest)

cell_ids_train <- filter(df_train_test, train_test == "train") %>% .$cell_id
cell_ids_test <- filter(df_train_test, train_test == "test") %>% .$cell_id

sce <- zellkonverter::readH5AD(args$input_h5ad)
assay(sce, 'logcounts') <- assay(sce, 'X')


sce <- sce[, colnames(sce) %in% c(cell_ids_train, cell_ids_test)]

data <- as.data.frame(t(logcounts(sce)))
data$cell_id <- rownames(data)

data <- inner_join(data, df_annot, by = "cell_id")

data_train <- filter(data, cell_id %in% cell_ids_train)
data_test <- filter(data, cell_id %in% cell_ids_test)

stopifnot(length(intersect(data_train$cell_id, data_test$cell_id)) == 0)

train_dir <- file.path(tempdir(), "train")
test_dir <- file.path(tempdir(), "test")

dir.create(train_dir)
dir.create(test_dir)
write.table(data_train,file = file.path(train_dir, "train.csv"),col.names = FALSE,row.names = FALSE,sep = ',')
write.table(data_test,file = file.path(test_dir, "test.csv"),col.names = FALSE,row.names = FALSE,sep = ',')

source(file.path(args$cytofLDA_path, 'CyTOF_LDAtrain.R'))
source(file.path(args$cytofLDA_path, 'CyTOF_LDApredict.R'))

LDA.Model <- CyTOF_LDAtrain(TrainingSamplesExt = train_dir,TrainingLabelsExt = '',mode = 'CSV',
                            RelevantMarkers =  c(3:40),
                            LabelIndex = 42, 
                            Transformation = FALSE)

Predictions <- CyTOF_LDApredict(LDA.Model,TestingSamplesExt = test_dir, mode = 'CSV', RejectionThreshold = 0)

Predictions <- unlist(Predictions)

# print(length(Predictions))
# print(dim(data_test))

df_output <- tibble(
  cell_id = data_test$cell_id,
  cell_type = Predictions,
  annotator = args$annotator,
  cohort = args$cohort,
  method = args$method
)

write_tsv(df_output, args$output_assignments)
