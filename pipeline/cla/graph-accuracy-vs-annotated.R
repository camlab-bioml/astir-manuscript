
suppressPackageStartupMessages({
  library(tidyverse)
  library(yardstick)
  library(matrixStats)
})

#' We need:
#' (1) train test directory
#' (2) cohort name
#' (3) directory with annotation tsvs
#' (4) astir assignment csv
#' (5) directory to cytofLDA results
#' (6) directory to other workflow results
#' (7) output plot
#' (8) output tsv of accuracies


devtools::load_all("taproom")

acc_wrap <- function(tt) {
  bind_rows(
  kap(tt, cell_type_annotated, cell_type_predicted),
  f_meas(tt, cell_type_annotated, cell_type_predicted),
  precision(tt, cell_type_annotated, cell_type_predicted),
  recall(tt, cell_type_annotated, cell_type_predicted)
  )
}

cohort <- snakemake@params[['cohort']]


df_train_test <- read_tsv(snakemake@input[['traintest']])

cell_ids_test <- filter(df_train_test, train_test == "test") %>% .$cell_id

# Annotations -------------------------------------------------------------

df_annot <- snakemake@input[['annotations']] %>%
  map_dfr(read_tsv)

df_annot <- bind_rows(df_annot1, df_annot2) %>% 
  rename(cell_type_annotated = cell_type)

cell_types <- unique(df_annot$cell_type_annotated)

df_annot <- filter(df_annot, cell_type_annotated != "Unclear")


# Astir -------------------------------------------------------------------

get_astir_assignments <- function(df_astir, threshold, cell_ids) {
  df <- filter(df_astir, X1 %in% cell_ids)
  astir_cell_types <- select(df, -X1) %>%
    taproom::get_celltypes(thresh=threshold)
  df_types <- tibble(
    cell_id = df$X1,
    cell_type_predicted = astir_cell_types
  )
}


df_astir <- read_csv(snakemake@input[['astir_assignments']])

astir_types_default <- get_astir_assignments(df_astir, 0.5, cell_ids_test)
astir_types_high_confidence <- get_astir_assignments(df_astir, 0.95, cell_ids_test)

df_astir_default <-   inner_join(df_annot, astir_types_default) %>% 
  # filter(cell_type_annotated != "Unclear",
  #        cell_type_predicted != "Unclear") %>% 
  mutate(
    cell_type_annotated = factor(cell_type_annotated, levels=cell_types),
    cell_type_predicted = factor(cell_type_predicted, levels=cell_types)
  ) %>% 
  group_by(annotator_test) %>% 
  do(
    acc_wrap(.)
  ) %>% 
  ungroup() %>% 
  mutate(method = "Astir default",
         annotator_train = "None")

df_astir_high_confidence <- inner_join(df_annot, astir_types_high_confidence) %>% 
 filter(cell_type_predicted != "Unknown") %>% 
  mutate(
    cell_type_annotated = factor(cell_type_annotated, levels=cell_types),
    cell_type_predicted = factor(cell_type_predicted, levels=cell_types)
  ) %>% 
  group_by(annotator_test) %>% 
  do(
    acc_wrap(.)
  ) %>% 
  ungroup() %>% 
  mutate(method = "Astir high confidence",
         annotator_train = "None")



# CytofLDA ----------------------------------------------------------------

types_cytoflda <- dir(snakemake@params[['cytofLDA_path']], pattern=paste0("annotations_cytofLDA_",cohort), full.names=TRUE) %>% 
  map_dfr(read_tsv) %>% 
  rename(annotator_train = annotator, cell_type_predicted = cell_type)


df_cytoflda <- inner_join(df_annot, types_cytoflda) %>% 
  # filter(cell_type_annotated != "Unclear",
  #        cell_type_predicted != "Unclear") %>% 
  mutate(
    cell_type_annotated = factor(cell_type_annotated, levels=cell_types),
    cell_type_predicted = factor(cell_type_predicted, levels=cell_types)
  ) %>% 
  group_by(annotator_train, annotator_test) %>% 
  do(
    acc_wrap(.)
  ) %>% 
  ungroup() %>% 
  mutate(method = "CytofLDA")


# Alternative workflows ---------------------------------------------------

df_other <- dir(snakemake@params['other_workflow_path'],
    pattern=cohort,
    full.names=TRUE) %>% 
  map_dfr(read_csv)

## REMOVE ME
# df_other$cell_type <- sample(df_other$cell_type)

df_other <- mutate(df_other, 
                   cell_type_predicted = cell_type,
                   cell_id = id)

df_other <- inner_join(df_annot, df_other)

# df_other <- mutate(df_other, method_long = paste0(method, "-", params))

df_other_acc <- df_other %>% mutate(
  cell_type_annotated = factor(cell_type_annotated, levels=cell_types),
  cell_type_predicted = factor(cell_type_predicted, levels=cell_types)
) %>% 
  group_by(annotator_test, method, params) %>% 
  do(
    acc_wrap(.)
  ) %>% 
  ungroup() %>% 
  mutate(annotator_train = "None")



# Overall plot ------------------------------------------------------------

df_plot <- bind_rows(
  df_astir_default,
  df_astir_high_confidence,
  df_cytoflda,
  df_other_acc
)

df_plot$annotator_test <- paste("Test annotator: ", df_plot$annotator_test)

ggplot(df_plot, aes(x = method, y = .estimate, fill=annotator_train)) +
  geom_bar(stat='identity', position = "dodge2") +
  geom_boxplot() +
  facet_grid(.metric ~ annotator_test)

ggsave(snakemake@output[['plot']])

write_tsv(df_plot, snakemake@output[['tsv']])





