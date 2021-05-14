
library(tidyverse)
library(yardstick)
library(yaml)
library(matrixStats)


source("pipeline/cla/helpers.R")

#' We need:
#' (1) train test 
#' (2) cohort name
#' (3) directory with annotation tsvs
#' (4) astir assignment csv
#' (5) directory to cytofLDA results
#' (6) directory to other workflow results
#' (7) output plot
#' (8) output tsv of accuracies
#' (9) coarse-fine mapping file


cohort <- snakemake@params[['cohort']]

devtools::load_all(snakemake@params[['taproom_path']])


coarse_mapping <- read_yaml(snakemake@input[['coarse_fine_mapping']])
coarse_mapping <- unlist(Biobase::reverseSplit(coarse_mapping))

coarse_mapping['Unclear'] <- 'Unclear'

remap_to_coarse <- function(dfr, mapping = coarse_mapping, input_column = 'cell_type_predicted', output_column = 'cell_type_predicted') {
  dfr[[output_column]] <- plyr::mapvalues(
    dfr[[input_column]],
    from=names(mapping),
    to=mapping
  )
  dfr
}

coarse_mapping['Unclear'] <- 'Unclear'


cell_types <- unique(coarse_mapping)


# Annotations -------------------------------------------------------------

cat("\n Reading previous annotations \n")

df_cluster <- read_csv(snakemake@input[['clustering']]) %>% 
  rename(cell_type_annotated=cell_type)

df_cluster <- filter(df_cluster, cell_type_annotated != "Apoptotic")

print(df_cluster)

df_train_test <- read_tsv(snakemake@input[['traintest']])

cell_ids_train <- filter(df_train_test, train_test == "train") %>% .$cell_id


df_cluster <- filter(df_cluster, !(cell_id %in% cell_ids_train))

# Astir -------------------------------------------------------------------

get_astir_assignments <- function(df_astir, threshold, cell_ids) {
  df <- filter(df_astir, X1 %in% cell_ids)
  astir_cell_types <- select(df, -X1) %>%
    taproom::get_celltypes(thresh=threshold)
  astir_cell_types[astir_cell_types == "Other"] <- "Unclear"
  astir_cell_types[astir_cell_types == "Unknown"] <- "Unclear"
  df_types <- tibble(
    cell_id = df$X1,
    cell_type_predicted = astir_cell_types
  )
}

cat("\n Reading ASTIR \n")

df_astir <- read_csv(snakemake@input[['astir_assignments']])
cell_ids_test <- setdiff(df_astir$X1, cell_ids_train)

astir_types_default <- get_astir_assignments(df_astir, 0.5, cell_ids_test)
astir_types_high_confidence <- get_astir_assignments(df_astir, 0.95, cell_ids_test)

astir_types_default <- remap_to_coarse(astir_types_default)
astir_types_high_confidence <- remap_to_coarse(astir_types_high_confidence)

df_astir_default <- inner_join(df_cluster, astir_types_default) %>% 
  mutate(
    cell_type_annotated = factor(cell_type_annotated, levels=cell_types),
    cell_type_predicted = factor(cell_type_predicted, levels=cell_types)
  ) %>% 
  do(
    acc_wrap(.)
  ) %>% 
  ungroup() %>% 
  mutate(method = "Astir",
         annotator_train = "None")

df_astir_high_confidence <- inner_join(df_cluster, astir_types_high_confidence) %>% 
 filter(cell_type_predicted != "Unknown") %>% 
  mutate(
    cell_type_annotated = factor(cell_type_annotated, levels=cell_types),
    cell_type_predicted = factor(cell_type_predicted, levels=cell_types)
  ) %>% 
  do(
    acc_wrap(.)
  ) %>% 
  ungroup() %>% 
  mutate(method = "Astir high confidence",
         annotator_train = "None")



# CytofLDA ----------------------------------------------------------------

cat("\n Reading cytofLDA \n")

types_cytoflda <- dir(snakemake@params[['cytofLDA_path']], pattern=paste0("annotations_cytofLDA_",cohort), full.names=TRUE) %>% 
  map_dfr(read_tsv) %>% 
  rename(annotator_train = annotator, cell_type_predicted = cell_type)

types_cytoflda <- remap_to_coarse(types_cytoflda)


df_cytoflda <- inner_join(df_cluster, types_cytoflda) %>% 
  # filter(cell_type_annotated != "Unclear",
  #        cell_type_predicted != "Unclear") %>% 
  mutate(
    cell_type_annotated = factor(cell_type_annotated, levels=cell_types),
    cell_type_predicted = factor(cell_type_predicted, levels=cell_types)
  ) %>% 
  group_by(annotator_train) %>% 
  do(
    acc_wrap(.)
  ) %>% 
  ungroup() %>% 
  mutate(method = "CytofLDA")


# ACDC --------------------------------------------------------------------

cat("\n Reading ACDC \n")


types_acdc <- dir(snakemake@params[['acdc_path']], 
                      pattern=paste0("annotations_acdc_*.*",cohort), full.names=TRUE) %>% 
  map_dfr(read_tsv) %>% 
  rename(annotator_train = annotator, cell_type_predicted = cell_type)

types_acdc$cell_type_predicted[types_acdc$cell_type_predicted == "unknown"] <- "Unclear"

types_acdc <- remap_to_coarse(types_acdc)


df_acdc <- inner_join(df_cluster, types_acdc) %>% 
  mutate(
    cell_type_annotated = factor(cell_type_annotated, levels=cell_types),
    cell_type_predicted = factor(cell_type_predicted, levels=cell_types)
  ) %>% 
  group_by(method) %>% 
  do(
    acc_wrap(.)
  ) %>% 
  ungroup() %>%
  mutate(annotator_train = "None")


# Alternative workflows ---------------------------------------------------

cat("\n Reading Other \n")

df_other <- dir(snakemake@params[['other_workflow_path']],
    pattern=cohort,
    full.names=TRUE) %>% 
  map_dfr(read_csv)

print(head(df_other))

df_other <- gather(df_other, annotation_method, cell_type_predicted, -(id:params))

df_other <- rename(df_other, cell_id = id)

df_other$method <- paste0(df_other$method, "_", df_other$annotation_method)

df_other <- select(df_other, -annotation_method)


df_other <- remap_to_coarse(df_other)

df_other <- inner_join(df_cluster, df_other, by="cell_id")

df_other_acc <- df_other %>% mutate(
  cell_type_annotated = factor(cell_type_annotated, levels=cell_types),
  cell_type_predicted = factor(cell_type_predicted, levels=cell_types)
) %>%
  group_by(method, params) %>%
  do(
    acc_wrap(.)
  ) %>%
  ungroup() %>%
  mutate(annotator_train = "None")



# Overall plot ------------------------------------------------------------


df_plot <- bind_rows(
  df_astir_default,
  # df_astir_high_confidence,
  df_cytoflda,
  df_other_acc,
  df_acdc
)

df_plot$annotator_test <- paste("Test annotator: ", df_plot$annotator_test)



df_plot <- mutate(df_plot, method = case_when(
  grepl("Astir", method) ~ "Astir",
  method == "acdc-absent" ~ "ACDC_absent",
  method == "acdc-no-consider" ~ "ACDC_no_consider",
  TRUE ~ method
))

df_plot <- df_plot[!is.nan(df_plot$.estimate),]

df_plot <- mutate(df_plot, .metric = case_when(
  .metric == "f_meas" ~ "F-measure",
  .metric == "kap" ~ "Cohen's\nkappa",
  .metric == "bal_accuracy" ~ "Balanced\naccuracy",
  .metric == "mcc" ~ "MCC",
  TRUE ~ stringr::str_to_title(.metric)
))

df_plot <- mutate(df_plot, method_type = case_when(
  grepl("LDA", method) ~ "Supervised",
  grepl("Astir|ACDC", method) ~ "Unsupervised",
  TRUE ~ "Cluster & interpret"
))

method_cols <- c(
  "Supervised"="#AA3939",
  "Unsupervised"="#882D61",
  "Cluster & interpret"="#AA6C39"
)

df_plot <- group_by(df_plot, .metric) %>%
  mutate(norm_estimate = (.estimate - mean(.estimate)) / sd(.estimate) )

ggplot(df_plot, aes(x = forcats::fct_reorder(method, norm_estimate), y = .estimate, fill = method_type)) +
  geom_bar(stat='identity', position = "dodge2") +
  geom_boxplot() +
  facet_grid(.metric ~ .)  +
  astir_paper_theme() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "Method", y = "Estimate") +
  scale_fill_manual(values=method_cols) +
  theme(legend.title = element_blank())

df_plot$cohort <- snakemake@params[['cohort']]

ggsave(snakemake@output[['plot']], width=10, height=10)

write_tsv(df_plot, snakemake@output[['tsv']])




