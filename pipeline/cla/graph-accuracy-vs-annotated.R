
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

print(snakemake@input)


devtools::load_all(snakemake@params[['taproom_path']])

acc_wrap <- function(tt) {
  # annotated <- unique(as.character(tt$cell_type_annotated))
  cell_types <- unique(intersect(tt$cell_type_annotated, tt$cell_type_predicted))
  
  tt$cell_type_annotated <- factor(tt$cell_type_annotated, levels = cell_types)
  tt$cell_type_predicted <- factor(tt$cell_type_predicted, levels = cell_types)

  bind_rows(
     tryCatch({kap(tt, cell_type_annotated, cell_type_predicted)}, error=function(e) NULL),
    tryCatch({bal_accuracy(tt, cell_type_annotated, cell_type_predicted)}, error=function(e) NULL),
     tryCatch({mcc(tt, cell_type_annotated, cell_type_predicted)}, error=function(e) NULL)
  )
}

cohort <- snakemake@params[['cohort']]




# Annotations -------------------------------------------------------------

# df_annot <- dir("annotation/annotations/basel/", pattern=paste0(cohort, "*.*annotation"),
#     full.names = TRUE) %>% 

print("Reading annotations:")
cat(snakemake@input[['annotations']])

df_annot <- snakemake@input[['annotations']] %>% 
  map_dfr(read_tsv) %>% 
  rename(cell_type_annotated = cell_type,
         annotator_test = annotator)

df_train_test <- read_tsv(snakemake@input[['traintest']])

df_train_test <- df_train_test[df_train_test$cell_id %in% df_annot$cell_id,]

cell_ids_test <- filter(df_train_test, train_test == "test") %>% .$cell_id

df_annot <- filter(df_annot, cell_id %in% cell_ids_test)

cell_types <- filter(df_annot, cell_id %in% cell_ids_test) %>% 
  .$cell_type_annotated %>% unique()

print("Finished reading annotations")




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


astir_types_default <- get_astir_assignments(df_astir, 0.5, cell_ids_test)


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

astir_types_high_confidence <- get_astir_assignments(df_astir, 0.95, cell_ids_test)

df_astir_high_confidence <-   inner_join(df_annot, astir_types_high_confidence) %>% 
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
  mutate(method = "Astir high confidence",
         annotator_train = "None")



# CytofLDA ----------------------------------------------------------------

cat("\n Reading cytofLDA \n")


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


# ACDC ----------------------------------------------------------------

cat("\n Reading ACDC \n")


types_acdc <- dir(snakemake@params[['acdc_path']], 
                      pattern=paste0("annotations_acdc_*.*",cohort), full.names=TRUE) %>% 
  map_dfr(read_tsv) %>% 
  rename(annotator_train = annotator, cell_type_predicted = cell_type)

types_acdc$cell_type_predicted[types_acdc$cell_type_predicted == "unknown"] <- "Unclear"


df_acdc <- inner_join(df_annot, types_acdc) %>% 
  mutate(
    cell_type_annotated = factor(cell_type_annotated, levels=cell_types),
    cell_type_predicted = factor(cell_type_predicted, levels=cell_types)
  ) %>% 
  group_by(annotator_train, annotator_test, method) %>% 
  do(
    acc_wrap(.)
  ) %>% 
  ungroup() #%>% 
  # mutate(method = "ACDC")


# Alternative workflows ---------------------------------------------------

cat("\n Reading other \n")

df_other <- dir(snakemake@params[['other_workflow_path']],
    pattern=cohort,
    full.names=TRUE) %>% 
  map_dfr(read_csv)


df_other <- gather(df_other, annotation_method, cell_type_predicted, -(id:params))

df_other <- rename(df_other, cell_id = id)

df_other$method <- paste0(df_other$method, "_", df_other$annotation_method)

df_other <- select(df_other, -annotation_method)

df_other <- inner_join(df_annot, df_other)

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

# save.image("~/deleteme.rds  ")

# Overall plot ------------------------------------------------------------


# Overall plot ------------------------------------------------------------

df_plot <- bind_rows(
  df_astir_default,
  df_cytoflda,
  df_other_acc,
  df_acdc
)

df_plot$annotator_test <- paste("Test annotator: ", df_plot$annotator_test)

df_plot <- mutate(df_plot, .metric = case_when(
  .metric == "f_meas" ~ "F-measure",
  .metric == "kap" ~ "Cohen's\nkappa",
  .metric == "bal_accuracy" ~ "Balanced\naccuracy",
  .metric == "mcc" ~ "MCC",
  TRUE ~ stringr::str_to_title(.metric)
))

df_plot <- mutate(df_plot, .metric = case_when(
  .metric == "f_meas" ~ "F-measure",
  .metric == "kap" ~ "Cohen's\nkappa",
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


df_plot <- df_plot[!is.nan(df_plot$.estimate),]

df_plot <- group_by(df_plot, .metric) %>%
  mutate(norm_estimate = (.estimate - mean(.estimate, na.rm=T)) / sd(.estimate, na.rm=T) )

fill_cols <- c("None"="grey50",
               "Annotator-1"=scales::muted('blue'),
               "Annotator-2"=scales::muted('red'))

ggplot(df_plot, aes(x = forcats::fct_reorder(method, norm_estimate), y = .estimate, fill = method_type)) +
  geom_bar(stat='identity', position = "dodge2") +
  geom_boxplot() +
  facet_grid(.metric ~ annotator_test, scales="free_y") +
  astir_paper_theme() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "Method", y = "Estimate") +
  scale_fill_manual(values=method_cols) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(size=8))

df_plot$cohort <- snakemake@params[['cohort']]

ggsave(snakemake@output[['plot']], width=10, height=10)

write_tsv(df_plot, snakemake@output[['tsv']])





