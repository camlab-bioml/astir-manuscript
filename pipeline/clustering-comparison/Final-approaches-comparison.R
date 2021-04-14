library(argparse)
library(cowplot)
library(ggpubr)
library(ComplexHeatmap)
library(GSVA)
library(tidyverse)
library(taproom)
library(DelayedArray)
library(SingleCellExperiment)
source("scripts/functions.R")

parser <- ArgumentParser()
parser$add_argument('--basel_astir_assignments', type = 'character')
parser$add_argument('--schapiro_astir_assignments', type = 'character')
parser$add_argument('--wagner_astir_assignments', type = 'character')
parser$add_argument('--zurich_astir_assignments', type = 'character')
parser$add_argument('--lin_astir_assignments', type = 'character')
parser$add_argument('--basel_sce', type = 'character')
parser$add_argument('--schapiro_sce', type = 'character')
parser$add_argument('--wagner_sce', type = 'character')
parser$add_argument('--zurich_sce', type = 'character')
parser$add_argument('--lin_sce', type = 'character')
parser$add_argument('--jackson_markers', type = 'character')
parser$add_argument('--schapiro_markers', type = 'character')
parser$add_argument('--wagner_markers', type = 'character')
parser$add_argument('--lin_markers', type = 'character')
parser$add_argument('--basel_files', type = 'character', nargs = '+')
parser$add_argument('--schapiro_files', type = 'character', nargs = '+')
parser$add_argument('--wagner_files', type = 'character', nargs = '+')
parser$add_argument('--lin_files', type = 'character', nargs = '+')
parser$add_argument('--zurich_files', type = 'character', nargs = '+')
parser$add_argument('--output_heatmap', type = 'character', nargs = '+')

args <- parser$parse_args()


# Define basic functions
read_astir_assignments <- function(path){
  read_csv(path) %>% 
  column_to_rownames("X1") %>% 
  mutate(cell_type = get_celltypes(.)) %>% 
  rownames_to_column("id") %>% 
  select(id, cell_type) %>% 
  dplyr::rename("cluster" = "cell_type")
}

create_expression_mat <- function(sce, astir_assignment){
  expression <- logcounts(sce) %>% 
    t() %>% 
    as_tibble() %>% 
     mutate(id = colnames(sce))
  
  expression <- expression %>% 
    left_join(select(astir_assignment, id, cluster)) %>% 
    select(-id)
}

# Read in Astir assignments
basel_astir <- read_astir_assignments(args$basel_astir_assignments)
schapiro_astir <- read_astir_assignments(args$schapiro_astir_assignments)
wagner_astir <- read_astir_assignments(args$wagner_astir_assignments)
zurich_astir <- read_astir_assignments(args$zurich_astir_assignments)
lin_astir <- read_astir_assignments(args$lin_astir_assignments)


# Read in sces
basel_sce <- readRDS(args$basel_sce)
schapiro_sce <- readRDS(args$schapiro_sce)
wagner_sce <- readRDS(args$wagner_sce)
zurich_sce <- readRDS(args$zurich_sce)
lin_sce <- readRDS(args$lin_sce)


# Create expression matrices
basel_expression <- create_expression_mat(basel_sce, basel_astir)
schapiro_expression <- create_expression_mat(schapiro_sce, schapiro_astir)
wagner_expression <- create_expression_mat(wagner_sce, wagner_astir)
zurich_expression <- create_expression_mat(zurich_sce, zurich_astir)
lin_expression <- create_expression_mat(lin_sce, lin_astir)

# read in markers
basel_markers <- zurich_markers <- read_markers(args$jackson_markers)
schapiro_markers <- read_markers(args$schapiro_markers)
wagner_markers <- read_markers(args$wagner_markers)
lin_markers <- read_markers(args$lin_markers)

# Cluster assignments
create_counts <- function(expression, markers, cohort){
  GSVA_Astir_celltypes <- expression %>% 
    filter(!cluster %in% c("Other", "Unknown")) %>% 
    assign_clusters(markers)
  
  astir_counts <- GSVA_Astir_celltypes$assignment %>% 
    group_by(cell_type) %>% 
    tally() %>% 
    mutate(cohort = cohort, method = "Astir")
  astir_counts <- astir_counts[,c("cohort", "method", "cell_type", "n")]
}

basel_counts <- create_counts(basel_expression, basel_markers, "Basel")
schapiro_counts <- create_counts(schapiro_expression, schapiro_markers, "Schapiro")
wagner_counts <- create_counts(wagner_expression, wagner_markers, "Wagner")
zurich_counts <- create_counts(zurich_expression, zurich_markers, "Zurich")
lin_counts <- create_counts(lin_expression, lin_markers, "Lin")



# Read in other methods assignments
read_in_cohort <- function(files, cohort){
    data <- lapply(files, read_csv) %>% 
                bind_rows() %>% 
                mutate(cohort = cohort)
    
    data <- data %>%
        select(-c(id, Manual_cell_type)) %>% 
        distinct(.keep_all = TRUE)

    data <- data %>%
        mutate(params = str_replace(params, "_", " ")) %>% 
        mutate(params = str_replace(params, "_", "-")) %>% 
        mutate(method = paste0(method, "-", params)) %>% 
        select(-params) %>%
        dplyr::rename("cell_type" = "GSVA_cell_type")
    
    data
}
# basel <- lapply(args$basel_files, read_csv) %>% 
#   bind_rows() %>% 
#   mutate(cohort = "Basel")
basel <- read_in_cohort(args$basel_files, "Basel")
head(basel)
# schapiro <- lapply(args$schapiro_files, read_csv) %>% 
#   bind_rows() %>% 
#   mutate(cohort = "Schapiro")
schapiro <- read_in_cohort(args$schapiro_files, "Schapiro")
head(schapiro)
# wagner <- lapply(args$wagner_files, read_csv) %>% 
#   bind_rows() %>% 
#   mutate(cohort = "Wagner")
wagner <- read_in_cohort(args$wagner_files, "Wagner")
head(wagner)
# zurich <- lapply(args$zurich_files, read_csv) %>% 
#   bind_rows() %>% 
#   mutate(cohort = "Zurich")
zurich <- read_in_cohort(args$zurich_files, "Zurich")
head(zurich)
# lin <- lapply(args$lin_files, read_csv) %>% 
#   bind_rows() %>% 
#   mutate(cohort = "Lin")
lin <- read_in_cohort(args$lin_files, "Lin")
head(lin)


# all_cohorts <- bind_rows(basel, schapiro, wagner, zurich, lin) %>% 
#   mutate(params = str_replace(params, "_", " ")) %>% 
#   mutate(params = str_replace(params, "_", "-")) %>% 
#   mutate(method = paste0(method, "-", params)) %>% 
#   select(-params)
all_cohorts <- bind_rows(basel, schapiro, wagner, zurich, lin)


# All counts for other methods
all_counts <- all_cohorts %>% 
#   select(-id) %>% 
#   distinct(.keep_all = TRUE) %>% 
  dplyr::group_by(cohort, method, cell_type) %>% 
  tally() %>% 
  ungroup()

all_counts$n[is.na(all_counts$cell_type)] <- NA

# Add other method counts and astir counts
all_counts <- bind_rows(all_counts, basel_counts, schapiro_counts, wagner_counts,
                        zurich_counts, lin_counts)



all_scores <- all_counts %>% 
  ungroup() %>%
  mutate(score = case_when(
    n == 1 ~ 1,
    n != 1 ~ -1
  )) %>% group_by(cohort, method) %>% 
  dplyr::summarize(score = sum(score)) %>% 
  ungroup()

# Get maximum clusters
max.clusters <- all_counts %>% 
  filter(!is.na(n)) %>% 
  pull(n) %>% 
  max()

if(max.clusters > 10){
  max.clusters <- 10
}

# Create colour palette
gradient <- colorRampPalette(c("#FFEBEE", "#F50025"))
pal <- c("#FFDA1F", "#00D262", 
         gradient(max.clusters - 1))
names(pal) <- c(0:(max.clusters - 1), "10+")

tooHigh <- as.character(c(0:9))



all_scores_wide <- all_scores %>% 
  pivot_wider(values_from = "score", names_from = "cohort") %>% 
  column_to_rownames("method") %>% 
  mutate(mean_score = rowMeans(., na.rm = TRUE))

plottingOrder = all_scores_wide %>% 
    arrange(desc(mean_score)) %>%
    rownames()
  
plot_eval_heatmap <- function(counts_df, scores_df, plottingOrder, select_cohort){
  cohort_counts <- filter(all_counts, cohort == "Basel")
  cohort_scores <- filter(all_scores, cohort == "Basel") %>% 
    select(-cohort)
  # filter dataframes to select required data
  cohort_counts <- filter(counts_df, cohort == select_cohort)
  cohort_scores <- filter(scores_df, cohort == select_cohort) %>% 
    select(-cohort)
  
  # Create matrices to plot
  # pivot wider such that cell types can be column names
  counts_wide <- cohort_counts %>% 
    pivot_wider(names_from = "cell_type", values_from = "n", values_fill = 0)
  
  if("NA" %in% colnames(counts_wide)){
    # Get the number of cell types for this cohort
    no_of_cell_types <- na.omit(cohort_counts$cell_type) %>% 
      unique() %>% 
      length()
    
    # Set all cluster numbers to NA for those methods where only one cluster was identified
    counts_wide[is.na(counts_wide$`NA`), 3:(3 + no_of_cell_types)] <- NA
    
    cohort_mat <- counts_wide %>% 
    select(-c(`NA`, cohort)) %>% 
    column_to_rownames("method") %>% 
    as.matrix()
  }else{
    cohort_mat <- counts_wide %>% 
      select(-cohort) %>% 
      column_to_rownames("method") %>% 
      as.matrix()
  }
  
  # Replace anything above 10 with 10+
  cohort_mat[!(cohort_mat %in% tooHigh) & !is.na(cohort_mat)] <- "10+"

  barplot_score <- data.frame(cohort_scores[match(rownames(cohort_mat), cohort_scores$method),])

  # Create annotation
  cohort.ha = rowAnnotation(Score = anno_barplot(barplot_score$score, 
                                                gp = gpar(fill = "lightgrey")),
                     width = unit(2, "cm"))
  
  # Create heatmap
  cohort.hm <- Heatmap(cohort_mat, col = pal,
                      row_order = plottingOrder,
                      right_annotation = cohort.ha,
                      column_title = select_cohort,
                      row_names_side = "left",
                      row_names_max_width = unit(7, "cm"),
                      name = "# clusters \nassigned to \ncell type")
  
  
  cohort.hm
}



basel.hm <- plot_eval_heatmap(counts_df = all_counts, 
                              scores_df = all_scores, 
                              plottingOrder = plottingOrder, 
                              select_cohort = "Basel")

schapiro.hm <- plot_eval_heatmap(counts_df = all_counts, 
                                scores_df = all_scores, 
                                plottingOrder = plottingOrder, 
                                select_cohort = "Schapiro")

wagner.hm <- plot_eval_heatmap(counts_df = all_counts, 
                               scores_df = all_scores, 
                               plottingOrder = plottingOrder, 
                               select_cohort = "Wagner")

zurich.hm <- plot_eval_heatmap(counts_df = all_counts, 
                               scores_df = all_scores, 
                               plottingOrder = plottingOrder, 
                               select_cohort = "Zurich")

lin.hm <- plot_eval_heatmap(counts_df = all_counts, 
                            scores_df = all_scores, 
                            plottingOrder = plottingOrder, 
                            select_cohort = "Lin")



#basel.hm + schapiro.hm + wagner.hm + zurich.hm + lin.hm



pdf(file = args$output_heatmap, width = 14, height = 7)
  hm_list = basel.hm + schapiro.hm + wagner.hm + zurich.hm + lin.hm
  draw(hm_list, ht_gap = unit(0.5, "cm"))
dev.off()

