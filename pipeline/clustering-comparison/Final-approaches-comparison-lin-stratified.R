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
parser$add_argument('--bladder_astir_assignments', type = 'character')
parser$add_argument('--breast_astir_assignments', type = 'character')
parser$add_argument('--GI_stomach_astir_assignments', type = 'character')
parser$add_argument('--GI_colorectal_astir_assignments', type = 'character')
parser$add_argument('--kidney_astir_assignments', type = 'character')
parser$add_argument('--liver_astir_assignments', type = 'character')
parser$add_argument('--lung_astir_assignments', type = 'character')
parser$add_argument('--lymph_astir_assignments', type = 'character')
parser$add_argument('--ovary_astir_assignments', type = 'character')
parser$add_argument('--pancreas_astir_assignments', type = 'character')
parser$add_argument('--prostate_astir_assignments', type = 'character')
parser$add_argument('--skin_astir_assignments', type = 'character')
parser$add_argument('--uterus_astir_assignments', type = 'character')

parser$add_argument('--bladder_sce', type = 'character')
parser$add_argument('--breast_sce', type = 'character')
parser$add_argument('--GI_stomach_sce', type = 'character')
parser$add_argument('--GI_colorectal_sce', type = 'character')
parser$add_argument('--kidney_sce', type = 'character')
parser$add_argument('--liver_sce', type = 'character')
parser$add_argument('--lung_sce', type = 'character')
parser$add_argument('--lymph_sce', type = 'character')
parser$add_argument('--ovary_sce', type = 'character')
parser$add_argument('--pancreas_sce', type = 'character')
parser$add_argument('--prostate_sce', type = 'character')
parser$add_argument('--skin_sce', type = 'character')
parser$add_argument('--uterus_sce', type = 'character')

parser$add_argument('--bladder_files', type = 'character', nargs = '+')
parser$add_argument('--breast_files', type = 'character', nargs = '+')
parser$add_argument('--GI_stomach_files', type = 'character', nargs = '+')
parser$add_argument('--GI_colorectal_files', type = 'character', nargs = '+')
parser$add_argument('--kidney_files', type = 'character', nargs = '+')
parser$add_argument('--liver_files', type = 'character', nargs = '+')
parser$add_argument('--lung_files', type = 'character', nargs = '+')
parser$add_argument('--lymph_files', type = 'character', nargs = '+')
parser$add_argument('--ovary_files', type = 'character', nargs = '+')
parser$add_argument('--pancreas_files', type = 'character', nargs = '+')
parser$add_argument('--prostate_files', type = 'character', nargs = '+')
parser$add_argument('--skin_files', type = 'character', nargs = '+')
parser$add_argument('--uterus_files', type = 'character', nargs = '+')

parser$add_argument('--bladder_absent_acdc', type = 'character', nargs = '+')
parser$add_argument('--breast_absent_acdc', type = 'character', nargs = '+')
parser$add_argument('--GI_stomach_absent_acdc', type = 'character', nargs = '+')
parser$add_argument('--GI_colorectal_absent_acdc', type = 'character', nargs = '+')
parser$add_argument('--kidney_absent_acdc', type = 'character', nargs = '+')
parser$add_argument('--liver_absent_acdc', type = 'character', nargs = '+')
parser$add_argument('--lung_absent_acdc', type = 'character', nargs = '+')
parser$add_argument('--lymph_absent_acdc', type = 'character', nargs = '+')
parser$add_argument('--ovary_absent_acdc', type = 'character', nargs = '+')
parser$add_argument('--pancreas_absent_acdc', type = 'character', nargs = '+')
parser$add_argument('--prostate_absent_acdc', type = 'character', nargs = '+')
parser$add_argument('--skin_absent_acdc', type = 'character', nargs = '+')
parser$add_argument('--uterus_absent_acdc', type = 'character', nargs = '+')

parser$add_argument('--bladder_no_consider_acdc', type = 'character', nargs = '+')
parser$add_argument('--breast_no_consider_acdc', type = 'character', nargs = '+')
parser$add_argument('--GI_stomach_no_consider_acdc', type = 'character', nargs = '+')
parser$add_argument('--GI_colorectal_no_consider_acdc', type = 'character', nargs = '+')
parser$add_argument('--kidney_no_consider_acdc', type = 'character', nargs = '+')
parser$add_argument('--liver_no_consider_acdc', type = 'character', nargs = '+')
parser$add_argument('--lung_no_consider_acdc', type = 'character', nargs = '+')
parser$add_argument('--lymph_no_consider_acdc', type = 'character', nargs = '+')
parser$add_argument('--ovary_no_consider_acdc', type = 'character', nargs = '+')
parser$add_argument('--pancreas_no_consider_acdc', type = 'character', nargs = '+')
parser$add_argument('--prostate_no_consider_acdc', type = 'character', nargs = '+')
parser$add_argument('--skin_no_consider_acdc', type = 'character', nargs = '+')
parser$add_argument('--uterus_no_consider_acdc', type = 'character', nargs = '+')

parser$add_argument('--lin_markers', type = 'character')
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

# Cluster assignments
create_counts <- function(expression, markers, cohort, method){
  GSVA_Astir_celltypes <- expression %>% 
    filter(!cluster %in% c("Other", "Unknown")) %>% 
    assign_clusters(markers)
  
  astir_counts <- manual_cluster_assignment(GSVA_Astir_celltypes$expression, markers) %>%
    group_by(Manual_cell_type) %>%
    tally() %>% 
    mutate(cohort = cohort, method = method) %>%
    dplyr::rename("cell_type" = "Manual_cell_type")

  # astir_counts <- GSVA_Astir_celltypes$assignment %>% 
  #   group_by(GSVA_cell_type) %>% 
  #   tally() %>% 
  #   mutate(cohort = cohort, method = "Astir") %>%
  #   dplyr::rename("cell_type" = "GSVA_cell_type")
  astir_counts <- astir_counts[,c("cohort", "method", "cell_type", "n")]
}

# Read in other methods assignments
read_in_cohort <- function(files, cohort){
    data <- lapply(files, read_csv) %>% 
                bind_rows() %>% 
                mutate(cohort = cohort)
    
    data <- data %>%
        select(-c(id, GSVA_cell_type)) %>% 
        distinct(.keep_all = TRUE)

    data <- data %>%
        mutate(params = str_replace(params, "_", " ")) %>% 
        mutate(params = str_replace(params, "_", "-")) %>% 
        mutate(method = paste0(method, "-", params)) %>% 
        select(-params) %>%
        dplyr::rename("cell_type" = "Manual_cell_type")
    
    data
}


plot_eval_heatmap <- function(counts_df, scores_df, plottingOrder, select_cohort){
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
                     width = unit(1.5, "cm"))
  
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

# Read in Astir assignments
bladder_astir <- read_astir_assignments(args$bladder_astir_assignments)
breast_astir <- read_astir_assignments(args$breast_astir_assignments)
GI_Stomach_astir <- read_astir_assignments(args$GI_stomach_astir_assignments)
GI_colorectal_astir <- read_astir_assignments(args$GI_colorectal_astir_assignments)
Kidney_astir <- read_astir_assignments(args$kidney_astir_assignments)
Liver_astir <- read_astir_assignments(args$liver_astir_assignments)
Lung_astir <- read_astir_assignments(args$lung_astir_assignments)
Lymph_astir <- read_astir_assignments(args$lymph_astir_assignments)
Ovary_astir <- read_astir_assignments(args$ovary_astir_assignments)
Pancreas_astir <- read_astir_assignments(args$pancreas_astir_assignments)
Prostate_astir <- read_astir_assignments(args$prostate_astir_assignments)
Skin_astir <- read_astir_assignments(args$skin_astir_assignments)
Uterus_astir <- read_astir_assignments(args$uterus_astir_assignments)


# Read in sces
bladder_sce <- readRDS(args$bladder_sce)
breast_sce <- readRDS(args$breast_sce)
GI_stomach_sce <- readRDS(args$GI_stomach_sce)
GI_colorectal_sce <- readRDS(args$GI_colorectal_sce)
kidney_sce <- readRDS(args$kidney_sce)
liver_sce <- readRDS(args$liver_sce)
lung_sce <- readRDS(args$lung_sce)
lymph_sce <- readRDS(args$lymph_sce)
ovary_sce <- readRDS(args$ovary_sce)
pancreas_sce <- readRDS(args$pancreas_sce)
prostate_sce <- readRDS(args$prostate_sce)
skin_sce <- readRDS(args$skin_sce)
uterus_sce <- readRDS(args$uterus_sce)

# Create expression matrices
bladder_expression <- create_expression_mat(bladder_sce, bladder_astir)
breast_expression <- create_expression_mat(breast_sce, breast_astir)
GI_stomach_expression <- create_expression_mat(GI_stomach_sce, GI_Stomach_astir)
GI_colorectal_expression <- create_expression_mat(GI_colorectal_sce, GI_colorectal_astir)
kidney_expression <- create_expression_mat(kidney_sce, Kidney_astir)
liver_expression <- create_expression_mat(liver_sce, Liver_astir)
lung_expression <- create_expression_mat(lung_sce, Lung_astir)
lymph_expression <- create_expression_mat(lymph_sce, Lymph_astir)
ovary_expression <- create_expression_mat(ovary_sce, Ovary_astir)
pancreas_expression <- create_expression_mat(pancreas_sce, Pancreas_astir)
prostate_expression <- create_expression_mat(prostate_sce, Prostate_astir)
skin_expression <- create_expression_mat(skin_sce, Skin_astir)
uterus_expression <- create_expression_mat(uterus_sce, Uterus_astir)

# read in markers
lin_markers <- read_markers(args$lin_markers)


bladder_counts <- create_counts(bladder_expression, lin_markers, "Bladder", "Astir")
breast_counts <- create_counts(breast_expression, lin_markers, "Breast", "Astir")
GI_stomach_counts <- create_counts(GI_stomach_expression, lin_markers, "GI stomach", "Astir")
GI_colorectal_counts <- create_counts(GI_colorectal_expression, lin_markers, "GI colorectal", "Astir")
kidney_counts <- create_counts(kidney_expression, lin_markers, "Kidney", "Astir")
liver_counts <- create_counts(liver_expression, lin_markers, "Liver", "Astir")
lung_counts <- create_counts(lung_expression, lin_markers, "Lung", "Astir")
lymph_counts <- create_counts(lymph_expression, lin_markers, "Lymph", "Astir")
ovary_counts <- create_counts(ovary_expression, lin_markers, "Ovary", "Astir")
pancreas_counts <- create_counts(pancreas_expression, lin_markers, "Pancreas", "Astir")
prostate_counts <- create_counts(prostate_expression, lin_markers, "Prostate", "Astir")
skin_counts <- create_counts(skin_expression, lin_markers, "Skin", "Astir")
uterus_counts <- create_counts(uterus_expression, lin_markers, "Uterus", "Astir")


bladder <- read_in_cohort(args$bladder_files, "Bladder")
breast <- read_in_cohort(args$breast_files, "Breast")
GI_stomach <- read_in_cohort(args$GI_stomach_files, "GI stomach")
GI_colorectal <- read_in_cohort(args$GI_colorectal_files, "GI colorectal")
kidney <- read_in_cohort(args$kidney_files, "Kidney")
liver <- read_in_cohort(args$liver_files, "Liver")
lung <- read_in_cohort(args$lung_files, "Lung")
lymph <- read_in_cohort(args$lymph_files, "Lymph")
ovary <- read_in_cohort(args$ovary_files, "Ovary")
pancreas <- read_in_cohort(args$pancreas_files, "Pancreas")
prostate <- read_in_cohort(args$prostate_files, "Prostate")
skin <- read_in_cohort(args$skin_files, "Skin")
uterus <- read_in_cohort(args$uterus_files, "Uterus")


acdc_counts <- function(file, sce, markers, cohort, method){
  acdc_assignment <- read_tsv(file) %>% 
    dplyr::rename("cluster" = "cell_type") %>%
    dplyr::rename("id" = "cell_id")

  expression <- create_expression_mat(sce, acdc_assignment)

  create_counts(expression, markers, cohort, method)
}

bladder_no_consider <- acdc_counts(args$bladder_no_consider_acdc, bladder_sce, lin_markers, "Bladder", "ACDC no-consider")
bladder_absent <- acdc_counts(args$bladder_absent_acdc, bladder_sce, lin_markers, "Bladder", "ACDC absent")
breast_no_consider <- acdc_counts(args$breast_no_consider_acdc, breast_sce, lin_markers, "Breast", "ACDC no-consider")
breast_absent <- acdc_counts(args$breast_absent_acdc, breast_sce, lin_markers, "Breast", "ACDC absent")

GI_stomach_no_consider <- acdc_counts(args$GI_stomach_no_consider_acdc, GI_stomach_sce, lin_markers, "GI stomach", "ACDC no-consider")
GI_stomach_absent <- acdc_counts(args$GI_stomach_absent_acdc, GI_stomach_sce, lin_markers, "GI stomach", "ACDC absent")

GI_colorectal_no_consider <- acdc_counts(args$GI_colorectal_no_consider_acdc, GI_colorectal_sce, lin_markers, "GI colorectal", "ACDC no-consider")
GI_colorectal_absent <- acdc_counts(args$GI_colorectal_absent_acdc, GI_colorectal_sce, lin_markers, "GI colorectal", "ACDC absent")

kidney_no_consider <- acdc_counts(args$kidney_no_consider_acdc, kidney_sce, lin_markers, "Kidney", "ACDC no-consider")
kidney_absent <- acdc_counts(args$kidney_absent_acdc, kidney_sce, lin_markers, "Kidney", "ACDC absent")

liver_no_consider <- acdc_counts(args$liver_no_consider_acdc, liver_sce, lin_markers, "Liver", "ACDC no-consider")
liver_absent <- acdc_counts(args$liver_absent_acdc, liver_sce, lin_markers, "Liver", "ACDC absent")

lung_no_consider <- acdc_counts(args$lung_no_consider_acdc, lung_sce, lin_markers, "Lung", "ACDC no-consider")
lung_absent <- acdc_counts(args$lung_absent_acdc, lung_sce, lin_markers, "Lung", "ACDC absent")

lymph_no_consider <- acdc_counts(args$lymph_no_consider_acdc, lymph_sce, lin_markers, "Lymph", "ACDC no-consider")
lymph_absent <- acdc_counts(args$lymph_absent_acdc, lymph_sce, lin_markers, "Lymph", "ACDC absent")

ovary_no_consider <- acdc_counts(args$ovary_no_consider_acdc, ovary_sce, lin_markers, "Ovary", "ACDC no-consider")
ovary_absent <- acdc_counts(args$ovary_absent_acdc, ovary_sce, lin_markers, "Ovary", "ACDC absent")

pancreas_no_consider <- acdc_counts(args$pancreas_no_consider_acdc, pancreas_sce, lin_markers, "Pancreas", "ACDC no-consider")
pancreas_absent <- acdc_counts(args$pancreas_absent_acdc, pancreas_sce, lin_markers, "Pancreas", "ACDC absent")

prostate_no_consider <- acdc_counts(args$prostate_no_consider_acdc, prostate_sce, lin_markers, "Prostate", "ACDC no-consider")
prostate_absent <- acdc_counts(args$prostate_absent_acdc, prostate_sce, lin_markers, "Prostate", "ACDC absent")

skin_no_consider <- acdc_counts(args$skin_no_consider_acdc, skin_sce, lin_markers, "Skin", "ACDC no-consider")
skin_absent <- acdc_counts(args$skin_absent_acdc, skin_sce, lin_markers, "Skin", "ACDC absent")

uterus_no_consider <- acdc_counts(args$uterus_no_consider_acdc, uterus_sce, lin_markers, "Uterus", "ACDC no-consider")
uterus_absent <- acdc_counts(args$uterus_absent_acdc, uterus_sce, lin_markers, "Uterus", "ACDC absent")


acdc_count <- bind_rows(bladder_no_consider, bladder_absent,
                        breast_no_consider, breast_absent,
                        GI_stomach_no_consider, GI_stomach_absent,
                        GI_colorectal_no_consider, GI_colorectal_absent,
                        kidney_no_consider, kidney_absent,
                        liver_no_consider, liver_absent,
                        lung_no_consider, lung_absent,
                        lymph_no_consider, lymph_absent,
                        ovary_no_consider, ovary_absent,
                        pancreas_no_consider, pancreas_absent,
                        prostate_no_consider, prostate_absent,
                        skin_no_consider, skin_absent,
                        uterus_no_consider, uterus_absent)


all_cohorts <- bind_rows(bladder, breast, GI_stomach, GI_colorectal, 
                        kidney, liver, lung, lymph, ovary, pancreas,
                        prostate, skin, uterus)

# All counts for other methods
all_counts <- all_cohorts %>% 
  dplyr::group_by(cohort, method, cell_type) %>% 
  tally() %>% 
  ungroup()

all_counts$n[is.na(all_counts$cell_type)] <- NA

# Add other method counts and astir counts
all_counts <- bind_rows(all_counts, bladder_counts, breast_counts, 
                        GI_stomach_counts, GI_colorectal_counts,
                        kidney_counts, liver_counts, lung_counts,
                        lymph_counts, ovary_counts, pancreas_counts,
                        prostate_counts, skin_counts, uterus_counts, acdc_count)


# Calculate scores for all methods across all cohorts
# First I need a list of cell types for each cohort
cohort_markers <- list(list("Lin", lin_markers$cell_types))

methods <- unique(all_counts$method)
cell_types <- names(lin_markers$cell_types)

cohort_cell_types <- expand.grid(list(methods, cell_types, unique(all_counts$cohort)))
colnames(cohort_cell_types) <- c("method", "cell_type", "cohort")

# Join with actual scores so that those cell types for which no cluster was found
# Can also be taken into account
all_scores <- full_join(all_counts, cohort_cell_types)
all_scores$n[is.na(all_scores$n)] <- 0

# Calculate scores
all_scores <- all_scores %>% 
  ungroup() %>% 
  mutate(score = case_when(
    n == 1 ~ 1,
    n != 1 ~ -1
  )) %>% group_by(cohort, method) %>% 
  dplyr::summarise(score = sum(score)) %>% 
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

bladder.hm <- plot_eval_heatmap(counts_df = all_counts, 
                              scores_df = all_scores, 
                              plottingOrder = plottingOrder, 
                              select_cohort = "Bladder")

breast.hm <- plot_eval_heatmap(counts_df = all_counts, 
                              scores_df = all_scores, 
                              plottingOrder = plottingOrder, 
                              select_cohort = "Breast")

GI_stomach.hm <- plot_eval_heatmap(counts_df = all_counts, 
                              scores_df = all_scores, 
                              plottingOrder = plottingOrder, 
                              select_cohort = "GI stomach")

GI_colorectal.hm <- plot_eval_heatmap(counts_df = all_counts, 
                              scores_df = all_scores, 
                              plottingOrder = plottingOrder, 
                              select_cohort = "GI colorectal")

kidney.hm <- plot_eval_heatmap(counts_df = all_counts, 
                              scores_df = all_scores, 
                              plottingOrder = plottingOrder, 
                              select_cohort = "Kidney")

liver.hm <- plot_eval_heatmap(counts_df = all_counts, 
                              scores_df = all_scores, 
                              plottingOrder = plottingOrder, 
                              select_cohort = "Liver")

lung.hm <- plot_eval_heatmap(counts_df = all_counts, 
                              scores_df = all_scores, 
                              plottingOrder = plottingOrder, 
                              select_cohort = "Lung")

lymph.hm <- plot_eval_heatmap(counts_df = all_counts, 
                              scores_df = all_scores, 
                              plottingOrder = plottingOrder, 
                              select_cohort = "Lymph")

ovary.hm <- plot_eval_heatmap(counts_df = all_counts, 
                              scores_df = all_scores, 
                              plottingOrder = plottingOrder, 
                              select_cohort = "Ovary")

pancreas.hm <- plot_eval_heatmap(counts_df = all_counts, 
                              scores_df = all_scores, 
                              plottingOrder = plottingOrder, 
                              select_cohort = "Pancreas")

prostate.hm <- plot_eval_heatmap(counts_df = all_counts, 
                              scores_df = all_scores, 
                              plottingOrder = plottingOrder, 
                              select_cohort = "Prostate")

skin.hm <- plot_eval_heatmap(counts_df = all_counts, 
                              scores_df = all_scores, 
                              plottingOrder = plottingOrder, 
                              select_cohort = "Skin")

uterus.hm <- plot_eval_heatmap(counts_df = all_counts, 
                              scores_df = all_scores, 
                              plottingOrder = plottingOrder, 
                              select_cohort = "Uterus")

pdf(file = args$output_heatmap, width = 27, height = 8)
  hm_list = bladder.hm + breast.hm + GI_stomach.hm + GI_colorectal.hm + kidney.hm +
            liver.hm + lung.hm + lymph.hm + ovary.hm + pancreas.hm + prostate.hm + skin.hm + 
            uterus.hm
  draw(hm_list, ht_gap = unit(0.5, "cm"))
dev.off()