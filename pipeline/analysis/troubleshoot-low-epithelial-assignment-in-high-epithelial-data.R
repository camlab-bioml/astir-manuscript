library(tidyverse)
library(ggalluvial)
library(devtools)
library(SingleCellExperiment)
library(broom)
library(ComplexHeatmap)
source("scripts/functions.R")
devtools::load_all("../taproom/")

# Read in phenograph and astir assignments
phenograph_assignments <- lapply(snakemake@input[['phenograph_assignments']], read_csv) %>%
  bind_rows()
#phenograph_path <- "output/squirrel/results/epithelial_overclustering/"
#phenograph_files <- dir(phenograph_path, pattern = "Epithelial_overclustering_Phenograph_clusters-")
#phenograph_assignments <- lapply(paste0(phenograph_path, phenograph_files), read_csv) %>% 
#  bind_rows

#"output/squirrel/astir_assignments/basel_astir_assignments.csv")
astir_assignments <- read_csv(snakemake@input[['astir_assignments']]) %>%
  column_to_rownames("X1") %>% 
  mutate(cell_type = get_celltypes(.)) %>% 
  select(cell_type) %>% 
  rownames_to_column("id") %>% 
  dplyr::rename("astir_cell_type" = "cell_type")

markers <- read_markers(snakemake@input[['markers']])#"markers/jackson-2020-markers-v4.yml")
lineage_markers <- markers$cell_types[c("Epithelial (luminal)", "T cells", "Macrophage", "Stromal")] %>% 
  unlist() %>% unique()
names(lineage_markers) <- NULL

# Read in basel expression data
basel_imc <- readRDS(snakemake@input[['sce']])#"output/squirrel/sces/basel_sce.rds")

# process data
phenograph_assignments <- phenograph_assignments %>% 
  select(-GSVA_cell_type) %>%
  dplyr::rename("phenograph_cell_type" = "Manual_cell_type") %>% 
  pivot_longer(-c("id", "phenograph_cell_type", "method", "params", "percent_epithelial"),
               values_to = "cluster", names_to = "k") %>% 
  drop_na(cluster) %>% 
  left_join(astir_assignments)

# Subset to only the k20
all_k20 <- phenograph_assignments %>% 
  filter(params == "all_markers_k20")

### Alluvial plot 
pdf(snakemake@output[['alluvial']], width = 10, height = 12)
  all_k20 %>% 
    select(percent_epithelial, astir_cell_type, phenograph_cell_type, cluster) %>% 
    group_by(percent_epithelial, astir_cell_type, phenograph_cell_type, cluster) %>% 
    tally() %>% 
    mutate(percent_epithelial = paste0("Percent epithelial: ", percent_epithelial)) %>% 
    ggplot(aes(y = n, axis1 = astir_cell_type, axis2 = cluster, axis3 = phenograph_cell_type)) +
      geom_alluvium(aes(fill = astir_cell_type)) +
      geom_stratum(width = 1/5, fill = "grey", color = "black") +
      geom_text(stat = "stratum", aes(label = after_stat(
        ifelse(as.numeric(stratum) %in% c(1,2,3,4), 
              as.character(stratum), 
              NA)))) +
      labs(y = "Cells", fill = "Astir cell types") +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_discrete(limits = c("Astir\ncell type", "Phenograph\ncluster", "Phenograph\nZ-score\ncell type"), 
                      expand = c(.2, .3)) +
      facet_wrap(~percent_epithelial, ncol = 2) +
      astir_paper_theme() +
      scale_fill_manual(values = jackson_basel_colours()) +
      theme(axis.ticks = element_blank(),
            axis.title = element_text(size = 14),
            strip.text.x = element_text(size = 14))
dev.off()



##### Differential analysis
# Required functions
allDuplicated <- function(vec){
  front <- duplicated(vec)
  back <- duplicated(vec, fromLast = TRUE)
  all_dup <- front + back > 0
  return(all_dup)
}

rem_low_cell_cores <- function(df){
  number_of_proteins <- length(unique(df$Protein))
  
  cores_keep <- df %>% 
    group_by(core, phenograph_cell_type) %>% 
    tally() %>% 
    filter(n > number_of_proteins) %>% 
    pull(core)
  
  cores_keep <- cores_keep[allDuplicated(cores_keep)]
  
  df %>% 
    filter(core %in% cores_keep)
}

plot_lm <- function(df, title = ""){
  df %>% 
    ggplot(aes(x = Protein, y = estimate)) +
    geom_boxplot() +
    ggtitle(title) +
    astir_paper_theme() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_hline(yintercept = 0, color = "red")
}


# Select data with k = 20 and percent epithelial = 30% and 99%
all_k20_30 <- all_k20 %>% 
  filter(percent_epithelial == 0.3)
all_k20_30_cells <- all_k20_30 %>% 
  pull(id)

all_k20_99 <- all_k20 %>% 
  filter(percent_epithelial == 0.99)
all_k20_99_cells <- all_k20_99 %>% 
  pull(id)


# Filter out expression
imc_k20_30_sce <- basel_imc[, all_k20_30_cells]
imc_k20_30_sce$phenograph_cell_type <- all_k20_30[all_k20_30$id == colnames(imc_k20_30_sce),]$phenograph_cell_type

imc_k20_99_sce <- basel_imc[, all_k20_99_cells]
imc_k20_99_sce$phenograph_cell_type <- all_k20_99[all_k20_99$id == colnames(imc_k20_99_sce),]$phenograph_cell_type


get_expression <- function(imc){
  logcounts(imc) %>% t() %>% 
    as.data.frame() %>% 
    select(lineage_markers) %>%
    mutate(phenograph_cell_type = imc$phenograph_cell_type) %>% 
    mutate(core = imc$core) %>% 
    pivot_longer(cols = -c(phenograph_cell_type, core), values_to = "Expression", names_to = "Protein") %>% 
    mutate(Expression = as.numeric(Expression))
}

lm_cell_type_specific <- function(imc_df, cell_type_1, cell_type_2){
  imc_df %>% 
    filter(phenograph_cell_type == cell_type_1 | phenograph_cell_type == cell_type_2) %>% 
    rem_low_cell_cores() %>% 
    mutate(phenograph_cell_type = factor(phenograph_cell_type, levels = c(cell_type_1, cell_type_2))) %>% 
    group_by(core, Protein) %>% 
    do(tidy(lm(Expression ~ phenograph_cell_type, data = .))) %>% 
    filter(term != "(Intercept)")
}

expression_30 <- get_expression(imc_k20_30_sce)
expression_99 <- get_expression(imc_k20_99_sce)
  
  
epi_macro_99 <- lm_cell_type_specific(expression_99, "Macrophage", "Epithelial (luminal)") %>% 
  mutate(percent_epithelial = 0.99, comparison = "Macrophage")
epi_macro_30 <- lm_cell_type_specific(expression_30, "Macrophage", "Epithelial (luminal)") %>% 
  mutate(percent_epithelial = 0.3, comparison = "Macrophage")

epi_stromal_99 <- lm_cell_type_specific(expression_99, "Stromal", "Epithelial (luminal)") %>% 
  mutate(percent_epithelial = 0.99, comparison = "Stromal")
epi_stromal_30 <- lm_cell_type_specific(expression_30, "Stromal", "Epithelial (luminal)") %>% 
  mutate(percent_epithelial = 0.3, comparison = "Stromal")

epi_T_cells_99 <- lm_cell_type_specific(expression_99, "T cells", "Epithelial (luminal)") %>% 
  mutate(percent_epithelial = 0.99, comparison = "T cells")
epi_T_cells_30 <- lm_cell_type_specific(expression_30, "T cells", "Epithelial (luminal)") %>% 
  mutate(percent_epithelial = 0.3, comparison = "T cells")

combined_lms <- rbind(epi_macro_30, epi_macro_99,
                      epi_stromal_30, epi_stromal_99,
                      epi_T_cells_30, epi_T_cells_99)

pdf(snakemake@output[['lm']], width = 7, height = 5)
  combined_lms %>% 
    mutate(percent_epithelial = paste("Percent epithelial:", percent_epithelial)) %>% 
    ggplot(aes(x = Protein, y = estimate, fill = comparison)) +
    geom_boxplot() +
    labs(title = "Epithelial (luminal) vs Macrophages, Stromal & T cells", fill = "Comparison") +
    scale_fill_manual(values = jackson_basel_colours()) +
    astir_paper_theme() +
    labs(y = "Estimate") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_hline(yintercept = 0, color = "red") +
    facet_wrap(~percent_epithelial, ncol = 1)
dev.off()





# Heatmaps
nonScaledHeatmap <- function(sce, scale_to = 1, max_out = NA){
  lc <- t(as.matrix(assay(sce, "logcounts")))
  
  if(is.na(max_out)){
    lc <- apply(lc, 2, function(x){((x-min(x)) / (max(x) - min(x))) * scale_to})
  }else{
    lc[lc > max_out] <- max_out
    lc[lc < -0] <- -0
  }
  
  celltypes <- colData(sce)[["phenograph_cell_type"]]
  
  celltype_annot <- HeatmapAnnotation(`Phenograph Z-score cell type` = celltypes, 
                                      which="column",
                                      col = list(`Phenograph Z-score cell type` = jackson_basel_colours()))  
  
  type_exprs <- Heatmap(t(lc), 
                        name = "Expression",
                        column_title = "Cell",
                        col=viridis(100),
                        top_annotation = celltype_annot,
                        show_column_names = FALSE,
                        column_order = order(celltypes))
  type_exprs
}

pdf(snakemake@output[['scaled_30']], height = 4, width = 7)
  nonScaledHeatmap(imc_k20_30_sce[lineage_markers,])
dev.off()

pdf(snakemake@output[['max_4_30']], height = 4, width = 7)
  nonScaledHeatmap(imc_k20_30_sce[lineage_markers,], max_out = 4)
dev.off()

pdf(snakemake@output[['scaled_99']], height = 4, width = 7)
  nonScaledHeatmap(imc_k20_99_sce[lineage_markers,])
dev.off()

pdf(snakemake@output[['max_4_99']], height = 4, width = 7)
  nonScaledHeatmap(imc_k20_99_sce[lineage_markers,], max_out = 4)
dev.off()

