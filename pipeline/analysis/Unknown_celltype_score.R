#!/usr/local/bin/Rscript
library(viridis)
library(ggplot2)
library(devtools)
devtools::load_all("~/taproom/")

args <- commandArgs(trailingOnly = TRUE)
cells <- readRDS(args[1])
cells <- readRDS("../imc-prototype-analysis/raw_data/zurich1/zurich1_subset_sce.rds")
#types <- read_csv(args[2]) %>% 
types <- read_csv("../imc-prototype-analysis/raw_data/zurich1/zurich1_subset_assignments_type.csv") %>% 
  column_to_rownames("X1")

cohort <- args[3]
output_dir <- args[4]

# Assign cell types
types$cell_type <- get_celltypes(types)
cells$cell_type <- types[colnames(cells), ]$cell_type

# Get expression - unknown cells
unknown.cells <- cells[, cells$cell_type == "Unknown"]
unknown.expr <- t(logcounts(unknown.cells)) %>% as.data.frame()

## Assigned cells
assigned.cells <- cells[, cells$cell_type != "Unknown"]
assigned.expr <- t(logcounts(assigned.cells)) %>% as.data.frame()


# Normalize expression data
thresh <- 2

unknown.expr <- scale(unknown.expr)
unknown.expr[unknown.expr > thresh] <- thresh
unknown.expr[unknown.expr < -thresh] <- -thresh

assigned.expr <- scale(assigned.expr)
assigned.expr[assigned.expr > thresh] <- thresh
assigned.expr[assigned.expr < -thresh] <- -thresh

unknown.expr <- as.data.frame(unknown.expr)
assigned.expr <- as.data.frame(assigned.expr)

# Create list of mutually exclusive genes
leuk <- c("CD3", "CD20", "CD45", "CD68")
epi <- c("E-Cadherin", 
         "pan Cytokeratin",
         "Cytokeratin 7",
         "Cytokeratin 8/18",
         "Cytokeratin 19",
         "Cytokeratin 5",
         "Cytokeratin 14")

gene_pairs_mutual_df <- expand.grid(leuk, epi)

gene_pairs <- lapply(seq_len(nrow(gene_pairs_mutual_df)), function(i) {
  c(as.character(gene_pairs_mutual_df[i,1]), 
    as.character(gene_pairs_mutual_df[i,2]))
})


# Calculate scores for each cell - unknown cells
tmp_score <- lapply(gene_pairs, function(x){
  unknown.expr[x[1]] * unknown.expr[x[2]]
}) 

unknown_score <- Reduce("+",tmp_score)
unknown_score$Cell <- "Unknown"

# assigned cells
tmp_score <- lapply(gene_pairs, function(x){
  assigned.expr[x[1]] * assigned.expr[x[2]]
}) 

assigned_score <- Reduce("+", tmp_score)
assigned_score$Cell <- "Assigned"

# Combine all scores
total <- rbind(unknown_score, assigned_score) %>% as.data.frame()
colnames(total) <- c("score", "Cell")

# Plot scores as histograms
pdf(paste0(output_dir, "Unknown_celltype_mutual_exclusivity_score_hist_", cohort, ".pdf"))
ggplot(total, aes(x = score, color = Cell, fill = Cell)) +
  geom_histogram(position = "identity", bins = 100) +
  scale_fill_viridis(discrete = T) +
  scale_color_viridis(discrete = T) +
  astir_paper_theme()
dev.off()

pdf(paste0(output_dir, "Unknown_celltype_mutual_exclusivity_score_qqplot_", cohort, ".pdf"))
qqplot(unknown_score$CD3, assigned_score$CD3, xlab = "Unknown score", ylab = "Assigned score")
abline(0, 1, col="red")
dev.off()
