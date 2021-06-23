#!/usr/local/bin/Rscript
library(SingleCellExperiment)
library(tidyverse)
library(scater)
library(ggplot2)
library(Hmisc)
library(devtools)
library(stringr)
library(reshape2)
devtools::load_all("../taproom/")

### [FUNCTIONS] ####
cors <- function(df) { 
  # turn all three matrices (r, n, and P into a data frame)
  M <- Hmisc::rcorr(as.matrix(df))
  # return the three data frames in a list return(Mdf)
  Mdf <- map(M, ~data.frame(.x, check.names =  FALSE))
  Mdf
}

select_tri <- function(mat, triangle){
  if(triangle == "lower"){
    mat[lower.tri(mat)] <- NA
  }else if(triangle == "upper"){
    mat[upper.tri(mat)] <- NA
  }
  as.matrix(mat)
}

formated_cors <- function(mat, triangle){
  mat %>% cors() %>% 
    map(~select_tri(., triangle)) %>% 
    melt(na.rm = TRUE) %>% 
    as_tibble %>% 
    pivot_wider(names_from = L1, values_from = "value") %>% 
    mutate(sig_level = case_when(
      P < corrected.sig[1] ~ "***",
      P < corrected.sig[2] & P > corrected.sig[1] ~ "**",
      P < corrected.sig[3] & P > corrected.sig[2] ~ "*",
      TRUE ~ NA_character_))
}

### [READ IN DATA] ####
args <- commandArgs(trailingOnly = TRUE)
cells <- args[1]
#cells <- "output/v6/schapiro_subset/schapiro_subset_sce.rds"
type <- args[2]
#type <- "output/v6/schapiro_subset/schapiro_subset_assignments_type.csv"
#state <- args[3]
markers <- read_markers(args[3])
#markers <- read_markers("markers/schapiro-markers.yml")
cohort <- args[4]
#cohort <- "schapiro"
output_dir <- args[5]

thresh <- 0.5

# read in data 
sce <- readRDS(cells)
type <- read_csv(type)

type$cell_type <- get_celltypes(select(type, -X1)) 
type <- select(type, X1, cell_type) %>%
  column_to_rownames("X1")
colData(sce)["cell_type"] <- type[colnames(sce),]


### [SUBSET DATA] #####
unknown.sce <- sce[, sce$cell_type == "Unknown" | sce$cell_type == "Other"]
unknown.sce <- unknown.sce[unique(unlist(markers$cell_types)),] %>% 
  logcounts() %>% t()
assigned.sce <- sce[, sce$cell_type != "Unknown" & sce$cell_type != "Other"]
assigned.sce <- assigned.sce[unique(unlist(markers$cell_types)),] %>% 
  logcounts() %>% t()

# Create plotting order
if(cohort == "basel" | cohort == "zurich1"){
  protein.order <- c("CD45", "CD20", "CD3", "CD68", "E-Cadherin", "pan Cytokeratin", 
                     "Cytokeratin 5", "Cytokeratin 14", "Cytokeratin 7",
                     "Cytokeratin 8/18", "Cytokeratin 19", "vWF", 
                     "Vimentin", "Fibronectin", "SMA")
  if(cohort == "basel"){
    title <- "Basel"
  }else{
    title <- "Zurich"
  }
}else if(cohort == "wagner"){
  protein.order <- c("SMA", "Vimentin", "CD45", "CD3", "CD68", "CD31", 
                     "CD49f", "ECadherin", "EpCAM", "panK", "K14", "K5", "K7", "K8K18")
  correct_names <- c("SMA", "Vimentin", "CD45", "CD3", "CD68", "CD31", 
                     "CD49f", "E-Cadherin", "EpCAM", "pan Cytokeratin", "Cytokeratin 14", 
                     "Cytokeratin 5", "Cytokeratin 7", "Cytokeratin 8/18")
  title <- "Wagner"
}else if(cohort == "schapiro"){
  protein.order <- c("CD68", "E-cadherin", "EpCAM", "Cytokeratin7", "Cytokeratin8-18",
                     "Vimentin", "Fibronectin")
  correct_names <- c("CD68", "E-Cadherin", "EpCAM", "Cytokeratin 7", "Cytokeratin 8/18",
                     "Vimentin", "Fibronectin")
  title <- "Schapiro"
}else if(cohort == "lin_cycif"){
  protein.order <- c("Vimentin", "Ecad", "Keratin", "Ki67", "PDL1", "CD45", "CD4")
  correct_names <- c("Vimentin", "E-Cadherin", "Keratin", "Ki67", "PDL1", "CD45", "CD4")

  title <- "Lin"
}

unknown.mat <- unknown.sce[, protein.order]
assigned.mat <- assigned.sce[, protein.order]

if(cohort == "wagner" | cohort == "schapiro" | cohort == "lin_cycif"){
  colnames(unknown.mat) <- correct_names
  colnames(assigned.mat) <- correct_names
  protein.order <- correct_names
}


penalty <- (length(protein.order) * length(protein.order)) / 2 - length(protein.order)
sig.levels <- c(0.001, 0.01, 0.05)
corrected.sig <- sig.levels / penalty

unknown.cor <- formated_cors(unknown.mat, "lower")
assigned.cor <- formated_cors(assigned.mat, "upper")

correlations <- rbind(unknown.cor, assigned.cor)

x_labels <- levels(correlations$Var2)
y_labels <- levels(correlations$Var1)

pdf(paste0(output_dir, "expression_correlation_", cohort, ".pdf"), height = 9.5, width = 9)
correlations %>% 
  ggplot(aes(as.numeric(Var2), as.numeric(Var1), fill = r)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low="#6D9EC1",mid="white",high="#E46726", 
                       midpoint = 0, limits=c(-1,1), name = "Pearson\nCorrelation  ") +
  labs(y = "Unknown and Other cell types", x = "Assigned cell types") +
  astir_paper_theme() + 
  scale_x_continuous(position = "bottom", expand=c(0,0), sec.axis = dup_axis(),
                     breaks = 1:length(x_labels), labels = x_labels) +
  scale_y_continuous(position = "left", expand=c(0,0), sec.axis = dup_axis(),
                     breaks = 1:length(y_labels), labels = y_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 24, hjust = 1),
        axis.text.y = element_text(size = 24),
        axis.title=element_text(size=30),
        plot.title = element_text(size=40),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.x.bottom = element_blank(), # remove titles
        axis.title.y.left = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.key.width = unit(2,"cm"),
        legend.title=element_text(size=25),
        legend.text=element_text(size=22))
  coord_fixed()
dev.off()

