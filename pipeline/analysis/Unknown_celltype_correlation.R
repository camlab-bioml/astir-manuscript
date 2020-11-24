#!/usr/local/bin/Rscript

library(SingleCellExperiment)
library(tidyverse)
library(scater)
library(ggplot2)
library(corrplot)
library(Hmisc)
library(ggcorrplot)
library(devtools)
devtools::load_all("~/taproom/")

### [FUNCTIONS] ####
get_celltypes <- function(prob_mat, thresh = 0.5) {
  if(is.data.frame(prob_mat)) {
    prob_mat <- as.matrix(prob_mat)
  }
  celltypes <- apply(prob_mat, 1, function(x) {
    colnames(prob_mat)[which.max(x)]
  })
  
  celltypes[rowMaxs(prob_mat) < thresh] <- "Unknown"

  if("Epithelial (basal)" %in% colnames(prob_mat)){
    basal <- prob_mat[, "Epithelial (basal)"]
  }else{
    basal <- rep(0, nrow(prob_mat))
  }

  if("Epithelial (other)" %in% colnames(prob_mat)){
    other <- prob_mat[, "Epithelial (other)"]
  }else{
    other <- rep(0, nrow(prob_mat))
  }

  if("Epithelial (luminal)" %in% colnames(prob_mat)){
    luminal <- prob_mat[, "Epithelial (luminal)"]
  }else{
    luminal <- rep(0, nrow(prob_mat))
  }
  
  celltypes[(basal + other + luminal > thresh) & 
              celltypes == "Unknown"] <- "Epithelial (indeterminate)"
  celltypes
}
assignIdentity <- function(raw.sce, types, states, dimReduct = F){
  #### REad in data
  if(is.character(raw.sce)){
    sce <- readRDS(raw.sce)
    
    ### Make sure wagner has id field
    if(any(colnames(colData(sce)) == "id") == F & grepl("wagner", raw.sce)){
      id <- sce %>% colData() %>% as.data.frame() %>% rownames()
      sce$id <- id
    }
  }else{
    sce <- raw.sce
  }
  
  #### ASSIGN CELL TYPES
  types <- read_csv(types)
  
  types_mat <- select(types, -X1) %>% 
    as.data.frame()
  rownames(types_mat) <- types$X1
  
  assignments <- get_celltypes(types_mat, thresh = 0.5) %>% as.data.frame()
  colnames(assignments) <- "cell_type"
  assignments$id <- rownames(assignments)
  
  #### CELL STATE
  states <- read_csv(states)
  states_df <- select(states, -X1) %>% 
    as.data.frame()
  
  rownames(states_df) <- states$X1
  
  colnames(states_df) <- colnames(states_df) %>% 
    str_replace_all(" ", "_") %>% 
    str_replace_all("\\(", ".") %>% 
    str_replace_all("\\)", ".") %>% 
    str_replace("HALLMARK_", "")
  
  stateNames <- colnames(states_df) #%>% 
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  states_df <- apply(states_df, 2, range01) %>% 
    as.data.frame()
  states_df$id <- rownames(states_df)
  
  
  typeState <- inner_join(assignments, states_df) %>% 
    pivot_longer(is.numeric, names_to = "state", values_to = "activation") %>% 
    group_by(cell_type, state) %>% 
    dplyr::mutate(activation.med = median(activation)) %>% 
    dplyr::mutate(fac = factor(ifelse(activation < activation.med, "lo", "hi"))) %>% 
    select(-activation.med) %>% 
    mutate(pathway_activation = paste0(state, "_", fac)) %>% 
    pivot_wider(names_from = state, values_from = c(activation, fac, pathway_activation), 
                names_glue = "{state}.{.value}") %>% 
    unite(type_summary, ends_with("pathway_activation")) %>%
    mutate(type_state = paste0(cell_type, "_", type_summary))
  
  
  if("B cell" %in% levels(typeState$cell_type) &
     "T cell" %in% levels(typeState$cell_type)){
    typeState$cell_type <- as.character(typeState$cell_type)
    
    typeState <- typeState %>% 
      mutate(cell_type = replace(cell_type, cell_type == "B cell", "B cells")) %>% 
      mutate(cell_type = replace(cell_type, cell_type == "T cell", "T cells"))
    
    typeState$cell_type <- as.factor(typeState$cell_type)
  }
  
  cellID <- typeState$id
  typeState <- typeState %>% 
    select(-id) %>% 
    as.data.frame()
  
  rownames(typeState) <- cellID
  
  newColNames <- typeState %>% colnames()
  
  colData(sce)[newColNames] <- typeState[colnames(sce), ]
  
  ### [PCA & UMAP PROJECTIONS] #####
  if(dimReduct == T){
    sce <- runPCA(sce, ncomponents = 10)
    sce <- runUMAP(sce)
  }
  
  return(list(sce = sce, pathways = paste0(stateNames, ".activation")))
}

### [READ IN DATA] ####
args <- commandArgs(trailingOnly = TRUE)
cells <- args[1]
type <- args[2]
state <- args[3]
markers <- read_markers(args[4])
cohort <- args[5]
output_dir <- args[6]

sce <- assignIdentity(cells, type, state)$sce

### [SUBSET DATA] #####
unknown.sce <- sce[, sce$cell_type == "Unknown" | sce$cell_type == "Other"]
unknown.sce <- unknown.sce[unique(unlist(markers$cell_types)),]
assigned.sce <- sce[, sce$cell_type != "Unknown" & sce$cell_type != "Other"]
assigned.sce <- assigned.sce[unique(unlist(markers$cell_types)),]

### [ANALYSIS & PLOTTING] #####
# Create plotting order
if(cohort == "basel" | cohort == "zurich1"){
  protein.order <- c("CD45", "CD20", "CD3", "E-Cadherin", "pan Cytokeratin", 
                    "Cytokeratin 5", "Cytokeratin 14", "Cytokeratin 7",
                    "Cytokeratin 8/18", "Cytokeratin 19", "CD68", "vWF", 
                    "Vimentin", "Fibronectin")
}else if(cohort == "wagner"){
  protein.order <- c("CD45", "CD3", "CD68", "CD24", "CD31", 
                     "CD49f", "ECadherin", "EpCAM", "panK", "K14", "K5", 
                     "SMA", "K7", "K8K18", "FAP")
}else if(cohort == "schapiro"){
  protein.order <- c("E-cadherin", "EpCAM", "Cytokeratin7", "Cytokeratin8-18", "CD68",
                     "Vimentin", "Fibronectin")
}

# Bonferroni correction
penalty <- (length(protein.order) * length(protein.order)) / 2 - length(protein.order)
sig.levels <- c(0.001, 0.01, 0.05)
corrected.sig <- sig.levels / penalty

# Create ordered matrices
unknown.mat <- t(logcounts(unknown.sce))[, protein.order]
assigned.mat <- t(logcounts(assigned.sce))[, protein.order]

# Create correlations
unknown.cor <- Hmisc::rcorr(unknown.mat)
assigned.cor <- Hmisc::rcorr(assigned.mat)


# Create plots
pdf(paste0(output_dir, "expression_correlation_unknown_cells_", cohort, ".pdf"))
ggcorrplot(unknown.cor$r, hc.order = F, type = "lower",
           outline.col = "white", p.mat = unknown.cor$P,
           sig.level = 0.05, insig = "blank",
           ggtheme = taproom::astir_paper_theme,
           title = "Unknown cells",
           colors = c("#6D9EC1", "white", "#E46726"))
dev.off()

pdf(paste0(output_dir, "expression_correlation_assigned_cells_", cohort, ".pdf"))
ggcorrplot(assigned.cor$r, hc.order = F, type = "lower",
           outline.col = "white", p.mat = assigned.cor$P,
           sig.level = 0.05, insig = "blank",
           ggtheme = taproom::astir_paper_theme,
           title = "Assigned cells",
           colors = c("#6D9EC1", "white", "#E46726"))
dev.off()