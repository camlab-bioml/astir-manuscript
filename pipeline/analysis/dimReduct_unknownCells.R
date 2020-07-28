#!/usr/local/bin/Rscript
library(SingleCellExperiment)
library(tidyverse)
library(scater)
library(rsvd)
library(DelayedArray)
library(ggplot2)

### define functions ####
get_celltypes <- function(prob_mat, thresh = 0.7) {
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
  
  assignments <- get_celltypes(types_mat) %>% as.data.frame()
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

### Actual code ######

source('~/FIt-SNE/fast_tsne.R', chdir=T)

args <- commandArgs(trailingOnly = TRUE)

basel_cells <- args[1]
basel_types <- args[2]
basel_states <- args[3]

zurich1_cells <- args[4]
zurich1_types <- args[5]
zurich1_states <- args[6]

wagner_cells <- args[7]
wagner_types <- args[8]
wagner_states <- args[9]
output_dir <- args[10]


### [GET SCE OBJECTS] #####
basel <- assignIdentity(basel_cells, basel_types, basel_states)$sce

basel.cellType <- basel$cell_type
colData(basel) <- NULL
basel$cell_type <- basel.cellType
basel$cohort <- "Basel"


zurich1 <- assignIdentity(zurich1_cells, zurich1_types, zurich1_states)$sce

zurich1.cellType <- zurich1$cell_type
colData(zurich1) <- NULL
zurich1$cell_type <- zurich1.cellType
zurich1$cohort <- "Zurich1"

wagner <- assignIdentity(wagner_cells, wagner_types, wagner_states)$sce

wagner.cellType <- wagner$cell_type
colData(wagner) <- NULL
wagner$cell_type <- wagner.cellType
wagner$cohort <- "Wagner"

### [RENAME MARKERS] #####
wagner.names <- wagner %>% rownames()
wagner.names[6] <- "Cytokeratin 8/18"
wagner.names[13] <- "Her2"
wagner.names[28] <- "E-Cadherin"
wagner.names[29] <- "Ki-67"
rownames(wagner) <- wagner.names

### [COMBINE & START CLUSTERING] #####
clustering.markers <- Reduce(intersect, list(basel %>% rownames(),
                                             #schapiro %>% rownames(),
                                             zurich1 %>% rownames(),
                                             wagner %>% rownames()))

sce <- cbind(basel[clustering.markers, ],
             #schapiro[clustering.markers, ],
             wagner[clustering.markers, ],
             zurich1[clustering.markers, ])

expression <- logcounts(sce) %>% 
  t() %>% 
  as.data.frame()

expression$cell_type <- sce[, rownames(expression)]$cell_type
expression$cohort <- as.factor(sce[, rownames(expression)]$cohort)

#expression <- expression %>% as.matrix()

tsne <- fftRtsne(expression[, 1:length(clustering.markers)] %>% as.matrix()) %>% 
  as.data.frame()
colnames(tsne) <- c("tSNE Dim1", "tSNE Dim2")

tsne$cell_type <- expression$cell_type
tsne$cohort <- expression$cohort

library(devtools)
devtools::load_all("~/taproom/")

### [SAVE PLOT AS PDF] #####
# save as png
png(paste0(output_dir, "tSNE_cellType_unknownCells.png"), width = 650)
ggplot(tsne, aes(x = `tSNE Dim1`, y = `tSNE Dim2`)) +
  geom_point(aes(color = cell_type)) +
  scale_color_manual(values = jackson_basel_colours(), name = "Cell types") +
  astir_paper_theme()
dev.off()
