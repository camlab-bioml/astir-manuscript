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

if(cohort == "lin_cycif"){
  thresh <- 0.7
}else{
  thresh <- 0.7
}

# read in data 
sce <- readRDS(cells)
type <- read_csv(type)

type$cell_type <- get_celltypes(select(type, -X1), thresh) 
type <- select(type, X1, cell_type) %>%
  column_to_rownames("X1")
colData(sce)["cell_type"] <- type[colnames(sce),]

#sce <- assignIdentity(cells, type, state)$sce

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
  protein.order <- c("SMA", "FAP", "Vimentin", "CD45", "CD3", "CD68", "CD24", "CD31", 
                     "CD49f", "ECadherin", "EpCAM", "panK", "K14", "K5", "K7", "K8K18")
  correct_names <- c("SMA", "FAP", "Vimentin", "CD45", "CD3", "CD68", "CD24", "CD31", 
                     "CD49f", "E-Cadherin", "EpCAM", "pan Cytokeratin", "Cytokeratin 14", 
                     "Cytokeratin 5", "Cytokeratin 7", "Cytokeratin 8/18")
  title <- "Wagner"
}else if(cohort == "schapiro"){
  protein.order <- c("CD68", "E-cadherin", "EpCAM", "Cytokeratin7", "Cytokeratin8-18",
                     "Vimentin", "Fibronectin")
  correct_names <- c("CD68", "E-cadherin", "EpCAM", "Cytokeratin 7", "Cytokeratin 8/18",
                     "Vimentin", "Fibronectin")
  title <- "Schapiro"
}else if(cohort == "lin_cycif"){
  protein.order <- c("Vimentin", "E_Cadherin", "Keratin")
  correct_names <- c("Vimentin", "E Cadherin", "Keratin")

  title <- "Lin"
}

unknown.mat <- unknown.sce[, protein.order]
assigned.mat <- assigned.sce[, protein.order]

if(cohort == "wagner" | cohort == "schapiro" | cohort == "lin-cycif"){
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

pdf(paste0(output_dir, "expression_correlation_", cohort, ".pdf"), height = 10, width = 9)
correlations %>% 
  ggplot(aes(as.numeric(Var2), as.numeric(Var1), fill = r)) +
  geom_tile(color = "white") +
  ggtitle(title) +
  scale_fill_gradient2(low="#6D9EC1",mid="white",high="#E46726", 
                       midpoint = 0, limits=c(-1,1), name = "Pearson\nCorrelation  ") +
  labs(y = "Unknown and Other cell types", x = "Assigned cell types") +
  astir_paper_theme() + 
  scale_x_continuous(position = "bottom", expand=c(0,0), sec.axis = dup_axis(),
                     breaks = 1:length(x_labels), labels = x_labels) +
  scale_y_continuous(position = "left", expand=c(0,0), sec.axis = dup_axis(),
                     breaks = 1:length(y_labels), labels = y_labels) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 20, hjust = 1),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size=20),
        plot.title = element_text(size=26),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.x.bottom = element_blank(), # remove titles
        axis.title.y.left = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1,"cm")) +
  coord_fixed()
dev.off()

