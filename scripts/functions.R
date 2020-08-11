combineRDS <- function(metadata, path){
  # Testing combining rds files
  # First, select all the cores that are non-tumour for subsequent exclusion
  nonTumourCores <- metadata %>% 
    filter(diseasestatus == "non-tumor") %>% 
    pull(core)
  
  # List all RDS files
  RDSfiles <- dir(path, pattern = ".rds", full.names = T)
  
  RDSfiles <- RDSfiles[-grep("Liver", RDSfiles)] # remove liver samples
  
  # Remove non tumour cores
  RDSfiles <- RDSfiles[-grep(paste(nonTumourCores, collapse="|"), RDSfiles)]
  
  listSCE <- lapply(RDSfiles, readRDS)
  sce <- do.call('cbind', listSCE)
  
  sce
}

createSCE <- function(path){
  RDSfiles <- dir(path, pattern = ".rds", full.names = T)
  
  listSCE <- lapply(RDSfiles, readRDS)
  sce <- do.call('cbind', listSCE)
  
  sce
}

assignIdentity <- function(raw.sce, types, states, dimReduct = F){
  #### REad in data
  # raw.sce <- "output/v4/zurich1_subset/zurich1_subset_sce.rds"
  # types <- "output/v4/zurich1_subset/zurich1_subset_assignments_type.csv"
  # states <- "output/v4/zurich1_subset/zurich1_subset_assignments_state.csv"
  
  if(is.character(raw.sce)){
    sce <- readRDS(raw.sce)
    
    ### Make sure wagner has id field
    if(any(colnames(colData(sce)) == "id") == F){
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
  
  assignments <- taproom::get_celltypes(types_mat) %>% as.data.frame()
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


create.alluvial <- function(df, method) {
  # Ceate alluvial plotting function
  # Define baseline plot
  plot <- ggplot(df, aes(y = df[, 3], axis = df[, 1], axis2 = df[, 2])) +
    geom_alluvium(aes(fill = df[, 1])) +
    geom_stratum(width = 1/20, fill = "grey", color = "black") +
    ylab("Cells") +
    labs(fill = "Cell Types") +
    scale_x_discrete(limits = c("Astir cell type", paste(method, "cluster")), 
                     expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    astir_paper_theme() +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.background = element_blank())
  
  if (length(unique(df[, 1])) < 13) {
    # For few starting states
    plot <- plot + scale_fill_brewer(type = "qual", palette = "Set1") +
      scale_fill_manual(values = jackson_basel_colours())
  }else {
    # There are too many categories in the starting column
    plot <- plot +
      theme(legend.position = "none")
  }
  
  plot
}
