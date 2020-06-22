### [PURPOSE OF THIS SCRIPT] #####
# This script reads in the original sc proteomics data and cell type and state 
# are then assigned. For cell states, in addition to the continuous value
# (between 0 and 1), an additional column is added classifying each cell as 
# either high or low for that state. This happens in a cell type specific manner.
# Finally, cell identity is a columnt that is composed of the cell type, followed
# by boolean classifications for each state, e.g: Macrophage_RTK-hi_proliferation-lo...

# This script returns a single cell experiment object and can be called by any 
# subsequent script.

assignIdentity2 <- function(raw.sce, types, states){
  #### REad in data
  # raw.sce <- "output/v4/zurich1_subset/zurich1_subset_sce.rds"
  # types <- "output/v4/zurich1_subset/zurich1_subset_assignments_type.csv"
  # states <- "output/v4/zurich1_subset/zurich1_subset_assignments_state.csv"
  sce <- readRDS(raw.sce)
  
  ### Make sure wagner has id field
  if(any(colnames(colData(sce)) == "id") == F & grepl("wagner", raw.sce)){
    id <- sce %>% colData() %>% as.data.frame() %>% rownames()
    sce$id <- id
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
  #   str_replace_all(" ", "_") %>% 
  #   str_replace_all("\\(", ".") %>% 
  #   str_replace_all("\\)", ".")
  # 
  # stateAbbrev <- str_replace(stateNames, "HALLMARK_", "")
  # stateMedians <- paste0(stateAbbrev, ".med")
  # stateFactors <- paste0(stateAbbrev, ".fac")
  # 
  # colnames(states_df) <- stateNames
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  states_df <- apply(states_df, 2, range01) %>% 
    as.data.frame()
  states_df$id <- rownames(states_df)
  
  
  typeState <- left_join(assignments, states_df) %>% 
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
  sce <- runPCA(sce, ncomponents = 10)
  sce <- runUMAP(sce)
  
  return(list(sce = sce, pathways = paste0(stateNames, ".activation")))
}



assignIdentity <- function(raw.sce, types, states){
  ### [DATA WRANGLING] #####
  # read in data
  #raw.sce <- "raw_data/wagner/wagner_subset_sce.rds"
  #types <- "raw_data/wagner/wagner_subset_assignments_type.csv"
  #states <- "raw_data/wagner/wagner_subset_assignments_state.csv"
  sce <- readRDS(raw.sce)
  
  if(any(colnames(colData(sce)) == "id") == F & grepl("wagner", raw.sce)){
    id <- sce %>% colData() %>% as.data.frame() %>% rownames()
    sce$id <- id
  }
  
  # ASSIGN CELL TYPES
  # ASSIGN CELL TYPES
  threshold <- 0.7
  df <- read_csv(types)
  
  assignments <- select(df, -X1) %>% as.matrix()
  
  # Add cell type assignments to df tbl
  df$cell_type <- taproom::get_celltypes(assignments)
  
  # Make sure the order of the cells is the same as in the sce object
  df <- df %>% 
    arrange(factor(X1, levels = sce$id))
  
  sce$cell_type <- df$cell_type
  
  
  # ASSIGN CONTINUOUS CELL STATE VALUES
  states <- read_csv(states) %>% 
    arrange(factor(X1, levels = sce$id)) %>% 
    dplyr::select(-X1) 
  
  
  #states_mat <- select(states, -X1) %>% 
  #  as.matrix()
  
  #states <- apply(states, 2, taproom:::winsorize_one, c(0.05, 0.95))
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  states <- apply(states, 2, range01)
  
  states <- as_tibble(states)
  
  stateNames <- states %>% # will become states_mat
    colnames() %>% 
    str_replace_all(" ", "_") %>% 
    str_replace_all("\\(", ".") %>% 
    str_replace_all("\\)", ".")

  stateAbbrev <- substr(stateNames, 1, 4)
  stateMedians <- paste0(stateAbbrev, ".med")
  stateFactors <- paste0(stateAbbrev, ".fac")
  
  #rownames(states_mat) <- states$X1
  #reducedDim(sce, 'states') <- states_mat[colnames(sce), ] # <- this will throw error
  
  #Add cell states
  sce.names <- sce %>% colData() %>% colnames()
  
  sce$state1 <- pull(states, 1)
  sce$state2 <- pull(states, 2)
  sce$state3 <- pull(states, 3)
  sce$state4 <- pull(states, 4)
  
  sce.names <- c(sce.names, stateNames)
  colnames(colData(sce)) <- sce.names
  
  # Create boolean classifications of cell states
  typeState <- sce %>% colData() %>% as.data.frame() %>% 
    select(cell_type, stateNames)
  
  # Calculate the median for each state (cell type specific)
  typeState <- typeState %>% dplyr::group_by(cell_type) %>% 
    dplyr::mutate(!!stateMedians[1] := median(.data[[stateNames[1]]]),
                 !!stateMedians[2] := median(.data[[stateNames[2]]]),
                 !!stateMedians[3] := median(.data[[stateNames[3]]]),
                 !!stateMedians[4] := median(.data[[stateNames[4]]])) %>% 
    ungroup()
  
  # Create factor variables for each state: hi vs lo
  typeState <- typeState %>% 
    dplyr::mutate(
      !!stateFactors[1] := factor(ifelse(.data[[stateNames[1]]] < 
                                           .data[[stateMedians[1]]], "lo", "hi")),
      !!stateFactors[2] := factor(ifelse(.data[[stateNames[2]]] < 
                                          .data[[stateMedians[2]]], "lo", "hi")),
      !!stateFactors[3] := factor(ifelse(.data[[stateNames[3]]] < 
                                          .data[[stateMedians[3]]], "lo", "hi")),
      !!stateFactors[4] := factor(ifelse(.data[[stateNames[4]]] < 
                                          .data[[stateMedians[4]]], "lo", "hi")))
  
  
  typeState <- typeState %>% 
    dplyr::mutate(type_state = paste(cell_type, 
                              stateAbbrev[1], .data[[stateFactors[1]]], 
                              stateAbbrev[2], .data[[stateFactors[2]]], 
                              stateAbbrev[3], .data[[stateFactors[3]]],
                              stateAbbrev[4], .data[[stateFactors[4]]], sep = "_"),
           state_summary = paste(stateAbbrev[1], .data[[stateFactors[1]]],
                                 stateAbbrev[2], .data[[stateFactors[2]]],
                                 stateAbbrev[3], .data[[stateFactors[3]]],
                                 stateAbbrev[4], .data[[stateFactors[4]]], sep = "_")) %>% 
    select(stateFactors, type_state, state_summary)
  
  # Add factor states & type state data to main sce object
  colData(sce) <- cbind(colData(sce), typeState)
  
  
  # REMOVE OTHER AND UNKNOWN CELL TYPES
  OtherCells <- which(sce$cell_type == "Other")
  UnknownCells <- which(sce$cell_type == "Unknown")
  
  sce <- sce[, -c(OtherCells, UnknownCells)]
  #colnames(colData(sce)) <- sce %>% colData() %>% as.data.frame() %>% colnames()
  
  
  ### [PCA & UMAP PROJECTIONS] #####
  sce <- runPCA(sce, ncomponents = 10)
  sce <- runUMAP(sce)
  
  #t <- list(sce = sce, pathways = stateNames)
  
  return(list(sce = sce, pathways = stateNames)) 
}



create.alluvial <- function(df, title) {
  # Ceate alluvial plotting function
  # Define baseline plot
  plot <- ggplot(df, aes(y = df[, 3], axis = df[, 1], axis2 = df[, 2])) +
    geom_alluvium(aes(fill = df[, 1])) +
    geom_stratum(width = 1/20, fill = "grey", color = "black") +
    ggtitle(title)
  
  
  if (length(unique(df[, 1])) < 13) {
    # For few starting states
    plot <- plot + scale_fill_brewer(type = "qual", palette = "Set1") +
      scale_fill_manual(values = jackson_basel_colours()) +
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            panel.background = element_blank(),
            axis.title.y = element_blank())
  }else {
    # There are too many categories in the starting column
    plot <- plot +
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            panel.background = element_blank(),
            legend.position = "none",
            axis.title.y = element_blank())
  }
  
  plot
}
