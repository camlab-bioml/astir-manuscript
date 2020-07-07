
suppressPackageStartupMessages({
  library(tidyverse)
  library(GSVA)
  library(argparse)
  library(yaml)
})


parser <- ArgumentParser(description = "Benchmark GSVA timings")

parser$add_argument('--input_dir', type = 'character')
parser$add_argument('--markers', type = 'character')
parser$add_argument('--n_cells', type = 'integer')
parser$add_argument('--method', type = 'character')
parser$add_argument('--dataset', type = 'character')
parser$add_argument('--output_file', type = 'character')
parser$add_argument('--max_epochs', type='integer')
parser$add_argument('--batch_size', type='integer')
parser$add_argument('--learning_rate', type='double')
parser$add_argument('--n_initial_epochs', type='integer')

args <- parser$parse_args()

markers <- read_yaml(args$markers)

csv_files <- dir(args$input_dir, pattern='csv', full.names = TRUE)
csv_files <- sample(csv_files) # randomize


expr_csv <- read_csv(csv_files[1])
i = 2

while(i < length(csv_files) && nrow(expr_csv) < args$n_cells) {
  expr_csv <- rbind(
    expr_csv,
    read_csv(csv_files[i])
  )
  i <- i+1
}

expr_csv <- sample_n(expr_csv, min(args$n_cells, nrow(expr_csv)))

g <- st <- expr_mat <- NULL

if("X" %in% names(expr_csv)) {
  expr_csv <- select(expr_csv, -X)
}

if("Y" %in% names(expr_csv)) {
  expr_csv <- select(expr_csv, -Y)
}

expr_mat <- select(expr_csv, -X1) %>% 
  as.matrix() %>% 
  t()

pathways_to_test <- markers$cell_states
pathways_to_test <- pathways_to_test[sapply(pathways_to_test, length) > 2]

heldout_proteins <- lapply(
  pathways_to_test,
  sample,
  size = 1
)


modified_markers <- lapply(names(pathways_to_test), function(pathway) {
  setdiff(
    pathways_to_test[[pathway]],
    heldout_proteins[[pathway]]
  )
})


# don't test proteins if there's any overlap
to_remove <- sapply(as.vector(heldout_proteins), function(hp) {
  hp %in% unlist(modified_markers)
})

heldout_proteins <- heldout_proteins[!to_remove]

lens <- sapply(markers$cell_types, length)
if(any(lens == 1)) {
  for(i in which(lens == 01)) {
    markers$cell_types[[i]] <- list(markers$cell_types[[i]])
  }
}

names(modified_markers) <- names(pathways_to_test)

if(args$method != "astir") {
  st <- system.time({
    g <- gsva(expr_mat, 
              modified_markers,
              method = args$method)
  })
} else {
  csv_file <- tempfile(tmpdir='tmp')
  output_file <- tempfile(tmpdir='tmp')
  expr_csv %>% 
    write_csv(csv_file)
  
  yaml_file <- tempfile(tmpdir='tmp')
  markers$cell_states <- modified_markers
  write_yaml(markers, yaml_file)

  cmd <- paste("/home/ltri/campbell/kcampbel/.conda/envs/imc/bin/python",
              "/home/ltri/campbell/kcampbel/.conda/envs/imc/bin/astir",
              "state",
              csv_file,
              yaml_file,
              output_file,
              "--max_epochs", args$max_epochs,
              "--batch_size", args$batch_size,
              "--learning_rate", args$learning_rate,
              "--n_init_epochs", args$n_initial_epochs,
              "--n_init 20"
  )


  st <- system.time({
    system(cmd)
  })
  
  output <- read_csv(output_file)
  output <- select(output, -X1)
  g <- t(as.matrix(output))
  
} 

get_stat <- function(pathway) {
# p <- names(modified_markers)[1]
  expression <- expr_mat[heldout_proteins[[pathway]],  ]
  pathway_score <- g[pathway,]
  
  fit <- lm(expression ~ pathway_score)
  s <- summary(fit)
  
  
  list(
    statistic = s$coefficients[2,3],
    coefficient = s$coefficients[2,1],
    correlation = cor(expression, pathway_score)
  )
}
stats <- lapply(
  names(heldout_proteins),
  get_stat
)

statistics <- sapply(stats, `[[`, 'statistic')
coefs <- sapply(stats, `[[`, 'coefficient')
correlations <- sapply(stats, `[[`, 'correlation')

save.image("deleteme.rdata")

directions <- sapply(names(heldout_proteins), function(mm) {
  proteins <- modified_markers[[mm]]
  apply(expr_mat[proteins,], 1, function(y) {
    cor(y, g[mm,])
  }) %>% 
    mean()
})

df_res <- tibble(
  method = args$method,
  dataset = args$dataset,
  n_cells = args$n_cells,
  statistics = statistics,
  coefficients = coefs,
  correlations = correlations,
  directions = sign(directions)
)

# 
# cors <- sapply(names(markers$cell_states), function(p) {
#   c1 <- t(g[p,,drop=FALSE])
#   c2 <- t(expr_mat[markers$cell_states[[ p ]], ])
#   (cor(c1,c2))
# }, simplify = FALSE)
# 
# cors <- unlist(cors)



write_csv(df_res, args$output_file)
