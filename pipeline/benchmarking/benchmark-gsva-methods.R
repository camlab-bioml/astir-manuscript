
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

expr_csv <- sample_n(expr_csv, args$n_cells)

g <- st <- expr_mat <- NULL

if(args$method == "astir") {
  csv_file <- tempfile()
  output_file <- tempfile()
  select(expr_csv, -X, -Y) %>% 
    write_csv(csv_file)
  
  expr_mat <- select(expr_csv, -X, -Y, -X1) %>% 
    as.matrix() %>% 
    t()
  
  cmd <- paste("astir state",
              csv_file,
              args$markers,
              output_file,
              "--max_epochs", args$max_epochs,
              "--batch_size", args$batch_size,
              "--learning_rate", args$learning_rate,
              "--n_init_epochs", args$n_initial_epochs
  )

  st <- system.time({
    system(cmd)
  })
  
  output <- read_csv(output_file)
  output <- select(output, -X1)
  g <- t(as.matrix(output))
  
} else {

  expr_mat <- select(expr_csv, -X1, -X, -Y) %>% 
    as.matrix() %>% 
    t()
  
  st <- system.time({
  g <- gsva(expr_mat, 
            markers$cell_states,
            method = args$method)
  })
}



cors <- sapply(names(markers$cell_states), function(p) {
  c1 <- t(g[p,,drop=FALSE])
  c2 <- t(expr_mat[markers$cell_states[[ p ]], ])
  (cor(c1,c2))
}, simplify = FALSE)

cors <- unlist(cors)

df_res <- tibble(
  method = args$method,
  dataset = args$dataset,
  n_cells = args$n_cells,
  elapsed = st['elapsed'],
  mean_correlation = mean(cors)
)

write_csv(df_res, args$output_file)