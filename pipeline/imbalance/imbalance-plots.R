

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

count <- dplyr::count
rename <- dplyr::rename

devtools::load_all("../taproom")


# Astir -------------------------------------------------------------------



files <- dir(snakemake@params[['input_dir_astir']], pattern = "astir_imbalance", full.names = TRUE)

read_astir <- function(f) {
  p_epithelial <- gsub(".csv", "", strsplit(f, "_")[[1]][3], fixed=TRUE)
  df <- read_csv(f)
  df$X1 <- NULL
  celltypes <- get_celltypes(df)
  df2 <- tibble(celltype = celltypes)
  df_count <- count(df2, celltype)
  df_count <- mutate(
    df_count,
    n_cluster = 1 * (n>0),
    prop_celltype = n / sum(n)
  )
  df_count <- 
    mutate(df_count,
           p_epithelial = p_epithelial,
           method = "Astir",
           workflow = "Astir")
  df_count
}

df_astir <- map_dfr(files, read_astir)


# Everything else ---------------------------------------------------------


files <- dir(snakemake@params[['input_dir_other']],
             pattern = "^Epithelial",
             full.names = TRUE)

read_other <- function(f) {
  df <- read_csv(f)
  names(df)[2] <- "cluster"
  df <- rename(df, 
               p_epithelial = percent_epithelial,
               celltype = `cell type`)
  df <- mutate(df, workflow = paste0(method, " - ", params))
  
  df_count <- count(df, method, workflow, p_epithelial, celltype, cluster) %>% 
    count(method, workflow, p_epithelial, celltype) %>% 
    rename(n_cluster = n) 
  
  df_count2 <- count(df, method, workflow, p_epithelial, celltype)
  
  df_count <- inner_join(df_count, df_count2,
                         by = c("method", "workflow", "p_epithelial", "celltype"))
  
  df_count <- mutate(
    df_count,
    prop_celltype = n / sum(n)
  )
  df_count
}

df_other <- map_dfr(files, read_other)

df_astir$p_epithelial <- as.numeric(df_astir$p_epithelial)
df <- bind_rows(df_astir, df_other)

df_mean <- group_by(df, celltype, p_epithelial, method) %>% 
  summarize(mean_prop_celltype = mean(prop_celltype),
            mean_n_cluster = mean(n_cluster)) %>% 
  ungroup() %>% 
  filter(grepl("Epithelial", celltype))

filter(df, grepl("Epithelial", celltype)) %>% 
  ggplot(aes(x = 100 * p_epithelial, colour = method)) +
  geom_line(aes(group = workflow, y = 100 * prop_celltype), size = 1.2, alpha = 0.3) +
  geom_point(data = df_mean, aes(y = 100 * mean_prop_celltype)) +
  geom_line(data = df_mean, aes(y = 100 * mean_prop_celltype), size = 1.5) +
  facet_wrap(~ method) +
  astir_paper_theme() +
  scale_colour_brewer(palette = "Set1", guide = FALSE) +
  labs(
    x = "% cells epithelial in sample",
    y = "% cells assigned to epithelial"
  )

n_cell_plot <- last_plot()
  

filter(df, grepl("Epithelial", celltype)) %>% 
  ggplot(aes(x = 100 * p_epithelial, colour = method)) +
  # geom_point() + 
  geom_line(aes(group = workflow, y = n_cluster), size = 1.2, alpha = 0.3) +
  geom_point(data = df_mean, aes(y = mean_n_cluster)) +
  geom_line(data = df_mean, aes(y = mean_n_cluster), size = 1.5) +
  facet_wrap(~ method) +
  scale_y_log10()+
  astir_paper_theme() +
  scale_colour_brewer(palette = "Set1", guide = FALSE) +
  labs(
    x = "% cells epithelial in sample",
    y = "# clusters corresponding to epithelial"
  )

prop_cell_plot <- last_plot()

plot_grid(
  n_cell_plot, 
  prop_cell_plot,
  ncol = 1
)

save.image("tpm.rds")

pdf(snakemake@output[['pdf']], width=4,height=5.2)
print(last_plot())
dev.off()
