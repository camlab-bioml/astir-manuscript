
library(tidyverse)
library(glue)
library(forcats)


devtools::load_all("../taproom")

schapiro_alt_mask_samples <- c('Cy1x5_32', 'Cy1x6_33', 'Cy1x8_35')
schapiro_users <- c('Catena', 'Jackson', 'Schulz')

dir_astir <- snakemake@params[['dir_astir']]
dir_other <- snakemake@params[['dir_other']]

get_df_for_sample <- function(sample) {
  path_for_sample <- glue("{dir_astir}/assignments_{sample}_{{user}}.csv")
  
  read_user <- function(user, path) {
    input_path <- glue(path)
    df <- read_csv(input_path)
    cell_ids <- df$X1
    df$X1 <- NULL
    celltypes <- taproom::get_celltypes(df, 0.5)
    tibble(
      cell_id = cell_ids,
      cell_type = celltypes,
      user = user
    )
  }
  
  df_all <- map_dfr(schapiro_users, read_user, path_for_sample)
  df_all$sample <- sample
  df_all
}

df <- map_dfr(schapiro_alt_mask_samples, get_df_for_sample)

ggplot(df, aes(x = cell_type, fill = user)) +
  geom_bar(position="dodge") +
  facet_wrap(~ sample)

# stop("Done")


# Investigate Cy1x6_33 ----------------------------------------------------
# 
# expression_dir <- "output/squirrel/schapiro_processed_alt_mask/"
# files <- dir(expression_dir, full.names=TRUE, pattern="Cy1x6_33")
# files <- files[grepl("csv", files)]
# 
# df_catena <- read_csv(files[1]) %>% 
#   mutate(user = "Catena")
# 
# df_jackson <- read_csv(files[2]) %>% 
#   mutate(user = "Jackson") 
# 
# df_exprs <- bind_rows(df_catena, df_jackson)
#   
# 
# df_tidy <- select(df_exprs, -X1) %>% 
#   gather(protein, expression, -user)
# 
# ggplot(df_tidy, aes(x = expression, fill = user)) +
#   geom_histogram() +
#   facet_wrap(~ protein, scales= "free")
# 
# 
# user <- "Catena"
# jackson_33_assignments <- read_csv(glue("output/squirrel/schapiro_astir_assignments_alt_mask/assignments_Cy1x6_33_{user}.csv"))
# cell_ids <- jackson_33_assignments$X1
# jackson_33_assignments$X1 <- NULL
# celltypes <- taproom::get_celltypes(jackson_33_assignments, 0.5)
# names(celltypes) <- cell_ids
# sce <- readRDS(glue("output/squirrel/schapiro_processed_alt_mask/Cy1x6_33_{user}.rds"))
# sce$cell_type <- celltypes[colnames(sce)]
# 
# createHeatmap(sce)

# stop("Done")

# Clustering part ---------------------------------------------------------

relevant_output <- dir(dir_other,full.names = TRUE)

clustering_outputs <- map_dfr(relevant_output, read_csv)
clustering_outputs$user <- sapply(strsplit(clustering_outputs$origin,' - '), `[`, 2)
clustering_outputs$sample <- sapply(strsplit(clustering_outputs$origin,' - '), `[`, 1)

clustering_outputs$method2 <- paste0(clustering_outputs$method, "-", clustering_outputs$params)

ggplot(clustering_outputs, aes(x = `cell type`, fill = user)) +
  geom_bar(position="dodge") +
  facet_grid(sample ~ method2)


# Comparison --------------------------------------------------------------

co <- dplyr::count(clustering_outputs, method2, user, sample, `cell type`) %>% 
  dplyr::rename(method = method2, cell_type = `cell type`)
ast <- dplyr::count(df, cell_type, user, sample) %>% 
  mutate(method = "Astir") %>% 
  filter(!grepl("Other|Unknown", cell_type))

df_all <- rbind(
  co,
  ast[,names(co)]
)

df_all <- spread(df_all, cell_type, n, fill = 0) %>% 
  # mutate(`Epithelial (all)` = `Epithelial (luminal)` + `Epithelial (other)`) %>% 
  gather(cell_type, n, -method, -user, -sample)

df_comp <- group_by(df_all, method, sample, cell_type) %>% 
  summarize(sd_n = sd(n), mean_n = median(n)) %>% 
  ungroup()

ggplot(df_comp, aes(x = fct_reorder(method, sd_n / mean_n), y = sd_n)) +
  geom_boxplot() +
  # facet_wrap(~ cell_type, scales = "free_y") +
  coord_flip()

ggplot(df_comp, aes(x = mean_n, y = sd_n, colour = method)) +
  geom_point()

df_freq <- group_by(df_all, method, user, sample) %>% 
  mutate(freq = n / sum(n)) %>% 
  ungroup()

df_comp2 <- group_by(df_freq, method, sample, cell_type) %>% 
  summarize(sd_n = sd(freq), mean_n = median(freq)) %>% 
  ungroup()

ggplot(df_comp2, aes(x = fct_reorder(method, sd_n), y = sd_n)) +
  geom_boxplot() +
  # facet_wrap(~ cell_type, scales = "free_y") +
  coord_flip()

df_rank <- group_by(df_comp2, sample, cell_type) %>% 
  mutate(rank = rank(sd_n)) %>% 
  ungroup() %>% 
  mutate(situation = paste0(sample, " - ", cell_type))

label_method <- function(method) {
  case_when(
    grepl("Astir", method) ~ "Astir",
    grepl("Phenograph", method) ~ "Phenograph",
    grepl("FlowSOM", method) ~ "FlowSOM",
    grepl("ClusterX", method) ~ "ClusterX"
  )
}

df_rank <- mutate(df_rank, overall_method = label_method(method))

df_rank_mean <- group_by(df_rank, overall_method, cell_type, sample) %>% 
  summarize(mean_rank = mean(rank)) %>% 
  ungroup()

ggplot(df_rank_mean, aes(x = cell_type, y = mean_rank, colour = overall_method, group = overall_method)) +
  geom_line(data = df_rank, 
            aes(x = cell_type, y = rank, colour = overall_method, group = method),
            alpha = 0.3,
            size = 1) +
  geom_line(size = 1.5) +
  geom_point(size = 2) +
  astir_paper_theme() +
  facet_wrap(~ sample, scales = "free_y") +
  scale_colour_brewer(palette = "Set1", name = "Method") +
  labs(x = "Cell type", y = "Method rank\n(s.d. of cell type proportions)") +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_x_discrete(expand = c(0.1,0.1))

pdf(snakemake@output[['pdf']], width = 10, height = 4.5)
print(last_plot())
dev.off()
