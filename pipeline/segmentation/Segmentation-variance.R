library(tidyverse)
library(glue)
library(devtools)
library(DelayedArray)
library(ggpubr)
library(RColorBrewer)
library(colorspace)
library(taproom)
source("scripts/functions.R")

### List all the files
acdc_files <- dir("output/phoenix/results/alternative_masks/", #snakemake@params[['dir_other']], 
                  full.names = TRUE, 
                  pattern = "ACDC")
other_files <- dir("output/phoenix/results/alternative_masks/", #snakemake@params[['dir_other']], 
                   full.names = TRUE,
                   pattern = "Alternative_masks")
astir_files <- dir("output/phoenix/schapiro_astir_assignments_alt_mask/", #snakemake@params[['dir_astir']], 
                   full.names = TRUE,
                   pattern = "assignments")

### Read in ACDC & other
acdc <- map_dfr(acdc_files, read_tsv)
other <- map_dfr(other_files, read_csv)
other$seed[is.na(other$seed)] <- 0


### Define all necessary variables to read in astir & read in
schapiro_alt_mask_samples <- c('Cy1x5_32', 'Cy1x6_33', 'Cy1x8_35')
schapiro_users <- c('Catena', 'Jackson', 'Schulz')
iteration <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
dir_astir <- "output/phoenix/schapiro_astir_assignments_alt_mask/"#snakemake@params[['dir_astir']]

get_df_for_sample <- function(sample, iteration) {
  path_for_sample <- glue("{dir_astir}/assignments_{sample}_{{user}}_{iteration}.csv")
  
  read_user <- function(user, path) {
    input_path <- glue(path)
    df <- read_csv(input_path)
    cell_ids <- df$X1
    df$X1 <- NULL
    celltypes <- taproom::get_celltypes(df, 0.5)
    tibble(
      cell_id = cell_ids,
      cell_type = celltypes,
      user = user,
      iteration = iteration
    )
  }
  
  df_all <- map_dfr(schapiro_users, read_user, path_for_sample)
  df_all$sample <- sample
  df_all
}

df <- lapply(schapiro_alt_mask_samples, function(x){
  lapply(iteration, function(k){
    get_df_for_sample(x, k)
  }) %>% bind_rows()
}) %>% bind_rows()



### Some data cleanup
df$method <- "Astir"

acdc <- dplyr::rename(acdc, "iteration" = "seed")
acdc$iteration <- as.character(acdc$iteration)

other <- dplyr::rename(other, "cell_id" = "id", "cell_type" = "Manual_cell_type", "iteration" = "seed")
other$user <- sapply(strsplit(other$core,' - '), `[`, 2)
other$sample <- sapply(strsplit(other$core,' - '), `[`, 1)
other$method <- paste(other$method, other$params)
other$iteration <- as.character(other$iteration)



### Bind everything together
seg <- bind_rows(df,
          acdc[,c("cell_id", "cell_type", "user", "iteration", "sample", "method")],
          other[,c("cell_id", "cell_type", "user", "iteration", "sample", "method")]) %>% 
  filter(!grepl("k_500", method))

seg_clean <- filter(seg, !grepl("Unknown|Other|unknown", cell_type))


#### Start plotting

# Get number of cells (cell type specific)
cell_counts <- seg_clean %>% 
  group_by(method, sample, user, iteration, cell_type) %>% 
  tally()

# Get total number of cells per sample and segmenter
total_cell_counts <- seg_clean %>% 
  filter(method == "Astir" & iteration == 0) %>% 
  group_by(sample, user) %>% 
  tally() %>% 
  dplyr::rename("total_cells" = "n")
  
# Combine the two
cell_counts <- left_join(cell_counts, total_cell_counts)

# Get cell type proportion & calculate mean proportion across seeds
cell_proportion_mean <- cell_counts %>% 
  mutate(fraction = n / total_cells) %>% 
  select(-n, -total_cells) %>% 
  ungroup() %>% 
  group_by(method, sample, user, cell_type) %>% 
  summarize(mean_segmenter_prop = mean(fraction)) %>% 
  ungroup()

# Calculate standard deviations across users
cell_proportion_sd <- cell_proportion_mean %>% 
  group_by(method, sample, cell_type) %>% 
  summarize(sd = sd(mean_segmenter_prop)) %>% 
  drop_na(sd) %>% 
  ungroup()

# Remove any method that is NA
cell_proportion_sd_no_na <- cell_proportion_sd %>% 
  pivot_wider(names_from = method, values_from = sd) %>% View
  select_if(~ !any(is.na(.))) %>% 
  pivot_longer(-c(sample, cell_type), names_to = "method", values_to = "sd")

# Calculate average standard deviation across algorithms, samples and cell types
cell_proportion_meanSD <- cell_proportion_sd_no_na %>% 
  mutate(algorithm = case_when(grepl("acdc", method) ~ "ACDC",
                               grepl("ClusterX", method) ~ "ClusterX",
                               grepl("FlowSOM", method) ~ "FlowSOM",
                               grepl("Phenograph", method) ~ "Phenograph",
                               TRUE ~ method)) %>% 
  group_by(algorithm, sample, cell_type) %>% 
  summarize(mean_sd = mean(sd))


cols <- brewer.pal(n = 5, name = "Set1")[c(2,4,5,3)]
# Get plotting order
plot_order <- cell_proportion_meanSD %>% 
  ungroup() %>% 
  select(-c(cell_type, sample)) %>% 
  group_by(algorithm) %>% 
  summarise(mean_sd = mean(mean_sd)) %>% 
  arrange(.$mean_sd) %>% pull(algorithm)


# print("plot order")
# print(plot_order)
# plot_order <- rowMeans(plot_order) %>% sort()

cell_proportion_meanSD$algorithm <- factor(cell_proportion_meanSD$algorithm, 
                                           levels = plot_order)



pdf(snakemake@output[['pdf']], width = 2.5, height = 4.5)
cell_proportion_meanSD %>% 
  ggplot(aes(x = algorithm, y = mean_sd, color = algorithm, fill = algorithm)) +
  geom_boxplot(lwd = 1.5, fatten = 0.8) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = lighten(cols, amount = 0.9)) +
  labs(x = "Method", y = "Mean standard deviation in cell type\nproportion between segmentations") +
  astir_paper_theme() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none")
dev.off()




## Scatterplot
algorithmic_sd <- cell_counts %>% 
  mutate(fraction = n / total_cells) %>% 
  select(-n, -total_cells) %>% 
  ungroup() %>% 
  group_by(method, sample, user, cell_type) %>% 
  summarize(algorithmic_sd = sd(fraction)) %>% 
  ungroup()

clean_algorithmic_sd <- algorithmic_sd %>% 
  filter(!grepl("ClusterX|Pheno", method))


test <- cell_proportion_sd_no_na %>% 
  mutate(algorithm = case_when(grepl("acdc", method) ~ "ACDC",
                               grepl("ClusterX", method) ~ "ClusterX",
                               grepl("FlowSOM", method) ~ "FlowSOM",
                               grepl("Phenograph", method) ~ "Phenograph",
                               TRUE ~ method))


pdf(snakemake@output[['suppl_pdf']], width = 4, height = 7)
inner_join(clean_algorithmic_sd, 
           cell_proportion_sd_no_na) %>% 
  mutate(algorithm = case_when(grepl("acdc", method) ~ "ACDC",
                               grepl("ClusterX", method) ~ "ClusterX",
                               grepl("FlowSOM", method) ~ "FlowSOM",
                               grepl("Phenograph", method) ~ "Phenograph",
                               TRUE ~ method)) %>% 
  ggplot(aes(x = sd, y = algorithmic_sd)) +
  geom_point(size = 3) +
  labs(x = "Mean standard deviation in cell type\nproportion between segmentations",
       y = "Standard deviation in cell type\nproportions between seeds") +
  astir_paper_theme() +
  theme(strip.text = element_text(size = 18)) +
  facet_wrap(~algorithm)
dev.off()