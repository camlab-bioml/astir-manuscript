library(tidyverse)
library(glue)
library(devtools)
library(DelayedArray)
library(ggpubr)
library(RColorBrewer)
library(colorspace)
devtools::load_all("../taproom/")

### List all the files
acdc_files <- dir("output/cardinal/results/alternative_masks/", full.names = TRUE, 
                  pattern = "ACDC")
other_files <- dir("output/cardinal/results/alternative_masks/", full.names = TRUE,
                   pattern = "Alternative_masks")
astir_files <- dir("output/cardinal/schapiro_astir_assignments_alt_mask/", full.names = TRUE,
                   pattern = "assignments")

### Read in ACDC & other
acdc <- map_dfr(acdc_files, read_tsv)
other <- map_dfr(other_files, read_csv)
other$seed[is.na(other$seed)] <- 0


### Define all necessary variables to read in astir & read in
schapiro_alt_mask_samples <- c('Cy1x5_32', 'Cy1x6_33', 'Cy1x8_35')
schapiro_users <- c('Catena', 'Jackson', 'Schulz')
iteration <- c('0', '1', '2', '3', '4')
dir_astir <- "output/cardinal/schapiro_astir_assignments_alt_mask/"

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
  pivot_wider(names_from = method, values_from = sd) %>% 
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

display.brewer.pal(n = 5, name = "Set1")
cols <- brewer.pal(n = 5, name = "Set1")[c(2,4,5,3)]
# Get plotting order
plot_order <- cell_proportion_meanSD %>% 
  ungroup() %>% 
  select(-c(cell_type, sample)) %>% 
  group_by(algorithm) %>% 
  summarise(mean_sd = mean(sd)) %>% 
  arrange(.$mean_sd) %>% pull(algorithm)


plot_order <- rowMeans(order) %>% sort()

cell_proportion_meanSD$algorithm <- factor(cell_proportion_meanSD$algorithm, levels = plot_order)


cell_proportion_meanSD %>% 
  ggplot(aes(x = algorithm, y = sd, color = algorithm, fill = algorithm)) +
  geom_boxplot(lwd = 1.5, fatten = 0.8) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = lighten(cols, amount = 0.9)) +
  labs(x = "Method", y = "Mean standard deviation in cell type\nproportion between segmentations") +
  astir_paper_theme() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none")




## Scatterplot
# 
# within_segmenter_sd <- seg_clean %>%
#   group_by(cell_type, sample, method, user, iteration) %>% 
#   tally() %>% 
#   ungroup() %>% 
#   group_by(cell_type, sample, method, user) %>% 
#   mutate(within_segmenter_sd = sd(n)) %>% 
#   mutate(norm_within_segmenter_sd = within_segmenter_sd / n)
#   
# seg_clean %>% 
#   group_by(cell_type, sample, method, user, iteration) %>% 
#   tally() %>% 
#   ungroup() %>% 
#   
# 
# left_join(select(within_segmenter_sd,
#                  -within_segmenter_sd), sd_data) %>% 
#   mutate(Algorithm = case_when(grepl("acdc", method) ~ "ACDC",
#                                grepl("FlowSOM", method) ~ "FlowSOM",
#                                TRUE ~ method)) %>% 
#   filter(!grepl("Phenograph|ClusterX", method)) %>% 
#   ggplot(aes(x = norm_within_segmenter_sd, y = norm_sd)) +
#   geom_point() +
#   labs(x = "Within segmenter standard deviation (normalized by cell count)",
#        y = "Across segmenter standard deviation\n(normalized by cell count)") + 
#   facet_wrap(~Algorithm) +
#   astir_paper_theme()
# 
# 
# 
# 
# 
# 
# df %>% 
#   group_by(sample, iteration, user, cell_type) %>% 
#   tally() %>% 
#   ggplot(aes(x = cell_type, y = n)) +
#   geom_bar(stat = 'identity')+
#   facet_grid(sample ~ user + iteration) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# df %>% 
#   pivot_wider(names_from = cell_type, values_from = iteration)
# 
# 
# 
# 
# 
# 
# 
# ### Average number of cells
# mean_cells <- seg %>% 
#   group_by(sample, method, user, iteration, cell_type) %>% 
#   tally() %>% 
#   ungroup() %>% 
#   group_by(sample, method, cell_type) %>% 
#   summarize(mean = mean(n)) %>% 
#   ungroup()
# 
# 
# mean_is_assigned <- seg %>% 
#   mutate(is_assigned = (!grepl("Unknown|unknown|Other", cell_type))) %>% 
#   group_by(sample, method, iteration, user) %>% 
#   summarize(is_assigned = count(is_assigned == TRUE)) %>% 
#   ungroup() %>% 
#   group_by(sample, method) %>% 
#   summarize(mean_is_assigned = mean(is_assigned))
# 
# left_join(mean_is_assigned, total_cells) %>% 
#   mutate(fraction_assigned = mean_is_assigned / mean_cell_no) %>% 
#   ggplot(aes(x = method, y = fraction_assigned)) +
#   geom_bar(stat = 'identity') +
#   facet_wrap(~sample) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
# 
# 
# seg %>% 
#   group_by(sample, method, user, iteration, cell_type) %>% 
#   tally() %>% 
#   filter(iteration == 0 & user == "Catena")
# 
# 
# 
# total_cells <- df %>% 
#   filter(iteration == 0) %>% 
#   group_by(sample, user) %>% 
#   tally() %>% 
#   ungroup() %>% 
#   group_by(sample) %>% 
#   mutate(mean_cell_no = mean(n)) %>% 
#   select(sample, mean_cell_no) %>% 
#   distinct()
# 
# 
# 
# mean_cells %>% 
#   filter(!grepl("Unknown|unknown|Other", cell_type)) %>% 
#   group_by(sample, method) %>% 
#   summarize(assigned_average_numeber = sum(mean)) %>% 
#   left_join(total_cells) %>% 
#   mutate(fraction_assigned_cells = assigned_average_numeber / mean_cell_no) %>% 
#   ggplot(aes(x = method, y = fraction_assigned_cells)) +
#   geom_bar(stat = 'identity') +
#   facet_wrap(~sample, ncol = 3) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
# 
# 
# 
# 
# 
# 
# 
# acdc_clean <- acdc %>% filter(cell_type != "unknown")
# 
# 
# acdc_clean %>% 
#   select(-c(cell_id, cohort)) %>% 
#   mutate(sample = paste(sample, " - ", method)) %>% 
#   group_by(sample, user, cell_type, iteration) %>% 
#   tally() %>% #View() 
#   ggplot(aes(x = cell_type, y = n, color = user)) +
#   geom_point() +
#   labs(y = "number of cells") +
#   #scale_fill_discrete() +
#   facet_wrap(~sample, ncol = 5) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  