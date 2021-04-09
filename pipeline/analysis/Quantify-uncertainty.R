library(SingleCellExperiment)
library(tidyverse)
library(devtools)
library(ggpubr)

devtools::load_all("../taproom/")

zurich <- readRDS("output/squirrel/sces/zurich1_sce.rds")
zurich_assignments <- read_csv("output/squirrel/astir_assignments/zurich1_astir_assignments.csv")
zurich_meta <- read_csv("data-raw/jackson-2020/SingleCell_and_Metadata/ZurichTMA/Zuri_PatientMetadata.csv")

# List all cores used for supplementary figure
selected_cores <- c("SP43_137_X9Y3", "slide_37_Cy1x8", "slide_66_Cy4x1", 
                    "SP41_126_X14Y7", "SP42_185_X7Y5", "slide_33_Cy7x3", 
                    "SP43_102_X1Y4", "SP43_270_X6Y8", "slide_21_Cy5x8", 
                    "SP43_172_X9Y2", "SP42_189_X5Y5", "slide_25_By15x1",
                    "SP42_59_X3Y9", "SP42_64_X14Y4", "SP42_171_X2Y7", 
                    "SP41_191_X15Y7", "slide_52_Cy3x1", "SP42_127_X2Y8", 
                    "SP41_133_X3Y6", "slide_35_By15x5")

# Get the maximum probability for each cell
max_prob <- zurich_assignments %>% 
  select(-Other) %>% 
  column_to_rownames("X1") %>% 
  as.matrix() %>% 
  rowMax()

max_prob_df <- tibble(cell = zurich_assignments$X1, 
                      max_prob = max_prob)

# Get all cells on cores of interest & merge with prob & slide size
plotted <- zurich[,grepl(paste(selected_cores, collapse = "|"), colnames(zurich))] %>% 
  colData() %>% 
  as.data.frame() %>% 
  select(Location_Center_X, Location_Center_Y) %>% 
  rownames_to_column("cell") %>% 
  left_join(max_prob_df) %>% 
  mutate(core = sub("_[^_]+$", "", cell)) %>% 
  left_join(select(zurich_meta, core, Height_FullStack, Width_FullStack))

all_zurich <- zurich %>% 
  colData %>% 
  as.data.frame() %>% 
  select(Location_Center_X, Location_Center_Y) %>% 
  rownames_to_column("cell") %>% 
  left_join(max_prob_df) %>% 
  mutate(core = sub("_[^_]+$", "", cell)) %>% 
  left_join(select(zurich_meta, core, Height_FullStack, Width_FullStack))


assign_margin_cells <- function(df){
  df %>% 
    group_by(core) %>% 
    mutate(core_area = Height_FullStack * Width_FullStack) %>% 
    mutate(inner_height = sqrt((9 * core_area * Height_FullStack) / (10 * Width_FullStack))) %>% 
    mutate(inner_width = sqrt((9 * core_area * Width_FullStack) / (10 * Height_FullStack))) %>% 
    mutate(x_margin = (Width_FullStack - inner_width) / 2) %>% 
    mutate(y_margin = (Height_FullStack - inner_height) / 2) %>% 
    mutate(location = ifelse(Location_Center_X > x_margin & 
                               Location_Center_X  < (inner_width + x_margin) &
                               Location_Center_Y > y_margin &
                               Location_Center_Y < (inner_height + y_margin), "inner", "outer"))
}

plot_max_prob_violin <- function(df){
  df %>% 
    ggplot(aes(x = location, y = max_prob)) +
    geom_violin() +
    labs(y = "Maximum probability assigned to each cell") +
    astir_paper_theme() + 
    stat_compare_means(method = "t.test", label.x = 1.4, label.y = 1.01)
}

classified_cells <- assign_margin_cells(plotted)
all_cells <- assign_margin_cells(all_zurich)
  
plot_max_prob_violin(classified_cells)
plot_max_prob_violin(all_cells)
  



rand_sample <- unique(all_zurich$core)[1:12]


all_cells %>% 
  filter(core %in% rand_sample) %>%
  ggplot(aes(x = Location_Center_X, y = Location_Center_Y, color = location)) +
  geom_point()+
  facet_wrap(~core, scales = "free")


cells_per_core <- all_cells %>% 
  group_by(core) %>% 
  tally()

summary(cells_per_core$n)

over_500_cores <- cells_per_core %>% 
  filter(n > 500) %>% 
  pull(core)

all_cells %>% 
  filter(core %in% over_500_cores) %>% 
  plot_max_prob_violin()





average_max_prob_per_core <- all_cells %>% 
  group_by(core) %>% 
  summarize(mean = mean(max_prob))




low_prob_cores <- average_max_prob_per_core %>% 
  filter(mean < mean(.$mean)) %>% 
  pull(core)

all_cells %>% 
  filter(core %in% low_prob_cores) %>% 
  plot_max_prob_violin()





ratio <- all_cells %>% 
  group_by(core, location) %>% 
  summarize(mean = mean(max_prob)) %>% 
  pivot_wider(names_from = "location", values_from = "mean") %>% 
  mutate(inner_outer_ration = inner/outer) %>% 
  left_join(cells_per_core)


ratio %>% 
  ggplot(aes(x = n, y = inner_outer_ration)) +
  geom_point() +
  #labs(y = "Mean probability assigned to each core") +
  astir_paper_theme()
