library(tidyverse)

input_dir <- "../../output/v1/benchmarking/geneset/"

csv_files <- dir(input_dir, full.names = TRUE)


df <- map_dfr(csv_files, read_csv)
print(length(csv_files))

df <- group_by(df, method, n_cells) %>% 
  summarise(elapsed = mean(elapsed),
            mean_correlation = mean(mean_correlation)) %>% 
  ungroup()

ggplot(df, aes(x = n_cells, y = elapsed, colour = method)) + 
  geom_line() +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

ggplot(df, aes(x = n_cells, y = abs(mean_correlation), colour = method)) + 
  geom_line() +
  geom_point() +
  scale_x_log10()
