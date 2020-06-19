library(tidyverse)
library(cowplot)

input_dir <- "../../output/v4/benchmarking/geneset/"

csv_files <- dir(input_dir, full.names = TRUE)



df <- map_dfr(csv_files, read_csv)

xs <- sort(unique(df$n_cells))

dfg <- group_by(df, method, dataset, n_cells) %>% 
  summarise(elapsed = mean(elapsed),
            mean_correlation = mean(mean_correlation)) %>% 
  ungroup()

plt1 <- ggplot(dfg, aes(x = n_cells, y = elapsed, colour = method)) + 
  geom_line() +
  geom_point() +
  scale_x_log10(breaks=xs) +
  scale_y_log10() +
  facet_wrap(~ dataset)

plt2 <- ggplot(dfg, aes(x = n_cells, y = abs(mean_correlation), colour = method)) + 
  geom_line() +
  geom_point() +
  scale_x_log10(breaks=xs) +
  facet_wrap(~ dataset)

plot_grid(
  plt1,
  plt2,
  ncol = 1
)

