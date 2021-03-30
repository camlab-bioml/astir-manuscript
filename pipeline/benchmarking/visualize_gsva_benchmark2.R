library(tidyverse)
library(cowplot)

theme_set(theme_cowplot(font_size=11))

input_dir <- "output/v4/benchmarking/geneset/"

csv_files <- dir(input_dir, full.names = TRUE)



df <- map_dfr(sample(csv_files), read_csv)

ggplot(df, aes(x = factor(n_cells), y = directions * (statistics), fill = method)) +
  facet_wrap(~ dataset) +
  geom_boxplot() +
  scale_y_log10()

# ggplot(df, aes(x = factor(n_cells), y = directions * (correlations), fill = method)) +
#   facet_wrap(~ dataset) +
#   geom_boxplot()



# plt2 <- ggplot(df, aes(x = factor(n_cells), y = abs(coefficients), fill = method)) +
#   facet_wrap(~ dataset) +
#   geom_boxplot() + 
#   ylim(-50, 100) + 
#   scale_y_log10()

# plot_grid(
#   plt1,
#   plt2,
#   ncol = 1
# )

stop("done")

# xs <- sort(unique(df$n_cells))
# 
# dfg <- group_by(df, method, dataset, n_cells) %>% 
#   summarise(elapsed = mean(elapsed),
#             mean_correlation = mean(mean_correlation)) %>% 
#   ungroup()
# 
# plt1 <- ggplot(dfg, aes(x = n_cells, y = elapsed, colour = method)) + 
#   geom_line() +
#   geom_point() +
#   scale_x_log10(breaks=xs) +
#   scale_y_log10() +
#   facet_wrap(~ dataset)
# 
# plt2 <- ggplot(dfg, aes(x = n_cells, y = abs(mean_correlation), colour = method)) + 
#   geom_line() +
#   geom_point() +
#   scale_x_log10(breaks=xs) +
#   facet_wrap(~ dataset)
# 
# plot_grid(
#   plt1,
#   plt2,
#   ncol = 1
# )

