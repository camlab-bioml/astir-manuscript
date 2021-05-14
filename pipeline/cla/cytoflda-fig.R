library(tidyverse)

df <- read_tsv(snakemake@input[[1]])

df <- filter(df, grepl("LDA", method))

ggplot(df, aes(x = annotator_train, y = .estimate, colour=annotator_test)) +
  geom_line(aes(group=annotator_test)) +
  geom_point() +
  facet_grid(.metric ~ cohort, scales = "free_y") +
  cowplot::theme_cowplot() +
  theme(legend.title = element_blank(),
        legend.position = "top")+
  labs(x = "Train annotator", y = "Estimate") +
  scale_colour_brewer(palette = "Dark2")

ggsave(snakemake@output[[1]], width = 10, height = 6)
