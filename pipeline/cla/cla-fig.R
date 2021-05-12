
library(tidyverse)
library(forcats)
library(taproom)
library(ggbeeswarm)

# We need:
# (1) output figs tables directory
# (2) 3x output PDF files

output_fig_dir <- snakemake@params[['output_fig_dir']]

astir_paper_theme <- function() {
  cowplot::theme_cowplot(font_size = 12) +
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(face="bold"),
          legend.text=element_text(size=12),
          legend.title=element_text(face = "bold", size=14))
}

method_cols <- c(
  "Supervised"="#AA3939",
  "Unsupervised"="#882D61",
  "Cluster & interpret"="#AA6C39"
)

# Cell level accuracy fig -------------------------------------------------

cla_files <- dir(output_fig_dir, pattern="annotation*.*tsv", full.names=TRUE)

df_cla <- map_dfr(cla_files, read_tsv)

df_cla$method <- gsub("Manual", "z_score", df_cla$method, fixed=TRUE)

df_cla <- filter(df_cla, !grepl("high", method))

df_cla <- mutate(df_cla, method = case_when(
    method == "acdc-no-consider" ~ "ACDC_no_consider",
    method == "acdc-absent" ~ "ACDC_absent",
    TRUE ~ method
  )) %>% 
  mutate(method_type = case_when(
  grepl("LDA", method) ~ "Supervised",
  grepl("Astir|ACDC", method) ~ "Unsupervised",
  TRUE ~ "Cluster & interpret"
))

df_cla <- mutate(df_cla, .metric = case_when(
  .metric == "Spec" ~ "Specificity",
  .metric == "Sens" ~ "Sensitivity",
  .metric == "Mcc" ~ "MCC",
  TRUE ~ .metric
))

df_cla <- mutate(df_cla, cohort = str_to_title(cohort))

df_cla$cohort <- gsub("1", "", df_cla$cohort)

df_cla$method <- gsub("_cell_type", "", df_cla$method)
df_cla$method <- gsub(" default", "", df_cla$method)

df_cla <- drop_na(df_cla, .estimate)

tmp <- group_by(df_cla, method) %>% 
  summarize(mean_est = mean(.estimate))
method_ordering <- arrange(tmp, mean_est) %>% .$method

df_cla_sum <- group_by(df_cla, .metric, method, method_type, cohort) %>% 
  summarize(mean_value = mean(.estimate), min_value = max(0, min(.estimate))) %>% 
  ungroup() %>% 
  group_by(.metric) %>% 
  mutate(min_value = min(min_value)) %>% 
  ungroup()

# ggplot(df_cla, aes(x = factor(method,  levels = method_ordering), y = .estimate, fill = method_type)) +
#   # geom_bar(stat='summary', alpha = 0.8) +
#   # geom_boxplot() +
#   geom_segment(data = df_cla_sum, aes(x = factor(method,  levels = method_ordering), xend=method, y=min_value, yend=mean_value, colour=method_type),
#                size=5) +
#   geom_quasirandom(shape=21, size=2.5) + 
#   facet_grid(.metric ~ cohort, scales = "free_y")  +
#   astir_paper_theme() +
#   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
#   labs(x = "Method", y = "Estimate",
#        title = "Comparison to manually annotated cell types") +
#   scale_fill_manual(values=method_cols, name = "Method type") +
#   scale_colour_manual(values=method_cols, name = "Method type") +
#   theme(legend.title = element_text(size=10))+
#   stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
#                geom = "crossbar", width = 0.9) +
#   theme(legend.position = "top",
#         panel.background = element_rect(fill='grey95'))


ggplot(df_cla, aes(x = factor(method,  levels = method_ordering), y = .estimate, fill = method_type)) +
  geom_bar(stat='summary', alpha = 0.8) +
  # geom_boxplot() +
  # geom_segment(data = df_cla_sum, aes(x = factor(method,  levels = method_ordering), xend=method, y=min_value, yend=mean_value, colour=method_type),
               # size=5) +
  geom_quasirandom(shape=21, size=2.5) + 
  facet_grid(.metric ~ cohort, scales = "free_y")  +
  astir_paper_theme() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  labs(x = "Method", y = "Estimate",
       title = "Comparison to manually annotated cell types") +
  scale_fill_manual(values=method_cols, name = "Method type") +
  scale_colour_manual(values=method_cols, name = "Method type") +
  theme(legend.title = element_text(size=10))+
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width = 0.9) +
  theme(legend.position = "top",
        panel.background = element_rect(fill='grey95'))

ggsave(snakemake@output[['annotation']], width=8, height = 8)

# df_cla2 <- filter(df_cla, method_type != "Supervised", !grepl("high", method))
# 
# df_cla2 <- group_by(df_cla2, cohort, .metric, annotator_test) %>% 
#   mutate(is_best = .estimate == max(.estimate, na.rm=TRUE),
#          ratio_from_best = max(.estimate, na.rm=TRUE) / .estimate) %>% 
#   ungroup()
# 
# ggplot(df_cla2, aes(x = method, y = ratio_from_best)) +
#   geom_boxplot() +
#   coord_flip()
# 
# df_cla_count <- dplyr::filter(df_cla2, is_best == TRUE) %>% 
#   dplyr::count(method)
# 
# ggplot(df_cla_count, aes(x = method, y = n)) +
#   geom_bar(stat = 'identity')







# Cluster level accuracy --------------------------------------------------

cluster_files <- dir(output_fig_dir, pattern="cluster*.*tsv", full.names=TRUE)
cluster_files <- c(cluster_files, file.path(output_fig_dir, "wagner/cla_cluster_wagner.tsv"))
cluster_files <- c(cluster_files, file.path(output_fig_dir, "lin-cycif/cla_cluster_lin-cycif.tsv"))

df_clus <- map_dfr(cluster_files, read_tsv)

df_clus$method <- gsub("Manual", "z_score", df_clus$method, fixed=TRUE)

df_clus <- filter(df_clus, !grepl("high", method))

df_clus <- mutate(df_clus, method = case_when(
  method == "acdc-no-consider" ~ "ACDC_no_consider",
  method == "acdc-absent" ~ "ACDC_absent",
  TRUE ~ method
)) %>% 
  mutate(method_type = case_when(
    grepl("LDA", method) ~ "Supervised",
    grepl("Astir|ACDC", method) ~ "Unsupervised",
    TRUE ~ "Cluster & interpret"
  ))

df_clus <- mutate(df_clus, cohort = str_to_title(cohort))

df_clus$cohort <- gsub("1", "", df_clus$cohort)

df_clus$.metric[df_clus$.metric == "Bal_accuracy"] <- "Balanced\naccuracy"

# df_clus$method <- paste0(df_clus$method, "-", df_clus$params)

df_clus$method <- gsub("_cell_type", "", df_clus$method)
df_clus$method <- gsub(" default", "", df_clus$method)

df_clus <- drop_na(df_clus, .estimate)

tmp <- group_by(df_clus, method) %>% 
  summarize(mean_est = mean(.estimate))
method_ordering <- arrange(tmp, mean_est) %>% .$method

df_clus_sum <- group_by(df_clus, .metric, method, method_type, cohort) %>% 
  summarize(mean_value = mean(.estimate), min_value = max(0, min(.estimate))) %>% 
  ungroup() %>% 
  group_by(.metric) %>% 
  mutate(min_value = min(min_value)) %>% 
  ungroup()

# ggplot(df_clus, aes(x = factor(method,  levels = method_ordering), y = .estimate, fill = method_type)) +
#   geom_bar(stat='summary', alpha = 0.8) +
#   # geom_boxplot() +
#   # geom_segment(data = df_clus_sum, aes(x = factor(method,  levels = method_ordering), xend=method, y=min_value, yend=mean_value, colour=method_type),
#                # size=5) +
#   geom_quasirandom(shape=21, size=2.5) + 
#   facet_grid(.metric ~ cohort, scales = "free_y")  +
#   astir_paper_theme() +
#   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
#   labs(x = "Method", y = "Estimate",
#        title = "Comparison to annotated clusters") +
#   scale_fill_manual(values=method_cols, name = "Method type") +
#   scale_colour_manual(values=method_cols, name = "Method type") +
#   theme(legend.title = element_text(size=10))+
#   stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
#                geom = "crossbar", width = 0.9) +
#   theme(legend.position = "top",
#         panel.background = element_rect(fill='grey95'))



ggplot(df_clus, aes(x = factor(method,  levels = method_ordering), y = .estimate, fill = method_type)) +
  geom_bar(stat='summary', alpha = 0.8) +
  # geom_boxplot() +
  # geom_segment(data = df_clus_sum, aes(x = factor(method,  levels = method_ordering), xend=method, y=min_value, yend=mean_value, colour=method_type),
  #              size=5) +
  geom_quasirandom(shape=21, size=2.5) + 
  facet_grid(.metric ~ cohort, scales = "free_y")  +
  astir_paper_theme() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  labs(x = "Method", y = "Estimate",
       title = "Comparison to annotated clusters") +
  scale_fill_manual(values=method_cols, name = "Method type") +
  scale_colour_manual(values=method_cols, name = "Method type") +
  theme(legend.title = element_text(size=10))+
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width = 0.9) +
  theme(legend.position = "top",
        panel.background = element_rect(fill='grey95'))

ggsave(snakemake@output[['clustering']], width=8, height = 8)



# Debug phenograph --------------------------------------------------------

df_pg <- filter(df_cla, grepl("Phenograph", method), grepl("z_score", method))

df_pg$k <- as.numeric(sapply(strsplit(df_pg$params, "_"), `[`, 4))
df_pg$markers <- sapply(strsplit(df_pg$params, "_"), function(s) paste0(s[1:2], collapse="_"))


ggplot(df_pg, aes(x = k, y = .estimate, colour=markers, linetype=annotator_test)) +
  geom_point() +
  geom_line() +
  facet_grid(.metric ~ cohort, scales = "free_y") +
  scale_x_log10() +
  astir_paper_theme()

ggsave(snakemake@output[['phenograph']], width=10,height=5)

