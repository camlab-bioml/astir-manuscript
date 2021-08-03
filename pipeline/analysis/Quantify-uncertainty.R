library(SingleCellExperiment)
library(tidyverse)
library(devtools)
library(ggpubr)

devtools::load_all("../taproom/")


# Read in data
zurich <- readRDS("output/phoenix/sces/zurich1_sce.rds")
zurich_assignments <- read_csv("output/phoenix/astir_assignments/zurich1_astir_assignments.csv")

# Get the maximum probability for each cell
max_prob <- zurich_assignments %>% 
  select(-Other) %>% 
  column_to_rownames("X1") %>% 
  as.matrix() %>% 
  rowMax()

max_prob_df <- tibble(cell = zurich_assignments$X1, 
                      max_prob = max_prob)


# Select Ay16
Ay16 <- zurich[,grepl("Ay16", colnames(zurich))]

# Get the probabilites associated with cells on Ay16
Ay16_prob <- max_prob_df[grepl("Ay16", max_prob_df$cell),]
names <-  Ay16_prob %>% 
  pull(cell)

Ay16_prob$core <- sapply(strsplit(names, "_"), function(x) x[4])

get_p <- function(x){
  p <- t.test(x1_4, x, alternative = "less")$p.value
  
  if(p == 0 | p < 2.2e-16){
    p <- "p < 2.2e-16"
  } else{
    p <- signif(p, 3)
  }
  
  p
}

x1_4 <- filter(Ay16_prob, grepl("Ay16x1|Ay16x2|Ay16x3|Ay16x4", core)) %>% 
  pull(max_prob)
x_5 <- filter(Ay16_prob, core == "Ay16x5") %>% 
  pull(max_prob)
x_6 <- filter(Ay16_prob, core == "Ay16x6") %>% 
  pull(max_prob)
x_7 <- filter(Ay16_prob, core == "Ay16x7") %>% 
  pull(max_prob)
x_8 <- filter(Ay16_prob, core == "Ay16x8") %>% 
  pull(max_prob)

x5.p <- get_p(x_5)
x6.p <- get_p(x_6)
x7.p <- get_p(x_7)
x8.p <- get_p(x_8)

# Plot probability
pdf("output/phoenix/figures/Staining_Astir_probability_Ay16.pdf", width = 7, height = 5.2)
Ay16_prob %>% 
  ggplot(aes(x = core, y = max_prob, fill = core)) +
  geom_boxplot() +
  # x1-4
  annotate("segment", x = "Ay16x1", xend = "Ay16x4", y = 0.5, yend = 0.5) +
  # x5
  annotate("segment", x = 2.5, xend = 2.5, y = 0.5, yend = 0.12) +
  
  annotate("segment", x = 2.5, xend = "Ay16x5", y = 0.33, yend = 0.33) +
  annotate("segment", x = "Ay16x5", xend = "Ay16x5", y = 0.33, yend = 0.35) +
  annotate("text", x = 3.75, y = 0.36, label = x5.p, hjust = 0.5) +
  
  # x6
  annotate("segment", x = 2.5, xend = "Ay16x6", y = 0.26, yend = 0.26) +
  annotate("segment", x = "Ay16x6", xend = "Ay16x6", y = 0.26, yend = 0.28) +
  annotate("text", x = 4, y = 0.29, label = x6.p, hjust = 0.5) +
  
  # x7
  annotate("segment", x = 2.5, xend = "Ay16x7", y = 0.19, yend = 0.19) +
  annotate("segment", x = "Ay16x7", xend = "Ay16x7", y = 0.19, yend = 0.21) +
  annotate("text", x = 4.75, y = 0.22, label = x7.p, hjust = 0.5) +
  
  # x8
  annotate("segment", x = 2.5, xend = "Ay16x8", y = 0.12, yend = 0.12) +
  annotate("segment", x = "Ay16x8", xend = "Ay16x8", y = 0.12, yend = 0.14) +
  annotate("text", x = 5.25, y = 0.15, label = x8.p, hjust = 0.5) +
  
  ylim(0,1) +
  labs(x = "Sample", y = "Maximum probability") +
  scale_fill_brewer(palette = "Blues") +
  astir_paper_theme() +
  theme(legend.position = "None",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# Get raw imc values (row sums over all channels)
raw_imc <- assays(Ay16)$raw_imc %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("cell")

raw_imc$total_raw_signal <- select(raw_imc, -cell) %>% rowSums()
raw_imc$core <- sapply(strsplit(raw_imc$cell, "_"), function(x) x[4])



# Plot raw imc values
pdf("output/phoenix/figures/Staining_raw_signal_Ay16.pdf", width = 7, height = 5.2)

x1_4 <- filter(raw_imc, grepl("Ay16x1|Ay16x2|Ay16x3|Ay16x4", core)) %>% 
  pull(total_raw_signal)
x_5 <- filter(raw_imc, core == "Ay16x5") %>% 
  pull(total_raw_signal)
x_6 <- filter(raw_imc, core == "Ay16x6") %>% 
  pull(total_raw_signal)
x_7 <- filter(raw_imc, core == "Ay16x7") %>% 
  pull(total_raw_signal)
x_8 <- filter(raw_imc, core == "Ay16x8") %>% 
  pull(total_raw_signal)


x5.p <- get_p(x_5)
x6.p <- get_p(x_6)
x7.p <- get_p(x_7)
x8.p <- get_p(x_8)


raw_imc %>% 
  ggplot(aes(x = core, y = total_raw_signal, fill = core)) +
  geom_boxplot() +
  # x1-4
  annotate("segment", x = "Ay16x1", xend = "Ay16x4", y = 100, yend = 100) +

  # # x5
  annotate("segment", x = 2.5, xend = 2.5, y = 100, yend = 460) +

  annotate("segment", x = 2.5, xend = "Ay16x5", y = 250, yend = 250) +
  annotate("segment", x = "Ay16x5", xend = "Ay16x5", y = 250, yend = 240) +
  annotate("text", x = 3.75, y = 265, label = x5.p, hjust = 0.5) +
 
  # # x6
  annotate("segment", x = 2.5, xend = "Ay16x6", y = 320, yend = 320) +
  annotate("segment", x = "Ay16x6", xend = "Ay16x6", y = 320, yend = 310) +
  annotate("text", x = 4, y = 335, label = x6.p, hjust = 0.5) +

  # # x7
  annotate("segment", x = 2.5, xend = "Ay16x7", y = 390, yend = 390) +
  annotate("segment", x = "Ay16x7", xend = "Ay16x7", y = 390, yend = 380) +
  annotate("text", x = 4.75, y = 405, label = x7.p, hjust = 0.5) +

  # # x8
  annotate("segment", x = 2.5, xend = "Ay16x8", y = 460, yend = 460) +
  annotate("segment", x = "Ay16x8", xend = "Ay16x8", y = 460, yend = 450) +
  annotate("text", x = 5.25, y = 475, label = x8.p, hjust = 0.5) +
  

  labs(x = "Sample", y = "Total raw signal") +
  scale_fill_brewer(palette = "Blues") +
  astir_paper_theme() +
  theme(legend.position = "None",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))



dev.off()




### Within patient comparison
Ay <- zurich[,grepl("Ay16x1|Ay16x2|Ay15x7|Ay15x8", colnames(zurich))]

raw_imc <- assays(Ay)$raw_imc %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("cell")

raw_imc$total_raw_signal <- select(raw_imc, -cell) %>% rowSums()
raw_imc$core <- sapply(strsplit(raw_imc$cell, "_"), function(x) x[4])


pdf("output/phoenix/figures/Staining_Within_patient_comparison.pdf", heigh = 5.2, width = 4)
raw_imc %>% 
  ggplot(aes(x = core, y = total_raw_signal, fill = core)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("Ay15x7", "Ay15x8"),
                                        c("Ay15x7", "Ay16x1"),
                                        c("Ay15x7", "Ay16x2"),
                                        c("Ay15x8", "Ay16x1"),
                                        c("Ay15x8", "Ay16x2"),
                                        c("Ay16x1", "Ay16x2")),
                     method.args = list(alternative = "greater"),
                     tip.length = 0) +
  labs(x = "Sample", y = "Total raw signal") +
  scale_fill_brewer(palette = "Blues") +
  astir_paper_theme() +
  theme(legend.position = "None",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# Get the probabilites associated with cells on Ay16
Ay16_prob <- max_prob_df[grepl("Ay16x1|Ay16x2|Ay15x7|Ay15x8", max_prob_df$cell),]
names <-  Ay16_prob %>% 
  pull(cell)

Ay16_prob$core <- sapply(strsplit(names, "_"), function(x) x[4])

pdf("output/phoenix/figures/Astir_max_prob_within_patient_comparison.pdf", width = 4, height = 5.2)
Ay16_prob %>% 
  ggplot(aes(x = core, y = max_prob, fill = core)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("Ay15x7", "Ay15x8"),
                                        c("Ay15x7", "Ay16x1"),
                                        c("Ay15x7", "Ay16x2"),
                                        c("Ay15x8", "Ay16x1"),
                                        c("Ay15x8", "Ay16x2"),
                                        c("Ay16x1", "Ay16x2")),
                     method.args = list(alternative = "greater"),
                     label.y = c(0.1, 0.17, 0.24, 0.31, 0.38, 0.45),
                     tip.length = 0) +
  ylim(0,1) +
  labs(x = "Sample", y = "Maximum probability") +
  scale_fill_brewer(palette = "Blues") +
  astir_paper_theme() +
  theme(legend.position = "None",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()