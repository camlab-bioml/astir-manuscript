#!/usr/local/bin/Rscript

library(tidyverse)
library(devtools)
devtools::load_all("~/taproom/")

### [READ IN ARGUMENTS] #####
args <- commandArgs(trailingOnly = TRUE)

state <- args[1] %>% read_csv()
type <- args[2] %>% read_csv()
cohort <- args[3]
output_dir <- args[4]

### [PROCESS DATA] #####
cells <- type %>% 
  mutate(cell_type = get_celltypes(select(., -X1))) %>% 
  select(X1, cell_type) %>% 
  rename(core = X1) %>% 
  mutate(core = sub("_[^_]+$", "", core))

coreFreqs <- cells %>% 
  group_by(core) %>% 
  count(cell_type)

# create plotting order
coreOrder <- coreFreqs %>% 
  filter(grepl("Epithelial", cell_type)) %>% 
  group_by(core) %>% 
  summarise(total = sum(n)) %>% 
  arrange(desc(total)) %>% 
  pull(core)

# change plotting order
coreFreqs$core <- factor(coreFreqs$core, levels = coreOrder)
  
# Plot
pdf(paste0(output_dir, "core_cell_identity_barchart_", cohort, ".pdf"))
ggplot(coreFreqs, aes(x = core, y = n, fill = cell_type)) + 
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = jackson_basel_colours()) +
  astir_paper_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank())
dev.off()
