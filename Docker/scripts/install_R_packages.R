packages <- c(
    'tidyverse',
    'devtools',
    'cowplot',
    'ggalluvial',
    'circlize',
    'broom',
    'viridis',
    'rsvd',
    'wesanderson',
    'Hmisc',
    'stringr',
    'reshape2',
    'yaml',
    'pheatmap',
    'factoextra',
    'dichromat',
    'ggpubr',
    'SingleCellExperiment',
    'scater',
    'ComplexHeatmap',
    'GSVA',
    'DelayedArray'
)

ip <- installed.packages()
packages <- packages[!(packages %in% rownames(ip))]

BiocManager::install(packages)