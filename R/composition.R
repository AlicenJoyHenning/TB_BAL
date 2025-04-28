# Script for analysis sample compositions
library(dplyr)

# Retrieve data
integrated <- readRDS("~/Projects/TB_BAL_data/data/integrated/annotated.rds")

info <- integrated@meta.data
info$sample <- sub(".*_", "", info$orig.ident)
info$origin_celltype <- paste0(info$orig.ident, "_", info$annotation)
info$sample_celltype <- paste0(info$sample, "_", info$annotation)

write.csv(info, 
          './data/annotated_object.csv',
          quote = FALSE)


# Count cells per sample and annotation
celltype_counts_normalized <- info %>%
  count(sample, annotation) %>%
  left_join(info %>% count(sample, name = "total_cells"), by = "sample") %>%
  mutate(proportion = n / total_cells)

# Subset for BAL and PBMCs 
BAL <- subset(info, grepl("BAL", info$orig.ident))
PBMC <- subset(info, grepl("PBMC", info$orig.ident))

# Correct normalization
BAL_normalized <- BAL %>%
  count(sample, annotation) %>%
  left_join(BAL %>% count(sample, name = "total_cells"), by = "sample") %>%
  mutate(proportion = n / total_cells)
BAL_normalized$total_cells <- NULL
BAL_normalized$n <- NULL
BAL_normalized

PBMC_normalized <- PBMC %>%
  count(sample, annotation) %>%
  left_join(PBMC %>% count(sample, name = "total_cells"), by = "sample") %>%
  mutate(proportion = n / total_cells)
PBMC_normalized$total_cells <- NULL
PBMC_normalized$n <- NULL
PBMC_normalized


