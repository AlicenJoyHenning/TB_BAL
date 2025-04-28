# Script for analysis sample compositions
library(dplyr)
library(ggplot2)

# Retrieve data ----
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

# Data manipulation -----
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

# Plot ----
colours <- c("T" = "#A6CEE3", 
             "B" = "#B2DE89",
             "ncMono" = "#FB9A99",
             "cMono" = "#33A02C", 
             "myeloid" = "#1D78B4")

# Prepare the data
celltype_counts_normalized <- celltype_counts_normalized %>%
  mutate(annotation = factor(annotation, levels = c("B", "T", "ncMono", "cMono", "myeloid"))) %>%
  group_by(sample) %>%
  mutate(proportion = n / sum(n))  # Normalize proportions

# Create the stacked bar plot
plot <- ggplot(celltype_counts_normalized, 
               aes(x = proportion, y = factor(sample), fill = annotation)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colours) +
  labs(fill = "") +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        legend.position = "bottom")  

# Save the plot as an SVG
ggsave("./plots/composition/PBMC_BAL_proportion.svg", plot, width = 8, height = 4)

# Repeat for isolated 

# Create the stacked bar plot
BAL_plot <- ggplot(BAL_normalized, 
               aes(x = proportion, y = factor(sample), fill = annotation)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colours) +
  labs(fill = "") +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        legend.position = "bottom")  

# Save the plot as an SVG
ggsave("./plots/composition/BAL_proportion.svg", BAL_plot, width = 8, height = 4)

PBMC_plot <- ggplot(PBMC_normalized, 
                   aes(x = proportion, y = factor(sample), fill = annotation)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colours) +
  labs(fill = "") +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        legend.position = "bottom")  

# Save the plot as an SVG
ggsave("./plots/composition/PBMC_proportion.svg", PBMC_plot, width = 8, height = 4)
