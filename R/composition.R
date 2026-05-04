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
# BAL <- subset(info, grepl("BAL", info$orig.ident))
PBMC <- subset(info, grepl("PBMC", info$orig.ident))
BAL <- readRDS( "~/Projects/TB_BAL_data/data/integrated/TB_BAL_lymphocytes_labelled.rds")
# Save for plot 

# Correct normalization
BAL_normalized <- BAL@meta.data %>%
  count(orig.ident, groups) %>%
  mutate(proportion = n / sum(n))

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
             "myeloid" = "#1D78B4",
             "Myeloid" = "#1D78B4"
            )


# Prepare the data
celltype_counts_normalized <- celltype_counts_normalized %>%
  mutate(annotation = factor(annotation, levels = c("B", "T", "ncMono", "cMono", "myeloid"))) %>%
  group_by(sample) %>%
  mutate(proportion = n / sum(n))  # Normalize proportions


celltype_counts_normalized <- BAL@meta.data %>%

  mutate(groups = factor(groups, levels = c("B", "T",  "Myeloid"))) %>%
  group_by(orig.ident) %>%
  count(groups) %>%
  mutate(proportion = n / sum(n))

# Create the stacked bar plot
celltype_probs <- BAL@meta.data %>%
  mutate(groups = factor(groups, levels = c("B", "T", "Myeloid"))) %>%
  group_by(orig.ident) %>%
  count(groups) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# Create the stacked bar plot
plot <- ggplot(celltype_probs, 
       aes(x = proportion, y = orig.ident, fill = groups)) +
  geom_bar(stat = "identity") +
  # Adds percentage labels on the axes
  scale_x_continuous(labels = scales::percent_format()) + 
  scale_fill_manual(values = colours) +
  labs(fill = "Cell Group", x = "Proportion of Sample", y = "Sample ID") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title.y = element_blank() # Matches your previous style
  )

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
