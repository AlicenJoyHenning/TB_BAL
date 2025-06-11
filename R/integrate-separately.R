library(Seurat)
library(ggplot2)

# PBMC ----
data_dir <- "~/Projects/TB_BAL_data/data/preprocessed/"
PBMC_1098 <- readRDS(paste0(data_dir, "PBMC_1098.rds"))
PBMC_1376 <- readRDS(paste0(data_dir, "PBMC_1376.rds"))
PBMC_1483 <- readRDS(paste0(data_dir, "PBMC_1483.rds"))
PBMC_1484 <- readRDS(paste0(data_dir, "PBMC_1484.rds"))
PBMC_1523 <- readRDS(paste0(data_dir, "PBMC_1523.rds"))
PBMC_1566 <- readRDS(paste0(data_dir, "PBMC_1566.rds"))
PBMC_1676 <- readRDS(paste0(data_dir, "PBMC_1676.rds"))

rds_list <- list(PBMC_1098, PBMC_1376, PBMC_1483, PBMC_1484, 
                  PBMC_1523, PBMC_1566, PBMC_1676)

cells <- c(
  unique(rds_list[[1]]$orig.ident), unique(rds_list[[2]]$orig.ident), 
  unique(rds_list[[3]]$orig.ident), unique(rds_list[[4]]$orig.ident), 
  unique(rds_list[[5]]$orig.ident), unique(rds_list[[6]]$orig.ident), 
  unique(rds_list[[7]]$orig.ident)
)
  

# Merge data
TB_PBMC_dataset <- merge(
  x = rds_list[[1]], 
  y = c(
    rds_list[[2]], rds_list[[3]], rds_list[[4]],rds_list[[5]], 
    rds_list[[6]], rds_list[[7]]
  ),
  add.cell.ids = cells
)

# Integrate (Ensure v5)
class(TB_PBMC_dataset[['RNA']]) 
packageVersion("Seurat")

# Normalise counts 
TB_PBMC_dataset <- NormalizeData(TB_PBMC_dataset) %>%
  FindVariableFeatures(selection.method = "mean.var.plot") %>% 
  ScaleData() %>% 
  RunPCA()

TB_PBMC_dataset <- IntegrateLayers(
  TB_PBMC_dataset, 
  method = CCAIntegration,
  assay = "RNA",
  orig.reduction = "pca", 
  k.weight = 66,
  new.reduction = "integrated")

TB_PBMC_dataset[["RNA"]] <- JoinLayers(TB_PBMC_dataset[["RNA"]])

TB_PBMC_dataset <- FindNeighbors(TB_PBMC_dataset, dims = 1:10) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:10)


# BAL ----
data_dir <- "~/Projects/TB_BAL_data/data/preprocessed/"
BAL_1098 <- readRDS(paste0(data_dir, "BAL_1098.rds"))
BAL_1523 <- readRDS(paste0(data_dir, "BAL_1523.rds"))

rds_list <- list(BAL_1098, BAL_1523)

cells <- c(
  unique(rds_list[[1]]$orig.ident), 
  unique(rds_list[[2]]$orig.ident))


# Merge data
TB_BAL_dataset <- merge(
  x = rds_list[[1]], 
  y = c(rds_list[[2]]),
  add.cell.ids = cells
)

# Integrate (Ensure v5)
class(TB_BAL_dataset[['RNA']]) 
packageVersion("Seurat")

# Normalise counts 
TB_BAL_dataset <- NormalizeData(TB_BAL_dataset) %>%
  FindVariableFeatures(selection.method = "mean.var.plot") %>% 
  ScaleData() %>% 
  RunPCA()

TB_BAL_dataset <- IntegrateLayers(
  TB_BAL_dataset, 
  method = CCAIntegration,
  assay = "RNA",
  orig.reduction = "pca", 
  k.weight = 66,
  new.reduction = "integrated")

TB_BAL_dataset[["RNA"]] <- JoinLayers(TB_BAL_dataset[["RNA"]])

TB_BAL_dataset <- FindNeighbors(TB_BAL_dataset, dims = 1:10) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:10)


# Transfer annotations ----
integrated <- readRDS("~/Projects/TB_BAL_data/data/integrated/annotated.rds")
annotations <- data.frame(
  annotation = integrated$annotation)

TB_PBMC_dataset$annotations <- annotations$annotation[match(
  rownames(TB_PBMC_dataset@meta.data), rownames(annotations)
)]

TB_BAL_dataset$annotations <- annotations$annotation[match(
  rownames(TB_BAL_dataset@meta.data), rownames(annotations)
)]


# Plot ----
colours <- c("T" = "#A6CEE3", 
             "B" = "#B2DE89",
             "ncMono" = "#FB9A99",
             "cMono" = "#33A02C", 
             "myeloid" = "#1D78B4")


plot <- DimPlot(TB_BAL_dataset, 
        pt.size = 1,
        cols = colours,
        reduction = "umap",
        group.by = "annotations") +
  NoAxes() +
  labs(title = "Integrated BAL samples") + 
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1))

ggsave(plot = plot, 
       filename = "~/Projects/TB_BAL/plots/integrated/BAL_UMAP.png",
       width = 5,         
       height = 4,         
       dpi = 300 
)

features <- FeaturePlot(TB_BAL_dataset, 
                        features = c("MS4A1", "CD3E", "NKG7", 
                                     "MARCO", "FCGR3A", "CD14"), 
                        cols = c("#E1E1E1", "#1D78B4"))

# Remove axes and legends from all plots
features <- features & 
  theme_void() &  
  theme(
    legend.position = "none", 
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8)
    ) 

BAL <- plot | features 

ggsave(plot = BAL, 
       filename = "~/Projects/TB_BAL/plots/integrated/BAL_markers.png",
       width = 8,         
       height = 6,         
       dpi = 300 
)


# PBMC ---

plot <- DimPlot(TB_PBMC_dataset, 
                pt.size = 1,
                cols = colours,
                reduction = "umap",
                group.by = "annotations") +
  NoAxes() +
  labs(title = "Integrated BAL samples") + 
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1))

ggsave(plot = plot, 
       filename = "~/Projects/TB_BAL/plots/integrated/PBMC_UMAP.png",
       width = 5,         
       height = 4,         
       dpi = 300 
)

features <- FeaturePlot(TB_PBMC_dataset, 
                        features = c("MS4A1", "CD3E", "NKG7", 
                                     "MARCO", "FCGR3A", "CD14"), 
                        cols = c("#E1E1E1", "#1D78B4"))

# Remove axes and legends from all plots
features <- features & 
  theme_void() &  
  theme(
    legend.position = "none", 
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.8)
  ) 

PBMC <- plot | features 

ggsave(plot = PBMC, 
       filename = "~/Projects/TB_BAL/plots/integrated/PBMC_markers.png",
       width = 8,         
       height = 6,         
       dpi = 300 
)
