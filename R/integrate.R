# Integrate the individually pre-processed objects 
# https://satijalab.org/seurat/articles/seurat5_integration

library(Seurat)

# Retrieve data 
data_dir <- "~/Projects/TB_BAL_data/data/preprocessed/"
rds_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE)
rds_list <- setNames(
  lapply(rds_files, readRDS),
  tools::file_path_sans_ext(basename(rds_files))
)

cells <- c(
  unique(rds_list[[1]]$orig.ident), unique(rds_list[[2]]$orig.ident), 
  unique(rds_list[[3]]$orig.ident), unique(rds_list[[4]]$orig.ident), 
  unique(rds_list[[5]]$orig.ident), unique(rds_list[[6]]$orig.ident), 
  unique(rds_list[[7]]$orig.ident), unique(rds_list[[8]]$orig.ident), 
  unique(rds_list[[9]]$orig.ident)
)

# Merge data
TB_dataset <- merge(
  x = rds_list[[1]], 
  y = c(
    rds_list[[2]], rds_list[[3]], rds_list[[4]],rds_list[[5]], 
    rds_list[[6]], rds_list[[7]], rds_list[[8]], rds_list[[9]]
  ),
  add.cell.ids = cells
)

# Integrate (Ensure v5)
class(TB_dataset[['RNA']]) 
packageVersion("Seurat")

# Normalise counts 
TB_dataset <- NormalizeData(TB_dataset) %>%
  FindVariableFeatures(selection.method = "mean.var.plot") %>% 
  ScaleData() %>% 
  RunPCA()

TB_dataset <- IntegrateLayers(
  TB_dataset, 
  method = CCAIntegration,
  assay = "RNA",
  orig.reduction = "pca", 
  k.weight = 66,
  new.reduction = "integrated")

TB_dataset[["RNA"]] <- JoinLayers(TB_dataset[["RNA"]])

TB_dataset <- FindNeighbors(TB_dataset, dims = 1:10) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:10)

# Save & visualise
plot_markers(TB_dataset, "small-subset-integrated", 
             output = "~/Projects/TB_BAL_data/data/")
integrated <- readRDS("~/Projects/TB_BAL_data/data/Integrated.rds")

# Annotate 
DimPlot(integrated, reduction = "umap", label = TRUE)
integrated$annotation <- ifelse(integrated$seurat_clusters %in% c(4, 5, 11), 
                                "T", "-")
integrated$annotation <- ifelse(integrated$seurat_clusters %in% c(9, 17), 
                                "B", integrated$annotation)
integrated$annotation <- ifelse(integrated$seurat_clusters == "14", 
                                "ncMono", integrated$annotation)
integrated$annotation <- ifelse(integrated$seurat_clusters == "12", 
                                "cMono", integrated$annotation)

Idents(integrated) <- "annotation"
markers <- Seurat::FindMarkers(integrated,
                               ident.1 = "unknown",
                               ident.2 = NULL)

markers <- subset(markers, p_val_adj <= 0.05 & avg_log2FC >= 2)
markers <- markers[order(-markers$avg_log2FC), ]
top10_genes <- rownames(head(markers, 10))


DimPlot(integrated, reduction = "umap", group.by = "orig.ident")

FeaturePlot(integrated, reduction = "umap",
            features = c("FCGR3A", "CD14", "MARCO", "CD68"), 
            label = TRUE, repel = TRUE)  


saveRDS(integrated, "~/Projects/TB_BAL_data/data/integrated/annotated.rds")

colours <- RColorBrewer::brewer.pal(12, 'Paired')
colours_full <- c(colours, colours)
plot_clusters <- DimPlot(integrated,
                         pt.size = 1,
                         cols = colours_full) +
  NoAxes() +
  labs(title = "Integrated BAL and PBMC of TB dataset") + 
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1))

annotated_clusters <- DimPlot(integrated,
                         pt.size = 1,
                         cols = colours_full) +
  NoAxes() +
  labs(title = "Integrated BAL and PBMC samples") + 
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1))
ggsave(plot = annotated_clusters, 
       filename = "~/Projects/TB_BAL/plots/annotated.svg",
       width = 7,         
       height = 6,         
       dpi = 300 
)

plot_clusters <- DimPlot(integrated,
                         pt.size = 1,
                         group.by = "orig.ident",
                         cols = colours_full) +
  NoAxes() +
  labs(title = "Batch effect from low quality BAL samples") + 
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1))
ggsave(plot = plot_clusters, 
       filename = "~/Projects/TB_BAL/plots/integrated/batch_effect.png",
       width = 7,         
       height = 6,         
       dpi = 300 
)
