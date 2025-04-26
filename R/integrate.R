# Integrate the individually pre-processed objects 
# https://satijalab.org/seurat/articles/seurat5_integration

library(Seurat)

# Retrieve data 
data_dir <- "./data/preprocessed"
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
  unique(rds_list[[9]]$orig.ident), unique(rds_list[[10]]$orig.ident), 
  unique(rds_list[[11]]$orig.ident), unique(rds_list[[12]]$orig.ident),
  unique(rds_list[[13]]$orig.ident)
)

# Merge data
TB_dataset <- merge(
  x = rds_list[[1]], 
  y = c(
    rds_list[[2]], rds_list[[3]], rds_list[[4]],rds_list[[5]], 
    rds_list[[6]], rds_list[[7]], rds_list[[8]], rds_list[[9]], 
    rds_list[[10]], rds_list[[11]], rds_list[[12]], rds_list[[13]] 
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
  new.reduction = "integrated")

TB_dataset[["RNA"]] <- JoinLayers(TB_dataset[["RNA"]])

TB_dataset <- FindNeighbors(TB_dataset, dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)

# Save & visualise
saveRDS(TB_dataset, "./data/integrated.rds")
plot_markers(TB_dataset, "Integrated BAL and PBMC samples from TB patients")

?DimPlot()
DimPlot(TB_dataset, reduction = "umap", group.by = "orig.ident")
IFeaturePlot(TB_dataset, reduction = "integrated", feature = "PTPRC")
  
