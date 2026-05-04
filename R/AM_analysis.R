# Script for analysis sample compositions
library(dplyr)
library(ggplot2)
library(Seurat)

# Retrieve data ----
BAL <- readRDS("~/Projects/TB_BAL_data/data/integrated/TB_BAL_dataset.rds")


# Visualise to check right data
DimPlot(BAL, label = TRUE, repel = TRUE)

# Isolate AMs ----
FeaturePlot(BAL, features = c("MARCO")) 
BAL_AM <- subset(BAL, BAL$seurat_clusters %in% c(0, 1, 2))
DimPlot(BAL, group.by = "orig.ident")

# Basic differentiation ----
FeaturePlot(BAL,
    
    features = c(
    "FABP4", "MARCO"
#"CD169", "CD163", "IFNA", "IL10", "IL6", "TNF", "TGFb")
))


# BAL integration ------

# Samples included based on valid barcode % > 90 & % mapped > 74 
data_dir <- "~/Projects/TB_BAL_data/data/preprocessed/"
BAL_1098 <- readRDS(paste0(data_dir, "BAL_1098.rds"))
BAL_1376 <- readRDS(paste0(data_dir, "BAL_1376.rds"))
BAL_1483 <- readRDS(paste0(data_dir, "BAL_1483.rds"))
BAL_1566 <- readRDS(paste0(data_dir, "BAL_1566.rds"))
BAL_1676 <- readRDS(paste0(data_dir, "BAL_1676.rds"))

rds_list <- list(BAL_1098, BAL_1376, BAL_1483, BAL_1566, BAL_1676)

cells <- c(
  unique(rds_list[[1]]$orig.ident), unique(rds_list[[2]]$orig.ident), 
  unique(rds_list[[3]]$orig.ident), unique(rds_list[[4]]$orig.ident), 
  unique(rds_list[[5]]$orig.ident)
)
  
# Merge data
TB_BAL_dataset <- merge(
  x = rds_list[[1]], 
  y = c(
    rds_list[[2]], rds_list[[3]], rds_list[[4]],rds_list[[5]]
  ),
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

TB_BAL_integrated <- IntegrateLayers(
  TB_BAL_dataset, 
  method = RPCAIntegration, # CCAIntegration
  assay = "RNA",
  orig.reduction = "pca", 
  k.weight = 66, # Reduced from 100 due to small cell number
  new.reduction = "integrated")

# Counts separate in v5, must join
TB_BAL_integrated[["RNA"]] <- JoinLayers(TB_BAL_integrated[["RNA"]])
TB_BAL_integrated # check, should be only one data, counts, scale.data 

# Tried with dims 1:10, not seeing good separation between myeloid and lymphoid cells 
TB_BAL_integrated <- FindNeighbors(TB_BAL_integrated, reduction = "integrated", dims = 1:10) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:10)

# Visualise integration
DimPlot(TB_BAL_integrated) # general clusters? Not really 
DimPlot(TB_BAL_integrated, group.by = "orig.ident") # samples separate, not integrated


# Alveolar Macrophages (AMs) signals ----

# Mature
FeaturePlot(TB_BAL_integrated, features = c(
    "MARCO", "PPARG", # all except small cluster at top - # 10 seurat cluster
    "FABP4", # all, except # 10 and extending 
    "SIGLEC1"
))

# Early
FeaturePlot(TB_BAL_integrated, features = c(
    "S100A8", # all same, not informative
    "VCAN", "FCN1" # highlight the same few cells 
))

# Other cell types or artifacts 
library(ggplot2)
TB_BAL_integrated$mt.percent <- PercentageFeatureSet(TB_BAL_integrated, pattern = "^MT-")

ggplot(TB_BAL_integrated@meta.data, aes(y = mt.percent, x = nFeature_RNA)) + 
  geom_point() 


# Running DamageDetection
library(DamageDetective)

# Finding ideal ribosomal penalty
penalty <- select_penalty(TB_BAL_integrated@assays$RNA$counts)

# Running detection with penalty 
TB_BAL_damage_free <- detect_damage(
    count_matrix = TB_BAL_integrated@assays$RNA$counts, 
    ribosome_penalty = penalty
)

# transfer labels & filter
TB_BAL_integrated$DamageDetective <- TB_BAL_damage_free$output$pANN
TB_BAL_filtered <- subset(TB_BAL_integrated, DamageDetective == "cell")

# Run clustering again 
TB_BAL_filtered <- FindNeighbors(TB_BAL_filtered, reduction = "integrated", dims = 1:10) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:10)


DimPlot(TB_BAL_filtered, label = TRUE, repel = TRUE) # Better separation of clusters
# Cluster 8 = t cells 
# Minimal B cells ~ cluster 10 at edge of cluster 3

FeaturePlot(TB_BAL_filtered, features = c(
 "PTPRC", # all immune?
 "CD3E",  # T?
 "MS4A1" # B?
))

FeaturePlot(TB_BAL_filtered, features = c(
    "MARCO", "PPARG"
))

# Remaining clusters all seem to be AM but show a distinct batch effect: 
# - BAL_1376 , BAL_1483 , BAL_1566, BAL_1676  cluster strongly on their own, BAL 1098 spread out

# First isolating the Alveolar Macrophages (AMs)
TB_BAL_AM <- subset(TB_BAL_filtered, !TB_BAL_filtered$seurat_clusters %in% c(8,10)) # Removing T & B

# Recluster 

TB_BAL_AM <- NormalizeData(TB_BAL_AM) %>%
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()

TB_BAL_AM <- FindNeighbors(TB_BAL_AM, reduction = "integrated", dims = 1:10) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:10)

DimPlot(TB_BAL_AM)
DimPlot(TB_BAL_AM, group.by = "orig.ident")
# Still strong batch effect for 1676 and 1376

# Testing if the main genes dividing 1676 and 1376 are technical. 
# Renaming clusters for testing 
TB_BAL_AM$batch <- ifelse(TB_BAL_AM$orig.ident %in% c("BAL_1676"), "BAL_1676", 
    ifelse(TB_BAL_AM$orig.ident %in% c("BAL_1376"), "BAL_1376", "combined"))

TB_BAL_AM <- SetIdent(TB_BAL_AM, value = "batch")

# Find top diff. exp genes separating 1676 from the rest
markers_1676 <- FindMarkers(TB_BAL_AM, 
                            ident.1 = "BAL_1676", 
                            min.pct = 0.25, 
                            logfc.threshold = 0.25)
head(markers_1676[markers_1676$avg_log2FC > 0, ], n = 10)
# Definitely techincal noise 

# Find markers for BAL_1376 vs all others 
markers_1376 <- FindMarkers(TB_BAL_AM, 
                            ident.1 = "BAL_1376", 
                            min.pct = 0.25, 
                            logfc.threshold = 0.25)
markers_1376 <- markers_1376[order(markers_1376$avg_log2FC, decreasing = TRUE), ]
head(markers_1376[markers_1376$avg_log2FC > 0, ], n = 10)
tail(markers_1376[markers_1376$avg_log2FC < 0, ], n = 10)



# Finding that 1376 is true biological signal difference!
FeaturePlot(TB_BAL_AM, features = c(
    # "MARCO", "ITGAM" # high in early 
    "C8B", "IFI27", # high in mature
))


# differentiation ie CD206, CD169, CD163, IFNA, IL10, IL6, TNF, TGFb?
FeaturePlot(TB_BAL_AM, features = c(
    # "MRC1", # high all 
 # "SIGLEC1" # low all 
    #"CD163" # low but slightly differen 
   # "IFNAR1" # low but slightly different
   # "IL10" , # None 
   # "IL6", # None 
    #"TNF" # Very low everywhere 
    "TGFB1" # medium everywhere 
))

# Other signals 
# From 1376
FeaturePlot(TB_BAL_AM, features = c(
    "CYP1B1", "ITGAM", "KLF2", "EMP1" , #only 1376
    "FN1", # much higher in 1376
    
))
DimPlot(TB_BAL_AM, group.by = "orig.ident")
# Likely early recruited population of alveolar macrophages in 1376

# What is cluster 8 being separated by?
TB_BAL_AM <- SetIdent(TB_BAL_AM, value = "seurat_clusters")

# Find top diff. exp genes separating cluster 8 from the rest
markers_cluster_8 <- FindMarkers(TB_BAL_AM, 
                            ident.1 = 8, 
                            min.pct = 0.25, 
                            logfc.threshold = 0.25)
head(markers_cluster_8[markers_cluster_8$avg_log2FC > 0, ], n = 10)
# Cluster 8 seems to be clear proliferating population of AMs!

FeaturePlot(TB_BAL_AM, features = c(
    "MKI67", "TYMS", "TOP2A", "CDK1" # all high in cluster 8
))

# Confirming cluster 8 with cell cycle scoring (built in Seurat lists)
s_genes <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes

TB_BAL_AM <- CellCycleScoring(
  object = TB_BAL_AM, 
  s.features = s_genes, 
  g2m.features = g2m_genes, 
  set.ident = TRUE
)

DimPlot(TB_BAL_AM, reduction = "umap", group.by = "Phase", pt.size = 2)

FeaturePlot(TB_BAL_AM, pt.size = 2,
    features = c(
    #"ITGAM" # "MARCO"
    "EMP1", 
    "CYP1B1" # definitiely only early
))


# For plots for analysis 
# Define groups for the clusters & show dot plot why
TB_BAL_AM$celltype <- ifelse(TB_BAL_AM$seurat_clusters == 8, "Proliferating", "-")
TB_BAL_AM$celltype <- ifelse(TB_BAL_AM$batch == "BAL_1376", "Early", TB_BAL_AM$celltype )
TB_BAL_AM$celltype <- ifelse(TB_BAL_AM$celltype == "-", "Mature", TB_BAL_AM$celltype )

# See if division makes sense 
DimPlot(TB_BAL_AM, group.by = "celltype", pt.size = 2) # Okay but some early & mature cells in wrong cluster 

# By cell cycle phase
TB_BAL_AM$celltype <- ifelse((TB_BAL_AM$Phase == "S") & (TB_BAL_AM$seurat_clusters == 8), "Proliferating", TB_BAL_AM$celltype)

# Make column for expression of NB genes, which is the main gene separating early vs mature
TB_BAL_AM$ITGAM <- TB_BAL_AM@assays$RNA$data["ITGAM", ]
TB_BAL_AM$CYP1B1 <- TB_BAL_AM@assays$RNA$data["CYP1B1", ]
TB_BAL_AM$EMP1 <- TB_BAL_AM@assays$RNA$data["EMP1", ]


# see distribution of ITGAM expression
ggplot(TB_BAL_AM@meta.data, aes(x = ITGAM, y = nFeature_RNA)) + geom_point() 
# clear group of cells with no ITGAM expression (aka zero)

# Targetting populations 
TB_BAL_AM$celltype <- ifelse((TB_BAL_AM$ITGAM == 0) & (TB_BAL_AM$celltype == "Early"), "Mature", TB_BAL_AM$celltype)
DimPlot(TB_BAL_AM, group.by = "celltype", pt.size = 2)

TB_BAL_AM$celltype <- ifelse(TB_BAL_AM$batch == "BAL_1376", "Early", TB_BAL_AM$celltype )

TB_BAL_AM$celltype <- ifelse((TB_BAL_AM$CYP1B1 > 0.5), "Early", TB_BAL_AM$celltype)
TB_BAL_AM$celltype <- ifelse((TB_BAL_AM$EMP1 > 1) & (TB_BAL_AM$celltype != "Proliferating"), "Early", TB_BAL_AM$celltype)

DimPlot(TB_BAL_AM, group.by = "celltype", pt.size = 2) # Better separation of early vs mature now



# Plotting -----

# Dot plot of top genes separating early vs mature vs proliferating 
TB_BAL_AM$celltype <- factor(TB_BAL_AM$celltype, levels = c("Early", "Proliferating", "Mature"))

DotPlot(
    TB_BAL_AM, 
    features = c("ITGAM", "CYP1B1", "EMP1", "MKI67", "TYMS", "TOP2A"), 
    group.by = "celltype"
) + coord_flip()

