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
BAL_1523 <- readRDS(paste0(data_dir, "BAL_1523.rds"))
BAL_1566 <- readRDS(paste0(data_dir, "BAL_1566.rds"))
BAL_1676 <- readRDS(paste0(data_dir, "BAL_1676.rds"))
PBMC_1523 <- readRDS(paste0(data_dir, "PBMC_1523.rds")) # Helps to separate out cell types

rds_list <- list(BAL_1098, BAL_1376, BAL_1483, BAL_1523, BAL_1566, BAL_1676, PBMC_1523)
cells <- c(
  unique(rds_list[[1]]$orig.ident), unique(rds_list[[2]]$orig.ident), 
  unique(rds_list[[3]]$orig.ident), unique(rds_list[[4]]$orig.ident), 
  unique(rds_list[[5]]$orig.ident), unique(rds_list[[6]]$orig.ident), 
  unique(rds_list[[7]]$orig.ident)
)
  
# Merge data
TB_BAL_dataset <- merge(
  x = rds_list[[1]], 
  y = c(
    rds_list[[2]], rds_list[[3]], rds_list[[4]],rds_list[[5]], rds_list[[6]], rds_list[[7]]
  ),
  add.cell.ids = cells
)

# Integrate (Ensure v5)
class(TB_BAL_dataset[['RNA']]) 
packageVersion("Seurat")

# Normalise counts 
TB_BAL_dataset <- NormalizeData(TB_BAL_dataset) %>%
  FindVariableFeatures(selection.method = "mean.var.plot") 

# Identify junk (Pseudogenes, Ribosomal, Mitochondrial, and ENSG IDs)
junk_features <- grep("^RP[L|S]|^MT-|^ENSG|P[0-9]$|^HLA-DPA2$|^HLA-DPB2$", 
                      VariableFeatures(TB_BAL_dataset), value = TRUE)

# Update VariableFeatures to exclude the junk
VariableFeatures(TB_BAL_dataset) <- setdiff(VariableFeatures(TB_BAL_dataset), junk_features)

TB_BAL_dataset <- ScaleData(TB_BAL_dataset) %>% 
  RunPCA()

TB_BAL_integrated <- IntegrateLayers(
  TB_BAL_dataset, 
  method = RPCAIntegration, # CCAIntegration
  assay = "RNA",
  orig.reduction = "pca", 
  # k.weight = 66, # Reduced from 100 due to small cell number
  new.reduction = "integrated")

# Counts separate in v5, must join
TB_BAL_integrated[["RNA"]] <- JoinLayers(TB_BAL_integrated[["RNA"]])
TB_BAL_integrated # check, should be only one data, counts, scale.data 

# Tried with dims 1:10, not seeing good separation between myeloid and lymphoid cells 
TB_BAL_integrated <- FindNeighbors(TB_BAL_integrated, reduction = "integrated", dims = 1:10) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:10)

# Visualise integration
DimPlot(TB_BAL_integrated) # general clusters 
DimPlot(TB_BAL_integrated, group.by = "orig.ident") # samples separate, not integrated


# Alveolar Macrophages (AMs) signals ----

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
DimPlot(TB_BAL_filtered, group.by = "orig.ident") # Still strong batch effect, but better than before

FeaturePlot(TB_BAL_filtered, features = c(
 "PTPRC", # all immune?
 "CD3E",  # T?
 "MS4A1" # B?
))

# Cluster 8 = t cells 
# Minimal B cells ~ cluster 10 at edge of cluster 3

# Mature
FeaturePlot(TB_BAL_filtered, features = c(
    "MARCO", "PPARG", # all except small cluster at top - # 10 seurat cluster
    "FABP4", # all, except # 10 and extending 
    "SIGLEC1"
))

# Early
FeaturePlot(TB_BAL_filtered, features = c(
    "S100A8", # all same, not informative
    "VCAN", "FCN1" # highlight the same few cells 
))

# Remaining clusters all seem to be AM but show a distinct batch effect: 
# - BAL_1376 , BAL_1483 , BAL_1566, BAL_1676  cluster strongly on their own, BAL 1098 spread out

# First isolating the Alveolar Macrophages (AMs)
# Remove PBMCs 
TB_BAL_filtered <- subset(TB_BAL_filtered, !TB_BAL_filtered$orig.ident == "PBMC_1523") 
DimPlot(TB_BAL_filtered)

DimPlot(TB_BAL_filtered, group.by = "seurat_clusters", label = TRUE, repel = TRUE) 
FeaturePlot(TB_BAL_filtered, features = c(
    "CD3E"
   # "MS4A1"
))

# converting to loupe for finer labelling of t and b cells (not big enough for precise clusters)
create_loupe_from_seurat(TB_BAL_filtered, 
output_dir = "/home/alicen/Projects/TB_BAL_data/data/LoupeR",
output_name = "TB_BAL_lympocytes_unlabelled")

lymph_labels <- read.csv("/home/alicen/Projects/TB_BAL_data/data/LoupeR/BAL_lymphocytes_labelled.csv")

# Match barcodes to tranfer labels 
TB_BAL_filtered$groups <- lymph_labels$Celltypes[match(colnames(TB_BAL_filtered), lymph_labels$Barcode)]
TB_BAL_filtered$groups <- ifelse(TB_BAL_filtered$groups == "Unknown", "Myeloid", TB_BAL_filtered$groups)
DimPlot(TB_BAL_filtered, group.by = "groups", pt.size = 2)

saveRDS(TB_BAL_filtered, "~/Projects/TB_BAL_data/data/integrated/TB_BAL_lymphocytes_labelled.rds")
# Save for plot 

# Continue with isolation
TB_BAL_AM$MS4A1 <- TB_BAL_AM@assays$RNA$data["MS4A1", ]
TB_BAL_AM <- subset(TB_BAL_AM, !TB_BAL_AM$seurat_clusters %in% c(8)) # Removing  B
TB_BAL_AM <- subset(TB_BAL_AM, TB_BAL_AM$MS4A1 == 0) # Removing  B

DimPlot(TB_BAL_AM, group.by = "seurat_clusters", label = TRUE, repel = TRUE) # T and B removed
DimPlot(TB_BAL_AM, group.by = "orig.ident") # Nice and combined

# reclustering to see if i can get better clusters
TB_BAL_AM <- FindNeighbors(TB_BAL_AM, reduction = "integrated", dims = 1:10) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:10)
DimPlot(TB_BAL_AM, group.by = "orig.ident")

# Testing if the main genes dividing 1376 are technical. 
# Renaming clusters for testing 
TB_BAL_AM$batch <- ifelse(TB_BAL_AM$orig.ident %in% c("BAL_1376"), "BAL_1376", "combined")

TB_BAL_AM <- SetIdent(TB_BAL_AM, value = "batch")

# Find markers for BAL_1376 vs all others 
markers_1376 <- FindMarkers(TB_BAL_AM, 
                            ident.1 = "BAL_1376", 
                            min.pct = 0.25, 
                            logfc.threshold = 0.25)
markers_1376 <- markers_1376[order(markers_1376$avg_log2FC, decreasing = TRUE), ]
head(markers_1376[markers_1376$avg_log2FC > 0, ], n = 10)
tail(markers_1376[markers_1376$avg_log2FC < 0, ], n = 10)


# Find markers 
TB_BAL_AM <- SetIdent(TB_BAL_AM, value = "celltype")


markers_mature <- FindMarkers(TB_BAL_AM, 
                            ident.1 = "mature", 
                            min.pct = 0.25, 
                            logfc.threshold = 0.25)
markers_mature <- markers_mature[order(markers_mature$avg_log2FC, decreasing = TRUE), ]
head(markers_mature[markers_mature$avg_log2FC > 0, ], n = 10)

markers_proliferating <- FindMarkers(TB_BAL_AM, 
                            ident.1 = "proliferating", 
                            min.pct = 0.25, 
                            logfc.threshold = 0.25)
markers_proliferating <- markers_proliferating[order(markers_proliferating$avg_log2FC, decreasing = TRUE), ]
head(markers_proliferating[markers_proliferating$avg_log2FC > 0, ], n = 10)

# Finding: 1376 is true biological signal difference!
FeaturePlot(TB_BAL_AM, features = c(
    # "CD36", "ITGAM" # high in early 
    # "C8B", "IFI27", # high in mature
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

# For plots for analysis -----
DimPlot(TB_BAL_AM, reduction = "umap", group.by = "seurat_clusters", pt.size = 2, label = TRUE, repel = TRUE)

# Define groups for the clusters & show dot plot why
TB_BAL_AM$celltype <- ifelse(TB_BAL_AM$seurat_clusters == 8, "Proliferating", "-")
DimPlot(TB_BAL_AM, group.by = "celltype", pt.size = 2) # clear proliferating cluster but early vs mature not clear

# Using Loupe for wrappin tool 
remotes::install_github("10XGenomics/loupeR")
loupeR::setup()
library("loupeR")
create_loupe_from_seurat(TB_BAL_AM, 
output_dir = "/home/alicen/Projects/TB_BAL_data/data/LoupeR",
output_name = "TB_BAL_AM_unlabelled")

# Read in the manual labels 
labels <- read.csv("/home/alicen/Projects/TB_BAL_data/data/LoupeR/BAL_labels.csv")

# Match barcodes to tranfer labels 
TB_BAL_AM$celltype <- labels$Proliferating[match(colnames(TB_BAL_AM), labels$Barcode)]
TB_BAL_AM$celltype <- ifelse(is.na(TB_BAL_AM$celltype), "mature", TB_BAL_AM$celltype)

DimPlot(TB_BAL_AM, group.by = "celltype", pt.size = 2)

saveRDS(TB_BAL_AM, "~/Projects/TB_BAL_data/data/integrated/TB_BAL_AM.rds")

# Plotting -----

# Dot plot of top genes separating early vs mature vs proliferating 
TB_BAL_AM$celltype <- factor(TB_BAL_AM$celltype, levels = c("mature", "proliferating", "early"))

markers_dot <- DotPlot(cols = c( "lightgrey", "#1D78B4"),
    TB_BAL_AM, 
    features = c( "MARCO", "THBS1","FABP4", "ITGAM", "CYP1B1", "EMP1", "MKI67", "TYMS", "TOP2A"), 
    group.by = "celltype"
) + 
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_blank())


ggsave("./plots/composition/markers_dot.svg", markers_dot, width = 9, height = 4)


# Basc UMAP 
umap_plot <- DimPlot(TB_BAL_AM, reduction = "umap", group.by = "celltype", pt.size = 2) + 
  theme(
    legend.position = "bottom",
    plot.title = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.ticks = element_blank(),
    axis.text = element_blank()) 

ggsave("./plots/AM_clusters_umap.png", umap_plot , width = 4, height = 4)

# Feature Plots 
# Generate the initial multi-panel plot
feature <- FeaturePlot(TB_BAL_AM, 
                       pt.size = 2,
                       features = c("MARCO", "CYP1B1", "MKI67"), 
                       order = TRUE, 
                       ncol = 3,
                       cols = c("lightgrey", "#1D78B4"))

# Iterate through each panel to apply the uniform style
for (i in 1:length(feature)) {
  feature[[i]] <- feature[[i]] + 
    NoAxes() + 
    theme(
      plot.title = element_text(face = "italic", hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = "bottom" # Bottom usually looks cleaner for a horizontal row
    )
}

ggsave("./plots/AM_clusters_features.svg", feature, width = 9, height = 4)


# Proportions library(dplyr)


celltype_probs <- TB_BAL_AM@meta.data %>%
  # Ensure we are counting celltype within each sample
  group_by(orig.ident) %>%
  count(celltype) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

celltype_probs <- celltype_probs %>%
  mutate(orig.ident = reorder(orig.ident, proportion, 
                              FUN = function(x) sum(x[celltype_probs$celltype == "Mature AM"])))

plot <- ggplot(celltype_probs, 
               aes(x = proportion, y = orig.ident, fill = celltype)) + # Changed 'groups' to 'celltype'
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = scales::percent_format(), expand = c(0,0)) + 
  scale_fill_manual(values = c("early" = "#619CFF", 
                               "proliferating" = "#00BA38", 
                               "mature" = "#E26B63")) +
  labs(fill = "", x = "", y = NULL) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank() # Clean vertical look
  )

# Save the plot as an SVG
ggsave("./plots/BAL_subtype_proportion.svg", plot, width = 8, height = 4)


# Differentiation expression plots -----

low_feature <- FeaturePlot(TB_BAL_AM, 
                       pt.size = 2,
                       features = c("MRC1", "TGFB1", "SIGLEC1", "CD163", "IFNAR1", "IFNAR2", "IL6", "IL10"), 
                       order = TRUE, 
                       cols = c("lightgrey", "#1D78B4"))

# Iterate through each panel to apply the uniform style
for (i in 1:length(low_feature)) {
  low_feature[[i]] <- low_feature[[i]] + 
    NoAxes() + 
    theme(
      plot.title = element_text(face = "italic", hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = "bottom" # Bottom usually looks cleaner for a horizontal row
    )
}


ggsave("./plots/differentiation_features.svg", low_feature, width = 8, height = 10)
