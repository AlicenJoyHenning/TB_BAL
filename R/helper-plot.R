library(Seurat)
library(ggplot2)


plot_markers <- function(counts, project_name){
  
  counts <- NormalizeData(counts, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(verbose = FALSE) %>%
    RunUMAP(dims = 1:10, verbose = FALSE)
  
  markers_plot <- FeaturePlot(
    counts, 
    features = c("EPCAM", "PDGFRA", "PTPRC", "MS4A1", "PECAM1",
                 "CD3E","NKG7", "CD14", "MARCO"), 
    cols = c("lightgrey", "#6ab5ba")) +
    theme(
      axis.text.x = element_text(angle = 90),
      axis.title.x = element_blank(),
      plot.caption = element_blank()
    )
  
  for (plot in 1:length(markers_plot)) {
    markers_plot[[plot]] <- markers_plot[[plot]] + NoLegend() + NoAxes() + theme(panel.border = element_rect(colour = "black", fill=NA, linewidth =1))
  }
  
  colours <- RColorBrewer::brewer.pal(12, 'Paired')
  colours_full <- c(colours, colours)
  
  plot_clusters <- DimPlot(counts,
                           pt.size = 1,
                           label = TRUE,
                           label.size = 6,
                           cols = colours_full) +
    labs(title = paste(project_name, " clustering and marker visualisation", sep = "")) + 
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(1), face = "bold"),
      panel.border = element_rect(colour = "black", fill=NA, linewidth =1))
  
  UMAP_markers <- plot_clusters + wrap_plots(markers_plot) + plot_layout(ncol = 2)
  
  ggsave(plot = UMAP_markers, 
         filename = paste0("./plots/", project_name, ".png"), 
         width = 12,         
         height = 6,         
         dpi = 300 
        )
}
