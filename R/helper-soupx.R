# Ambient RNA correction 
# Retrieve raw and filtered alignment output to estimate contamination profile 
library(SoupX)
library(Seurat)

# Function
run_soupX <- function(filtered_path, raw_path, project_name){
  
  tod <- suppressWarnings(Read10X(raw_path))
  toc <- suppressWarnings(Read10X(filtered_path))
  check <- suppressWarnings(CreateSeuratObject(
    counts = toc,
    project = project_name,
    min.cells = 1))  # since so few cells were sequenced, matrices are small & all genes can be retained 
  
  if (length(Cells(check)) >= 1000) {
    
    message("Ambient RNA correction...")
    
    sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
      
    # Estimate the soup
    sc <- estimateSoup(sc)
      
    # clustering (required) of filtered matrix
    seurat_soup <- suppressWarnings(CreateSeuratObject(toc, min.cells = 0))
    seurat_soup <- suppressWarnings(SCTransform(seurat_soup, verbose = FALSE) %>%
                                        RunPCA(verbose = FALSE) %>%
                                        RunUMAP(dims = 1:30, verbose = FALSE) %>%
                                        FindNeighbors(dims = 1:30, verbose = FALSE) %>%
                                        FindClusters(verbose = FALSE))
      
    # add clusters to the channel
    meta.data <- seurat_soup@meta.data
    umap.embedding <- seurat_soup@reductions$umap@cell.embeddings
    sc <- suppressWarnings(setClusters(sc, setNames(meta.data$seurat_clusters, rownames(meta.data))))
    sc <- suppressWarnings(setDR(sc, umap.embedding, c("UMAP_1", "UMAP_2")))
      
    # with defined clusters, run the main Soup X function to calculate ambient RNA profile
    sc <- autoEstCont(sc, verbose = FALSE, doPlot = FALSE, forceAccept = TRUE)
      
    # Output integer matrix of soup-corrected reads (unzipped output)
    adj.matrix <- suppressWarnings(adjustCounts(sc, roundToInt = T))
      
    # Output results in seurat object
    seurat <- suppressWarnings(CreateSeuratObject(counts = adj.matrix,
                                                    min.cells = 1, # genes expressed in at least 1 cell to be retained 
                                                  project = project_name))
  } else {
    
    # Create seurat object straight from filtered output
    seurat <- Read10X(filtered_path)
    seurat <- suppressWarnings(CreateSeuratObject(
      counts = seurat,
      project = project_name,
      min.cells = 1))
    
  }
  
  return(seurat)
}

# Test 
# test <- run_soupX(
#   filtered_path = "./zipped/BAL_1098/BAL_1098Solo.out/Gene/filtered/",
#   raw_path =  "./zipped/BAL_1098/BAL_1098Solo.out/Gene/raw/",
#   project_name = "BAL_1098"
# )
# dim(test)
# [1] 23621  1500

  
