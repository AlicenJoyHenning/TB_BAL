# Retain true, single cells using DoubletFinder 
library(DoubletFinder)

run_DF <- function(counts) {
  # Prepare object (required)
  counts_DF <- NormalizeData(counts, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:10, verbose = FALSE) %>%
    FindClusters(verbose = FALSE)
  
  # Open a connection to a temporary file
  tmp_file <- tempfile()
  tmp_conn <- file(tmp_file, open = "wt")
  
  # Redirect output to the temporary file
  sink(tmp_conn)
  sink(tmp_conn, type = "message")
  
  # Ensure sinks are reset and file is closed/deleted even if an error occurs
  on.exit({
    sink(NULL)
    sink(NULL, type = "message")
    close(tmp_conn)
    if (file.exists(tmp_file)) file.remove(tmp_file)
  }, add = TRUE)
  
  # pK Identification (no ground-truth)
  sweep.res.list <- paramSweep(counts_DF, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  pK <- bcmvn %>%
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  
  # Homotypic Doublet Proportion Estimate
  annotations <- counts_DF@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075 * nrow(counts_DF@meta.data))  # adjust as needed
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  nExp_poi.adj <- nExp_poi.adj - (nExp_poi.adj * (1 / 7))
  
  # Run DoubletFinder
  counts_DF <- doubletFinder(
    counts_DF,
    PCs = 1:10,
    pN = 0.25,
    pK = pK,
    nExp = nExp_poi.adj,
    reuse.pANN = FALSE
  )
  
  # Extract DoubletFinder classifications
  df_col_name <- grep("^DF\\.classifications", 
                      colnames(counts_DF@meta.data), value = TRUE)
  doublet_classifications <- counts_DF@meta.data[[df_col_name]]
  counts$DF <- doublet_classifications
  
  # Subset to singlets
  counts <- subset(counts, DF == "Singlet")
  
  return(counts)
}
