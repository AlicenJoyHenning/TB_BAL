# Filter RBCs 

filter_rbc <- function(counts){
  # Estimate red blood cell contamination 
  hemo_gene <- c("HBA1", "HBA2", "HBB")
  
  counts$rbc <- PercentageFeatureSet(
    object = counts,
    features = intersect(hemo_gene, rownames(counts@assays$RNA)),
    assay = "RNA"
  )
  
  counts <- subset(counts, rbc <= 5)
  
  return(counts)
}