# devtools::install_github(
#   "AlicenJoyHenning/DamageDetective",
#   build_vignettes = TRUE
# )
library(DamageDetective)
library(Seurat)

run_DD <- function(counts){

  filtered_matrix <- detect_damage(
    count_matrix = counts@assays$RNA$counts,
    ribosome_penalty = 1,
    seed = 7
  )
  
  counts$DamageDetective <- ifelse(filtered_matrix$output$DamageDetective > 0.5, 
                                   "damaged", "cell")
  counts <- subset(counts, DamageDetective == "cell")
  
  
  # Save sample
  output_path <- paste0("./plots/", project_name, ".png")
  ggsave(plot = filtered_matrix$plot,  
         filename = output_path, 
         width = 6,         
         height = 4,         
         dpi = 300)
  
  return( counts)
}
