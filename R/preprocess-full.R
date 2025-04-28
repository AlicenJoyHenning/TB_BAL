# Pre-process single cell data

# Load libraries ---- 
packages <- c(
  "DoubletFinder", "DropletQC", "DropletUtils", "Matrix", "Seurat",
  "SoupX", "cowplot", "dplyr", "ggplot2", "gridExtra", "knitr",
  "multtest", "patchwork", "tidyr"
)

lapply(packages, library, character.only = TRUE)

# Load helper scripts ----
source("./R/helper-soupx.R")
source("./R/helper-rbc.R")
source("./R/helper-DoubletFinder.R")
source("./R/helper-DamageDetective.R")
source("./R/helper-plot.R")

# Core function ----
preprocess <- function(
    filtered_path, 
    raw_path, 
    project_name, 
    output = "~/Projects/TB_BAL_data/data/"
){
  # Ambient RNA correction 
  counts <- run_soupX(filtered_path, raw_path, project_name)
  
  # Filter red blood cell contamination 
  counts <- filter_rbc(counts)
  
  # Filter dead 
  counts <- run_DD(counts)
  
  # Filter doublets 
  counts <- run_DF(counts)
  
  # View sample
  plot_markers(counts, project_name, output)
  
  # Ensure correct label 
  counts$orig.ident <- project_name
  
  # Save sample
  output_path <- paste0(output, project_name, ".rds")
  saveRDS(counts, output_path)
  
  return()
}

# Run the samples -----

samples <- data.frame(
  BAL_1098 = c("~/Projects/TB_BAL_data/zipped/BAL_1098/BAL_1098Solo.out/Gene/filtered", 
               "~/Projects/TB_BAL_data/zipped/BAL_1098/BAL_1098Solo.out/Gene/raw"),
  BAL_1376 =  c("~/Projects/TB_BAL_data/zipped/BAL_1376/BAL_1376Solo.out/Gene/filtered", 
                "~/Projects/TB_BAL_data/zipped/BAL_1376/BAL_1376Solo.out/Gene/raw"),
  BAL_1483 = c("~/Projects/TB_BAL_data/zipped/BAL_1483/BAL_1483Solo.out/Gene/filtered", 
               "~/Projects/TB_BAL_data/zipped/BAL_1483/BAL_1483Solo.out/Gene/raw"),
  BAL_1523 = c("~/Projects/TB_BAL_data/zipped/BAL_1523/BAL_1523Solo.out/Gene/filtered", 
               "~/Projects/TB_BAL_data/zipped/BAL_1523/BAL_1523Solo.out/Gene/raw"),
  BAL_1566 = c("~/Projects/TB_BAL_data/zipped/BAL_1566/BAL_1566Solo.out/Gene/filtered", 
               "~/Projects/TB_BAL_data/zipped/BAL_1566/BAL_1566Solo.out/Gene/raw"),
  BAL_1676 = c("~/Projects/TB_BAL_data/zipped/BAL_1676/BAL_1676Solo.out/Gene/filtered", 
               "~/Projects/TB_BAL_data/zipped/BAL_1676/BAL_1676Solo.out/Gene/raw"),
  PBMC_1098 = c("~/Projects/TB_BAL_data/zipped/PBMC_1098/PBMC_1098Solo.out/Gene/filtered", 
                "~/Projects/TB_BAL_data/zipped/PBMC_1098/PBMC_1098Solo.out/Gene/raw"),
  PBMC_1376 = c("~/Projects/TB_BAL_data/zipped/PBMC_1376/PBMC_1376Solo.out/Gene/filtered", 
                "~/Projects/TB_BAL_data/zipped/PBMC_1376/PBMC_1376Solo.out/Gene/raw"), 
  PBMC_1483 = c("~/Projects/TB_BAL_data/zipped/PBMC_1483/PBMC_1483Solo.out/Gene/filtered", 
                "~/Projects/TB_BAL_data/zipped/PBMC_1483/PBMC_1483Solo.out/Gene/raw"), 
  PBMC_1484 = c("~/Projects/TB_BAL_data/zipped/PBMC_1484/PBMC_1484Solo.out/Gene/filtered", 
                "~/Projects/TB_BAL_data/zipped/PBMC_1484/PBMC_1484Solo.out/Gene/raw"), 
  PBMC_1523 = c("~/Projects/TB_BAL_data/zipped/PBMC_1523/PBMC_1523Solo.out/Gene/filtered", 
                "~/Projects/TB_BAL_data/zipped/PBMC_1523/PBMC_1523Solo.out/Gene/raw"), 
  PBMC_1566 = c("~/Projects/TB_BAL_data/zipped/PBMC_1566/PBMC_1566Solo.out/Gene/filtered", 
                "~/Projects/TB_BAL_data/zipped/PBMC_1566/PBMC_1566Solo.out/Gene/raw"), 
  PBMC_1676 = c("~/Projects/TB_BAL_data/zipped/PBMC_1676/PBMC_1676Solo.out/Gene/filtered", 
                "~/Projects/TB_BAL_data/zipped/PBMC_1676/PBMC_1676Solo.out/Gene/raw")
)


for (sample in names(samples)){
  message("Running ", sample, "...")
  # Define inputs
  filtered_path <- samples[[sample]][[1]]
  raw_path <- samples[[sample]][[2]]
  project_name <- sample
  
  # Run the function
  preprocess(filtered_path, raw_path, project_name)
}



