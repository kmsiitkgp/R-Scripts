#!/usr/bin/env Rscript

# Read and store variables from CLI
cli <- commandArgs(trailingOnly = TRUE) 
args <- strsplit(cli, "=", fixed = TRUE)

for (e in args){
  argname <- e[1]
  argval <- e[2]
  assign(argname, argval)
}

# NOTE: All variables and functions are defined within the file below
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")

#******************************************************************************#
#                       REMOVE MIXED POPULATION (PRIMARY)                      #
#******************************************************************************#

# Identify and annotate mixed population clusters using UMAP of module scores
# NOTE: Use Module scores plotted as UMAP, Violin plots to decide if a cluster
# is contaminant. If values on violin plots are less than 0.1, ignore and dont
# consider as contaminant.

# Use Harmony as reduction since it is lenient as compared to rpca
if (proj=="scRNASeq_BBN_C57B6"){
  reduc <- "Harmony"
  Mixed <- list("Epithelial" = c(19,18,25,28,6,31), 
                "Fibroblasts" = c(27,29,30,10),              
                "Myeloid" = c(23,31,22,25),                     
                "Lymphoid" = c(25,17,20))
}

if (proj=="scRNASeq_BBN_Rag"){
  reduc <- "Harmony"
  Mixed <- list("Epithelial" = c(5,51,30,13,33,24,38,23,45), 
                "Fibroblasts" = c(34,16,14,22,31),              
                "Myeloid" = c(39,37,15,36,34,39),                     
                "Lymphoid" = c())
}

if (proj=="scRNASeq_Chen"){
  reduc <- "Harmony"
  Mixed <- list("Epithelial" = c(37,21,36,18), 
                "Fibroblasts" = c(16,13,14),              
                "Myeloid" = c(15,11,17),                     
                "Lymphoid" = c(18,9,30))
}

if (proj=="scRNASeq_GSE222315"){
  reduc <- "Harmony"
  Mixed <- list("Epithelial" = c(14,15,22,24), 
                "Fibroblasts" = c(19,10),              
                "Myeloid" = c(15,18,1,11),                     
                "Lymphoid" = c(28,34)) 
}

if (proj=="scRNASeq_HRA003620"){
  reduc <- "Harmony"
  Mixed <- list("Epithelial" = c(23,26,31,5), 
                "Fibroblasts" = c(15,14,11),              
                "Myeloid" = c(10,18,27,22),                     
                "Lymphoid" = c(16,29,30)) 
}

if (proj=="scRNASeq_Jinfen"){
  reduc <- "Harmony"
  Mixed <- list("Epithelial" = c(5), 
                "Fibroblasts" = c(),              
                "Myeloid" = c(12,19),                     
                "Lymphoid" = c(18)) 
}

# Remove mixed population clusters
for (celltype in c("Epithelial", "Fibroblasts", "Myeloid", "Lymphoid")){
  
  if (lengths(Mixed)[[celltype]] != 0){
    res <- 1.4
    remove_mixed_clusters(res, reduc, celltype, Mixed)
  }
}

#******************************************************************************#
#                       PERFORM A FINAL SUBTYPE ANALYSIS                       #
#******************************************************************************#

for (celltype in c("Epithelial", "Fibroblasts", "Myeloid", "Lymphoid")){
  
  # Load the integrated seurat object
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                      dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
  
  filt <- prep_data(integrated_seurat, celltype)
  sct <- sctransform_data(filt)
  
  kweight <- sct@meta.data %>% dplyr::count(Sample) %>% dplyr::select(n) %>% min()/2
  kweight <- min(kweight, 100) 
  integ <- integrate_data(sct, kweight)
  integ <- cluster_data(integ, celltype)
  res <- 1.4
  integ <- add_module_scores(res, celltype, "All Markers")
  #plot_pre_integration(sct)
  
  for (reduc in c("CCA", "RPCA", "Harmony", "JointPCA")){
    res <- 1.4
    plot_conserved_modules(res, reduc, celltype, "All Markers")
  }
  
  res <- 1.4
  reduc <- "Harmony"
  idents <- paste0("cluster.", res, ".", base::tolower(reduc))
  plot_post_integration(res, reduc, idents, celltype)
}

#******************************************************************************#