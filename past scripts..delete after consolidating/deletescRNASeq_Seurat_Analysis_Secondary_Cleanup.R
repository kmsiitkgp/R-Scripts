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
#                     REMOVE MIXED POPULATION (SECONDARY)                      #
#******************************************************************************#

# Identify and annotate mixed population clusters using UMAP of module scores
# NOTE: Use Module scores plotted as UMAP, Violin plots to decide if a cluster
# is contaminant. If values on violin plots are less than 0.1, ignore and dont
# consider as contaminant. 

if (proj=="scRNASeq_BBN_C57B6"){
  reduc <- "Harmony"
  Mixed <- list("Epithelial" = c(), 
                "Fibroblasts" = c(27),              
                "Myeloid" = c(),                     
                "Lymphoid" = c())
}

if (proj=="scRNASeq_BBN_Rag"){
  reduc <- "rpca"
  Mixed <- list("Epithelial" = c(), 
                "Fibroblasts" = c(),            
                "Myeloid" = c(),                    
                "Lymphoid" = c()) 
}

if (proj=="scRNASeq_Chen"){
  Mixed <- list("Epithelial" = c(), 
                "Fibroblasts" = c(),              
                "Myeloid" = c(),                     
                "Lymphoid" = c()) 
}

if (proj=="scRNASeq_Jinfen"){
  Mixed <- list("Epithelial" = c(1), 
                "Fibroblasts" = c(),              
                "Myeloid" = c(),                     
                "Lymphoid" = c()) 
}

if (proj=="scRNASeq_HRA003620"){
  Mixed <- list("Epithelial" = c(), 
                "Fibroblasts" = c(),              
                "Myeloid" = c(),                     
                "Lymphoid" = c(31)) 
}

# Remove mixed population clusters
for (celltype in c("Epithelial", "Fibroblasts", "Myeloid", "Lymphoid")){
  
  if (lengths(Mixed)[[celltype]] != 0){
    res <- "1.4"
    reduc <- "Harmony"
    remove_mixed_clusters(res, reduc, celltype, Mixed)
    plot_conserved_modules(res, reduc, celltype, "All Markers")
  }

  res <- "0.4"
  reduc <- "Harmony"
  get_markers(res, reduc, celltype)
}