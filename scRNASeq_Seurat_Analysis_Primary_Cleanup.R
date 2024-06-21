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

# Almost always, some fibroblasts/myeloid/lymphoid etc cells get clustered 
# within an epithelial cluster and vice versa. Such contaminants usually tend
# to form distinct tiny clusters when we plot UMAPs of epithelial cells. 
# Traditionally, we would identify these clusters based on module scores in UMAP
# and remove them. However, we might have to run a second round of cleanup to 
# remove such contaminants. MOreover, one has to manually go through the UMAPs 
# and mark the contaminating clusters for each celltype. This is extremely time
# consuming.

# An alternative faster, better and accurate strategy is to identify 
# contaminants based on UCell scores. For instance, cells that have higher 
# epithelial module score than fibroblast module score within a fibroblast 
# cluster are most likely contaminants. Hence, we can mark such cells and 
# remove them. This doesnt need a second round of cleanup as UCell scores are 
# constant for each cell.

# TRADITIONAL METHOD USING UMAP (NOT RECOMMENDED)
# # Use Harmony as reduction since it is lenient as compared to rpca
# if (proj=="scRNASeq_BBN_C57B6"){
#   reduc <- "Harmony"
#   Mixed <- list("Epithelial" = c(19,18,25,28,6,31),
#                 "Fibroblasts" = c(27,29,30,10),
#                 "Myeloid" = c(23,31,22,25),
#                 "Lymphoid" = c(25,17,20))
# }
# 
# if (proj=="scRNASeq_BBN_Rag"){
#   reduc <- "Harmony"
#   Mixed <- list("Epithelial" = c(5,51,30,13,33,24,38,23,45),
#                 "Fibroblasts" = c(34,16,14,22,31),
#                 "Myeloid" = c(39,37,15,36,34,39),
#                 "Lymphoid" = c())
# }
# 
# if (proj=="scRNASeq_Chen"){
#   reduc <- "Harmony"
#   Mixed <- list("Epithelial" = c(37,21,36,18),
#                 "Fibroblasts" = c(16,13,14),
#                 "Myeloid" = c(15,11,17),
#                 "Lymphoid" = c(18,9,30))
# }
# 
# if (proj=="scRNASeq_GSE222315"){
#   reduc <- "Harmony"
#   Mixed <- list("Epithelial" = c(14,15,22,24),
#                 "Fibroblasts" = c(19,10),
#                 "Myeloid" = c(15,18,1,11),
#                 "Lymphoid" = c(28,34))
# }
# 
# if (proj=="scRNASeq_HRA003620"){
#   reduc <- "Harmony"
#   Mixed <- list("Epithelial" = c(23,26,31,5),
#                 "Fibroblasts" = c(15,14,11),
#                 "Myeloid" = c(10,18,27,22),
#                 "Lymphoid" = c(16,29,30))
# }
# 
# if (proj=="scRNASeq_Jinfen"){
#   reduc <- "Harmony"
#   Mixed <- list("Epithelial" = c(5),
#                 "Fibroblasts" = c(),
#                 "Myeloid" = c(12,19),
#                 "Lymphoid" = c(18))
# }
# 
# # Remove mixed population clusters
# for (celltype in c("Epithelial", "Fibroblasts", "Myeloid", "Lymphoid")){
# 
#   if (lengths(Mixed)[[celltype]] != 0){
#     res <- 1.4
#     remove_mixed_clusters(res, reduc, celltype, Mixed)
#   }
# }

# cell_class has annotations based on UCell scores while cell_type has 
# annotations based on UMAP. 
# We can identify all cells where cell_class != cell_type and mark them as 
# contaminants. However, this could be errorneous.
# Alternatively, identify all cells that have similar module scores for more 
# than one cell type. These cells could be doublets.

df <- integrated_seurat@meta.data %>% 
  dplyr::select(contains("UCell"))

# Find max of each column
m <- apply(X=df,MARGIN=2,FUN=max)

# Divide by max of each column
df <- sweep(x=df, MARGIN=2, STATS=m, FUN="/")


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