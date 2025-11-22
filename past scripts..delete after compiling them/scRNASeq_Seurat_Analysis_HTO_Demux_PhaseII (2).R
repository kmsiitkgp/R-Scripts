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

# NOTE: Run Phase I first. Then, check diagnostic plots and run Phase II.

#******************************************************************************#
#                                   PHASE II                                   #
#******************************************************************************#

# Create a list of samples which will be added to each barcode.
# Since folder names correspond to sample name, we just use list.files()
samples <- list.files(path = hto_matrix_path)

# NOTE: STOP!!! Check UMAP diagnostic plots from previous step before proceeding
# to next step. Based on diagnostic plots, subset singlets from appropriate 
# seurat object based on appropriate method.

for (s in samples){
  
  #******************EXTRACT SINGLETS & SAVE AS SEURAT OBJECT******************#
  
  # IMPORTANT: Define which margin and method to use for each sample
  if (s %in% c("B1","B3","B4")){
    margin <- "1"
    method <- "MULTI_classification.global"
  } else{
    margin <- "2"
    method <- "MULTI_classification.global"
  }
  cat("\n", s, ":", margin, ":", method, "\n")
  
  # Read demuxed seurat object
  sample.seurat <- readRDS(paste0(demux_results, "margin", margin, "/", s, ".rds"))
  Idents(sample.seurat) <- method
  singlet <- base::subset(x = sample.seurat,
                          idents = c("Singlet"))
  
  # IMPORTANT: Define sample specific subsetting to remove BAD HTOs
  if (s %in% c("B1","B4")){
    singlet <- base::subset(x = singlet,
                            MULTI_ID == "HTO-E", 
                            invert = TRUE)
  }
  
  # Assign proper idents based on demux method used
  y <- dplyr::if_else(method == "MULTI_classification.global", "MULTI_ID", "hash.ID")
  
  # Add a new column "HTO_Final" that contains the finalized classification
  singlet@meta.data <- singlet@meta.data %>% dplyr::mutate(HTO_Final = get(y))
  
  cat("\nNumber of singlets:", nrow(singlet@meta.data), "\n")
  
  # Save the singlets with appropriate sample name
  saveRDS(singlet, file=paste0(demux_results, "singlets/", s, ".rds"))
  
  #***************************DIAGNOSTICS ON SINGLETS**************************#
  
  # Cluster and visualize cells using the usual scRNA-seq workflow, and examine
  # for the potential presence of batch effects. Since RNA assay was normalized, 
  # we don't need to normalize again. Even if use normalize, results will be same
  # as logNormalize is done within cells.
  
  # Select the top 1000 most variable features
  singlet <- Seurat::FindVariableFeatures(object = singlet,
                                          assay = "RNA",
                                          selection.method = "vst",
                                          nfeatures = 2000)
  
  # Scaling RNA data, we only scale the variable features here for efficiency
  singlet <- Seurat::ScaleData(object = singlet,
                               features = VariableFeatures(singlet),
                               assay = "RNA")
  
  # Run PCA
  singlet <- Seurat::RunPCA(object = singlet,
                            assay = "RNA",
                            features = VariableFeatures(singlet))
  
  # We select the top 10 PCs for clustering and UMAP
  singlet <- Seurat::FindNeighbors(object = singlet,
                                   reduction = "pca",
                                   dims = 1:40)
  
  singlet <- Seurat::FindClusters(object = singlet,
                                  resolution = 0.6,
                                  verbose = FALSE)
  
  singlet <- Seurat::RunUMAP(object = singlet,
                             reduction = "pca",
                             dims = 1:40)
  
  Idents(sample.seurat) <- y
  
  # Plot UMAP
  p1 <- Seurat::DimPlot(object = singlet,
                  group.by = y)
  
  ggsave(filename = paste0(s, "_", y, "_Singlet_UMAP.tiff"),
         plot = p1,
         device = "jpeg",
         path = diagnostics_path,
         scale = 1,
         width = 15,
         height = 7.5,
         units = c("in"),
         dpi = 600,
         limitsize = TRUE,
         bg = "white")
}

#******************************************************************************#