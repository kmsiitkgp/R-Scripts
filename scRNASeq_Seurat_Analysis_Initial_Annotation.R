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
#                          STEP 4: CLUSTER ANNOTATION                          #
#******************************************************************************#

# NOTE: If you assign same cell type to multiple clusters, verify that the
# clusters indeed are similar. This can be done by checking the clustering at a
# lower resolution. At lower resolution, if the clusters you assigned the same
# cell type are merged, then, your cell type identification is correct.
# In this example, cluster 15 and cluster 25 at resolution 1.4 seemed similar.
# Indeed, the 2 clusters are merged into cluster 5 at resolution 1.2 or
# cluster 0 at resolution 1

# (i) DO NOT USE symbols like "/" in cluster names as it will create errors when
# saving files based on cluster names.
# (ii) You can use Seurat::RenameIdents() and then add cell type identity to
# metadata using integrated_seurat$cell_type <- Idents(integrated_seurat1).

# Easier alternative is to use custom code below which has multiple benefits:
# (a) Using RenameIdents() permanently changes Idents. So, if we want to change
# annotation, we have to read the Seurat object again. To avoid this, we can
# store the annotated data in new object "integrated_seurat1" but this increases
# memory used by R
# (b) If we save the annotated "integrated_seurat1" as rds, it will take another
# 10GB space in cluster.
# By using custom code, we aren't altering the idents of the seurat object.
# This allows to rename the clusters freely on the same seurat object.
# Also, we can save simply overwrite "integrated_seurat" as rds saving space

# Renaming using RenameIdents() is NOT recommended
# integrated_seurat <- RenameIdents(object = integrated_seurat, "0" = "")

if (proj=="scRNASeq_BBN_C57B6"){
  # rpca works better especially for fibroblasts. In harmony, there are patches
  # of fibroblasts in middle of epithelial clusters
  reduc <- "rpca"
  clusters <- list("Fibroblasts" = c(4,12,22,8,39,11,18,27,16,24,40,21,42,50),
                   "Myofibroblasts" = c(35),
                   "Epithelial" = c(0,1,5,7,15,25,30,48,46,49),
                   "Myeloid - Macrophages" = c(3,20,29,32,36),
                   "Myeloid - DCs" = c(38),
                   "Myeloid - MDSC" = c(),
                   "Myeloid - Granulocytes" = c(),
                   "Myeloid - Mast" = c(9,14,17,19,26),
                   "Lymphoid - B" = c(2,47),
                   "Lymphoid - Plasma" = c(44),
                   "Lymphoid - T" = c(10,13,23,28),
                   "Lymphoid - NK" = c(33),
                   "Endothelial" = c(6,43,45,31),
                   "Endothelial - Lymphatic" = c(37),
                   "Neurons" = c(41),
                   "Myocytes" = c(),
                   "Erythrocytes" = c(34),
                   "Unclassified" = c())
}

if (proj=="scRNASeq_BBN_Rag"){
  # rpca works better especially for epithelial. In harmony, there are patches 
  # of epithelial cells in many other non-epithelial clusters.
  reduc <- "rpca"
  clusters <- list("Fibroblasts" = c(3,9,20,39,31,26),
                   "Myofibroblasts" = c(30),
                   "Epithelial" = c(18,28,0,1,2,4,5,6,7,12,13,14,15,16,21,22,25,29,35,37,41),
                   "Myeloid - Macrophages" = c(8,10,34,19,32),
                   "Myeloid - DCs" = c(33),
                   "Myeloid - MDSC" = c(),
                   "Myeloid - Granulocytes" = c(),
                   "Myeloid - Mast" = c(11,17,24,38),
                   "Lymphoid - B" = c(),
                   "Lymphoid - Plasma" = c(),
                   "Lymphoid - T" = c(),
                   "Lymphoid - NK" = c(),
                   "Endothelial" = c(23,27,36),
                   "Endothelial - Lymphatic" = c(40),
                   "Neurons" = c(),
                   "Myocytes" = c(),
                   "Erythrocytes" = c(),
                   "Unclassified" = c())
}

if (proj=="scRNASeq_Chen"){
  # rpca works better especially for epithelial. In harmony, there are patches 
  # of epithelial cells in many other non-epithelial clusters.
  reduc <- "rpca"
  clusters <- list("Fibroblasts" = c(14,19),
                   "Myofibroblasts" = c(20,27,28),
                   "Epithelial" = c(0,3,4,7,9,11,12,13,17,24,25,26,30,32,39,43),
                   "Myeloid - Macrophages" = c(22,31),
                   "Myeloid - DCs" = c(),
                   "Myeloid - MDSC" = c(),
                   "Myeloid - Granulocytes" = c(21),
                   "Myeloid - Mast" = c(33),
                   "Lymphoid - B" = c(18),
                   "Lymphoid - Plasma" = c(34,38,29),
                   "Lymphoid - T" = c(2,5,8,10,37,40,36),
                   "Lymphoid - NK" = c(15,16),
                   "Endothelial" = c(1,6,23,35,42),
                   "Endothelial - Lymphatic" = c(),
                   "Neurons" = c(41),
                   "Myocytes" = c(),
                   "Erythrocytes" = c(),
                   "Unclassified" = c())
}

if (proj=="scRNASeq_GSE222315"){
  # rpca works better especially for epithelial. In harmony, there are patches 
  # of epithelial cells in many other non-epithelial clusters.
  reduc <- "rpca"
  clusters <- list("Fibroblasts" = c(9,27),
                   "Myofibroblasts" = c(2,23,48,49),
                   "Epithelial" = c(4,6,7,12,15,17,18,19,22,26,28,30,32,35,39,43,45,46,24,38),
                   "Myeloid - Macrophages" = c(13,20),
                   "Myeloid - DCs" = c(37),
                   "Myeloid - MDSC" = c(),
                   "Myeloid - Granulocytes" = c(41),
                   "Myeloid - Mast" = c(10,52),
                   "Lymphoid - B" = c(0,34,54,36,44,47),
                   "Lymphoid - Plasma" = c(31,53),
                   "Lymphoid - T" = c(5,8,11,42,25,40),
                   "Lymphoid - NK" = c(1,14,16,29,33),
                   "Endothelial" = c(3,50,51,21),
                   "Endothelial - Lymphatic" = c(),
                   "Neurons" = c(),
                   "Myocytes" = c(),
                   "Erythrocytes" = c(),
                   "Unclassified" = c())
}

if (proj=="scRNASeq_HRA003620"){
  # rpca works better especially for epithelial. In harmony, there are patches 
  # of epithelial cells in many other non-epithelial clusters.
  reduc <- "rpca"
  clusters <- list("Fibroblasts" = c(10,18),
                   "Myofibroblasts" = c(37),
                   "Epithelial" = c(0,2,4,9,11,16,19,21,33,28),
                   "Myeloid - Macrophages" = c(13,15,23),
                   "Myeloid - DCs" = c(35),
                   "Myeloid - MDSC" = c(),
                   "Myeloid - Granulocytes" = c(),
                   "Myeloid - Mast" = c(20,22,25),
                   "Lymphoid - B" = c(7,32),
                   "Lymphoid - Plasma" = c(30),
                   "Lymphoid - T" = c(8,24,31,5,1,6,3,14,12,29,38,17,26),
                   "Lymphoid - NK" = c(),
                   "Endothelial" = c(27,36),
                   "Endothelial - Lymphatic" = c(34),
                   "Neurons" = c(),
                   "Myocytes" = c(),
                   "Erythrocytes" = c(),
                   "Unclassified" = c())
}

if (proj=="scRNASeq_Jinfen"){
  reduc <- "rpca"
  clusters <- list("Fibroblasts" = c(26),
                   "Myofibroblasts" = c(),
                   "Epithelial" = c(19),
                   "Myeloid - Macrophages" = c(0,8,29,24,22),
                   "Myeloid - DCs" = c(28),
                   "Myeloid - MDSC" = c(),
                   "Myeloid - Granulocytes" = c(),
                   "Myeloid - Mast" = c(2,6,9,14,16,20,21),
                   "Lymphoid - B" = c(27),
                   "Lymphoid - Plasma" = c(),
                   "Lymphoid - T" = c(5,7,17,23,25,1,3,4,10,11,12,13,30),
                   "Lymphoid - NK" = c(),
                   "Endothelial" = c(31),
                   "Endothelial - Lymphatic" = c(),
                   "Neurons" = c(),
                   "Myocytes" = c(),
                   "Erythrocytes" = c(),
                   "Unclassified" = c(15,18,32))
}

if (proj=="scRNASeq_Simon"){
  # rpca works better especially for epithelial. In harmony, there are patches 
  # of epithelial cells in many other non-epithelial clusters.
  reduc <- "rpca"
  clusters <- list("Fibroblasts" = c(13,21,27),
                   "Myofibroblasts" = c(12,19,11,18,28,29,5),
                   "Epithelial" = c(0,2,4,6,7,9,10,24,25,26,33,34),
                   "Myeloid - Macrophages" = c(),
                   "Myeloid - DCs" = c(),
                   "Myeloid - MDSC" = c(),
                   "Myeloid - Granulocytes" = c(31),
                   "Myeloid - Mast" = c(),
                   "Lymphoid - B" = c(),
                   "Lymphoid - Plasma" = c(20),
                   "Lymphoid - T" = c(14),
                   "Lymphoid - NK" = c(),
                   "Endothelial" = c(22),
                   "Endothelial - Lymphatic" = c(35),
                   "Neurons" = c(),
                   "Myocytes" = c(),
                   "Erythrocytes" = c(),
                   "Unclassified" = c(1,3,8,15,16,17,23,30,32))
}

if (proj=="scRNASeq_GSE217093"){
  # rpca works better especially for epithelial. In harmony, there are patches 
  # of epithelial cells in many other non-epithelial clusters.
  reduc <- "harmony"
  clusters <- list("Fibroblasts" = c(6,12,27,25,28,38,45,47,54),
                   "Myofibroblasts" = c(35),
                   "Epithelial" = c(0,1,2,3,4,5,7,11,13,14,15,19,20,22,26,30,31,33,39,40,41,44,46,48,51,52,53,55,56,59,61),
                   "Myeloid - Macrophages" = c(16,17,32,34),
                   "Myeloid - DCs" = c(50),
                   "Myeloid - MDSC" = c(),
                   "Myeloid - Granulocytes" = c(),
                   "Myeloid - Mast" = c(8,18,29,43,49),
                   "Lymphoid - B" = c(37),
                   "Lymphoid - Plasma" = c(36,42),
                   "Lymphoid - T" = c(9,57),
                   "Lymphoid - NK" = c(24),
                   "Endothelial" = c(10,21,23,58,60),
                   "Endothelial - Lymphatic" = c(),
                   "Neurons" = c(),
                   "Myocytes" = c(),
                   "Erythrocytes" = c(62),
                   "Unclassified" = c())
}

res <- 1.4
# reduc is defined for each proj within the if statement
celltype <- NULL
integ <- annotate_data_umap(res, reduc, celltype, clusters)
integ <- annotate_data_score(integ, celltype)
  
#******************************************************************************#
#       STEP 5: PERFORM INITIAL SUBTYPE ANALYSIS ON REQUIRED CELL TYPES        #
#******************************************************************************#

# NOTE: We have defined "Myeloid - MDSCs, Myeloid - Mast Cells etc..so, use
# "Myeloid" in for loop, NOT "Myeloid Cells"

# NOTE: The intensity of module scores using the same set of features on a 
# celltype specific Seurat object will be reduced as compared to the intial 
# Seurat object containing all celltypes. This is because most of the features 
# have high expression in a celltype specific seurat object which affects the 
# binning of genes in AddModuleScore(). Conversely, contaminating cells are 
# expected to have strong expression in celltype specific Seurat objects.

# An alternative to AddModuleScore() is AddModuleScore_UCell() which is not
# affected by subsetting seurat objects as the score is calculated differently.

# Load the full seurat object
integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn.rds"))

for (celltype in c("Epithelial", "Fibroblasts", "Myeloid", "Lymphoid")){
  
  # prep_data() filters celltypes based on cell_class column of metadata which 
  # is based on UCell scoring, not based on cluster based annotation we 
  # performed in the previous step
  filt <- prep_data(integrated_seurat, celltype)
  sct <- sctransform_data(filt)
  
  kweight <- sct@meta.data %>% dplyr::count(Sample) %>% dplyr::select(n) %>% min()/2
  kweight <- min(kweight, 100) 
  integ <- integrate_data(sct, kweight)
  integ <- cluster_data(integ, celltype)
  integ <- add_module_scores(integ, "All Markers")
  save_data(integ, celltype)
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