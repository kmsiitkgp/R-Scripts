proj <- "scRNASeq_Krizia"
species <- "Mus musculus"
assay <- "RNA"

clusters <- list("Hepatocytes"         = c(),
                 "Pancreatic.Acinar"   = c(),
                 "Pancreatic.Islet"    = c(),
                 "B.Plasma"            = c(11),
                 "T.NK"                = c(7),
                 "Fibroblasts"         = c(3,4),
                 "Macrophages"         = c(6,8),
                 "Mast"                = c(2),
                 "Platelets"           = c(15),  
                 "Dendritic"           = c(9),
                 "Endothelial"         = c(5,13),
                 "Lymph.Endothelial"   = c(),
                 "Myocytes.Myofibroblasts"  = c(10),
                 "CAFs"                = c(),
                 "Epithelial"          = c(1,12),
                 "Neurons"             = c(),
                 "Epithelial.I"        = c(),
                 "Epithelial.II"       = c(),
                 "Epithelial.III"      = c(),
                 "Epithelial.IV"       = c(),
                 "Unclassified"        = c(14))

subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
reduction <- "umap.harmony"
color.col <- "Cell.Type"
filename <- "UMAP.Cell.Types.Group"
output_path <- proj.params$seurat_dir
split.col <- "Condition"
plot_umap(subset_seurat, reduction, color.col, filename, output_path, split.col)