proj <- "scRNASeq_Vera"
species <- "Homo sapiens"
assay <- "RNA"

clusters <- list("Hepatocytes"         = c(),
                 "Pancreatic.Acinar"   = c(),
                 "Pancreatic.Islet"    = c(),
                 "B.Plasma"            = c(22),
                 "T.NK"                = c(3),
                 "Fibroblasts"         = c(5,12,30),
                 "Macrophages"         = c(4,13),
                 "Mast"                = c(14),
                 "Dendritic"           = c(),
                 "Endothelial"         = c(1,17,25),
                 "Lymph.Endothelial"   = c(),
                 "Myocytes.Myofibroblasts"  = c(6,23,27),
                 "CAFs"                = c(),
                 "Epithelial"          = c(),
                 "Neurons"             = c(),
                 "Epithelial.I"        = c(8,7,18,26),
                 "Epithelial.II"       = c(10,11,19,28,35),
                 "Epithelial.III"      = c(),
                 "Epithelial.IV"       = c(),
                 "Unclassified"        = c(2,9,15,16,20,21,24,29,31,32,33,34))

integ.ann <- annotate_manual_sc_sp(integ.final, clusters, resolution, reduction, suffix, seurat_results)

# Visualize annotations
subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
reduction <- "Harmony"
color.col <- "Cell.Type"
filename <- "UMAP.Cell.Types"
plot_umap(subset_seurat, reduction, color.col, filename, proj.params$seurat_dir)

# Visualize markers
idents <- "Sample"
features <- c("GH1") #,"IGF1")
split.col <- "Cell.Type"
filename <- "test"
plot_dot_plot(integ.ann, idents, features, filename, proj.params$seurat_dir, split.col)

integ.ann@meta.data <- integ.ann@meta.data %>%
  dplyr::mutate(Groups = paste(Disease, Condition, Tissue, sep=".")) %>% 
  dplyr::mutate(Groups = gsub(pattern=" ", "", Groups))

idents <- "Groups"
features <- c("GH1") #,"IGF1")
split.col <- "Cell.Type"
filename <- "test"
plot_dot_plot(integ.ann, idents, features, filename, proj.params$seurat_dir, split.col)

idents <- "Patient"
features <- c("GH1") #,"IGF1")
split.col <- "Cell.Type"
filename <- "test"
plot_dot_plot(integ.ann, idents, features, filename, proj.params$seurat_dir, split.col)
  
# Visualize features in feature plot
subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
reduction <- "Harmony"
features <- c("IGF1", rownames(integ.ann@assays$RNA)[grepl("^GH.*", rownames(integ.ann@assays$RNA))])  
features <- c(rownames(integ.ann@assays$RNA)[grepl("^TGFB.*", rownames(integ.ann@assays$RNA))])
split.col <- "Groups" #"Patient"
for (feature in features){
  filename <- feature
  plot_features(subset_seurat, feature, reduction, filename, proj.params$seurat_dir, split.col)
}



# clusters <- list("Hepatocytes"         = c(),
#                  "Pancreatic.Acinar"   = c(),
#                  "Pancreatic.Islet"    = c(),
#                  "B.Plasma"            = c(18),
#                  "T.NK"                = c(1),
#                  "Fibroblasts"         = c(4,14,30),
#                  "Macrophages"         = c(5,21),
#                  "Mast"                = c(15,32),
#                  "Dendritic"           = c(),
#                  "Endothelial"         = c(2,17,23),
#                  "Lymph.Endothelial"   = c(),
#                  "Myocytes.Myofibroblasts"  = c(6,19,28),
#                  "CAFs"                = c(),
#                  "Epithelial"          = c(),
#                  "Neurons"             = c(),
#                  "Epithelial.I"        = c(8,11,16,25,26),
#                  "Epithelial.II"       = c(7,10,27,34),
#                  "Epithelial.III"      = c(),
#                  "Epithelial.IV"       = c(),
#                  "Unclassified"        = c(3,9,12,13,20,22,24,29,31,33))
# integ.ann.133K <- annotate_manual_sc_sp(integ.final.133K, clusters, resolution, reduction, suffix, seurat_results)


for (g in unique(integ.final@meta.data$Groups)){
  obj <- subset(integ.final, Groups == g)
  
  
  
  p_list <- c(p_list, list(p))
}

