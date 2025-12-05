#!/usr/bin/env Rscript

# Read and store variables from command line interface (CLI)
cli <- base::commandArgs(trailingOnly = TRUE) 
args <- base::strsplit(x = cli, split = "=", fixed = TRUE)

for (e in args){
  argname <- e[1]
  argval <- e[2]
  assign(argname, argval)
}

# Run the Custom_Functions.R script
path1 <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/Custom_Functions.R"
path2 <- "/hpc/home/kailasamms/projects/scRNASeq/Custom_Functions.R"
if (file.exists(path1)) {
  source(path1)
} else if (file.exists(path2)) {
  source(path2)
}

# ---- PROJECT SET UP ----

# parent_dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data"
# gmt_dir    <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets"
# scripts_dir <- NULL

parent_dir  <- "/hpc/home/kailasamms/scratch"
gmt_dir     <- "/hpc/home/kailasamms/projects/GSEA_genesets"
scripts_dir <- "/hpc/home/kailasamms/projects/scRNASeq"

contrasts <- c()
deseq2.override <- list()
heatmap.override <- list()
volcano.override <- list()

proj.params <- setup_project(proj             = proj,
                             species          = species,  #"Mus musculus", "Homo sapiens"
                             contrasts        = contrasts,
                             parent_dir       = parent_dir,
                             gmt_dir          = gmt_dir,
                             scripts_dir      = scripts_dir,
                             deseq2.override  = deseq2.override,
                             heatmap.override = heatmap.override,
                             volcano.override = volcano.override)

# ---- Single cell Workflow ----

assay <- "RNA"
seurat_list <- list() 
raw_metadata <- data.frame()
reference_samples <- NULL  # change if you want reference based integration
samples <- list.files(path = proj.params$filt_matrix_dir)

# Loop through each sample specified in the 'samples' vector
for (sample in samples){

  cat("\n--- Processing Sample:", sample, " ---\n")

  # Load cellranger results (raw UMI matrix/barcodes)
  s.obj <- load_cellranger(sample = sample,
                           matrix_dir = proj.params$raw_matrix_dir)

  # Identify empty droplets using DropletUtils & Cell Ranger's filtered output
  s.obj <- classify_dropletutils(sample_seurat = s.obj)
  s.obj <- classify_cellranger(sample_seurat = s.obj,
                               filt_matrix_dir = proj.params$filt_matrix_dir)

  # Calculate standard QC metrics and identify 'High Quality' barcodes
  s.obj <- calc_qc_metrics(sample_seurat = s.obj,
                           assay = assay)

  # Identify singlets/doublets from 'High Quality' barcodes using DoubletFinder and scDblFinder
  s.obj <- classify_doubletfinder(sample_seurat = s.obj)
  s.obj <- classify_scdblfinder(sample_seurat = s.obj)

  # Retain 'High Quality Singlets'
  result <- filter_singlets(sample_seurat = s.obj)

  # Capture metadata from the raw object for QC visualization
  raw_metadata <- dplyr::bind_rows(raw_metadata, result$raw_metadata)

  # Store the filtered Seurat object in a list, named by sample
  seurat_list[[sample]] <- result$sample_seurat
}

# Generate QC plots
plot_qc(raw_metadata = raw_metadata,
        output_dir = proj.params$seurat_dir)

# Merge all filtered, high-quality Seurat objects from 'seurat_list' into one object
# This also typically loads and integrates external metadata using 'metafile'.
filt.obj <- merge_filtered(seurat_list = seurat_list,
                           assay = assay,
                           meta_file = proj.params$metafile,
                           output_dir = proj.params$seurat_dir)

# Run SCTransform on the merged object
sct.obj <- run_sctransform(filtered_seurat = filt.obj,
                           assay = assay,
                           s_genes = proj.params$cell_cycle$S,
                           g2m_genes = proj.params$cell_cycle$G2M)

# Integrate to remove technical batch effects between samples
integ.obj <- integrate_sct_data(sct_seurat = sct.obj,
                                assay = assay)

# Perform clustering on the integrated data
integ.clust <- cluster_integrated_data(integrated_seurat = integ.obj,
                                       assay = assay)

# Filter out sparse clusters often associated with low-quality cells or debris
integ.final <- remove_sparse_clusters(integrated_seurat = integ.clust,
                                      assay = assay,
                                      output_dir = proj.params$seurat_dir)

# Calculate the optimal clustering res
optimal_res <- calc_optimal_resolution(integrated_seurat = integ.final,
                                       reduction = "Harmony",
                                       output_dir = proj.params$seurat_dir)

# Find differential expressed genes (markers)
resolutions <- unique(c(0.2, 0.4, 0.6, 0.8, 1, optimal_res))
for (res in resolutions){
  identify_markers(integrated_seurat = integ.final,
                   res = res,
                   reduction = "Harmony",
                   output_dir = proj.params$seurat_dir)
}

# Plot metrics
plot_metrics_post_integration(integrated_seurat = integ.final,
                              output_dir = proj.params$diagnostics_dir)

# Calculate module scores
integ.final <- calc_module_scores(integrated_seurat = integ.final,
                                  reduction = "Harmony",
                                  marker_file = proj.params$markerfile,
                                  output_dir = proj.params$seurat_dir)

# Annotate the clusters based on predefined markers
integ.ann <- annotate_clusters(integrated_seurat = integ.final,
                               reduction = "Harmony",
                               assay = assay,
                               output_dir = proj.params$seurat_dir)

##integ.ann <- readRDS(file.path(proj.params$seurat_dir, "integrated_seurat_ann.rds"))

# Run FindAllMarkers on annotated clusters
identify_markers(integrated_seurat = integ.ann,
                 res = NULL,
                 reduction = NULL,
                 cluster_col = "CellType",
                 output_dir = proj.params$seurat_dir)


# clusters <- list("Hepatocytes"              = c(),
#                  "Pancreatic.Acinar"        = c(),
#                  "Pancreatic.Islet"         = c(),
#                  "B.Plasma"                 = c(),
#                  "T.NK"                     = c(),
#                  "Fibroblasts"              = c(),
#                  "Macrophages"              = c(),
#                  "Mast"                     = c(),
#                  "Platelets"                = c(),
#                  "Dendritic"                = c(),
#                  "Endothelial"              = c(),
#                  "Lymph.Endothelial"        = c(),
#                  "Myocytes.Myofibroblasts"  = c(),
#                  "Neurons"                  = c(),
#                  "Epithelial"               = c(),
#                  "Met.Panc"                 = c(),
#                  "Primary.Panc"             = c(),
#                  "Ductal.Cells"             = c(),
#                  "Unclassified"             = c())
# integ.ann <- annotate_manual_sc_sp(integ.final, assay, clusters, res, reduction, proj.params$seurat_dir)
# 
# # Visualize annotations in dim plot
# subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
# reduction <- "umap.harmony"
# color.col <- "Cell.Type"
# filename <- "UMAP_Annotated"
# output_path <- proj.params$seurat_dir
# split.col <- NULL
# plot_umap(subset_seurat, reduction, color.col, filename, output_path, split.col)
# 
# # Visualize markers in dot plot
# subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
# assay <- "Spatial.016um"
# ident.1 <- "Cell.Type"
# ident.2 <- NULL
# features <- find_top_features(subset_seurat , assay, proj.params)
# #features <- c("GH1")
# filename <- "Dot_Plot_Markers"
# output_path <- proj.params$seurat_dir
# split.col <- NULL
# plot_dot_custom(subset_seurat, assay, proj.params, ident.1, ident.2, features, filename, output_path, split.col)
# 
# # Visualize clusters from each column metadata
# coi <- c("Sample", "Patient", "Condition", "Treatment", "Disease")
# subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
# reduction <- "umap.harmony"
# color.col <- "Cell.Type"
# output_path <- proj.params$seurat_dir
# for (split.col in coi){
#   if(!is.null(subset_seurat@meta.data[[split.col]] %>% unique())){
#     filename <- paste0("UMAP_", split.col)
#     plot_umap(subset_seurat, reduction, color.col, filename, output_path, split.col)
#   }
# }
# 
# # Tabulate population frequencies
# split.cols <- c("Sample", "Patient", "Condition", "Treatment", "Disease")
# subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
# output_path <- proj.params$seurat_dir
# tabulate_frequency(subset_seurat, split.cols, output_path)
# 
# # [OPTIONAL] Visualize features in feature plot
# subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
# reduction <- "umap.harmony"
# features <- c("IGF1", rownames(integ.ann@assays$RNA)[grepl("^GH.*", rownames(integ.ann@assays$RNA))])  
# filename <- "Feature_Plot"
# output_path <- proj.params$seurat_dir
# split.col <- NULL #"Groups" #"Patient"
# plot_features(subset_seurat, features, reduction, filename, output_path, split.col)
# 
# # ---- Pseudobulk Workflow ----
# 
# assay <- "RNA"
# integ.ann <- readRDS(file.path(proj.params$seurat_dir, paste0("integrated.seurat.", assay, ".ann.rds")))
# subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
# comparison.col <- "Comparisons"
# subset_seurat@meta.data <- subset_seurat@meta.data %>% dplyr::mutate(Comparisons = paste0(Treatment, ".", Cell.Type))
# p <- prep_pseudobulk(subset_seurat, comparison.col, assay)
# 
# contrasts <- c("Tumor.Fibroblasts-Normal.Fibroblasts",
#                "Tumor.T.NK-Normal.T.NK",
#                "Tumor.Epithelial-Normal.Epithelial",
#                "Tumor.Endothelial-Normal.Endothelial",
#                "Tumor.Myocytes.Myofibroblasts-Normal.Myocytes.Myofibroblasts",
#                "Tumor.Mast-Normal.Mast",
#                "Tumor.Macrophages-Normal.Macrophages",
#                "Tumor.Dendritic-Normal.Dendritic",
#                "Tumor.B.Plasma-Normal.B.Plasma")  
# proj.params$contrast <- contrasts
# 
# main_analysis(proj_dir, p$meta_data, p$read_data, proj.params)
# 
# # ---- Seurat DEG Workflow ----
# 
# subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
# comparison.col <- "Treatment"
# subset.col <- "Cell.Type"
# assay <- "RNA"
# output_path <- proj.params$deseq2_dir
# degs_df <- find_degs_seurat(subset_seurat, comparison.col, output_path, celltype.col, assay)
