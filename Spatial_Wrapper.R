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

parent.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data"
gmt.dir    <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets"

parent.dir  <- "/hpc/home/kailasamms/scratch"
gmt.dir     <- "/hpc/home/kailasamms/projects/GSEA_genesets"
scripts.dir <- "/hpc/home/kailasamms/projects/scRNASeq"

proj.params <- setup_project(
  proj = proj,
  species = species,
  contrasts = contrasts,
  parent.dir = parent.dir,
  gmt.dir = gmt.dir,
  deseq2.override = deseq2.override,
  heatmap.override = heatmap.override,
  volcano.override = volcano.override
)

# ---- Spatial Workflow ----

raw_metadata <- data.frame(Cell = c(""))
samples <- list.files(path = proj.params$filt_matrix_dir)
reference.samples <- NULL
bins <- c(8,16)

for (b in bins){
  
  spatial_samples <- c()
  assay <- dplyr::case_when(b == 2 ~ "Spatial.002um",
                            b == 8 ~ "Spatial.008um",
                            b == 16 ~ "Spatial.016um")
  for (s in samples){
    
    spatial_sample <- paste0(s, ".", assay)
    spatial_samples <- c(spatial_samples, spatial_sample)
    
    s.obj <- read_spaceranger(s, b, proj.params$filt_matrix_dir)
    # DoubletFinder and scDblFinder are designed to detect doublets based on 
    # expression profile similarities per droplet. In spatial transcriptomics
    # (Visium, Slide-seq, etc.), each “spot” or “capture location” contains RNA
    # from multiple cells by design.So, we cannot use them.
    s.obj <- calc_qc_metrics_sc_sp(s.obj, assay)
    s.obj <- mark_low_quality_sc_sp(s.obj)
    s.obj <- filter_singlets_sc_sp(s.obj)
    assign(paste0(s, ".", assay), s.obj)
  }
  filt.obj <- merge_filtered_sc_sp(spatial_samples, assay, proj.params, proj.params$seurat_dir)
}

# If integrating multiple samples, use below codes
# NOTE: Run workflow for each bin size separately
assay <- "Spatial.016um"  #"Spatial.008um"
filt.obj <- readRDS(file.path(proj.params$seurat_dir, paste0("filtered.seurat.", assay, ".rds")))

# (Optional) Remove bad samples 
#good_samples <- c("A06", "A07", "A10", "A11", "A12", "B07", "B10", "C07", "C09", "C10", "D08", "D09", "D10", "D12", "E08", "E09")
good_samples <- filt.obj@meta.data %>% 
  dplyr::count(Patient) %>% 
  dplyr::filter(n >= 500) %>%
  dplyr::pull(Patient)
filt.obj <- subset(filt.obj, grepl(pattern=paste0(good_samples, collapse="|"), x = filt.obj@meta.data$Patient))

sct.obj <- sctransform_sc_sp(filt.obj, assay, proj.params$seurat_dir)
integ <- integrate_sc_sp(sct.obj, assay, reference.samples, proj.params$seurat_dir)
integ.clust <- cluster_sc_sp(integ, assay, proj.params$seurat_dir)
integ.final <- remove_sparse_clusters_sc_sp(integ.clust, assay, proj.params$seurat_dir)
plot_metrics_post_integration_sc_sp(integ.final, assay, proj.params$diagnostics_dir)

reduction <- "Harmony"
for (res in c(0.4, 0.6, 0.8, 1, 1.2, 1.4)){
  identify_markers_sc_sp(integ.final, assay, res, reduction, proj.params$seurat_dir)
}
resolution <- 1
integ.final <- calc_module_scores(integ.final, assay, proj.params, proj.params$seurat_dir)
plot_module_scores(integ.final, assay, resolution, reduction, proj.params, proj.params$seurat_dir)

clusters <- list("Hepatocytes"              = c(4,9,11,12,15,20,23),
                 "Pancreatic.Acinar"        = c(3,24),
                 "Pancreatic.Islet"         = c(26),
                 "B.Plasma"                 = c(13),
                 "T.NK"                     = c(),
                 "Fibroblasts"              = c(7,8,16,19,22),
                 "Macrophages"              = c(14),
                 "Mast"                     = c(),
                 "Platelets"                = c(),
                 "Dendritic"                = c(),
                 "Endothelial"              = c(),
                 "Lymph.Endothelial"        = c(),
                 "Myocytes.Myofibroblasts"  = c(17),
                 "Neurons"                  = c(),
                 "Epithelial"               = c(5,10,2,21,25),
                 "Met.Panc"                 = c(1),
                 "Primary.Panc"             = c(6,18),
                 "Ductal.Cells"             = c(),
                 "Unclassified"             = c())
integ.ann <- annotate_manual_sc_sp(integ.final, assay, clusters, resolution, reduction, proj.params$seurat_dir)

# Visualize annotations in dim plot
subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
reduction <- "umap.harmony"
color.col <- "Cell.Type"
filename <- "UMAP_Annotated"
output_path <- proj.params$seurat_dir
split.col <- NULL
plot_umap(subset_seurat, reduction, color.col, filename, output_path, split.col)

# Visualize markers in dot plot
subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
assay <- "Spatial.016um"
ident.1 <- "Cell.Type"
ident.2 <- NULL
features <- find_top_features(subset_seurat , assay, proj.params)
filename <- "Dot_Plot_Markers"
output_path <- proj.params$seurat_dir
split.col <- NULL
plot_dot_custom(subset_seurat, assay, proj.params, ident.1, ident.2, features, filename, output_path, split.col)

# Visualize clusters from each column metadata
coi <- c("Sample", "Patient", "Condition", "Treatment", "Disease")
subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
reduction <- "umap.harmony"
color.col <- "Cell.Type"
output_path <- proj.params$seurat_dir
for (split.col in coi){
  if(!is.null(subset_seurat@meta.data[[split.col]] %>% unique())){
    filename <- paste0("UMAP_", split.col)
    plot_umap(subset_seurat, reduction, color.col, filename, output_path, split.col)
  }
}

# Tabulate population frequencies
split.cols <- c("Sample", "Patient", "Condition", "Treatment", "Disease")
subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
output_path <- proj.params$seurat_dir
tabulate_frequency(subset_seurat, split.cols, output_path)

# [OPTIONAL] Plot spatial map from each column metadata
coi <- c("Sample", "Patient", "Condition", "Treatment", "Disease")
subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
for (split.col in coi){
  if(!is.null(subset_seurat@meta.data[[split.col]] %>% unique())){
    
    filename <- paste0("Spatial_Plot_", split.col)
    p <- Seurat::SpatialDimPlot(object = subset_seurat, 
                                group.by = split.col,
                                image.scale="lowres",   # "hires" hides show H&E image
                                pt.size.factor = 4)  
    
    # Save combined plot
    ggsave(filename  = file.path(proj.params$seurat_dir, paste0(filename,".pdf")),
           plot      = p,
           limitsize = FALSE,
           bg        = "white")
  }
}
  
# [OPTIONAL] View Spatial annotation of slides
assay <- "Spatial.016um"
output_path <- proj.params$seurat_dir
plot.seurat <- Seurat::SplitObject(object = integ.ann,
                                   split.by = "Sample")
for(x in 1:length(plot.seurat)){

  # These limits were identified earlier based on repeated plotting
  x1 <- c(750, 940, 1250, 1565, 1885, 2225)
  y1 <- c(400, 750, 1065, 1385, 1700, 1900)
  x2 <- c(925, 1250, 1560, 1900, 2250)
  y2 <- c(375,  575,  875, 1175, 1475, 1775)

  plot_spatial_map(plot.seurat[[x]], x1, y1, x2, y2, assay, output_path)
}

# ---- Pseudobulk Workflow ----

assay <- "Spatial.016um"
integ.ann <- readRDS(file.path(proj.params$seurat_dir, paste0("integrated.seurat.", assay, ".ann.rds")))
subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
comparison.col <- "Comparisons"
subset_seurat@meta.data <- subset_seurat@meta.data %>% dplyr::mutate(Comparisons = paste0(Treatment, ".", Cell.Type))

# IMPORTANT: In spatial analysis, Sample corresponds to Slide.
# However, we want to aggregate counts at level of patients, not slide.
subset_seurat@meta.data <- subset_seurat@meta.data %>%
  dplyr::rename(Slide = Sample, 
                Sample = Patient)

p <- prep_pseudobulk(subset_seurat, comparison.col, assay)

meta_data <- p$meta_data
read_data <- p$read_data
trial <- FALSE
main_analysis(meta_data, read_data, proj.params, trial)

# ---- Seurat DEG Workflow ----

subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
comparison.col <- "Treatment"
subset.col <- "Cell.Type"
assay <- "Spatial.016um"
output_path <- proj.params$deseq2_dir
degs_df <- find_degs_seurat(subset_seurat, comparison.col, output_path, celltype.col, assay)

# ---- USE SECTION BELOW TO DETERMINE CO-ORDINATES FOR MULTI-SAMPLE SLIDES ----

# NOTE: Sometimes, multiple samples are present in each slide. So, we need to
# determine Xmin, Xmax, Ymin, Ymax for each sample in order to include them
# in the extra_metada.

# Sample metadata
Normal.Liver <- c("E06", "B07", "C07", "A09", "B09", "C10", "D10")
Normal.Pancreas <- c("A06", "D08", "E08", "A11", "B12", "C12", "D13")
Panc.Cancer <- c("A08", "B08", "C08", "C09", "D09", "E09", "E10", "A12", "A13", "B13", "C13")
Met.Panc.Cancer <- c("B06", "C06", "D06", "A07", "D07", "E07", "A10", "B10", "B11", "C11", "D11", "D12")

# Add sample metadata based to each seurat object corresponding to specific bin size
for (s in slides){
  for (b in bins){
    assay <- dplyr::case_when(b == 2 ~ "Spatial.002um",
                              b == 8 ~ "Spatial.008um",
                              b == 16 ~ "Spatial.016um")
    
    sample.seurat <- get(paste0(s, ".", assay))
    
    # These limits were identified earlier based on repeated plotting
    x1 <- c(750, 940, 1250, 1565, 1885, 2225)
    y1 <- c(400, 750, 1065, 1385, 1700, 1900)
    x2 <- c(925, 1250, 1560, 1900, 2250)
    y2 <- c(375,  575,  875, 1175, 1475, 1775)
  }
}

# ---- NOTES ----

# NOTE: Normalization is key step in single and spatial analysis. The authors
# of Seurat indicate SCTransform() MUST be done for each sample individually but
# this shouldnt be followed blindly. The main objective of SCTransform() is to
# normalize raw gene expression data while "removing technical noise" and 
# stabilizing variance across cells or spots.

# Consider, you have a visiumHD slide with multiple samples and each sample
# comes from a different tissue. Now, you have 2 such VisiumHD slides and you
# sequence each slide separately. Now, should we do SCTransform() at sample 
# level? or tissue level? or slide level?

# Since Sequencing was done per slide, slide ~ technical batch and you want to 
# model and normalize technical noise within each batch, SCTransform() MUST be
# done at slide level. 

# Since each sample is from a different tissue (e.g., liver vs kidney), their
# expression profiles are inherently different. If you normalize per-sample or 
# per-tissue, you WILL remove real biological differences.

# Decide which reduction gives (i) clear separation of clusters 
# (ii) clusters cells together i.e. B cells are not split between clusters.
# In all cases so far, RPCA and Harmony seemed to work best.