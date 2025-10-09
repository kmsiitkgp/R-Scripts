#!/usr/bin/env Rscript

# Read and store variables from command line interface (CLI)
cli <- base::commandArgs(trailingOnly = TRUE) 
args <- base::strsplit(x = cli, split = "=", fixed = TRUE)

for (e in args){
  argname <- e[1]
  argval <- e[2]
  assign(argname, argval)
}

#source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/Custom_Functions.R")
source("/hpc/home/kailasamms/projects/scRNASeq/Custom_Functions.R")

# ---- Project Directories ----

proj.dir     <- file.path("/hpc/home/kailasamms/scratch", proj)
scripts.dir  <- "/hpc/home/kailasamms/projects/scRNASeq"

proj.params <- list(
  meta.file               = file.path(scripts.dir, paste0(proj, "_Metadata.xlsx")),
  cell.cycle.marker.file  = file.path(scripts.dir, "Cell_Cycle_Markers.xlsx"),
  cell.type.marker.file   = file.path(scripts.dir, "Cell_Type_Markers.xlsx"),
  raw_matrix_dir          = file.path(proj.dir, "raw_feature_bc_matrix"),
  filt_matrix_dir         = file.path(proj.dir, "filt_feature_bc_matrix"),
  hto_matrix_dir          = file.path(proj.dir, "raw_hto_bc_matrix"),
  diagnostics_dir         = file.path(proj.dir, "diagnostics"),
  seurat_dir              = file.path(proj.dir, "seurat"),
  deseq2_dir              = file.path(proj.dir, "deseq2"),
  demux_dir               = file.path(proj.dir, "demux"),
  pyscenic_dir            = file.path(proj.dir, "pyscenic"),
  scvelo_dir              = file.path(proj.dir, "scvelo"),
  velocyto_dir            = file.path(proj.dir, "velocyto"),
  cellphonedb_dir         = file.path(proj.dir, "cellphonedb"),
  cellchat_dir            = file.path(proj.dir, "cellchat"),
  dorothea_dir            = file.path(proj.dir, "dorothea"),
  
  proj                    = proj,                   # Project name
  species                 = species,                # Species name (e.g. "Mus musculus", "Homo sapiens")
  
  counts_dir              = file.path(proj.dir, "counts"),     # Directory containing count files
  deseq2_dir              = file.path(proj.dir, "deseq2"),     # Directory to store DESeq2 results
  contrast                = contrasts,                         # Vector of contrasts for DE analysis
  deseq2.design           = "Comparisons",                     # DESeq2 design formula or column name
  deseq2.lfc.cutoff       = 0,                                 # Log fold change cutoff for significance
  deseq2.padj.cutoff      = 0.1,                               # Adjusted p-value cutoff for significance
  deseq2.batch.correct    = FALSE,                             # Boolean, whether to apply batch correction
  
  heatmap.force.log       = FALSE,                       # Force log transform on heatmap data (default FALSE, auto detect)
  heatmap.col.ann         = c("Treatment"),              # Columns from metadata used as column annotation
  heatmap.row.ann         = NULL,                        # Columns from metadata_row for row annotation
  heatmap.col.gaps        = NULL,                        # Columns to define gaps in heatmap columns
  heatmap.row.gaps        = NULL,                        # Columns to define gaps in heatmap rows
  heatmap.col.cluster     = c("Treatment"),              # Clustering for columns ("all", "alphabetical", or metadata column)
  heatmap.row.cluster     = "all",                       # Clustering for rows ("all", "alphabetical", or metadata column)
  heatmap.palette         = "rdbu",                      # Color palette for heatmap ("rdbu" or "vrds")
  heatmap.ann.palette     = "discrete",                  # Annotation palette type ("discrete" or "sequential")
  heatmap.border.color    = NA,                          # Border color of heatmap cells, NA for no border
  heatmap.show.expr.legend= TRUE,                        # Show expression legend on heatmap (set FALSE if overlapping annotations)
  heatmap.title           = NA,                          # Title for heatmap (default NA = no title)
  heatmap.format          = "tiff",                      # Output format for heatmap ("tiff", "jpeg", etc.)
  
  volcano.lfc.cutoff      = 0.58,                         # Log fold change cutoff for volcano plot
  volcano.padj.cutoff     = 0.05,                         # Adjusted p-value cutoff for volcano plot
  volcano.color           = "vrds",                       # Color palette for volcano plot ("vrds", etc.)
  volcano.label.genes     = NULL                          # Optional vector of genes to label on volcano plot
)

# ---- Global Options ----

options(future.globals.maxSize = 1e15)            
options(Seurat.object.assay.version = "v5")       
set.seed(1234)

# ---- Load Marker Files ----

cell_type_markers  <- openxlsx::read.xlsx(xlsxFile = proj.params$cell.type.marker.file)
cell_cycle_markers <- openxlsx::read.xlsx(xlsxFile = proj.params$cell.cycle.marker.file)

# ---- Cell Cycle Genes ----

cell_cycle_genes <- c(
  cell_cycle_markers$Human_Gene, 
  cell_cycle_markers$Mouse_Gene)

s_genes <- c(
  cell_cycle_markers$Human_Gene[cell_cycle_markers$Phase == "S"],
  cell_cycle_markers$Mouse_Gene[cell_cycle_markers$Phase == "S"])

g2m_genes <- c(
  cell_cycle_markers$Human_Gene[cell_cycle_markers$Phase == "G2/M"],
  cell_cycle_markers$Mouse_Gene[cell_cycle_markers$Phase == "G2/M"])

# ---- Single cell Workflow ----

raw_metadata <- data.frame(Cell = c(""))
samples <- list.files(path = proj.params$filt_matrix_dir)
reference.samples <- NULL

for (s in samples){
  
  s.obj <- read_cellranger(s, proj.params$raw_matrix_dir)
  s.obj <- mark_empty_droplets_dropletutils(s.obj)
  s.obj <- mark_empty_droplets_cellranger(s.obj, proj.params$filt_matrix_dir)
  s.obj <- doublet_finder(s.obj)
  s.obj <- scdbl_finder(s.obj)
  s.obj <- calc_qc_metrics_sc_sp(s.obj, "RNA")
  s.obj <- mark_low_quality_sc_sp(s.obj)
  raw_metadata <- generate_plotdata(s.obj, raw_metadata)
  s.obj <- filter_singlets_sc_sp(s.obj)
  assign(s, s.obj)
}

plot_qc(raw_metadata, proj.params$seurat_dir)

assay <- "RNA"
filt.obj <- merge_filtered_sc_sp(samples, assay, proj.params, proj.params$seurat_dir)
# filt.obj <- readRDS(file.path(proj.params$seurat_dir, paste0("filtered.seurat.", assay, ".rds")))
# scDblFinder identified ~16K doublets marked as singlets by DoubletFinder
# filt.obj <- subset(filt.obj, scDblFinder == "Singlet")
sct.obj <- sctransform_sc_sp(filt.obj, assay, proj.params$seurat_dir)
integ <- integrate_sc_sp(sct.obj, assay, reference.samples, proj.params$seurat_dir)
integ.clust <- cluster_sc_sp(integ, assay, proj.params$seurat_dir)
integ.final <- remove_sparse_clusters_sc_sp(integ.clust, assay, proj.params$seurat_dir)
plot_metrics_post_integration_sc_sp(integ.final, assay, proj.params$diagnostics_dir)

reduction <- "Harmony"
for (res in c(0.4, 0.6, 0.8, 1, 1.2, 1.4)){
  identify_markers_sc_sp(integ.final, assay, res, reduction, proj.params$seurat_dir)
}

resolution <- 0.8
integ.final <- calc_module_scores(integ.final, assay, proj.params, proj.params$seurat_dir)
plot_module_scores(integ.final, assay, resolution, reduction, proj.params, proj.params$seurat_dir)

clusters <- list("Hepatocytes"              = c(),
                 "Pancreatic.Acinar"        = c(),
                 "Pancreatic.Islet"         = c(),
                 "B.Plasma"                 = c(),
                 "T.NK"                     = c(),
                 "Fibroblasts"              = c(),
                 "Macrophages"              = c(),
                 "Mast"                     = c(),
                 "Platelets"                = c(),
                 "Dendritic"                = c(),
                 "Endothelial"              = c(),
                 "Lymph.Endothelial"        = c(),
                 "Myocytes.Myofibroblasts"  = c(),
                 "Neurons"                  = c(),
                 "Epithelial"               = c(),
                 "Met.Panc"                 = c(),
                 "Primary.Panc"             = c(),
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
#features <- c("GH1")
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

# [OPTIONAL] Visualize features in feature plot
subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
reduction <- "umap.harmony"
features <- c("IGF1", rownames(integ.ann@assays$RNA)[grepl("^GH.*", rownames(integ.ann@assays$RNA))])  
filename <- "Feature_Plot"
output_path <- proj.params$seurat_dir
split.col <- NULL #"Groups" #"Patient"
plot_features(subset_seurat, features, reduction, filename, output_path, split.col)

# ---- Pseudobulk Workflow ----

assay <- "RNA"
integ.ann <- readRDS(file.path(proj.params$seurat_dir, paste0("integrated.seurat.", assay, ".ann.rds")))
subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
comparison.col <- "Comparisons"
subset_seurat@meta.data <- subset_seurat@meta.data %>% dplyr::mutate(Comparisons = paste0(Treatment, ".", Cell.Type))
p <- prep_pseudobulk(subset_seurat, comparison.col, assay)

contrasts <- c("Tumor.Fibroblasts-Normal.Fibroblasts",
               "Tumor.T.NK-Normal.T.NK",
               "Tumor.Epithelial-Normal.Epithelial",
               "Tumor.Endothelial-Normal.Endothelial",
               "Tumor.Myocytes.Myofibroblasts-Normal.Myocytes.Myofibroblasts",
               "Tumor.Mast-Normal.Mast",
               "Tumor.Macrophages-Normal.Macrophages",
               "Tumor.Dendritic-Normal.Dendritic",
               "Tumor.B.Plasma-Normal.B.Plasma")  
proj.params$contrast <- contrasts

main_analysis(proj.dir, p$meta_data, p$read_data, proj.params)

# ---- Seurat DEG Workflow ----

subset_seurat <- subset(integ.ann, (Cell.Type == "Unclassified"), invert=TRUE)
comparison.col <- "Treatment"
subset.col <- "Cell.Type"
assay <- "RNA"
output_path <- proj.params$deseq2_dir
degs_df <- find_degs_seurat(subset_seurat, comparison.col, output_path, celltype.col, assay)
