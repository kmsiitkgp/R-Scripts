#!/usr/bin/env Rscript

# Read and store variables from command line interface (CLI)
cli <- base::commandArgs(trailingOnly = TRUE) 
args <- base::strsplit(x = cli, split = "=", fixed = TRUE)

for (e in args){
  argname <- e[1]
  argval <- e[2]
  assign(argname, argval)
}

# NOTE: All variables and functions are defined within the file below
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables_new.R")

#******************************************************************************#

# Create a list of sample names which will be added to each barcode.
# Since folder names correspond to sample name, we just use list.files()
samples <- list.files(path = raw_matrix_path) 

# Create empty dataframe to store raw metadata
raw_metadata <- data.frame(Cell = c("")) #Sample = as.factor(1), nUMIs = c(0)

for (s in samples){
  s.obj <- read_cellranger(s, raw_matrix_path)
  s.obj <- mark_emptydroplets_dropletutils(s.obj)
  s.obj <- mark_emptydroplets_cellranger(s.obj)
  s.obj <- doublet_finder(s.obj)
  s.obj <- scdbl_finder(s.obj)
  s.obj <- calc_qc_metrics(s.obj)
  s.obj <- mark_low_quality(s.obj)
  raw_metadata <- generate_plotdata(s.obj, raw_metadata)
  s.obj <- filter_singlets(s.obj)
  assign(s, s.obj)
}

plot_qc(raw_metadata, seurat_results)
#xpectr::suppress_mw(generate_whitelist(filt, seurat_results))
filt              <- format_filtered(samples, seurat_results)
sct               <- sctransform_singlecell(filt, seurat_results)
reference.samples <- NULL
kweight <- min(sct@meta.data %>% dplyr::count(Sample) %>% dplyr::select(n) %>% min()/2, 100) 
integ             <- integrate_singlecell(sct, reference.samples, kweight, seurat_results)
clust             <- cluster_singlecell(integ, seurat_results)
final             <- remove_sparse_clusters(clust, seurat_results)
suffix <- "Full"
plot_metrics_post_integration(final, suffix, diagnostics_path)
resolution <- 0.8
reduction <- "Harmony"
species <- "Mus musculus"
identify_markers(final, resolution, reduction, species, suffix, seurat_results)

# Decide which reduction gives (i) clear separation of clusters 
# (ii) clusters cells together i.e. B cells are not split between clusters.
# In all cases so far, RPCA and Harmony seemed to work best.