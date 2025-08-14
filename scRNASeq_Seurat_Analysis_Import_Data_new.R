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
source("/hpc/home/kailasamms/projects/scRNASeq/Custom.Functions.R.Analysis.R")

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
  s.obj <- filter_singlets_sc.sp(s.obj)
  assign(s, s.obj)
}

plot_qc(raw_metadata, seurat_results)
#xpectr::suppress_mw(generate_whitelist(filt, seurat_results))

# Merge all samples
# NOTE: seurat objects MUST have been loaded into R prior to this step
# Import any other meta data associated with data set
# NOTE: This xslx file should have column named "Unique_ID" whose values matches 
# with  column "Unique_ID" of seurat object's metadata.
extra_metadata <-  openxlsx::read.xlsx(xlsxFile = paste0(scripts_path, metafile)) %>%
  dplyr::select(everything(), -Comments)
filt.obj <- merge_filtered_sc.sp(samples, "RNA", extra_metadata, seurat_results)

# Perform SCTransformation
sct.obj  <- sctransform_sc.sp(filt.obj, "RNA", seurat_results)

# Integrate all samples
reference.samples <- NULL
kweight <- min(sct@meta.data %>% dplyr::count(Sample) %>% dplyr::select(n) %>% min()/2, 100) 
integ <- integrate_sc.sp(sct.obj, "RNA", reference.samples, kweight, seurat_results)
integ.clust <- cluster_sc.sp(integ, "RNA", seurat_results)
integ.final <- remove_sparse_clusters_sc.sp(integ.clust, "RNA", seurat_results)

suffix <- "Full"
resolution <- 0.8
reduction <- "Harmony"
identify_markers_sc.sp(integ.final, "RNA", resolution, reduction, suffix, seurat_results)
plot_metrics_post_integration_sc.sp(integ.final, suffix, diagnostics_path)

# Decide which reduction gives (i) clear separation of clusters 
# (ii) clusters cells together i.e. B cells are not split between clusters.
# In all cases so far, RPCA and Harmony seemed to work best.