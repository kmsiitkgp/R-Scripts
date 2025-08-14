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

# Create a list of slide names which will be added to each barcode.
# Since folder names correspond to slide name, we just use list.files()
slides <- list.files(path = filt_matrix_path) 
bins <- c(8,16)

for (s in slides){
  for (b in bins){
    assay <- dplyr::case_when(b == 2 ~ "Spatial.002um",
                              b == 8 ~ "Spatial.008um",
                              b == 16 ~ "Spatial.016um")
    
    s.obj <- read_spaceranger(s, b, filt_matrix_path)
    s.obj <- calc_qc_metrics_sc.sp(s.obj, assay)
    s.obj <- mark_low_quality_sc.sp(s.obj)
    s.obj <- filter_singlets_sc.sp(s.obj)
    assign(paste0(s, ".", assay), s.obj)
  }
}

# Save the filtered seurat object for each bin size separately
for (b in bins){
  
  # Define the variables "spatial_samples" where the seurat objects are stored
  spatial_samples <- c()
  for (s in slides){
    assay <- dplyr::case_when(b == 2 ~ "Spatial.002um",
                              b == 8 ~ "Spatial.008um",
                              b == 16 ~ "Spatial.016um")
    filename <- paste0(s, ".", assay)
    spatial_samples <- c(spatial_samples, filename)
  }
  
  # Merge all samples belonging to specific bin size
  # NOTE: seurat objects MUST have been loaded into R prior to this step
  extra_metadata <- read.xlsx(paste0(scripts_path, proj, "_Metadata.xlsx"))
  filt.obj <- merge_filtered_sc.sp(spatial_samples, assay, extra_metadata, seurat_results)
}

# If integrating multiple samples, use below codes
# NOTE: Run workflow for each bin size separately
for (assay in c("Spatial.008um", "Spatial.016um")){
  
  filt.obj <- readRDS(paste0(seurat_results, "filtered.seurat.", assay, ".rds"))
  
  # Remove bad samples
  good_samples <- c("A06", "A07", "A10", "A11", "A12", "B07", "B10", "C07", "C09", "C10", "D08", "D09", "D10", "D12", "E08", "E09")
  filt.obj <- subset(filt.obj, grepl(pattern=paste0(good_samples, collapse="|"), x=filt.obj@meta.data$Patient))
  filt.obj <- subset(filt.obj, subset = (nUMIs>=50 & nGenes>=25))
  
  # Perform SCTransformation at sample level, not slide level
  sctransform.var <- "Sample"
  sct.obj <- sctransform_sc.sp(filt.obj, assay, sctransform.var, seurat_results)
  
  # Integrate all samples belonging to specific bin size
  reference.samples <- NULL
  kweight <- min(sct.obj@meta.data %>% dplyr::count(Sample) %>% dplyr::select(n) %>% min()/2, 100) 
  integ <- integrate_sc.sp(sct.obj, assay, reference.samples, kweight, seurat_results)
  
  # Perform clustering at various resolution using different reductions
  integ.clust <- cluster_sc.sp(integ, assay, seurat_results)
  
  # Remove sparse clusters and save the results
  integ.final <- remove_sparse_clusters_sc.sp(integ.clust, assay, seurat_results)
  
  # Plot metrics post integration (DO this ONLY after adding sample info)
  plot_metrics_post_integration_sc.sp(integ.final, assay, diagnostics_path)
  
  # Find markers
  for (res in c(0.4, 0.8, 1, 1.4)){  
    resolution <- res
    reduction <- "harmony"
    identify_markers_sc.sp(integ.final, assay, resolution, reduction, assay, seurat_results)
  }
  
  # Plot spatial map of samples and groups
  Seurat::SpatialDimPlot(object = integ.final, 
                         group.by = "Patient",
                         image.scale="lowres",   # "hires" hides show H&E image
                         pt.size.factor = 4)
  
  # plot.seurat <- Seurat::SplitObject(object = integ.final,  
  #                                    split.by = "Sample")
  # for(x in 1:length(plot.seurat)){
  #   
  #   # These limits were identified earlier based on repeated plotting
  #   x1 <- c(750, 940, 1250, 1565, 1885, 2225)
  #   y1 <- c(400, 750, 1065, 1385, 1700, 1900)
  #   x2 <- c(925, 1250, 1560, 1900, 2250)
  #   y2 <- c(375,  575,  875, 1175, 1475, 1775)
  #   
  #   plot_spatial_map(plot.seurat[[x]], x1, y1, x2, y2, assay, diagnostics_path)
  # }
}

# Visualize by UMAP the contribution of each patient based on treatment
for (assay in c("Spatial.008um", "Spatial.016um")){
  
  integ.final <- readRDS(paste0(seurat_results, "integrated.seurat.", assay, ".rds"))
  
  # File names, reductions, splits  for each of the figures
  groups <- unique(integ.final@meta.data$Treatment)
  filenames <- paste0("UMAP.Treatment.", groups, ".", assay)
  reductions <- rep(x="umap.harmony", times=length(filenames))
  splits <- rep(x="Patient", times=length(filenames))
  
  for (i in 1:length(groups)){ 
    
    plot.seurat <- subset(integ.final, Treatment == groups[i])
    
    plot.seurat <- Seurat::SplitObject(object = plot.seurat,
                                       split.by = splits[i])
    
    purrr::map(.x = c(1:length(plot.seurat)),
               .f = function(x){  
                 Idents(plot.seurat[[x]]) <- "cluster.0.4.harmony"
                 Seurat::DimPlot(object = plot.seurat[[x]],
                                 reduction = reductions[i],
                                 group.by = "cluster.0.4.harmony",
                                 pt.size = 0.1,
                                 order = TRUE,  # plot positive cells above negative cells
                                 label = TRUE,
                                 raster = FALSE,
                                 combine = TRUE) +
                   NoLegend() +
                   my_theme + 
                   coord_cartesian(xlim = c(floor(min(integ.final@reductions[[reductions[i]]]@cell.embeddings[,1])), 
                                            ceiling(max(integ.final@reductions[[reductions[i]]]@cell.embeddings[,1]))), 
                                   ylim = c(floor(min(integ.final@reductions[[reductions[i]]]@cell.embeddings[,2])), 
                                            ceiling(max(integ.final@reductions[[reductions[i]]]@cell.embeddings[,2])))) +
                   ggplot2::labs(title = names(plot.seurat)[x]) 
               }) %>% cowplot::plot_grid(plotlist=.,
                                         align="hv",
                                         axis="tblr",
                                         nrow=ceiling(sqrt(length(plot.seurat))),
                                         ncol=floor(sqrt(length(plot.seurat))),
                                         rel_widths=1,
                                         rel_heights=1,
                                         greedy=TRUE,
                                         byrow=TRUE)
    # Save the plot
    ggplot2::ggsave(filename = paste0(filenames[i], ".tiff"),
                    plot = last_plot(),
                    device = "jpeg",
                    path = diagnostics_path,
                    scale = 1,
                    width = 4*floor(sqrt(length(plot.seurat))),
                    height = 4*ceiling(sqrt(length(plot.seurat))),
                    units = c("in"),
                    dpi = 600,
                    limitsize = TRUE,
                    bg = "white")
  }
}

# Decide which reduction gives (i) clear separation of clusters 
# (ii) clusters cells together i.e. B cells are not split between clusters.
# In all cases so far, RPCA and Harmony seemed to work best.

# Annotate clusters
# 16um bin, resolution 1.4
clusters <- list("Hepatocytes"         = c(2,6,9,11,13,15,19,23),
                 "Pancreatic.Acinar"   = c(1,20),
                 "Pancreatic.Islet"    = c(28),
                 "B.Plasma"            = c(14),
                 "T.NK"                = c(),
                 "Fibroblasts"         = c(4,16,22,25,27),
                 "Macrophages"         = c(7),
                 "Dendritic"           = c(),
                 "Endothelial"         = c(),
                 "Lymph.Endothelial"   = c(),
                 "Myocytes"            = c(),
                 "CAFs"                = c(17),
                 "Epithelial"          = c(),
                 "Neurons"             = c(),
                 "Epithelial.I"        = c(3,18,24),
                 "Epithelial.II"       = c(8,10),
                 "Epithelial.III"      = c(5,12,21),
                 "Epithelial.IV"       = c(26),
                 "Unclassified"        = c())

assay <- "Spatial.016um"
resolution <- 1.4
reduction <- "harmony"
integ.final <- readRDS(paste0(seurat_results, "integrated.seurat.", assay, ".rds"))
integ.annotated <- annotate_manual_sc.sp(integ.final, clusters, resolution, reduction, assay, seurat_results)


# QC [Project Specific]
# I notice there are few hepatocyptes in normal pancreas etc. So, we remove
# such wrongly annotated cells
metadata <- integ.annotated@meta.data
metadata <- metadata %>%
  dplyr::filter(!(Cell.Type == "Hepatocytes" & Treatment %in% c("Normal.Pancreas", "Tumor.Pancreas") | 
                    Cell.Type == "Pancreatic.Acinar" & Treatment %in% c("Normal.Liver") |
                    Cell.Type ==  "CAFs" & Treatment %in% c("Normal.Liver", "Normal.Pancreas")))
integ.annotated@meta.data <- metadata


# Plot final UMAP and also UMAP split by treatment
assay <- "Spatial.016um"
filenames <- c(paste0(c("UMAP.Treatment.", "UMAP.Final."), assay))
reductions <- rep(x="umap.harmony", times=length(filenames))
splits <- c("Treatment", NA)
integ.annotated <- readRDS(paste0(seurat_results, "integrated.seurat.", assay, ".ann.rds"))
for (i in 1:length(filenames)){ 
  
  # Keep color palette consistent
  c <- integ.annotated@meta.data$Cell.Type %>% unique()
  cols <- my_palette[1:length(c)]
  names(cols) <- c
  
  
  plot.seurat <- integ.annotated
  if (!is.na(splits[i])){
    plot.seurat <- Seurat::SplitObject(object = plot.seurat,
                                       split.by = splits[i])
  } else{
    plot.seurat <- list(Full = plot.seurat)
  }
  
  purrr::map(.x = c(1:length(plot.seurat)),
             .f = function(x){  
               Idents(plot.seurat[[x]]) <- "Cell.Type"
               Seurat::DimPlot(object = plot.seurat[[x]],
                               reduction = reductions[i],
                               group.by = "Cell.Type",
                               cols = cols,
                               pt.size = 0.1,
                               order = TRUE,  # plot positive cells above negative cells
                               label = FALSE,
                               raster = FALSE,
                               combine = TRUE) +
                 #NoLegend() +
                 my_theme + 
                 coord_cartesian(xlim = c(floor(min(integ.final@reductions[[reductions[i]]]@cell.embeddings[,1])), 
                                          ceiling(max(integ.final@reductions[[reductions[i]]]@cell.embeddings[,1]))), 
                                 ylim = c(floor(min(integ.final@reductions[[reductions[i]]]@cell.embeddings[,2])), 
                                          ceiling(max(integ.final@reductions[[reductions[i]]]@cell.embeddings[,2])))) +
                 ggplot2::labs(title = names(plot.seurat)[x]) 
             }) %>% cowplot::plot_grid(plotlist=.,
                                       align="hv",
                                       axis="tblr",
                                       nrow=ceiling(sqrt(length(plot.seurat))),
                                       ncol=floor(sqrt(length(plot.seurat))),
                                       rel_widths=1,
                                       rel_heights=1,
                                       greedy=TRUE,
                                       byrow=TRUE)
  # Save the plot
  ggplot2::ggsave(filename = paste0(filenames[i], ".tiff"),
                  plot = last_plot(),
                  device = "jpeg",
                  path = diagnostics_path,
                  scale = 1,
                  width = 7*floor(sqrt(length(plot.seurat))),
                  height = 7*ceiling(sqrt(length(plot.seurat))),
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = "white")
}

#****************DE ANALYSIS & HEATMAP USING ALL GROUPS*********************#

# We need to identify set of genes varying across groups before plotting heatmap.
# If we have just 2 groups, we can identify DEGs using Wald test with contrasts
# in DESeq2. To identify DEGs across multiple groups, we need to use LRT in
# DESeq2. Although, we can do pairwise DE analysis between all groups, it is 
# NOT RECOMMENDED.

assay <- "Spatial.016um"
integ.annotated <- readRDS(paste0(seurat_results, "integrated.seurat.", assay, ".ann.rds"))

################# Generate metadata from seurat object
meta_data <- integ.annotated@meta.data %>%
  dplyr::filter(stringr::str_detect(string=Cell.Type, pattern="Epithelial")) %>%
  dplyr::select(Cell.Type, Treatment, Patient) %>%
  dplyr::mutate(Sample.ID = paste(Cell.Type, Treatment, Patient, sep=".")) %>% #Patient,
  dplyr::add_count(Sample.ID) %>%
  dplyr::filter(n>=150) %>%
  dplyr::distinct_at(.vars="Sample.ID", .keep_all = TRUE)

################## Define DEG parameters
groups <- meta_data %>% dplyr::pull(Cell.Type) %>% unique()
# groups_combn <- utils::combn(x=groups, m=2)
# contrasts <- apply(X=x, MARGIN=2, FUN=function(x){paste(x, collapse="-")})
DEG.params  <- list(contrast    = c("Epithelial.III-Epithelial.II"),
                    design      = "Cell.Type",
                    design.ref  = c("Cell.Type:Epithelial.I"),
                    lfc.cutoff  = 0,
                    padj.cutoff = 0.1,
                    deseq2.batch.correct = FALSE,
                    proj        = "Spatial_Bhowmick",
                    species     = "Homo sapiens")

degs_seurat <- integ.annotated
degs_seurat@meta.data <- degs_seurat@meta.data %>%
  dplyr::mutate(Sample.ID = paste(Cell.Type, Treatment, Patient, sep="."))

degs_seurat <- subset(degs_seurat, 
                      Sample.ID %in% meta_data$Sample.ID)

################# Generate read data from seurat object
# NOTE: The read data will have "the reads of all cells belonging to a single
# sample" merged together in each column. First, create a list of samples
samples <- degs_seurat@meta.data %>%
  dplyr::select(Sample.ID) %>%
  unlist(., use.names=FALSE) %>%
  unique()

# Second, create an empty dataframe with rows=genes and columns=samples
read_data <- data.frame(matrix(NA, nrow=nrow(degs_seurat@assays[[assay]]$counts), ncol=length(samples)))
rownames(read_data) <- rownames(degs_seurat@assays[[assay]]$counts)
colnames(read_data) <- samples

# Thirdly, we will add row-wise, the counts of each gene for each sample
for(i in samples){
  
  # Create a list of cells for each sample
  cells_subset <- degs_seurat@meta.data %>% dplyr::filter(Sample.ID == i) %>% rownames()
  
  # Use data.frame to convert "." in sparse matrix to "0"
  subset <- data.frame(degs_seurat@assays[[assay]]$counts[,cells_subset])
  read_data[,i]  <- rowSums(subset)
}

read_data <- read_data[rowSums(read_data) != 0,]
read_data <- read_data %>%
  tibble::rownames_to_column("SYMBOL")

################## Format metadata & readdata
meta_data <- prep_metadata(meta_data, read_data)
read_data <- prep_readdata(read_data, meta_data)
l <- check_data(read_data, meta_data)
meta_data <- l[[2]]
read_data <- l[[1]]
# norm_counts <- norm_counts_DESeq2(meta_data, read_data, seurat_results)
# norm_counts <- norm_counts[, -c(2,3)]

################## Perform all pairwise DE analysis

# Perform DESeq2() using in-built batch modelling
# Create DESeq2 object with appropriate variables in design
dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                      colData=meta_data, 
                                      design=~1)
design(dds) <- as.formula(paste0("~", DEG.params$design))
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~1)
res <- results(dds_lrt)
summary(res)

# Get vst transformed counts for plotting heatmap
vst <- DESeq2::vst(dds, blind=FALSE) 
vst_counts <- SummarizedExperiment::assay(vst) %>%
  data.frame() %>%
  tibble::rownames_to_column("SYMBOL")

# Define heatmap parameters
heatmap.params <- list(anno.row       = NULL,        
                       anno.column    = c("Cell.Type", "Treatment", "Patient"),
                       row.split      = NA,     
                       col.split      = c("Cell.Type"), 
                       row.cluster    = c("all"),  # c("alphabetical", "group", "all")
                       col.cluster    = c("group"),  # c("alphabetical", "group", "all")
                       discrete_panel = TRUE, 
                       log.transform  = FALSE,  #vst_counts are already log transformed
                       scale          = TRUE,
                       border_color   = "white",
                       bar_width      = NA,              # NA , 5
                       bar_height     = NA,              # NA , 5
                       width          = 5,              # NA
                       height         = 5,              # NA 
                       angle          = 90,
                       matrix_color   = "rdbu",          # c("vrds", "rdbu")
                       expr_legend    = TRUE,  # set FALSE if it overlaps with annotation legends
                       file_format    = "tiff")

metadata_column <- meta_data
metadata_row <- NULL
plot_genes <- res %>% 
  data.frame() %>% 
  dplyr::filter(padj < 0.05, abs(log2FoldChange) >= 1) %>% 
  rownames()
disp_genes <- ""
file_suffix <- ""

plot_heatmap(vst_counts, metadata_column, metadata_row, heatmap.params,
             plot_genes, disp_genes, file_suffix, seurat_results)

assay <- "Spatial.016um"
Idents(degs_seurat) <- "Cell.Type"
DefaultAssay(degs_seurat) <- assay
de_markers <-  Seurat::FindAllMarkers(object=degs_seurat,
                                        assay=assay,
                                        logfc.threshold=0.25,
                                        test.use="wilcox",
                                        slot="data",
                                        min.pct=0.1,
                                        min.diff.pct= -Inf,
                                        only.pos=TRUE)

# # Find top markers for each major cluster
# top_markers <- de_markers %>%
#   tibble::rownames_to_column("SYMBOL") %>%
#   dplyr::filter(avg_log2FC >= 0.58 & p_val_adj < 0.05) %>%
#   dplyr::group_by(cluster) %>%
#   dplyr::arrange(desc(avg_log2FC)) %>%  #desc(ratio)
#   dplyr::slice_head(n=200) %>%
#   ungroup()
# plot_genes <- top_markers$SYMBOL

# Find top markers for each major cluster
de_markers <- de_markers %>%
  tibble::rownames_to_column("SYMBOL") %>%
  dplyr::filter(p_val_adj < 0.05)
plot_genes <- de_markers$SYMBOL

plot_heatmap(vst_counts, metadata_column, metadata_row, heatmap.params,
             plot_genes, disp_genes, file_suffix, seurat_results)

### Pseudo bulk DESeq2

# Define DEG parameters
DEG.params  <- list(contrast    = c("Epithelial.III-Epithelial.II"),
                    design      = "Cell.Type",
                    design.ref  = c("Cell.Type:Epithelial.I"),
                    lfc.cutoff  = 0,
                    padj.cutoff = 0.1,
                    deseq2.batch.correct = FALSE,
                    proj        = "Spatial_Bhowmick",
                    species     = "Homo sapiens")


# Run DESeq2 workflow
annotations <- get_annotations()
meta_data <- prep_metadata(meta_data, read_data)
read_data <- prep_readdata(read_data, meta_data)
l <- check_data(read_data, meta_data)
meta_data <- l[[2]]
read_data <- l[[1]]
formula_string <- paste0("~", DEG.params$design)
dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                      colData=meta_data, 
                                      design=~1)
design(dds) <- as.formula(formula_string)
n <- 1
approach <- ""
dds.list <- run_deseq2(dds, meta_data, DEG.params, n, approach, seurat_results)



#*****************************************************************************#


#*****[USE SECTION BELOW TO DETERMINE CO-ORDINATES FOR MULTI-SAMPLE SLIDES]****#

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
