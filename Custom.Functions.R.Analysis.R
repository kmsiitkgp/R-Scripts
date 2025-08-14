#!/usr/bin/env Rscript

#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

pkgs <- c("BiocManager", "remotes", "AnnotationHub", "ensembldb", "org.Hs.eg.db",
          "org.Mm.eg.db", "fgsea", "clusterProfiler", "progeny", "dorothea", 
          "viper", "DESeq2", "sva", "GSVA", "RcisTarget", "glmGamPoi", "Seurat",
          "harmony", "hdf5r", "arrow", "leidenAlg", "scCustomize", "reticulate", 
          "ashr", "infercnv", 
          "UCell", "scDblFinder", "DropletUtils", "batchelor", 
          "ClusterFoldSimilarity", "CellChat", 
          "Banksy", "SeuratWrappers", "presto", "DoubletFinder", "SeuratData", 
          "oligo",
          "oligoData", "illuminaHumanv4.db", "hgu133plus2.db", "GEOquery", 
          "affy", "lumi", "openxlsx", "dplyr", "tibble", "stringr", "purrr", 
          "ggplot2", "ggplotify", "ggrepel", "ggpubr", "ggfortify", "ggridges",
          "ggbeeswarm", "pheatmap", "VennDiagram", "survival", "survminer", 
          "UpSetR", "umap", "plot3D", "cowplot", "viridis", "RColorBrewer", 
          "colorspace", "enrichplot", "ComplexHeatmap")

for (pkg in pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    message(paste("Loaded", pkg))
  } else {
    message(paste("Package", pkg, "is not installed — skipping"))
  }
}

# NOTE:  survminer handles %++% while dplyr handles %>%

#******************************************************************************#
#                       DEFINE GLOBAL OPTIONS & VARIABLES                      #
#******************************************************************************#

# ------------------------------
# General Configuration
# ------------------------------

options(future.globals.maxSize = 1e15)            # Increase object size limit
options(Seurat.object.assay.version = "v5")       # Set Seurat object version
# options(scipen = 999)                           # Avoid scientific notation (not recommended)

# ------------------------------
# Project Metadata & File Paths
# ------------------------------

# NOTE: 'proj' must be defined by the parent script
# Metadata Excel file should be named <proj>_Metadata.xlsx
# Must contain a column "Unique_ID" matching Seurat object metadata

metafile <- paste0(proj, "_Metadata.xlsx")

# Define paths relative to project
scripts_path        <- "/hpc/home/kailasamms/projects/scRNASeq/"
parent_path         <- paste0("/hpc/home/kailasamms/scratch/", proj, "/")

filt_matrix_path    <- paste0(parent_path, "filt_feature_bc_matrix/")
raw_matrix_path     <- paste0(parent_path, "raw_feature_bc_matrix/")
hto_matrix_path     <- paste0(parent_path, "raw_hto_bc_matrix/")
diagnostics_path    <- paste0(parent_path, "diagnostics/")
demux_results       <- paste0(parent_path, "results_demux/")
seurat_results      <- paste0(parent_path, "results_seurat/")
pyscenic_results    <- paste0(parent_path, "results_pyscenic/")
scvelo_results      <- paste0(parent_path, "results_scvelo/")
velocyto_results    <- paste0(parent_path, "results_velocyto/")
cellphonedb_results <- paste0(parent_path, "results_cellphonedb/")
cellchat_results    <- paste0(parent_path, "results_cellchat/")
dorothea_results    <- paste0(parent_path, "results_dorothea/")

# ------------------------------
# Cell Cycle Gene Markers
# ------------------------------

# Load cell cycle genes
cell_cycle_markers <- openxlsx::read.xlsx(xlsxFile = paste0(scripts_path, "Cell_Cycle_Markers.xlsx"))

cell_cycle_genes <- c(cell_cycle_markers$Human_Gene, cell_cycle_markers$Mouse_Gene)

s_genes <- c(cell_cycle_markers$Human_Gene[cell_cycle_markers$Phase == "S"],
             cell_cycle_markers$Mouse_Gene[cell_cycle_markers$Phase == "S"])

g2m_genes <- c(cell_cycle_markers$Human_Gene[cell_cycle_markers$Phase == "G2/M"],
               cell_cycle_markers$Mouse_Gene[cell_cycle_markers$Phase == "G2/M"])

# ------------------------------
# Global ggplot Theme
# ------------------------------

my_theme <- ggplot2::theme(
  plot.title    = element_text(family = "sans", face = "bold",  colour = "black", size = 15, hjust = 0.5),
  legend.title  = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0,   vjust = 1, angle = 0),
  axis.title.x  = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0.5, vjust = 0, angle = 0),
  axis.title.y  = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0.5, vjust = 1, angle = 90),
  legend.text   = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5),
  axis.text.x   = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5, vjust = 0.5, angle = 45),
  axis.text.y   = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5, vjust = 0.5, angle = 0))
# strip.text.x=element_text(family="sans", face="bold",  colour="black", size=10, hjust=0.5),
# legend.background=element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
# legend.position="right",
# legend.justification="left",
# legend.direction="vertical",
# legend.key.height=unit(0.5, 'cm'),
# legend.key.width =unit(0.5, 'cm'), 
# legend.text.align=0)

# ------------------------------
# Custom Color Palette
# ------------------------------

# Base color palette (up to 15 classes)
my_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#8C564B",
                "#E377C2", "#BCBD22", "#17BECF", "#FFC61E", "#762A83",
                "#333333", "#FF1F5B", "#B8E80C", "#9b19f5", "#DC0AB4")

scanpy_default_102 <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
  "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
  "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
  "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#e6550d", "#31a354",
  "#756bb1", "#636363", "#9e9ac8", "#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4",
  "#fed9a6", "#ffffcc", "#e5d8bd", "#fddaec", "#f2f2f2", "#a6cee3", "#1f78b4",
  "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
  "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072",
  "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5",
  "#ffed6f", "#a65628", "#ffff99", "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
  "#66a61e", "#e6ab02", "#a6761d", "#666666", "#e41a1c", "#377eb8", "#4daf4a",
  "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999", "#66c2a5",
  "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3",
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d",
  "#666666")

# Add transparency variants
my_palette <- c(my_palette,
                colorspace::adjust_transparency(col = my_palette, alpha = 0.2),
                colorspace::adjust_transparency(col = my_palette, alpha = 0.4),
                colorspace::adjust_transparency(col = my_palette, alpha = 0.6),
                colorspace::adjust_transparency(col = my_palette, alpha = 0.8))

# Optional: Alternate palette (commented)
# my_palette <- c("#000000", "#D9D9D9", "#003C30", "#beb9db", "#1D91C0", "#A6CEE3",
#                 "#50E991", "#A6D854", "#74C476", "#C7E9B4", "#00bfa0", "#E5F5F9",
#                 "#EDBF33", "#E6D800", "#FFF7BC", "#ffee65", "#C7EAE5", "#67001F",
#                 "#CB181D", "#FD8D3C", "#FC9272", "#EF3B2C", "#F16913", "#6A51A3",
#                 "#762A83", "#D4B9DA", "#0bb4ff", "#E60049", "#AE017E", "#DF65B0",
#                 "#FDCCE5", "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
#                 "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5", "#C51B7D",
#                 "#DE77AE", "#7F7F7F", "#9467BD")

#******************************************************************************#
#               SINGLE CELL & SPATIAL ANALYSIS RELATED FUNCTIONS               #       
#******************************************************************************#

#' Import Single cell transcriptomic data from Cell Ranger output
#'
#' Reads a Cell Ranger feature-barcode matrix and creates a Seurat object.
#' 
#' @param sample Character. Sample name.
#' @param path Character. Path to directory containing the feature-barcode matrix.
#'
#' @return A Seurat object containing raw/filtered count data.
#' @export
read_cellranger <- function(sample, path){
  
  # Construct the full data directory path
  data_dir <- base::file.path(path, sample)
  
  # Check if path exists
  if (!dir.exists(data_dir)) {
    stop("The directory '", data_dir, "' does not exist!")
  }
  
  # Read the feature-barcode matrix with gene symbols
  # gene.column = 1 → Ensembl IDs
  # gene.column = 2 → Gene symbols (required for mito/ribo/heme ratios)
  counts <- tryCatch({
    Seurat::Read10X(data.dir = paste0(path, sample),
                    gene.column = 2,
                    cell.column = 1,
                    unique.features = TRUE,
                    strip.suffix = FALSE)
  }, error = function(e) {
    stop("Failed to read 10X matrix from '", data_dir, "': ", e$message)
  })
  
  # Create Seurat object with raw matrix
  # Set min.cells = 0 and min.features = 0 to retain all barcodes for EmptyDrops()
  sample_seurat <- SeuratObject::CreateSeuratObject(counts = counts,
                                                    project = sample,
                                                    assay = "RNA",
                                                    names.field = 1,
                                                    names.delim = "_",
                                                    meta.data = NULL,
                                                    min.cells = 0,
                                                    min.features = 0)
  
  # Logging matrix type based on path
  if (grepl("raw", path, ignore.case = TRUE)){
    message("Raw Feature Barcode Matrix imported for '", sample, "'.")
  } else if (grepl("filt", path, ignore.case = TRUE)){
    message("Filtered Feature Barcode Matrix imported for '", sample, "'.")
  } else {
    message("Feature-barcode matrix (raw/filtered) imported for '", sample, "'.")
  }
  
  return(invisible(sample_seurat))
}

#' Import Spatial Transcriptomics data from SpaceRanger output
#'
#' Reads a Space Ranger feature-barcode matrix and creates a Seurat object.
#' 
#' @param sample Character. Sample name.
#' @param bin Numeric. Bin size in micrometers (2,8,16).
#' @param path Character. Path to directory containing the feature-barcode matrix.
#'
#' @return A Seurat object with Spatial assay containing filtered count data
#' @export
read_spaceranger <- function(sample, bin, path){
  
  # Construct the full data directory path
  data_dir <- base::file.path(path, sample)
  
  # Check if path exists
  if (!dir.exists(data_dir)) {
    stop("The directory '", data_dir, "' does not exist!")
  }
  
  # Define expected matrix file path
  matrix_file <- base::file.path(data_dir, "filtered_feature_bc_matrix.h5")
  
  # Check that matrix file exists
  if (!file.exists(matrix_file)) {
    stop("Matrix file 'filtered_feature_bc_matrix.h5' not found in: ", data_dir)
  }
  # Create Seurat object with filtered matrix
  sample_seurat <- Seurat::Load10X_Spatial(data.dir = data_dir,
                                           filename = "filtered_feature_bc_matrix.h5",
                                           assay = "Spatial",
                                           slice = sample,
                                           bin.size = bin,
                                           filter.matrix = TRUE,
                                           to.upper = FALSE,
                                           image = NULL)
  
  # Add orig.ident manually as Load10X_Spatial() doesn't set project/sample name
  sample_seurat@meta.data <- sample_seurat@meta.data %>% 
    dplyr::mutate(orig.ident = sample)
  
  message("Filtered Barcode Matrix imported for", sample, "at bin size", bin, "µm.")
  
  return(invisible(sample_seurat))
}

#' Identify empty droplets using DropletUtils
#'
#' @param sample_seurat Seurat object of a single sample from raw count matrix.
#'
#' @return Seurat object with 'DropletUtils' column ("Non-Empty Droplet"/"Empty Droplet") added to metadata.
#' @export
mark_empty_droplets_dropletutils <- function(sample_seurat){ 
  
  # Convert to SingleCellExperiment object
  sce <- Seurat::as.SingleCellExperiment(x = sample_seurat)
  
  # Set reproducible seed for emptyDrops
  set.seed(100)
  
  # NOTE: If FDR > 0.05 for some droplets AND Limited == TRUE, it indicates that
  # with more iterations, the FDR of these droplets can be reduced.
  
  # Iteratively run emptyDrops until limited droplets with high FDR are resolved
  niters <- 10000
  n_improve <- 1  # trigger loop
  
  while (n_improve > 0){
    e.out <- DropletUtils::emptyDrops(m = SingleCellExperiment::counts(sce),
                                      niters = niters)
    
    df <- as.data.frame(e.out)
    n_improve <- nrow(df %>%
                        dplyr::filter(Limited == TRUE, FDR > 0.05))
    message("emptyDrops check: ", n_improve, " droplets need more iterations (FDR > 0.05 & Limited=TRUE); niters = ", niters)
    niters <- niters + 10000
  }
  
  # Identify true (non-empty) cells with FDR ≤ 0.05
  true_cells <- df %>%
    dplyr::filter(FDR <= 0.05) %>% 
    rownames()
  
  # Annotate metadata with droplet classification
  sample_seurat@meta.data <- sample_seurat@meta.data %>%
    dplyr::mutate(Cell = rownames(.)) %>%
    dplyr::mutate(DropletUtils = dplyr::case_when(Cell %in% true_cells ~ "Non-Empty Droplet",
                                                  TRUE ~ "Empty Droplet"))
  
  # Log output
  ident <- as.character(unique(sample_seurat@meta.data$orig.ident))
  message("DropletUtils classification complete for sample: '", ident, "'")
  
  return(invisible(sample_seurat))
}


#' Mark empty droplets based on CellRanger-filtered matrix
#'
#' @param sample_seurat Seurat object of a single sample after mark_empty_droplets_dropletutils()
#' @param filt_matrix_path Character. Path to the filtered matrix directory
#'
#' @return Seurat object with 'CellRanger' column ("Non-Empty Droplet"/"Empty Droplet") added to metadata
#' @export
mark_empty_droplets_cellranger <- function(sample_seurat, filt_matrix_path){
  
  # Get sample name
  sample <- sample_seurat@meta.data$orig.ident %>% unique() %>% as.character()
  
  # Read the filtered barcode-feature matrix output of CellRanger
  sample_seurat_filt <- read_cellranger(sample, filt_matrix_path)
  
  # Mark cells absent in filtered barcode-feature matrix as empty droplets
  sample_seurat@meta.data <- sample_seurat@meta.data %>%
    dplyr::mutate(Cell = rownames(.)) %>%
    dplyr::mutate(CellRanger = dplyr::case_when(Cell %in% colnames(sample_seurat_filt) ~ "Non-Empty Droplet",
                                                TRUE ~ "Empty Droplet"))
  
  # Log output
  message("CellRanger empty droplets identified for sample: '", sample, "'")
  
  return(sample_seurat)
}

#' Identify doublets using DoubletFinder
#'
#' @param sample_seurat Seurat object of a single sample after mark_empty_droplets_cellranger()
#'
#' @return Seurat object with 'DoubletFinder' column ("Singlet"/"Doublet") added to metadata
#' @export
doublet_finder <- function(sample_seurat){
  
  # Filter out empty droplets before doublet identification
  subset_seurat <- subset(x = sample_seurat,
                          subset = (DropletUtils == "Empty Droplet" & CellRanger == "Empty Droplet"),
                          invert = TRUE)
  
  # Preprocess
  subset_seurat <- subset_seurat %>% Seurat::NormalizeData() %>% 
    Seurat::FindVariableFeatures() %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA()
  
  # Find significant PCs
  stdev_pc <- subset_seurat@reductions$pca@stdev
  percent_stdev_pc <- (stdev_pc / sum(stdev_pc)) * 100
  cumulative_stdev_pc <- cumsum(percent_stdev_pc)
  
  pc1 <- which(cumulative_stdev_pc > 90 & percent_stdev_pc < 5)[1]
  pc2 <- sort(which((percent_stdev_pc[1:(length(percent_stdev_pc)-1)] - 
                       percent_stdev_pc[2:length(percent_stdev_pc)]) > 0.1),
              decreasing = TRUE)[1] + 1
  n_pcs <- min(pc1, pc2)
  
  # pK Identification (no prior info)
  # Introduces artificial doublets in varying proportions into real dataset,
  # preprocesses the data and calculates proportion of artificial nearest 
  # neighbors. Output is a list of proportions of artificial nearest neighbors
  # for varying combinations of pK and pN. Optimal pK is the max of bimodality
  # coefficient (BCmvn) distribution
  
  # Run UMAP, neighbors, clustering for DoubletFinder paramSweep
  subset_seurat <- subset_seurat %>%
    Seurat::RunUMAP(dims = 1:n_pcs) %>%
    Seurat::FindNeighbors(dims = 1:n_pcs) %>%
    Seurat::FindClusters(resolution = 0.1)
  
  # Find optimal pK for DoubletFinder
  sweep.res <- DoubletFinder::paramSweep(subset_seurat, PCs = 1:n_pcs, sct = FALSE)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res, GT=FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  
  optimal_pK <- bcmvn %>% 
    dplyr::slice_max(order_by = BCmetric, n=1) %>%
    dplyr::pull(pK) %>%
    as.numeric()
  #optimal_pK <- as.numeric(as.character(optimal_pK[[1]]))
  
  # Estimate expected homotypic doublets based on 10X multiplet rate
  # https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled
  # From the above link, we can see that the multiplet rate is 8*10^-6 per cell
  multiplet_rate <- 8e-6
  n_cells <- nrow(subset_seurat@meta.data)
  n_exp <- round(multiplet_rate * n_cells)
  
  # Adjust for homotypic doublets
  homotypic.prop <- DoubletFinder::modelHomotypic(subset_seurat@meta.data$seurat_clusters)
  n_exp_adj <- round(n_exp * (1 - homotypic.prop))
  
  # Run DoubletFinder
  subset_seurat <- DoubletFinder::doubletFinder(seu = subset_seurat,
                                                PCs = 1:n_pcs,
                                                pN = 0.25,  #default
                                                pK = optimal_pK,
                                                nExp = nExp_poi.adj)
  
  # Rename classification column to 'DoubletFinder'
  colnames(subset_seurat@meta.data)[grepl(pattern = "DF.classifications", x = colnames(subset_seurat@meta.data))] <- "DoubletFinder"
  
  # Prepare classification dataframe for merging
  doublet_df <- subset_seurat@meta.data %>% 
    dplyr::select(DoubletFinder) %>%
    tibble::rownames_to_column(var = "Cell") %>%
    dplyr::mutate(DoubletFinder = stringr::str_to_title(DoubletFinder))
  
  # Merge classification back to original sample metadata
  sample_metadata <- sample_seurat@meta.data %>%
    tibble::rownames_to_column(var = "Cell") %>%
    dplyr::left_join(doublet_df, by=c("Cell"="Cell")) %>%
    dplyr::mutate(DoubletFinder = dplyr::case_when(is.na(DoubletFinder) ~ "Empty Droplet",
                                                   TRUE ~ DoubletFinder)) %>%
    dplyr::mutate(index = Cell) %>%
    tibble::column_to_rownames(var = "index")
  
  # Assign updated metadata
  sample_seurat@meta.data <- sample_metadata
  
  # Log output
  ident <- as.character(unique(sample_seurat@meta.data$orig.ident))
  message("DoubletFinder doublets identified for sample: '", ident, "'")
  
  return(invisible(sample_seurat))
}

#' Identify doublets using scDblFinder
#'
#' @param sample_seurat Seurat object of a single sample after doublet_finder()
#'
#' @return Seurat object with 'scDblFinder' column ("Singlet"/"Doublet") added to metadata
#' @export
scdbl_finder <- function(sample_seurat){
  
  # Filter out empty droplets before doublet identification
  subset_seurat <- subset(x = sample_seurat,
                          subset = (DropletUtils == "Empty Droplet" & CellRanger == "Empty Droplet"),
                          invert = TRUE)
  
  # Convert to SingleCellExperiment
  sce <- Seurat::as.SingleCellExperiment(x = subset_seurat)
  
  # Run scDblFinder
  scDbl <- scDblFinder::scDblFinder(sce = sce, 
                                    clusters = NULL,
                                    samples = NULL,
                                    dbr = NULL)
  
  # Extract classifications
  dbl_df <- scDbl@colData@listData %>% 
    data.frame() %>% 
    dplyr::rename(scDblFinder = scDblFinder.class) %>% 
    dplyr::mutate(Cell = scDbl@colData@rownames) %>%
    dplyr::mutate(scDblFinder = stringr::str_to_title(scDblFinder)) %>%
    dplyr::select(Cell, scDblFinder)
  
  # Merge into original metadata
  sample_metadata <- sample_seurat@meta.data %>%
    tibble::rownames_to_column(var = "Cell") %>%
    dplyr::left_join(dbl_df, by=c("Cell"="Cell")) %>%
    dplyr::mutate(scDblFinder = dplyr::case_when(is.na(scDblFinder) ~ "Empty Droplet",
                                                 TRUE ~ scDblFinder)) %>%
    dplyr::mutate(index = Cell) %>%
    tibble::column_to_rownames(var = "index")
  
  # Assign updated metadata
  sample_seurat@meta.data <- meta_updated
  
  # Log output
  ident <- as.character(unique(sample_seurat@meta.data$orig.ident))
  message("scDblFinder doublets identified for sample: '", ident, "'")
  
  return(invisible(sample_seurat))
}

#' Calculate cell-level QC metrics for Single Cell/Spatial data
#'
#' @param sample_seurat Seurat object of a single sample after scdbl_finder()
#' @param assay Character. Assay name to compute QC metrics on (e.g. "RNA", "Spatial.008um")
#'
#' @return Seurat object with QC metrics added to metadata
#' @export
calc_qc_metrics_sc_sp <- function(sample_seurat, assay){
  
  # Compute mitochondrial percent
  sample_seurat <- Seurat::PercentageFeatureSet(object = sample_seurat,
                                                pattern = "^[Mm][Tt]-",
                                                col.name = "MitoPercent",
                                                assay = assay)
  
  # Compute ribosomal percent
  sample_seurat <- Seurat::PercentageFeatureSet(object = sample_seurat,
                                                pattern = "^[Rr][Pp][SsLl]", 
                                                col.name = "RiboPercent",
                                                assay = assay)
  
  # Compute hemoglobin percent
  sample_seurat <- Seurat::PercentageFeatureSet(object = sample_seurat,
                                                pattern = "^[Hh][Bb][AaBb]-",                                                 
                                                col.name = "HemePercent",
                                                assay = assay)

  # Rename columns to be more intuitive and add the QC metrics:
  # (i)    Cell      : unique identifiers corresponding to each cell i.e. barcodes
  # (ii)   Sample    : sample names
  # (iii)  nUMIs     : number of transcripts per cell
  # (iv)   nGenes    : number of genes per cell
  # (v)    nHTO_UMIs : number of HTO reads per cell
  # (vi)   nHTOs     : number of HTO types per cell
  # (vii)  MitoRatio : MitoPercent/100
  # (viii) RiboRatio : RiboPercent/100 
  # (ix)   HemeRatio : HemePercent/100
  # (x)    Novelty   : log ratio of genes per UMI
 
  sample_metadata <- sample_seurat@meta.data %>% 
    tibble::rownames_to_column(var = "barcode") %>%
    dplyr::mutate(Cell = paste0(orig.ident, "_", barcode),
                  Sample = orig.ident,
                  nUMIs = get(paste0("nCount_", assay)),
                  nGenes = get(paste0("nFeature_", assay)),
                  MitoRatio = MitoPercent / 100,
                  RiboRatio = RiboPercent / 100,
                  HemeRatio = HemePercent / 100,
                  Novelty = log10(nGenes) / log10(nUMIs))
  
  # Handle HTO metadata if present
  hto_cols <- c("nCount_HTO", "nFeature_HTO", "HTO_Final")
  if (any(hto_cols %in% colnames(sample_metadata))){
    sample_metadata <- sample_metadata %>% 
      dplyr::rename(nHTO_UMIs = nCount_HTO,
                    nHTOs = nFeature_HTO)
  }
  
  # Add spatial coordinates if available
  images <- names(sample_seurat@images)
  if (length(sample_seurat@images) > 0){
    df_coords <- data.frame()
    for (n in 1:length(images)){
      image <-  gsub(pattern=".002|0.008|.016|um", replacement="", x=images[n])
      df <- data.frame("barcode" = sample_seurat@images[[n]]@boundaries$centroids@cells, 
                       X = sample_seurat@images[[n]]@boundaries$centroids@coords[,1], 
                       Y = sample_seurat@images[[n]]@boundaries$centroids@coords[,2])
      
      df_coords <- dplyr::bind_rows(df_coords, df)
    }
    
    sample_metadata <- sample_metadata %>% 
      dplyr::left_join(df_coords, by=c("barcode"="barcode")) 
  }
  
  # Select columns to keep
  keep_cols <- c("Cell", "Sample", "nUMIs", "nGenes", "nHTO_UMIs", "nHTOs", 
                 "HTO_Final", "MitoRatio", "RiboRatio", "HemeRatio", "Novelty",
                 "DropletUtils", "CellRanger", "DoubletFinder", "scDblFinder", 
                 "X", "Y", "barcode")
  keep_cols <- base::intersect(keep_cols, colnames(sample_metadata))
  
  # Restore rownames and subset columns
  sample_metadata <- sample_metadata %>%
    dplyr::select(all_of(keep_cols)) %>%
    tibble::column_to_rownames(var = "barcode")
  
  # Assign updated metadata
  sample_seurat@meta.data <- sample_metadata
  
  # Log output
  ident <- as.character(unique(sample_seurat@meta.data$Sample))
  message("Cell-level QC metrics calculated for sample: '", ident, "'")
  
  return(invisible(sample_seurat))
}

#' Identify low quality cells for Single Cell/Spatial data
#'
#' @param sample_seurat Seurat object of a single sample after QC metric calculation
#'
#' @return Seurat object with 'QC' column ("Singlet"/"Empty Droplet"/"Doublet"/"Low Quality") added to metadata
#' @export
mark_low_quality_sc_sp <- function(sample_seurat){
  
  # Define lenient cutoff thresholds for QC metrics
  if(length(names(sample_seurat@images)) == 0){
    
    # Single-cell cutoffs
    gene_cutoff <- 250
    umi_cutoff <- 500
    mito_cutoff <- 0.2
    ribo_cutoff <- 0.05
    novelty_cutoff <- 0.8
    
    sample_seurat@meta.data <- sample_seurat@meta.data %>% 
      dplyr::mutate(QC = dplyr::case_when((DropletUtils == "Empty Droplet" & CellRanger == "Empty Droplet") ~ "Empty Droplet",
                                          (DoubletFinder == "Doublet" & scDblFinder == "Doublet") ~ "Doublet",
                                          (nGenes >= gene_cutoff & 
                                             nUMIs >= umi_cutoff & 
                                             MitoRatio <= mito_cutoff & 
                                             Novelty >= novelty_cutoff) ~ "Singlet", #RiboRatio >= ribo_cutoff &
                                          TRUE ~ "Low Quality"))
  } else {
    
    # Spatial cutoffs (more lenient)
    gene_cutoff <- 10   # reduced from 250 used for single cell
    umi_cutoff <- 20    # reduced from 500 used for single cell
    mito_cutoff <- 0.2
    ribo_cutoff <- 0.05
    novelty_cutoff <- 0.8
    
    sample_seurat@meta.data <- sample_seurat@meta.data %>% 
      dplyr::mutate(QC = dplyr::case_when((nGenes >= gene_cutoff & 
                                             nUMIs >= umi_cutoff & 
                                             MitoRatio <= mito_cutoff & 
                                             Novelty >= novelty_cutoff) ~ "Singlet", #RiboRatio >= ribo_cutoff &
                                          TRUE ~ "Low Quality"))
  }
  
  # Log output
  ident <- as.character(unique(sample_seurat@meta.data$Sample))
  message("Good quality singlets identified for sample: '", ident, "'")
  
  return(invisible(sample_seurat))
}

#' Append raw metadata from a single sample for plotting purposes
#'
#' @param sample_seurat Seurat object of a single sample
#' @param raw_metadata Data frame containing accumulated raw metadata from previous samples
#'
#' @return Data frame with raw metadata of the current sample appended
#' @export
generate_plotdata <- function(sample_seurat, raw_metadata){
  
  # Append metadata from current sample and filter out any rows with missing 
  # Sample info. This will be used for QC plots later
  raw_metadata <- dplyr::bind_rows(raw_metadata, sample_seurat@meta.data) %>%
    dplyr::filter(!is.na(Sample))
  
  # Log output
  ident <- as.character(unique(sample_seurat@meta.data$Sample))
  message("Raw metadata (for plotting) appended for sample: '", ident, "'")
  
  return(invisible(raw_metadata))
}

#' Filter High-Quality Singlets for Single-cell/Spatial data
#'
#' @param sample_seurat Seurat object of a single sample after mark_low_quality_sc_sp()
#'
#' @return A Seurat object containing only cells labeled as "Singlet" in the QC column.
#' 
#' @export
filter_singlets_sc_sp <- function(sample_seurat){
  
  # Check for required column
  if (!"QC" %in% colnames(sample_seurat@meta.data)) {
    stop("Missing 'QC' column in metadata. Please run `mark_low_quality_sc_sp()` first.")
  }
  
  # Subset to retain only singlets
  sample_seurat <- base::subset(x = sample_seurat, subset = QC == "Singlet")

  # Log output
  ident <- as.character(unique(sample_seurat@meta.data$Sample))
  message("Retained high-quality singlets for sample(s): ", ident)
  
  return(invisible(sample_seurat))
}

#' Merge and Annotate Filtered Seurat Objects (Single-cell or Spatial)
#'
#' This function merges a list of filtered Seurat objects (from single-cell or spatial data),
#' optionally enriches metadata using an external table, and saves the final object to disk.
#'
#' @param samples A character vector of Seurat object names (must exist in the current environment).
#' @param assay A character string indicating the assay type ("RNA", "Spatial.008um", etc.).
#' @param extra_metadata A data frame containing additional metadata to join into the merged Seurat object.
#' @param output_path A character string specifying where to save the resulting `.rds` file.
#'
#' @return A merged and metadata-annotated Seurat object.
#'
#' @export
merge_filtered_sc_sp <- function(samples, assay, extra_metadata, output_path){
  
  # Create a merged seurat object after all the above filtering
  # NOTE: Samples will have the same barcodes. To keep track of cell identities
  # (i.e. barcodes) coming from each sample after merging, we add a prefix
  # (i.e. sample name) to each barcode using "add.cell.ids"
  samples.seurat <- lapply(samples, get)
  
  # Merge seurat objects with unique prefixes
  cell_prefixes <- gsub(pattern=".Spatial.*", replacement="", x=samples)
  merged.seurat <- base::merge(x = samples.seurat[[1]],   #get(paste0(samples[1])
                                 y = samples.seurat[-1],    #lapply(paste0(samples[2:length(samples)]), get)
                                 add.cell.ids = cell_prefixes,
                                 merge.data = FALSE)
  
  # Remove HTO assay if present to avoid complications during integration
  if ("HTO" %in% Assays(merged.seurat)){
    merged.seurat[["HTO"]] <- NULL
  }
  
  # Add extra metadata if any to seurat object
  if (nrow(extra_metadata) > 0){
    meta_data <- merged.seurat@meta.data
    
    if (assay == "RNA"){
      meta_data <- meta_data %>%
        dplyr::mutate(Unique_ID = dplyr::case_when(!is.na(HTO_Final) ~ paste0(Sample, "_", HTO_Final), 
                                                   is.na(HTO_Final) ~ paste0(Sample))) %>%
        dplyr::left_join(extra_metadata, by = ("Unique_ID" = "Unique_ID"))
    } else {
      meta_data <- meta_data %>%
        dplyr::left_join(extra_metadata, by = c("Sample" = "Slide")) %>% 
        dplyr::filter(dplyr::between(X, Xmin, Xmax) & dplyr::between(Y, Ymin, Ymax))
    }
    
    # Add row names before replacing metadata in Seurat object as dplyr::left_join() will remove row names.
    if (!"Cell" %in% colnames(meta_data)) stop("Column 'Cell' containing rownames is not found in metadata after dplyr join.")
    merged.seurat@meta.data <- meta_data %>%
      dplyr::mutate(index = Cell) %>%
      tibble::column_to_rownames(var = "index")  # VERY VERY IMPORTANT
  }
  
  # Save merged seurat object
  if (assay == "RNA"){
    filename <- file.path(output_path, "filtered.seurat.rds")
  } else{  # For Spatial.008um and Spatial.016um assays
    filename <- file.path(output_path, paste0("filtered.seurat.", assay, ".rds"))
  }
  saveRDS(merged_seurat, file = filename)
  
  # Log output
  message("Filtered Seurat object saved to: ", filename)
  
  return(merged.seurat)
}

#' Generate Quality Control (QC) Plots for Single-Cell or Spatial Metadata
#'
#' Generates standard QC plots such as cell counts, UMIs, genes, mitochondrial and ribosomal ratios, novelty scores,
#' and combined scatter plots, and saves them as PDF files.
#'
#' @param raw_metadata A data frame containing cell-level QC metadata with columns:
#'   `Sample`, `QC`, `nUMIs`, `nGenes`, `MitoRatio`, `RiboRatio`, `Novelty`.
#' @param output_path A string path to the directory where QC plots will be saved.
#'
#' @return Generates and saves multiple QC plots to `output_path`. No object is returned.
#'
#' @export
plot_qc <- function(raw_metadata, output_path){
  
  ### Helper function definitions
  qc_levels <- c("Empty Droplet", "Doublet", "Low Quality", "Singlet")
  fill_colors <- c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
                   "Doublet" = "#1F77B4", "Low Quality" = "#D62728")
  
  ### Visualize cell counts per sample
  cell_qc <- function(meta){
    
   df <- meta %>%
      dplyr::count(Sample, QC) %>%
      data.frame() %>%
      dplyr::mutate(QC = factor(QC, levels = qc_levels))
    
    p <- ggplot(data = df, aes(x = Sample, y = n, fill = QC)) + 
      # position = "dodge" for grouped; "stack" for stacked
      # stat = "identity" if y axis defined; "count" if y axis determined based on X axis frequency
      geom_bar(stat = "identity", position = position_dodge(0.9), , drop = FALSE) +             
      theme_classic() +               # display with x and y axis lines and no gridlines
      my_theme +
      labs(x = "Sample", y = "Cell Counts", title = "Number of Cells") +
      coord_cartesian(ylim = c(1,10000000), clip = "off", expand = FALSE) +
      scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000, 1000000)) +
      scale_fill_manual(values = fill_colors) +
      geom_text(aes(label = n, ymin = 0.1, ymax = 1), 
                position = position_dodge(width = 0.9), y = 0.1, hjust = 0, angle = 90)
    #geom_text(stat ="count", aes(label = after_stat(count)), y = 0, hjust = 0, angle = 90)
    
    return(p)
  }
  
  ### Visualize nUMIs, nGenes, MitoRatio, RiboRatio, Novelty per sample
  violin_qc <- function(meta, yvar, ylab, title, cutoff = NULL, ylog = TRUE, ylim = NULL) {
    df <- meta %>%
      dplyr::mutate(QC = factor(QC, levels = qc_levels))
    
    p <- ggplot(data = df, aes(x = Sample, y = .data[[yvar]], fill = QC)) +
      geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) +
      geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) +
      theme_classic() + my_theme +
      labs(x = "Sample", y = ylab, title = title) +
      scale_fill_manual(values = fill_colors)
    
    if (!is.null(cutoff)) p <- p + geom_hline(yintercept = cutoff, linetype = 2)
    if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim, clip = "off")
    if (ylog) p <- p + scale_y_log10()
    
    return(p)
  }
  
  ### Visualize number of genes/cell, number of UMIs/cell & MitoRatio together.
  # Bottom left quadrant : Poor quality cells with low genes & UMIs per cell 
  # Top right quadrant   : Good quality cells with high genes & UMIs per cell
  # Bottom right quadrant: Cells with low genes but high UMIs per cell. These 
  # could be dying cells or population of low complexity cells (i.e erythrocytes)
  gene_umi_mito_qc <- function(meta){
    
    umi_cutoff <- 500
    gene_cutoff <- 250
    df <- meta %>%
      dplyr::mutate(QC = factor(QC, levels = qc_levels))
    
    ggplot(data = df, aes(x = nUMIs, y = nGenes, color = MitoRatio)) +
      geom_point(alpha = 0.5) +
      theme_classic() + 
      my_theme + 
      labs(x = "Number of UMIs", y = "Number of Genes",	 title = "UMIs vs Genes (Colored by MitoRatio)") +
      coord_cartesian(xlim = c(1, 1000000), ylim = c(1, 20000), clip = "off") +
      scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) + 
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000)) + 
      scale_color_viridis(option = "D", limits = c(0, 1)) + 		# limits sets max and min values of gradient
      #facet_wrap(.~Sample, nrow = 4) +   #split the plot by X-axis label
      facet_wrap(Sample ~ QC, ncol = 4) +
      geom_vline(xintercept = umi_cutoff, linetype = "dashed") +    	#draw a vertical line at x=500 i.e.UMIs cutoff
      geom_hline(yintercept = gene_cutoff, linetype = "dashed") +    #draw a horizontal line at y =250 i.e. Genes cutoff
      stat_smooth(method=lm, color="yellow", se = FALSE)
  }
  
  # === Plot generators
  plot_list <- list(
    Cell_Counts                = cell_qc,
    UMI_Distribution           = function(x) violin_qc(x, "nUMIs", "Number of UMIs", "UMI Distribution", 500, TRUE, c(1, 1e6)),
    Gene_Distribution          = function(x) violin_qc(x, "nGenes", "Number of Genes", "Gene Distribution", 250, TRUE, c(1, 30000)),
    MitoRatio_Distribution     = function(x) violin_qc(x, "MitoRatio", "MitoRatio", "MitoRatio Distribution", 0.2, TRUE, c(1e-5, 1)),
    RiboRatio_Distribution     = function(x) violin_qc(x, "RiboRatio", "RiboRatio", "RiboRatio Distribution", 0.05, TRUE, c(1e-4, 1)),
    Novelty_Score_Distribution = function(x) violin_qc(x, "Novelty", "Novelty", "Novelty Score Distribution", 0.8, FALSE, c(0.3, 1)),
    Genes_UMI_MitoRatio_Distribution = gene_umi_mito_qc)
  
  # === Generate and save plots
  for (plot_name in names(plot_list)) {
    p <- plot_list[[plot_name]](raw_metadata)   #  p <- get(funcs[i])(raw_metadata)
    ggplot2::ggsave(filename = paste0("QC_", plot_name, ".pdf"),
                    plot = p,
                    device = "pdf",
                    path = output_path,
                    width = 11,
                    height = 8,
                    dpi = 600,
                    units = "in")
  }
  
  # Log output
  message("QC plots generated and saved to: ", output_path)
}
  
#' Perform SCTransformation and Dimensional Reduction on Single Cell or Spatial Data
#'
#' Performs normalization, cell cycle scoring, SCTransform normalization, 
#' variable feature filtering, and PCA dimensionality reduction on both raw and
#' SCT assays. It handles both single-cell and spatial transcriptomics data.
#'
#' @param filtered.seurat A Seurat object that has been filtered for quality control.
#' @param assay Character. Name of the assay (e.g., "RNA", "Spatial").
#' @param sctransform.var Character. Metadata column to split samples (e.g., "Sample").
#' @param output_path Character. Path for output files (currently unused).
#' @param s_genes Character vector of S phase genes.
#' @param g2m_genes Character vector of G2M phase genes.
#'
#' @return A Seurat object with SCTransform normalization and PCA reductions on both the RNA and SCT assays.
#'
#' @details 
#' - Performs log-normalization before cell cycle scoring.
#' - Joins split layers before cell cycle scoring.
#' - Calculates cell cycle difference score (G2M - S).
#' - Applies SCTransform per sample (split by `sctransform.var`).
#' - Regresses out cell cycle and mitochondrial content.
#' - Filters variable features to exclude ribosomal, mitochondrial, and predicted genes.
#' - Performs PCA on both RNA and SCT assays using their respective variable features.
#' @export
sctransform_sc_sp <- function(filtered.seurat, assay, sctransform.var, output_path){
  
  stopifnot(is(filtered.seurat, "Seurat"))
  if (!assay %in% names(filtered.seurat@assays)) {
    stop("Assay '", assay, "' not found in Seurat object.")
  }
  if (!sctransform.var %in% colnames(filtered.seurat@meta.data)) {
    stop("Metadata column '", sctransform.var, "' not found in Seurat object.")
  }
  
  message("Normalizing data before cell cycle scoring...")
  filtered.seurat <- Seurat::NormalizeData(object = filtered.seurat,
                                           assay = assay,
                                           normalization.method = "LogNormalize",
                                           scale.factor = 10000,
                                           margin = 1,
                                           verbose = TRUE)
  
  message("Joining layers for cell cycle scoring...")
  # NOTE: CellCycleScoring uses a single data layer (i.e. log norm counts) but
  # currently, data layer for each sample is stored separately. Join them first.
  filtered.seurat@assays[[assay]] <- SeuratObject::JoinLayers(filtered.seurat@assays[[assay]])
  
  message("Scoring cell cycle...")
  filtered.seurat <- Seurat::CellCycleScoring(object = filtered.seurat,
                                              s.features = intersect(s_genes,rownames(filtered.seurat@assays[[assay]]@features)),
                                              g2m.features = intersect(g2m_genes, rownames(filtered.seurat@assays[[assay]]@features)),
                                              ctrl = NULL)
  
  # Calculate CC.Score to regress out the difference between G2M & S scores
  # https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette
  filtered.seurat$CC.Score <- filtered.seurat$G2M.Score-filtered.seurat$S.Score
  
  message("Splitting assay by sample before SCTransform...")
  # NOTE: All cells within same batch MUST be analyzed together
  filtered.seurat@assays[[assay]] <- base::split(x = filtered.seurat@assays[[assay]],
                                                 f = filtered.seurat@meta.data[[sctransform.var]])
  
  message("Running SCTransform with CC.Score and MitoRatio regression...") 
  # https://github.com/satijalab/seurat/issues/7342
  sct.seurat <- Seurat::SCTransform(object = filtered.seurat,
                                    assay = assay,
                                    new.assay.name = "SCT",
                                    do.correct.umi = TRUE,
                                    ncells = 5000,
                                    variable.features.n = 3000,
                                    vars.to.regress = c("CC.Score","MitoRatio"),
                                    do.scale = FALSE,
                                    do.center = TRUE,
                                    vst.flavor = "v2",
                                    return.only.var.genes = TRUE,
                                    verbose = TRUE)
  
  message("Filtering out ribosomal, mitochondrial, RIKEN, predicted genes from variable features...")
  # NOTE: Do this to so PCA, UMAP and clustering are not influenced by these genes.
  var_f <- sct.seurat@assays[["SCT"]]@var.features
  var_f <- var_f[!grepl(pattern = "^[Rr][Pp][SsLl]|R[Ii][Kk]$|^[Mm][Tt]-|^G[Mm][0-9.]+$", 
                        x = var_f)]
  sct.seurat@assays[["SCT"]]@var.features <- var_f
  cat("\nFinal number of variable features:", length(var_f), "\n")
  
  message("Scaling and running PCA on original assay...")
  # NOTE: Needed for scVI integration
  sct.seurat <- Seurat::ScaleData(object = sct.seurat, 
                                  assay = assay,
                                  features = VariableFeatures(sct.seurat))
  sct.seurat <- Seurat::RunPCA(object = sct.seurat,
                               assay = assay,
                               features = Seurat::VariableFeatures(sct.seurat),
                               reduction.name = paste0(assay, ".pca"),
                               reduction.key = "PC_")
  
  message("Running PCA on SCT assay...")
  sct.seurat <- Seurat::RunPCA(object = sct.seurat,
                               assay = "SCT",
                               features = Seurat::VariableFeatures(sct.seurat),
                               reduction.name = "sct.pca",
                               reduction.key = "PC_")
  
  # # Create .rds object for sct seurat object
  # if (assay == "RNA"){
  #   saveRDS(sct.seurat, file=paste0(output_path, "sct.seurat.rds"))
  # } else{
  #   # For Spatial.008um and Spatial.016um assays
  #   saveRDS(sct.seurat, file=paste0(output_path, "sct.seurat.", assay, ".rds"))
  # }
  
  message("SCTransform workflow completed successfully.")
  return(invisible(sct.seurat))
}

#' Perform Integration on Single-cell or Spatial Seurat Object
#'
#' Applies multiple integration strategies (CCA, RPCA, Harmony, JointPCA) using SCT-normalized PCA embeddings.
#'
#' @param sct.seurat A Seurat object after SCTransform and PCA.
#' @param assay Character. Assay name (e.g. "RNA" or "Spatial").
#' @param reference.samples Character vector. Sample IDs to use as integration references.
#' @param kweight Numeric. `k.weight` used by RPCA and similar integration methods.
#' @param output_path Character. Optional path to save results (not currently used).
#'
#' @return Invisibly returns the integrated Seurat object.
#' @export
integrate_sc_sp <- function(sct.seurat, assay, reference.samples, kweight, output_path){
  
  stopifnot(is(sct.seurat, "Seurat"))
  if (!assay %in% names(sct.seurat@assays)) {
    stop("Assay '", assay, "' not found in Seurat object.")
  }
  if (!"sct.pca" %in% names(sct.seurat@reductions)) {
    stop("Reduction 'sct.pca' is missing. Please run PCA on SCT assay before integration.")
  }
  
  message("Reference samples: ", paste(reference.samples, collapse = ", "))
  message("k.weight: ", kweight)
 
  
  integrated.seurat <- sct.seurat
  
  # NOTE: While Seurat::CCAIntegration(), Seurat::RPCAIntegration() &
  # Seurat::JointPCAIntegration() all have dims=1:30 as default, only 
  # Seurat::JointPCAIntegration() gives error when integrated.seurat has fewer 
  # than 30 dims. So, we explicitly specify this.
  # Determine the maximum number of dimensions available for integration
  max_dims <- min(30, ncol(integrated.seurat@reductions$sct.pca@cell.embeddings))
  
  integration.methods <- c("CCA", "RPCA", "Harmony", "JointPCA")
  
  for (r in integration.methods){
    
    message("Running ", method, " integration...")
    
    reduction.name <- paste0("integrated.", tolower(method))
    
    integrated.seurat <- Seurat::IntegrateLayers(object = integrated.seurat,
                                                 method = paste0(method, "Integration"),
                                                 normalization.method = "SCT",
                                                 orig.reduction = "sct.pca", 
                                                 new.reduction =  reduction.name,
                                                 reference = reference.samples,
                                                 k.weight = kweight,    # for RPCA
                                                 dims = 1:max_dims,
                                                 verbose = TRUE)
  }
  
  # NOTE: scVI needs raw counts. Vignette also uses it on RNA assay
  # NOTE: We use variable features of SCT assay for integration.
  # NOTE: We use pca reduction from RNA assay (derived using variable features of SCT assay)
  # FastMNN throws error "Error in checkBatchConsistency(batches, cells.in.columns = TRUE)"
  # for (r in c("scVI", "FastMNN")){
  #   DefaultAssay(integrated.seurat) <- assay
  #   integrated.seurat <- Seurat::IntegrateLayers(object = integrated.seurat,
  #                                                method = paste0(r, "Integration"),
  #                                                normalization.method = "LogNormalize",
  #                                                orig.reduction = paste0(assay, ".pca"),
  #                                                features = integrated.seurat@assays$SCT@var.features,
  #                                                new.reduction = paste0("integrated.", base::tolower(r)),
  #                                                reference = ref_samples,
  #                                                k.weight = kweight,                                    # for RPCA
  #                                                conda_env = "/hpc/home/kailasamms/miniconda3/envs/R",  # for scVI
  #                                                verbose = FALSE)
  # }
  
  # Merge layers after integration
  integrated.seurat@assays[[assay]] <- SeuratObject::JoinLayers(integrated.seurat@assays[[assay]])
  
  # Optional RDS saving block (currently commented out)
  # rds_path <- if (assay == "RNA") {
  #   file.path(output_path, "integrated.seurat.rds")
  # } else {
  #   file.path(output_path, paste0("integrated.seurat.", assay, ".rds"))
  # }
  # saveRDS(integrated.seurat, file = rds_path)
  
  message("Integration completed.")
  return(invisible(integrated.seurat))
}

### Perform clustering [Single cell + Spatial]
# Input is seurat object after integration
# Output is seurat object after clustering
cluster_sc_sp <- function(integrated.seurat, assay, output_path){
  
  #***************STEP 8A: FIND NEAREST NEIGHBORS FOR EVERY CELL***************#
  
  # Determine the K-nearest neighbor graph
  for (r in c("CCA", "RPCA", "Harmony", "JointPCA")){
    max_dims <- min(40, ncol(integrated.seurat@reductions[[paste0("integrated.", base::tolower(r))]]@cell.embeddings))
    integrated.seurat <- Seurat::FindNeighbors(object = integrated.seurat,
                                               reduction = paste0("integrated.", base::tolower(r)),
                                               dims = 1:max_dims,
                                               k.param = 30,
                                               graph.name = c(paste0("graph_nn.", base::tolower(r)),
                                                              paste0("graph_snn.", base::tolower(r))))
  }
  
  #**********STEP 8B: SEPARATE CELLS INTO CLUSTERS BASED ON SNN GRAPH**********#
  
  # Determine the clusters for various resolutions
  for (r in c("CCA", "RPCA", "Harmony", "JointPCA")){
    for (res in c(0.4, 0.6, 0.8, 1, 1.2, 1.4)){
      integrated.seurat <- Seurat::FindClusters(object = integrated.seurat,
                                                resolution = res,
                                                graph.name = paste0("graph_snn.", base::tolower(r)),
                                                cluster.name = paste0("cluster.", res, ".", base::tolower(r)),
                                                modularity.fxn = 1,
                                                algorithm = 4)     #4=Leiden is best
    }
  }
  
  #**********STEP 8C: PERFORM DIMENSIONAL REDUCTION FOR VISUALIZATION**********#
  
  # Run UMAP
  for (r in c("CCA", "RPCA", "Harmony", "JointPCA")){ 
    max_dims <- min(40, ncol(integrated.seurat@reductions[[paste0("integrated.", base::tolower(r))]]@cell.embeddings))
    integrated.seurat <- Seurat::RunUMAP(object = integrated.seurat,
                                         dims = 1:max_dims,
                                         n.neighbors = 30L,
                                         reduction = paste0("integrated.", base::tolower(r)),
                                         reduction.name = paste0("umap.", base::tolower(r)))
  }
  
  # # Create .rds object for integrated seurat object
  # if (assay == "RNA"){
  #   saveRDS(integrated.seurat, file=paste0(output_path, "integrated.seurat.rds"))
  # } else{
  #   # For Spatial.008um and Spatial.016um assays
  #   saveRDS(integrated.seurat, file=paste0(output_path, "integrated.seurat.", assay, ".rds"))
  # }
  
  cat("Clustering completed", "\n")
  return(integrated.seurat)
}

### Remove sparse clusters [Single cell + Spatial]
# Input is seurat object after clustering
# Output is seurat object after removing clusters with less than 5 cells
remove_sparse_clusters_sc_sp <- function(integrated.seurat, assay, output_path){
  
  # Get all available resolutions at different reductions
  col_id <- colnames(integrated.seurat@meta.data %>% 
                       dplyr::select(starts_with("cluster.")))
  
  sparse_cells <- c()
  for (id in col_id){
    # Identify clusters which have less than 5 cells
    sparse_clusters <- integrated.seurat@meta.data %>%
      dplyr::count(get(id)) %>%
      dplyr::filter(n <=5) %>%
      dplyr::select(identity(1)) %>%
      unlist(.,use.names=FALSE) %>%
      as.character() %>%
      as.numeric()
    
    print(sparse_clusters)
    
    # Identify the cells in these clusters
    cells <- integrated.seurat@meta.data %>%
      dplyr::filter(get(id) %in% sparse_clusters) %>%
      dplyr::select(Cell) %>%
      unlist(.,use.names=FALSE)
    
    # Create a list of cells identified in sparse clusters at all resolutions 
    # and reductions
    sparse_cells <- c(sparse_cells, cells)
  }
  
  # Remove sparse_cells
  integrated.seurat <- subset(x=integrated.seurat,
                              subset = (Cell %in% unique(sparse_cells)),
                              invert=TRUE)
  
  cat("\nCells removed:", length(unique(sparse_cells)), "\n")
  
  # Create .rds object for integrated seurat object
  if (assay == "RNA"){
    saveRDS(integrated.seurat, file=paste0(output_path, "integrated.seurat.rds"))
  } else{
    # For Spatial.008um and Spatial.016um assays
    saveRDS(integrated.seurat, file=paste0(output_path, "integrated.seurat.", assay, ".rds"))
  }
  
  cat("Integrated seurat object saved after removing sparse clusters (below 5 cells)", "\n")
  return(integrated.seurat) 
}

### Plot metrics post integration [Single cell + Spatial]
# Input is seurat object after clustering & removal of sparse clusters, suffix to be addded to filename, and save path
# Output is (i) series of UMAPs at resolution 0.8 using Harmony reduction
# (ii) series of UMAPs at different resolution using every available reduction 
plot_metrics_post_integration_sc_sp <- function(integrated.seurat, assay, output_path){
  
  #                   Filename                       = c(reduction,             split.by,idents)  
  plot.params <- list(Pre.Integration.PCA.           = c("sct.pca",            "Sample", "cluster.0.8.harmony"),
                      Post.Integration.PCA.          = c("integrated.harmony", "Sample", "cluster.0.8.harmony"),
                      UMAP.Sample.                   = c("umap.harmony",       "Sample", "cluster.0.8.harmony"),
                      UMAP.Phase.                    = c("umap.harmony",       "Phase",  "cluster.0.8.harmony"),
                      UMAP.All.Resolutions.CCA       = c("umap.cca",            NA,      "All"),
                      UMAP.All.Resolutions.RPCA.     = c("umap.rpca",           NA,      "All"),
                      UMAP.All.Resolutions.JointPCA. = c("umap.jointpca",       NA,      "All"),
                      UMAP.All.Resolutions.Harmony.  = c("umap.harmony",        NA,      "All"),
                      UMAP.Singlets.Doublets.        = c("umap.harmony",        NA,      "cluster.0.8.harmony"),
                      UMAP.Numerical.Metrics.        = c("umap.harmony",        NA,      "cluster.0.8.harmony"))
  
  for (i in 1:length(plot.params)){
    
    if((sum(plot.params[[i]][2] %in% colnames(integrated.seurat@meta.data)) > 0) & 
       (sum(plot.params[[i]][1] %in% names(integrated.seurat@reductions)) > 0) &
       (sum(plot.params[[i]][3] %in% colnames(integrated.seurat@meta.data)) > 0)){
      
      plot.seurat <- Seurat::SplitObject(object = integrated.seurat,
                                         split.by = plot.params[[i]][2])
      
      purrr::map(.x = c(1:length(plot.seurat)),
                 .f = function(x){  
                   Idents(plot.seurat[[x]]) <- plot.params[[i]][3]   # define name of cluster
                   Seurat::DimPlot(object = plot.seurat[[x]],
                                   reduction = plot.params[[i]][1],  # define reduction to use
                                   group.by = plot.params[[i]][3],   # define color of cluster
                                   pt.size = 0.1,
                                   order = TRUE,  # plot positive cells above negative cells
                                   label = TRUE,
                                   raster = FALSE,
                                   combine = TRUE) +
                     NoLegend() +
                     my_theme + 
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
    }
    
    else if(grepl(pattern="All.Resolutions", x=names(plot.params)[i]) &
            (sum(plot.params[[i]][1] %in% names(integrated.seurat@reductions)) > 0)){
      
      purrr::map(.x = c(0.4, 0.6, 0.8, 1, 1.2, 1.4),
                 .f = function(x){  
                   idents <- paste0("cluster.", x, gsub(pattern="umap", replacement="", x=plot.params[[i]][1]))
                   Idents(integrated.seurat) <- idents
                   Seurat::DimPlot(object = integrated.seurat,
                                   reduction = plot.params[[i]][1],
                                   group.by = idents,
                                   pt.size = 0.1,
                                   order = TRUE,  # plot positive cells above negative cells
                                   label = TRUE,
                                   raster = FALSE,
                                   combine = TRUE) +
                     NoLegend() +
                     my_theme + 
                     ggplot2::labs(title = base::gsub(pattern="cluster.", replacement="", x=idents))
                 }) %>% cowplot::plot_grid(plotlist=.,
                                           align="hv",
                                           axis="tblr",
                                           nrow=2,
                                           ncol=3,
                                           rel_widths=1,
                                           rel_heights=1,
                                           greedy=TRUE,
                                           byrow=TRUE)
    }
    
    else if ((sum(c("DropletUtils", "CellRanger", "DoubletFinder", "scDblFinder", "QC") %in% colnames(integrated.seurat@meta.data)) == 5) &
             (sum(plot.params[[i]][1] %in% names(integrated.seurat@reductions)) > 0) &
             (sum(plot.params[[i]][3] %in% colnames(integrated.seurat@meta.data)) > 0)){
      
      purrr::map(.x = c("DropletUtils", "CellRanger", "DoubletFinder", "scDblFinder", "QC"),
                 .f = function(x){ 
                   Idents(integrated.seurat) <- plot.params[[i]][3]
                   Seurat::DimPlot(object = integrated.seurat,
                                   reduction =  plot.params[[i]][1],
                                   group.by = x,
                                   pt.size = 0.1,
                                   order = c("Doublet"),  # plot doublets on above rest of cells
                                   label = FALSE,
                                   raster = FALSE,
                                   combine = TRUE) +
                     #NoLegend() +
                     my_theme + 
                     ggplot2::labs(title = x)
                 }) %>% cowplot::plot_grid(plotlist=.,
                                           align="hv",
                                           axis="tblr",
                                           nrow=2,
                                           ncol=3,
                                           rel_widths=1,
                                           rel_heights=1,
                                           greedy=TRUE,
                                           byrow=TRUE)
    }
    
    else if ((sum(c("nUMIs", "nGenes", "S.Score", "G2M.Score", "CC.Score", "MitoRatio") %in% colnames(integrated.seurat@meta.data)) == 6) &
             (sum(plot.params[[i]][1] %in% names(integrated.seurat@reductions)) > 0) &
             (sum(plot.params[[i]][3] %in% colnames(integrated.seurat@meta.data)) > 0)){
      
      purrr::map(.x = c("nUMIs", "nGenes", "S.Score", "G2M.Score", "CC.Score", "MitoRatio"),
                 .f = function(x){ 
                   Idents(integrated.seurat) <- plot.params[[i]][3]
                   Seurat::FeaturePlot(object = integrated.seurat,
                                       reduction = plot.params[[i]][1],
                                       features = x,
                                       min.cutoff='q10',
                                       pt.size = 0.1,
                                       order = TRUE, 
                                       label = FALSE,
                                       raster = FALSE,
                                       combine = TRUE) +
                     #NoLegend() +
                     my_theme + 
                     ggplot2::labs(title = x)
                 }) %>% cowplot::plot_grid(plotlist=.,
                                           align="hv",
                                           axis="tblr",
                                           nrow=2,
                                           ncol=3,
                                           rel_widths=1,
                                           rel_heights=1,
                                           greedy=TRUE,
                                           byrow=TRUE)
    }
    else {
      ggplot() + geom_blank()
      cat("Invalid params:", names(plot.params)[i], ":", plot.params[[i]][1], ":",  plot.params[[i]][2], ":", plot.params[[i]][3], "\n")
    }
    
    # Save the plot
    ggplot2::ggsave(filename = paste0(names(plot.params)[i], assay, ".tiff"),
                    plot = last_plot(),
                    device = "jpeg",
                    path = output_path,
                    scale = 1,
                    width = dplyr::case_when(!is.na(plot.params[[i]][2]) ~ 4*floor(sqrt(length(plot.seurat))),
                                             TRUE ~ 4*3),
                    height = dplyr::case_when(!is.na(plot.params[[i]][2]) ~ 4*ceiling(sqrt(length(plot.seurat))),
                                              TRUE ~ 4*2),
                    units = c("in"),
                    dpi = 600,
                    limitsize = TRUE,
                    bg = "white")
  }
}

### Identify markers for each cluster [Single cell + Spatial]
# Input is seurat object after clustering & removal of sparse clusters
identify_markers_sc_sp <- function(integrated.seurat, assay, resolution, reduction, suffix, output_path){
  
  # Set default assay
  DefaultAssay(integrated.seurat) <- assay
  
  # Change active.ident
  idents <- paste0("cluster.", resolution, ".", base::tolower(reduction))
  Idents(object=integrated.seurat) <- idents
  
  # Find ALL markers
  all_markers <- Seurat::FindAllMarkers(object=integrated.seurat,
                                        assay=assay,
                                        logfc.threshold=0.25,
                                        test.use="wilcox",
                                        slot="data",
                                        min.pct=0.1,
                                        min.diff.pct= -Inf,
                                        only.pos=TRUE)
  
  # Get annotations from ENSEMBL
  annotations_list <- get_annotations()
  if (length(intersect(annotations_list[[1]]$SYMBOL, all_markers$gene)) > 
      length(intersect(annotations_list[[2]]$SYMBOL, all_markers$gene))){
    annotations <- annotations_list[[1]]
  } else {
    annotations <- annotations_list[[2]]
  }
  
  # Add gene descriptions
  all_markers <- all_markers %>% 
    dplyr::mutate(pct.1=dplyr::if_else(pct.1 == 0, 0.001, pct.1),
                  pct.2=dplyr::if_else(pct.2 == 0, 0.001, pct.2),
                  ratio=pct.1/pct.2) %>%
    dplyr::left_join(y=annotations, by=c("gene"="ENSEMBL_SYMBOL")) %>%
    dplyr::relocate(cluster, gene, CHR, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, ratio, DESCRIPTION) %>%
    dplyr::distinct(cluster, gene, avg_log2FC, pct.1, pct.2, .keep_all = TRUE)
  
  # Find top markers for each major cluster
  top_markers <- all_markers %>%
    dplyr::filter(avg_log2FC >= 0.58 & p_val_adj < 0.05) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%  #desc(ratio)
    dplyr::slice_head(n=30) %>%
    ungroup()
  
  # Create a matrix of markers in heatmap format
  mat <- all_markers %>% 
    dplyr::filter(p_val_adj < 0.05) %>%
    tidyr::pivot_wider(id_cols=cluster, names_from = gene, values_from = avg_log2FC, values_fill = 0.0) %>%
    tibble::column_to_rownames("cluster") %>%
    scale()   # column wise scaling so each gene has mean 0, stdev = 1
  
  row_dist <- stats::dist(x=mat, method = "euclidean", diag = TRUE, upper = TRUE)  # distance between rows based on columns
  row_clust <- hclust(row_dist)                                                    # clustering based on distance calculated
  row_order <- rownames(mat[row_clust$order,])
  col_dist <- stats::dist(x=t(mat), method = "euclidean", diag = TRUE, upper = TRUE)
  col_clust <- hclust(col_dist)
  col_order <- colnames(mat[,col_clust$order])
  mat <- mat[row_order, col_order]
  
  mat[mat < 0] <- 0 
  mat <- mat %>% 
    t() %>%
    data.frame() %>% 
    mutate(across(where(is.numeric), ~round(., 2)))
  
  # Save all the markers
  filename <- paste0(proj, ".Markers.All.", idents, ".", suffix,".xlsx")
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb=wb, sheetName="Matrix")
  openxlsx::writeData(wb=wb, sheet="Matrix", x=mat, rowNames = TRUE)
  openxlsx::addWorksheet(wb=wb, sheetName="All_Markers")
  openxlsx::writeData(wb=wb, sheet="All_Markers", x=all_markers)
  openxlsx::addWorksheet(wb=wb, sheetName="Top_Markers")
  openxlsx::writeData(wb=wb, sheet="Top_Markers", x=top_markers)
  openxlsx::saveWorkbook(wb=wb, file=paste0(output_path, filename), overwrite=TRUE)
} 

### Calculate module scores of gene sets
# Input is seurat object after clustering & removal of sparse clusters & 
# xlsx marker file with following format
#   Tcell    Bcell
#   CD4      BANK1
#   CD8A
calc_module_scores <- function(integrated.seurat, marker.file.with.path){
  
  # Read marker file
  marker_df <- read.xlsx(file = marker.file.with.path)
  
  # Set default assay
  DefaultAssay(integrated.seurat) <- "RNA"
  
  # Iterate through each celltype and plot its module scores
  for (i in 1:ncol(marker_df)){
    
    features <- marker_df[,i] %>% unlist(use.names=FALSE)
    features <- rownames(integrated.seurat@assays$RNA$data)[tolower(rownames(integrated.seurat@assays$RNA$data)) %in% 
                                                              tolower(features)]
    features <- list(sort(features))
    
    # Calculate module scores
    if (length(features) > 0){
      integrated.seurat <- Seurat::AddModuleScore(object=integrated.seurat,
                                                  features=features,
                                                  assay="RNA",
                                                  slot="data",
                                                  name=make.names(colnames(marker_df)[i]))
      
      names(features) <- make.names(colnames(marker_df)[i])
      integrated.seurat <- UCell::AddModuleScore_UCell(obj=integrated.seurat,
                                                       features=features,
                                                       assay="RNA",
                                                       slot="data",
                                                       name="_UCell")
    }
  }
  
  return(integrated.seurat)
}

### Annotate based on clusters variable defined by user
annotate_manual_sc_sp <- function(integrated.seurat, clusters, resolution, reduction, suffix, output_path){
  
  idents <- paste0("cluster.", resolution, ".", base::tolower(reduction))
  
  # Make sure you have assigned all clusters to one of the cell types
  # NOTE: "integrated_snn_res.1.4" etc are stored as factors. 
  # So, use as.character() and then as.numeric() to get accurate cluster values
  list_1 <- integrated.seurat@meta.data %>% 
    dplyr::count(get(idents)) %>% 
    dplyr::select(identity(1)) %>% 
    unlist(use.names=FALSE) %>% 
    as.character() %>% 
    as.numeric() %>% 
    sort()
  
  list_2 <- clusters %>% 
    unlist(., use.names=FALSE) %>% 
    sort()
  
  # Proceed with annotation ONLY if all clusters have been annotated
  if (identical(list_1, list_2)){
    print("All Clusters have been annotated")
    
    # Extract metadata from Seurat object, assign appropriate resolution to
    # seurat_clusters column and add Cell.Type, Cell.Subtype columns
    data <- integrated.seurat@meta.data %>% 
      dplyr::mutate(seurat_clusters=get(idents),
                    Cell.Type=NA, Cell.Subtype=NA)
    
    # Assign cell type based on cluster numbers within seurat_clusters column
    for (j in 1:nrow(data)){
      for (i in 1:length(clusters)){
        if (as.numeric(as.character(data$seurat_clusters[j])) %in% clusters[[i]]){
          data$Cell.Type[j] <- names(clusters[i])
        }
      }
    }
    
    # Check summary of cell counts
    print(data %>% dplyr::count(Cell.Type) %>% dplyr::arrange(n))
    cat("\n")
    
    # Import the metadata into Seurat object
    integrated.seurat@meta.data <- data
  } else {
    cat("\nYou missed annotating these clusters:\t", setdiff(list_1, list_2))
    cat("\nThese clusters are not present in data:\t", setdiff(list_2, list_1))
    cat("\nThese clusters have duplicate annotation:\t", list_2[duplicated(list_2)])
  }
  
  saveRDS(integrated.seurat, paste0(seurat_results, "integrated.seurat.", suffix, ".ann.rds"))
  
  return(integrated.seurat)
}  

#******************************************************************************#
#                             MARKER IDENTIFICATION                            #
#******************************************************************************#

find_markers_across_datasets <- function(datasets){
  
  # Rather than trying to identify clusters that are similar across datasets 
  # and then trying to find marker genes, we try to find which genes are 
  # co-expressed frequently across datasets and then assign cell types based on 
  # such co-expressed gene sets. Since our approach is focused on identifying 
  # markers that are frequently co-expressed across all clusters from all 
  # datasets, it is robust to:
  # (i) sequencing depth (low UMIs dataset vs high UMIs dataset)
  # (ii) cell composition (immune enriched vs whole tumor dataset)
  # (iii) experiment type (single cell vs single nuclei dataset)
  # (iv) purity (high ambient RNA vs low ambient RNA dataset)
  # (v) clustering resolution
  
  datasets <- c("scRNASeq_BBN_C57BL6", "scRNASeq_BBN_Rag", "scRNASeq_GSE217093", 
                "scRNASeq_Jinfen", "scRNASeq_Jyoti", "scRNASeq_GSE164557",
                "scRNASeq_Chen", "scRNASeq_GSE222315", "scRNASeq_HRA003620")
  
  # Get human to mouse ortholog mapping
  ortho <- get_orthologs()
  
  # Create empty dataframe to store markers from all datasets
  markers <- data.frame(cluster = c(""))
  
  # Read markers identified at resolution Harmony 0.8 from each dataset
  for (proj in datasets){
    
    seurat_results <- paste0("/hpc/home/kailasamms/scratch/", proj, "/results_seurat/")
    
    # Genes present repeatedly due to DESCRIPTION column. Remove them.
    df <- read.xlsx(paste0(seurat_results, proj, ".Markers.All.cluster.0.8.harmony.Full.xlsx")) %>%
      dplyr::mutate(Proj = gsub(pattern="scRNASeq_", replacement="", x=proj),
                    Proj.cluster = paste0(Proj, ".", cluster),
                    Avg_Expr = avg_log2FC*pct.1) %>%
      dplyr::distinct_at(c("gene", "cluster", "Proj", "p_val_adj", "avg_log2FC", "pct.1", "pct.2"), .keep_all = TRUE) %>%
      dplyr::select(Proj, cluster, Proj.cluster, gene, p_val_adj, avg_log2FC, pct.1, pct.2, ratio, Avg_Expr)
    
    # Merge markers from all datasets
    markers <- dplyr::bind_rows(markers, df) %>%
      dplyr::filter(!is.na(Proj))
    cat(nrow(markers), "\n")
  }
  
  # Add Human orthologs after removing poor markers
  markers <- markers %>%
    dplyr::filter(p_val_adj <= 0.05) %>% #, pct.1 >= 0.4, ratio >= 2) %>%
    dplyr::left_join(ortho, by=c("gene"="Mouse")) %>%
    dplyr::mutate(avg_log2FC = round(avg_log2FC, 2),
                  ratio = round(ratio, 2),
                  Avg_Expr = round( Avg_Expr, 2),
                  p_val_adj = round(p_val_adj, 2),
                  SYMBOL = dplyr::case_when(is.na(Human) ~ gene,
                                            TRUE ~ Human))
  
  # Get top 100 markers based on avg_log2FC, ratio, Avg_Expr, pct.1 for each cluster
  markers_log2FC <- markers %>%
    dplyr::group_by(Proj.cluster) %>%
    dplyr::slice_max(n=100, order_by = avg_log2FC)
  
  markers_ratio <- markers %>%
    dplyr::group_by(Proj.cluster) %>%
    dplyr::slice_max(n=100, order_by = ratio) %>%
    dplyr::filter(ratio > 1)
  
  markers_pct1 <- markers %>%
    dplyr::group_by(Proj.cluster) %>%
    dplyr::slice_max(n=100, order_by = pct.1) %>%
    # if you exclude this filter, you will get specific but sparsely expressed genes
    dplyr::filter(pct.1 >= 0.4)  
  
  markers_expr <- markers %>%
    dplyr::group_by(Proj.cluster) %>%
    dplyr::slice_max(n=100, order_by = Avg_Expr)
  
  # Merge top 100 markers for each cluster and remove duplicates
  markers_top <- dplyr::bind_rows(markers_log2FC, markers_ratio, 
                                  markers_pct1, markers_expr) %>%
    dplyr::distinct_all(.keep_all = TRUE)
  
  # Get all possible combinations of markers
  marker.pair.df <- data.frame()
  count <- 0
  for (i in unique(markers_top$Proj.cluster)){
    
    # Get markers from each proj and cluster
    markers.subset <- markers_top %>%
      dplyr::filter(Proj.cluster == i)
    
    count <- count+1
    cat(count, ":", nrow(markers.subset), "\t")
    
    if (nrow(markers.subset) >= 2){
      df <- utils::combn(x=markers.subset$SYMBOL, m=2) %>%
        #df <- mixtools::perm(n=length(markers.subset$SYMBOL), r=2, v=markers.subset$SYMBOL)
        t() %>%
        data.frame()
      
      marker.pair.df <- dplyr::bind_rows(marker.pair.df, df)
    }
  }  
  
  # Count all possible combinations of markers from all datasets
  marker.pair.df <- marker.pair.df %>%
    dplyr::rename(PairA = identity(1), PairB = identity(2)) %>%
    dplyr::add_count(PairA, PairB) %>%
    dplyr::distinct_at(c("PairA", "PairB"), .keep_all = TRUE) %>%
    dplyr::rename(n_clusters = n) %>%
    dplyr::filter(n_clusters >=3)   # remove combinations not observed in even 3 clusters
  
  # Calculate overlapping number of genes between PairA and PairB
  marker.pair.df$n_common <- 0  # number of genes commonly coexpressed between A & B
  marker.pair.df$n_PairA <- 0   # number of genes coexpressed with A
  marker.pair.df$n_PairB <- 0   # number of genes coexpressed with B
  marker.pair.df$n_ratio <- 0   # n_common/(n_PairA+n_PairB-n_Common)
  for (i in 1:nrow(marker.pair.df)){
    
    # Find all genes coexpressed with PairA
    coexp_A1 <- marker.pair.df %>% 
      dplyr::filter(PairA == marker.pair.df$PairA[i]) %>% 
      dplyr::select(PairB) %>% unlist(use.names=FALSE)
    coexp_A2 <- marker.pair.df %>% 
      dplyr::filter(PairB ==  marker.pair.df$PairA[i]) %>% 
      dplyr::select(PairA) %>% unlist(use.names=FALSE)
    coexp_A <- unique(coexp_A1, coexp_A2)
    
    # Find all genes coexpressed with PairB
    coexp_B1 <- marker.pair.df %>% 
      dplyr::filter(PairA ==  marker.pair.df$PairB[i]) %>% 
      dplyr::select(PairB) %>% unlist(use.names=FALSE)
    coexp_B2 <- marker.pair.df %>% 
      dplyr::filter(PairB == marker.pair.df$PairB[i]) %>% 
      dplyr::select(PairA) %>% unlist(use.names=FALSE)
    coexp_B <- unique(coexp_B1, coexp_B2)
    
    # Calculate stats
    marker.pair.df$n_common[i] <- length(intersect(coexp_A, coexp_B))
    marker.pair.df$n_PairA[i] <- length(coexp_A)
    marker.pair.df$n_PairB[i] <- length(coexp_B)
    marker.pair.df$n_ratio[i] <- marker.pair.df$n_common[i]/(marker.pair.df$n_PairA[i]+marker.pair.df$n_PairB[i]-marker.pair.df$n_common[i])
    
    cat(i, "\t")
  }
  
  # Remove all pairs that have poor overlap (n_ratio < 0.5) 
  top.marker.pair.df <- marker.pair.df %>% 
    dplyr::filter(n_ratio >= 0.5)
  
  final.markers <- list(T.NK.cell        = c("CD3D"),
                        B.Plasma.cell    = c("CD79A"),
                        Erythrocyte      = c("HBB"),
                        Mast.cell        = c("KIT"),   # tissue resident granule producing cell
                        #Granulocyte      = c(),        # blood resident granule producing cell (Basophil, Eosinophil, Neutrophil)                     
                        Monocyte         = c("GOS2"),  # blood resident phagocyte
                        Macrophage       = c("C1QA"),  # tissue resident phagocyte
                        #Dendritic.cell   = c(),
                        Endothelial.cell = c("VWF"),
                        Myocyte          = c("MYH11"),
                        Neurons          = c("KCNA1"))
  
  # Fill the marker list
  for (i in 1:length(final.markers)){
    
    final.markers[[i]] <- c(final.markers[[i]], marker.pair.df %>% 
                              dplyr::filter(PairA == final.markers[[i]]) %>% 
                              dplyr::select(PairB) %>%
                              unlist(use.names=FALSE))
  }
  
  # Convert to dataframe
  max_l <- max(lengths(final.markers)) 
  final.markers.df <- lapply(X=final.markers, FUN=function(x){c(x, base::rep(x="", times=max_l-length(x)))})
  
  # Save the clustered similarity matrix
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Markers")
  openxlsx::writeData(wb, sheet = "Markers", x = final.markers.df, rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName = "Datasets")
  openxlsx::writeData(wb, sheet = "Datasets", x = data.frame(datasets), rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName = "All")
  openxlsx::writeData(wb, sheet = "All", x = marker.pair.df, rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName = "Top")
  openxlsx::writeData(wb, sheet = "Top", x = top.marker.pair.df, rowNames = FALSE)
  openxlsx::saveWorkbook(wb = wb, file = "Markers.Compiled.xlsx", overwrite = TRUE)
}

plot_markers_across_datasets <- function(datasets, markers){
  
  # Read integrated seurat object for each dataset
  for (d in datasets){
    seurat_results <- paste0("/hpc/home/kailasamms/scratch/", d, "/results_seurat/")
    integrated.seurat <- readRDS(paste0(seurat_results, "integrated.seurat.rds"))
    Idents(integrated.seurat) <- "cluster.0.8.harmony"
    assign(d, integrated.seurat)
  }
  
  #*****************Plot UMAP at Harmony 0.8 for each dataset******************#
  
  purrr::map(.x = datasets,
             .f = function(x){
               obj <- get(x)
               Seurat::DimPlot(object=obj,
                               reduction="umap.harmony",
                               cols=my_palette,
                               pt.size=0.2,
                               label.size=1,
                               order = TRUE,  # plot doublets on above rest of cells
                               label = TRUE,
                               raster = FALSE,
                               combine = TRUE) +
                 NoLegend() +
                 ggplot2::labs(title = x,  x="UMAP_1", y="UMAP_2") +
                 my_theme}) %>%
    cowplot::plot_grid(plotlist=.,
                       align="hv",
                       axis="tblr",
                       #nrow=,
                       ncol=3,
                       rel_widths=1,
                       rel_heights=1,
                       greedy=TRUE,
                       byrow=TRUE)
  
  # Save the plot
  ggplot2::ggsave(filename = "Marker_UMAPs.tiff",
                  plot = last_plot(),
                  device = "jpeg",
                  #path = ,
                  scale = 1,
                  width = 4*3,
                  height = 4*ceiling(length(datasets)/3),
                  units = c("in"),
                  dpi = 300,
                  limitsize = TRUE,
                  bg = "white")
  
  #*****************Plot UMAP of identified markers for each dataset******************#
  
  # Plot these putative markers on UMAP
  for (f in markers){
    
    purrr::map(.x = datasets,
               .f = function(x){
                 obj <- get(x)
                 Seurat::FeaturePlot(object = obj,
                                     features = intersect(rownames(obj@assays$RNA$counts), c(f, stringr::str_to_title(f))),
                                     reduction = "umap.harmony",
                                     cols = c("grey", viridis(n = 10, option = "C", direction = -1)),
                                     pt.size = 0.2,
                                     label.size = 1,
                                     min.cutoff='q10',
                                     order = TRUE,  # plot doublets on above rest of cells
                                     label = TRUE,
                                     raster = FALSE,
                                     combine = TRUE) +
                   NoLegend() +
                   # scale_colour_gradientn(colours=rev(brewer.pal(n=11, name="RdBu"))[5:11])
                   ggplot2::labs(title = x, x="UMAP_1", y="UMAP_2") +
                   my_theme}) %>% cowplot::plot_grid(plotlist=.,
                                                     align="hv",
                                                     axis="tblr",
                                                     #nrow=,
                                                     ncol=3,
                                                     rel_widths=1,
                                                     rel_heights=1,
                                                     greedy=TRUE,
                                                     byrow=TRUE)
    
    ggplot2::ggsave(filename = paste0(f, ".tiff"),
                    plot = last_plot(),
                    device = "jpeg",
                    #path = ,
                    scale = 1,
                    width = 4*3,
                    height = 4*ceiling(length(datasets)/3),
                    units = c("in"),
                    dpi = 300,
                    limitsize = TRUE,
                    bg = "white")
    
  }
}

#******************************************************************************#
#                               GENE ANNOTATIONS                               #
#******************************************************************************#

### Get human-mouse orthologs
# Output is a dataframe with columns DB.Class.Key, Human, Mouse
# NOTE: Mouse genes (H2-Q10,H2-Q8,H2-Q7,H2-Q6,H2-Q4,H2-Q2,H2-Q1,H2-T23,H2-K1,
# H2-D1) are othologs of the same human gene HLA-A. Similarly, human genes 
# ZNG1A, ZNG1B, ZNG1C, ZNG1E, ZNG1F) are orthologs of the same mouse gene (Zng1).
# So, we make them unique as well syntactically valid (hyphens repalced with .)
# to avoid errors in data analysis.
get_orthologs <- function(){
  
  # This website has a list of human mouse orthologs
  df <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  # Get human genes
  df_h <- df %>% 
    dplyr::select(DB.Class.Key, Common.Organism.Name, Symbol) %>%
    dplyr::filter(Common.Organism.Name == "human") %>%
    dplyr::select(DB.Class.Key, Symbol) %>%
    dplyr::distinct_at(c("DB.Class.Key","Symbol"),.keep_all = TRUE) %>%
    dplyr::rename(Human = Symbol)
  
  # Get mouse genes 
  df_m <- df %>% 
    dplyr::select(DB.Class.Key, Common.Organism.Name, Symbol) %>%
    dplyr::filter(Common.Organism.Name == "mouse, laboratory") %>%
    dplyr::select(DB.Class.Key, Symbol) %>%
    dplyr::distinct_at(c("DB.Class.Key","Symbol"),.keep_all = TRUE) %>%
    dplyr::rename(Mouse = Symbol)
  
  # Get human-mouse orthologs & remove genes that dont have orthologs
  df_h_m <- dplyr::full_join(df_h, df_m, by=c("DB.Class.Key"="DB.Class.Key")) %>%
    base::replace(is.na(.), "None") %>%
    dplyr::filter(Human != "None", Mouse != "None")
  #df_h_m[is.na(df_h_m)] <- "None"
  
  # Similar orthologs (mouse and human gene names are identical)
  conf_h_m <- df_h_m %>% dplyr::filter(Human == base::toupper(Mouse))
  
  # Dissimilar orthologs (mouse and human gene names are NOT identical)
  fix_h_m <- df_h_m %>% 
    dplyr::filter(!(Human %in% conf_h_m$Human)) %>%
    dplyr::filter(!(Mouse %in% conf_h_m$Mouse))
  
  # Merge
  df <- dplyr::bind_rows(conf_h_m, fix_h_m) %>%
    dplyr::mutate(Human = make.names(Human, unique=TRUE),
                  Mouse = make.names(Mouse, unique=TRUE))
  
  # Save the excel file  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb=wb, sheetName="Homologs")
  openxlsx::writeData(wb=wb, sheet="Homologs", x=df)
  openxlsx::saveWorkbook(wb=wb, file="Human.Mouse.Homologs.xlsx", overwrite=TRUE)
  
  return (df)
  
  # # Get mouse and human annotations
  # df <- get_annotations()
  # h <- df[[1]]
  # m <- df[[2]]
  # 
  # # Identify identical genes
  # common <- base::intersect(h$SYMBOL, base::toupper(m$SYMBOL))
  # common_h <- h %>% dplyr::filter(h$SYMBOL %in% common)
  # common_m <- m %>% dplyr::filter(base::toupper(m$SYMBOL) %in% common)
  # 
  # # Identify non-common genes
  # uniq_h <- h %>% dplyr::filter(!(h$SYMBOL %in% common))
  # uniq_m <- m %>% dplyr::filter(!(base::toupper(m$SYMBOL) %in% common))
  # 
  # # Identify similar genes (TIME CONSUMING. SO, save excel and use it in future)
  # # Eg: CD3 and CD3D etc
  # similar_h <- c()
  # similar_m <- c()
  # for (p in unique(uniq_h$SYMBOL)){
  #   
  #   if (sum(grepl(pattern=paste0(p, "."), x=toupper(unique(uniq_m$SYMBOL)))) > 0){
  #     similar_h <- c(similar_h, p)
  #     similar_m <- c(similar_m, unique(uniq_m$SYMBOL)[grepl(pattern=paste0(p, "."), x=toupper(unique(uniq_m$SYMBOL)))][1])
  #                  
  #     #cat(p, ":", unique(uniq_m$SYMBOL)[grepl(pattern=paste0(p, "."), x=toupper(unique(uniq_m$SYMBOL)))][1], "\n")
  #   }
  # }
  # 
  # # Generate final dataframe of human-mouse homolog
  # common_df <- data.frame(HUMAN = unique(sort(common_h$SYMBOL)), 
  #                  MOUSE = unique(sort(common_m$SYMBOL)))
  # similar_df <- data.frame(HUMAN = similar_h,
  #                          MOUSE = similar_m)
  # df1 <- dplyr::bind_rows(common_df, similar_df)
}

#******************************************************************************#
#                      SPATIAL ANALYSIS RELATED FUNCTIONS                      #       
#******************************************************************************#

# We need to run SCT on each assay. But for cell cycle scoring, NA is present
# in barcodes from other assays like 002um and 016um. So, it throws error.
# Therefore, create a separate seurat object for each bin size so that each
# seurat object has only one assay equivalent to "RNA" assay of single cell data

### Perform clustering
# Input is seurat object after SCTransformation. There is NO integration step.
# Output is seurat object after clustering
cluster_spatial <- function(integrated.seurat, output_path){
  
  #***************STEP 8A: FIND NEAREST NEIGHBORS FOR EVERY CELL***************#
  
  # Determine the K-nearest neighbor graph
  integrated.seurat <- Seurat::FindNeighbors(object = integrated.seurat,
                                             reduction = paste0(assay, ".sct.pca"),
                                             dims = 1:40,
                                             k.param = 30,
                                             graph.name = c(paste0(assay, ".graph_nn"),
                                                            paste0(assay, ".graph_snn")))
  
  #**********STEP 8B: SEPARATE CELLS INTO CLUSTERS BASED ON SNN GRAPH**********#
  
  # Determine the clusters for various resolutions
  for (res in c(0.4, 0.6, 0.8, 1, 1.2, 1.4)){
    integrated.seurat <- Seurat::FindClusters(object = integrated.seurat,
                                              resolution = res,
                                              graph.name = paste0(assay, ".graph_snn"),
                                              cluster.name = paste0(assay, ".cluster.", res),
                                              modularity.fxn = 1,
                                              algorithm = 4)     #4=Leiden is best
  }
  
  #**********STEP 8C: PERFORM DIMENSIONAL REDUCTION FOR VISUALIZATION**********#
  
  # Run UMAP
  integrated.seurat <- Seurat::RunUMAP(object = integrated.seurat,
                                       dims = 1:40,
                                       n.neighbors = 30L,
                                       reduction = paste0(assay, ".sct.pca"),
                                       reduction.name = paste0(assay, ".umap"))
  
  # # Create .rds object for integrated seurat object
  # saveRDS(integrated.seurat, file=paste0(output_path, "integrated.seurat.rds"))
  
  cat("Clustering completed", "\n")
  return(integrated.seurat)
}

### Plot location of samples on spatial map
# Input is seurat object of a single slide with columns X, Y, Sample, Group
plot_spatial_map <- function(plot.seurat, x1, y1, x2, y2, suffix, output_path){
  
  Sample <- unique(plot.seurat@meta.data$Sample)
  
  # Get an idea of co-ordinates
  xmin <- (min(plot.seurat@meta.data$X)%/%100-1)*100
  xmax <- (max(plot.seurat@meta.data$X)%/%100+1)*100
  ymin <- (min(plot.seurat@meta.data$Y)%/%100-1)*100
  ymax <- (max(plot.seurat@meta.data$Y)%/%100+1)*100
  
  if (Sample == "TMA1-A1"){
    # Flip on X axis to get correct orientation for TMA1-A1
    p1 <- ggplot(plot.seurat@meta.data, mapping=aes(x=X, y=Y, color = Sample)) +
      geom_point() +
      scale_x_reverse(limits =c(xmax, xmin), breaks=seq(xmax, xmin, -100), position="top") +
      scale_y_continuous(limits =c(ymin, ymax), breaks=seq(ymin, ymax, 100), position="left") +
      geom_vline(xintercept = x1) +
      geom_hline(yintercept = y1) +
      scale_color_manual(values=my_palette)
    
    p2 <- ggplot(plot.seurat@meta.data, mapping=aes(x=X, y=Y, color = Group)) +
      geom_point() +
      scale_x_reverse(limits =c(xmax, xmin), breaks=seq(xmax, xmin, -100), position="top") +
      scale_y_continuous(limits =c(ymin, ymax), breaks=seq(ymin, ymax, 100), position="left") +
      geom_vline(xintercept = x1) +
      geom_hline(yintercept = y1) +
      scale_color_manual(values=my_palette)
  } else if (Sample == "TMA1-D1"){
    # Flip on Y axis to get correct orientation for TMA1-D1
    p1 <- ggplot(plot.seurat@meta.data, mapping=aes(x=X, y=Y, color = Sample)) +
      geom_point() +
      scale_x_continuous(limits =c(xmin, xmax), breaks=seq(xmin, xmax, 100), position="top") +
      scale_y_reverse(limits =c(ymax, ymin), breaks=seq(ymax, ymin, -100), position="left") +
      geom_vline(xintercept = x2) +
      geom_hline(yintercept = y2) +
      scale_color_manual(values=my_palette)
    
    p2 <- ggplot(plot.seurat@meta.data, mapping=aes(x=X, y=Y, color = Group)) +
      geom_point() +
      scale_x_continuous(limits =c(xmin, xmax), breaks=seq(xmin, xmax, 100), position="top") +
      scale_y_reverse(limits =c(ymax, ymin), breaks=seq(ymax, ymin, -100), position="left") +
      geom_vline(xintercept = x2) +
      geom_hline(yintercept = y2) +
      scale_color_manual(values=my_palette)
  }
  
  # Save the plot
  ggplot2::ggsave(filename=paste0(Sample, ".", suffix, ".Sample.Map.tiff"),
                  plot=p1+p2,
                  device="jpeg",
                  path=output_path,
                  width=25,
                  height=8.5,
                  units=c("in"),
                  dpi=600,
                  limitsize=TRUE,
                  bg="white")
}  

#******************************************************************************#
#                             DEPRECATED FUNCTIONS                             #
#******************************************************************************#

# DEPRECATED (used during Seurat v3)
v3_sctransform_singlecell <- function(filtered.seurat){
  
  # Seurat v5 stores counts of each sample in separate layers. Merge them.
  filtered.seurat@assays$RNA <- SeuratObject::JoinLayers(filtered.seurat@assays$RNA)
  
  # Split each sample into a seurat object to get a list of seurat object
  split.seurat <- Seurat::SplitObject(object = filtered.seurat,
                                      split.by = "Sample")
  
  # Remove samples with less than 50 cells so that RunPCA() doesnt give error
  split.seurat <- split.seurat[names(split.seurat)[sapply(split.seurat, ncol) > 50]]
  
  for (i in 1:length(split.seurat)){
    
    # Normalize the data
    split.seurat[[i]] <- Seurat::NormalizeData(object = split.seurat[[i]],
                                               assay = "RNA",
                                               normalization.method = "LogNormalize",
                                               scale.factor = 10000,
                                               margin = 1,
                                               verbose = TRUE)
    
    # Perform cell cycle scoring
    split.seurat[[i]]  <- Seurat::CellCycleScoring(object = split.seurat[[i]],
                                                   s.features = intersect(s_genes,rownames(split.seurat[[i]]@assays$RNA@features)),
                                                   g2m.features = intersect(g2m_genes, rownames(split.seurat[[i]]@assays$RNA@features)),
                                                   ctrl = NULL)
    
    split.seurat[[i]]$CC.Score <- split.seurat[[i]]$G2M.Score-split.seurat[[i]]$S.Score
    
    # SCTransform() is better than FindVariableFeatures() & ScaleData()
    # split.seurat[[i]] <- Seurat::FindVariableFeatures(object = split.seurat[[i]],
    #                                                   assay = "RNA",
    #                                                   selection.method = "vst",
    #                                                   nfeatures = 2000)
    # split.seurat[[i]] <- Seurat::ScaleData(object = split.seurat[[i]],
    #                                        features = NULL,
    #                                        assay = "RNA",
    #                                        vars.to.regress = NULL)
    
    # Perform scaling & variable feature identification usign SCTransform()
    split.seurat[[i]] <- Seurat::SCTransform(object =  split.seurat[[i]],
                                             assay = "RNA",
                                             new.assay.name = "SCT",
                                             do.correct.umi = TRUE,
                                             ncells = 5000,
                                             variable.features.n = 3000,
                                             vars.to.regress = c("CC.Score","MitoRatio"),
                                             do.scale = FALSE,
                                             do.center = TRUE,
                                             vst.flavor = "v2",
                                             return.only.var.genes = TRUE)
    
    
    # Remove ribosomal, Riken, predicted and mitochondrial genes from
    # VariableFeatures so that PCA, UMAP and hence clustering are not affected
    var_f <- split.seurat[[i]]@assays$SCT@var.features
    var_f <- var_f[!grepl(pattern = "^[Rr][Pp][SsLl]|R[Ii][Kk]$|^[Mm][Tt]-|^G[Mm][0-9.]+$", 
                          x=var_f)]
    
    split.seurat[[i]]@assays$SCT@var.features <- var_f
    cat("\nFinal number of variable features:", length(var_f), "\n")
    
    # Perform dimensional reduction using PCA on SCT assay variable features
    split.seurat[[i]] <- Seurat::RunPCA(object = split.seurat[[i]],
                                        assay = "SCT",
                                        features = NULL,
                                        ndims.print = 1,
                                        nfeatures.print = 1,
                                        reduction.name = "pca",
                                        reduction.key = "PC_")
    
    # Perform dimensional reduction using UMAP on PCA dimensions
    split.seurat[[i]] <- Seurat::RunUMAP(object = split.seurat[[i]],
                                         dims = 1:40,
                                         reduction = "pca",
                                         reduction.name = "umap",
                                         reduction.key = "UMAP_")
    
  }
  
  return(split.seurat)
}

# DEPRECATED (used during Seurat v3)
v3_integrate_singlecell <- function(sct, ref_samples){
  
  # NOTE: In v3, Harmony integration was not possible as FindIntegrationAnchors()
  # doesnt support reduction="harmony". Moreover, integration using rpca, cca, 
  # jpca couldnt be stored in same object since FindIntegrationAnchors() output
  # varies for each method
  
  #***STEP 7B: SELECT 3000 MOST VARIABLE GENES TO USE FOR INTEGRATING THE DATA***#
  integ_features <- Seurat::SelectIntegrationFeatures(object.list=split.seurat,
                                                      nfeatures=3000,
                                                      assay=NULL) #c("SCT", "SCT"),
  
  
  #******************STEP 7C: FIND RESIDUALS FOR MISSING GENES*******************#
  split.seurat <- Seurat::PrepSCTIntegration(object.list=split.seurat,
                                             assay="SCT",
                                             anchor.features=integ_features)
  
  #******STEP 7D: FIND COMMON ANCHORS BETWEEN SAMPLES TO INTEGRATE THE DATA******#
  integ_anchors <- Seurat::FindIntegrationAnchors(object.list=split.seurat,
                                                  reference=ref_samples,
                                                  anchor.features=integ_features,
                                                  scale=TRUE,
                                                  normalization.method="SCT",
                                                  sct.clip.range=NULL,
                                                  reduction="rpca", #"cca", "jpca", "rlsi"
                                                  l2.norm=TRUE,
                                                  dims=1:30)
  
  #******STEP 7E: FIND OPTIMUM k.weight FOR USE IN Seurat::IntegrateData()*****#
  # Find minimum anchors between 2 datasets
  kweight1 <- as.data.frame(integ_anchors@anchors) %>%
    dplyr::group_by(dataset1, dataset2) %>%
    distinct_at("cell1", .keep_all=TRUE) %>%
    dplyr::summarize(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(n) %>%
    unlist(use.names=FALSE) %>%
    min()
  
  # Find half of number of cells in sample with least cell count
  kweight2 <- filtered.seurat@meta.data %>%
    dplyr::count(Sample) %>%
    dplyr::filter(n >=50) %>%
    dplyr::select(n) %>%
    unlist(use.names=FALSE) %>%
    min()
  
  kweight2 <- floor(kweight2/2)
  
  kweight <- base::min(kweight1, kweight2, 100)
  dplyr::if_else(kweight >= 100, 100, kweight)
  cat("\n", celltype, "\tkweight1:", kweight1, "\tkweight2:", kweight2, "\tkweight:",  kweight, "\n")
  
  # NOTE: Integration will not fail anymore. If it fails, identify the 2
  # datasets that are involved in the error and use kweight=number of anchors
  # for these 2 datasets.
  cat("\nNumber of unique anchors between datasets\n")
  print(as.data.frame(integ_anchors@anchors) %>%
          dplyr::group_by(dataset1, dataset2) %>%
          distinct_at("cell1", .keep_all=TRUE) %>%
          dplyr::summarize(n=n()), n=1000)
  
  #************************STEP 7F: INTEGRATE THE DATA*************************#
  # NOTE: weight.reduction=NULL means new PCA will be calculated & used to
  # calculate anchor weights
  integrated.seurat.rpca <- Seurat::IntegrateData(anchorset=integ_anchors.rpca,
                                                  new.assay.name="integrated",
                                                  normalization.method="SCT",
                                                  features=NULL,
                                                  features.to.integrate=NULL,
                                                  dims=1:30,
                                                  k.weight=kweight, #default is 100
                                                  weight.reduction=NULL,
                                                  sd.weight=1)
  
  #**STEP 7G: RUN PCA USING 3000 INTEGRATION FEATURES & UMAP USING FIRST 40 PCs**#
  integrated.seurat <- Seurat::RunPCA(object=integrated.seurat,
                                      assay="integrated",
                                      features=NULL)
  
  integrated.seurat <- Seurat::RunUMAP(object=integrated.seurat,
                                       dims=1:40,
                                       reduction="pca")
  return(integrated.seurat)
}

# DEPRECATED (used during Seurat v3)
### Generate whitelist for CITESeq
# Input is filtered seurat object
# Output is a list of csv files - one per batch containing valid barcodes
v3_generate_whitelist <- function(filtered.seurat, output_path){
  
  # Extract barcodes and split by "_"
  bc <- filtered.seurat@meta.data$Cell
  
  # Adjust this based on how your samples are named
  # NOTE: There will be multiple samples within each batch
  barcodes <- data.frame(stringr::str_split_fixed(bc, "_", 2)) %>%
    dplyr::rename(Batch = identity(1), Barcodes = identity(2)) %>%
    dplyr::mutate(Barcodes = stringr::str_replace(Barcodes, "-1", ""),
                  Batch = gsub(pattern="-GEX.*", replacement="", x=Batch))
  
  # Remove duplicate barcodes within each batch
  barcodes <- barcodes %>% 
    dplyr::group_by(Batch) %>% 
    dplyr::distinct_at("Barcodes", .keep_all=TRUE) %>% 
    as.data.frame()
  
  # Check how many barcodes are present in each batch
  barcodes %>% dplyr::group_by(Batch) %>% dplyr::count()
  
  # Save barcodes from each batch to individual csv files
  for (i in unique(barcodes$Batch)){
    whitelist <- barcodes %>%
      dplyr::filter(Batch == i) %>%
      dplyr::select(Barcodes)
    
    write.table(x = whitelist,
                file = paste0(scripts_path, proj, "_", i, "_whitelist.csv"),
                row.names = FALSE,
                col.names = FALSE)
  }
}

# DEPRECATED (used during Seurat v3)
# Input is path to folder containing h5ad
# Output is a raw seurat object
v3_read_h5ad <- function(input_path){
  
  # Load h5ad (useful if analyzing collaborator data in h5ad format)
  SeuratDisk::Convert(source = paste0(input_path, proj, ".h5ad"),
                      dest = "h5seurat",
                      assay="RNA",
                      overwrite = FALSE)
  
  raw_seurat <- SeuratDisk::LoadH5Seurat(file = paste0(input_path, proj, ".h5seurat"))
  
  return(raw_seurat)
}

### Import data from output of citeseq
# Input is path to demux_results folder
# Ouput is seurat object of each sample
v3_read_citeseq <- function(input_path){
  
  # Create a list of samples that have been demultiplexed already
  files <- list.files(path = paste0(input_path, "singlets/"),
                      full.names = FALSE)
  samples <- gsub(pattern="\\..*", replacement="", x=files)
  
  # Loop through each of the individual object in demux directory & import data
  for (i in 1:length(files)){
    
    # Read the seurat object containing demultiplexed singlets
    sample.seurat <- readRDS(file = paste0(demux_results, "singlets/", files[i]))
    
    # Assign the seurat object to its corresponding variable
    assign(samples[i], sample.seurat)
  }
}

SpatialFeaturePlotBlend <- function(cells_obj, column_1, column_2, combine=TRUE){
  
  # Convert decimal number to hexadecimal. Pad with 0s if only a single
  # character following conversion.
  as_hex <- function(num) {
    hex_str <- as.character(as.hexmode(num))
    if (nchar(hex_str) == 1) {
      hex_str <- paste0("0", hex_str)
    }
    
    return(hex_str)
  }
  
  metadata_to_hexadecimal <- function(in_dat) {
    apply(in_dat, 2,
          function(x) {
            # Make minimum 0
            x - min(x)
          }) %>%
      apply(2,
            function(x) {
              # Constrain to range [0, 255]
              round(255 * (x / max(x)))
            }) %>%
      apply(1,
            function(x) {
              # Convert to hexadecimal codes
              toupper(paste0("#", as_hex(x[1]), as_hex(x[2]), "00"))
            })
  }
  
  blend_plot_theme <- theme(legend.position="none",
                            plot.title=element_text(hjust=0.5))
  
  plot_list <- lapply(c(column_1, column_2),
                      function(column) {
                        max_color <- if_else(column == column_1,
                                            "#FF0000", "#00FF00")
                        SpatialFeaturePlot(cells_obj, column) +
                          scale_fill_gradient(low="#000000",
                                              high=max_color) +
                          ggtitle(column) +
                          blend_plot_theme
                      })
  
  dat <- FetchData(cells_obj, c(column_1, column_2))
  colors <- as.matrix(dat) %>% metadata_to_hexadecimal()
  
  new_md_column <- paste0(column_1, "_vs_", column_2)
  cells_obj[[new_md_column]] <- colors
  names(colors) <- as.character(colors)
  
  plot_list[[3]] <- SpatialDimPlot(cells_obj, new_md_column, cols=colors) +
    ggtitle(paste0(column_1, "_", column_2)) +
    blend_plot_theme
  
  side_length <- 100
  legend_grid <- expand.grid(seq(from=min(dat[, column_1]),
                                 to=max(dat[, column_1]),
                                 length.out=side_length),
                             seq(from=min(dat[, column_2]),
                                 to=max(dat[, column_2]),
                                 length.out=side_length))
  colnames(legend_grid) <- c(column_1, column_2)
  legend_colors <- metadata_to_hexadecimal(legend_grid)
  legend_grid$color <- legend_colors
  names(legend_colors) <- legend_colors
  
  legend <- ggplot(legend_grid,
                   aes(x=.data[[column_1]], y=.data[[column_2]],
                       color=color)) +
    geom_point(shape=15, size=1.9) +
    scale_color_manual(values=legend_colors) +
    coord_cartesian(expand=FALSE) +
    theme(legend.position="none", aspect.ratio=1,
          panel.background=element_blank())
  
  plot_list[[4]] <- wrap_plots(ggplot() + theme_void(), legend,
                               ggplot() + theme_void(), ncol=1,
                               heights=c(0.2, 0.6, 0.2))
  
  if (combine == FALSE) {
    return(plot_list)
  } else {
    p <- wrap_plots(plot_list, nrow=1,
                    widths=c(0.28, 0.28, 0.28, 0.16))
    return(p)
  }
}

#******************************************************************************#
#                      RNA SEQ ANALYSIS RELATED FUNCTIONS                      #       
#******************************************************************************#

### Read individual count files.txt and merge all counts into "Read_data.xlsx"
# Input is path to counts file, project and method (STAR, HTSeq)
compile_raw_counts <- function(count_dir, proj, method){
  
  data_path <- gsub(pattern = paste0("counts/", method, "_raw_counts/"), 
                    replacement = "", x = count_dir)
  
  # Create a list of all txt files within folder that will be analyzed
  files <- list.files(path=count_dir)
  
  # Create an empty dataframe with 0's
  read_data <- data.frame(0)
  
  # Create the reads table 
  for (i in 1:length(files)){
    
    # Read the txt file
    temp_file <- read.table(file=paste0(count_dir, files[i]), header=FALSE, sep="\t")  
    
    if (method == "HTSEQ"){
      # Remove last 5 rows in HTSeq count output: __no_feature,  __ambiguous, 
      # __too_low_aQual, __not_aligned, __alignment_not_unique
      temp_file <- temp_file[1:(nrow(temp_file)-5),]
    } else if (method == "STAR"){
      # Remove top 4 rows in STAR count output:  
      # N_unmapped, N_multimapping, N_noFeature, N_ambiguous
      temp_file <- temp_file[5:nrow(temp_file),]
    }
    
    # The 1st Column will have Ensembl ids. 
    # For HTSeq, gene counts may be in 2nd or 3rd column. 
    # For STAR, gene counts are in 2nd (unstranded), 3rd (+), 4th (-) column. 
    # Append appropriate column.
    if (method == "HTSEQ" & sum(temp_file[2], na.rm=TRUE) == 0 & sum(temp_file[3], na.rm=TRUE) > 0){
      read_data <- bind_cols(read_data, temp_file[,3])
    } else if (method == "HTSEQ" &  sum(temp_file[2], na.rm=TRUE) > 0 & sum(temp_file[3], na.rm=TRUE) ==0){
      read_data <- bind_cols(read_data, temp_file[,2])                  
    } else if (method == "STAR" & abs((sum(temp_file[2])/sum(temp_file[3])) - (sum(temp_file[2])/sum(temp_file[4]))) < 2){
      print("Unstranded")
      read_data <- bind_cols(read_data, temp_file[,2])                  
    } else if (method == "STAR" & sum(temp_file[3]) > 3*sum(temp_file[4])){
      print("Pos stranded")
      read_data <- bind_cols(read_data, temp_file[,3])                  
    } else if (method == "STAR" & sum(temp_file[4]) > 3*sum(temp_file[3])){
      print("Neg stranded")
      read_data <- bind_cols(read_data, temp_file[,4])                  
    } else{
      print("Error: Gene counts NOT PRESENT in either column 2 or 3 of count file")
    }
    
    # Rename the column names to sample names
    colnames(read_data)[i+1] <- gsub(pattern="\\..*$|ReadsPerGene.out.tab", replacement="", x=files[i])
  }
  
  # Check if all count files have same order of genes in the rows so that files can be merged together
  temp_file <- read.table(file=paste0(count_dir, files[1]), header=FALSE, sep="\t")
  gene_list <- temp_file[, 1]
  for (i in 1:length(files)){
    
    temp_file <- read.table(file=paste0(count_dir, files[i]), header=FALSE, sep="\t")
    genes <- temp_file[, 1]
    
    if (!identical(gene_list, genes)){ 
      print("Gene order is different between the count files")
    }
  }
  
  if (method == "STAR"){
    gene_list <- gene_list[5:length(gene_list)]
  } else if(method == "HTSEQ"){
    gene_list <-gene_list[1:(length(gene_list)-5)]
  }
  
  # Add gene names to 1st column
  read_data[, 1] <- gene_list
  colnames(read_data)[1] <- "SYMBOL"
  colnames(read_data) <- gsub(pattern="ReadsPerGene", replacement="", x=colnames(read_data))
  
  # Remove genes with 0 counts in all samples
  read_data <- read_data[rowSums(read_data[,-1]) != 0,]
  
  # Remove samples with 0 counts in all genes
  read_data <- read_data[,colSums(read_data[,-1]) != 0]
  
  # Save the results as xlsx file
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Raw_counts")
  openxlsx::writeData(wb, sheet="Raw_counts", x=read_data)
  openxlsx::saveWorkbook(wb, file=paste0(data_path, proj, ".raw.counts.xlsx"),
                         overwrite=TRUE)
  
  return(read_data)
}

# Example: If samples were isolated on different days & prepared using different
# kits, "Batch" column must have values like "1_Ribo", "1_Poly-A", "2_Ribo", etc
# NOTE: Make sure there are no white spaces in the Target and Reference columns
# in excel file. R will change any white space (" ") to dot ("."). So, replace 
# white space (" ") with underscore ("_") before importing into R.

# DESeq2 automatically removes genes that have 0 counts but doesnt remove 
# samples that have 0 counts for all genes (unlikely scenario). So, remove
# such samples, else, the geometric mean will be 0 for all genes and DESeq2 
# will halt execution.

# Built-in PCA plot in DESeq2
# Use getMethod("plotPCA", "DESeqTransform") to understand how DESeq2 makes 
# PCA plot, it uses top 500 genes with highest variance and uses scale=FALSE
# pcaData <- plotPCA(object = vsd, intgroup = "Sample.ID", returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))

# Calculate Principal components
# (i) PCA is performed across columns. So, have variables i.e. genes on columns.
# (ii) "Error cannot rescale a constant/zero column to unit variance" is due  
# to genes that have 0 expression across all samples. Remove them.
# (iii) All NA must be replaced with 0
# (iv) prcomp() is better than princomp()

# Normalized counts are influenced by sizeFactors.
# sizeFactors are affected by number of samples (all samples vs subset of samples)
# sizeFactors are NOT affected by design formula.
# sizeFactors MUST be estimated first before normalization.
# Normalized counts from dds object are NOT batch corrected. We do this below.
# https://www.biostars.org/p/490181/

# design doesnt affect size factors. Hence, normalized counts are not affected by design
# but vst counts are affected by design blind=TRUE vs blind=FALSE

# AnnotationHub has SYMBOL-ENSEMBL_ID info ONLY.
# AnnotationDbi has SYMBOL-ENSEMBL_ID as well as SYMBOL-ENTREZ_ID info.
# hubCache(AnnotationHub()) to find location where cache is stored and delete
# it and start fresh if you get errors like "Error: failed to load resource"
# NOTE: mapIds can ONLY retrieve one of "EMSEMBL/SYMBOL/GENETYPE" at a time
# mapping <- AnnotationDbi::mapIds(x=org.Hs.eg.db, 
#                                  keys=keys(org.Hs.eg.db),
#                                  keytype="ENTREZID", 
#                                  column="SYMBOL") %>%
#   as.data.frame(do.call(cbind, list(.))) %>%
#   tibble::rownames_to_column("ENTREZID") %>%
#   dplyr::rename(ENTREZID=identity(1), SYMBOL=identity(2))

# Two types of analysis: 
# (i) Gene Set Enrichment Analysis (GSEA)
# (ii) Over Representation Analysis (ORA)
# GSEA uses all DEGs ordered by fold change or other parameter
# ORA uses ONLY significant DEGs and ignores fold change etc

# NOTE: Genes MUST be ranked i.e. sorted in descending fold change. You can 
# also rank based on log2FC & p value like: sign(df$log2fc)*(-log10(df$pval)))

# NOTE: Genes MUST be stored in list format, not as a dataframe.

# NOTE: No NA MUST be present in SYMBOL column. Else, fgsea::collapsePathways()
# will give "Error in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  : 
# Not all stats values are finite numbers"
# You can figure out using table(is.na(names(DEGs_list))), 
# is.infinite(DEGs_df$log2FoldChange) or sapply(DEGs_df, class) to make sure
# log2FoldChange and padj are numeric

# NOTE: If your excel file has "inf" in padj or log2FoldChange columns, the
# column will be read into R as character column instead of numeric. So, remove
# text values from log2FoldChange and padj columns

# Define score type in fgseaMultilevel() based on fold change values.
# NOTE: Use "pos", if you are ONLY interested in activated pathways.
# NOTE: Use "neg", if you are ONLY interested in inhibited pathways. 
# NOTE: Else, use "std" for both activated & inhibited pathways. 

# NOTE: If you run multiple gene sets like C5 and C2 together, padj will not 
# be significant as there will be too many multiple comparisons. So, run
# each gene set separately and merge results

# If you ordered your gene list based on fold change, then +ve NES indicates
# that the genes in this gene set are mostly at the top of your gene list
# (hence, most of them are upregulated) and -ve NES indicates that the genes
# in this gene set are mostly at the bottom of your gene list (hence, most 
# of them are downregulated)

# NOTE: Output of fgsea is a data.table & data.frame. 
# "leadingEdge" column is a list of genes. 
# So, DO NOT FORCE the output of fgsea to a dataframe as this will lead to 
# data loss from "leadingEdge" column & affect plotting using fgsea::plotEnrichment()

# NOTE: DO NOT USE labels for defining colors due to reasons below. 
# RECOMMEND using a named vector.
# NOTE: If using labels, sort labels in alphabetical order and then assign 
# color because R by default will arrange the labels in alphabetical order 
# first and then match them to colors indicated in values vector and then 
# color the plot. The coloring in the legend is however dependent on the 
# order of labels vector and values vector. To understand, create a plot first 
# using the sort and then without the sort(). 

### combatseq Batch correction:
# NOTE: This batch correction of known factors is done on raw counts
# The batch corrected raw reads are used in DESeq2
batch_correct_combat <- function(meta_data, read_data, formula_string){ 
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~1)
  dds <- DESeq2::estimateSizeFactors(dds) 
  dat  <- DESeq2::counts(dds, normalized = TRUE)

  # Full model matrix with the variable of interest
  mod  <- stats::model.matrix(as.formula(formula_string), colData(dds))
  
  if (length(unique(as.numeric(meta_data$Batch))) > 1){
    # Instead of using group & full_mod=TRUE, use covar_mod
    read_data_combat <- sva::ComBat_seq(counts=as.matrix(read_data), 
                                        batch=as.numeric(meta_data$Batch), 
                                        #group=as.numeric(as.factor(meta_data$Condition)),
                                        #full_mod = TRUE,
                                        group = NULL,
                                        covar_mod = mod)
  } else{
    read_data_combat <- read_data
  }
  return(read_data_combat) 
}

### svaseq Batch correction:
# NOTE: This batch correction of unknown factors is done on normalized counts
# NOTE: svaseq() can find n number of surrogate variables. If we model for all 
# of them there could be over correction. Hence, we limit batch correction to
# only the top 3 surrogate variables.
# Here, we just create a new object sva_dds with sva design variables
batch_correct_sva <- function(meta_data, read_data, formula_string){
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~1)
  dds <- DESeq2::estimateSizeFactors(dds) 
  dat  <- DESeq2::counts(dds, normalized = TRUE)
  
  # Remove lowly expressed genes before finding surrogate variables
  idx  <- rowMeans(dat) > 1
  dat  <- dat[idx, ]
  
  # Full model matrix with the variable of interest #Reduced/null model matrix with only an intercept term
  mod  <- stats::model.matrix(as.formula(formula_string), colData(dds))
  # Reduced/null model matrix with only an intercept term
  mod0 <- stats::model.matrix(~ 1, colData(dds))
  # Estimate all surrogate variables
  svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0)
  
  # If there are more than 3 SV, estimate upto 3 SVs
  # svseq$sv values will change based on n.sv
  if (svseq$n.sv > 3){
    svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0, n.sv=3)
  }
  
  # Add these SVs as columns to the DESeqDataSet and then add them to the design
  # ddssva$SV1 <- svseq$sv[,1]
  # ddssva$SV2 <- svseq$sv[,2]
  # design(ddssva) <- ~ SV1 + SV2 + id
  
  ddssva <- dds
  for (i in 1:ncol(svseq$sv)){
    var <- paste0("SV",i)
    ddssva[[var]] <- svseq$sv[,i]
  }
  
  design(ddssva) <- as.formula(paste0("~", 
                                      paste0("SV", seq(1:ncol(svseq$sv)), collapse = "+"), 
                                      "+", 
                                      gsub(pattern="~", replacement="",x=formula_string)))
  
  return(ddssva)
}

norm_counts_combat <- function(meta_data, read_data_combat, output_path){
  
  # design doesnt affect size factors. Hence, normalized counts are not affected by design
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data_combat,
                                        colData=meta_data, 
                                        design=~ 1)
  # Estimate sizefactors
  dds <- DESeq2::estimateSizeFactors(object = dds)
  
  # EXtract normalized counts from dds object
  normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
  
  normalized_counts <- normalized_counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID")
  
  # Add gene names
  normalized_counts <- add_annotation(normalized_counts, annotations)
  
  # Save the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Norm_counts")
  openxlsx::writeData(wb, sheet="Norm_counts", x=normalized_counts, rowNames=FALSE)
  openxlsx::saveWorkbook(wb, file=paste0(output_path, "Normalized.counts.combat.xlsx"), 
                         overwrite=TRUE)
  
  return(normalized_counts)
}

# svaseq corrected normalized counts:
# NOTE: ddssva object from svaseq_batch has the top 2 surrogate variables that 
# will be used in DESeq2() but the normalized counts from ddssva object are NOT 
# batch corrected. We do this below.  https://www.biostars.org/p/121489/
# NOTE: Lowly expressed genes are removed before finding surrogate variables.
# So, number of genes is lower than number of DESeq2 normalized counts excel.
norm_counts_sva <- function(meta_data, read_data, formula_string, output_path){
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~1)
  dds <- DESeq2::estimateSizeFactors(dds) 
  dat  <- DESeq2::counts(dds, normalized = TRUE)
  
  # Remove lowly expressed genes before finding surrogate variables
  idx  <- rowMeans(dat) > 1
  dat  <- dat[idx, ]
  
  # Full model matrix with the variable of interest #Reduced/null model matrix with only an intercept term
  mod  <- stats::model.matrix(as.formula(formula_string), colData(dds))
  # Reduced/null model matrix with only an intercept term
  mod0 <- stats::model.matrix(~ 1, colData(dds))
  # Estimate all surrogate variables
  svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0)
  
  # If there are more than 3 SV, estimate upto 3 SVs
  # svseq$sv values will change based on n.sv
  if (svseq$n.sv > 3){
    svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0, n.sv=3)
  }
  
  # %*% indicates Matrix multiplication
  X <- base::cbind(mod, svseq$sv)
  Hat <- base::solve(t(X) %*% X) %*% t(X)
  beta <- (Hat %*% t(dat))
  P <- ncol(mod)
  corrected_data <- dat - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])
  
  normalized_counts <- corrected_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID")
  
  normalized_counts <- add_annotation(normalized_counts, annotations) 
  
  # Save batch corrected normalized counts for entire dataset
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Norm_counts_batch_corrected")
  openxlsx::writeData(wb, sheet="Norm_counts_batch_corrected", x=normalized_counts_batch, rowNames=FALSE)
  openxlsx::saveWorkbook(wb, file=paste0(output_path, "Normalized.counts.SVA.xlsx"), 
                         overwrite=TRUE)
}

### Perform DESeq2
run_deseq2 <- function(dds, meta_data, DEG.params, n, approach, output_path){  
  # betaPrior: default=FALSE, shrunken LFCs are obtained later using lfcShrink
  # Define all possible coefficient vectors
  # Get the model matrix
  mod_mat <- model.matrix(design(dds), colData(dds))
  
  # Identify the design factors
  design.factors <- str_split(string=DEG.params$design, pattern="\\+|:|\\*")[[1]] %>% unique()
  
  # Create a replicate of meta_data with all possible groups that could be compared based on the design
  df <- colData(dds) %>% 
    data.frame() %>% 
    tidyr::unite(col = "Groups", design.factors, sep=".")
  
  # Get all possible coefficent vectors
  groups <- df$Groups %>% unique()
  for (i in groups){
    val <- colMeans(mod_mat[df$Groups == i,])
    assign(x = i, value = val)
  }
  
  # Define the contrast
  contrast <- base::eval(expr = base::parse(text = DEG.params$contrast[n]))
  
  # # NOTE: coeff MUST match one of columns in resultsNames(dds)
  # # Define contrast and coeff [NOT RECOMMENDED]
  # DE_levels <- meta_data %>%
  #   dplyr::select(all_of(DEG.params$Variable[n])) %>%
  #   unlist(use.names=FALSE) %>%
  #   unique() %>%
  #   as.vector()
  # 
  # contrast <- c(DEG.params$Variable[n],
  #               DE_levels[DE_levels==DEG.params$Target[n]],
  #               DE_levels[DE_levels==DEG.params$Reference[n]])
  
  #coeff <- paste(contrast[1], contrast[2], "vs", contrast[3], sep='_')

  # Perform lfcshrinkage to account for variability between replicates
  # For ashr, if res is provided, then coef and contrast are ignored.
  # lfcshrinkage will not change the number of DEGs and affects only logFC
  
}

#******************************************************************************#
#                      PATHWAY ANALYSIS RELATED FUNCTIONS                      #
#******************************************************************************#

plot_ora <- function(ora_results, file_suffix, output_path){
  
  # Remove under score from pathway names
  ora_results <- ora_results %>%
    dplyr::mutate(Description = gsub(pattern="_", replacement=" ", x=Description))
  #Description = gsub(pattern="GOBP |GOCC |GOMF |REACTOME ", replacement="GOBP:|GOCC:|GOMF:|REACTOME:",x=Description))
  
  # Plot bar plots
  # (NOT RECOMMEDNED since it gives info only on enrichment ratio & pvalue)
  ggplot2::ggplot(data = ora_results,
                  aes(x = k.K,
                      y = reorder(Description, k.K),
                      #fill = p.adjust,
                      fill = pvalue)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Enrichment Ratio",
                  y = "",
                  fill = "padj",
                  title = "Pathways") +
    geom_text(aes(label = Count), position = position_dodge(width = 1), hjust = -0.5) +
    #ggplot2::coord_cartesian(xlim = c(0, max(abs(enriched_result$k.K)))) +
    ggplot2::theme(aspect.ratio = 2,
                   plot.title =   element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                   plot.caption = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0),
                   axis.title.x = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                   axis.title.y = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                   axis.text.x =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                   axis.text.y =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 1),
                   legend.title = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0, vjust = 1),
                   legend.text =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                   #legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
                   legend.position = "bottom", #"right",
                   legend.justification = "center",
                   legend.direction = "horizontal", #"vertical",
                   legend.key.height = unit(0.5, 'cm'),     #unit(0.75, 'cm'),
                   legend.key.width = unit(1.25, 'cm')) +   #unit(0.5, 'cm'),
    viridis::scale_fill_viridis(option = "viridis", limits = c(0, 0.05))
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("ORA_Bar_Plot_", file_suffix, ".tiff"),
                  plot = last_plot(),
                  device = "jpeg",
                  path = output_path,
                  scale = 1,
                  width = 12,
                  height = 12,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
  # Plot dot plots (RECOMMENDED)
  ggplot2::ggplot(data = ora_results,
                  aes(x = k.K,
                      y = reorder(Description, k.K),
                      #label = Upstream.Regulator,
                      #color = p.adjust,
                      color = pvalue,
                      size = Count)) +
    ggplot2::geom_point() +
    ggplot2::scale_size_continuous(range = c(10,20)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Enrichment Ratio",
                  y = "",
                  title = "",
                  fill = "padj") +
    ggplot2::theme(aspect.ratio = 2,
                   plot.title =   element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                   plot.caption = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0),
                   axis.title.x = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                   axis.title.y = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                   axis.text.x =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                   axis.text.y =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 1),
                   legend.title = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0, vjust = 1),
                   legend.text =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                   #legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
                   legend.position = "right",
                   legend.justification = "left",
                   legend.direction = "vertical",
                   legend.key.height = unit(0.5, 'cm'),
                   legend.key.width = unit(0.5, 'cm')) +
    viridis::scale_color_viridis(option = "viridis") +
    ggplot2::scale_size(breaks = sapply(as.vector(quantile(c(min(ora_results$Count), max(ora_results$Count)))), floor))
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("ORA_Dot_Plot_", file_suffix, ".tiff"),
                  plot = last_plot(),
                  device = "jpeg",
                  path = output_path,
                  scale = 1,
                  width = 12,
                  height = 8,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
}

#******************************************************************************#
#                                 VENN DIAGRAM                                 #
#******************************************************************************#

#' Plot Venn Diagram and Save Overlapping Genes
#'
#' This function takes a data frame of gene lists (one list per column upto 4 columns),
#' plots a Venn diagram based on the gene overlaps, saves the diagram as a TIFF file,
#' and exports an Excel workbook containing the overlapping genes and the original data.
#'
#' @param data A data frame where each column is a list/vector of gene names. Each column represents a group.
#'             NAs are automatically removed from each column.
#' @param path A character string specifying the directory path where output files will be saved.
#'             The directory must exist.
#' @param suffix A character string used to label output files and plot titles.
#'
#' @return Invisibly returns NULL. Side effects: saves Venn diagram TIFF and Excel file to `path`.
#'
#' @details
#' The function supports 1 to 4 columns in `data`.
#' It wraps column names for better display on the plot and applies color palettes
#' depending on the number of groups.
#' Overlapping genes are calculated and saved in an Excel workbook with two sheets:
#' one for overlaps and one with the input data.
#' @export
plot_venn <- function(data, path, suffix){
  
  # Input validation
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (!dir.exists(path)) stop("`path` directory does not exist.")
  if (!is.character(suffix) || length(suffix) != 1) stop("`suffix` must be a single string.")
  
  plot_title <- suffix
  ncol <- ncol(data)
  
  # Clean column names (replace _ and . with space)
  colnames(data) <- stringr::str_replace_all(colnames(data), c("_" = " ", "\\." = " "))
  
  # Set cat.pos, cat.dist, cex and palette based on number of columns
  if (ncol == 4){
    pos <- c(330, 15, 330, 15)
    dist <- c(0.27, 0.25, 0.15, 0.13)
    cex = 2
    palette1 <- c("#C8E7F5", "#00008C", "#F6D2E0", "#E75480")        
  } else if (ncol == 3){
    pos <- c(0, 0, 180)
    dist <- c(0.1, 0.1, 0.1)
    cex = 2
    palette1 <- c("#C8E7F5", "#F6D2E0", "#db6d00")     
  } else if (ncol == 2){
    pos <- c(0, 0)
    dist <- c(0.05, 0.05)
    cex = 2.75
    palette1 <- c("#C8E7F5", "#db6d00")                             
  } else if (ncol == 1){
    pos <- c(0)
    dist <- c(0.1)
    cex = 2.75
    palette1 <- c("#F6D2E0")                                         
  } else {
    stop("`data` must have between 1 and 4 columns.")
  }
  
  # Create a dataframe to store the wrapped column names
  annotation <- data.frame(Labels = stringr::str_wrap(colnames(data), width = 10))
  
  # Convert the data frame to a named list (removing NAs)
  genes <- base::vector(mode = "list", length = ncol(data))
  names(genes) <- annotation$Labels
  
  for (i in 1:ncol(data)){
    
    # remove NA values and create a list of genes for each label
    genes[[i]] <- data[!is.na(data[i]),i]
  }
  
  # Plot the venn diagram
  VennDiagram::venn.diagram(x = genes,
                            main = plot_title, 
                            category.names = annotation$Labels,
                            filename = file.path(path, paste0("Venn_Diagram_", suffix, ".tiff")),
                            output = TRUE,
                            scaled =FALSE,
                            imagetype = "tiff",
                            height = 11, 
                            width = 11,
                            units = "in",
                            resolution = 600,
                            compression = "lzw",
                            #amount of white space around Venn Diagram in grid units
                            margin = 0.3,    
                            
                            # Formatting the shapes of venn diagram
                            lwd = 1.5,                 #thickness of line
                            lty = 1,                   #type of line
                            col = "black",             #color of line
                            
                            # Formatting numbers inside venn diagram
                            cex = cex,                 #font size (2 or 2.75)
                            fontface = "bold",         #font style
                            fontfamily = "sans",       #font type
                            
                            # Formatting title of venn diagram
                            main.cex = 2,              #font size
                            main.fontface = "bold",    #font style
                            main.fontfamily = "sans",  #font type
                            main.col = "black",        #font color
                            
                            # Formatting category of venn diagram
                            cat.cex = 2,               #font size
                            cat.fontface = "bold",     #font style
                            cat.fontfamily = "sans",   #font type
                            cat.col = palette1,  #"black",
                            
                            # Formatting colors of venn diagram
                            fill = palette1,
                            alpha = rep(0.5, ncol), #0.5=50% transparency, 1=0% transparency
                            #cat.default.pos = "outer",    
                            
                            cat.pos = pos,    
                            cat.dist = dist, 
                            disable.logging = TRUE,
                            ext.text = TRUE)
  
  #******************************************************************************#
  #                          SAVE THE OVERLAPPING GENES                          #
  #******************************************************************************#
  
  # Save the list of overlapping genes. NOTE: You need to manually figure out
  # which genes belong to which overlap based on number of genes overlapping
  
  # Calculate overlaps
  overlap <- VennDiagram::calculate.overlap(x = genes)
  
  # Identify maximum number of genes present in any overlap
  max_len = max(lengths(overlap))
  
  # Create an dataframe of size length(overlap), max with NAs
  results = data.frame(matrix("", nrow = max_len, ncol = length(overlap)))
  rownames(results) <- paste0("Gene#", seq(max_len))
  colnames(results) <- paste0("Intersection#", seq(length(overlap)))
  
  # Populate the dataframe with gene names
  for (i in 1:length(overlap)){
    if (length(overlap[[i]]) > 0){
      for (j in 1:length(overlap[[i]])){
        results[j,i] <- overlap[[i]][j]
        #results[[j,i]] <- overlap[[i]][j]
      }
    }
  }
  
  # Save results to Excel
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Output")
  openxlsx::writeData(wb, sheet = "Output", x = results, rowNames = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "Input")
  openxlsx::writeData(wb, sheet = "Input", x = data, rowNames = FALSE)
  openxlsx::saveWorkbook(wb, file = file.path(path, paste0("Overlap_", suffix, ".xlsx")),
                         overwrite = TRUE,  returnValue = FALSE)
  
  invisible(NULL)
}

#******************************************************************************#
#                                    HEATMAP                                   #
#******************************************************************************#

#' Plot Heatmap of Normalized Gene Expression
#'
#' Generates a heatmap from normalized gene expression counts with optional 
#' row and column annotations, clustering, scaling, and color customization.
#' Saves the processed matrix and heatmap plot to specified output.
#' 
#' @param norm_counts A data frame containing normalized expression counts.
#'   Must include a column named \code{SYMBOL} with gene identifiers 
#'   (ENSEMBL_ID, ENTREZ_ID, or SYMBOL). Column names (except \code{SYMBOL}) 
#'   correspond to sample IDs.
#'   If genes are duplicated, only the row with the highest total expression 
#'   across samples is retained.
#' @param metadata_column A data frame of sample metadata. Must include a 
#'   \code{Sample.ID} column matching sample names in \code{norm_counts}, and 
#'   any columns specified in \code{heatmap.params$anno.column} for annotation.
#' @param metadata_row A data frame of gene metadata. Must include a \code{SYMBOL} 
#'   column matching genes in \code{norm_counts}, and any columns specified in 
#'   \code{heatmap.params$anno.row} for annotation.
#' @param heatmap.params A list of parameters controlling heatmap behavior:
#'   \describe{
#'     \item{log.transform}{Logical. Whether to log2-transform the counts (default FALSE).}
#'     \item{scale}{Logical. Whether to scale gene expression across samples (default FALSE).}
#'     \item{anno.column}{Character vector of metadata column names to use as column annotations.}
#'     \item{anno.row}{Character vector of metadata column names to use as row annotations.}
#'     \item{col.split}{Column name in \code{metadata_column} to define column gaps/clusters.}
#'     \item{row.split}{Column name in \code{metadata_row} to define row gaps/clusters.}
#'     \item{row.cluster}{One of \code{"group"}, \code{"alphabetical"}, or \code{"all"} for row ordering/clustering.}
#'     \item{col.cluster}{One of \code{"group"}, \code{"alphabetical"}, or \code{"all"} for column ordering/clustering.}
#'     \item{matrix_color}{Name of color palette function for heatmap colors (e.g., \code{viridis_pal}).}
#'     \item{border_color}{Color for heatmap cell borders.}
#'     \item{bar_width, bar_height}{Numeric. Width and height of heatmap cells.}
#'     \item{expr_legend}{Logical. Whether to show expression legend.}
#'     \item{angle}{Numeric or character. Rotation angle for column labels.}
#'     \item{width, height}{Numeric. Dimensions of the output heatmap file.}
#'     \item{file_format}{Character. Output file format (e.g., \code{"pdf"}, \code{"png"}).}
#'     \item{discrete_panel}{Logical. Whether to use discrete color panels for annotations.}
#'   }
#' @param plot_genes Character vector of gene identifiers to include in the heatmap.
#' @param disp_genes Character vector of gene names to be labeled on the heatmap.
#' @param file_suffix Character string appended to output file names.
#' @param output_path Character string specifying the directory where output files 
#'   (heatmap images and data matrices) are saved.
#' 
#' @details
#' - The function filters \code{norm_counts} to genes in \code{plot_genes} (case-insensitive).
#' - Genes with zero expression across all samples are removed.
#' - Duplicate genes retain only the row with the highest summed expression.
#' - The data matrix is optionally log-transformed and scaled.
#' - Annotations for samples (columns) and genes (rows) are added from \code{metadata_column} and \code{metadata_row}.
#' - Row and column clustering/order can be controlled via \code{heatmap.params}.
#' - The heatmap is saved as a file, and the processed matrix is saved as an Excel file.
#' 
#' @return
#' Invisibly returns \code{NULL}. Side effects include saving the heatmap plot 
#' and data matrix to disk.
#' @export
plot_heatmap <- function(norm_counts, metadata_column, metadata_row, heatmap.params,
                         plot_genes, disp_genes, file_suffix, output_path){
  
  
  # NOTE: If you are generating heatmaps using gene counts, you MUST USE
  # vst transformed counts. These are already log transformed. So, just scale
  # each gene across samples
  
  #****************************************************************************#
  # Format matrix for heatmap: filter genes, handle duplicates, transform, scale
  #****************************************************************************#
  
  mat <- norm_counts %>%
    # Keep only genes that need to be plotted
    dplyr::filter(base::toupper(SYMBOL) %in% base::toupper(plot_genes)) %>%
    # Replace NA with 0
    base::replace(is.na(.), 0) %>%
    # If there are duplicated genes, keep only row for highest expressing copy
    dplyr::mutate(n = rowSums(.[,-1])) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(n) %>%
    dplyr::ungroup() %>%
    # Remove genes with 0 expression in all samples
    dplyr::filter(n != 0) %>%
    dplyr::select(everything(), -n) %>%
    # Move gene names to rownames
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("SYMBOL") %>%
    # Make sure all columns are numeric
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
  # mat[, unlist(lapply(mat, is.numeric))]    #alternative to mutate(across())
  
  # Make rownames proper
  rownames(mat) <- make.names(rownames(mat), unique = TRUE)
  colnames(mat) <- make.names(colnames(mat), unique = TRUE)
  
  # Perform log transform if needed. Count data is usually skewed right or left.
  # So, it will mostly be red or blue. log transform to make it less skewed.
  if (heatmap.params$log.transform == TRUE){
    mat <- log(1+mat, base = 2)
  }
  
  # Currently, genes are in rows and samples are in columns.
  # We need to scale each gene (in rows) across samples (in columns)
  # However, scale() can ONLY perform scaling within each column.
  # So, we transpose the dataframe first, so that genes are in columns & samples
  # are in rows & then perform scaling.
  if (heatmap.params$scale == TRUE){
    mat <- mat %>% t() %>% scale() %>% t()
  }
  
  # After scaling, NA values could have been introduced. So, replace NA with 0
  mat[is.na(mat)] <- 0
  
  # Keep only samples common to metadata and mat
  metadata_column <- metadata_column %>% 
    dplyr::filter(make.names(Sample.ID) %in% make.names(colnames(mat)))
  
  # Arrange samples in mat in the same order as in metadata_column.
  # NOTE: This is important because in the next step we assign rownames to
  # col_annotation assuming identical sample order between metadata_column & mat
  #mat <- mat[,metadata_column[,columns]]
  
  #****************************************************************************#
  # Define column and row annotations based on metadata and heatmap params
  #****************************************************************************#
  
  # Define column annotation for samples
  if(length(heatmap.params$anno.column) > 0){
    col_annotation <- metadata_column %>%
      dplyr::select(Sample.ID, all_of(heatmap.params$anno.column)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("Sample.ID")
    rownames(col_annotation) <- make.names(rownames(col_annotation))
  } else {
    col_annotation <- NA
  }
  
  #gtools::invalid(heatmap.params$anno.column))
  
  # Define row annotation for genes
  if(length(heatmap.params$anno.row) > 0){
    row_annotation <- metadata_row %>%
      dplyr::select(SYMBOL, all_of(heatmap.params$anno.row)) %>%
      dplyr::mutate(SYMBOL = make.names(SYMBOL, unique=TRUE)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("SYMBOL") %>%
      as.data.frame()
    rownames(row_annotation) <- make.names(rownames(row_annotation))
  } else {
    # If you define as NA, etc, it will count as 1 variable during annotation_colors calculation
    row_annotation <- NULL  
  }
  
  #****************************************************************************#
  # Generate annotation colors (column + row) based on unique factor levels
  #****************************************************************************#
  
  # NOTE: There is ONLY 1 parameter in pheatmap() i.e. annotation_colors to 
  # define colors for both row and column annotations
  
  # The colors should be specified in the following format:
  # $Sample
  # FB1       FB2       FB3       FB4       FB5        FC       MB1
  # "#BF812D" "#35978F" "#C51B7D" "#7FBC41" "#762A83" "#E08214" "#542788"
  # $Sex
  # Female      Male
  # "#9E0142" "#E41A1C"
  
  # This is an example of how this needs to be specified
  # ann_colors <- list(Column_Group1 = c(`Immune Depleted` = "#CB181D", `Immune Enriched` = "#A6D854"),
  #                    Row_Group1 = c(`Pro-tumor` = "white", `Anti-tumor` = "white"))
  
  ann_colors <- list()
  colors <- c("#E08214", "#762A83", "#C51B7D", "#7FBC41", "#35978F",
              "#BF812D", "#542788", "#D6604D", "#4393C3", "#878787",
              "#E41A1C", "#F781BF", "#4DAF4A", "#FFFFBF", "#377EB8",
              "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#999999",
              "#66C2A5", "#FC8D62", "#000000", "#9E0142", "#1A1A1A",
              "#74006F", "#FFC606", "#F6D2E0", "#C8E7F5")
  
  # Convert col_annotation and row_annotation to list for easier analysis
  col_list <- base::lapply(X = as.list(col_annotation), FUN = unique)
  row_list <- base::lapply(X = as.list(row_annotation), FUN = unique)
  ann_list <- c(col_list, row_list)
  
  # Define total number of variables to be colored in row and column annotation
  ann_total <- length(ann_list)
  
  # Go through each variable and specify colors for each element within variable
  if (ann_total > 1){
    start_n <- 1
    end_n <- 0
    for (i in 1:ann_total){
      
      elements <- ann_list[[i]]
      palette_n <- length(elements)
      # Create a vector of transparency values. Remove alpha=0 as it is always white
      alphas <- base::seq(from=0, to=1, by=1/palette_n)
      alphas <- base::setdiff(alphas, 0)
      
      if (palette_n > 1 & heatmap.params$discrete_panel == FALSE){
        
        # Create a color palette with different transparencies
        palette_colors <- colorspace::adjust_transparency(col = colors[i], alpha = alphas[1])
        for (n in 2:palette_n){
          palette_colors <- c(palette_colors,
                              colorspace::adjust_transparency(col = colors[i], alpha = alphas[n]))
        }
      } else if (palette_n == 1 & heatmap.params$discrete_panel == FALSE) {
        palette_colors <- colorspace::adjust_transparency(col = colors[i], alpha = alphas[1])
      } else{
        end_n <- start_n + palette_n - 1
        palette_colors <- colors[start_n:end_n]
        start_n <- end_n+1
      }
      # Sort the elements within each variable
      names(palette_colors) <- sort(elements)
      
      # Merge all annotation colors  
      ann_colors <- c(ann_colors, list(palette_colors))
    }
    # Specify variable names to annotation colors
    names(ann_colors) <- names(ann_list)
    
  } else if (ann_total == 1){
    elements <- ann_list[[1]]
    palette_n <- length(elements)
    palette_colors <- colors[1:palette_n]
    names(palette_colors) <- sort(elements)
    ann_colors <- c(ann_colors, list(palette_colors))
    names(ann_colors) <- names(ann_list)
  } else {
    cat("There are no row or column annotations")
  }
  
  #ann_colors_col$Score[[1]] <- "#0c2c84"
  #ann_colors_col$Score[[2]] <- "#d73027"
  
  #****************************************************************************#
  # Define heatmap color scale and breaks
  #****************************************************************************#
  
  vrds <- viridis_pal()(100)
  rdbu <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
  
  # NOTE: breaks correspond to numerical ranges for color palette's bins
  # i.e. 0 to length(my_palette)
  
  # Define breaks
  if(max(mat) == 0){
    breaks <- c(seq(from = floor(min(mat)), to = 0, length.out = 100))
    my_palette <- my_palette[1:50]
  } else if (min(mat) == 0){
    breaks <- c(seq(from = 0, to = ceiling(max(mat)), length.out = 100))
    my_palette <- my_palette[50:100]
  } else if(min(mat) < -3 | max(mat) > 3){
    breaks <- c(seq(-3, 0, length.out = 50), seq(3/100, 3, length.out = 50))
    my_palette <- my_palette[1:100]
  } else{
    breaks <- c(seq(from = floor(min(mat)), to = 0, length.out = 50), seq(from = max(mat)/100, to = ceiling(max(mat)), length.out = 50))
    my_palette <- my_palette[1:100]
  }
  
  #****************************************************************************#
  # Define gaps for rows and columns based on grouping splits
  #****************************************************************************#
  
  # NOTE: count() automatically arranges by alphabetical order. So, 
  
  # Define gaps for columns
  if (!gtools::invalid(heatmap.params$col.split)){
    gaps_col <- col_annotation %>% 
      dplyr::count(get(heatmap.params$col.split)) %>%
      dplyr::mutate(n = cumsum(n)) %>% 
      dplyr::pull(n) %>%
      .[. < ncol(mat)]
  } else {
    gaps_col <- NULL
  }
  
  # Define gaps for rows
  if (!gtools::invalid(heatmap.params$row.split)){
    gaps_row <- row_annotation %>% 
      dplyr::count(get(heatmap.params$row.split)) %>%
      dplyr::mutate(n = cumsum(n)) %>% 
      dplyr::pull(n) %>%
      .[. < nrow(mat)]
  } else {
    gaps_row <- NULL
  }
  
  #   # Count() automatically arranges by alphabetical order unfortunately
  #   # So, we do this manually.
  #   element_names <- c()
  #   element_counts <- c()
  #   c <- 0
  #   for (i in 1:nrow(col_annotation)){
  #     if (!(col_annotation[,gap_columns][i] %in% element_names)){
  #       element_names <- c(element_names, col_annotation[,gap_columns][i])
  #       element_counts <- c(element_counts, c)
  #       c <- 1
  #     } else{
  #       c <- c+1
  #     }
  #   }
  #   element_counts <- c(element_counts,c)
  
  #****************************************************************************#
  # Define ordering for rows and columns according to clustering parameters
  #****************************************************************************#
  
  # reorder_axis <- function(annotation, mat_sub, split_col, cluster_method, axis = c("row","col")){
  #   axis <- match.arg(axis)
  #   elements <- annotation %>% dplyr::select(all_of(split_col)) %>% unique() %>% unlist() %>% sort()
  #   order_vec <- c()
  #   for(g in elements){
  #     items <- rownames(annotation)[annotation[[split_col]] == g]
  #     sub_mat <- if(axis == "row") mat_sub[items, , drop=FALSE] else mat_sub[, items, drop=FALSE]
  #     if(cluster_method == "group" && length(items) > 1){
  #       hc <- if(axis == "row") hclust(dist(sub_mat)) else hclust(dist(t(sub_mat)))
  #       ord <- if(axis == "row") rownames(sub_mat)[hc$order] else colnames(sub_mat)[hc$order]
  #       order_vec <- c(order_vec, ord)
  #     } else if(cluster_method == "group" && length(items) == 1){
  #       order_vec <- c(order_vec, items)
  #     } else if(cluster_method == "alphabetical"){
  #       order_vec <- c(order_vec, sort(items))
  #     } else if(cluster_method == "all"){
  #       order_vec <- c(order_vec, items)
  #       message(sprintf("No %s clustering performed; %s.split will be affected by %s.cluster", axis, axis, axis))
  #     } else {
  #       stop(sprintf("%s.cluster must be one of 'group', 'all', or 'alphabetical'", axis))
  #     }
  #   }
  #   order_vec
  # }
  # 
  # # Determine row order
  # if (!gtools::invalid(heatmap.params$row.split)) {
  #   row_order <- reorder_axis(row_annotation, mat, heatmap.params$row.split, heatmap.params$row.cluster, axis="row")
  # } else {
  #   if (heatmap.params$row.cluster == "alphabetical") {
  #     row_order <- sort(rownames(mat))
  #   } else if (heatmap.params$row.cluster == "all") {
  #     row_order <- rownames(mat)[hclust(dist(mat))$order]
  #   } else if (heatmap.params$row.cluster == "group") {
  #     row_order <- reorder_axis(row_annotation, mat, heatmap.params$anno.row[1], "group", axis="row")
  #   } else {
  #     stop("row.cluster must be 'group', 'all' or 'alphabetical'")
  #   }
  # }
  # 
  # # Determine column order
  # if (!gtools::invalid(heatmap.params$col.split)) {
  #   col_order <- reorder_axis(col_annotation, mat, heatmap.params$col.split, heatmap.params$col.cluster, axis="col")
  # } else {
  #   if (heatmap.params$col.cluster == "alphabetical") {
  #     col_order <- sort(colnames(mat))
  #   } else if (heatmap.params$col.cluster == "all") {
  #     col_order <- colnames(mat)[hclust(dist(t(mat)))$order]
  #   } else if (heatmap.params$col.cluster == "group") {
  #     col_order <- reorder_axis(col_annotation, mat, heatmap.params$anno.column[1], "group", axis="col")
  #   } else {
  #     stop("col.cluster must be 'group', 'all' or 'alphabetical'")
  #   }
  # }
  
  if (!gtools::invalid(heatmap.params$row.split)){
    
    # Important to sort so row_elements is similar to gaps_row <- row_annotation %>% dplyr::count(get(heatmap.params$anno.row[1])) %>% dplyr::mutate(n = cumsum(n))
    # Else gaps will be wrong sometimes
    row_elements <- row_annotation %>% 
      dplyr::select(all_of(heatmap.params$row.split)) %>% 
      unique() %>% 
      unlist(use.names=FALSE) %>% 
      sort()
    
    row_order <- c()
    for (g in row_elements){
      items <- rownames(row_annotation)[which(row_annotation %>% dplyr::select(all_of(heatmap.params$row.split)) == g)]
      temp_mat <- mat[items,]
      
      if (heatmap.params$row.cluster == "group" & length(items) > 1){
        rowclust <- hclust(dist(temp_mat))
        row_order <- c(row_order, rownames(temp_mat[rowclust$order,]))
      } else if(heatmap.params$row.cluster == "group" & length(items) == 1){
        row_order <- c(row_order, items)
      }else if (heatmap.params$row.cluster == "alphabetical"){
        row_order <- c(row_order, sort(rownames(temp_mat)))
      } else if (heatmap.params$row.cluster == "all"){
        row_order <- c(row_order, rownames(temp_mat))
        cat("No row clustering performed. row.split will be affected by row.cluster")
      } else {
        cat("row.cluster must be either 'group', 'all' or 'alphabetical'")
      }
    }
  } else if (gtools::invalid(heatmap.params$row.split)){
    
    if (heatmap.params$row.cluster == "alphabetical"){
      row_order <- sort(rownames(mat))
    } else if (heatmap.params$row.cluster == "all"){
      rowclust <- hclust(dist(mat))
      row_order <- rownames(mat[rowclust$order,])
    } else if (heatmap.params$row.cluster == "group"){
      
      # Important to sort so row_elements is similar to gaps_row <- row_annotation %>% dplyr::count(get(heatmap.params$anno.row[1])) %>% dplyr::mutate(n = cumsum(n))
      # Else gaps will be wrong sometimes
      row_elements <- row_annotation %>%
        dplyr::select(all_of(heatmap.params$anno.row[1])) %>%
        unique() %>%
        unlist(use.names=FALSE) %>%
        sort()
      
      row_order <- c()
      for (g in row_elements){
        items <- rownames(row_annotation)[which(row_annotation %>% dplyr::select(all_of(heatmap.params$anno.row[1])) == g)]
        temp_mat <- mat[items,]
        rowclust <- hclust(dist(temp_mat))
        row_order <- c(row_order, rownames(temp_mat[rowclust$order,]))
      }
    } else {
      cat("row.cluster must be either 'group', 'all' or 'alphabetical'")
    }
  }
  
  if (!gtools::invalid(heatmap.params$col.split)){
    
    # Important to sort so col_elements is similar to 
    # gaps_col <- col_annotation %>% dplyr::count(get(heatmap.params$anno.column[1])) %>% dplyr::mutate(n = cumsum(n))
    # Else gaps will be wrong sometimes
    col_elements <- col_annotation %>% 
      dplyr::select(all_of(heatmap.params$col.split)) %>%
      unique() %>% 
      unlist(use.names=FALSE) %>% 
      sort()
    
    col_order <- c()
    for (g in col_elements){
      items <- rownames(col_annotation)[which(col_annotation %>% dplyr::select(all_of(heatmap.params$col.split)) == g)]
      temp_mat <- mat[,items]
      
      if (heatmap.params$col.cluster == "group" & length(items) > 1){
        colclust <- hclust(dist(t(temp_mat)))
        col_order <- c(col_order, colnames(temp_mat[,colclust$order]))
      } else if(heatmap.params$col.cluster == "group" & length(items) == 1){
        col_order <- c(col_order, items)
      } else if (heatmap.params$col.cluster == "alphabetical"){
        col_order <- c(col_order, sort(colnames(temp_mat)))
      } else if (heatmap.params$col.cluster == "all"){
        col_order <- c(col_order, colnames(temp_mat))
        cat("No col clustering performed. col.split will be affected by col.cluster")
      } else {
        cat("col.cluster must be either 'group', 'all' or 'alphabetical'")
      }
    }
  } else if (gtools::invalid(heatmap.params$col.split)){
    
    if (heatmap.params$col.cluster == "alphabetical"){
      col_order <- sort(colnames(mat))
    } else  if(heatmap.params$col.cluster == "all"){
      colclust <- hclust(dist(t(mat)))
      col_order <- colnames(mat[,colclust$order])
    } else if (heatmap.params$col.cluster == "group"){
      
      # Important to sort so col_elements is similar to gaps_col <- col_annotation %>% dplyr::count(get(heatmap.params$anno.column[1])) %>% dplyr::mutate(n = cumsum(n))
      # Else gaps will be wrong sometimes
      col_elements <- col_annotation %>% 
        dplyr::select(all_of(heatmap.params$anno.col[1])) %>%
        unique() %>% 
        unlist(use.names=FALSE) %>% 
        sort()
      
      col_order <- c()
      for (g in col_elements){
        items <- rownames(col_annotation)[which(col_annotation %>% dplyr::select(all_of(heatmap.params$anno.col[1])) == g)]
        temp_mat <- mat[,items]
        colclust <- hclust(dist(t(temp_mat)))
        col_order <- c(col_order, colnames(temp_mat[,colclust$order]))
      }
    } else{
      cat("col.cluster must be either 'group', 'all' or 'alphabetical'")
    }
  }
  
  # Arrange the rows and columns
  reordered <- mat[row_order, col_order]
  
  #****************************************************************************#
  # Prepare labels for display (gene names for rows, sample names for columns)
  #****************************************************************************#
  
  display_col <- colnames(reordered)
  display_row <- dplyr::if_else(rownames(reordered) %in% disp_genes, rownames(reordered), "")
 
  #****************************************************************************#
  # Save heatmap matrix to Excel file
  #****************************************************************************#
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Heatmap_matrix")
  if (ncol(reordered) > 500){
    openxlsx::writeData(wb, sheet = "Heatmap_matrix", x = t(reordered), rowNames = TRUE)
  } else{
    openxlsx::writeData(wb, sheet = "Heatmap_matrix", x = reordered, rowNames = TRUE)
  }
  openxlsx::saveWorkbook(wb, file = paste0(output_path, "Heatmap_matrix_", file_suffix, ".xlsx"),
                         overwrite = TRUE)
  
  #****************************************************************************#
  # Plot heatmap with pheatmap
  #****************************************************************************#
  pheatmap::pheatmap(mat                      = as.matrix(reordered),
                     color                    = get(heatmap.params$matrix_color),
                     breaks                   = breaks,
                     border_color             = heatmap.params$border_color,
                     cellwidth                = heatmap.params$bar_width,
                     cellheight               = heatmap.params$bar_height,
                     scale                    = "none",
                     cluster_rows             = FALSE,
                     cluster_cols             = FALSE,,
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method        = "complete",
                     legend                   = heatmap.params$expr_legend, 
                     legend_breaks            = NA,
                     legend_labels            = NA,
                     annotation_row           = row_annotation,
                     annotation_col           = col_annotation,
                     annotation_colors        = ann_colors,
                     annotation_legend        = TRUE,
                     annotation_names_row     = FALSE,
                     annotation_names_col     = FALSE, # set to TRUE if more than 1
                     show_rownames            = dplyr::if_else(length(disp_genes) < 80, TRUE, FALSE, missing = NULL),
                     show_colnames            = dplyr::if_else(length(unique(display_col)) < 50, TRUE, FALSE, missing = NULL),
                     fontsize                 = 5,
                     fontsize_row             = 5,
                     fontsize_col             = 5,
                     gaps_row                 = gaps_row,
                     gaps_col                 = gaps_col,
                     angle_col                = heatmap.params$angle, #"90",
                     fontsize_number          = 0.8*fontsize,
                     labels_row               = display_row,
                     labels_col               = display_col,
                     width                    = heatmap.params$width,
                     height                   = heatmap.params$height,
                     filename                 = paste0(output_path, "Heatmap_", file_suffix, ".", heatmap.params$file_format))
  
  invisible(NULL)
}

#******************************************************************************#
#                       SURVIVAL CURVE RELATED FUNCTIONS                       #
#******************************************************************************#

# Read this paper for survival analyis
# https://doi.org/10.1093/jncimonographs/lgu024

# NOTE:  Output of prep_expr_df is df
#log_norm_counts is matrix with SYMBOLS as rownames
prep_expr_df <- function(log_norm_counts, meta_data, plot_genes, survival_params){
  
  # Merge expression data with survival data
  if (survival_params$gene_sig_score == TRUE){
    
    # Calculate gene signature score
    expr_df <- as.data.frame(advanced_Z(plot_genes, log_norm_counts))
    
    expr_df <- expr_df %>%
      data.frame() %>%
      dplyr::rename(combined.exp = identity(1)) %>%
      tibble::rownames_to_column("Sample.ID") %>%
      dplyr::inner_join(meta_data, by=c("Sample.ID"="Sample.ID")) %>%
      dplyr::select(Sample.ID, combined.exp, Time, Status)
  } else {
    expr_df <- log_norm_counts %>%
      t() %>%
      data.frame() %>%
      tibble::rownames_to_column("Sample.ID") %>%
      dplyr::inner_join(meta_data, by=c("Sample.ID"="Sample.ID")) %>%
      dplyr::select(Sample.ID, all_of(plot_genes), Time, Status)
  }
  
  # Add split_by column to expr_df to define groups in order to calculate multiple_cutoff
  if (!is.na(survival_params$split_by)){
    expr_df <- expr_df %>% 
      dplyr::left_join(meta_data %>% dplyr::select(Sample.ID, survival_params$split_by),
                       by=c("Sample.ID"="Sample.ID"))
  }
  
  return(expr_df)
}

# NOTE:  Output of calc_cutoffs is list(df,ls)
# If plotting by Sex, make sure to create column "model" based on which lines will be split
calc_cutoffs <- function(df, gene, group, survival_params){
  
  # Identify upper & lower cutoffs based on stratify_criteria
  #*************************Split samples by median**************************#
  if(survival_params$stratify_criteria == "m"){
    quartiles <- stats::quantile(x = df[[gene]], 
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    
    cutoff_lower_end <- quartiles[[3]]
    cutoff_upper_end <- quartiles[[3]]
    cutoff_middle <- "NA"
  }
  
  #****************Split samples into top and bottom tertiles****************#
  else if(survival_params$stratify_criteria == "t"){
    tertiles <- stats::quantile(x = df[[gene]],
                                probs = c(0, 0.33, 0.66, 1),
                                na.rm=TRUE)
    
    cutoff_lower_end <- tertiles[[2]]
    cutoff_upper_end <- tertiles[[3]]
    cutoff_middle <- "NA"
  }
  
  #***************Split samples into top and bottom quartiles****************#
  else if(survival_params$stratify_criteria == "q"){
    quartiles <- stats::quantile(x = df[[gene]], 
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    
    cutoff_lower_end <- quartiles[[2]]
    cutoff_upper_end <- quartiles[[4]]
    cutoff_middle <- quartiles[[3]]
  }
  
  #*********************Split expression range by thirds*********************#
  else if(survival_params$stratify_criteria == "th"){
    quartiles <- stats::quantile(x = df[[gene]], 
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    iqr <- stats::IQR(x = df[[gene]],
                      na.rm=TRUE)
    
    # Normal range of expression values lie between cutoff_lower & cutoff_upper
    cutoff_upper <- quartiles[[4]]+1.5*iqr
    cutoff_lower <- dplyr::if_else(quartiles[[1]]-1.5*iqr > 0, quartiles[[1]]-1.5*iqr, 0)
    
    # Based on normal range of expression, identify onethird & twothird cutoff
    cutoff_lower_end <- cutoff_lower + (cutoff_upper-cutoff_lower)/3
    cutoff_upper_end <- cutoff_lower + (cutoff_upper-cutoff_lower)*2/3
    cutoff_middle <- "NA"
  }
  
  #***************Split expression range using optimal cutoff****************#
  else if(survival_params$stratify_criteria == "o"){
    
    # Sometimes quartiles will look like: 
    # 0%       25%      50%      75%     100% 
    # 0.000000 0.000000 0.000000 0.000000 3.495493 
    # In such cases, surv_cutpoint() will fail. So, we add extra if() here.
    quartiles <- stats::quantile(x = df[[gene]], 
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    
    if (quartiles[[4]] > quartiles[[2]]){
      res.cut <- survminer::surv_cutpoint(data = df,
                                          time = "Time",
                                          event = "Status",
                                          variables = gene)
      
      cutoff_lower_end <- res.cut$cutpoint$cutpoint
      cutoff_upper_end <- res.cut$cutpoint$cutpoint
      cutoff_middle <- "NA"
    } else{
      #cat("Surv cutpoint unable to detect optimum cutoff")
      cutoff_lower_end <- "NA"
      cutoff_upper_end <- "NA"
      cutoff_middle <- "NA"
    }
  }
  
  # Categorize the sample based on above cutoffs
  if (survival_params$plot_all_bins == TRUE & survival_params$stratify_criteria == "q"){
    df <- df %>% 
      dplyr::mutate(Expression = dplyr::case_when(get(gene) > cutoff_upper_end ~ "HIGH",
                                                  get(gene) <= cutoff_lower_end ~ "LOW",
                                                  get(gene) <= cutoff_middle ~ "MED_LOW",
                                                  TRUE ~ "MED_HIGH"))
    
  } else if (survival_params$plot_all_bins == TRUE) {
    df <- df %>% 
      dplyr::mutate(Expression = dplyr::case_when(get(gene) > cutoff_upper_end ~ "HIGH",
                                                  get(gene) <= cutoff_lower_end ~ "LOW",
                                                  TRUE ~ "MID"))
    
  } else if (survival_params$stratify_criteria == "none") {
    #When plotting by Sex, Treatment response, we dont use expression data.
    df <- df %>% 
      dplyr::mutate(Expression = model)
    cutoff_lower_end <- NA
    cutoff_upper_end <- NA
    cutoff_middle <- NA
    
  } else {
    df <- df %>% 
      dplyr::mutate(Expression = dplyr::case_when(get(gene) > cutoff_upper_end ~ "HIGH",
                                                  get(gene) <= cutoff_lower_end ~ "LOW",
                                                  TRUE ~ "MID")) %>%
      dplyr::filter(Expression != "MID")
  }
  
  # # Print the cutoffs
  # cat("\nGene         :", gene)
  # cat("\nGroup        :", group)
  # cat("\nLower cutoff :", cutoff_lower_end)
  # cat("\nUpper cutoff :", cutoff_upper_end)
  # cat("\nMiddle cutoff:", cutoff_middle)
  
  # Create a list to store cutoff values
  ls <- list("group" = c(), 
             "gene" = c(), 
             "lower" = c(), 
             "upper" = c(), 
             "middle" = c())
  
  ls$group <- c(group)
  ls$gene <- c(gene)
  ls$lower <- c(cutoff_lower_end)
  ls$upper <- c(cutoff_upper_end)
  ls$middle <- c(cutoff_middle)
  
  # Return the df and the cutoffs
  return(list(df, ls))
}

# NOTE:  Output of calc_surv_stats is list. 
# It also generate survival plot with risk table
calc_surv_stats <- function(df, group, prefix, output_path){
  
  # If all samples belong to one group (like LOW or HIGH or males or female),
  # then quit the function as comparison cannot be done
  if (nrow(df %>% dplyr::count(model)) > 1){
    
    # Create a survival object where Alive = 0, Dead = 1
    surv_object <- survival::Surv(time = df$Time,
                                  event = df$Status,
                                  type = "right",
                                  origin = 0)
    
    # Create a formula for plotting survival curve
    surv_formula <- surv_object ~ model
    
    # Create a fit for survival curve.
    # NOTE: survival::survfit() gives error in ggsurvplot(). Use survminer::surv_fit()
    surv_curve <- survminer::surv_fit(formula = surv_formula,
                                      data = df,
                                      type = "kaplan-meier",
                                      group.by = NULL,
                                      match.fd = FALSE)
    
    # # Check summary of the survival curve with time duration of our interest
    # cat("\nRange of survival (months):", range(df$Time, na.rm=TRUE), "\n")
    # base::summary(surv_curve, times = base::seq(from = floor(range(df$Time, na.rm=TRUE)[[1]]),
    #                                                 to = ceiling(range(df$Time, na.rm=TRUE)[[2]]),
    #                                                 by = 3))
    
    # Create a Cox model for the survival curve and calculate stats
    cox_model <- survival::coxph(formula = surv_formula,
                                 data = df)
    #print(summary(cox_model))
    cat("\n")
    
    # Calculate HR, 95% CI for HR, p-val
    # NOTE: Variable mentioned in summary(cox_model) is numerator in h1(t)/h0(t).
    # The reference variable h0(t) will not be mentioned in co-efficients.
    # Make sure this is not the reference level i.e. low expression. If this is
    # the reference, then we need to reverse the HR ratio, legend labels
    #print(names(cox_model$coefficients))  
    
    # If samples belong to more than 2 groups (like LOW, MID, HIGH), then we 
    # cannot have survival stats. So, we set them to 0.
    if (nrow(df %>% dplyr::count(model)) == 2){
      # Store HR and CI
      if (stringr::str_detect(names(cox_model$coefficients), survival_params$reference)){
        HR <- round(exp(-cox_model$coefficients[[1]]), 2)
        CI <- round(exp(-confint(cox_model)), 2)
        CI_1 <- CI[1]
        CI[1] <- CI[2]
        CI[2] <- CI_1
      } else {
        HR <- round(exp(cox_model$coefficients[[1]]),2)
        CI <- round(exp(confint(cox_model)),2)
      }
      
      # Store pvalues by different methods
      pvals <- c()
      for (test.method in c("survdiff", "1", "n", "sqrtN", "S1","S2", "FH_p=1_q=1")){
        
        # Some of the methods give error. So, we catch them and skip
        if (sum(str_detect(string = class(tryCatch(survminer::surv_pvalue(fit = surv_curve,
                                                                          method = test.method,
                                                                          test.for.trend = FALSE,
                                                                          combine = FALSE), error = function(e) e)),
                           pattern = "error")) == 0){
          p_val <- survminer::surv_pvalue(fit = surv_curve,
                                          method = test.method,
                                          test.for.trend = FALSE,
                                          combine = FALSE)
          pvals <- c(pvals, p_val[[2]])
        } else{
          pvals <- c(pvals, 1)
          print(test.method)
        } 
      } 
    } else {
      HR <- 0
      CI <- c(0, 0)
      pvals <- c(0,0,0,0,0,0,0) # 7 pvalues for each method
    }
    
    # Plot the survival curve using survminer::ggsurvplot() instead of base::plot()
    # ggsurvplot() produces a list of ggplot objects: a survival curve and a risk table
    # Saving it using cowplot() first and then using ggsave() works nicely as
    # compared to saving directly using ggsave()
    if(survival_params$plot_curve == TRUE){
      
      # Plot the survival curve
      legend_label <- df %>% 
        dplyr::count(model) %>% 
        dplyr::select(model) %>% 
        unlist(.,use.names=FALSE)
      
      # We identify proper breaks based on max duration of the dataset
      # We want a maximum of 10 timepoint intervals that are multiples of 12
      max_time <- max(df$Time,na.rm=TRUE)
      n <- floor(max_time/10/12)*12
      if(max_time %/% n <= 10){
        breaks <- n
      } else{
        breaks <- n+12
      }
      
      surv_plot <- survminer::ggsurvplot(fit = surv_curve,
                                         pval = FALSE,
                                         palette = survival_params$color_palette,
                                         linetype = "solid",
                                         size = 1.5,                       # thickness of line
                                         
                                         # Format the legend
                                         legend  = "top",                  # position of legend
                                         legend.title = survival_params$legend_title,
                                         legend.labs = survival_params$legend_label,
                                         
                                         # Format the axes
                                         break.time.by = breaks,           # break X axis in time intervals of 12 months
                                         xlab = "Time (Months)",           # customize X axis label
                                         ylab = "Survival Probability",    # customize Y axis label
                                         title = dplyr::if_else(gene == "combined.exp", "", gene),
                                         
                                         # Format confidence intervals
                                         conf.int = survival_params$conf_interval,
                                         #conf.int.fill = ?,               # color to fill confidence interval
                                         conf.int.style = "ribbon",        # confidence interval style
                                         conf.int.alpha = 0.3,             # confidence fill color transparency
                                         
                                         # Format the risk table
                                         risk.table = survival_params$plot_risk_table,
                                         risk.table.title = "Number at risk",
                                         risk.table.y.text.col = TRUE,     # color of risk table text annotations
                                         risk.table.pos = "out",           # draw risk table outside survival plot
                                         
                                         # Format the censor points
                                         censor = TRUE,
                                         censor.shape = '|',
                                         censor.size = 5)
      
      surv_plot$table <- surv_plot$table + 
        coord_cartesian(x=c(0,ceiling(max_time/breaks)*breaks), clip = "off")
      
      surv_plot$plot <- surv_plot$plot + 
        coord_cartesian(x=c(0,ceiling(max_time/breaks)*breaks), clip = "off")
      
      # Plot p and HR value
      method_plot <- "log-rank"
      p_plot <- pvals[1]  
      
      grob1 <- grobTree(textGrob(label = paste0("p = ", formatC(p_plot, format = "e", digits = 1),
                                                "\nHR = ", round(HR,1), " [", round(CI[1],1), ", ", round(CI[2],1), "]",
                                                "\nMethod = ", method_plot),
                                 x = 0.50,
                                 y = 0.90,
                                 hjust = 0,
                                 gp = grid::gpar(fontfamily="Times", fontface="bold", col="black", fontsize=10)))
      
      # Add p values and HR values to plot
      surv_plot$plot <- surv_plot$plot %++%
        ggplot2::annotation_custom(grob1)
      
      cowplot::plot_grid(plotlist = surv_plot,
                         align = "hv",
                         axis = "tblr",
                         nrow = 2,  
                         ncol = 1, 
                         rel_widths = 1,
                         rel_heights = c(1,0.45),
                         labels = NULL,
                         label_size = 14,
                         label_fontfamily = NULL,
                         label_fontface = "bold",
                         label_colour = NULL)
      
      f_name <- paste0(prefix, "_", group, "_", survival_params$stratify_criteria, ".pdf")
      f_name <- gsub("/", "-", x=f_name)
      
      # Save the plot
      ggplot2::ggsave(filename = f_name,
                      plot = last_plot(),
                      device = "pdf",
                      path = output_path,
                      width = 7,
                      height = 7,
                      units = c("in"),
                      dpi = 300,
                      limitsize = TRUE,
                      bg = NULL)
    }
  } else {
    HR <- 0
    CI <- c(0, 0)
    pvals <- c(0,0,0,0,0,0,0) # 7 pvalues for each method
  }
  
  # Create a list to store survival stats
  ls <- list("group" = c(), 
             "HR" = c(), 
             "CI_lower" = c(), 
             "CI_upper" = c(), 
             "pvalue" =c(), 
             "logrank" = c(), 
             "reg_logrank.late" = c(), 
             "Gehan_Breslow.early" = c(),
             "Tarone_Ware.early" = c(), 
             "Peto_Peto.early" = c(),  
             "modified_Peto_Peto" = c(), 
             "Fleming_Harrington" = c())
  
  ls$group               <- c(group)
  ls$HR                  <- c(HR)
  ls$CI_lower            <- c(CI[1])
  ls$CI_upper            <- c(CI[2])
  ls$logrank             <- c(pvals[1])
  ls$reg_logrank.late    <- c(pvals[2])
  ls$Gehan_Breslow.early <- c(pvals[3])
  ls$Tarone_Ware.early   <- c(pvals[4])
  ls$Peto_Peto.early     <- c(pvals[5])
  ls$modified_Peto_Peto  <- c(pvals[6])
  ls$Fleming_Harrington  <-c(pvals[7])
  
  return(ls)
}

# NOTE: Output of plot_survival is list(df,ls)
plot_survival <- function(expr_df, gene, survival_params, prefix, output_path){
  
  # Create an empty dataframe to store expr_df and classification from calculate_cutoffs()
  survival_data <- data.frame(model = " ")
  
  # Create a list to store results of calculate_cutoffs, surv_plot et
  stats <- list("gene" = c(), 
                "group" = c(),
                "lower_cutoff" = c(),
                "middle_cutoff" = c(),
                "upper_cutoff" = c(),
                "HR" = c(), 
                "CI_lower" = c(), 
                "CI_upper" = c(), 
                "logrank" = c(), 
                "reg_logrank.late" = c(), 
                "Gehan_Breslow.early" = c(),
                "Tarone_Ware.early" = c(), 
                "Peto_Peto.early" = c(),  
                "modified_Peto_Peto" = c(), 
                "Fleming_Harrington" = c())
  
  # Create a list of groups for multiple_cutoff calculation 
  if (survival_params$multiple_cutoff == TRUE & !is.na(survival_params$split_by)){
    cutoff_groups <- expr_df %>% 
      dplyr::add_count(get(survival_params$split_by)) %>%
      dplyr::filter(n>2) %>%
      dplyr::select(all_of(survival_params$split_by)) %>% 
      unlist(use.names=FALSE) %>% 
      unique()
  } else if (survival_params$multiple_cutoff == TRUE & is.na(survival_params$split_by)){
    print("Please define survival_params$split_by to calculate multiple cutoffs")
  } else {
    cutoff_groups <- c(NA)
  }
  
  # STEP 1: Calculate cutoffs
  # If cutoffs need to be calculated for each group, subset the expr_df and pass
  # it to calculate_cutoffs(). Else, pass entire expr_df to calculate_cutoffs()
  for (group in cutoff_groups){
    
    # Subset the expr_df for each group to calculate cutoffs
    if (survival_params$multiple_cutoff == TRUE & !is.na(survival_params$split_by)){
      df <- expr_df %>% dplyr::filter(get(survival_params$split_by) == group)
    } else if (survival_params$multiple_cutoff == TRUE & is.na(survival_params$split_by)){
      cat("\n 'split_by' variable is undefined but 'multiple_cutoff' is set to TRUE")
    } else{
      df <- expr_df
    }
    
    # Calculate cutoffs for each group
    mat <- calc_cutoffs(df, gene, group, survival_params)
    
    ##### Save the data from output of calculate_cutoffs()
    survival_data       <- dplyr::bind_rows(survival_data, mat[[1]])
    stats$gene          <- c(stats$gene,          mat[[2]]$gene)
    #stats$group         <- c(stats$group,         mat[[2]]$group)
    stats$lower_cutoff  <- c(stats$lower_cutoff,  mat[[2]]$lower)
    stats$middle_cutoff <- c(stats$middle_cutoff, mat[[2]]$middle)
    stats$upper_cutoff  <- c(stats$upper_cutoff,  mat[[2]]$upper)
  }
  
  # Populate the model variable by concatenating "Expression" and "split_by"
  if (!is.na(survival_params$split_by)){
    survival_data <- survival_data %>%
      dplyr::mutate(model = paste0(Expression, "_", get(survival_params$split_by))) %>%
      dplyr::filter(!is.na(Sample.ID))
  } else {
    survival_data <- survival_data %>%
      dplyr::mutate(model = Expression) %>%
      dplyr::filter(!is.na(Sample.ID))
  }
  
  # Create a list of groups for plotting survival curves 
  if (!is.na(survival_params$split_by)){
    plot_groups <- expr_df %>% 
      dplyr::add_count(get(survival_params$split_by)) %>%
      dplyr::filter(n>2) %>%
      dplyr::select(all_of(survival_params$split_by)) %>% 
      unlist(use.names=FALSE) %>% 
      unique()
  } else {
    plot_groups <- c(NA)
  }
  
  # STEP 2: Calculate survival stats
  # If each group has to be plotted in separate plots, subset the survival_data
  # and pass it to calc_surv_stats(). Else, pass entire survival_data to 
  # calc_surv_stats().
  for (group in plot_groups){
    
    # Subset the survival_data for each group to generate separate plots
    if (survival_params$split_plot == TRUE & !is.na(survival_params$split_by)){
      df <- survival_data %>% dplyr::filter(get(survival_params$split_by) == group)
    } else if (survival_params$split_plot == TRUE & is.na(survival_params$split_by)){
      cat("\n 'split_by' variable is undefined but 'split_plot' is set to TRUE")
    } else{
      df <- survival_data
    }
    
    # Calculate survival stats for each group
    cox_stats <- calc_surv_stats(df, group, prefix, output_path)
    
    ##### Save the data from output of calc_surv_stats()
    stats$group               <- c(stats$group,               cox_stats$group)
    stats$HR                  <- c(stats$HR,                  cox_stats$HR)
    stats$CI_lower            <- c(stats$CI_lower,            cox_stats$CI_lower)
    stats$CI_upper            <- c(stats$CI_upper,            cox_stats$CI_upper)
    stats$logrank             <- c(stats$logrank,             cox_stats$logrank)
    stats$reg_logrank.late    <- c(stats$reg_logrank.late,    cox_stats$reg_logrank.late)
    stats$Gehan_Breslow.early <- c(stats$Gehan_Breslow.early, cox_stats$Gehan_Breslow.early)
    stats$Tarone_Ware.early   <- c(stats$Tarone_Ware.early,   cox_stats$Tarone_Ware.early)
    stats$Peto_Peto.early     <- c(stats$Peto_Peto.early,     cox_stats$Peto_Peto.early)
    stats$modified_Peto_Peto  <- c(stats$modified_Peto_Peto,  cox_stats$modified_Peto_Peto)
    stats$Fleming_Harrington  <- c(stats$Fleming_Harrington,  cox_stats$Fleming_Harrington)
  }
  return(list(survival_data, stats))
}

# We need to calculate a combined expression value of all genes in the gene
# signature for each sample. Next, we use surv_cutpoint() to classify samples 
# into high and low groups and plot survival curves.
# There are several approaches to calculate the combined expression value. While
# normal z-scaling is logical, https://doi.org/10.1186/gb-2006-7-10-r93 
# recommends a more advanced z-scaling which is implemented below. Refer
# "Measuring gene set expression" in "Materials & stratify_criteria" section.
advanced_Z <- function(gset, eset) {
  
  # Compute z-score for gene set of interest
  ix <- is.element(toupper(rownames(eset)), toupper(gset))
  cat(sum(ix), "\n")
  
  if (sum(ix)>0){
    avg_gset <- base::apply(X=eset[ix,], MARGIN=2, FUN=mean, na.rm=TRUE)
    avg_all <- base::apply(X=eset, MARGIN=2, FUN=mean, na.rm=TRUE)
    sd_all <- base::apply(X=eset, MARGIN=2, FUN=sd, na.rm=TRUE)
    z <- (avg_gset - avg_all)*sqrt(sum(ix))/sd_all
  } else{
    z <- NA
  }
  
  return(z)
}

normal_Z <- function(gset, eset) {
  
  # Compute z-score for gene set of interest
  eset <- eset[gset,]
  a <- t(scale(t(eset)))
  z <- colSums(a, na.rm=TRUE) 
  
  return(z)
}

#******************************************************************************#
#                       CRISPR ANALYSIS RELATED FUNCTIONS                      #
#******************************************************************************#

# Dataframe as input
# Column "Gene" MUST be present
# Column Gene MUST have control sgRNAs labelled as "none" and/or "safe"
calc_t_score <- function(data){
  
  # Create a dataframe of control sgRNAs
  data_control <- data %>%
    dplyr::filter(Gene %in% c("none", "safe"))
  
  median_ctrl <- median(data_control$LFC, na.rm=TRUE)
  sd_ctrl <- sd(data_control$LFC, na.rm=TRUE)
  
  # Normalize to control sgRNAs
  data <- data %>%
    dplyr::mutate(pZ = (LFC-median_ctrl)/sd_ctrl)
  
  data_control <- data_control %>%
    dplyr::mutate(pZ = (LFC-median_ctrl)/sd_ctrl)
  
  U_ctrl <- median(data_control$pZ)
  Var_ctrl <- var(data_control$pZ)
  N_ctrl <- mean((data %>% dplyr::count(Gene))$n)
  # Nctrl is the average number of sgRNAs per gene in a given screen
  
  data <- data %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(U_gene = median(pZ),
                  Var_gene = var(pZ),
                  N_gene = n(),
                  U_ctrl = U_ctrl,
                  Var_ctrl = Var_ctrl,
                  N_ctrl = N_ctrl,
                  S_gene = (Var_gene*(N_gene-1)) + (Var_ctrl*(N_ctrl-1)),
                  t_score = (U_gene - U_ctrl)/sqrt(S_gene/N_gene + S_gene/N_ctrl),
                  Abs_t_score = abs(t_score)) %>%
    dplyr::select(Gene, U_gene, Var_gene, N_gene, U_ctrl, Var_ctrl, N_ctrl, S_gene, t_score, Abs_t_score) %>%
    dplyr::distinct_at("Gene", .keep_all = TRUE)
  
  return(data)
}

# data is a dataframe output of calc_t_score
# save_path is location to save file
plot_t_score <- function(data, save_path, disp_genes, suffix){
  
  y_cutoff <- sort(data$Abs_t_score, decreasing = TRUE)[100]
  xmin <- floor(min(data$U_gene))
  xmax <- ceiling(max(data$U_gene))
  ymin <- 0
  ymax <- max(data$Abs_t_score)
  
  color_breaks <- c(-20,0,20)
  p <- ggplot2::ggplot(data = data,
                       aes(x = U_gene, 
                           y = Abs_t_score,
                           size = Abs_t_score,
                           #color = pz,
                           fill = U_gene)) +
    # Plot dot plot
    ggplot2::geom_point(col="black", 
                        shape=21,
                        stroke=0.5,
                        position=position_jitter(h=0.01,w=0.01)) +
    # Define the theme of plot
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "U_gene") +
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax), clip = "off") +
    #scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
    #scale_y_continuous(breaks = seq(0, 5, by = 1)) +
    ggplot2::guides(size = "none",
                    fill = guide_colourbar(theme = theme(legend.key.width  = unit(0.75, "lines"),
                                                         legend.key.height = unit(10, "lines"),
                                                         legend.ticks = element_blank(),
                                                         legend.frame = element_rect(colour = "Black",
                                                                                     linewidth = 0.5)))) +
    # Define the color of the dots
    ggplot2::scale_fill_viridis_c(option="turbo", limits =c(-5,3))
  
  if (length(disp_genes) > 0){
    p <- p + ggrepel::geom_text_repel(data = data %>% dplyr::filter(Gene %in% disp_genes),
                                      mapping = aes(label = Gene),
                                      size = 5,
                                      show.legend = FALSE,
                                      direction = "both",   #"y"
                                      box.padding = 2.5,      # increases line length somehow
                                      point.padding = 0.1,  # distance around point = dist between line and point
                                      max.overlaps = nrow(data),
                                      position = position_quasirandom())
    
    
  }
  #geom_hline(yintercept= y_cutoff, linetype ="dotted")
  
  # scale_fill_gradientn(colors=c("#007ba7", "Black","#FFFF00"), 
  #                      limits=c(-20, 20), 
  #                      values=c(0, scales::rescale(color_breaks, from = range(color_breaks)), 1))
  #scale_fill_gradient2(low="#007ba7", mid="Black", high="Yellow", midpoint = 0, limits=c(-5, 2))
  #scale_fill_continuous_diverging(palette = "Tofino")
  
  ggplot2::ggsave(filename = paste0(suffix, ".jpg"),
                  plot = p,
                  device = "tiff",
                  path = save_path,
                  width = 7,
                  height = 7,
                  units = c("in"),
                  dpi = 300,
                  limitsize = TRUE,
                  bg = NULL)
}

#******************************************************************************#
#                              CALCULATE padj AND log2FC
#******************************************************************************#

# Function to calculate pval and log2FoldChange
# norm_counts is matrix with log2 transformed values or non-log transformed values
# DO NOT use log10 transformed values
calc_stats <- function(norm_counts, metadata, Target, Reference, log2_transformed_already){
  
  # Perform t.test
  SYMBOL <- c()
  expt <- c()
  control <- c()
  pval <- c()
  for (j in 1:nrow(norm_counts)){
    
    data <- norm_counts[j,] %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Sample") %>%
      dplyr::left_join(metadata, by = c("Sample" = "Sample"))
    colnames(data) <- c("Sample", "values", "Condition")
    
    # Use F.test() to determine if variances are equal or not.
    # NOTE: p < 0.05 => unequal variance
    
    # NOTE: If Iso <- c(NA, NA, 18) and IP <- c(17, NA, 18),
    # var.test and t.test with throw error since Iso has only 1 value.
    
    # NOTE: If Iso <- c(1,1,1) and IP <- c(1,1,1),
    # t.test will throw error "data are essentially constant".
    
    # NOTE: If Iso <- c(0,0,0) and IP <- c(1,1,1),
    # t.test will throw error "data are essentially constant".
    
    if (sum(!is.na(data[data$Condition == Reference, ]$values)) > 1 &
        sum(!is.na(data[data$Condition == Target, ]$values)) > 1 & 
        (length(unique(data[data$Condition == Reference, ]$values)) + 
         length(unique(data[data$Condition == Target, ]$values)) > 2)){
      #   f_test <- var.test(formula = values ~ Condition, 
      #                      data = data,
      #                      alternative = "two.sided")
      #   
      #   if(!is.na(f_test$p.value)){
      # if (f_test$p.value < 0.05){
      # t_test <- t.test(formula = values ~ Condition,
      #                  data = data,
      #                  alternative = "two.sided",
      #                  var.equal = FALSE)
      # }
      
      # # Remove outliers
      # ref_data <- data[data$Condition == Reference, ]$values
      # low_ref <- quantile(ref_data, na.rm=TRUE)[2] - 1.5*IQR(ref_data, na.rm=TRUE)
      # high_ref <- quantile(ref_data, na.rm=TRUE)[3] + 1.5*IQR(ref_data, na.rm=TRUE)
      # 
      # target_data <- data[data$Condition == Target, ]$values
      # low_tar <- quantile(target_data, na.rm=TRUE)[2] - 1.5*IQR(target_data, na.rm=TRUE)
      # high_tar <- quantile(target_data, na.rm=TRUE)[3] + 1.5*IQR(target_data, na.rm=TRUE)
      # 
      # data <- data %>%
      #   dplyr::filter(!(Condition == Reference & (values > high_ref | values < low_ref))) %>%
      #   dplyr::filter(!(Condition == Target & (values > high_tar | values < low_tar)))
      
      
      # Calculate p values, mean expression
      t_test <- stats::t.test(formula = values ~ Condition, 
                              data = data,
                              alternative = "two.sided",
                              var.equal = FALSE)
      
      if (grepl(Reference, names(t_test$estimate[1]))){
        SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[1]])
        expt <- c(expt, t_test$estimate[[2]])
      } else if (grepl(Reference, names(t_test$estimate[2]))) {
        SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[2]])
        expt <- c(expt, t_test$estimate[[1]])
      }
    } else{
      SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
      pval <- c(pval, 1) # Note: DO NOT SET to NA. It will increase padj.
      control <- c(control, mean(data[data$Condition == Reference, ]$values, na.rm=TRUE))
      expt <- c(expt, mean(data[data$Condition == Target, ]$values, na.rm=TRUE))
    }
  }
  
  stats_df <- data.frame(SYMBOL, expt, control, pval)
  stats_df$padj <- stats::p.adjust(p = stats_df$pval, method = "fdr", n = length(stats_df$pval))
  if(log2_transformed_already){
    stats_df$log2FoldChange <- stats_df$expt - stats_df$control  # if data is already log transformed
  }else{
    stats_df$log2FoldChange <- log(stats_df$expt/stats_df$control, base=2)
  }
  
  result <- norm_counts %>% 
    tibble::rownames_to_column(var = "SYMBOL") %>% 
    dplyr::left_join(stats_df, by = c("SYMBOL" = "SYMBOL"))
  
  return(result)
}

flattenCorrMatrix_pmatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  df <- data.frame(row = rownames(cormat)[row(cormat)[ut]],
                   column = rownames(cormat)[col(cormat)[ut]],
                   cor  =(cormat)[ut],
                   p = pmat[ut])
  
  return(df)
}

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  df <- data.frame(row = rownames(cormat)[row(cormat)[ut]],
                   column = rownames(cormat)[col(cormat)[ut]],
                   cor  =(cormat)[ut])
  return(df)
}

# Compare two dataframe and output similarity
compare_deg_results <- function(df1, df2, file_suffix, output_path){
  
  df1 <- df1 %>% 
    dplyr::select(SYMBOL, ENSEMBL_ID, ENSEMBL_SYMBOL, log2FoldChange, padj) %>%
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::mutate(padj = round(padj, 2), log2FoldChange = round(log2FoldChange, 2)) 
  
  df2 <- df2 %>% 
    dplyr::select(SYMBOL, ENSEMBL_ID, ENSEMBL_SYMBOL, log2FoldChange, padj) %>%
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::mutate(padj = round(padj, 2), log2FoldChange = round(log2FoldChange, 2)) 
  
  merged_df <- dplyr::full_join(df1, df2,by=c("ENSEMBL_ID"="ENSEMBL_ID")) %>%
    dplyr::mutate(SYMBOL = dplyr::case_when(!is.na(SYMBOL.x) ~ SYMBOL.x,
                                            !is.na(SYMBOL.y) ~ SYMBOL.y,
                                            TRUE ~ ENSEMBL_ID)) %>%
    dplyr::select(SYMBOL, ENSEMBL_ID, log2FoldChange.x,  log2FoldChange.y, padj.x, padj.y) %>%
    dplyr::mutate(Group = dplyr::case_when(padj.x <= 0.05 & padj.y <= 0.05 & log2FoldChange.x >= 0.58 & log2FoldChange.y >= 0.58 ~ "Up in both",
                                           padj.x <= 0.05 & padj.y <= 0.05 & log2FoldChange.x <= -0.58 & log2FoldChange.y <= -0.58 ~ "Down in both",
                                           padj.x <= 0.05 & log2FoldChange.x >= 0.58 & (padj.y > 0.05 | is.na(padj.y)) ~ "Up x",
                                           padj.x <= 0.05 & log2FoldChange.x <= -0.58 & (padj.y > 0.05 | is.na(padj.y)) ~ "Down x",
                                           (padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y >= 0.58 ~ "Up y",
                                           (padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y <= -0.58 ~ "Down y",
                                           TRUE ~ "Not Significant in both"))
  
  up.x.y <- nrow(merged_df %>% dplyr::filter(log2FoldChange.x >= 0.58, log2FoldChange.y >= 0.58, padj.x <=0.05, padj.y <= 0.05))
  down.x.y <- nrow(merged_df %>% dplyr::filter(log2FoldChange.x <= -0.58, log2FoldChange.y <= -0.58, padj.x <=0.05, padj.y <= 0.05))
  up.x <- nrow(merged_df %>% dplyr::filter(padj.x <= 0.05 & log2FoldChange.x >= 0.58 & (padj.y > 0.05 | is.na(padj.y))))
  down.x <- nrow(merged_df %>% dplyr::filter(padj.x <= 0.05 & log2FoldChange.x <= -0.58 & (padj.y > 0.05 | is.na(padj.y))))
  up.y <- nrow(merged_df %>% dplyr::filter((padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y >= 0.58))
  down.y <- nrow(merged_df %>% dplyr::filter((padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y <= -0.58))               
  
  
  ggplot2::ggplot(data = merged_df, 
                  mapping = aes(x=log2FoldChange.x, y = log2FoldChange.y, color = Group)) +
    geom_point(size=0.75) +
    theme_classic() +
    coord_cartesian(clip = "off") +
    geom_hline(yintercept = c(-0.58, 0.58), color="black", linetype="dashed") +
    geom_vline(xintercept = c(-0.58, 0.58), color="black", linetype="dashed") +
    annotate(geom="text", label=up.x.y, x=10, y=10, col="black", size=5) +
    annotate(geom="text", label=down.x.y , x=-10, y=-10, col="black", size=5) +
    annotate(geom="text", label=up.x, x=10, y=0, col="black", size=5) +
    annotate(geom="text", label=down.x, x=-10, y=0, col="black", size=5) +
    annotate(geom="text", label=up.y, x=0, y=10, col="black", size=5) +
    annotate(geom="text", label=down.y, x=0, y=-10, col="black", size=5)
  
  ggplot2::ggsave(filename = paste0(file_suffix, ".jpg"),
                  plot = last_plot(),
                  device = "tiff",
                  path = output_path,
                  width = 7,
                  height = 7,
                  units = c("in"),
                  dpi = 300,
                  limitsize = TRUE,
                  bg = NULL)
}

#******************************************************************************#
#                         MICROARRAY RELATED FUNCTIONS                         #                       
#******************************************************************************#

# Values should be raw i.e. untransformed. DO NOT use log transformed values etc.
# DO NOT replace NA with 0 etc. Leave NA as they are.
# NA values will be imputed/replaced with average of non-NA values.
# If all values are NA, they will be set to 0.
# NO duplicated genes MUST be present.
impute_with_mean <- function(raw_counts){
  
  # Replace NA with average
  for (j in 1:nrow(raw_counts)){
    
    data <- raw_counts[j,] %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Sample") %>%
      dplyr::left_join(metadata, by = c("Sample" = "Sample"))
    colnames(data) <- c("Sample", "values", "Condition")
    
    # If all values for Reference == NA or Target == NA, set them to 0. 
    # Replace NA with average wherever possible.
    data <- data %>% 
      dplyr::mutate(values = as.numeric(values)) %>%
      dplyr::group_by(Condition) %>% 
      dplyr::mutate(average = mean(values, na.rm=TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(average = dplyr::if_else(is.na(average), 0, average),
                    values = dplyr::if_else(is.na(values), average, values))
    
    data <- data %>% 
      dplyr::select(Sample, values) %>% 
      tibble::column_to_rownames("Sample") %>% 
      t()
    
    if(all(colnames(data) == colnames(raw_counts))){
      raw_counts[j,] <- data[,colnames(raw_counts[j,])] %>% unlist(use.names=FALSE)
    }
  }
  
  imputed_counts <- raw_counts %>% 
    dplyr::mutate(across(.cols=everything(), .fns = as.numeric))
  
  return(imputed_counts)
}

# Perform normalization before imputation so counts can be compared across samples
# NOTE: meta_data MUST have column "Sample" that matches with column names of
# raw_counts 
# NOTE: meta_data MUST have column "Condition" that defines the groups
quantile_norm <- function(raw_counts, meta_data, quant_norm){
  
  # Perform quantile normalization
  if (quant_norm == TRUE){
    
    # https://doi.org/10.1038/s41598-020-72664-6
    # The above paper recommends to quantile normalize within each group rather
    # than whole dataset
    quant_norm_counts <- data.frame(matrix(data=NA, nrow=nrow(raw_counts), ncol=0))
    
    # Subset raw_counts for each group
    for (c in unique(meta_data$Condition)){
      
      samples <- meta_data %>% 
        dplyr::filter(Condition %in% c) %>% 
        dplyr::select(Sample) %>% 
        unlist(use.names=FALSE)
      
      counts <- raw_counts[,samples]
      counts <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(counts)))
      rownames(counts) <- rownames(raw_counts)
      colnames(counts) <- samples
      
      quant_norm_counts <- dplyr::bind_cols(quant_norm_counts, counts)
    }
  } else {
    quant_norm_counts <- raw_counts
  }
  
  return(quant_norm_counts)
}

#******************************************************************************#
#                          NOTES ON PATHWAY ANALYSIS                           #
#******************************************************************************#

# Enrichment analysis takes differential data from every measured gene and looks
# for pathways displaying significantly coordinated shifts in those values.
# https://www.biostars.org/p/12182/
# https://www.biostars.org/p/17628/
# https://www.reddit.com/r/bioinformatics/comments/11o7mrv/gene_set_enrichment_analysis_do_you_separate_out/?rdt=59356
# https://groups.google.com/g/gsea-help/c/oXsBOAUYnH4
# https://www.biostars.org/p/132575/
# https://support.bioconductor.org/p/85681/

#******************************************************************************#
#                      NOTES ON RNA SEQ BATCH CORRECTION                       #
#******************************************************************************#

# NOTE: It is RECOMMENDED to perform batch correction ONLY if you know the batch
# information for all the samples in meta_data.
# https://support.bioconductor.org/p/76099/
# https://support.bioconductor.org/p/9149116/
# https://support.bioconductor.org/p/133222/#133225/
# https://support.bioconductor.org/p/125386/#125387
# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#using-sva-with-deseq2

# There are multiple approaches to batch correction:
# (i) Modelling batch effect using DESeq2 design [RECOMMENDED]
# Add Batch to design of DESEq2 like design=~Batch + id
# (ii) Using combatseq to remove batch effects
# Use combatseq to get batch corrected read counts and then use the corrected 
# counts as input for DESeq2() and a design without batch variable (design=~ id)
# (iii) Using sva to identify hidden batch variables and model them in DeSeq2
# Use sva to find upto 2 surrogate variables that could be causing batch effects
# and include them in the design of DESeq2 like design=~SV1 + SV2 + id

# NOTE: I have tested (i) and (ii) and found that 
# (a) Almost all DEGs identified as significant in (i) are present in (ii) 
# (b) All significant DEGs with -ve log2FC from (i) also have -ve log2FC in (ii)
# (c) All significant DEGs with +ve log2FC from (i) also have +ve log2FC in (ii)
# (d) (ii) identifies many other significant DEGs missing in (i)
# (e) log2FC from (ii) differs from (i) for all genes but difference is minimal
# for most of the significant DEGs

# NOTE: I have tested (i) and (iii) and found that
# (a) Majority of significant DEGs in (iii) match with (i) but there are many
# significant DEGs in (i) absent in (iii) and vice versa
# (b) Significant DEGs with -ve log2FC from (i) also have -ve log2FC in (iii)
# (c) Significant DEGs with +ve log2FC from (i) also have +ve log2FC in (iii) 
# (d) (iii) identifies many other significant DEGs missing in (i) but (iii)
# fails to identify many DEGs identified by (i)
# (e) log2FC from (iii) differs from (i) for MOST genes to a great extent,
# HOWEVER, DDR2KO samples had log2FC(DDR2) of -0.7 in (iii) but only -0.37 in 
# (i) and (ii). So, sva might infact be removing batch effects and detecting 
# true biological effects !!! 

#******************************************************************************#
#                         NOTES ON SINGLE CELL ANALYSIS                        #
#******************************************************************************#

# https://ggplot2.tidyverse.org/reference/guide_colourbar.html
# https://stackoverflow.com/questions/56777529/how-to-pass-bash-variable-into-r-script
# https://github.com/satijalab/seurat/issues/4082
# Major changes in Seurat v5 
# https://satijalab.org/seurat/articles/seurat5_integration
# https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html
# https://satijalab.org/seurat/articles/integration_mapping
# https://satijalab.org/seurat/articles/integration_introduction
# https://satijalab.org/seurat/reference/integratelayers
# https://rdrr.io/r/base/split.html
# https://satijalab.org/seurat/reference/aggregateexpression
# https://github.com/satijalab/seurat/issues/2101

# # Indicate if data is from human or mice. We will adjust gene names accordingly.
# species <- dplyr::case_when(proj %in% c("scRNASeq_Chen",
#                                         "scRNASeq_Simon",
#                                         "visium_GSE171351",
#                                         "scRNASeq_HRA003620",
#                                         "scRNASeq_GSE222315") ~ "Homo sapiens",
#                             TRUE ~ "Mus musculus")

#******************************************************************************#
#        NORMALIZE DATA, IDENTIFY HIGHLY VARIABLE FEATURES, SCALE DATA,        # 
#                PERFORM DIMENSIONAL REDUCTION USING PCA & UMAP                #
#******************************************************************************#

# Use the sctransform method as a more accurate method of normalizing, 
# estimating the variance of the filtered data, and identifying the most 
# variable genes. By default, sctransform accounts for cellular sequencing 
# depth (i.e. nUMIs). Also, we can regress out variation from cell cycle genes
# and mitochondrial genes if needed. 

# Refer https://satijalab.org/seurat/articles/sctransform_vignette.html
# The residuals (normalized values) are stored in pbmc[["SCT"]]@scale.data and 
# used directly as input to PCA. Please note that this matrix is non-sparse, and
# can therefore take up a lot of memory if stored for all genes. To save memory,
# we store these values only for variable genes, by setting the 
# return.only.var.genes=TRUE by default in the SCTransform().

# To assist with visualization and interpretation, we also convert Pearson 
# residuals back to ‘corrected’ UMI counts. You can interpret these as the UMI 
# counts we would expect to observe if all cells were sequenced to the same depth.
# The ‘corrected’ UMI counts are stored in pbmc[["SCT"]]@counts. 

# The log-normalized versions of these corrected counts are stored in 
# pbmc[["SCT"]]@data, which are very helpful for visualization.

# You can use the corrected log-normalized counts for differential expression
# and integration. However, in principle, it would be most optimal to perform
# these calculations directly on the residuals (stored in the scale.data slot) 
# themselves.

#******************************************************************************#
#               PREPARE THE DATA FOR INTEGRATION & INTEGRATE DATA              #
#******************************************************************************#

# As you see from the UMAP, the cells cluster differently in each sample. 
# To find the same cell population (say macrophages) between 2 samples,
# it is necessary for both samples to have similar clustering pattern in UMAP.
# So, we have to integrate the samples. 

# The goal of integration is to ensure that cell types of one 
# condition/dataset align with the same cell types of the other 
# conditions/datasets (e.g. macrophages in control samples align with 
# macrophages in stimulated condition).

# To integrate, we will use the shared highly variable genes from each 
# condition identified using SCTransform, then, we will "integrate" or 
# "harmonize" the conditions to overlay cells that are similar or have a 
# "common set of biological features" between groups. 

# STEP 7A: DECLARE REFERENCE SAMPLES FOR INTEGRATING THE DATA
# STEP 7B: SELECT 3000 MOST VARIABLE GENES TO USE FOR INTEGRATING THE DATA
# STEP 7C: FIND RESIDUALS FOR MISSING GENES
# Each sample has different 3000 most variable genes. Gene X which is most 
# variable among cells of "sample A" may not be one of the top 3000 most 
# variable genes in "sample B". PrepSCTIntegration() will calculate Pearson 
# residuals for missing genes so that all samples have the same 3000 genes

# STEP 7D: FIND COMMON ANCHORS BETWEEN SAMPLES TO INTEGRATE THE DATA
# NOTE: Data must be scaled & PCA must have been run before doing cca or rpca
# in this step. cca is computationally intensive if more than 2 samples are 
# integrated. In such cases, use "rpca". Also, using reference based integration
# is faster.

# (i) Perform canonical correlation analysis (CCA):
# CCA identifies shared sources of variation between the conditions/groups. It
# is a form of PCA, in that it identifies the greatest sources of variation in
# the data, but only if it is shared or conserved across the conditions/groups
# (using the 3000 most variant genes from each sample). This step roughly aligns
# the cells using the greatest shared sources of variation.

# NOTE: The shared highly variable genes are used because they are the most 
# likely to represent those genes distinguishing the different cell types 
# present.

# (ii) Identify anchors or mutual nearest neighbors (MNNs) across datasets 
# (sometimes incorrect anchors are identified): MNNs can be thought of as 
# 'best buddies'. For each cell in one condition:   
# (a) The cell's closest neighbor in the other condition is identified based on
# gene expression values - it's 'best buddy'.
# (b) The reciprocal analysis is performed, and if the two cells are 'best 
# buddies' in both directions, then those cells will be marked as anchors to 
# 'anchor' the two datasets together.

# NOTE: The difference in expression values between cells in an MNN pair 
# provides an estimate of the batch effect, which is made more precise by 
# averaging across many such pairs. A correction vector is obtained and applied
# to the expression values to perform batch correction."
# 
# (iii) Filter anchors to remove incorrect anchors:
# Assess the similarity between anchor pairs by the overlap in their local 
# neighborhoods (incorrect anchors will have low scores)

# STEP 7E: FIND OPTIMUM k.weight FOR USE IN Seurat::IntegrateData()
# k.weight MUST be less than number of anchors. Else, error will be thrown.

# STEP 7F: INTEGRATE THE DATA
# Use anchors and corresponding scores to transform the cell expression values,
# allowing for the integration of the conditions/datasets (different samples, 
# conditions, datasets, modalities)

# NOTE: Transformation of each cell uses a weighted average of the two cells of 
# each anchor across anchors of the datasets. Weights determined by cell 
# similarity score (distance between cell and k nearest anchors) and anchor 
# scores, so cells in the same neighborhood should have similar correction values.

# If cell types are present in one dataset, but not the other, then the cells 
# will still appear as a separate sample-specific cluster.

# STEP 7G: RUN PCA USING 3000 INTEGRATION FEATURES & UMAP USING FIRST 40 PCs
# You need to run PCA and UMAP after integration in order to visualize correctly
# because IntegrateData() uses a different set of 3000 variable genes. So, new
# PCs will need to be calculated.
# Note: If you used SCTransform() before integration, you don't need to run 
# ScaleData() after integration. However, if you ONLY used NormalizeData() 
# before integration, you need to use ScaleData() after integration.

#************************STEP 7B: INTEGRATE THE DATA*************************#

# NOTE: The work of SelectIntegrationFeatures(), PrepSCTIntegration(), 
# FindIntegrationAnchors() and IntegrateData() are done by IntegrateLayers().
# Additionally, a new reduction which is equivalent of RunPCA() is also 
# created after integration.

# NOTE: RPCA needs proper kweight. Else, it throws error. I have not yet found
# a way to calculate optimal kweight unlike seurat v3. If script gives error
# regarding kweight, use the kweight it recommends in the error and re-run.

#******************************************************************************#
#                 CLUSTER THE CELLS & REMOVE SCARCE CLUSTERS                   #
#******************************************************************************#

# FindNeighbors() uses the user indicated "reduction" to calculate the k-nearest
# neighbors and construct the SNN graph.
# FindClusters() then performs graph-based clustering on the SNN graph. 

# NOTE: It is recommended to adjust k.param of FindNeighbors() [default=20] to 
# the same value as n.neighbors of UMAP() [default=30] 
# https://github.com/satijalab/seurat/issues/2152

#**************************STEP 8C: MERGE ALL LAYERS*************************#

# Once integrative analysis is complete, you can rejoin the layers - which 
# collapses the individual datasets together and recreates the original 
# counts and data layers. You will need to do this before performing any 
# differential expression analysis. However, you can always resplit the 
# layers in case you would like to reperform integrative analysis.

#******************************************************************************#


# ### Visualize the number of UMIs per cell
# umi_qc <- function(meta){
#   
#   umi_cutoff <- 500
#   # metadata <- meta %>%
#   #   dplyr::mutate(QC = factor(QC, levels = qc_levels))
#   
#   # p <- ggplot(data = metadata, aes(x = Sample, y = nUMIs, fill = QC)) +
#     # geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) + 
#     # geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) +
#     # theme_classic() +  
#     # my_theme +
#     #labs(x = "Sample", y = "Number of UMIs", title = "Distribution of UMIs") +
#     #coord_cartesian(ylim = c(1,1000000), clip = "off") +
#     #scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) +
#     # scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
#     #                              "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
#     #geom_hline(yintercept = umi_cutoff, linetype = 2)
#   
#   return(p)
# }

# ### Visualize the number of genes per cell
# gene_qc <- function(meta){
#   
#   gene_cutoff <- 250
#   # metadata <- meta %>%
#   #   dplyr::mutate(QC = factor(QC, levels = c("Empty Droplet","Doublet", "Low Quality", "Singlet")))
#   
#   # p <- ggplot(data = metadata, aes(x = Sample, y = nGenes, fill = QC)) +
#   #   geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) + 
#   #   geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) +
#   #   theme_classic() + 
#   #   my_theme + 
#   #  labs(x = "Sample", y = "Number of Genes", title = "Distribution of Genes") +
#     coord_cartesian(ylim = c(1, 30000), clip = "off") +
#     scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000)) + 
#     # scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
#     #                              "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
#     #geom_hline(yintercept = gene_cutoff, linetype = 2)
#   
#   return(p)
# }

# # Visualize the MitoRatio of each cell
# mito_qc <- function(meta){
#   
#   mito_cutoff <- 0.2
#   # metadata <- meta %>%
#   #   dplyr::mutate(QC = factor(QC, levels = c("Empty Droplet","Doublet", "Low Quality", "Singlet")))
#   
#   # p <- ggplot(data = metadata, aes(x = Sample, y = MitoRatio, fill = QC)) +
#   #   geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) + 
#   #   geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) + 
#   #   theme_classic() + 
#   #   my_theme +
#     #labs(x = "Sample", y = "MitoRatio", title = "Distribution of MitoRatio") +
#     coord_cartesian(ylim = c(0.00001, 1), clip = "off") +
#     scale_y_log10(breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1)) + 
#     # scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
#     #                              "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
#     #geom_hline(yintercept = mito_cutoff, linetype = 2)
#   
#   return(p)
# }

# ### Visualize the RiboRatio of each cell
# ribo_qc <- function(meta){
#   
#   ribo_cutoff <- 0.05
#   # metadata <- meta %>%
#   #   dplyr::mutate(QC = factor(QC, levels = c("Empty Droplet","Doublet", "Low Quality", "Singlet")))
#   
#   # p <- ggplot(data = metadata, aes(x = Sample, y = RiboRatio, fill = QC)) +
#   #   geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) + 
#   #   geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) + 
#   #   theme_classic() + 
#   #   my_theme +
#    # labs(x = "Sample", y = "RiboRatio", title = "Distribution of RiboRatio") +
#     coord_cartesian(ylim = c(0.0001, 1), clip = "off") +
#     scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) + 
#     # scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
#     #                              "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
#     #geom_hline(yintercept = ribo_cutoff, linetype = 2)
#   
#   return(p)
# }

# ### Visualize the novelty or complexity of each cell
# novelty_qc <- function(meta){
#   
#   novelty_cutoff <- 0.8
#   # metadata <- meta %>%
#   #   dplyr::mutate(QC = factor(QC, levels = c("Empty Droplet","Doublet", "Low Quality", "Singlet")))
#   
#   # p <- ggplot(data = metadata, aes(x = Sample, y = Novelty, fill = QC)) +
#   #   geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) + 
#   #   geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) + 
#   #   theme_classic() + 
#   #   my_theme +
#     # labs(x = "Sample", y = "Novelty", title = "Distribution of Novelty Score") +
#     coord_cartesian(ylim = c(0.3, 1), clip = "off") +
#     scale_y_log10(breaks = c(0.3, 1)) + 
#     # scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
#     #                              "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
#     #geom_hline(yintercept = novelty_cutoff, linetype = 2)
#   
#   return(p)
# }

# # Plot all QC metrics before and after QC
# funcs <- c("cell_qc", "umi_qc", "gene_qc", "mito_qc", "ribo_qc", "novelty_qc",
#            "gene_umi_mito_qc")
# 
# filenames <- c("Cell_Counts", "UMI_Distribution", "Gene_Distribution",
#                "MitoRatio_Distribution", "RiboRatio_Distribution", 
#                "Novelty_Score_Distribution", "Genes_UMI_MitoRatio_Distribution")
# 
#for (i in 1:length(funcs)){

# # Plot QC metrics
# purrr::map(.x = c("raw_metadata"), .f = get(funcs[i])) %>% 
#   cowplot::plot_grid(plotlist = ., align = "hv", axis = "tblr", nrow = 1, ncol = 1)

#p <- get(funcs[i])(raw_metadata)

#   # Save the plot
#   ggplot2::ggsave(filename = paste0("QC_", filenames[i], ".pdf"),
#                   plot = p, #last_plot(),
#                   device = "pdf",
#                   path = output_path,
#                   #scale = 1,
#                   width = 11,
#                   height = 8,
#                   units = c("in"),	 
#                   dpi = 600,
#                   limitsize = TRUE,
#                   bg = NULL)
# }

#cat("QC plots generated\n")  
#}