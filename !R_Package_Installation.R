#!/usr/bin/env Rscript

.libPaths(new = "/hpc/home/kailasamms/R/x86_64-pc-linux-gnu-library/4.1")
.libPaths()
chooseCRANmirror(ind=1)

# NOTE: Install all packages from submit node NOT computing node. 

# NOTE: Some packages need hdf5r package which in turn needs HDF5 library files
# for installation. Others like SSGSEA need gcc compiler. Rcpp which is needed 
# by most R packages doesnt get installed if you skip "module load hdf5/1.8.18"
# So, before installing do the following:
# conda activate R
# module load hdf5/1.8.18
# module load gcc/11.1.0
# module load zlib/1.3
# module load libtiff/4.2.0
# module load cmake/3.19.5        #3.2.3 doesnt work
# R

# NOTE: If you get "libxml not found", "Cannot find xml2-config",
# "ERROR: configuration failed for package ‘XML’", then, 
# Exit R, then "conda activate R", then "conda install -c conda-forge r-xml"
# IMPORTANT: DO NOT UPDATE XML even if outdated. 

#******************************************************************************#
#                          INSTALL NECESSARY PACKAGES                          #
#******************************************************************************#

### Non-CRAN repository managing packages
install.packages(pkgs = c("BiocManager", "remotes"),
                 repos ='http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

# ------------------------------
# Data Analysis Packages
# ------------------------------

## Gene Annotation
library(AnnotationHub)        # Access annotation resources
library(ensembldb)            # Gene model annotation
library(org.Hs.eg.db)         # Human gene database
library(org.Mm.eg.db)         # Mouse gene database

## Microarray Analysis
library(affy)                 # Affymetrix array processing
library(illuminaHumanv4.db)   # Illumina array annotation
library(lumi)                 # Illumina expression analysis
library(limma)                # Linear modeling for arrays

## RNA-seq & Differential Expression
library(DESeq2)               # RNA-seq differential expression
library(sva)                  # Batch effect correction

## Gene Set Enrichment & Pathway Analysis
library(fgsea)                # Fast GSEA
library(clusterProfiler)      # Enrichment analysis tools
library(enrichplot)           # Enrichment result visualization
library(progeny)              # Pathway activity estimation
library(dorothea)             # TF-target interaction network
library(viper)                # TF activity inference

# ------------------------------
# Data Wrangling Packages
# ------------------------------

## Core Data Manipulation
library(dplyr)                # Data manipulation grammar
library(tibble)               # Tidy data frames
library(purrr)                # Functional programming tools
library(stringr)              # String manipulation

## Excel File Handling
library(openxlsx)             # Read/write Excel files

# ------------------------------
# Basic Plotting Packages
# ------------------------------

## Core Plotting Tools
library(ggplot2)              # Grammar of graphics
library(ggpubr)               # Publication-ready plots
library(cowplot)              # Combine multiple plots

## Labeling & Positioning
library(ggrepel)              # Smart label placement
library(ggbeeswarm)           # Bee swarm plots

## Color Palettes & Aesthetics
library(RColorBrewer)         # Color palettes
library(viridis)              # Colorblind-friendly scales
library(colorspace)           # Color manipulation tools

# ------------------------------
# Specialized Plotting Packages
# ------------------------------

## Heatmaps & Ridge Plots
library(pheatmap)             # Pretty heatmaps
library(ggridges)             # Ridge (joy) plots

## Dimensional & 3D Plots
library(plot3D)               # 3D plotting

## Set Visualization
library(VennDiagram)          # Venn diagrams
library(UpSetR)               # UpSet intersection plots

## Survival Analysis
library(survival)             # Survival modeling
library(survminer)            # Survival plot enhancements

## Single-Cell Customization
library(scCustomize)          # Customize Seurat plots

# ------------------------------
# Single-Cell Analysis Packages
# ------------------------------

## Core Seurat Ecosystem
library(Seurat)               # Single-cell analysis toolkit
library(SeuratData)           # Load Seurat datasets
library(SeuratWrappers)       # Extra Seurat functions
library(patchwork)            # Arrange ggplots easily

## Dimensionality Reduction & Integration
library(harmony)              # Batch correction integration
library(umap)                 # UMAP dimension reduction

## Cluster & Cell Type Analysis
library(Banksy)               # Label transfer and mapping
library(UCell)                # Signature scoring method
library(ClusterFoldSimilarity)# Cluster similarity across datasets

## Quality Control
library(DropletUtils)         # Empty droplet detection
library(DoubletFinder)        # Detect doublets in Seurat
library(scDblFinder)          # Another doublet detection method

## Cell-Cell Communication
library(CellChat)             # Infer cell-cell signaling

## Utilities
library(xpectr)               # Suppress warnings and messages

# ------------------------------
# Commented/Optional Libraries
# ------------------------------

# library(ChIPQC)           # ChIP-seq quality control
# library(pathview)         # Pathway visualization tool
# library(hrbrthemes)       # Modern ggplot themes
# library(infercnv)         # CNV from single-cell RNA
# library(SeuratDisk)       # HDF5 (.h5ad) conversion
# library(SCopeLoomR)       # Read .loom format files
# library(SCENIC)           # SC regulon analysis
# library(Augur)            # Cell-type prioritization
# library(ktplots)          # CellPhoneDB visualizations






### Annotation packages
BiocManager::install(pkgs = c("AnnotationHub", "ensembldb", 
                              "org.Hs.eg.db", "org.Mm.eg.db"),
                     force = FALSE,
                     INSTALL_opts = '--no-lock')

### Functional data analysis packages
BiocManager::install(pkgs = c("fgsea", "clusterProfiler", "decoupleR", 
                              "progeny", "dorothea", 
                              "viper", "DESeq2", "sva", "GSVA", "glmGamPoi", 
                              "ashr"),
                     force = FALSE,
                     INSTALL_opts = '--no-lock')

### Single cell data analysis packages
Sys.setenv("ARROW_R_DEV"=TRUE, "LIBARROW_BINARY"=FALSE,
           "ARROW_WITH_ZSTD"="ON", "ARROW_DEPENDENCY_SOURCE"="BUNDLED")
install.packages(pkgs = c("Seurat", "harmony", "hdf5r", "arrow", "leidenAlg",
                          "scCustomize", "reticulate"),
                 repos ='http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

BiocManager::install(pkgs = c("infercnv", "UCell", "scDblFinder", "RcisTarget",
                              "DropletUtils", "batchelor", "ClusterFoldSimilarity"),
                     force = FALSE,
                     INSTALL_opts = '--no-lock')

# Sometimes connection is timed out when connecting to GitHub. So, use a proxy.
library("httr")
httr::set_config(use_proxy("8.8.8.8", port = 8080))
remotes::install_github(repo = c("jinworks/CellChat", 
                                 "prabhakarlab/Banksy@devel", 
                                 "satijalab/seurat-wrappers", # needed for Banksy
                                 "immunogenomics/presto",     # needed for faster FindMarkers()
                                 "chris-mcginnis-ucsf/DoubletFinder",
                                 "satijalab/seurat-data"),
                        force = FALSE,
                        INSTALL_opts = '--no-lock')

### Microarray data analysis packages
BiocManager::install(pkgs = c("oligo", "oligoData", "illuminaHumanv4.db",
                              "hgu133plus2.db", "GEOquery", "affy", "lumi"),
                     force = FALSE,
                     INSTALL_opts = '--no-lock')

### Data wrangling packages
install.packages(pkgs = c("openxlsx", "dplyr", "tibble", "stringr", "purrr"),
                 repos ='http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

### Basic and Specialized Graph plotting packages
install.packages(pkgs = c("ggplot2", "ggplotify", "ggrepel", "ggpubr", 
                          "ggfortify", "ggridges", "ggbeeswarm", "pheatmap", 
                          "VennDiagram", "survival", "survminer", "UpSetR", 
                          "umap", "plot3D", "cowplot", "viridis", 
                          "RColorBrewer", "colorspace"),
                 repos = 'http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

BiocManager::install(pkgs = c("enrichplot", "ComplexHeatmap"),
                     force = FALSE,
                     INSTALL_opts = '--no-lock')

# List outdated packages
utils::old.packages(lib.loc = .libPaths(),
             repos = 'http://cran.us.r-project.org')

# Update outdated packages
utils::update.packages(lib.loc = .libPaths(),
                    repos = 'http://cran.us.r-project.org')

# List packages that couldn't be installed
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
         
# Display packages that couldn't be installed
cat("These pacakges have NOT yet been installed:", sort(pkgs[!(pkgs %in% installed.packages()[,1])]), sep="\n")


#BiocManager::install(pkgs = c("BiocNeighbors", ),
# "mojaveazure/seurat-disk",
# "aertslab/SCENIC",
# "aertslab/SCopeLoomR",
# "neurorestore/Augur"
# "ChIPQC",  "metap", "SeuratDisk", "DoubletFinder", "Augur", ,
# "SCENIC", "SCopeLoomR", "ktplots", , "wordcloud",, 
# "pathview", "reticulate", "ashr", "multtest",  
