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
# module load cmake/3.19.5    #3.2.3 doesnt work  
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

### Annotation packages
BiocManager::install(pkgs = c("AnnotationHub", "ensembldb", 
                              "org.Hs.eg.db", "org.Mm.eg.db"),
                     force = FALSE,
                     INSTALL_opts = '--no-lock')

### Functional data analysis packages
BiocManager::install(pkgs = c("fgsea", "clusterProfiler", "progeny", "dorothea", 
                              "viper", "DESeq2", "sva", "GSVA"),
                     force = FALSE,
                     INSTALL_opts = '--no-lock')

### Single cell data analysis packages
Sys.setenv("ARROW_R_DEV"=TRUE, "LIBARROW_BINARY"=FALSE,
           "ARROW_WITH_ZSTD"="ON", "ARROW_DEPENDENCY_SOURCE"="BUNDLED")
install.packages(pkgs = c("Seurat", "harmony", "hdf5r", "arrow"),
                 repos ='http://cran.us.r-project.org',
                 force = FALSE,
                 INSTALL_opts = '--no-lock')

BiocManager::install(pkgs = c("infercnv"),
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

BiocManager::install(pkgs = c("enrichplot"),
                     force = FALSE,
                     INSTALL_opts = '--no-lock')

# Sometimes connection is timed out when connecting to GitHub. So, use a proxy.
library("httr")
httr::set_config(use_proxy("8.8.8.8", port = 8080))
remotes::install_github(repo = c("nicolerg/ssGSEA2"),
                        
                        INSTALL_opts = '--no-lock')

# List outdated packages
utils::old.packages(lib.loc = .libPaths(),
             repos = 'http://cran.us.r-project.org')

# Update outdated packages
utils::update.packages(lib.loc = .libPaths(),
                    repos = 'http://cran.us.r-project.org')

# List packages that couldn't be installed
pkgs <- c("BiocManager", "remotes", "AnnotationHub", "ensembldb", 
          "org.Hs.eg.db", "org.Mm.eg.db", "fgsea", "clusterProfiler", "progeny",
          "dorothea", "viper", "DESeq2", "sva", "GSVA", "Seurat", "harmony", 
          "hdf5r", "arrow", "infercnv", "oligo", "oligoData", "illuminaHumanv4.db", 
          "hgu133plus2.db", "GEOquery", "affy", "lumi", "openxlsx", "dplyr", 
          "tibble", "stringr", "purrr", "ggplot2", "ggplotify", "ggrepel", "ggpubr", 
          "ggfortify", "ggridges", "ggbeeswarm", "pheatmap", "VennDiagram", 
          "survival", "survminer", "UpSetR", "umap", "plot3D", "cowplot", 
          "viridis", "RColorBrewer", "colorspace", "enrichplot", "ssGSEA2")
         
# Display packages that couldn't be installed
cat("These pacakges have NOT yet been installed:", sort(pkgs[!(pkgs %in% installed.packages()[,1])]), sep="\n")

# # arrow package is key dependency for SCENIC
# # Install using compute node. It might fail using submit node.
# #BiocNeighbors is key dependency for CellChat


# # "BiocNeighbors" is key dependency for CellChat package
# # "RcisTarget" is key dependency for SCENIC package
# BiocManager::install(pkgs = c("BiocNeighbors", "RcisTarget", "UCell", 
#                               "glmGamPoi", "scCustomize"),
#                      force = FALSE,
#                      INSTALL_opts = '--no-lock')

# "mojaveazure/seurat-disk",
# "aertslab/SCENIC",
# "aertslab/SCopeLoomR",
# "sqjin/CellChat",
# "stephens999/ashr",
# "chris-mcginnis-ucsf/DoubletFinder",
# "neurorestore/Augur"
# "ChIPQC",  "metap", "SeuratDisk", "DoubletFinder", "Augur", "glmGamPoi",
# "SCENIC", "SCopeLoomR", "ktplots", "CellChat", "wordcloud", "UCell", 
# "pathview", "reticulate", "ashr", "multtest", "RcisTarget", "BiocNeighbors", 
