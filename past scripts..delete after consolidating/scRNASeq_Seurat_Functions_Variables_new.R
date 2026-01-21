#!/usr/bin/env Rscript

#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

# Data analysis packages
#library("TCGAbiolinks")         # Needed for TCGA data analysis
library("ensembldb")            # Needed for annotating genes
library("AnnotationHub")        # Needed for annotating genes
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("affy")                 # Needed for micro-array analysis
library("lumi")                 # Needed for micro-array analysis
library("illuminaHumanv4.db")  # Needed for micro-array analysis
library("limma")                # Needed for micro-array analysis
#library("ChIPQC")               # Needed for ChIP analysis
library("fgsea")                # Needed for GSEA analysis
library("enrichplot")           # Needed for GSEA analysis
library("clusterProfiler")      # Needed for GSEA analysis
#library("pathview")            # Needed for GSEA analysis
library("DESeq2")               # Needed for Differential Expression analysis
library("progeny")
library("dorothea")
library("viper")
#library("infercnv")
library("sva")

# Data wrangling packages
library("openxlsx")             # Needed for reading, writing xlsx files
library("dplyr")                # Needed for data wrangling
library("tibble")               # Needed for data wrangling
library("stringr")              # Needed for data wrangling
library("purrr")                # Needed for data wrangling

# Graph plotting packages
library("ggplot2")              # Needed for making graphs
library("cowplot")              # Needed for merging multiple graphs
library("viridis")              # Needed for nice graph coloring
library("RColorBrewer")         # Needed for nice graph coloring
library("ggrepel")              # Needed for making graphs prettier
library("ggpubr")              # Needed for adding p values to graphs
library("ggbeeswarm")           # Needed for proper positioning of labels in scatter plots
library("colorspace")

# Specialized Graph plotting packages
library("pheatmap")             # Needed for making heatmap
library("ggridges")             # Needed for making ridgeplots
#library("hrbrthemes")           # Needed for modern look of plots
library("VennDiagram")          # Needed for making Venn diagram
library("survival")             # Needed for making survival curves
library("survminer")            # Needed for making survival curves, to handle %++% in survival function
library("scCustomize")          # Needed fro customizing Seurat plots

# Single cell analysis packages
library("Seurat")               # Needed for single cell analysis
library("SeuratData")
library("SeuratWrappers")
library("Banksy")
library("UCell")
#library("SeuratDisk")           # Needed for reading h5ad files
#library("SCopeLoomR")           # Needed for reading loom files
library("harmony")              # Needed for single cell analysis
#library("SCENIC")              # Needed for SCENIC analysis
library("DropletUtils")         # Needed for identifying empty droplets
library("DoubletFinder")        # Needed for identifying doublets
library("scDblFinder")          # Needed for identifying doublets
#library("Augur")
#library("ktplots")             # Needed for plotting cellphonedb results
library("CellChat")
library("patchwork")

#
library("xpectr")             # Suppress warnings , messages

#******************************************************************************#
#                       DEFINE GLOBAL OPTIONS & VARIABLES                      #
#******************************************************************************#

# Change default limit for allowable object sizes within R 
options(future.globals.maxSize=1e15)
options(Seurat.object.assay.version = "v5")
options(scipen=999)                         # disables scientific notation (e.g., 1e+05)

# NOTE: proj variable is henceforth defined within scRNASeq wrapper and will be
# read by the Rscript calling this R file.

# Choose xlsx file with metadata for the project. 
# NOTE: This xlsx file MUST be named in the format: <proj>_Metadata.xlsx.
# NOTE: This xslx file should have column named "Unique_ID" whose values matches 
# with  column "Unique_ID" of seurat object's metadata.
metafile <- paste0(proj, "_Metadata.xlsx")

# Define directory paths
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

# Create a list of S and G2M markers
cell_cycle_markers <- openxlsx::read.xlsx(xlsxFile = paste0(scripts_path, "Cell_Cycle_Markers.xlsx"))

cell_cycle_genes <- c(cell_cycle_markers$Human_Gene, 
                      cell_cycle_markers$Mouse_Gene)
s_genes <- c(cell_cycle_markers$Human_Gene[which(cell_cycle_markers$Phase=="S")], 
             cell_cycle_markers$Mouse_Gene[which(cell_cycle_markers$Phase=="S")])
g2m_genes <- c(cell_cycle_markers$Human_Gene[which(cell_cycle_markers$Phase=="G2/M")],
               cell_cycle_markers$Mouse_Gene[which(cell_cycle_markers$Phase=="G2/M")])

# Define axis font etc to use in all plots
my_theme <- ggplot2::theme(plot.title=  element_text(family="sans", face="bold",  colour="black", size=15, hjust=0.5),
                           legend.title=element_text(family="sans", face="bold",  colour="black", size=12, hjust=0,   vjust=1,   angle=0),
                           axis.title.x=element_text(family="sans", face="bold",  colour="black", size=12, hjust=0.5, vjust=0,   angle=0),
                           axis.title.y=element_text(family="sans", face="bold",  colour="black", size=12, hjust=0.5, vjust=1,   angle=90),
                           legend.text= element_text(family="sans", face="plain", colour="black", size=10, hjust=0.5),
                           axis.text.x= element_text(family="sans", face="plain", colour="black", size=10, hjust=0.5, vjust=0.5, angle=45),
                           axis.text.y= element_text(family="sans", face="plain", colour="black", size=10, hjust=0.5, vjust=0.5, angle=0))
# strip.text.x=element_text(family="sans", face="bold",  colour="black", size=10, hjust=0.5),
# legend.background=element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
# legend.position="right",
# legend.justification="left",
# legend.direction="vertical",
# legend.key.height=unit(0.5, 'cm'),
# legend.key.width =unit(0.5, 'cm'), 
# legend.text.align=0)

# Assign colors for UMAP. The current palette supports up to 33 cell types
my_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#8C564B", 
                "#E377C2", "#BCBD22", "#17BECF", "#FFC61E", "#762A83",
                "#333333", "#FF1F5B", "#B8E80C", "#9b19f5", "#DC0AB4")

my_palette <- c(my_palette, 
                colorspace::adjust_transparency(col=my_palette, alpha=0.2), 
                colorspace::adjust_transparency(col=my_palette, alpha=0.4),
                colorspace::adjust_transparency(col=my_palette, alpha=0.6),
                colorspace::adjust_transparency(col=my_palette, alpha=0.8))

# my_palette <- c("#000000","#D9D9D9","#003C30","#beb9db","#1D91C0","#A6CEE3",
#                 "#50E991","#A6D854","#74C476","#C7E9B4","#00bfa0","#E5F5F9",
#                 "#EDBF33","#E6D800","#FFF7BC","#ffee65","#C7EAE5","#67001F",
#                 "#CB181D","#FD8D3C","#FC9272","#EF3B2C","#F16913","#6A51A3",
#                 "#762A83","#D4B9DA","#0bb4ff","#E60049","#AE017E","#DF65B0",
#                 "#FDCCE5","#AEC7E8","#FFBB78","#98DF8A","#FF9896","#C5B0D5",
#                 "#C49C94","#F7B6D2","#C7C7C7","#DBDB8D","#9EDAE5","#C51B7D",
#                 "#DE77AE","#7F7F7F","#9467BD")

#******************************************************************************#
#                    SINGLE CELL ANALYSIS RELATED FUNCTIONS                    #       
#******************************************************************************#

### Import data from output of cellranger
# Input is sample name and path to feature barcode matrix
# Ouput is seurat object
read_cellranger <- function(sample, path){
  
  # Read the feature-barcode matrices into dgCMatrix object
  # NOTE: gene.column=1 imports Ensembl ids while gene.column=2 imports gene
  # symbols from features.tsv. We need gene symbols as we will calculate 
  # mitoratio, riboratio, hemeratio using gene names
  dgCMatrix <- Seurat::Read10X(data.dir = paste0(path, sample),
                               gene.column = 2,  
                               cell.column = 1,
                               unique.features = TRUE,
                               strip.suffix = FALSE)
  
  # Create a seurat object for each dgCMatrix object
  # Since EmptyDrops() is run on raw matrix, set min.cells=0 & min.features=0
  sample.seurat <- SeuratObject::CreateSeuratObject(counts = dgCMatrix,
                                                    project = sample,
                                                    assay = "RNA",
                                                    names.field = 1,
                                                    names.delim = "_",
                                                    meta.data = NULL,
                                                    min.cells = 0,
                                                    min.features = 0)
  if (grepl("raw", path)){
    cat(paste0("Raw Feature Barcode Matrix imported for '", sample, "'\n"))
  } else if (grepl("filt", path)){
    cat(paste0("Filtered Feature Barcode Matrix imported for '", sample, "'\n"))
  } else {
    cat("Raw or Filtered Barcode Matrix imported for '", sample, "'\n")
  }
  
  return(sample.seurat)
}

### Identify empty droplets using DropletUtils()
# Input is seurat object of a single sample after read_cellranger()  
# Output is seurat object with DropletUtils column (having Singlet/Empty Droplet) added to metadata
mark_emptydroplets_dropletutils <- function(sample.seurat){ 
  
  sce.sample.seurat <- Seurat::as.SingleCellExperiment(x = sample.seurat)
  
  # If FDR > 0.05 for some droplets AND Limited == TRUE, it indicates that with
  # more iterations, the FDR of these droplets can be reduced.
  set.seed(100)      # to obtain reproducible results
  n_improve <- 1
  niters <- 10000
  
  while (n_improve > 0){
    e.out <- DropletUtils::emptyDrops(m = SingleCellExperiment::counts(sce.sample.seurat),
                                      niters = niters)
    n_improve <- nrow(e.out %>% 
                        data.frame() %>% 
                        dplyr::filter(Limited == TRUE, FDR > 0.05))
    cat("n_improve:", n_improve, "\tniters:", niters, "\n")
    niters <- niters + 10000
  }
  
  true_cells <- e.out %>% 
    data.frame() %>% 
    dplyr::filter(FDR <= 0.05) %>% 
    rownames()
  
  # Mark cells as empty droplets
  sample.seurat@meta.data <- sample.seurat@meta.data %>%
    dplyr::mutate(Cell = rownames(.)) %>%
    dplyr::mutate(DropletUtils = dplyr::case_when(Cell %in% true_cells ~ "Non-Empty Droplet",
                                                  TRUE ~ "Empty Droplet"))
  
  cat(paste0("DropletUtils empty droplets identified for '", as.character(unique(sample.seurat@meta.data$orig.ident)), "'\n"))
  return(sample.seurat)
}

### Mark empty droplets identified using CellRanger
# Input is seurat object of a single sample after mark_emptydroplets_dropletutils() 
# Output is seurat object with CellRanger column (having Singlet/Empty Droplet) added to metadata
mark_emptydroplets_cellranger <- function(sample.seurat){
  
  # Read the filtered barcode-feature matrix output of cellranger
  sample <- s.obj@meta.data$orig.ident %>% unique() %>% as.character()
  sample.seurat.filt <- read_cellranger(sample, filt_matrix_path)
  
  # Mark cells absent in filtered barcode-feature matrix as empty droplets
  sample.seurat@meta.data <- sample.seurat@meta.data %>%
    dplyr::mutate(Cell = rownames(.)) %>%
    dplyr::mutate(CellRanger = dplyr::case_when(Cell %in% colnames(sample.seurat.filt) ~ "Non-Empty Droplet",
                                                TRUE ~ "Empty Droplet"))
  
  cat(paste0("CellRanger empty droplets identified for '", as.character(unique(sample.seurat@meta.data$orig.ident)), "'\n"))
  return(sample.seurat)
}

### Identify doublets using DoubletFinder()
# Input is seurat object of a single sample after mark_emptydroplets_cellranger() 
# Output is seurat object with DoubletFinder column (having Singlet/Doublet) added to metadata
doublet_finder <- function(sample.seurat){
  
  # Filter out empty droplets before doublet identification
  subset.seurat <- subset(x = sample.seurat,
                          subset = (DropletUtils == "Empty Droplet" & CellRanger == "Empty Droplet"),
                          invert = TRUE)
  
  # Preprocess each sample
  subset.seurat <- Seurat::NormalizeData(subset.seurat)
  subset.seurat <- Seurat::FindVariableFeatures(subset.seurat)
  subset.seurat <- Seurat::ScaleData(subset.seurat)
  subset.seurat <- Seurat::RunPCA(subset.seurat)
  
  # Find significant PCs
  stdev_pc <- subset.seurat@reductions$pca@stdev
  percent_stdev_pc <- (stdev_pc/sum(stdev_pc))*100
  cumulative_stdev_pc <- cumsum(percent_stdev_pc)
  pc1 <- which(cumulative_stdev_pc > 90 & percent_stdev_pc < 5)[1]
  pc2 <- sort(which((percent_stdev_pc[1:(length(percent_stdev_pc)-1)] - 
                       percent_stdev_pc[2:length(percent_stdev_pc)]) > 0.1),
              decreasing = TRUE)[1] + 1
  min_pc <- min(pc1, pc2)
  
  # pK Identification (no prior info)
  # Introduces artificial doublets in varying proportions into real dataset,
  # preprocesses the data and calculates proportion of artificial nearest 
  # neighbors. Output is a list of proportions of artificial nearest neighbors
  # for varying combinations of pK and pN. Optimal pK is the max of bimodality
  # coefficient (BCmvn) distribution
  subset.seurat <- Seurat::RunUMAP(subset.seurat, dims = 1:min_pc)
  subset.seurat <- Seurat::FindNeighbors(subset.seurat, dims = 1:min_pc)
  subset.seurat <- Seurat::FindClusters(subset.seurat, resolution = 0.1)
  sweep.res <- DoubletFinder::paramSweep(subset.seurat, PCs = 1:min_pc, sct = FALSE)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res, GT=FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  optimal_pK <- bcmvn %>% 
    dplyr::slice_max(order_by = BCmetric) %>%
    dplyr::select(pK)
  optimal_pK <- as.numeric(as.character(optimal_pK[[1]]))
  
  # pN
  default_pN <- 0.25
  
  # Homotypic doublet estimation
  # https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled
  # From the above link, we can see that the multiplet rate is 8*10^-6 per cell
  multiplet_rates_10X <- 8*10^-6*nrow(subset.seurat@meta.data)
  nExp_poi <- round(multiplet_rates_10X*nrow(subset.seurat@meta.data))
  annotations <- subset.seurat@meta.data$seurat_clusters
  homotypic.prop <- DoubletFinder::modelHomotypic(annotations)
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # # Convert to seurat v3 object
  # subset.seurat_v3 <- scCustomize::Convert_Assay(seurat_object = subset.seurat, convert_to = "V3")
  
  # Estimate doublets
  subset.seurat <- DoubletFinder::doubletFinder(seu = subset.seurat,
                                                PCs = 1:min_pc,
                                                pN = default_pN, 
                                                pK = optimal_pK,
                                                nExp = nExp_poi.adj)
  
  # Rename column name
  colnames(subset.seurat@meta.data)[grepl(pattern = "DF.classifications", x = colnames(subset.seurat@meta.data))] <- "DoubletFinder"
  
  # Extract Cell and DoubletFinder
  df <- subset.seurat@meta.data %>% 
    dplyr::mutate(Cell = rownames(.)) %>%
    dplyr::select(Cell, DoubletFinder) %>%
    dplyr::mutate(DoubletFinder = stringr::str_to_title(DoubletFinder))
  
  # Add DoubletFinder to metadata
  sample.seurat@meta.data <- sample.seurat@meta.data %>%
    dplyr::mutate(Cell = rownames(.)) %>%
    dplyr::left_join(df, by=c("Cell"="Cell")) %>%
    dplyr::mutate(DoubletFinder = dplyr::case_when(is.na(DoubletFinder) ~ "Empty Droplet",
                                                   TRUE ~ DoubletFinder)) %>%
    dplyr::mutate(index = Cell) %>%
    tibble::column_to_rownames(var = "index")
  
  cat(paste0("DoubletFinder doublets identified for '", as.character(unique(sample.seurat@meta.data$orig.ident)), "'\n"))
  return(sample.seurat)
}

### Identify doublets using scDblFinder()
# Input is seurat object of a single sample after doublet_finder()
# Output is seurat object with scDblFinder column (having Singlet/Doublet) added to metadata
scdbl_finder <- function(sample.seurat){
  
  # Filter out empty droplets before doublet identification
  subset.seurat <- subset(x = sample.seurat,
                          subset = (DropletUtils == "Empty Droplet" & CellRanger == "Empty Droplet"),
                          invert = TRUE)
  
  # Convert to single cell experiment object
  sce.subset.seurat <- Seurat::as.SingleCellExperiment(x = subset.seurat)
  
  # Estimate doublets
  scDbl <- scDblFinder::scDblFinder(sce = sce.subset.seurat, 
                                    clusters = NULL,
                                    samples = NULL,
                                    dbr = NULL)
  
  # Extract Cell and scDblFinder
  df <- scDbl@colData@listData %>% 
    data.frame() %>% 
    dplyr::rename(scDblFinder = scDblFinder.class) %>% 
    dplyr::mutate(Cell = scDbl@colData@rownames) %>%
    dplyr::select(Cell, scDblFinder) %>%
    dplyr::mutate(scDblFinder = stringr::str_to_title(scDblFinder))
  
  # Add scDblFinder to metadata
  sample.seurat@meta.data <- sample.seurat@meta.data %>%
    dplyr::mutate(Cell = rownames(.)) %>%
    dplyr::left_join(df, by=c("Cell"="Cell")) %>%
    dplyr::mutate(scDblFinder = dplyr::case_when(is.na(scDblFinder) ~ "Empty Droplet",
                                                 TRUE ~ scDblFinder)) %>%
    dplyr::mutate(index = Cell) %>%
    tibble::column_to_rownames(var = "index")
  
  cat(paste0("scDblFinder doublets identified for '", as.character(unique(sample.seurat@meta.data$orig.ident)), "'\n"))
  return(sample.seurat)
}

### Calculate cell-level QC metrics
# Input is seurat object of a single sample after scdbl_finder()
# Ouput is seurat object with QC metrics added to metadata
calc_qc_metrics <- function(sample.seurat){
  
  # Compute percent mito percent
  sample.seurat <- Seurat::PercentageFeatureSet(object = sample.seurat,
                                                pattern = "^[Mm][Tt]-",
                                                features = NULL,
                                                col.name = "MitoPercent",
                                                assay = "RNA")
  
  # Compute percent ribo percent
  sample.seurat <- Seurat::PercentageFeatureSet(object = sample.seurat,
                                                pattern = "^[Rr][Pp][SsLl]", 
                                                features = NULL,
                                                col.name = "RiboPercent",
                                                assay = "RNA")
  
  # Compute percent hemoglobin percent
  sample.seurat <- Seurat::PercentageFeatureSet(object = sample.seurat,
                                                pattern = "^[Hh][Bb][AaBb]-", 
                                                features = NULL,
                                                col.name = "HemePercent",
                                                assay = "RNA")
  
  # Extract metadata
  sample_metadata <- sample.seurat@meta.data
  
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
  sample_metadata <- sample_metadata %>% 
    dplyr::mutate(Cell = paste0(orig.ident, "_", rownames(sample_metadata)),
                  Sample = orig.ident,
                  nUMIs = nCount_RNA,
                  nGenes = nFeature_RNA,
                  MitoRatio = MitoPercent/100,
                  RiboRatio = RiboPercent/100,
                  HemeRatio = HemePercent/100,
                  Novelty = log10(nGenes)/log10(nUMIs), .keep="unused")
  
  # If HTO tag info is available in metadata of seurat object, rename 
  # nCount_HTO, nFeature_HTO and HTO_Final columns in metadata
  if (sum(colnames(sample.seurat@meta.data) %in% c("nCount_HTO", "nFeature_HTO", "HTO_Final")) > 0){
    sample_metadata <- sample_metadata %>% 
      dplyr::mutate(nHTO_UMIs = nCount_HTO,
                    nHTOs = nFeature_HTO,
                    HTO_Final = HTO_Final)
  } else {
    sample_metadata <- sample_metadata %>% 
      dplyr::mutate(nHTO_UMIs = 0,
                    nHTOs = 0,
                    HTO_Final = NA) 
  }
  
  # Replace metadata in raw Seurat object with updated column names
  sample.seurat@meta.data <- sample_metadata %>%
    dplyr::select(Cell, Sample, nUMIs, nGenes, nHTO_UMIs, nHTOs, HTO_Final, 
                  MitoRatio, RiboRatio, HemeRatio, Novelty, DropletUtils, 
                  CellRanger, DoubletFinder, scDblFinder)
  
  cat(paste0("Cell-level QC metrics calculated for '", as.character(unique(sample.seurat@meta.data$Sample)), "'\n"))
  return(sample.seurat)
}

### Identify low quality cells
# Input is seurat object of a single sample after calc_qc_metrics()
# Ouput is seurat object with column (having Singlet/Empty Droplet/Doublet/Low Quality) added to metadata
mark_low_quality <- function(sample.seurat){
  
  # These are "very lenient" hard-cutoffs to mark poor quality cells
  gene_cutoff <- 250
  umi_cutoff <- 500
  mito_cutoff <- 0.2
  ribo_cutoff <- 0.05
  novelty_cutoff <- 0.8
  
  # Mark the poor quality cells
  sample.seurat@meta.data <- sample.seurat@meta.data %>% 
    dplyr::mutate(QC = dplyr::case_when((DropletUtils == "Empty Droplet" & CellRanger == "Empty Droplet") ~ "Empty Droplet",
                                        (DoubletFinder == "Doublet" & scDblFinder == "Doublet") ~ "Doublet",
                                        (nGenes >= gene_cutoff & nUMIs >= umi_cutoff & MitoRatio <= mito_cutoff & Novelty >= novelty_cutoff) ~ "Singlet", #RiboRatio >= ribo_cutoff &
                                        TRUE ~ "Low Quality"))
  
  cat(paste0("Good quality singlets identified for '", as.character(unique(sample.seurat@meta.data$Sample)), "'\n"))
  return(sample.seurat)
}

### Append raw metadata from each sample for plotting purpose
# Input is a seurat object of a single sample and dataframe containing raw metadata after mark_low_quality() but before filter_singlets()
# Output is dataframe with raw metadata of current sample added
generate_plotdata <- function(sample.seurat, raw_metadata){
  
  # Append the raw metadata of each sample which will be used for QC plots later
  raw_metadata <- dplyr::bind_rows(raw_metadata, sample.seurat@meta.data) %>%
    dplyr::filter(!is.na(Sample))
  
  cat(paste0("Raw metadata (for plotting purpose) generated for '", as.character(unique(sample.seurat@meta.data$Sample)), "'\n"))
  return(raw_metadata)
}

### Remove low quality cells
# Input is seurat object of a single sample after mark_low_quality()
# Ouput is seurat object with ONLY singlets
filter_singlets <- function(sample.seurat){
  
  # Keep ONLY singlets
  sample.seurat <- base::subset(x = sample.seurat,
                                subset = (QC == "Singlet"))
  
  cat(paste0("Good quality singlets retained for '", as.character(unique(sample.seurat@meta.data$Sample)), "'\n"))
  return(sample.seurat)
}

### Add additional metadata from metafile, merge into single seurat object & save
# Input is a list of samples, path to metafile, metafile name & path to store the filtered seurat object
# Output is a single filtered seurat object
format_filtered <- function(samples, output_path){
  
  # Create a merged seurat object after all the above filtering
  # NOTE: Samples will have the same barcodes. To keep track of cell identities
  # (i.e. barcodes) coming from each sample after merging, we add a prefix
  # (i.e. sample name) to each barcode using "add.cell.ids"
  samples.seurat <- lapply(samples, get)
  filtered.seurat <- base::merge(x = samples.seurat[[1]],   #get(paste0(samples[1])
                                 y = samples.seurat[-1],    #lapply(paste0(samples[2:length(samples)]), get)
                                 add.cell.ids = samples,
                                 merge.data = FALSE)
  
  # Remove HTO assay if present to avoid complications during integration
  if (sum(Assays(filtered.seurat) %in% "HTO") > 0){
    filtered.seurat[["HTO"]] <- NULL
  }
  
  # Import any other meta data associated with data set
  # NOTE: This xslx file should have column named "Unique_ID" whose values matches 
  # with  column "Unique_ID" of seurat object's metadata.
  extra_metadata <-  openxlsx::read.xlsx(xlsxFile = paste0(scripts_path, metafile)) %>%
    dplyr::select(everything(), -Comments)
  
  # Merge imported metadata with existing metadata
  # NOTE: Add row names before replacing metadata in Seurat object as left_join 
  # will remove row names.
  filtered.seurat@meta.data <- filtered.seurat@meta.data %>%
    dplyr::mutate(Unique_ID = dplyr::case_when(!is.na(HTO_Final) ~ paste0(Sample, "_", HTO_Final), 
                                               is.na(HTO_Final) ~ paste0(Sample))) %>%
    dplyr::left_join(extra_metadata, by=("Unique_ID"="Unique_ID")) %>%
    dplyr::mutate(index = Cell) %>%
    tibble::column_to_rownames(var = "index")
  
  # Create .rds object for filtered seurat object
  saveRDS(filtered.seurat, file=paste0(output_path, "filtered.seurat.rds"))
  
  cat("Filtered seurat object saved\n")
  return(filtered.seurat)
}

### Generate QC plots
# Input is dataframe containing raw metadata & path to store the QC plots
# Output are plots
plot_qc <- function(raw_metadata, output_path){
  
  ### Visualize the number of cell counts per sample
  cell_qc <- function(meta){
    
    metadata <- meta %>%
      dplyr::count(Sample, QC) %>%
      data.frame() %>%
      dplyr::mutate(QC = factor(QC, levels = c("Empty Droplet","Doublet", "Low Quality", "Singlet")))
    
    p <- ggplot(data = metadata, aes(x = Sample, y = n, fill = QC)) + 
      # position = "dodge" for grouped; "stack" for stacked
      # stat = "identity" if y axis defined; "count" if y axis determined based on X axis frequency
      geom_bar(position = position_dodge(0.9), stat = "identity", drop = FALSE) +             
      theme_classic() +               # display with x and y axis lines and no gridlines
      my_theme +
      labs(x = "Sample", y = "Cell Counts", title = "Number of Cells") +
      coord_cartesian(ylim = c(1,10000000), clip = "off", expand = FALSE) +
      scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000, 1000000)) +
      scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
                                   "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
      geom_text(aes(label = n, ymin = 0.1, ymax = 1), 
                position = position_dodge(width = 0.9), y = 0.1, hjust = 0, angle = 90)
    #geom_text(stat ="count", aes(label = after_stat(count)), y = 0, hjust = 0, angle = 90)
    
    return(p)
  }
  
  ### Visualize the number of UMIs per cell
  umi_qc <- function(meta){
    
    umi_cutoff <- 500
    metadata <- meta %>%
      dplyr::mutate(QC = factor(QC, levels = c("Empty Droplet","Doublet", "Low Quality", "Singlet")))
    
    p <- ggplot(data = metadata, aes(x = Sample, y = nUMIs, fill = QC)) +
      geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) + 
      geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) +
      theme_classic() +  
      my_theme +
      labs(x = "Sample", y = "Number of UMIs", title = "Distribution of UMIs") +
      coord_cartesian(ylim = c(1,1000000), clip = "off") +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) +
      scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
                                   "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
      geom_hline(yintercept = umi_cutoff, linetype = 2)
    
    return(p)
  }
  
  ### Visualize the number of genes per cell
  gene_qc <- function(meta){
    
    gene_cutoff <- 250
    metadata <- meta %>%
      dplyr::mutate(QC = factor(QC, levels = c("Empty Droplet","Doublet", "Low Quality", "Singlet")))
    
    p <- ggplot(data = metadata, aes(x = Sample, y = nGenes, fill = QC)) +
      geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) + 
      geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) +
      theme_classic() + 
      my_theme + 
      labs(x = "Sample", y = "Number of Genes", title = "Distribution of Genes") +
      coord_cartesian(ylim = c(1, 30000), clip = "off") +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000)) + 
      scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
                                   "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
      geom_hline(yintercept = gene_cutoff, linetype = 2)
    
    return(p)
  }
  
  # Visualize the MitoRatio of each cell
  mito_qc <- function(meta){
    
    mito_cutoff <- 0.2
    metadata <- meta %>%
      dplyr::mutate(QC = factor(QC, levels = c("Empty Droplet","Doublet", "Low Quality", "Singlet")))
    
    p <- ggplot(data = metadata, aes(x = Sample, y = MitoRatio, fill = QC)) +
      geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) + 
      geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) + 
      theme_classic() + 
      my_theme +
      labs(x = "Sample", y = "MitoRatio", title = "Distribution of MitoRatio") +
      coord_cartesian(ylim = c(0.00001, 1), clip = "off") +
      scale_y_log10(breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1)) + 
      scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
                                   "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
      geom_hline(yintercept = mito_cutoff, linetype = 2)
    
    return(p)
  }
  
  ### Visualize the RiboRatio of each cell
  ribo_qc <- function(meta){
    
    ribo_cutoff <- 0.05
    metadata <- meta %>%
      dplyr::mutate(QC = factor(QC, levels = c("Empty Droplet","Doublet", "Low Quality", "Singlet")))
    
    p <- ggplot(data = metadata, aes(x = Sample, y = RiboRatio, fill = QC)) +
      geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) + 
      geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) + 
      theme_classic() + 
      my_theme +
      labs(x = "Sample", y = "RiboRatio", title = "Distribution of RiboRatio") +
      coord_cartesian(ylim = c(0.0001, 1), clip = "off") +
      scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) + 
      scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
                                   "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
      geom_hline(yintercept = ribo_cutoff, linetype = 2)
    
    return(p)
  }
  
  ### Visualize the novelty or complexity of each cell
  novelty_qc <- function(meta){
    
    novelty_cutoff <- 0.8
    metadata <- meta %>%
      dplyr::mutate(QC = factor(QC, levels = c("Empty Droplet","Doublet", "Low Quality", "Singlet")))
    
    p <- ggplot(data = metadata, aes(x = Sample, y = Novelty, fill = QC)) +
      geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) + 
      geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) + 
      theme_classic() + 
      my_theme +
      labs(x = "Sample", y = "Novelty", title = "Distribution of Novelty Score") +
      coord_cartesian(ylim = c(0.3, 1), clip = "off") +
      scale_y_log10(breaks = c(0.3, 1)) + 
      scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
                                   "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
      geom_hline(yintercept = novelty_cutoff, linetype = 2)
    
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
    
    metadata <- meta %>%
      dplyr::mutate(QC = factor(QC, levels = c("Empty Droplet","Doublet", "Low Quality", "Singlet")))
    
    p <- ggplot(data = metadata, aes(x = nUMIs, y = nGenes, color = MitoRatio)) +
      geom_point() +
      theme_classic() + 
      my_theme + 
      labs(x = "Number of UMIs", y = "Number of Genes",	 title = "Distribution of UMIs, Genes & MitoRatio") +
      coord_cartesian(xlim = c(1, 1000000), ylim = c(1, 20000), clip = "off") +
      scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) + 
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000)) + 
      #facet_wrap(.~Sample, nrow = 4) +   #split the plot by X-axis label
      facet_wrap(Sample ~ QC, ncol = 4) +
      stat_smooth(method=lm, color="yellow") +
      geom_vline(xintercept = umi_cutoff) +    	#draw a vertical line at x=500 i.e.UMIs cutoff
      geom_hline(yintercept = gene_cutoff) +    #draw a horizontal line at y =250 i.e. Genes cutoff
      scale_color_viridis(option = "D", limits = c(0, 1)) 		# limits sets max and min values of gradient 
    
    return(p)
  }
  
  # Plot all QC metrics before and after QC
  funcs <- c("cell_qc", "umi_qc", "gene_qc", "mito_qc", "ribo_qc", "novelty_qc",
             "gene_umi_mito_qc")
  
  filenames <- c("Cell_Counts", "UMI_Distribution", "Gene_Distribution",
                 "MitoRatio_Distribution", "RiboRatio_Distribution", 
                 "Novelty_Score_Distribution", "Genes_UMI_MitoRatio_Distribution")
  
  for (i in 1:length(funcs)){
    
    # # Plot QC metrics
    # purrr::map(.x = c("raw_metadata"), .f = get(funcs[i])) %>% 
    #   cowplot::plot_grid(plotlist = ., align = "hv", axis = "tblr", nrow = 1, ncol = 1)
    
    p <- get(funcs[i])(raw_metadata)
    
    # Save the plot
    ggplot2::ggsave(filename = paste0("QC_", filenames[i], ".pdf"),
                    plot = p, #last_plot(),
                    device = "pdf",
                    path = output_path,
                    #scale = 1,
                    width = 11,
                    height = 8,
                    units = c("in"),	 
                    dpi = 600,
                    limitsize = TRUE,
                    bg = NULL)
  }
  
  cat("QC plots generated\n")  
}

### Perform SCTransofrmation
# Input is filtered seurat object
# Output is seurat object after SCTransformation & PCA reduction
sctransform_singlecell <- function(filtered.seurat, output_path){
  
  # NOTE: In v3, we perform SCTransform on each sample.seurat object separately.
  # In v5, we perform SCTransform on each sample stored in layers in RNA assay.
  
  # Normalize the data before cell cycle scoring
  filtered.seurat <- Seurat::NormalizeData(object = filtered.seurat,
                                           assay = "RNA",
                                           normalization.method = "LogNormalize",
                                           scale.factor = 10000,
                                           margin = 1,
                                           verbose = TRUE)
  
  # CellCycleScoring uses a single data layer (i.e. log norm counts) 
  # However, data layer for each sample is stored separately. Merge them first.
  filtered.seurat@assays$RNA <- SeuratObject::JoinLayers(filtered.seurat@assays$RNA)
  
  # Perform cell cycle scoring
  filtered.seurat <- Seurat::CellCycleScoring(object = filtered.seurat,
                                              s.features = intersect(s_genes,rownames(filtered.seurat@assays$RNA@features)),
                                              g2m.features = intersect(g2m_genes, rownames(filtered.seurat@assays$RNA@features)),
                                              ctrl = NULL)
  
  # Regress out the difference between the G2M and S phase scores
  # https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette
  filtered.seurat$CC.Score <- filtered.seurat$G2M.Score-filtered.seurat$S.Score
  
  # SCTransform MUST be performed on each sample INDIVIDUALLY. Split them first.
  filtered.seurat@assays$RNA <- base::split(x = filtered.seurat@assays$RNA, 
                                            f = filtered.seurat$Sample)
  
  # Perform normalization, variable feature identification & scaling 
  # https://github.com/satijalab/seurat/issues/7342
  sct.seurat <- Seurat::SCTransform(object = filtered.seurat,
                                    assay = "RNA",
                                    new.assay.name = "SCT",
                                    do.correct.umi = TRUE,
                                    ncells = 5000,
                                    variable.features.n = 3000,
                                    vars.to.regress = c("CC.Score","MitoRatio"),
                                    do.scale = FALSE,
                                    do.center = TRUE,
                                    vst.flavor = "v2",
                                    return.only.var.genes = TRUE,
                                    verbose = FALSE)
  
  # Remove ribosomal, Riken, predicted, mitochondrial genes from VariableFeatures
  # so that PCA, UMAP and hence clustering are not influenced by these genes
  var_f <- sct.seurat@assays$SCT@var.features
  var_f <- var_f[!grepl(pattern = "^[Rr][Pp][SsLl]|R[Ii][Kk]$|^[Mm][Tt]-|^G[Mm][0-9.]+$", 
                        x = var_f)]
  
  # Replace variable features of SCT assay
  sct.seurat@assays$SCT@var.features <- var_f
  cat("\nFinal number of variable features:", length(var_f), "\n")
  
  # Scale data & run PCA on RNA assay (Needed for scVI integration)
  sct.seurat <- Seurat::ScaleData(object = sct.seurat, 
                                  assay = "RNA",
                                  features = VariableFeatures(sct.seurat))
  sct.seurat <- Seurat::RunPCA(object = sct.seurat,
                               assay = "RNA",
                               features = VariableFeatures(sct.seurat),
                               reduction.name = "rna.pca",
                               reduction.key = "PC_")
  
  # Perform dimensional reduction using PCA on SCT assay variable features
  sct.seurat <- Seurat::RunPCA(object = sct.seurat,
                               assay = "SCT",
                               features = VariableFeatures(sct.seurat),
                               reduction.name = "sct.pca",
                               reduction.key = "PC_")
  
  # Create .rds object for sct seurat object
  saveRDS(sct.seurat, file=paste0(output_path, "sct.seurat.rds"))
  
  cat("scTransform completed", "\n")
  return(sct.seurat)
}

### Perform Integration
# Input is seurat object after SCTransformation & PCA reduction
# Output is seurat object after integration
integrate_singlecell <- function(sct.seurat, reference.samples, kweight, output_path){
  
  cat("Reference.samples:", reference.samples, "\n")
  cat("kweight:", kweight, "\n")
  
  integrated.seurat <- sct.seurat
  for (r in c("CCA", "RPCA", "Harmony", "JointPCA")){
    integrated.seurat <- Seurat::IntegrateLayers(object = integrated.seurat,
                                                 method = paste0(r, "Integration"),
                                                 normalization.method = "SCT",
                                                 orig.reduction = "sct.pca", 
                                                 new.reduction = paste0("integrated.", base::tolower(r)),
                                                 reference = reference.samples,
                                                 k.weight = kweight,    # for RPCA
                                                 verbose = FALSE)
  }
  
  # NOTE: scVI needs raw counts. Vignette also uses it on RNA assay
  # NOTE: We use variable features of SCT assay for integration.
  # NOTE: We use pca reduction from RNA assay (derived using variable features of SCT assay)
  # FastMNN throws error "Error in checkBatchConsistency(batches, cells.in.columns = TRUE)"
  # for (r in c("scVI", "FastMNN")){
  #   DefaultAssay(integrated.seurat) <- "RNA"
  #   integrated.seurat <- Seurat::IntegrateLayers(object = integrated.seurat,
  #                                                method = paste0(r, "Integration"),
  #                                                normalization.method = "LogNormalize",
  #                                                orig.reduction = "rna.pca",
  #                                                features = integrated.seurat@assays$SCT@var.features,
  #                                                new.reduction = paste0("integrated.", base::tolower(r)),
  #                                                reference = ref_samples,
  #                                                k.weight = kweight,                                    # for RPCA
  #                                                conda_env = "/hpc/home/kailasamms/miniconda3/envs/R",  # for scVI
  #                                                verbose = FALSE)
  # }
  
  # Create .rds object for integrated seurat object
  saveRDS(integrated.seurat, file=paste0(output_path, "integrated.seurat.rds"))
  
  cat("Integration completed", "\n")
  return(integrated.seurat)
}

### Perform clustering
# Input is seurat object after integration
# Output is seurat object after clustering
cluster_singlecell <- function(integrated.seurat, output_path){
  
  #***************STEP 8A: FIND NEAREST NEIGHBORS FOR EVERY CELL***************#
  
  # Determine the K-nearest neighbor graph
  for (r in c("CCA", "RPCA", "Harmony", "JointPCA")){
    integrated.seurat <- Seurat::FindNeighbors(object = integrated.seurat,
                                               reduction = paste0("integrated.", base::tolower(r)),
                                               dims = 1:40,
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
                                                algorithm = 4,     #4=Leiden is best
                                                method = "matrix")
    }
  }
  
  #**********STEP 8C: PERFORM DIMENSIONAL REDUCTION FOR VISUALIZATION**********#
  
  # Run UMAP
  for (r in c("CCA", "RPCA", "Harmony", "JointPCA")){ 
    integrated.seurat <- Seurat::RunUMAP(object = integrated.seurat,
                                         dims = 1:40,
                                         n.neighbors = 30L,
                                         reduction = paste0("integrated.", base::tolower(r)),
                                         reduction.name = paste0("umap.", base::tolower(r)))
  }
  
  #**************************STEP 8D: MERGE ALL LAYERS*************************#
  
  integrated.seurat@assays$RNA <- SeuratObject::JoinLayers(integrated.seurat@assays$RNA)
  
  # Create .rds object for integrated seurat object
  saveRDS(integrated.seurat, file=paste0(output_path, "integrated.seurat.rds"))
  
  cat("Clustering completed", "\n")
  return(integrated.seurat)
}

### Remove sparse clusters
# Input is seurat object after clustering
# Output is seurat object after removing clusters with less than 5 cells
remove_sparse_clusters <- function(integrated.seurat, output_path){
  
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
  saveRDS(integrated.seurat, file=paste0(output_path, "integrated.seurat.rds"))
  
  cat("Integrated seurat object saved after removing sparse clusters (below 5 cells)", "\n")
  return(integrated.seurat) 
}

### Plot metrics post integration
# Input is seurat object after clustering & removal of sparse clusters
# Output is (i) series of UMAPs at resolution 0.8 using Harmony reduction
# (ii) series of UMAPs at different resolution using every available reduction 
plot_metrics_post_integration <- function(integrated.seurat, suffix, output_path){
  
  # File names for each of 10 figures
  filenames <- paste0(c("Pre.Integration.PCA.", "Post.Integration.PCA.", "UMAP.Sample.", "UMAP.Phase.", 
                        "UMAP.All.Resolutions.CCA.", "UMAP.All.Resolutions.RPCA.", 
                        "UMAP.All.Resolutions.JointPCA.", "UMAP.All.Resolutions.Harmony.",
                        "UMAP.Singlets.Doublets.", "UMAP.Numerical.Metrics."), suffix)
  
  # Reductions to be used for each of 10 figures
  reductions <- c("sct.pca", "integrated.harmony", "umap.harmony", "umap.harmony",
                  "umap.cca", "umap.rpca", "umap.jointpca", "umap.harmony",
                  "umap.harmony", "umap.harmony")
  
  # Variable on which seurat object needs to be split for each of 10 figures
  splits <- c("Sample", "Sample", "Sample", "Phase", NA, NA, NA, NA, NA, NA)
  
  for (i in 1:length(filenames)){ 
    
    reduction.parameter <-  reductions[i]
    
    # Plotting PCA (Pre, Post- integration), UMAP (Sample, Phase) for each sample at Harmony 0.8
    if (splits[i] %in% c("Sample", "Phase")){
      plot.seurat <- Seurat::SplitObject(object = integrated.seurat,
                                         split.by = splits[i])
      
      purrr::map(.x = c(1:length(plot.seurat)),
                 .f = function(x){  
                   Idents(plot.seurat[[x]]) <- "cluster.0.8.harmony"
                   Seurat::DimPlot(object = plot.seurat[[x]],
                                   reduction = reduction.parameter,
                                   group.by = "cluster.0.8.harmony",
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
                                           nrow=ceiling(length(plot.seurat)/3),
                                           ncol=3,
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
                      width = 4*3,
                      height = 4*ceiling(length(plot.seurat)/3),
                      units = c("in"),
                      dpi = 600,
                      limitsize = TRUE,
                      bg = "white")
    }
    
    # Plotting UMAP for CCA, RPCA, JointPCA & Harmony at each resolution
    else if (is.na(splits[i]) & i < 9){
      
      purrr::map(.x = c(0.4, 0.6, 0.8, 1, 1.2, 1.4),
                 .f = function(x){  
                   idents <- paste0("cluster.", x, gsub(pattern="umap", replacement="", x=reduction.parameter))
                   Idents(integrated.seurat) <- idents
                   Seurat::DimPlot(object = integrated.seurat,
                                   reduction = reduction.parameter,
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
      
      # Save the plot
      ggplot2::ggsave(filename = paste0(filenames[i], ".tiff"),
                      plot = last_plot(),
                      device = "jpeg",
                      path = diagnostics_path,
                      scale = 1,
                      width = 4*3,
                      height = 4*2,
                      units = c("in"),
                      dpi = 600,
                      limitsize = TRUE,
                      bg = "white")
    }
    
    # Plotting UMAP of singlets, doublets, low quality cells at Harmony 0.8
    else if (is.na(splits[i]) & i == 9){
      
      purrr::map(.x = c("DropletUtils", "CellRanger", "DoubletFinder", "scDblFinder", "QC"),
                 .f = function(x){ 
                   Idents(integrated.seurat) <- "cluster.0.8.harmony"
                   Seurat::DimPlot(object = integrated.seurat,
                                   reduction = reduction.parameter,
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
      
      # Save the plot
      ggplot2::ggsave(filename = paste0(filenames[i], ".tiff"),
                      plot = last_plot(),
                      device = "jpeg",
                      path = diagnostics_path,
                      scale = 1,
                      width = 4*3,
                      height = 4*2,
                      units = c("in"),
                      dpi = 600,
                      limitsize = TRUE,
                      bg = "white")
    }
    
    # Plotting UMAP of numerical metrics at Harmony 0.8
    else{
      purrr::map(.x = c("nUMIs", "nGenes", "S.Score", "G2M.Score", "CC.Score", "MitoRatio"),
                 .f = function(x){ 
                   Idents(integrated.seurat) <- "cluster.0.8.harmony"
                   Seurat::FeaturePlot(object = integrated.seurat,
                                       features = x,
                                       reduction = reduction.parameter,
                                       pt.size = 0.1,
                                       min.cutoff='q10',
                                       order = TRUE,  # plot doublets on above rest of cells
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
      
      # Save the plot
      ggplot2::ggsave(filename = paste0(filenames[i], ".tiff"),
                      plot = last_plot(),
                      device = "jpeg",
                      path = diagnostics_path,
                      scale = 1,
                      width = 4*3,
                      height = 4*2,
                      units = c("in"),
                      dpi = 600,
                      limitsize = TRUE,
                      bg = "white")
      
    }
  }
}

### Identify markers for each cluster
# Input is seurat object after clustering & removal of sparse clusters
identify_markers <- function(integrated.seurat, resolution, reduction, suffix, output_path){
  
  # Set default assay
  DefaultAssay(integrated.seurat) <- "RNA"
  
  # Change active.ident
  idents <- paste0("cluster.", resolution, ".", base::tolower(reduction))
  Idents(object=integrated.seurat) <- idents
  
  # Find ALL markers
  all_markers <- Seurat::FindAllMarkers(object=integrated.seurat,
                                        assay="RNA",
                                        features=NULL,
                                        logfc.threshold=0.25,
                                        test.use="wilcox",
                                        slot="data",
                                        min.pct=0.1,
                                        min.diff.pct=0.1,
                                        node=NULL,
                                        verbose=TRUE,
                                        only.pos=TRUE,
                                        max.cells.per.ident=Inf,
                                        random.seed=1,
                                        latent.vars=NULL,
                                        min.cells.feature=3,
                                        min.cells.group=1,
                                        pseudocount.use=1,
                                        mean.fxn=NULL,
                                        fc.name=NULL,
                                        base=2,
                                        return.thresh=0.01,
                                        densify=FALSE)
  
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
    dplyr::left_join(y=unique(annotations[, c("SYMBOL", "CHR", "DESCRIPTION")]), by=c("gene"="SYMBOL")) %>%
    dplyr::relocate(cluster, gene, CHR, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, ratio, DESCRIPTION)
  
  # Find top markers for each major cluster
  top_markers <- all_markers %>%
    dplyr::filter(avg_log2FC >= 0.58 & p_val_adj < 0.05) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%  #desc(ratio)
    dplyr::slice_head(n=30) %>%
    ungroup()
  
  # Save all the markers
  filename <- paste0(proj, ".Markers.All.", idents, ".", suffix,".xlsx")
  wb <- openxlsx::createWorkbook()
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

### Annotate clusters
# Input is seurat object after clustering & removal of sparse clusters
# Output is seurat object with cell annotations based on our marker file
# NOTE: Cluster annotation can be performed in multiple ways: 
# (i) manually rename each cluster (needs to be done for each dataset, time consuming)
# (ii) automatically using seurat and ucell module scores (RECOMMENDED)

# STEP1: Annotate each cell based on seurat and ucell module scores
# STEP2: Determine which celltypes have highest percentage within each cluster
# STEP3: Re-annotate all cells within each cluster to match STEP2 classification













#******************************************************************************#
#                               GENE ANNOTATIONS                               #
#******************************************************************************#

# This function returns 2 dataframes: first for human & second for mouse with 
# following columns: "ENSEMBL_ID", "ENTREZ_ID", "SYMBOL", "SYMBOL_ENTREZ", 
# "BIOTYPE", "BIOTYPE_ENTREZ", "START", "END", "CHR", "STRAND", "DESCRIPTION"

get_annotations <- function(){
  
  # AnnotationHub has SYMBOL-ENSEMBL_ID info ONLY.
  # AnnotationDbi has SYMBOL-ENSEMBL_ID as well as SYMBOL-ENTREZ_ID info.
  # hubCache(AnnotationHub()) to find location where cache is stored and delete
  # it and start fresh if you get errors like "Error: failed to load resource"
  
  for (species in c("Homo sapiens", "Mus musculus")){
    
    #**************************GET ENSEMBL ANNOTATIONS***************************#
    # Connect to AnnotationHub
    ah <- AnnotationHub::AnnotationHub()
    
    # Access the Ensembl database for organism
    ahDb <- AnnotationHub::query(x=ah,
                                 pattern=c(species, "EnsDb"),
                                 ignore.case=TRUE)
    
    # Acquire the latest annotation files
    id <- ahDb %>%
      mcols() %>%
      rownames() %>%
      tail(n=1)
    
    # Download the appropriate Ensembldb database
    edb <- ah[[id]]
    
    # Extract gene-level information from database
    ensembl <- ensembldb::genes(x=edb,
                                return.type="data.frame")
    
    # Select annotations of interest
    ensembl <- ensembl %>%
      dplyr::rename(ENSEMBL_ID=gene_id, SYMBOL=gene_name, 
                    BIOTYPE=gene_biotype, START=gene_seq_start, END=gene_seq_end, 
                    CHR=seq_name, STRAND=seq_strand, DESCRIPTION=description,
                    ENSEMBL_TRANSCRIPT = canonical_transcript) %>%
      dplyr::mutate(SYMBOL=dplyr::case_when(nchar(SYMBOL) == 0 ~ NA,
                                            TRUE ~ SYMBOL)) %>%
      dplyr::select(ENSEMBL_ID, ENSEMBL_TRANSCRIPT, SYMBOL, BIOTYPE, START, END, CHR, STRAND, DESCRIPTION)
    
    #***************************GET ENTREZ ANNOTATIONS***************************# 
    
    # NOTE: mapIds can ONLY retrieve one of "EMSEMBL/SYMBOL/GENETYPE" at a time
    # mapping <- AnnotationDbi::mapIds(x=org.Hs.eg.db, 
    #                                  keys=keys(org.Hs.eg.db),
    #                                  keytype="ENTREZID", 
    #                                  column="SYMBOL") %>%
    #   as.data.frame(do.call(cbind, list(.))) %>%
    #   tibble::rownames_to_column("ENTREZID") %>%
    #   dplyr::rename(ENTREZID=identity(1), SYMBOL=identity(2))
    
    if (species == "Homo sapiens"){
      entrez <- AnnotationDbi::select(x=org.Hs.eg.db, keys=keys(org.Hs.eg.db),
                                      columns=c("ENSEMBL", "SYMBOL","GENETYPE"))
    } else if (species == "Mus musculus"){
      entrez <- AnnotationDbi::select(x=org.Mm.eg.db, keys=keys(org.Mm.eg.db),
                                      columns=c("ENSEMBL", "SYMBOL","GENETYPE"))
    }
    
    colnames(entrez) <- c("ENTREZ_ID", "ENSEMBL_ID", "SYMBOL_ENTREZ", "BIOTYPE_ENTREZ")
    
    # Merge ensembl and entrez dataframes
    annotations <- dplyr::full_join(ensembl, entrez, by=c("ENSEMBL_ID"="ENSEMBL_ID")) %>%
      dplyr::select(ENSEMBL_ID, ENSEMBL_TRANSCRIPT, ENTREZ_ID, SYMBOL, 
                    SYMBOL_ENTREZ, BIOTYPE, BIOTYPE_ENTREZ, START, END, CHR, 
                    STRAND, DESCRIPTION)
    
    # Save to dataframe
    if (species == "Homo sapiens"){
      df1 <- annotations
    } else {
      df2 <- annotations
    }
  }
  
  #DBI::dbDisconnect(conn=ah)
  return(list(df1, df2))
}

#******************************************************************************#
#                      SPATIAL ANALYSIS RELATED FUNCTIONS                      #       
#******************************************************************************#

sctransform_spatial <- function(filtered_seurat){
  
  # Normalize, identify variable features, SCTransform each dataset independently
  sct <- Seurat::SCTransform(object = filtered_seurat,
                             assay = "Spatial",
                             new.assay.name = "SCT",
                             reference.SCT.model = NULL,
                             do.correct.umi = TRUE,
                             ncells = 5000,
                             residual.features = NULL,
                             variable.features.n = 3000,
                             vars.to.regress = NULL,
                             do.scale = FALSE,
                             do.center = TRUE,
                             return.only.var.genes = FALSE) #important
  
  # Remove ribosomal, Riken, predicted and mitochondrial genes from
  # VariableFeatures so that PCA, UMAP and hence clustering are not affected
  var_f <- sct@assays$SCT@var.features
  var_f <- var_f[!grepl(pattern = "^[Rr][Pp][SsLl]|R[Ii][Kk]$|^[Mm][Tt]-|^G[Mm][0-9.]+$", 
                        x = var_f)]
  
  sct@assays$SCT@var.features <- var_f
  cat("\nFinal number of variable features:", length(var_f), "\n")
  
  # Perform dimensional reduction using PCA on SCT assay variable features
  sct <- Seurat::RunPCA(object = sct,
                        assay = "SCT",
                        features = NULL,
                        ndims.print = 1,
                        nfeatures.print = 1,
                        reduction.name = "pca",
                        reduction.key = "PC_")
  
  # Perform dimensional reduction using UMAP on PCA dimensions
  sct <- Seurat::RunUMAP(object = sct,
                         dims = 1:40,
                         reduction = "pca",
                         reduction.name = "umap",
                         reduction.key = "UMAP_")
  
  return(sct)
}

cluster_spatial_data <- function(integrated_seurat){
  
  # Unlike single cell data, each spatial tissue is analyzed individually.
  # So, there is no integration involved using cca, rpca, harmony etc..
  
  #***************STEP 8A: FIND NEAREST NEIGHBORS FOR EVERY CELL***************#
  
  # Determine the K-nearest neighbor graph
  integrated_seurat <- Seurat::FindNeighbors(object=integrated_seurat,
                                             reduction="pca",
                                             dims=1:40,
                                             k.param =30)
  
  #**********STEP 8B: SEPARATE CELLS INTO CLUSTERS BASED ON SNN GRAPH**********#
  
  # Determine the clusters for various resolutions
  for (res in c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4)){
    integrated_seurat <- Seurat::FindClusters(object=integrated_seurat,
                                              resolution=res,
                                              modularity.fxn=1,
                                              algorithm=3,     #4=Leiden
                                              method="matrix")
  }
  
  return(integrated_seurat)
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
                        max_color <- ifelse(column == column_1,
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

#******************************************************************************#
#                                     NOTES                                    #
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
