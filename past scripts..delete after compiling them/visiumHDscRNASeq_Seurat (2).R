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
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")

# NOTE: In single cell experiment, we combine multiple samples into a single 
# seurat object and do analysis. In spatial experiment, we analyze each sample
# individually.

# NOTE: In visium v1 spatial transcriptomics, there are 4992 barcoded spots 
# which capture RNA from tissue placed above these spots. 
# So, unlike scRNASeq, each barcode corresponds to a spot, not a cell.
# https://www.youtube.com/watch?v=VwNk4d-0RJc
# https://kb.10xgenomics.com/hc/en-us/articles/360035848191-How-many-spots-are-within-a-single-capture-area-on-the-Visium-v1-Spatial-Gene-Expression-Slide

# NOTE: Reads are obtained from both background (non-tissue) and tissue specific
# areas and stored in raw_feature_bc_matrix and filtered_fearure_bc_matrix
# respectively. Since, we are ONLY interested in reads obtained from tissues,
# use ONLY filtered_feature_bc_matrix.

#******************************************************************************#
#                       STEP 1: SETUP THE SEURAT OBJECT                        #
#******************************************************************************#

#*********************IMPORTING GEX (GENE EXPRESSION) DATA*********************#

# raw_matrix_path
# ->Sample1
#   ->binned_outputs
#     ->square_002um
#       ->raw_feature_bc_matrix.h5
#       ->filtered_feature_bc_matrix.h5
#       ->spatial
#     ->square_008um
#       ->raw_feature_bc_matrix.h5
#       ->filtered_feature_bc_matrix.h5
#       ->spatial
#     ->square_016um
#       ->raw_feature_bc_matrix.h5
#       ->filtered_feature_bc_matrix.h5
#       ->spatial

# DIRECTORY STRUCTURE IS IMPORTANT:
# data.dir MUST have a H5 file specified by filename parameter as well as folder
# named "spatial" containing the image

# Create a list of sample names which will be added to each barcode.
# Since folder names correspond to sample name, we just use list.files()
samples <- list.files(path = filt_matrix_path)

# Loop through each of the individual folders in parent directory & import data
for(i in samples){
  
  sample.seurat <- Load10X_Spatial(data.dir = paste0(filt_matrix_path, i),
                                   filename = "filtered_feature_bc_matrix.h5",
                                   assay = "Spatial",
                                   slice = i,
                                   bin.size = c(2,8,16),
                                   filter.matrix = TRUE,
                                   to.upper = FALSE,
                                   image = NULL)
  
  # Since we use filtered_matrix which ONLY has reads from tissues, 
  # filter.matrix = TRUE or FALSE didnt show any difference in sample.seurat
  
  # [i think wrong] If filter.matrix is set to FALSE, all 4992 spots will be recorded in 
  # sample.seurat@images$B8@coordinates. Else, only spots that are over tissue,
  # will be recorded in sample.seurat@images$B8@coordinates.
  
  # Unlike Seurat::Read10X(), Seurat::Load10X_Spatial doesnt have ability to 
  # specify project parameter. So, we manually do it.
  sample.seurat@meta.data <- sample.seurat@meta.data %>% 
    dplyr::mutate(orig.ident = i)
  
  # Assign the seurat object to its corresponding variable
  assign(paste0(i, ".filt"), sample.seurat)
  cat("DATA IMPORTED FOR ", i, ".filt dataset\n")
}

#******************************************************************************#
#                           STEP 2: QUALITY CONTROL                            #
#******************************************************************************#

#**********************STEP 2A: CALCULATE ALL QC METRICS***********************#

# NOTE:  We are going to do the same QC on each sample. So, we can merge the 
# individual seurat objects into a single seurat object and calculate metrics on
# the merged seurat object. However, if the number of cells exceeds 1.5 million,
# then this step will fail as there are simply too many cells. So, it is better
# to calculate QC metrics on individual samples, perform filtering on individual
# samples and then merge the filtered seurat objects which have fewer cells.

# Initialize an empty dataframe where the class of each column resembles those 
# of raw_metadata. We will use this dataframe for making QC plots later.
raw_metadata <- data.frame(Cell = c(""), 
                           Sample = as.factor(1), 
                           nUMIs = c(0), 
                           nGenes = c(0), 
                           MitoRatio = c(0), 
                           RiboRatio = c(0), 
                           Novelty = c(0))

# Calculate QC metrics for each sample individually using raw matrices
for(i in paste0(samples, ".filt")){
  
  sample.seurat <- get(i)
  
  # NOTE: PercentageFeatureSet() ONLY calculates metrics for cells within the 
  # indicated assay. If no assay is indicated, it will use default assay. So, we
  # loop through each assay and calculate metrics
  
  # Get list of assays
  assay.list <- Assays(sample.seurat)
  
  for (assay in assay.list){
    # Compute percent mito percent
    sample.seurat <- Seurat::PercentageFeatureSet(object = sample.seurat,
                                                  pattern = "^[Mm][Tt]-",
                                                  features = NULL,
                                                  col.name = "MitoPercent",
                                                  assay = assay)
    
    # Compute percent ribo percent
    sample.seurat <- Seurat::PercentageFeatureSet(object = sample.seurat,
                                                  pattern = "^[Rr][Pp][SsLl]", 
                                                  features = NULL,
                                                  col.name = "RiboPercent",
                                                  assay = assay)
  }
  
  # Extract metadata
  sample_metadata <- sample.seurat@meta.data
  
  # Rename columns to be more intuitive and add the additional QC metrics:
  # (i)     Cell      : unique identifiers corresponding to each cell i.e. barcodes
  # (ii)    Sample    : sample names
  # (iii)   nUMIs     : number of transcripts per cell
  # (iv)    nGenes    : number of genes per cell
  # (v)     nHTO_UMIs : number of HTO reads per cell
  # (vi)    nHTOs     : number of HTO types per cell
  # (vii)   MitoRatio : MitoPercent/100
  # (viii)	RiboRatio : RiboPercent/100  
  # (ix)    Novelty   : log ratio of genes per UMI
  sample_metadata <- sample_metadata %>% 
    dplyr::mutate(Cell = paste0(orig.ident, "_", rownames(sample_metadata)),
                  Sample = orig.ident,
                  nUMIs = dplyr::case_when(!is.na(nFeature_Spatial.002um) & !is.na(nCount_Spatial.002um) ~ nCount_Spatial.002um,
                                           !is.na(nFeature_Spatial.008um) & !is.na(nCount_Spatial.008um) ~ nCount_Spatial.008um,
                                           !is.na(nFeature_Spatial.016um) & !is.na(nCount_Spatial.016um) ~ nCount_Spatial.016um,
                                           TRUE ~ NA),
                  nGenes = dplyr::case_when(!is.na(nFeature_Spatial.002um) & !is.na(nCount_Spatial.002um) ~ nFeature_Spatial.002um,
                                            !is.na(nFeature_Spatial.008um) & !is.na(nCount_Spatial.008um) ~ nFeature_Spatial.008um,
                                            !is.na(nFeature_Spatial.016um) & !is.na(nCount_Spatial.016um) ~ nFeature_Spatial.016um,
                                            TRUE ~ NA),
                  MitoRatio = MitoPercent/100,
                  RiboRatio = RiboPercent/100,
                  Novelty = log10(nGenes)/log10(nUMIs)) %>%
    dplyr::select(Cell, Sample, contains(c("nFeature", "nCount")), nUMIs, nGenes, MitoRatio, RiboRatio, Novelty)
  
  # Replace the metadata in raw Seurat object
  sample.seurat@meta.data <- sample_metadata
  
  # Append raw metadata of each seurat object which will be used for QC plots later
  raw_metadata <- dplyr::bind_rows(raw_metadata, sample_metadata)
  
  # Assign the seurat object to its corresponding variable
  assign(i, sample.seurat)
}

#******************************STEP 2B: PERFORM QC*****************************#

# In single cell analysis, we have to perform QC to remove bad cells in order to
# accurately annotate each cell. In spatial analysis, RNA from tissue is 
# captured within the spot below the tissue. Since sequencing depth is usually 
# way lower than single cell, we use almost all reads (lenient cutoffs) captured
# in the spot to annotate the spots.

# NOTE: Unlike PercentageFeatureSet(), subset() uses all the cells irrespective 
# of the default assay. 2um bins are too small to contain 100 UMIs or 50 genes 
# and the entire Spatial.002um assay will be removed. So, we use lenient cutoffs

# Perform QC for each sample individually
for(i in paste0(samples, ".filt")){
  
  sample.seurat <- get(i)
  
  # If default assay is 2um and no cells of default assay pass QC, then Seurat 
  # will give error. So, change default assay to 8um
  DefaultAssay(sample.seurat) <- "Spatial.008um"
  
  gene_cutoff <- 1           # reduced to 1 from 250 used for spatial
  umi_cutoff <- 1            # reduced to 1 from 500 used for spatial
  #mito_cutoff <- 0.2
  #ribo_cutoff <- 0.05
  #novelty_cutoff <- 0.8  	    # use 0.8 as starting point. Maximum 0.9
  
  sample.seurat <- base::subset(x = sample.seurat,
                                subset = (nGenes >= gene_cutoff) &
                                  (nUMIs >= umi_cutoff))
  #(MitoRatio <= mito_cutoff) &
  # (RiboRatio >= ribo_cutoff) &
  # (Novelty >= novelty_cutoff))
  
  # Assign the seurat object to its corresponding variable
  assign(i, sample.seurat)
}

# Create a merged Seurat object
# NOTE: Samples will have same barcodes. To keep track of cell identities 
# (i.e. barcodes) coming from each sample after merging, we add a prefix 
# (i.e. sample name) to each barcode using "add.cell.ids"

# IMPORTANT: We DO NOT merge samples in spatial experiments. We analyze each
# sample individually.

# NOTE: spatial objects have multiple assays. add.cell.ids modifies rownames
# ONLY for default assay. So, sample name is not added to cells from other
# assays. This creates problem in identifying common_bc. So, we exclude
# add.cell.ids

#****************************STEP 2C: SAVE THE DATA****************************#

for(i in paste0(samples, ".filt")){
  
  sample.seurat <- get(i)
  saveRDS(sample.seurat, file=paste0(seurat_results, i,".rds"))
}

#******************************************************************************#
#                       STEP 3: RUN THE STANDARD PIPELINE                      #
#******************************************************************************#

# Create workbook to save markers
wb <- openxlsx::createWorkbook()
for (i in samples){
  
  object <- get(paste0(i, ".filt"))
  
  # Run workflow on all assays
  for (assay in c("Spatial.002um", "Spatial.008um", "Spatial.016um")){
    
    DefaultAssay(object) <- assay
    object <- Seurat::NormalizeData(       object = object, assay = assay, normalization.method = "LogNormalize")
    object <- Seurat::FindVariableFeatures(object = object, assay = assay, nfeatures = 2000)
    object <- Seurat::ScaleData(           object = object, assay = assay, features = VariableFeatures(object))
    
    # Create a new 'sketch' assay using 50k cells
    sample_n <- 50000
    e <- "Error"
    while (class(e) == "character"){
      e <- tryCatch(Seurat::SketchData(object = object, assay = assay,
                                       ncells = sample_n,
                                       method = "LeverageScore",
                                       features = VariableFeatures(object),
                                       sketched.assay = paste0(assay, ".sketch")), error = function(msg){
                                         print("Reducing number of cells sampled")
                                         return("Error")})
      sample_n <- sample_n-5000
    }
    
    object <- e
    
    # Switch analysis to sketched cells
    DefaultAssay(object) <- paste0(assay, ".sketch")
    object <- Seurat::FindVariableFeatures(object = object, 
                                           assay = paste0(assay, ".sketch"), 
                                           nfeatures = 2000)
    
    object <- Seurat::ScaleData(object = object, 
                                assay = paste0(assay, ".sketch"), 
                                features = VariableFeatures(object))
    
    object <- Seurat::RunPCA(object = object, 
                             assay = paste0(assay, ".sketch"), 
                             reduction.name = paste0(assay, ".pca.sketch"))
    
    object <- Seurat::FindNeighbors(object = object, 
                                    assay = paste0(assay, ".sketch"), 
                                    reduction = paste0(assay, ".pca.sketch"), 
                                    dims = 1:50)
    
    object <- Seurat::FindClusters(object = object, 
                                   algorithm = 4,
                                   cluster.name = paste0(assay, ".cluster.sketch"), 
                                   resolution = 0.6)
    
    object <- Seurat::RunUMAP(object = object, 
                              reduction.name = paste0(assay, ".umap.sketch"), 
                              reduction = paste0(assay, ".pca.sketch"), 
                              dims = 1:50)
    
    # Project the cluster labels, dimensional reductions (PCA and UMAP) that we
    # learned from the 50,000 sketched cells to the entire dataset
    object <- Seurat::ProjectData(object = object, 
                                  assay = assay,
                                  sketched.assay = paste0(assay, ".sketch"),
                                  sketched.reduction = paste0(assay, ".pca.sketch"),
                                  full.reduction = paste0(assay, ".pca.full"),
                                  dims = 1:50,
                                  refdata = list(seurat_cluster.projected = paste0(assay, ".cluster.sketch")),
                                  umap.model = paste0(assay, ".umap.sketch"))
    
    object <- Seurat::RunUMAP(object = object, 
                              reduction.name = paste0(assay, ".umap.full"), 
                              reduction = paste0(assay, ".pca.full"), 
                              dims = 1:50)
    
    # Append assay name to the newly generated metadata columns so they dont get overwritten
    object@meta.data <- object@meta.data %>%
      dplyr::rename(!!rlang::sym(paste0(assay, ".cluster.full")) := "seurat_cluster.projected",
                    !!rlang::sym(paste0(assay, ".cluster.full.score")) := "seurat_cluster.projected.score")
    
    # We can visualize the clustering results for the sketched cells, as well as the
    # projected clustering results for the full dataset:
    DefaultAssay(object) <- paste0(assay, ".sketch")
    Idents(object) <- paste0(assay, ".cluster.sketch")
    p1 <- DimPlot(object, reduction = paste0(assay, ".umap.sketch"), label = T) + 
      ggtitle("Sketched clustering") + 
      theme(legend.position = "none")
    
    # switch to full dataset
    DefaultAssay(object) <- assay
    Idents(object) <- paste0(assay, ".cluster.full")
    p2 <- DimPlot(object, reduction = paste0(assay, ".umap.full"), label = T) + 
      ggtitle("Projected clustering (full dataset)") + 
      theme(legend.position = "none")
    
    # Save the plot
    ggplot2::ggsave(filename=paste0(assay, ".umap.tiff"),
                    plot=p1+p2,
                    device="jpeg",
                    path=diagnostics_path,
                    width=11,
                    height=8.5,
                    units=c("in"),
                    dpi=600,
                    limitsize=TRUE,
                    bg="white")
    
    # Find markers
    markers <- Seurat::FindAllMarkers(object = object, 
                                      assay = assay,
                                      only.pos = TRUE)
    
    top10 <- markers %>%
      dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
      dplyr::group_by(cluster) %>%
      dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
      slice_max(order_by=avg_log2FC*pct.1, n = 10) %>%
      ungroup() %>%
      data.frame() %>%
      dplyr::arrange(cluster, desc(avg_log2FC*pct.1))
    
    # Save markers to worksheet
    openxlsx::addWorksheet(wb, sheetName = paste0(i, assay))
    openxlsx::writeData(wb, sheet = paste0(i, assay), x = markers, rowNames = FALSE)
    openxlsx::addWorksheet(wb, sheetName = paste0(i, assay, "top10"))
    openxlsx::writeData(wb, sheet = paste0(i, assay, "top10"), x = top10, rowNames = FALSE)
    
  }
  
  # Save xlsx file with markers
  openxlsx::saveWorkbook(wb, file = paste0(seurat_results, i, "_Markers.xlsx"), overwrite = TRUE)
  
  # Save the object
  saveRDS(object, file=paste0(seurat_results, i,".rds"))
}
  
  # DefaultAssay(object) <- "Spatial.008um"
  # 
  # vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + 
  #   theme(axis.text = element_text(size = 4)) + 
  #   NoLegend()
  # count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + 
  #   theme(legend.position = "right")
  
  # # switch back to 8um
  # #object@assays$Spatial.008um@features@dimnames[[1]][25:30]
  # DefaultAssay(object) <- "Spatial.008um"
  # p2 <- SpatialFeaturePlot(object, features = "GAPDH") + ggtitle("Hpca expression (8um)")
  
  
  #   # Visualize the unsupervised clusters based on their spatial location. 
  #   SpatialDimPlot(object, label = T, repel = T, label.size = 4)
  #   
  #   # Plot the spatial location of different clusters individually.
  #   Idents(object) <- "seurat_cluster.projected"
  #   cells <- CellsByIdentities(object, idents = c(0, 4, 32, 34, 35))
  #   p <- SpatialDimPlot(object,
  #                       cells.highlight = cells[setdiff(names(cells), "NA")],
  #                       cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
  #   ) + NoLegend()
  #   p
  #   
  #   
  #   # visualize the top gene expression markers for each cluster:
  #   # Create downsampled object to make visualization either
  #   DefaultAssay(object) <- "Spatial.008um"
  #   Idents(object) <- "seurat_cluster.projected"
  #   object_subset <- subset(object, cells = Cells(object[["Spatial.008um"]]), downsample = 1000)
  #   
  #   
  #   
  #   object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
  #   p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
  #   p
  #   
  #   # Identifying spatially-defined tissue domains
  #   object <- Banksy::RunBanksy(object,
  #                               lambda = 0.8, verbose = TRUE,
  #                               assay = "Spatial.008um", slot = "data", features = "variable",
  #                               k_geom = 50)
  #   
  #   DefaultAssay(object) <- "BANKSY"
  #   object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(object), npcs = 30)
  #   object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:30)
  #   object <- FindClusters(object, cluster.name = "banksy_cluster", resolution = 0.5)
  #   
  #   Idents(object) <- "banksy_cluster"
  #   p <- SpatialDimPlot(object, group.by = "banksy_cluster", label = T, repel = T, label.size = 4)
  #   p
  #   
  #   # highlight the spatial location of each tissue domain individually
  #   banksy_cells <- CellsByIdentities(object)
  #   p <- SpatialDimPlot(object, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()
  #   p
  #   
  # }
  # 

  # 
  # # Load rds file of seurat objects
  # for (i in samples){
  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  #   sample.seurat <- sctransform_spatial_data(sample.seurat)
  #   sample.seurat <- cluster_spatial_data(sample.seurat)
  #   saveRDS(sample.seurat, file=paste0(seurat_results, i,".rds"))
  #   
  #   DefaultAssay(sample.seurat) <- "SCT"
  #   p1 <- DimPlot(sample.seurat, reduction = "umap", label = TRUE)
  #   p2 <- SpatialDimPlot(sample.seurat, label = TRUE, label.size = 3)
  #   p <- p1 + p2
  #   ggsave(filename = paste0(seurat_results, "Clusters_on_slide_", i, ".jpg"), 
  #          plot = p)
  # }
  # 
  # # # There is no integration like scRNA Seq. We analyse each slide individually.
  # # # Unfortunately, all images are stored within each sample. So, we remove 
  # # # unwanted images from each sample. If more than 1 image is present in each 
  # # # sample, SpatialDimPlot() will give error.
  # # integ_data <- sct_data
  # # for (i in 1:length(sct_data)){
  # #   integ_data[[i]] <- cluster_data(sct_data[[i]])
  # #   integ_data[[i]]@images <- integ_data[[i]]@images[names(integ_data[[i]]@images) == names(sct_data)[[i]]]
  # # }
  # 
  # # Color spots were GOI are present based on expression
  # for (i in samples){
  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  #   DefaultAssay(sample.seurat) <- "SCT"
  #   Seurat::SpatialFeaturePlot(object = sample.seurat, 
  #                              features = c("CD8A", "CD8B", "NPEPPS", "CDH12"),
  #                              ncol = 4,
  #                              slot = "data")
  #   ggsave(filename = paste0(seurat_results, "Feature_plot_", i, ".jpg"),
  #          plot = last_plot(),
  #          units = c("in"),
  #          width = 11,
  #          height = 8)
  #   
  #   # This can ONLY plot 2 genes at a time
  #   SpatialFeaturePlotBlend(sample.seurat, "CD8A", "NPEPPS")
  # }
  # 
  # # Color spots were CD8 and NPEPPS are expressed in same plot. 
  # # NOTE: This is NOT expression based. We just color the cells that express our 
  # # GOI in different colors.
  # # Use SCT assay to get counts, not Spatial assay as it has too much background.
  # for (i in samples){
  #   
  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  #   
  #   # Make sure the genes you want are present in the assay
  #   GOI <- intersect(c("CD8A", "CD8B"), rownames(sample.seurat@assays$SCT@data))
  #   cd8_df <- sample.seurat@assays$SCT$counts[GOI, ]
  #   
  #   if (length(GOI) > 1){
  #     cd8_cells <- colnames(cd8_df[,colSums(cd8_df) > 0])
  #   } else {
  #     cd8_cells <- names(cd8_df[(cd8_df>0)])
  #   }
  #   
  #   GOI <- intersect(c("NPEPPS"), rownames(sample.seurat@assays$SCT@data))
  #   npepps_df <- sample.seurat@assays$SCT$counts[c("NPEPPS"), ]
  #   npepps_cells <- names(npepps_df[(npepps_df>0)])
  #   
  #   SpatialPlot(object = sample.seurat, 
  #               cells.highlight = list(CD8 = cd8_cells, NPEPPS = npepps_cells), 
  #               cols.highlight =  c("green", "blue", "grey"),
  #               pt.size.factor = 2)
  #   
  #   ggsave(filename = paste0(seurat_results, "Location_plot_", i, ".jpg"),
  #          plot = last_plot(),
  #          units = c("in"),
  #          width = 11,
  #          height = 8)
  # }
  # 
  # # # Classify each spot as CD8+ or CD8-
  # # #for (i in samples){
  # # i <- "B8"
  # #   
  # #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  # #   #rownames(sample.seurat@meta.data) <- sample.seurat@meta.data$Cell
  # #   
  # #   df <- sample.seurat@assays$Spatial$counts[c("CD8A", "CD8B"), ]
  # #   pos_cells  <- colnames(df[,colSums(df) > 0])
  # #   neg_cells <- setdiff(colnames(df), pos_cells)
  # #   pos_cells <- paste0(i,"_", pos_cells)
  # #   neg_cells <- paste0(i,"_", neg_cells)
  # #   sample.seurat@meta.data <- sample.seurat@meta.data %>% 
  # #     dplyr::mutate(CD8_status = dplyr::case_when(Cell %in% pos_cells ~ "CD8_pos",
  # #                                                 TRUE ~ "CD8_neg"))
  # # }
  # 
  # 
  # i <- "B8"
  # sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  # 
  # a <- GetTissueCoordinates(sample.seurat)
  # df <- data.frame(sample.seurat@images$B8@coordinates)
  # 
  # ### CHECK THAT WE ARE SUCCESSFULLY ABLE TO IDENTIFY NEIGHBORS ###
  # # Find a cell in the center of the tissue
  # row_cell <- median(df$row)
  # col_cell <- median(df$col)
  # cell <- rownames(df %>% dplyr::filter(row==row_cell, col==col_cell))
  # 
  # # Find adjacent cells
  # cell1 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell-2)))
  # cell2 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell+2)))
  # cell3 <- rownames(df %>% dplyr::filter(row==(row_cell-2), col==col_cell))
  # cell4 <- rownames(df %>% dplyr::filter(row==(row_cell+2), col==col_cell))
  # cell5 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell-1)))
  # cell6 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell+1)))
  # cell7 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell-1)))
  # cell8 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell+1)))
  # 
  # SpatialPlot(object = sample.seurat, 
  #             cells.highlight = list(cell, unlist(lapply(paste0("cell", seq(1,8)), get), use.names=FALSE)), 
  #             #facet.highlight = TRUE,
  #             cols.highlight =  c("yellow", "red", "black"), 
  #             pt.size.factor = 2)
  # 
  # #******************************************************************************#
  # 
  # i <- "B8"
  # sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  # df <- data.frame(sample.seurat@images$B8@coordinates)
  # expr_df <- sample.seurat@assays$SCT$data[rownames(sample.seurat@assays$SCT$data) %in% c("NPEPPS", "CD8A", "CD8B"),]
  # 
  # # Identify cells that express NPEPPS
  # npepps <- expr_df["NPEPPS",]
  # npepps <- npepps[npepps > 0]
  # npepps <- names(npepps)
  # # Not all npepps cells have co-ordinates. Remove those that lack co-ordinates.
  # npepps <- intersect(npepps, rownames(df))
  # 
  # expr_df <- data.frame(expr_df)
  # # For every cell that expresses NPEPPS, calculate total expression of CD8A and 
  # # CD8B in 1st neighbor
  # t_expr <- c()
  # npepps_expr <- c()
  # for (i in npepps){
  #   row_cell <- df[rownames(df) == i,]$row
  #   col_cell <- df[rownames(df) == i,]$col
  #   
  #   # Find adjacent cells
  #   cell1 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell-2)))
  #   cell2 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell+2)))
  #   cell3 <- rownames(df %>% dplyr::filter(row==(row_cell-2), col==col_cell))
  #   cell4 <- rownames(df %>% dplyr::filter(row==(row_cell+2), col==col_cell))
  #   cell5 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell-1)))
  #   cell6 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell+1)))
  #   cell7 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell-1)))
  #   cell8 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell+1)))
  #   
  #   # SpatialPlot(object = sample.seurat, 
  #   #             cells.highlight = list(i, unlist(lapply(paste0("cell", seq(1,8)), get), use.names=FALSE)), 
  #   #             #facet.highlight = TRUE,
  #   #             cols.highlight =  c("yellow", "red", "black"), 
  #   #             pt.size.factor = 2)
  #   
  #   cells <- intersect(make.names(unlist(lapply(paste0("cell", seq(1,8)), get))), colnames(expr_df))
  #   print(cells)
  #   
  #   t <- expr_df %>% 
  #     dplyr::select(all_of(cells)) %>%
  #     tibble::rownames_to_column("Gene") %>%
  #     dplyr::filter(Gene != "NPEPPS") %>%
  #     tibble::column_to_rownames("Gene")
  #   
  #   n <- expr_df %>% 
  #     dplyr::select(all_of(cells)) %>%
  #     tibble::rownames_to_column("Gene") %>%
  #     dplyr::filter(Gene == "NPEPPS") %>%
  #     tibble::column_to_rownames("Gene")
  #   
  #   t_expr <- c(t_expr, sum(t))
  #   npepps_expr <- c(npepps_expr, sum(n))
  # }
  # 
  # df <- data.frame(npepps_expr, t_expr)
  # # Save batch corrected normalized counts for entire dataset
  # wb <- openxlsx::createWorkbook()
  # openxlsx::addWorksheet(wb, sheetName = "spatial")
  # openxlsx::writeData(wb, sheet = "spatial", x = df, rowNames = FALSE)
  # openxlsx::saveWorkbook(wb,
  #                        file = paste0(seurat_results, "Correlation.xlsx"),
  #                        overwrite = TRUE)
  # 

  # 
  # # pseudobulk the counts based on donor-condition-celltype
  # pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("stim", "donor_id", "seurat_annotations"))
  # 
  # # each 'cell' is a donor-condition-celltype pseudobulk profile
  # tail(Cells(pseudo_ifnb))
  # 
  # # Identification of Spatially Variable Features
  # 
  # # Seurat offers two workflows to identify molecular features that correlate 
  # # with spatial location within a tissue. The first is to perform differential 
  # # expression based on pre-annotated anatomical regions within the tissue, which
  # # may be determined either from unsupervised clustering or prior knowledge. 
  # # This strategy works will in this case, as the clusters above exhibit clear 
  # # spatial restriction.
  # 
  # de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
  # SpatialFeaturePlot(object = brain, 
  #                    features = rownames(de_markers)[1:3], 
  #                    alpha = c(0.1, 1), 
  #                    ncol = 3)
  # 
  # 
  # # An alternative approach, implemented in FindSpatiallyVariables(), is to 
  # # search for features exhibiting spatial patterning in the absence of 
  # # pre-annotation. The default method (method = 'markvariogram), is inspired by 
  # # the Trendsceek, which models spatial transcriptomics data as a mark point 
  # # process and computes a ‘variogram’, which identifies genes whose expression 
  # # level is dependent on their spatial location. More specifically, this process
  # # calculates gamma(r) values measuring the dependence between two spots a 
  # # certain “r” distance apart. By default, we use an r-value of ‘5’ in these 
  # # analyses, and only compute these values for variable genes (where variation 
  # # is calculated independently of spatial location) to save time.
  # 
  # 
  # brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
  #                                        selection.method = "moransi")
  # top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
  # SpatialFeaturePlot(brain, 
  #                    features = top.features, 
  #                    alpha = c(0.1, 1),
  #                    ncol = 3)
  # 
  # 
  # 
  # # Integration with single-cell data i.e. label transfer from single-cell data
  # # NOTE: At ~50um size, spots from the visium assay will encompass the expression 
  # # profiles of multiple cells. 
  # 
  # # Users may be interested to ‘deconvolute’ each of the spatial voxels to predict
  # # the underlying composition of cell types. We tested a wide variety of 
  # # deconvolution and integration methods, using a reference scRNA-seq dataset of
  # # ~14,000 adult mouse cortical cell taxonomy from the Allen Institute, 
  # # generated with the SMART-Seq2 protocol. 
  # 
  # # We consistently found superior performance using integration methods (as 
  # # opposed to deconvolution methods), likely because of substantially different 
  # # noise models that characterize spatial and single-cell datasets, and 
  # # integration methods are specifically designed to be robust to 
  # # these differences. 
  # 
  # # We therefore apply the ‘anchor’-based integration workflow introduced in 
  # # Seurat v3, that enables the probabilistic transfer of annotations from a 
  # # reference to a query set. 
  # 
  # # While many of the methods are conserved (both procedures begin by identifying 
  # # anchors), there are two important distinctions between data transfer and 
  # # integration:
  # # (i) In data transfer, Seurat does not correct or modify query expression data.
  # # (ii) In data transfer, Seurat has an option (set by default) to project the 
  # # PCA structure of a reference onto the query, instead of learning a joint 
  # # structure with CCA. We generally suggest using this option when projecting 
  # # data between scRNA-seq datasets.
  # 
  # # After finding anchors, we use the TransferData() function to classify the 
  # # query cells based on reference data (a vector of reference cell type labels). 
  # # TransferData() returns a matrix with predicted IDs and prediction scores, 
  # # which we can add to the query metadata.
  # 
  # # NOTE: Make sure reference seurat object is SCTransformed.
  # integrated_seurat <- readRDS("/hpc/home/kailasamms/scratch/scRNASeq_Chen/results_seurat/integrated_seurat_snn.rds")
  # 
  # # Seurat v3 vs Seurat v5 issues: Check if cells in graph are in same order
  # all(Cells(integrated_seurat@graphs$integrated_snn) == colnames(integrated_seurat))
  # # if FALSE, reorder the cells in all existing graphs
  # integrated_seurat@graphs$integrated_snn <- 
  #   integrated_seurat@graphs$integrated_snn[colnames(integrated_seurat), colnames(integrated_seurat)]
  # integrated_seurat@graphs$integrated_nn <- 
  #   integrated_seurat@graphs$integrated_nn[colnames(integrated_seurat), colnames(integrated_seurat)]
  # 
  # # Remove ambiguous cells before label transfer
  # integrated_seurat <- subset(integrated_seurat, 
  #                             sub_type == "Unclassified", 
  #                             invert=TRUE)
  # 
  # # NOTE: For FindTransferAnchors() to work, both reference and query MUST be
  # # SCtransformed()  https://github.com/satijalab/seurat/issues/3937
  # # So, run SCTransform() on RNA assay.
  # DefaultAssay(integrated_seurat) <- "integrated"
  # 
  # # Make sure all assays in Seuratv3 object has same order of cells
  # acells <- colnames(x = integrated_seurat[["integrated"]])
  # ocells <- colnames(x = integrated_seurat)
  # 
  # integrated_seurat@assays$RNA@counts <- integrated_seurat@assays$RNA@counts[,ocells]
  # integrated_seurat@assays$RNA@data <- integrated_seurat@assays$RNA@data[,ocells]
  # integrated_seurat@assays$integrated@data <- integrated_seurat@assays$integrated@data[,ocells]
  # integrated_seurat@assays$integrated@scale.data <- integrated_seurat@assays$integrated@scale.data[,ocells]
  # integrated_seurat@assays$SCT@counts <- integrated_seurat@assays$SCT@counts[,ocells]
  # integrated_seurat@assays$SCT@data <- integrated_seurat@assays$SCT@data[,ocells]
  # integrated_seurat@assays$SCT@scale.data <- integrated_seurat@assays$SCT@scale.data[,ocells]
  # 
  # # NOTE: If you subset the spatial object, re-run SCTranform() on it.
  # for (i in samples){
  #   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  #   DefaultAssay(sample.seurat) <- "SCT"
  #   
  #   # Run ONLY if you had subset the spatial data
  #   # sample.seurat <- SCTransform(sample.seurat, 
  #   #                              assay = "Spatial",
  #   #                              verbose = FALSE)
  #   # sample.seurat <- RunPCA(sample.seurat,
  #   #                         verbose = FALSE)
  #   
  #   # Find anchors between reference and query
  #   anchors <- FindTransferAnchors(reference = integrated_seurat,
  #                                  query = sample.seurat,
  #                                  normalization.method = "SCT",
  #                                  reference.assay = "integrated",
  #                                  query.assay = "SCT")
  #   
  #   predictions.assay <- TransferData(anchorset = anchors, 
  #                                     refdata = integrated_seurat$cell_type,
  #                                     prediction.assay = TRUE,
  #                                     weight.reduction = sample.seurat[["pca"]],
  #                                     dims = 1:30)
  #   sample.seurat[["predictions"]] <- predictions.assay
  #   
  #   DefaultAssay(sample.seurat) <- "predictions"
  #   
  #   SpatialFeaturePlot(sample.seurat, 
  #                      features = rownames(predictions.assay@data), 
  #                      pt.size.factor = 1.6,
  #                      #ncol = 4,
  #                      crop = FALSE)
  #   
  #   ggsave(filename = paste0(seurat_results, "Labeltransfer_plot_", i, ".jpg"),
  #          plot = last_plot(),
  #          units = c("in"),
  #          width = 11,
  #          height = 8)
  #   
  #   
  #   
  #   # sample.seurat <- FindSpatiallyVariableFeatures(sample.seurat,
  #   #                                                assay = "predictions",
  #   #                                                selection.method = "moransi",
  #   #                                                features = rownames(sample.seurat),
  #   #                                                r.metric = 5,
  #   #                                                slot = "data")
  #   
  #   # top.clusters <- head(SpatiallyVariableFeatures(sample.seurat, 
  #   #                                                selection.method = "moransi"), 4)
  #   # SpatialPlot(object = sample.seurat, 
  #   #             features = top.clusters,
  #   #             ncol = 2)
  # }