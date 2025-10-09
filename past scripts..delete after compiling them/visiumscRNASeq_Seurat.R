#!/usr/bin/env Rscript

# NOTE: All variables and functions are defined within the file below

source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")

# NOTE: In visium v1 spatial transcriptomics, there are 4992 barcoded spots 
# which capture RNA from tissue placed above these spots. 
# So, unlike scRNASeq, each barcode corresponds to a spot, not a cell.
# https://www.youtube.com/watch?v=VwNk4d-0RJc
# https://kb.10xgenomics.com/hc/en-us/articles/360035848191-How-many-spots-are-
# within-a-single-capture-area-on-the-Visium-v1-Spatial-Gene-Expression-Slide

#******************************************************************************#
#                       STEP 1: SETUP THE SEURAT OBJECT                        #
#******************************************************************************#

#************************IMPORTING DATA FROM h5AD FILE*************************#

# data.dir MUST have a folder named "spatial" with the image and H5 file 
# specified by filename parameter.

# Create a list of samples which will be added to each barcode.
# Since folder names correspond to sample name, we just use list.files()
samples <- list.files(path = feature_matrix_path)

# Loop through each of the individual folders in parent directory & import data
for(i in samples){
  
  sample.seurat <- Load10X_Spatial(data.dir = paste0(feature_matrix_path, i),
                                   filename = "raw_feature_bc_matrix.h5",
                                   assay = "Spatial",
                                   slice = i,
                                   filter.matrix = TRUE,
                                   to.upper = FALSE,
                                   image = NULL)
  # if filter.matrix is set to FALSE, all 4992 spots will be recorded in 
  #sample.seurat@images$B8@coordinates. Else, only spots that are over tissue,
  # will be recorded in sample.seurat@images$B8@coordinates.
  
  # Unlike Seurat::Read10X(), Seurat::Load10X_Spatial doesnt have ability to 
  # specify project parameter. So, we manually do it.
  sample.seurat@meta.data <- sample.seurat@meta.data %>% 
    dplyr::mutate(orig.ident = i)
  
  #print(dim(sample.seurat@meta.data))
  
  # Assign the seurat object to its corresponding variable
  assign(i, sample.seurat)
  
  # Explore the meta.data slot
  cat("\nFirst few rows of ", i, "\n")
  print(head(sample.seurat@meta.data))
  cat("\nLast few rows of ", i, "\n")
  print(tail(sample.seurat@meta.data))
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
# samples and then merge the filtered seurat objects.

# Initialize an empty dataframe where the class of each column resembles those 
# of raw_metadata. We will use this dataframe for making plots later.
raw_metadata <- data.frame(Cell = c(""), Sample = as.factor(1), 
                           nUMIs = c(0), nGenes = c(0), 
                           MitoRatio = c(0), RiboRatio = c(0), Novelty = c(0))

# Calculate QC metrics for each sample individually  
for(i in samples){
  
  sample.seurat <- get(i)
  
  # Compute percent mito percent
  sample.seurat <- Seurat::PercentageFeatureSet(object = sample.seurat,
                                                pattern = dplyr::if_else(species=="Homo sapiens", "MT-", "mt-"),
                                                features = NULL,
                                                col.name = "MitoPercent",
                                                assay = NULL)
  
  # Compute percent ribo percent
  sample.seurat <- Seurat::PercentageFeatureSet(object = sample.seurat,
                                                pattern = dplyr::if_else(species=="Homo sapiens", "^RP[SL]", "^Rp[sl]"),
                                                features = NULL,
                                                col.name = "RiboPercent",
                                                assay = NULL)
  
  # Extract metadata
  sample_metadata <- sample.seurat@meta.data
  
  # Rename columns to be more intuitive and add the additional QC metrics:
  # (i)     Cell      : Unique identifiers corresponding to each cell = barcodes
  # (ii)    Sample    : sample name
  # (iii)   nUMIs     : number of transcripts per cell
  # (iv)    nGenes    : number of genes per cell
  # (v)     nHTO_UMIs : number of HTO reads per cell
  # (vi)    nHTOs     : number of HTOs types per cell
  # (vii)   MitoRatio : MitoPercent/100
  # (viii)	RiboRatio : RiboPercent/100  
  # (ix)    Novelty   : log ratio of genes per UMI
  sample_metadata <- sample_metadata %>% 
    dplyr::transmute(Cell = paste0(orig.ident, "_", rownames(sample_metadata)),
                     Sample = orig.ident,
                     nUMIs = nCount_Spatial,
                     nGenes = nFeature_Spatial,
                     MitoRatio = MitoPercent/100,
                     RiboRatio = RiboPercent/100,
                     Novelty = log10(nGenes)/log10(nUMIs))
  
  # Replace the metadata in raw Seurat object
  sample.seurat@meta.data <- sample_metadata
  
  # Save raw metadata
  raw_metadata <- dplyr::bind_rows(raw_metadata, sample_metadata)
  
  # Assign the seurat object to its corresponding variable
  assign(i, sample.seurat)
}

#******************************STEP 2B: PERFORM QC*****************************#

# Perform QC for each sample individually
for(i in samples){
  
  sample.seurat <- get(i)
  
  gene_cutoff <- 50              # reduced to 50 from 250 used for scRNASeq
  umi_cutoff <- 500
  mito_cutoff <- 0.2
  ribo_cutoff <- 0.05
  novelty_cutoff <- 0.8  	        # use 0.8 as starting point. Maximum 0.9
  
  sample.seurat <- base::subset(x = sample.seurat,
                                subset = ((nGenes >= gene_cutoff) &
                                            (nUMIs >= umi_cutoff) &
                                            (MitoRatio <= mito_cutoff) &
                                            # (RiboRatio >= ribo_cutoff) &
                                            # (!is.na(filtered_seurat@meta.data$Patient)) &
                                            (Novelty >= novelty_cutoff)))
  
  # Assign the seurat object to its corresponding variable
  assign(i, sample.seurat)
}

# Create a merged Seurat object.
# NOTE: Samples will have same barcodes. To keep track of cell identities 
# (i.e.barcodes) coming from each sample after merging, we add a prefix 
# (i.e. sample name) to each barcode using add.cell.ids.
filtered_seurat <- base::merge(x = get(samples[1]),
                               y = lapply(samples[2:length(samples)], get),
                               add.cell.ids = samples,
                               merge.data = FALSE)

#****************************STEP 2C: SAVE THE DATA****************************#

# Create .rds object for filtered seurat object to load at any time
saveRDS(filtered_seurat, file=paste0(seurat_results, "filtered_seurat.rds"))

for(i in samples){
  saveRDS(get(i), file=paste0(seurat_results, i,".rds"))
}

#******************************************************************************#

#*****************STEP 2D: VISUALIZE DATA BEFORE AND AFTER QC******************#

# Remove dummy first row from raw_metadata
raw_metadata <- raw_metadata[-1,]
rownames(raw_metadata) <- raw_metadata$Cell

filtered_metadata <- filtered_seurat@meta.data

# Visualize the number of cell counts per sample
cell_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, fill=Sample)) + 
    geom_bar() +              
    theme_classic() +         #display with x and y axis lines and no gridlines
    labs(x = "Sample", y = "Number of Cells", title = stringr::str_wrap(paste0("Number of Cells ", tag), 30)) +
    #coord_cartesian(ylim = c(0, 20000)) +
    geom_text(stat="count", aes(label=after_stat(count)), vjust = -1) +
    my_theme
}

# Visualize the number of UMIs per cell
umi_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, y=nUMIs, fill=Sample)) +
    geom_violin() +        
    theme_classic() +       
    labs(x = "Sample", y = "Number of UMIs", title = stringr::str_wrap(paste0("Distribution of UMIs ", tag),30)) +
    coord_cartesian(ylim = c(100, 1000000)) +
    my_theme +        
    scale_y_log10(breaks = c(100, 1000, 10000, 100000, 1000000)) +  		     #display y axis in log scale
    geom_boxplot(width=0.1) + 
    geom_hline(yintercept = gene_cutoff, linetype = 2)
}

# Visualize the number of genes per cell
gene_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, y=nGenes, fill = Sample)) +
    geom_violin() +         
    theme_classic() +       
    labs(x = "Sample", y = "Number of Genes", title = stringr::str_wrap(paste0("Distribution of Genes ", tag),30)) +
    coord_cartesian(ylim = c(1, 30000)) +
    my_theme +        
    scale_y_log10() +    	
    geom_boxplot(width=0.1) + 
    geom_hline(yintercept = gene_cutoff, linetype = 2)
}

# Visualize the MitoRatio of each cell
mito_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, y=MitoRatio, fill = Sample)) +
    geom_violin() +         
    theme_classic() +       
    labs(x = "Sample", y = "MitoRatio", title = stringr::str_wrap(paste0("Distribution of MitoRatio ", tag),30)) +
    coord_cartesian(ylim = c(0, 1)) +
    my_theme +        
    geom_boxplot(width=0.1) + 
    geom_hline(yintercept = mito_cutoff, linetype = 2)
}

# Visualize the RiboRatio of each cell
ribo_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, y=RiboRatio, fill = Sample)) +
    geom_violin() +         
    theme_classic() +       
    labs(x = "Sample", y = "RiboRatio", title = stringr::str_wrap(paste0("Distribution of RiboRatio ", tag),30)) +
    coord_cartesian(ylim = c(0, 1)) +
    my_theme +        
    geom_boxplot(width=0.1) + 
    geom_hline(yintercept = ribo_cutoff, linetype = 2)
}

# Visualize the novelty or complexity of each cell
novelty_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=Sample, y=Novelty, fill = Sample)) +
    geom_violin() +     
    theme_classic() + 
    labs(x = "Sample", y = "Novelty Score", title = stringr::str_wrap(paste0("Distribution of Novelty Score ", tag),30)) +
    coord_cartesian(ylim = c(0, 1)) +
    my_theme +        
    geom_boxplot(width=0.1) + 
    geom_hline(yintercept = novelty_cutoff, linetype = 2)
}

# Visualize number of genes/cell, number of UMIs/cell & MitoRatio together.
# Bottom left quadrant : Poor quality cells with low genes & UMIs per cell 
# Top right quadrant   : Good quality cells with high genes & UMIs per cell
# Bottom right quadrant: Cells with low genes but high UMIs per cell. These 
# could be dying cells or population of low complexity cells (i.e erythrocytes)
gene_umi_mito_qc <- function(metadata, tag){
  
  ggplot(data = get(metadata), aes(x=nUMIs, y=nGenes, color = MitoRatio)) +
    geom_point() +
    theme_classic() + 
    labs(x = "Number of UMIs", y = "Number of Genes",	 title = paste0("Distribution of UMIs, Genes & MitoRatio ", tag)) +
    my_theme + 
    coord_cartesian(xlim = c(100, 1000000), ylim = c(100, 20000)) +
    scale_x_log10(breaks = c(100, 1000, 10000, 100000, 1000000)) + 
    scale_y_log10() + 
    facet_wrap(.~Sample, nrow = 4) +   #split the plot by X-axis label
    stat_smooth(method=lm, color="yellow") +
    geom_vline(xintercept = umi_cutoff) +    	#draw a vertical line at x=500 i.e.UMIs cutoff
    geom_hline(yintercept = gene_cutoff) +    #draw a horizontal line at y =250 i.e. Genes cutoff
    scale_color_viridis(option = "D", limits = c(0, 1)) 		# limits sets max and min values of gradient 
}

# Plot all QC metrics before and after QC
funcs <- c("cell_qc", "umi_qc", "gene_qc", "mito_qc", "ribo_qc", "novelty_qc",
           "gene_umi_mito_qc")

filenames <- c("Cell_Counts", "UMI_Distribution", "Gene_Distribution",
               "MitoRatio_Distribution", "RiboRatio_Distribution", 
               "Novelty_Score_Distribution", "Genes_UMI_MitoRatio_Distribution")

for (i in 1:length(funcs)){
  
  # Plot QC metrics
  purrr::map2(.x = c("raw_metadata", "filtered_metadata"),
              .y = c("Pre QC", "Post QC"),
              .f = get(funcs[i])) %>% 
    cowplot::plot_grid(plotlist = .,
                       align = "hv",
                       axis = "tblr",
                       nrow = 2,  
                       ncol = 1, 
                       rel_widths = 1,
                       rel_heights = 1,
                       labels = NULL,
                       label_size = 14,
                       label_fontfamily = NULL,
                       label_fontface = "bold",
                       label_colour = NULL,
                       label_x = 0,
                       label_y = 1,
                       hjust = -0.5,
                       vjust = 1.5,
                       #scale = 1,
                       greedy = TRUE,
                       byrow = TRUE,
                       cols = NULL,
                       rows = NULL)  
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("QC_", filenames[i], ".pdf"),
                  plot = last_plot(),
                  device = "pdf",
                  path = seurat_results,
                  #scale = 1,
                  #width = dplyr::if_else(i==7, 17,length(samples)/2),
                  height = dplyr::if_else(i==7, 22, 11),
                  units = c("in"),	 
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
}

#******************************************************************************#
#                       STEP 3: RUN THE STANDARD PIPELINE                      #
#******************************************************************************#

# Load rds file of seurat objects
for (i in samples){
  sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  sample.seurat <- sctransform_spatial_data(sample.seurat)
  sample.seurat <- cluster_spatial_data(sample.seurat)
  saveRDS(sample.seurat, file=paste0(seurat_results, i,".rds"))
  
  DefaultAssay(sample.seurat) <- "SCT"
  p1 <- DimPlot(sample.seurat, reduction = "umap", label = TRUE)
  p2 <- SpatialDimPlot(sample.seurat, label = TRUE, label.size = 3)
  p <- p1 + p2
  ggsave(filename = paste0(seurat_results, "Clusters_on_slide_", i, ".jpg"), 
         plot = p)
}

# # There is no integration like scRNA Seq. We analyse each slide individually.
# # Unfortunately, all images are stored within each sample. So, we remove 
# # unwanted images from each sample. If more than 1 image is present in each 
# # sample, SpatialDimPlot() will give error.
# integ_data <- sct_data
# for (i in 1:length(sct_data)){
#   integ_data[[i]] <- cluster_data(sct_data[[i]])
#   integ_data[[i]]@images <- integ_data[[i]]@images[names(integ_data[[i]]@images) == names(sct_data)[[i]]]
# }

# Color spots were GOI are present based on expression
for (i in samples){
  sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  DefaultAssay(sample.seurat) <- "SCT"
  Seurat::SpatialFeaturePlot(object = sample.seurat, 
                             features = c("CD8A", "CD8B", "NPEPPS", "CDH12"),
                             ncol = 4,
                             slot = "data")
  ggsave(filename = paste0(seurat_results, "Feature_plot_", i, ".jpg"),
         plot = last_plot(),
         units = c("in"),
         width = 11,
         height = 8)
  
  # This can ONLY plot 2 genes at a time
  SpatialFeaturePlotBlend(sample.seurat, "CD8A", "NPEPPS")
}

# Color spots were CD8 and NPEPPS are expressed in same plot. 
# NOTE: This is NOT expression based. We just color the cells that express our 
# GOI in different colors.
# Use SCT assay to get counts, not Spatial assay as it has too much background.
for (i in samples){
  
  sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  
  # Make sure the genes you want are present in the assay
  GOI <- intersect(c("CD8A", "CD8B"), rownames(sample.seurat@assays$SCT@data))
  cd8_df <- sample.seurat@assays$SCT$counts[GOI, ]
  
  if (length(GOI) > 1){
    cd8_cells <- colnames(cd8_df[,colSums(cd8_df) > 0])
  } else {
    cd8_cells <- names(cd8_df[(cd8_df>0)])
  }
  
  GOI <- intersect(c("NPEPPS"), rownames(sample.seurat@assays$SCT@data))
  npepps_df <- sample.seurat@assays$SCT$counts[c("NPEPPS"), ]
  npepps_cells <- names(npepps_df[(npepps_df>0)])
  
  SpatialPlot(object = sample.seurat, 
              cells.highlight = list(CD8 = cd8_cells, NPEPPS = npepps_cells), 
              cols.highlight =  c("green", "blue", "grey"),
              pt.size.factor = 2)
  
  ggsave(filename = paste0(seurat_results, "Location_plot_", i, ".jpg"),
         plot = last_plot(),
         units = c("in"),
         width = 11,
         height = 8)
}
  
  




# # Classify each spot as CD8+ or CD8-
# #for (i in samples){
# i <- "B8"
#   
#   sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
#   #rownames(sample.seurat@meta.data) <- sample.seurat@meta.data$Cell
#   
#   df <- sample.seurat@assays$Spatial$counts[c("CD8A", "CD8B"), ]
#   pos_cells  <- colnames(df[,colSums(df) > 0])
#   neg_cells <- setdiff(colnames(df), pos_cells)
#   pos_cells <- paste0(i,"_", pos_cells)
#   neg_cells <- paste0(i,"_", neg_cells)
#   sample.seurat@meta.data <- sample.seurat@meta.data %>% 
#     dplyr::mutate(CD8_status = dplyr::case_when(Cell %in% pos_cells ~ "CD8_pos",
#                                                 TRUE ~ "CD8_neg"))
# }


i <- "B8"
sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))

a <- GetTissueCoordinates(sample.seurat)
df <- data.frame(sample.seurat@images$B8@coordinates)

### CHECK THAT WE ARE SUCCESSFULLY ABLE TO IDENTIFY NEIGHBORS ###
# Find a cell in the center of the tissue
row_cell <- median(df$row)
col_cell <- median(df$col)
cell <- rownames(df %>% dplyr::filter(row==row_cell, col==col_cell))

# Find adjacent cells
cell1 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell-2)))
cell2 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell+2)))
cell3 <- rownames(df %>% dplyr::filter(row==(row_cell-2), col==col_cell))
cell4 <- rownames(df %>% dplyr::filter(row==(row_cell+2), col==col_cell))
cell5 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell-1)))
cell6 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell+1)))
cell7 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell-1)))
cell8 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell+1)))

SpatialPlot(object = sample.seurat, 
            cells.highlight = list(cell, unlist(lapply(paste0("cell", seq(1,8)), get), use.names=FALSE)), 
            #facet.highlight = TRUE,
            cols.highlight =  c("yellow", "red", "black"), 
            pt.size.factor = 2)

#******************************************************************************#

i <- "B8"
sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
df <- data.frame(sample.seurat@images$B8@coordinates)
expr_df <- sample.seurat@assays$SCT$data[rownames(sample.seurat@assays$SCT$data) %in% c("NPEPPS", "CD8A", "CD8B"),]

# Identify cells that express NPEPPS
npepps <- expr_df["NPEPPS",]
npepps <- npepps[npepps > 0]
npepps <- names(npepps)
# Not all npepps cells have co-ordinates. Remove those that lack co-ordinates.
npepps <- intersect(npepps, rownames(df))

expr_df <- data.frame(expr_df)
# For every cell that expresses NPEPPS, calculate total expression of CD8A and 
# CD8B in 1st neighbor
t_expr <- c()
npepps_expr <- c()
for (i in npepps){
  row_cell <- df[rownames(df) == i,]$row
  col_cell <- df[rownames(df) == i,]$col
  
  # Find adjacent cells
  cell1 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell-2)))
  cell2 <- rownames(df %>% dplyr::filter(row==row_cell, col==(col_cell+2)))
  cell3 <- rownames(df %>% dplyr::filter(row==(row_cell-2), col==col_cell))
  cell4 <- rownames(df %>% dplyr::filter(row==(row_cell+2), col==col_cell))
  cell5 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell-1)))
  cell6 <- rownames(df %>% dplyr::filter(row==(row_cell-1), col==(col_cell+1)))
  cell7 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell-1)))
  cell8 <- rownames(df %>% dplyr::filter(row==(row_cell+1), col==(col_cell+1)))
  
  # SpatialPlot(object = sample.seurat, 
  #             cells.highlight = list(i, unlist(lapply(paste0("cell", seq(1,8)), get), use.names=FALSE)), 
  #             #facet.highlight = TRUE,
  #             cols.highlight =  c("yellow", "red", "black"), 
  #             pt.size.factor = 2)
  
  cells <- intersect(make.names(unlist(lapply(paste0("cell", seq(1,8)), get))), colnames(expr_df))
  print(cells)
  
  t <- expr_df %>% 
    dplyr::select(all_of(cells)) %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::filter(Gene != "NPEPPS") %>%
    tibble::column_to_rownames("Gene")
  
  n <- expr_df %>% 
    dplyr::select(all_of(cells)) %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::filter(Gene == "NPEPPS") %>%
    tibble::column_to_rownames("Gene")
  
  t_expr <- c(t_expr, sum(t))
  npepps_expr <- c(npepps_expr, sum(n))
}

df <- data.frame(npepps_expr, t_expr)
# Save batch corrected normalized counts for entire dataset
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "spatial")
openxlsx::writeData(wb, sheet = "spatial", x = df, rowNames = FALSE)
openxlsx::saveWorkbook(wb,
                       file = paste0(seurat_results, "Correlation.xlsx"),
                       overwrite = TRUE)

















# pseudobulk the counts based on donor-condition-celltype
pseudo_ifnb <- AggregateExpression(ifnb, assays = "RNA", return.seurat = T, group.by = c("stim", "donor_id", "seurat_annotations"))

# each 'cell' is a donor-condition-celltype pseudobulk profile
tail(Cells(pseudo_ifnb))




































# Identification of Spatially Variable Features

# Seurat offers two workflows to identify molecular features that correlate 
# with spatial location within a tissue. The first is to perform differential 
# expression based on pre-annotated anatomical regions within the tissue, which
# may be determined either from unsupervised clustering or prior knowledge. 
# This strategy works will in this case, as the clusters above exhibit clear 
# spatial restriction.

de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, 
                   features = rownames(de_markers)[1:3], 
                   alpha = c(0.1, 1), 
                   ncol = 3)


# An alternative approach, implemented in FindSpatiallyVariables(), is to 
# search for features exhibiting spatial patterning in the absence of 
# pre-annotation. The default method (method = 'markvariogram), is inspired by 
# the Trendsceek, which models spatial transcriptomics data as a mark point 
# process and computes a ‘variogram’, which identifies genes whose expression 
# level is dependent on their spatial location. More specifically, this process
# calculates gamma(r) values measuring the dependence between two spots a 
# certain “r” distance apart. By default, we use an r-value of ‘5’ in these 
# analyses, and only compute these values for variable genes (where variation 
# is calculated independently of spatial location) to save time.


brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
                                       selection.method = "moransi")
top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
SpatialFeaturePlot(brain, 
                   features = top.features, 
                   alpha = c(0.1, 1),
                   ncol = 3)



# Integration with single-cell data i.e. label transfer from single-cell data
# NOTE: At ~50um size, spots from the visium assay will encompass the expression 
# profiles of multiple cells. 

# Users may be interested to ‘deconvolute’ each of the spatial voxels to predict
# the underlying composition of cell types. We tested a wide variety of 
# deconvolution and integration methods, using a reference scRNA-seq dataset of
# ~14,000 adult mouse cortical cell taxonomy from the Allen Institute, 
# generated with the SMART-Seq2 protocol. 

# We consistently found superior performance using integration methods (as 
# opposed to deconvolution methods), likely because of substantially different 
# noise models that characterize spatial and single-cell datasets, and 
# integration methods are specifically designed to be robust to 
# these differences. 

# We therefore apply the ‘anchor’-based integration workflow introduced in 
# Seurat v3, that enables the probabilistic transfer of annotations from a 
# reference to a query set. 

# While many of the methods are conserved (both procedures begin by identifying 
# anchors), there are two important distinctions between data transfer and 
# integration:
# (i) In data transfer, Seurat does not correct or modify query expression data.
# (ii) In data transfer, Seurat has an option (set by default) to project the 
# PCA structure of a reference onto the query, instead of learning a joint 
# structure with CCA. We generally suggest using this option when projecting 
# data between scRNA-seq datasets.

# After finding anchors, we use the TransferData() function to classify the 
# query cells based on reference data (a vector of reference cell type labels). 
# TransferData() returns a matrix with predicted IDs and prediction scores, 
# which we can add to the query metadata.

# NOTE: Make sure reference seurat object is SCTransformed.
integrated_seurat <- readRDS("/hpc/home/kailasamms/scratch/scRNASeq_Chen/results_seurat/integrated_seurat_snn.rds")

# Seurat v3 vs Seurat v5 issues: Check if cells in graph are in same order
all(Cells(integrated_seurat@graphs$integrated_snn) == colnames(integrated_seurat))
# if FALSE, reorder the cells in all existing graphs
integrated_seurat@graphs$integrated_snn <- 
  integrated_seurat@graphs$integrated_snn[colnames(integrated_seurat), colnames(integrated_seurat)]
integrated_seurat@graphs$integrated_nn <- 
  integrated_seurat@graphs$integrated_nn[colnames(integrated_seurat), colnames(integrated_seurat)]

# Remove ambiguous cells before label transfer
integrated_seurat <- subset(integrated_seurat, 
                            sub_type == "Unclassified", 
                            invert=TRUE)

# NOTE: For FindTransferAnchors() to work, both reference and query MUST be
# SCtransformed()  https://github.com/satijalab/seurat/issues/3937
# So, run SCTransform() on RNA assay.
DefaultAssay(integrated_seurat) <- "integrated"

# Make sure all assays in Seuratv3 object has same order of cells
acells <- colnames(x = integrated_seurat[["integrated"]])
ocells <- colnames(x = integrated_seurat)

integrated_seurat@assays$RNA@counts <- integrated_seurat@assays$RNA@counts[,ocells]
integrated_seurat@assays$RNA@data <- integrated_seurat@assays$RNA@data[,ocells]
integrated_seurat@assays$integrated@data <- integrated_seurat@assays$integrated@data[,ocells]
integrated_seurat@assays$integrated@scale.data <- integrated_seurat@assays$integrated@scale.data[,ocells]
integrated_seurat@assays$SCT@counts <- integrated_seurat@assays$SCT@counts[,ocells]
integrated_seurat@assays$SCT@data <- integrated_seurat@assays$SCT@data[,ocells]
integrated_seurat@assays$SCT@scale.data <- integrated_seurat@assays$SCT@scale.data[,ocells]

# NOTE: If you subset the spatial object, re-run SCTranform() on it.
for (i in samples){
  sample.seurat <- readRDS(paste0(seurat_results, i, ".rds"))
  DefaultAssay(sample.seurat) <- "SCT"
  
  # Run ONLY if you had subset the spatial data
  # sample.seurat <- SCTransform(sample.seurat, 
  #                              assay = "Spatial",
  #                              verbose = FALSE)
  # sample.seurat <- RunPCA(sample.seurat,
  #                         verbose = FALSE)
  
  # Find anchors between reference and query
  anchors <- FindTransferAnchors(reference = integrated_seurat,
                                 query = sample.seurat,
                                 normalization.method = "SCT",
                                 reference.assay = "integrated",
                                 query.assay = "SCT")
  
  predictions.assay <- TransferData(anchorset = anchors, 
                                    refdata = integrated_seurat$cell_type,
                                    prediction.assay = TRUE,
                                    weight.reduction = sample.seurat[["pca"]],
                                    dims = 1:30)
  sample.seurat[["predictions"]] <- predictions.assay
  
  DefaultAssay(sample.seurat) <- "predictions"
  
  SpatialFeaturePlot(sample.seurat, 
                     features = rownames(predictions.assay@data), 
                     pt.size.factor = 1.6,
                     #ncol = 4,
                     crop = FALSE)
  
  ggsave(filename = paste0(seurat_results, "Labeltransfer_plot_", i, ".jpg"),
         plot = last_plot(),
         units = c("in"),
         width = 11,
         height = 8)
  
  
  
  # sample.seurat <- FindSpatiallyVariableFeatures(sample.seurat,
  #                                                assay = "predictions",
  #                                                selection.method = "moransi",
  #                                                features = rownames(sample.seurat),
  #                                                r.metric = 5,
  #                                                slot = "data")
  
  # top.clusters <- head(SpatiallyVariableFeatures(sample.seurat, 
  #                                                selection.method = "moransi"), 4)
  # SpatialPlot(object = sample.seurat, 
  #             features = top.clusters,
  #             ncol = 2)
}




