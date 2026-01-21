#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

# Data wrangling packages   
library("openxlsx")
library("dplyr")
library("tibble")
library("stringr")
library("purrr")
library("Matrix")

# Graph plotting packages
library("ggplot2")
library("cowplot")
library("viridis")
library("RColorBrewer")
library("ggrepel")
library("ggpubr")
library("extrafont")
library("ggbeeswarm")

# Specialized Graph plotting packages
library("ComplexHeatmap")
library("pheatmap")

#******************************************************************************#
#                           DEFINE GLOBAL PARAMETERS                           #
#******************************************************************************#

parent_path <- "C:/Users/KailasammS/Desktop/"
results_path <- "C:/Users/KailasammS/Desktop/"

# Define any filename you want added to final file
file_suffix <- NULL

# Define cutoff for padj and log2FC
padj_cutoff <- 0.05
log2FC_cutoff <- 0.58

# Define aspect ratio
# NOTE: width of plot has been defined in plot_bar() as height*aspect_ratio.
# So, large the aspect ratio, wider the plot
aspect_ratio <- 1

# Assign colors for UMAP. The current palette supports up to 33 cell types
my_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", 
                "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                "#FFC61E", "#762A83", "#333333", "#FF1F5B", "#B8E80C",
                "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
                "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5",
                "#FFFFBF", "#C51B7D", "#DE77AE", "#F1B6DA", "#FDE0EF", 
                "#F7F7F7", "#E6F5D0", "#B8E186")

# Define axis font etc to use in all plots
my_theme <- ggplot2::theme(aspect.ratio = aspect_ratio,
                           plot.title =   element_text(family = "sans", face = "bold",  colour = "black", size = 15, hjust = 0.5),
                           axis.title.x = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0.5, vjust = 0,   angle = 0),
                           axis.title.y = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0.5, vjust = 1,   angle = 90),
                           legend.title = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0,   vjust = 1,   angle = 0),
                           axis.text.x =  element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5, vjust = 0.5, angle = 45),
                           axis.text.y =  element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5, vjust = 0.5, angle = 0),
                           legend.text =  element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5),
                           strip.text.x = element_text(family = "sans", face = "bold",  colour = "black", size = 10, hjust = 0.5),
                           #legend.background = element_rect(fill = "lightblue", size = 0.5, linetype = "solid", colour = "darkblue"),
                           legend.position = "right",
                           legend.justification = "left",
                           legend.direction = "vertical",
                           legend.key.height = unit(0.5, 'cm'),
                           legend.key.width  = unit(0.5, 'cm'), 
                           legend.text.align = 0)

#******************************************************************************#
#                         DEFINE PARAMETERS FOR GGPLOT2                        #             
#******************************************************************************#

plot_param <- list(
  
  # column to be plotted on x axis
  "data_x" = c("Sample"),                 #"k.K", "Count, "Sample"
  
  # Define title of x axis
  "title_x" = c(),  #"Encrichment Ratio", "Number of Genes"
  
  # column to be plotted on y axis
  "data_y" = c("nUMIs"), #"Description", "nUMIs"
  
  # Define title of y axis
  "title_y" = c("nUMIs"),   #"Percent Composition", "nUMIs"
  
  # column to be used for filling bar colors
  "data_fill" = c("Sample"), #"FDR", "Subgroups", "Sample"
  
  # column to be used for coloring the dots
  "data_color" = c(), #"FDR", "Subgroups"
  
  # column to be used for determining the size of dots
  "data_size" = c(), #"Count"
  
  # Define title of legend
  "title_legend_fill" = c("Sample"),   #"log10(padj)", "Subgroups", "Sample"
  "title_legend_color" = c(),             #"log10(padj)", "Subgroups"
  "title_legend_size" = c(),              #"Number of Genes"
  
  # Define title of plot
  "title_plot" = c("Distribution of UMIs") #"Distribution of UMIs"
)

#******************************************************************************#
#                      IMPORT DATA FOR BAR CHART/DOT PLOT                      #             
#******************************************************************************#

# NOTE: Bar plots gives information only on enrichment ratio & pvalue while dot
# plots give info on number of genes in each pathway additionally.

# NOTE: data MUST have ALL of the columns you define above
# (i) data_y ~ "Description" column containing pathway names
# (ii) data_x ~ "k.K" column containing ratio or total number of genes
# (iii) data_fill/data_color ~ "FDR" column containing log10(padj) values
# (iv) data_size ~ "Count" column containing total number of genes
data <- openxlsx::read.xlsx("C:/Users/KailasammS/Desktop/Data_metascape_DDR2KO.xlsx")

# data <- data %>%
#   dplyr::filter(!!rlang::sym(plot_param$data_fdr) < log10(padj_cutoff)) %>%
#   dplyr::mutate(!!rlang::sym(plot_param$data_pathway) := stringr::str_wrap(get(plot_param$data_pathway), width = 22)) %>%
#   dplyr::slice_min(order_by = get(plot_param$data_fdr), n = 10)

data <- data %>%
  dplyr::filter(FDR < log10(padj_cutoff)) %>%
  dplyr::mutate(Description = stringr::str_wrap(Description, width = 22)) %>%
  dplyr::slice_min(order_by = FDR, n = 10)

# Plot bar chart
plot_bar(data, plot_param)

# Plot dot chart
plot_dot(data, plot_param)

#******************************************************************************#
#                      IMPORT DATA FOR STACKED BAR CHART                       #             
#******************************************************************************#

# Define if you want to label percentages within stacked bar chart
label_percent <- "TRUE"

# Define if input data is already converted to % and can be be plotted as is
already_percent <- "FALSE"

# NOTE: data MUST have ALL of the columns you define above
# (i) data_x ~ "Sample" column containing sample names. This will be on X-axis
# (ii) Rest of columns MUST be the subtype names i.e. sub-categories within each
# sample which will be stacked in each bar which will automatically be named 
# "Subgroups

# Sample  Epithelial-Basal  Epithelial-Luminal
# FB1     350               900
# FB2     700               100

data <- openxlsx::read.xlsx("C:/Users/KailasammS/Desktop/y.xlsx")

# Plot stacked bar chart
plot_stackedbar(data, plot_param, label_percent)

# celltype <- "Epithelial"
 
# # Load the integrated seurat object
# integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
#                                     dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))
 
# # Remove unwanted cells and samples
# integrated_seurat <- subset(x = integrated_seurat,
#                             subset = (cell_class %in% c("Mixed", "Unclassified") | (Condition %in% c("Normal"))),
#                             invert = TRUE)
 
# # Calculate cells per cluster
# data <- integrated_seurat@meta.data %>% 
#   dplyr::group_by(Sample, sub_type) %>%
#   dplyr::count() %>%
#   tidyr::pivot_wider(id_cols = Sample, names_from = sub_type, values_from = n) %>%
#   dplyr::rename_with(make.names)

#******************************************************************************#
#                       IMPORT DATA FOR VIOLIN/BOX PLOT                        #             
#******************************************************************************#

# Define if y axis should be in log scale
log_scale_y <- TRUE

# Define if y axis should have an intercept
show_intercept_y <- TRUE

# Define the value for y axis intercept
yintercept_cutoff <- 5

# NOTE: data MUST have ALL of the columns you define above
# (i) data_x ~ "Sample" column containing sample names. This will be on X-axis
# (ii) data_y ~ "nUMIs" column containing y axis values.

# Sample  nUMIs   
# VC      4.2    
# VC      11.5    
# VC      7.3     
# CJ      5.8     
# CJ      6.4     

data <- openxlsx::read.xlsx("C:/Users/KailasammS/Desktop/x.xlsx")

# Plot violin plot
plot_violin(data, plot_param)

# Plot box plot
plot_box(data, plot_param)

#******************************************************************************#
#                        IMPORT DATA FOR HISTOGRAM PLOT                        #             
#******************************************************************************#

# Define bin_width for histogram
bin_width <- 1

# NOTE: data MUST have ALL of the columns you define above
# (i) data_x ~ "Sample" column containing sample names. This will be on X-axis
# (ii) data_y ~ "nUMIs" column containing y axis values.
data <- openxlsx::read.xlsx("C:/Users/KailasammS/Desktop/x.xlsx")

# Plot histogram
plot_histogram(data, plot_param)

#******************************************************************************#
#               IMPORT DATA & DEFINE PARAMETERS FOR VOLCANO PLOT               #             
#******************************************************************************#

# Import expression data without log2FC and pval
# NOTE: data_mat MUST have 
# (i) "SYMBOL" as rownames, not in 1st column
# (ii) rest of column names MUST match with metadata$Sample
# Values should be raw i.e. untransformed. DO NOT use log transformed values etc.
# DO NOT replace NA with 0 etc. Leave NA as they are.
# NA values will be imputed/replaced with average of non-NA values.
# If all values are NA, they will be set to 0.
# It also performs quantile normalization if specified for mass spec data

data_mat <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Data_Mass_Spec.xlsx")) %>%
  tibble::column_to_rownames(colnames(.)[1])
data_mat <- log(1+data_mat, base = 2)

# Define/import metadata
# NOTE: metadata MUST have 
# (i) "Sample" as 1st column containing sample names
# (ii) "Condition" as 2nd column which MUST match with "Target" and "Reference"
# variables defined above
# (iii) "Source" as 3rd column which MUST contain the celltype

metadata <- data.frame("Sample" = colnames(data_mat), 
                       "Condition" = c(rep(x = "Control", times = 3), rep(x = "KO", times = 3)),
                       "Source" = c(rep(x = "", times = 6)))

# Define the target and reference groups
# NOTE: These values MUST be present in "Condition" column of metadata
Target <- "KO"
Reference <- "Control"

# Calculate pval and log2FC if not already calculated 
volcano_df <- calc_stats(data_mat, metadata)

# (II) Alternatively, import expression data with log2FC and pval
# NOTE: volcano_df MUST have "Gene", "padj", "log2FC", "Source" columns
# NOTE: Output of calc_stats() can be used as input for plot_volcano()
volcano_df <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Results_metabolomics_Control_vs_KO.xlsx"))

# (III) Define any genes you want to mark in volcano plot
disp_genes <- c("Guanidinoacetate", "Phosphoserine")
disp_genes <- volcano_df %>%
  dplyr::filter(padj < 0.05 & log2FC > 0) %>%
  dplyr::select(Gene) %>%
  unlist(use.names = FALSE) %>%
  unique()

padj_cutoff <- 0.05
log2_cutoff <- 0
color_by <- "Direction"
same_color <- "FALSE"
draw_line <- "TRUE"
plot_param <- list(
  
  # column to be plotted on x axis
  "data_x" = c(),                 #"k.K", "Count, "Sample"
  
  # Define title of x axis
  "title_x" = c(),  #"Encrichment Ratio", "Number of Genes"
  
  # column to be plotted on y axis
  "data_y" = c(), #"Description", "nUMIs"
  
  # Define title of y axis
  "title_y" = c(),   #"Percent Composition", "nUMIs"
  
  # column to be used for filling bar colors
  "data_fill" = c(), #"FDR", "Subgroups", "Sample"
  
  # column to be used for coloring the dots
  "data_color" = c(), #"FDR", "Subgroups"
  
  # column to be used for determining the size of dots
  "data_size" = c(), #"Count"
  
  # Define title of legend
  "title_legend_fill" = c(),   #"log10(padj)", "Subgroups", "Sample"
  "title_legend_color" = c("log10(padj)"),             #"log10(padj)", "Subgroups"
  "title_legend_size" = c(1),              #"Number of Genes"
  
  # Define title of plot
  "title_plot" = c() #"Distribution of UMIs"
)

# Make volcano plots
plot_volcano(volcano_df, disp_genes,  file_suffix)

#******************************************************************************#
#                 IMPORT DATA & DEFINE PARAMETERS FOR HEATMAP                  #             
#******************************************************************************#

# Define if data is log transformed already
already_log <- FALSE

# Define if data is scaled already
already_scaled <- FALSE

# Define if you want to perform unsupervised row and column clustering
# NOTE: If set to FALSE, samples (columns) and genes(rows) will be arranged 
# in alphabetical order (default) in heatmap. If you want to arranged in 
# specific order, define below.
row_clustering <- TRUE    # Usually TRUE
col_clustering <- TRUE    # Usually FALSE

# Define if you want genes or samples to be arranged in alphabetical order
# NOTE: If set to FALSE, write the plot_genes in order you want in heatmap
# NOTE: If row_clustering == TRUE, then row_clustering_alphabetical is irrelevant
row_clustering_alphabetical <- FALSE
col_clustering_alphabetical <- FALSE

# Define if heatmap columns must have gaps and based on which column of metadata
gaps_in_col <- FALSE
gap_columns <- "Sample"   # Irrelevant if gaps_in_col is FALSE

# Define if heatmap rows must have gaps and based on which column of metadata_row
gaps_in_row <- FALSE
gap_rows <- "Pathway"    # Irrelevant if gaps_in_row is FALSE

# List annotations you want on heatmap
# NOTE: anno_columns MUST match one of the column names in metadata_column while
# anno_rows MUST match one of the column names in metadata_row
anno_columns <- c("Condition")
anno_rows <- c("Pathway")

# Define if you want to make color annotations of rows or columns or both
color_by_cols <- TRUE
color_by_rows <- FALSE

# Define colors for heatmap
my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)
# my_palette <- viridis(100)
# my_palette <- c(colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)[1:49], "#FFFFFF",
#                 colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100)[50:99])

# NOTE: normalized_counts is a dataframe with gene names in 1st column "SYMBOL"
# NOTE: normalized_counts will have genes in rows and samples in columns
# NOTE: metadata is a dataframe with sample names in 1st column "Sample"
# NOTE: metadata_row is a dataframe with gene names in 1st column "SYMBOL"

# NOTE: If there are duplicated gene names, plot_heatmap() will keep only 
# expression data for 1 copy of duplicated gene that has highest expression.
# NOTE: "Error in check.length("fill"):'gpar' element 'fill' must not be length 0"
# This error means sample_annotation dataframe doesnt match with columns of mat
# Try using mat = t(mat) to see if it fixes the error.
# NOTE: If you set scale = none, then you are plotting exact progeny scores
# NOTE: If you set scale = row, then colors do not represent progeny scores 
# (i.e. actual activity) but relative activity

# To plot heat map we need an input matrix of the following format:
#         sample1   sample2
# Gene1    x          z  
# Gene2    y          w
# The rownames of the input matrix will be gene names and colnames will be
# the sample names

# Also, you need to create an annotation matrix to pass to "annotation_col"
# parameter of pheatmap() in the following format:
#            Sex
# Sample1    Male
# Sample2    Female
# The rownames of the annotation matrix will be the sample names

# (I) Read expr data
normalized_counts <- openxlsx::read.xlsx(xlsxFile = 
                                           "C:/Users/KailasammS/Box/Saravana@cedars/10. Ongoing Projects/Boopati project/Yneg_Mass_Spec_Results.xlsx")

normalized_counts <- normalized_counts[, -1]
colnames(normalized_counts)[1] <- "SYMBOL"

# (II) Define metadata for columns if you want to annotate columns
metadata_column <- openxlsx::read.xlsx(xlsxFile = 
                                  "C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/Hany_Y/Metadata.xlsx")

metadata_column <- data.frame("Sample" = colnames(normalized_counts[,-1]), 
                       "Condition" = c("Yneg_IP", "Yneg_IP", "Yneg_iso", "Yneg_iso", "Ypos_IP", "Ypos_IP", "Ypos_iso", "Ypos_iso"))

# (III) Define metadata for rows if you want to annotate rows
metadata_row <- NULL

# (IV) Define genes to plot
plot_genes <- c("Vegfd", "Camk2b", "Mknk2", "Camk1d", "Ncoa1", "Prkd3", 
                "Pik3r5", "Crebbp", "Braf", "Ep300", "Egln3", "Nos2", "Hspa1b",
                "Vegfa", "Adm", "Ldhb", "Rras2", "Mmp2", "Mmp3")
plot_genes <- normalized_counts$SYMBOL

# (V) Genes to display in heatmap
disp_genes <- c("Aldoa", "Pgk1", "Pgam1")
disp_genes <- plot_genes

# (VI) Define which column of metadata to plot as columns in heatmap
columns <- "Sample"

file_suffix <- ""
results_path <- "C:/Users/KailasammS/Desktop/"
# Run function
plot_heatmap(normalized_counts, metadata_column, metadata_row, plot_genes, disp_genes, file_suffix, results_path)

################################################################################
