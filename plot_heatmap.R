#******************************************************************************#
#                                 HEATMAP PLOT                                 #
#******************************************************************************#

# Function to plot heatmap
# NOTE: The input normalized_counts MUST be a dataframe.  
# (i) First column of normalized_counts MUST be "SYMBOL". Duplicates ok.
# (ii) plot_genes is a vector of genes to be plotted
# (iii) disp_genes is a vector of plotted gene names that must be labelled
# (iv) metadata_column MUST contain the column names specified in anno_columns
# metadata_column MUST have rownames.
# (v) metadata_row MUST contain a column named "SYMBOL" containing gene names 
# without duplicated gene names. metadata_row MUST have rownames.
# (vi) If plotting gene signature score, name the column as "Score". The script
# has been adjusted to use same color as survival curves for HIGH LOW samples

# NOTE: Other variables that MUST be defined are:
# perform_log_transform; perform_scaling; row_clustering; col_clustering;
# row_clustering_alphabetical; col_clustering_alphabetical;
# gaps_in_row; gap_rows; gaps_in_col; gap_columns anno_columns; anno_rows; 
# my_palette; color_by_cols, color_by_rows 
# columns (this defines which column of metadata will be plotted as columns, 
# usually the column containing the sample names)
plot_heatmap <- function(normalized_counts, metadata_column, metadata_row, 
                         plot_genes, disp_genes, file_suffix, file_format, 
                         results_path, bar_width, bar_height, expr_legend){
  
  #****************************************************************************#
  # Format the matrix for heatmap
  #****************************************************************************#
  
  mat <- normalized_counts %>%
    # Keep only genes that need to be plotted
    dplyr::filter(str_to_upper(SYMBOL) %in% str_to_upper(plot_genes)) %>%
    # If there are duplicated genes, keep only data for highest expressing copy
    dplyr::mutate(n = rowSums(.[,-1])) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(n) %>%
    dplyr::ungroup() %>%
    # Duplicated genes with 0 expression in all samples still remain, remove them
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    dplyr::select(everything(), -n) %>%
    # Move gene names to rownames
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("SYMBOL") %>%
    # Make sure all values are numeric
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
  # mat[, unlist(lapply(mat, is.numeric))]    #alternative to mutate(across())
  
  # Perform log transform if needed. Count data is usually skewed right or left.
  # So, it will mostly be red or blue. log transform to make it less skewed.
  if (perform_log_transform == TRUE){
    mat <- log(1+mat, base = 2)
  }
  
  # Perform scaling for each gene across samples/cells if needed. Since scale() 
  # performs only column scaling, we transpose the dataframe first, so that 
  # genes are on columns & cells are on rows & then perform scaling.
  if (perform_scaling == TRUE){
    mat <- mat %>% t() %>% scale() %>% t()
  }
  
  # Replace NA values with 0
  mat[is.na(mat)] <- 0
  
  # Remove rows which have 0 in all samples
  mat <- mat[rowSums(mat) != 0,]
  
  # Keep only genes that need to be plotted
  mat <- mat[intersect(plot_genes, rownames(mat)), ]
  
  # Keep ONLY samples common in metadata_column and mat
  metadata_column <- metadata_column %>% 
    dplyr::filter(make.names(get(columns)) %in% make.names(colnames(mat)))
  
  # Arrange samples in mat in the same order as in metadata_column. 
  # NOTE: This is important because in the next step we assign rownames to 
  # col_annotation assuming identical sample order between metadata_column & mat
  mat <- mat[,metadata_column[,columns]]
  
  #****************************************************************************#
  # Define column and row annotations
  #****************************************************************************#
  
  if (gtools::invalid(anno_columns)){
    col_annotation <- NA
  } else {
    col_annotation <- metadata_column %>% dplyr::select(all_of(anno_columns))
    rownames(col_annotation) <- colnames(mat)
  }
  
  # Define row annotation for genes
  if (gtools::invalid(anno_rows)){
    row_annotation <- NA
  } else {
    row_annotation <- dplyr::left_join(x = mat %>% as.data.frame() %>% tibble::rownames_to_column("SYMBOL"),
                                       y = metadata_row,
                                       by=c("SYMBOL"="SYMBOL")) %>%
      dplyr::select(everything(), -colnames(mat)) %>%
      dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
      tibble::column_to_rownames("SYMBOL")
  }
  
  #****************************************************************************#
  # Define colors column annotation
  #****************************************************************************#
  
  ann_colors_col <- list()
  ann_colors_row <- list()
  colors <- c("#BF812D", "#35978F", "#C51B7D", "#7FBC41", "#762A83",
              "#E08214", "#542788", "#D6604D", "#4393C3", "#878787",
              "#E41A1C", "#F781BF", "#4DAF4A", "#FFFFBF", "#377EB8", 
              "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#999999",
              "#66C2A5", "#FC8D62", "#000000", "#9E0142", "#1A1A1A",
              "#74006F", "#FFC606", "#F6D2E0", "#C8E7F5")
  #colors <- c("#74006F", "#FFC606")    # Female: Purple, Male:Gold
  #colors <- c("#F6D2E0", "#C8E7F5")    # Female: Pink,   Male:Blue (BBN paper)
  
  if (color_by_cols == TRUE){
    for (y in 1:length(anno_columns)){
      
      # Find number of different elements in each column annotation
      # Eg: Sex has 2: Male, Female; Race has 3: White, Asian, African
      elements <- sort(unique(metadata_column[,anno_columns[y]]))
      palette_n <- length(elements)
      
      if (length(anno_columns) > 1){
        # Create a vector of transparency values. 
        alphas <- base::seq(from=0, to=1, by=1/palette_n)
        
        #Remove alpha=0 as it is always white
        alphas <- setdiff(alphas, 0)
        
        # Create a color palette with different transparencies 
        palette_colors <- colorspace::adjust_transparency(col = colors[y], alpha = alphas[1])
        for (n in 2:palette_n){
          palette_colors <- c(palette_colors, 
                              colorspace::adjust_transparency(col = colors[y], alpha = alphas[n]))
        }
      } else{
        palette_colors <- colors[1:palette_n]
      }
      palette_colors <- rev(palette_colors)
      names(palette_colors) <- rev(elements)  # reversed so order is "LOW","HIGH"
      ann_colors_col <- c(ann_colors_col, list(palette_colors))
    }
    names(ann_colors_col) <- gsub(pattern="\\.", replacement=" ", x=anno_columns)
  }
  
  ann_colors_col$Score[[1]] <- "#0c2c84"
  ann_colors_col$Score[[2]] <- "#d73027"
  
  if (color_by_rows == TRUE){
    for (y in 1:length(unique(metadata_row[,anno_rows]))){
      
      # Find number of different elements in each row annotation
      # Eg: Sex has 2: Male, Female; Race has 3: White, Asian, African
      elements <- sort(unique(unlist(mat[y,], use.names=FALSE)))
      palette_n <- length(elements)
      
      # Create a vector of transparency values. 
      alphas <- base::seq(from=0, to=1, by=1/palette_n)
      
      #Remove alpha=0 as it is always white
      alphas <- setdiff(alphas, 0)
      
      # Create a color palette with different transparencies 
      palette_colors <- colorspace::adjust_transparency(col = colors[y], alpha = alphas[1])
      for (n in 2:palette_n){
        palette_colors <- c(palette_colors, 
                            colorspace::adjust_transparency(col = colors[y], alpha = alphas[n]))
      }
      names(palette_colors) <- elements
      ann_colors_row <- c(ann_colors_row, list(palette_colors))
    }
    names(ann_colors_row) <- gsub(pattern="\\.", replacement=" ", x=rownames(mat))
  }  
  
  # $Sample
  # FB1       FB2       FB3       FB4       FB5        FC       MB1     
  # "#BF812D" "#35978F" "#C51B7D" "#7FBC41" "#762A83" "#E08214" "#542788"
  # $Sex
  # Female      Male
  # "#9E0142" "#E41A1C"
  # $Condition
  # Tumor    Normal
  # "#377EB8" "#4DAF4A"
  
  #****************************************************************************#
  # Determine breaks for heatmap color scale.
  # breaks correspond to numerical ranges for color palette's bins 
  # i.e. 0 to length(my_palette)
  #****************************************************************************#
  
  if(max(mat) == 0){
    breaks <- c(seq(from = min(mat), to = 0, length.out = 100))
    my_palette <- my_palette[1:50]
  } else if (min(mat) == 0){
    breaks <- c(seq(from = 0, to = max(mat), length.out = 100))
    my_palette <- my_palette[50:100]
  } else if(min(mat) < -3 | max(mat) > 3){
    breaks <- c(seq(-1.5, 0, length.out = 50), seq(1.5/100, 1.5, length.out = 50))
  } else{
    breaks <- c(seq(from = min(mat), to = 0, length.out = 50), seq(from = max(mat)/100, to = max(mat), length.out = 50))
  }
  
  #****************************************************************************#
  # Define vectors indicating positions where you want to have gaps in heatmap
  #****************************************************************************#
  
  if (gaps_in_row == TRUE){
    
    # Count() automatically arranges by alphabetical order unfortunately
    # So, we do this manually.
    element_names <- c()
    element_counts <- c()
    c <- 0
    for (i in 1:nrow(row_annotation)){
      if (!(row_annotation[,anno_rows][i] %in% element_names)){
        element_names <- c(element_names, row_annotation[,anno_rows][i])
        element_counts <- c(element_counts, c)
        c <- 1
      } else{
        c <- c+1
      }
    }
    element_counts <- c(element_counts,c)
    
    gaps_row <- data.frame("Description" = element_names, n = element_counts[-1]) %>%
      dplyr::mutate(n = cumsum(n)) %>%
      dplyr::select(n) %>%
      unlist(use.names = FALSE)
  } else{
    gaps_row <- NULL
  }
  
  if (gaps_in_col == TRUE){
    
    # Count() automatically arranges by alphabetical order unfortunately
    # So, we do this manually.
    element_names <- c()
    element_counts <- c()
    c <- 0
    for (i in 1:nrow(col_annotation)){
      if (!(col_annotation[,gap_columns][i] %in% element_names)){
        element_names <- c(element_names, col_annotation[,gap_columns][i])
        element_counts <- c(element_counts, c)
        c <- 1
      } else{
        c <- c+1
      }
    }
    element_counts <- c(element_counts,c)
    
    gaps_col <- data.frame("Description" = element_names, n = element_counts[-1]) %>%
      dplyr::mutate(n = cumsum(n)) %>%
      dplyr::select(n) %>%
      unlist(use.names = FALSE)
  } else {
    gaps_col <- NULL
  }
  
  #****************************************************************************#
  # Define how samples will be arranged in the heatmap
  #****************************************************************************#
  
  if (row_clustering_alphabetical == TRUE){
    mat <- mat[sort(rownames(mat)),]
  } 
  
  if (col_clustering_alphabetical == TRUE){
    mat <- mat[,sort(colnames(mat))]
  } 
  
  if(row_clustering == TRUE){
    # cluster and re-order rows
    rowclust <- hclust(dist(mat))
    reordered <- mat[rowclust$order,]
  } else{
    reordered <- mat
  }
  
  if(col_clustering == TRUE){
    # cluster and re-order columns
    colclust <- hclust(dist(t(mat)))
    reordered <- reordered[, colclust$order]
  } else{
    reordered <- reordered
  }
  
  # If you have 4 groups and want to cluster samples within each of these groups
  if(col_clustering_within_group == TRUE){
    element_names <- unique(col_annotation$Group)
    col_order <- c()
    for (g in element_names){
      temp_mat <- reordered[,rownames(col_annotation)[which(col_annotation == g)]] 
      colclust <- hclust(dist(t(temp_mat)))
      col_order <- c(col_order, colnames(temp_mat[,colclust$order]))
    }
  } else{
    col_order <- colnames(reordered)
  }
  
  if(row_clustering_within_group == TRUE){
    element_names <- unique(row_annotation$Group)
    row_order <- c()
    for (g in element_names){
      temp_mat <- reordered[rownames(row_annotation)[which(row_annotation == g)],]
      rowclust <- hclust(dist(temp_mat))
      row_order <- c(row_order,rownames(temp_mat[rowclust$order,]))
    }
  } else{
    row_order <- rownames(reordered)
  }
  reordered <- reordered[row_order,col_order]
  
  #****************************************************************************#
  # List genes and samples you want to display in the plot
  #****************************************************************************#
  
  display_col <- colnames(reordered)
  display_row <- data.frame("gene" = rownames(reordered)) %>%
    dplyr::mutate(gene = dplyr::case_when(gene %in% disp_genes ~ gene, 
                                          TRUE ~ "")) %>% 
    unlist(use.names=FALSE)
  
  if (ncol(reordered) > 500){
    
    # Save the clustered scores in xlsx
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "Heatmap_matrix")
    openxlsx::writeData(wb, sheet = "Heatmap_matrix", x = t(reordered), rowNames = TRUE)
    openxlsx::saveWorkbook(wb, file = paste0(results_path, "Heatmap_matrix", file_suffix, ".xlsx"), 
                           overwrite = TRUE)
  } else{
    # Save the clustered scores in xlsx
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "Heatmap_matrix")
    openxlsx::writeData(wb, sheet = "Heatmap_matrix", x = reordered, rowNames = TRUE)
    openxlsx::saveWorkbook(wb, file = paste0(results_path, "Heatmap_matrix", file_suffix, ".xlsx"), 
                           overwrite = TRUE)
  }
  
  #****************************************************************************#
  # Plot heatmap
  #****************************************************************************#
  pheatmap::pheatmap(mat = as.matrix(reordered),
                     color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
                     breaks = breaks, 
                     border_color = NA, #"grey90", "white"
                     cellwidth = bar_width,
                     cellheight = bar_height,
                     # cellwidth = dplyr::case_when(ncol(reordered) > 25 ~ NA,
                     #                              ncol(reordered) <= 25 ~ 6), 
                     # cellheight = dplyr::case_when(nrow(reordered) > 25 ~ NA,
                     #                               nrow(reordered) <= 25 ~ 6), 
                     scale = "none",   
                     cluster_rows = row_clustering,   #cluster the rows
                     cluster_cols = col_clustering,   #cluster the columns
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method = "complete",
                     legend = expr_legend,  # set FALSE if it overlaps with annotation legends
                     legend_breaks = NA,
                     legend_labels = NA, 
                     annotation_row = row_annotation,
                     annotation_col = col_annotation,
                     annotation_colors = ann_colors_col,
                     annotation_legend = TRUE,
                     annotation_names_row = FALSE,
                     annotation_names_col = FALSE, #set to TRUE if more than 1 
                     show_rownames = dplyr::if_else(length(disp_genes) < 80, TRUE, FALSE, missing = NULL), 
                     show_colnames = dplyr::if_else(length(unique(display_col)) < 50, TRUE, FALSE, missing = NULL),
                     fontsize = 5, 
                     fontsize_row = 5, 
                     fontsize_col = 5,
                     gaps_row = gaps_row,
                     gaps_col = gaps_col,
                     angle_col = "45",
                     fontsize_number = 0.8*fontsize, 
                     labels_row = display_row, 
                     labels_col = display_col,
                     #width = 11,
                     #height = 11,
                     filename = paste0(results_path, "Heatmap_", file_suffix, ".", file_format))
}
