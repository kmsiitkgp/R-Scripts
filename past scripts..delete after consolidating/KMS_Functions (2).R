# This R file contains all the user defined functions which can be imported
# to analyze bulk RNA Seq, single cell RNA Seq, plot graphs etc

#******************************************************************************#
#                                   BAR PLOT                                   #
#******************************************************************************#

# Function to plot bar chart
plot_bar <- function(data, plot_param){
  
  ggplot2::ggplot(data = data,
                  aes(x = !!rlang::sym(plot_param$data_x), 
                      y = reorder(stringr::str_to_upper(string = !!rlang::sym(plot_param$data_y)), desc(!!rlang::sym(plot_param$data_fill))),
                      fill = !!rlang::sym(plot_param$data_fill))) +
    
    # Plot bar plot
    ggplot2::geom_col(width = 0.75) +
    
    # Define the theme of plot
    ggplot2::theme_classic() +
    
    # Define the axis, plot headings
    ggplot2::labs(x = plot_param$title_x,
                  y = plot_param$title_y,
                  fill = plot_param$title_legend_fill,
                  color = plot_param$title_legend_color,
                  size = plot_param$title_legend_size,
                  title = plot_param$title_plot) +
    
    # Define x-axis start and end
    ggplot2::coord_cartesian(xlim = c(0, ceiling(max(abs(data[,plot_param$data_x]))))) +
    
    # Define the color of the bars
    viridis::scale_fill_viridis(discrete = FALSE, option = "D") +
    
    # Adjust font size, style
    my_theme +
    ggplot2::theme(axis.text.y = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 1))
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Bar_plot_", file_suffix, ".pdf"),
                  plot = last_plot(),
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = 2+nrow(data),
                  height = (2+nrow(data))*aspect_ratio,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
  # Create a workbook to store the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "Bar")
  openxlsx::writeData(wb = wb, sheet = "Bar", x = data)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "Bar_plot_", file_suffix, "_Results.xlsx"),
                         overwrite = TRUE)
  
}

#******************************************************************************#
#                                   DOT PLOT                                   #
#******************************************************************************#

# Function to plot dot plot
plot_dot <- function(data, plot_param){
  
  ggplot2::ggplot(data = data,
                  aes(x = !!rlang::sym(plot_param$data_x), 
                      y = reorder(stringr::str_to_upper(string = !!rlang::sym(plot_param$data_y)), desc(!!rlang::sym(plot_param$data_fill))),
                      size = !!rlang::sym(plot_param$data_size),
                      color = !!rlang::sym(plot_param$data_color))) +
    
    # Plot dot plot
    ggplot2::geom_point() +
    
    # Define the theme of plot
    ggplot2::theme_classic() +
    
    # Define the axis, plot headings
    ggplot2::labs(x = plot_param$title_x,
                  y = plot_param$title_y,
                  color = plot_param$title_legend_color,
                  size = plot_param$title_legend_size,
                  title = plot_param$title_plot) +
    
    # Define x-axis start and end
    ggplot2::coord_cartesian(xlim = c(0, ceiling(max(abs(data[,plot_param$data_x]))))) +
    
    # Define the color of the dots
    viridis::scale_color_viridis(discrete = FALSE, option = "D") +
    
    # Adjust font size, style
    my_theme +
    ggplot2::theme(axis.text.y = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 1))
  
  #ggplot2::scale_size(breaks = sapply(as.vector(quantile(c(min(data$Count), max(data$Count)))), floor))
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Dot_plot_", file_suffix, ".pdf"),
                  plot = last_plot(),
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = 2+nrow(data),
                  height = (2+nrow(data))*aspect_ratio,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
  # Create a workbook to store the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "Dot")
  openxlsx::writeData(wb = wb, sheet = "Dot", x = data)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "Dot_plot_", file_suffix, "_Results.xlsx"),
                         overwrite = TRUE)
}

#******************************************************************************#
#                               STACKED BAR PLOT                               #
#******************************************************************************#

# Function to plot stacked bar chart
plot_stackedbar <- function(data, plot_param, label_percent){
  
  if (already_percent == FALSE){
    # Calculate percent of each sub type for each sample
    data <- data %>%
      data.frame() %>%
      base::replace(is.na(.), 0) %>%
      dplyr::mutate(across(.cols = c(everything(), -plot_param$data_x), .fns = as.numeric)) %>%
      dplyr::mutate(across(.cols = c(everything(), -plot_param$data_x), .fns = function(x) x*100/rowSums(x = select_if(., is.numeric)))) %>%
      tidyr::pivot_longer(cols = !rlang::sym(plot_param$data_x), names_to = plot_param$data_fill, values_to = "Percent") %>%
      dplyr::mutate(n_gene = gsub(pattern="X", replacement="", x=n_gene))
  }
  
  # Plot stacked bar chart
  p <- ggplot2::ggplot(data = data, 
                       aes(x = !!rlang::sym(plot_param$data_x), 
                           y = Percent, 
                           fill = !!rlang::sym(plot_param$data_fill))) +
    
    # Plot stacked bar plot
    ggplot2::geom_col(position = "stack", width = 0.95) +
    
    # Define the theme of plot
    ggplot2::theme_classic() + 
    
    # Define the axis, plot headings
    ggplot2::labs(x = plot_param$title_x,
                  y = plot_param$title_y,
                  fill = plot_param$title_legend_fill,
                  title = plot_param$title_plot) +
    
    # Define the color of the bars
    ggplot2::scale_fill_manual(values = my_palette) +
    
    # Adjust font size, style
    my_theme
  
  if(label_percent == "TRUE"){
    p <- p +
      ggplot2::geom_text(aes(x = !!rlang::sym(plot_param$data_x), label = round(Percent,1)), 
                         position = position_stack(vjust = 0.5), 
                         fontface = "plain", colour = "white", size = 3, 
                         check_overlap = TRUE)
  }
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Stacked_Bar_", file_suffix, ".pdf"),
                  #plot = last_plot(),
                  plot = p,
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = dplyr::if_else((2+nrow(data)) < 10, (2+nrow(data)), 10),
                  height = dplyr::if_else((2+nrow(data))*aspect_ratio < 10, (2+nrow(data))*aspect_ratio, 10),
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
  # Create a workbook to store the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "Stacked_Bar")
  openxlsx::writeData(wb = wb, sheet = "Stacked_Bar", x = data)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "Stacked_Bar_", file_suffix, "_Results.xlsx"),
                         overwrite = TRUE)
  
}


#******************************************************************************#
#                                 VIOLIN PLOT                                  #
#******************************************************************************#

# Function to plot violin plot
plot_violin <- function(data, plot_param){
  
  p <- ggplot(data = data, 
              aes(x = !!rlang::sym(plot_param$data_x), 
                  y = !!rlang::sym(plot_param$data_y), 
                  fill = !!rlang::sym(plot_param$data_fill))) +
    
    # Plot violin plot with a small box plot within it
    geom_violin(trim=FALSE) + 
    geom_boxplot(width = 0.1) +
    
    # Define the theme of plot
    theme_classic() + 
    
    # Define the axis, plot headings
    labs(x = plot_param$title_x, 
         y = plot_param$title_y, 
         title = plot_param$title_plot)  +
    
    # Define y-axis start and end
    #coord_cartesian(ylim = c(0, ceiling(max(data[,plot_param$data_y])/10)*10)) +
    
    # Define the color of the violins
    ggplot2::scale_fill_manual(values = my_palette) +
    
    # Adjust font size, style
    my_theme 
  
  if (log_scale_y == "TRUE"){
    p <- p +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) 	#display y axis in log scale 
  }        
  
  if (show_intercept_y == "TRUE"){
    p <- p +
      geom_hline(yintercept = yintercept_cutoff, linetype = 2)
  }
  
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Violin_plot_", file_suffix, ".pdf"),
                  #plot = last_plot(),
                  plot = p,
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = 2+length(unique(data[,plot_param$data_x])),
                  height = (2+length(unique(data[,plot_param$data_x])))*aspect_ratio,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
}

#******************************************************************************#
#                                   BOX PLOT                                   #
#******************************************************************************#

# Function to plot box plot
plot_box <- function(data, plot_param){
  
  p <- ggplot(data = data, 
              aes(x = !!rlang::sym(plot_param$data_x), 
                  y = !!rlang::sym(plot_param$data_y), 
                  fill = !!rlang::sym(plot_param$data_fill))) +
    
    # Plot box plot
    geom_boxplot() +
    
    # Define the theme of plot
    theme_classic() + 
    
    # Define the axis, plot headings
    labs(x = plot_param$title_x, 
         y = plot_param$title_y, 
         title = plot_param$title_plot)  +
    
    # Define y-axis start and end
    #coord_cartesian(ylim = c(0, ceiling(max(data[,plot_param$data_y])/10)*10)) +
    
    # Define the color of the violins
    ggplot2::scale_fill_manual(values = my_palette) +
    
    # Adjust font size, style
    my_theme 
  
  if (log_scale_y == "TRUE"){
    p <- p +
      scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) 	#display y axis in log scale 
  }        
  
  if (show_intercept_y == "TRUE"){
    p <- p +
      geom_hline(yintercept = yintercept_cutoff, linetype = 2)
  }
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Box_plot_", file_suffix, ".pdf"),
                  #plot = last_plot(),
                  plot = p,
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = 2+length(unique(data[,plot_param$data_x])),
                  height = (2+length(unique(data[,plot_param$data_x])))*aspect_ratio,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
  
}

#******************************************************************************#
#                                HISTOGRAM PLOT                                #
#******************************************************************************#

# Function to plot histogram
plot_histogram <- function(data, plot_param){
  
  ggplot(data = data, 
         aes(x = !!rlang::sym(plot_param$data_x), 
             color = !!rlang::sym(plot_param$data_color),
             fill = !!rlang::sym(plot_param$data_fill))) +
    
    # Plot box plot
    geom_histogram(binwidth = bin_width) +
    
    # Define the theme of plot
    theme_classic() + 
    
    # Define the axis, plot headings
    labs(x = plot_param$title_x, 
         y = plot_param$title_y, 
         title = plot_param$title_plot)  +
    
    # Define y-axis start and end
    #coord_cartesian(ylim = c(0, ceiling(max(data[,plot_param$data_y])/10)*10)) +
    
    # Adjust font size, style
    my_theme 
  
  # Save the plot
  ggplot2::ggsave(filename = paste0("Histogram_plot_", file_suffix, ".pdf"),
                  #plot = last_plot(),
                  plot = p,
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  width = 2+length(unique(data[,plot_param$data_x])),
                  height = (2+length(unique(data[,plot_param$data_x])))*aspect_ratio,
                  units = c("in"),
                  dpi = 600,
                  limitsize = TRUE,
                  bg = NULL)
}

#******************************************************************************#
#                                 VOLCANO PLOT                                 #
#******************************************************************************#

# Function to calculate pval and log2FC
calc_stats <- function(data_mat, metadata, file_suffix){
  
  genes <- c()
  expt <- c()
  control <- c()
  pval <- c()
  
  for (j in 1:nrow(data_mat)){
    
    data <- data_mat[j,] %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Sample") %>%
      dplyr::left_join(metadata, by = c("Sample" = "Sample"))
    colnames(data) <- c("Sample", "values", "Condition")
    
    # If all values for Reference == NA or Target == NA, set them to 0. 
    # If at least 1 value is not NA, do not make any changes. 
    if (sum(data[data$Condition == Reference, ]$values, na.rm=TRUE) == 0){
      data[data$Condition == Reference, ]$values <- rep(x=0, times=sum(metadata$Condition == Reference))
    }
    if (sum(data[data$Condition == Target, ]$values, na.rm=TRUE) == 0){
      data[data$Condition == Target, ]$values <- rep(x=0, times=sum(metadata$Condition == Target))
    }
    
    # Remove rows which have NA in values.
    # NOTE: This step not needed as f.test/t.test automatically removes entries
    # which have NA
    #data <- data[!is.na(data$values),]
    
    # Use F.test() to determine if variances are equal or not.
    # NOTE: p < 0.05 => unequal variance
    # NOTE: If 3 values for Iso are "NA","NA","18" and 3 values for IP are 
    # "17", "NA", "18", var.test and t.test with throw error since Iso has only
    # 1 value.
    # NOTE: If 3 values for Iso are 1,1,1 and 3 values for IP are 1,1,1, t.test 
    # will throw error "data are essentially constant".
    # NOTE: If 3 values for Iso are 0,0,0 and 3 values for IP are 1,1,1, t.test 
    # will throw error "data are essentially constant".
    
    if (sum(!is.na(data[data$Condition == Reference, ]$values)) > 1 &
        sum(!is.na(data[data$Condition == Target, ]$values)) > 1 & 
        length(unique(data$values)) > 1 &
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
      
      # Calculate p values, mean expression
      t_test <- t.test(formula = values ~ Condition, 
                       data = data,
                       alternative = "two.sided",
                       var.equal = TRUE)
      
      if (grepl(Reference, names(t_test$estimate[1]))){
        genes <- c(genes, rownames(data_mat)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[1]])
        expt <- c(expt, t_test$estimate[[2]])
      } else if (grepl(Reference, names(t_test$estimate[2]))) {
        genes <- c(genes, rownames(data_mat)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[2]])
        expt <- c(expt, t_test$estimate[[1]])
      }
    } else{
      genes <- c(genes, rownames(data_mat)[j])
      pval <- c(pval, 1) # Note: DO NOT SET to NA. It will increase padj.
      control <- c(control, mean(data[data$Condition == Reference, ]$values, na.rm=TRUE))
      expt <- c(expt, mean(data[data$Condition == Target, ]$values, na.rm=TRUE))
    }
  }
  
  stats_df <- data.frame(genes, expt, control, pval)
  stats_df$padj <- stats::p.adjust(p = stats_df$pval, method = "fdr", n = length(stats_df$pval))
  stats_df$log2FC <- stats_df$expt - stats_df$control
  
  result <- data_mat %>% 
    tibble::rownames_to_column(var = "Gene") %>% 
    dplyr::left_join(stats_df, by = c("Gene" = "genes"))
  
  # Save the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "1")
  openxlsx::writeData(wb, sheet = "1", x = result, rowNames = FALSE)
  openxlsx::saveWorkbook(wb, 
                         file = paste0(results_path, "Volcano_Results_", file_suffix, ".xlsx"), 
                         overwrite = TRUE)
  
  return(result)
}

# Function to plot volcano plots
plot_volcano <- function(volcano_df, disp_genes, file_suffix){
  
  # Categorize the data points
  volcano_df <- volcano_df %>%
    dplyr::mutate(Direction = dplyr::case_when(padj < padj_cutoff & log2FC > log2_cutoff ~ paste0("Up in ", Source, "_", Target),
                                               padj < padj_cutoff & log2FC < -log2_cutoff ~ paste0("Up in ", Source, "_", Reference),
                                               TRUE ~ "Not Significant"),
                  padj = dplyr::case_when(is.na(padj) ~ 0,
                                          TRUE ~ padj),
                  Significance = dplyr::case_when(abs(log2FC) >= log2_cutoff & padj <= 0.05 & padj > 0.01 ~ "FDR < 0.05",
                                                  abs(log2FC) >= log2_cutoff & padj <= 0.01 & padj > 0.001 ~ "FDR < 0.01",
                                                  abs(log2FC) >= log2_cutoff & padj <= 0.001  ~ "FDR < 0.001",
                                                  TRUE ~ "Not Significant"))
  
  # Define the limits of the x-axis and y-axis
  x_limits <- as.vector(quantile(volcano_df$log2FC, probs = c(0,1)))
  if (is.infinite(max(x_limits))){
    x_limits <- as.vector(quantile(volcano_df$log2FC, probs = c(0,0.999)))
  }
  y_limits <- as.vector(quantile(-log10(volcano_df$padj), na.rm = TRUE, probs = c(0,1)))
  if (is.infinite(max(y_limits))){
    y_limits <- as.vector(quantile(-log10(volcano_df$padj), na.rm = TRUE, probs = c(0,0.999)))
  }
  
  # NOTE: If using labels, sort labels in alphabetical order and then assign 
  # color because R by default will arrange the labels in alphabetical order 
  # first and then match them to colors indicated in values vector and then 
  # color the plot. The coloring in the legend is however dependent on the 
  # order of labels vector and values vector. To understand, create plot first 
  # using the sort and then without the sort(). 
  
  # NOTE: DO NOT USE labels for defining colors due to reasons above. 
  # RECOMMEND using a named vector 
  #volcano_palette <- c("#808080", RColorBrewer::brewer.pal(11, "RdYlBu")[c(11)], RColorBrewer::brewer.pal(11, "RdYlBu")[c(1)])
  if (color_by == "Significance"){
    volcano_palette <- c(viridis(4)[4], viridis(4)[3], viridis(4)[2], viridis(4)[1])
    names(volcano_palette) <- c("Not Significant", "FDR < 0.05", "FDR < 0.01", "FDR < 0.001")
  } else if (color_by == "Direction" & same_color == "TRUE"){
    x <- unique(volcano_df$Direction)
    volcano_palette <- dplyr::case_when(grepl(Target, x) ~ "orange", 
                                        grepl(Reference, x) ~ "purple", 
                                        TRUE ~ "grey")
    names(volcano_palette) <- unique(volcano_df$Direction)
  } else if (color_by == "Direction" & same_color == "FALSE"){
    x <- unique(volcano_df$Direction)
    volcano_palette <- c("grey", my_palette[1:(length(x)-1)])
    names(volcano_palette) <- unique(volcano_df$Direction)
  }
  
  ggplot2::ggplot(data = volcano_df, 
                  aes(x = log2FC, 
                      y = -log10(padj), 
                      label = Gene,
                      color = get(color_by), 
                      shape = Direction,
                      size = 2)) +
    
    # Plot dot plot
    ggplot2::geom_point() +
    
    # Define the theme of plot
    ggplot2::theme_classic() +
    
    # Define the axis, plot headings
    ggplot2::labs(x = expression("log"[2]*"FC"),
                  y = expression("-log"[10]*"padj"),
                  color = plot_param$title_legend_color,
                  size = plot_param$title_legend_size,
                  shape = "Direction",
                  color = color_by,
                  title = paste0(Target, " vs ", Reference)) +   
    
    # Draw line to mark the cutoffs
    geom_vline(xintercept = c(-log2_cutoff,log2_cutoff), color = "black", linetype = "dotted", linewidth = 0.5) +
    geom_hline(yintercept = -log10(padj_cutoff), color = "black", linetype = "dotted", linewidth = 0.5) +
    
    # Define the axis tick marks
    scale_x_continuous(breaks = seq(-ceiling(max(abs(x_limits))), ceiling(max(abs(x_limits))), 1)) +
    #scale_y_continuous(breaks = seq(0, ceiling(y_limits[2]/10)*10, 20)) +
    
    # Define x-axis start and end
    coord_cartesian(xlim = c(-ceiling(max(abs(x_limits))), ceiling(max(abs(x_limits))))) +
    
    # Adjust size of symbols in legend. Since, legend is based on color, we use color = 
    guides(colour = guide_legend(override.aes = list(size = 3)),
           shape = guide_legend(override.aes = list(size = 3))) +
    
    # Define the color of the dots
    #scale_color_viridis_d()+
    scale_color_manual(values = volcano_palette) +
    
    # Adjust font size, style
    my_theme +
    
    # Add gene labels
    geom_text_repel(data = volcano_df %>% dplyr::filter(Gene %in% disp_genes, padj < padj_cutoff),
                    mapping = aes(label = Gene),
                    #size = 2,
                    force = 0.5,
                    point.size = 1,
                    angle = 0,
                    #vjust = 0,
                    #hjust = 0,
                    #direction = "y",
                    box.padding = 1,  # increases line length somehow
                    point.padding = 0.1,
                    max.overlaps = Inf,
                    xlim = c(NA, NA),
                    ylim = c(-Inf,NA),
                    min.segment.length = dplyr::if_else(draw_line == "TRUE", 0.2, Inf),
                    #arrow = arrow(length = unit(0.015, "npc")),
                    position = position_quasirandom())
  
  # Save the plot
  ggplot2::ggsave(filename =  paste0("Volcano_Plot",  file_suffix, ".pdf"),
                  plot = last_plot(),
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  #width = 8.5,
                  #height = 9,
                  units = c("in"),	 
                  dpi = 300,
                  limitsize = TRUE,
                  bg = "white")
}

#******************************************************************************#
#                                 HEATMAP PLOT                                 #
#******************************************************************************#

# Function to plot heatmap
# NOTE: First column of normalized_counts MUST be "SYMBOL"
# columns parameter defines which column of metadata will be plotted as columns
plot_heatmap <- function(normalized_counts, metadata_column, metadata_row, plot_genes, disp_genes, file_suffix){
  
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
    # Make sure all values are numerics
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
  # mat[, unlist(lapply(mat, is.numeric))]    #alternative to mutate(across())
  
  # Perform log transform if needed. Count data is usually skewed right or left.
  # So, it will mostly be red or blue. log transform to make it less skewed.
  if (already_log == FALSE){
    mat <- log(1+mat, base = 2)
  }
  
  # Perform scaling for each gene across samples/cells if needed. Since scale() 
  # performs only column scaling, we transpose the dataframe first, so that 
  # genes are on columns & cells are on rows & then perform scaling.
  if (already_scaled == FALSE){
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
  mat <- mat[, metadata_column[,columns]]
  
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
    row_annotation <- dplyr::left_join(mat %>% as.data.frame() %>% tibble::rownames_to_column("SYMBOL"),
                                       metadata_row,
                                       by=c("SYMBOL"="SYMBOL")) %>%
      dplyr::select(everything(), -colnames(mat)) %>%
      tibble::column_to_rownames("SYMBOL")
  }
  
  #****************************************************************************#
  # Define colors column annotation
  #****************************************************************************#
  
  groups <- unique(unlist(col_annotation, use.names=FALSE))
  
  colors <- c("#BF812D", "#35978F", "#C51B7D", "#7FBC41", "#762A83",
              "#E08214", "#542788", "#D6604D", "#4393C3", "#878787",
              "#1A1A1A", "#FFFFBF", "#377EB8", "#4DAF4A", "#984EA3",
              "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
              "#66C2A5", "#FC8D62", "#000000", "#9E0142", "#E41A1C")
  #colors <- c("#74006F", "#FFC606")    # Female: Purple, Male:Gold
  #colors <- c("#F6D2E0", "#C8E7F5")    # Female: Pink,   Male:Blue (BBN paper)
  colors <- colors[1:length(groups)]
  names(colors) <- groups
  
  ann_colors <- list(colors[1:lengths(lapply(col_annotation, unique))[[1]]])
  
  x <- 2
  while (x <= ncol(col_annotation)){ 
    ann_colors <- c(ann_colors, 
                    list(colors[(sum(lengths(ann_colors))+1):(sum(lengths(ann_colors))+lengths(lapply(col_annotation, unique))[[x]])]))
    x <- x+1
  }
  names(ann_colors) <- colnames(col_annotation)
  
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
  # Define how samples will be arranged in the heatmap
  #****************************************************************************#
  
  if (row_clustering_alphabetical == TRUE){
    mat <- mat[sort(rownames(mat)),]
  } 
  if (col_clustering_alphabetical == TRUE){
    mat <- mat[,sort(colnames(mat))]
  } 
  #****************************************************************************#
  # Determine breaks for heatmap color scale.
  # breaks correspond to numerical ranges for color palette's bins 
  # i.e. 0 to length(my_palette)
  #****************************************************************************#
  
  if(max(mat) == 0){
    breaks <- c(seq(from = min(mat), to = 0, length.out = 100))
    my_palette <- my_palette[1:100]
  } else if (min(mat) == 0){
    breaks <- c(seq(from = 0, to = max(mat), length.out = 100))
    my_palette <- my_palette[100:200]
  } else if(min(mat) < -3 | max(mat) > 3){
    breaks <- c(seq(-3, 0, length.out = 50), seq(3/100, 3, length.out = 50))
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
      if (!(col_annotation[,anno_columns][i] %in% element_names)){
        element_names <- c(element_names, col_annotation[,anno_columns][i])
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
  # Perform row and column clustering
  #****************************************************************************#
  
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
  
  #****************************************************************************#
  # List genes and samples you want to display in the plot
  #****************************************************************************#
  
  display_col <- colnames(reordered)
  display_row <- data.frame("gene" = rownames(reordered)) %>%
    dplyr::mutate(gene = dplyr::case_when(gene %in% disp_genes ~ gene, 
                                          TRUE ~ "")) %>% 
    unlist(use.names=FALSE)
  
  # Save the clustered scores in xlsx
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Heatmap_matrix")
  openxlsx::writeData(wb, sheet = "Heatmap_matrix", x = reordered, rowNames = TRUE)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "Heatmap_matrix", file_suffix, ".xlsx"), 
                         overwrite = TRUE)
  
  #****************************************************************************#
  # Plot heatmap
  #****************************************************************************#
  pheatmap::pheatmap(mat = as.matrix(reordered),
                     color = my_palette,
                     breaks = breaks, 
                     border_color = "grey90", #"white"
                     cellwidth = dplyr::case_when(ncol(reordered) > 25 ~ 1,
                                                  ncol(reordered) <= 25 ~ 6), 
                     cellheight = dplyr::case_when(nrow(reordered) > 25 ~ 1,
                                                   nrow(reordered) <= 25 ~ 6), 
                     scale = "none",   
                     cluster_rows = row_clustering,   #cluster the rows
                     cluster_cols = col_clustering,   #cluster the columns
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method = "complete",
                     legend = TRUE, 
                     legend_breaks = NA,
                     legend_labels = NA, 
                     annotation_row = row_annotation,
                     annotation_col = col_annotation,
                     annotation_colors = ann_colors,
                     annotation_legend = TRUE,
                     annotation_names_row = FALSE,
                     annotation_names_col = TRUE,
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
                     width = 8.5,
                     height = 11,
                     filename = paste0(results_path, "Heatmap", file_suffix, ".jpg"))
}
