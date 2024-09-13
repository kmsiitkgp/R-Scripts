#plot_pca

# Note: When you use autoplot, each dot represents the rows of read_data.
# So, if plotting read_data from DESeq2, use tranposed read_data like t(read_data)
# Also, remove genes which have zero expression in all samples to avoid error messages
# Also, use normalized counts not raw counts to replicate PCA plot from DESeq2

iris.pca <- prcomp(data, 
                   center = TRUE, 
                   scale. = TRUE) 


library(ggfortify) 
iris.pca.plot <- autoplot(iris.pca, label=TRUE)








# This R file contains all the user defined functions which can be imported
# to analyze bulk RNA Seq, single cell RNA Seq, plot graphs etc

# Assign colors for UMAP. The current palette supports up to 33 cell types
my_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", 
                "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                "#FFC61E", "#762A83", "#333333", "#FF1F5B", "#B8E80C",
                "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
                "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5",
                "#FFFFBF", "#C51B7D", "#DE77AE", "#F1B6DA", "#FDE0EF", 
                "#F7F7F7", "#E6F5D0", "#B8E186")

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
    ggplot2::scale_fill_manual(values = rep(my_palette, times = ceiling(length(get(plot_param$title_legend_fill))/length(my_palette)))) +
    
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

