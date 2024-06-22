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

#******************************************************************************#
#                                 VOLCANO PLOT                                 #
#******************************************************************************#

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
  
  # NOTE: DO NOT USE labels for defining colors due to reasons below. 
  # RECOMMEND using a named vector.
  # NOTE: If using labels, sort labels in alphabetical order and then assign 
  # color because R by default will arrange the labels in alphabetical order 
  # first and then match them to colors indicated in values vector and then 
  # color the plot. The coloring in the legend is however dependent on the 
  # order of labels vector and values vector. To understand, create a plot first 
  # using the sort and then without the sort(). 
  
  
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
                      #color = padj,
                      fill = log2FC,
                      size = abs(log2FC)*-log10(padj),
                      shape = Direction)) +
    
    # Plot dot plot
    ggplot2::geom_point(col="black", 
                        shape=21,
                        stroke=0.5,
                        position=position_jitter(h=0.01,w=0.01)) +
    
    # Define the theme of plot
    ggplot2::theme_classic() +
    
    ggplot2::guides(colour = guide_legend(override.aes = list(size = 3)),
                    shape = guide_legend(override.aes = list(size = 3)),
                    size = "none",
                    fill = guide_colourbar(theme = theme(legend.key.width  = unit(0.75, "lines"),
                                                         legend.key.height = unit(10, "lines"),
                                                         legend.ticks = element_blank(),
                                                         legend.frame = element_rect(colour = "Black",
                                                                                     linewidth = 1)))) +
    
    # Define the color of the dots
    scale_fill_gradient2(low="#007ba7", mid="Black", high="Yellow") +
    
    # Define the axis, plot headings
    ggplot2::labs(x = expression("log"[2]*"FC"),
                  y = expression("-log"[10]*"padj"),
                  color = "Direction",
                  fill = "log2FC",
                  size = "log2FC",
                  shape = "Direction",
                  color = color_by,
                  title = paste0(Target, " vs ", Reference)) +   
    
    # Draw line to mark the cutoffs
    geom_vline(xintercept = c(-log2_cutoff,log2_cutoff), color = "black", linetype = "dotted", linewidth = 0.5) +
    geom_hline(yintercept = -log10(padj_cutoff), color = "black", linetype = "dotted", linewidth = 0.5) +
    
    # Define x-axis start and end
    coord_cartesian(xlim = c(-ceiling(max(abs(x_limits))), ceiling(max(abs(x_limits))))) +
    
    # Adjust font size, style
    my_theme 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ggplot2::ggplot(data = volcano_df, 
                  aes(x = log2FC, 
                      y = -log10(padj), 
                      label = Gene,
                      color = get(color_by), 
                      shape = Direction)) +
    
    # Plot dot plot
    ggplot2::geom_point(size=0.2) +
    
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