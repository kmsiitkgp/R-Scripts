#!/usr/bin/env Rscript

#**********************STEP 6: MAKE PIE/STACKED BAR CHART**********************#

# Read proteome file if you already prepared it
proteome <- openxlsx::read.xlsx(xlsxFile = paste0(results_path, "Proteome/HPA_Proteome.xlsx"))

# Create list of cells corresponding to import file names
cells <- c("Epithelial_Cells", "Non-Epithelial_Cells", "Fibroblasts", "Myeloid_Cells","B_Cells", "T_Cells" )

# Plot pie chart
for (cell in cells[1:2]){
  
  # Read DEGs file
  DEGs <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "new figs BBN/DEGs_id_Male_vs_Female_", cell, ".xlsx"))
  # DEGs <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "DEGs_id_BBN_vs_Normal_Epithelial Cells.xlsx"))
  # DEGs <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "DEGs_id_BBN_vs_Normal_Fibroblasts.xlsx"))
  DEGs <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "/Fibroblast/", "DEGs_Fibroblast_proj.xlsx"))
  
  # If any of the column is missing column names, you will see
  # "Error in initialize(...):attempt to use zero-length variable name"
  DEGs <- DEGs %>% 
    dplyr::filter(padj < 0.05 & base::abs(log2FoldChange) >= 0.58) %>%
    dplyr::left_join(proteome, by = c("gene_name" = "mouse_gene_name")) %>%
    dplyr::mutate(Class = dplyr::if_else(is.na(Class), "Unknown", Class)) %>%
    dplyr::mutate(regulation = dplyr::if_else(log2FoldChange > 0, "up", "down", "NA"))
  
  # Check which of the secreted proteins are plasma proteins etc
  secreted_DEGs <- DEGs %>% 
    dplyr::filter(Class == "Secreted" | Class == "Secreted & Membrane" | Class == "Membrane") 
  
  secreted_DEGs <- dplyr::left_join(secreted_DEGs, 
                                    secreted_class %>% select("Ensembl", "X.Protein.class."),
                                    by=c("Ensembl"="Ensembl"))
  
  
  
  # # Plot pie chart for male and female
  # for (direction in c("up", "down")){
  # 
  #   # Calculate % of each cell type for sample being analyzed
  #   data <- DEGs %>%
  #     dplyr::filter(regulation == direction) %>%
  #     dplyr::count(class) %>%
  #     dplyr::mutate(percent = round(100*n/sum(n, na.rm=TRUE), digits = 0), label_percent = paste0(percent,"%")) %>%
  #     dplyr::arrange(class)
  #   
  #   # If you run ggplot without themevoid(), you will see 0 and 100% dont overlap.
  #   # We need to add the percentage value of 1st cell group  to get accurate label positions
  #   # data <- data %>%
  #   #   mutate(ypos = 100 + data$percent[1]-(cumsum(percent)-percent/2))
  #   
  #   ggplot2::ggplot(data = data, aes(x = "", y = percent, fill = class)) +
  #     geom_bar(stat = "identity", width = 3, color = "white") +
  #     coord_polar(theta = "y", start = 0, direction = -1) +
  #     #geom_label(aes(x = 1.6, y = ypos, label = label_percent), size = 5, label.size = NA, fill = NA) +
  #     #geom_label(aes(x = 2, label = label_percent), position = position_stack(vjust = 0.5), color = "black", size = 5, check_overlap = TRUE) +
  #     geom_text(aes(x=3, label=label_percent), position=position_stack(vjust=0.5), fontface="bold", colour="black", size=5, check_overlap=TRUE) +
  #     theme_void() +        #remove background, grid, numeric labels
  #     scale_fill_brewer(palette = "Set1",
  #                       aesthetics = "fill") +
  #     ggplot2::labs(title = dplyr::if_else(direction == "up", "Male", "Female"),
  #                   fill = "Protein Class",
  #                   x = "",
  #                   y = "") +
  #     ggplot2::theme(plot.title =   element_text(family="sans", face="bold",  colour="black", size=15, hjust = 0.5),
  #                    plot.caption = element_text(family="sans", face="bold",  colour="black", size=10, hjust = 0),
  #                    axis.title.x = element_text(family="sans", face="bold",  colour="black", size=15, hjust = 0.5),
  #                    axis.title.y = element_text(family="sans", face="bold",  colour="black", size=15, hjust = 0.5),
  #                    axis.text.x =  element_blank(),
  #                    axis.text.y =  element_blank(),
  #                    legend.title = element_text(family="sans", face="bold",  colour="black", size=15, hjust = 0),
  #                    legend.text =  element_text(family="sans", face="bold",  colour="black", size=12, hjust = 0.5),
  #                    strip.text.x = element_text(family="sans", face="bold",  colour="black", size=15, hjust = 0.5),
  #                    legend.position = "right",
  #                    legend.direction = "vertical",
  #                    legend.text.align = 0,
  #                    axis.line=element_blank(),
  #                    axis.ticks=element_blank())
  #   
  #   # Save the plot
  #   ggplot2::ggsave(filename = paste0("Pie_chart_Proteome_Class_", cell, "_", direction, "regulated.pdf"),
  #                   plot = last_plot(),
  #                   device = "pdf",
  #                   path = paste0(parent_path, "Proteome/"),
  #                   scale = 1,
  #                   #width = 11,
  #                   #height = 8.5,
  #                   units = c("in"),
  #                   dpi = 600,
  #                   limitsize = TRUE,
  #                   bg = NULL)
  # }
}

# Plot stacked bar chart
data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "BBN original figs/DEGs_id_Male_vs_Female_", cells[1], ".xlsx"))
data <- data %>% dplyr::mutate(cell_type = cells[1])

for (cell in cells[2:2]){
  
  # Read DEGs file
  DEGs <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "BBN original figs/DEGs_id_Male_vs_Female_", cell, ".xlsx"))
  DEGs <- DEGs %>% dplyr::mutate(cell_type = cell)

  
  # Merge all data to be plotted
  data <- rbind(data, DEGs)

  # Calculate percentages
  data <- data %>%
    dplyr::filter(padj < 0.05 & base::abs(log2FoldChange) >= 0.58) %>%
    dplyr::left_join(proteome, by = c("gene_name" = "mouse_gene_name")) %>%
    dplyr::mutate(class = dplyr::if_else(is.na(Class), "Unknown", Class)) %>%
    dplyr::mutate(regulation = dplyr::if_else(log2FoldChange > 0, "up", "down", "NA")) %>%
    dplyr::group_by(cell_type, regulation) %>%
    dplyr::count(class) %>%
    dplyr::mutate(percent = round(100*n/sum(n, na.rm=TRUE), digits = 2), label_percent = paste0(percent,"%")) %>%
    dplyr::arrange(class) %>%
    dplyr::mutate(plot.x=paste0(stringr::str_trunc(cell_type, 2, "right", ellipsis = "") ,"_", stringr::str_trunc(regulation, 2, "right", ellipsis = "")))


  # Plot stacked bar chart
  ggplot2::ggplot(data = data, aes(x = plot.x, y = percent, fill = class)) +
    geom_bar(stat = "identity", position = "stack", width = 0.95) +
    theme_classic() + 
    scale_fill_brewer(palette = "Set1",
                      aesthetics = "fill") +
    geom_text(aes(x=plot.x, label=n), position=position_stack(vjust=0.5), fontface="bold", colour="white", size=6, check_overlap=TRUE) +
    ggplot2::labs(title = "",
                  fill = "Protein Class",
                  x = "",
                  y = "Percent composition")
  
}