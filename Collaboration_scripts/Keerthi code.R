df <- read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/excel valued chemokine.xlsx")
df <- df %>% tibble::column_to_rownames("Treatment") %>% scale() %>% t()

pheatmap(df)

heatmap.params <- list(anno.row       = NULL,        # NULL, c("Group")
                       anno.column    = NULL,
                       row.split      = NA,     
                       col.split      = NA,
                       row.cluster    = c("all"),  # c("alphabetical", "group", "all")
                       col.cluster    = c("all"),  # c("alphabetical", "group", "all")
                       discrete_panel = FALSE, 
                       log.transform  = FALSE,
                       scale          = FALSE,
                       border_color   = "white",
                       bar_width      = NA,              # NA , 5
                       bar_height     = NA,              # NA , 5
                       width          = 5,              # NA
                       height         = 5,              # NA 
                       angle          = 90,
                       matrix_color   = "rdbu",          # c("vrds", "rdbu")
                       expr_legend    = TRUE,  # set FALSE if it overlaps with annotation legends
                       file_format    = "tiff")

norm_counts <- df %>% data.frame() %>% tibble::rownames_to_column("SYMBOL") 
metadata_column <- data.frame("Sample.ID" = colnames(norm_counts))
metadata_row <- NULL
plot_genes <- rownames(df)
disp_genes <- plot_genes
file_suffix <- ""
output_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/"
plot_heatmap(norm_counts, metadata_column, metadata_row, heatmap.params,
                         plot_genes, disp_genes, file_suffix, output_path)
