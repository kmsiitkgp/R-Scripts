library("Seurat")
library("dplyr")                
library("tibble")
library("openxlsx")

for (proj in c("scRNASeq_Koltsova_sn", "scRNASeq_Koltsova")){
  
  # path containing rds file
  seurat_results <- paste0("/hpc/home/kailasamms/scratch/", proj, "/results_seurat/")
  
  # Read the processed rds file
  integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn.rds"))
  
  # Subset Male patients
  integrated_seurat <- subset(x = integrated_seurat,
                              subset = Sex == "Male")
  # Plot the UMAP
  gg1 <- Seurat::DimPlot(object=integrated_seurat,
                         reduction="umap.harmony",
                         cols=my_palette,
                         label=FALSE,
                         group.by="cell_type",
                         split.by=NULL, #split,
                         shape.by=NULL,
                         pt.size=0.2,
                         label.size=5,
                         repel=FALSE,
                         raster=FALSE) +
    ggplot2::labs(fill="CLUSTERS",
                  x="UMAP_1",
                  y="UMAP_2") +
    my_theme
  
  ggplot2::ggsave(filename=paste0("UMAP_", proj, ".pdf"),
                  plot=gg1,
                  device="pdf")
  
  # Protein coding Y genes
  y_genes <- c("Ddx3y", "Eif2s3y", "Kdm5d","Uty","Rbmy","Sly","Sry","Uba1y",
               "Usp9y","Zfy1","Zfy2","H2al2b","H2al2c","Orly","Rbm31y","Srsy",
               "Ssty1","Ssty2")
  
  # Filter out Y genes absent in our data
  features <- intersect(y_genes, rownames(integrated_seurat@assays$RNA$data))
  
  # Identify cell types in our data
  celltypes <- unique(integrated_seurat@meta.data$cell_type)
  celltypes <- celltypes[!is.na(celltypes)]
  
  # Create an empty dataframe where LOY % for each celltype will be stored
  final_df <- data.frame(Sample="",n_gene=0,n=0, Percent=0, celltype="")
  
  # Iterate through each cell type and calculate LOY %
  for (c in celltypes){
    
    # Subset cells of specific celltype
    seurat_obj <- subset(x=integrated_seurat,
                         cell_type == c)
    
    # Extract expression info of Y genes
    df <- seurat_obj@assays$RNA$data
    df <- df[features,]
    df <- df[rowSums(df) != 0,]
    
    # Convert count matrix to binary format [1=Expressed, 0=Not expressed]
    df[df>0] <- 1
    
    # calculate LOY%
    gene_count_per_cell <- data.frame(counts = colSums(df)) %>%
      tibble::rownames_to_column("Cell") %>%
      dplyr::mutate(Cell = gsub(pattern= "_.*" , replacement="", x=Cell)) %>%
      dplyr::group_by(Cell, counts) %>%
      dplyr::summarise(n = n()) %>%
      dplyr::rename(Sample=Cell, n_gene=counts) %>%
      dplyr::mutate(Percent = 100*n/sum(n)) %>%
      dplyr::filter(n_gene == 0) %>%
      dplyr::mutate(celltype = c)
    
    # Merge LOY% from all celltypes
    final_df <- dplyr::bind_rows(final_df, gene_count_per_cell)
  }
  
  # Remove dummy row and format the dataframe
  final_df <- final_df[-1,]
  final_df <- final_df %>% 
    dplyr::select(Sample, Percent, celltype) %>% 
    tidyr::pivot_wider(id_cols=celltype, names_from=Sample, values_from=Percent, values_fill = NA)
  
  # Save the clustered matrix
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = paste0("LOY Percent_", proj))
  openxlsx::writeData(wb, sheet = paste0("LOY Percent_", proj), x = final_df, rowNames = FALSE)
}
openxlsx::saveWorkbook(wb, file = "LOY Percent.xlsx", overwrite = TRUE)