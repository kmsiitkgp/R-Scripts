# https://rpubs.com/bman/infercnvExampleOct222018
# NOTE: Your seurat object MUST have noth normal and tumor samples so that
# we can define the normal cells as reference for each cell type. If your 
# seurat object only has tumor samples, use the lymphoid and myeloid cells as
# normal cells

#******************************************************************************#
#                       STEP 1: GET ENSEMBL ANNOTATIONS                        #
#******************************************************************************#

get_annotations <- function(database){
  
  #**************************GET ENSEMBL ANNOTATIONS***************************#
  # Connect to AnnotationHub
  ah <- AnnotationHub()
  
  # Access the Ensembl database for organism
  ahDb <- AnnotationHub::query(ah,
                               pattern = c(species, "EnsDb"),
                               ignore.case = TRUE)
  
  # Acquire the latest annotation files
  id <- ahDb %>%
    mcols() %>%
    rownames() %>%
    tail(n = 1)
  
  # Download the appropriate Ensembldb database
  edb <- ah[[id]]
  
  # Extract gene-level information from database
  ensembl <- ensembldb::genes(x = edb,
                              return.type = "data.frame")
  
  # Select annotations of interest
  ensembl <- ensembl %>%
    dplyr::rename(chr = seq_name, SYMBOL = gene_name, ID = gene_id) %>%
    #dplyr::select(ID, SYMBOL, chr, description) %>%
    dplyr::distinct(ID, .keep_all = TRUE)
  
  #***************************GET ENTREZ ANNOTATIONS***************************# 
  
  # # Map Entrez id with "ENSEMBL" or "SYMBOL" whichever maps highest genes
  # # Here we map with SYMBOL
  # mapping <- AnnotationDbi::mapIds(x = org.Hs.eg.db, keys = read_data$ENTREZID,
  #                                  keytype = "ENTREZID", column = "SYMBOL")
  # 
  # mapping <- as.data.frame(do.call(cbind, list(mapping))) %>%
  #   tibble::rownames_to_column("ENTREZID")
  # 
  # colnames(mapping) <- c("ENTREZID","GENE")
  # 
  # # Merge read_data with mapped data
  # read_data <- read_data %>%
  #   dplyr::left_join(mapping, by=c("ENTREZID"="ENTREZID")) %>%
  #   dplyr::select(GENE, everything(), -ENTREZID)
  # 
  # # Check how many ENTREZIDs couldn't be mapped
  # unmapped <- nrow(read_data %>% dplyr::filter(is.na(.[[1]])))
  # cat("Genes in total :", nrow(read_data), "\t")
  # cat("Genes mapped   :", nrow(read_data)-unmapped, "\t")
  # cat("Genes unmapped :", unmapped, "\t")
  
  if (species == "Homo sapiens"){
    entrez <- AnnotationDbi::select(x = org.Hs.eg.db, keys = keys(org.Hs.eg.db),
                                    keytype = "ENTREZID", column = "SYMBOL")
  } else{
    entrez <- AnnotationDbi::select(x = org.Mm.eg.db, keys = keys(org.Mm.eg.db),
                                    keytype = "ENTREZID", column = "SYMBOL")
  }
  colnames(entrez) <- c("ID", "SYMBOL")
  
  # Choose which annotations to use
  if(database == "Ensembl"){
    annotations <- ensembl
  } else{
    annotations <- entrez
  }
  
  return(annotations)
}

#******************************************************************************#
#                             STEP 2: RUN INFERCNV                             #
#******************************************************************************#

# Load the integrated seurat object
celltype <- NULL
integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn", 
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))

integrated_seurat <- subset(x = integrated_seurat,
                            cell_class %in% c("Unclassified", "Mixed"),
                            invert = TRUE)

# Create gene order file
ann <- get_annotations(database)

gene_order <- ann %>%
  dplyr::rename(start = gene_seq_start, stop = gene_seq_end) %>% 
  dplyr::select(SYMBOL, chr, start, stop) %>%
  dplyr::filter(SYMBOL %in% rownames(integrated_seurat@assays$RNA@data)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  tibble::column_to_rownames("SYMBOL")

# Create annotation file
annotation <- integrated_seurat@meta.data %>% 
  dplyr::mutate(cell_class = paste0(Condition, "_", cell_class)) %>%
  dplyr::select(Cell, cell_class) %>% 
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("Cell")

# Define normal cells
ref <- unique(annotation$cell_class[grepl(pattern = "Normal", x = annotation$cell_class)])

# # If your seurat object doesnt have normal samples
# annotation <- integrated_seurat@meta.data %>% 
#   dplyr::select(Cell, cell_class) %>% 
#   tibble::remove_rownames() %>%
#   tibble::column_to_rownames("Cell")
# ref <- setdiff(unique(annotation$cell_class), c("Epithelial", "Fibroblasts"))

# Create infercnv object
# infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = "../inst/extdata/oligodendroglioma_expression_downsampled.counts.matrix.gz",
#                                                annotations_file="../inst/extdata/oligodendroglioma_annotations_downsampled.txt",
#                                                delim="\t",
#                                                gene_order_file="../inst/extdata/gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt",
#                                                ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)"))

infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = integrated_seurat@assays$RNA@counts,
                                               annotations_file = annotation,
                                               delim = "\t",
                                               gene_order_file = gene_order,
                                               ref_group_names = ref)

# Run with default parameters. This step does log normalization, anscombe 
# transformation etc etc...so use raw counts.
# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
infercnv_obj_default <- infercnv::run(infercnv_obj = infercnv_obj,
                                      cutoff = 0.1,  
                                      out_dir = seurat_results,
                                      cluster_by_groups = TRUE,
                                      denoise = TRUE,
                                      png_res= 600,
                                      no_prelim_plot = TRUE)

infercnv::plot_cnv(infercnv_obj = obj,
                   out_dir = seurat_results,
                   title = "inferCNV",
                   obs_title = "Observations (Cells)",
                   ref_title = "References (Cells)",
                   cluster_by_groups = TRUE,
                   cluster_references = TRUE,
                   plot_chr_scale = FALSE,
                   output_filename = "infercnv",
                   output_format = "pdf",
                   png_res = 300)