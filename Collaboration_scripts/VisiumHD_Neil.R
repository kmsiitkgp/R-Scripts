proj <- "VisiumHD_Bhowmick"
species <- "Homo sapiens"
contrasts <- c("Metastatic.Tumor.Pancreas.Met.Panc-Primary.Tumor.Pancreas.Primary.Panc",
               "Metastatic.Tumor.Pancreas.Macrophages-Normal.Liver.Macrophages") 

parent.dir <- "/hpc/home/kailasamms/scratch"
gmt.dir <- "/hpc/home/kailasamms/projects/GSEA_genesets"

# DESeq2 overrides
deseq2.override <- list(
  contrasts     = contrasts,
  #design        = "Comparisons",            # DESeq2 design formula or column name
  #lfc.cutoff    = 0,                        # Log fold change cutoff for significance
  #padj.cutoff   = 0.1,                      # Adjusted p-value cutoff for significance
  batch.correct = FALSE                     # Boolean, whether to apply batch correction
)

# Heatmap overrides
heatmap.override <- list(
  #force.log        = TRUE,                  # Force log transformation
  col.ann          = c("Treatment"),         # Column annotation
  #row.ann          = NULL,                  # Row annotation
  #col.gaps         = NULL,                  # Column gaps
  #row.gaps         = NULL,                  # Row gaps
  col.cluster      = c("Treatment"),         # Column clustering
  #row.cluster      = "all",                 # Row clustering
  #palette         = "rdbu",                 # Heatmap palette
  #ann.palette     = "discrete",             # Annotation palette
  #border.color    = NA,                     # Cell border color
  #show.expr.legend = TRUE,                  # Show expression legend
  #title           = "",                     # Heatmap title
  format           = "tiff"                  # Output file format
)

# Volcano plot overrides
volcano.override <- list(
  #lfc.cutoff  = 0.58,                         # Minimum log2 fold-change to highlight
  #padj.cutoff = 0.05,                      # Adjusted p-value cutoff
  #color       = "vrds",                    # Color palette
  #label.genes = c()                         # Genes to label on the plot
)



















for (n in seq_along(contrasts)) {
  
  
  meta_data <- read.xlsx(file.path(proj.params$proj.dir, paste0(proj, "_Metadata.xlsx")))
  DEGs_df <- openxlsx::read.xlsx(file.path(proj.params$proj.dir, "geomx/0.1, 0.1/lmer_results.xlsx"),
                                 sheet = contrasts[n])
  
  vst_counts <- openxlsx::read.xlsx(file.path(proj.params$proj.dir, "geomx/0.1, 0.1/Q3_Norm_Counts.xlsx")) %>%
    tibble::column_to_rownames("SYMBOL")
  
  contrast <- contrasts[n]
  contrast.dir <- proj.params$contrast.dir[n]
  
  if (!dir.exists(contrast.dir)) dir.create(contrast.dir, recursive = TRUE)
  
  output_path <- proj.params$deseq2.dir[n]
  plot_volcano_module(DEGs_df, proj.params, contrast, output_path)
  
  output_path <- proj.params$deseq2.dir[n]
  samples <- filter_samples_by_contrast(meta_data, contrast)
  samples <- intersect(colnames(vst_counts), samples)
  sig_genes <- DEGs_df %>% dplyr::filter(padj <= 0.05) %>% dplyr::pull(SYMBOL)
  norm_counts <- vst_counts[sig_genes, samples, drop = FALSE]
  metadata_col <- meta_data %>% dplyr::filter(Sample_ID %in% samples)
  metadata_row <- NULL
  disp_genes <- c()
  if (length(sig_genes) != 0){
    
    ph <- plot_heatmap_module(norm_counts, proj.params, metadata_col, metadata_row, disp_genes)
    
    pdf(file = file.path(output_path, "Heatmap.pdf"),
        width = 10, 
        height = 10)
    gridExtra::grid.arrange(grobs = list(ph$ph$gtable), ncol = 1)
    dev.off()
    
    save_xlsx(ph$mat, file.path(output_path, "Heatmap_Matrix.xlsx"), "Heatmap_matrix", row_names = TRUE)
  }
  
  # ---- Pathway Analysis ----
  
  # Perform pathway analysis
  output_path <- proj.params$pathway.dir[n]
  gsea_res_list <- pathway_analysis_module(DEGs_df, proj.params, output_path)
  
  # Plot Bar, Dot, Heatmap for top 10 Up & Down Pathways
  top_pathways_gsea <- gsea_res_list$consensus %>%
    dplyr::group_by(Collection, Consensus, Description) %>%
    dplyr::slice_min(order_by = match(method, c("FGSEA", "GSEA", "ORA")), n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Collection, Consensus) %>%
    dplyr::slice_max(order_by = abs(NES), n = 10, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  top_pathways_ora <- gsea_res_list$consensus %>%
    dplyr::filter(method == "ORA") %>%
    dplyr::group_by(Collection, Consensus) %>%
    dplyr::slice_min(order_by = padj, n = 10, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  output_path <- proj.params$pathway.dir[n]
  samples <- filter_samples_by_contrast(meta_data, contrast)
  samples <- intersect(colnames(vst_counts), samples)
  plot_pathway_module(top_pathways_gsea, vst_counts, meta_data, samples, "GSEA", output_path)
  plot_pathway_module(top_pathways_ora, vst_counts, meta_data, samples, "ORA", output_path)
  
  # ---- Transcription Factor (TF) Analysis ----
  
  # Perform TF analysis on DEGs using t-statistics
  t_stats_mat <- DEGs_df %>% 
    as.data.frame() %>%
    dplyr::mutate(t = -log10(padj) * log2FoldChange) %>%
    dplyr::filter(!is.na(t)) %>%
    dplyr::select(SYMBOL, t) %>%
    tibble::column_to_rownames("SYMBOL") %>%
    as.matrix()
  
  tf_res_degs <- tf_analysis_module(t_stats_mat, proj.params$species)
  
  # Perform TF analysis on vst counts 
  samples <- filter_samples_by_contrast(meta_data, contrast)
  samples <- intersect(colnames(vst_counts), samples)
  norm_counts_sub <- vst_counts[, samples] %>% as.matrix()
  
  tf_res_counts <- tf_analysis_module(norm_counts_sub, proj.params$species)
  
  # Plot Bar, Heatmap for top 20 Up and Down TFs
  contrast <- contrasts[n]
  output_path <- proj.params$tf.dir[n]
  samples <- filter_samples_by_contrast(meta_data, contrast)
  samples <- intersect(colnames(vst_counts), samples)
  plot_tf_module(tf_res_degs, contrast, meta_data, samples, output_path, n_tfs = 20)
  plot_tf_module(tf_res_counts, contrast, meta_data, samples, output_path, n_tfs = 20)
}
