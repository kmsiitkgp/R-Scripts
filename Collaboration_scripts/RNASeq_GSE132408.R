
proj <- "GSE132408"
species <- "Homo sapiens"
contrasts <- c("TET2KO-WT",
               "TET2KO.IFNy-WT.IFNy")

# DESeq2 overrides
deseq2.override <- list(
  contrasts     = contrasts ,               # Vector of contrasts for DE analysis
  design        = "Comparisons",            # DESeq2 design formula or column name
  #lfc_cutoff    = 0,                        # Log fold change cutoff for significance
  #padj_cutoff   = 0.1,                      # Adjusted p-value cutoff for significance
  batch_correct = FALSE                     # Boolean, whether to apply batch correction
)

# Heatmap overrides
heatmap.override <- list(
  col_annotations    = c("Condition", "Treatment"),        # NULL, 1 or more columns from metadata_col for column annotation
  #row_annotations    = NULL,        # NULL, 1 or more columns from metadata_row for row annotation
  #col_gap_by         = NULL,        # NULL, 1 column from metadata_col to define column gaps in heatmap
  #row_gap_by         = NULL,        # NULL, 1 column from metadata_row to define row gaps in heatmap
  col_cluster_by     = "allx",       # NULL, 1 column from metadata_col, "all", "alphabetical" for clustering columns
  row_cluster_by     = "all",       # NULL, 1 column from metadata_row, "all", "alphabetical" for clustering columns
  #plot_title         = NULL,        # NULL, Title for heatmap (default NULL i.e. no title)
  #heatmap_palette    = "rdbu",      # Color palette for heatmap matrix ("rdbu" or "vrds")
  #annotation_palette = "discrete",  # Color palette for heatmap annotation ("discrete" or "continuous")
  #border_color       = NA,          # Color of heatmap cell borders (default NA i.e. no border)
  #force_log          = FALSE,       # Force log transform (default FALSE i.e. auto detect)
  #show_expr_legend   = TRUE,        # Show expression legend (set FALSE if annotations overlap)
  save_plot          = TRUE,       # Save the heatmap plot as pdf (default FALSE i.e. no save)
  save_matrix        = TRUE        # Save the heatmap matrix as xlsx (default FALSE i.e. no save)
)

# Volcano plot overrides
volcano.override <- list(
  #top_n        = 5,                  # If label_genes = NULL, label top_n genes (default 5 genes)
  #lfc_cutoff   = 0.58,               # Log fold change cutoff (default 0.58)
  padj_cutoff  = 0.05                # Adjusted p-value cutoff (default 0.05)
)

