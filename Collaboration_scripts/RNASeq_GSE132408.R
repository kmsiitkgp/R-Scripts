# In this case, the replicates for each group are clustered together.
# So, the within-group variability is low. Hence, we do not subset the data

proj <- "GSE132408"
species <- "Homo sapiens"
contrasts <- c("TET2KO-WT",
               "TET2KO.IFNy-WT.IFNy")

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
  col.ann          = c("Condition", "Treatment"),                  # Column annotation
  #row.ann          = NULL,                  # Row annotation
  #col.gaps         = NULL,                  # Column gaps
  #row.gaps         = NULL,                  # Row gaps
  col.cluster      = c("Condition"),                 # Column clustering
  #row.cluster      = "all",                 # Row clustering
  #palette         = "rdbu",                # Heatmap palette
  #ann.palette     = "discrete",            # Annotation palette
  #border.color    = NA,                    # Cell border color
  #show.expr.legend = TRUE,                  # Show expression legend
  #title           = "",                    # Heatmap title
  format           = "tiff"                 # Output file format
)

# Volcano plot overrides
volcano.override <- list(
  #lfc.cutoff  = 0.58,                         # Minimum log2 fold-change to highlight
  #padj.cutoff = 0.05,                      # Adjusted p-value cutoff
  #color       = "vrds",                    # Color palette
  #label.genes = c()                         # Genes to label on the plot
)

