# Initially, I analyzed Sandrine and Supriya samples together. However, 
# not many DEGs were identified when comparing groups within Supriya dataset
# like TQuadRFL vs TQuadFC. However, if I analyze Supriya dataset separately
# from Sandrine dataset, I see lot of DEGs. As Supriya dataset was generated
# at a different time point than Sandrine dataset, I decided to perform their
# analysis separately.

# ---- User Override Project Directories & Parameters ----

proj <- "Sandrine_Quad"
species <- "Mus musculus"
contrasts <- c("Quadriceps.RF-Quadriceps.FC",
               "Quadriceps.RFL-Quadriceps.FC")

parent.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data"
gmt.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets"

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
  col.ann          = c("Treatment"),                  # Column annotation
  #row.ann          = NULL,                  # Row annotation
  #col.gaps         = NULL,                  # Column gaps
  #row.gaps         = NULL,                  # Row gaps
  col.cluster      = c("Treatment"),                 # Column clustering
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


# Remove bad samples [CASE by CASE BASIS]
read_data <- read_data %>% 
  dplyr::select(everything(), -c("SBQuadFc2", "SBQuadFc4"))






cGAS_STING_human <- c(
  # Sensors / Adapters
  "CGAS", "STING1", "TREX1", "SMPDL3A", "TAB1",
  # Kinases / Signaling
  "TBK1", "IKBKA", "IKBKB", "IKBKG", "JAK1", "TYK2",
  # Transcription Factors
  "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF7", "IRF9", "STAT6", "NFKB1", "NFKB2",
  # Inflammasome
  "NLRP1", "NLRP3", "NLRP6", "NLRP12", "ASC", "CASP1",
  # Type I IFNs & Receptors
  "IFNB1", "IFNA1", "IFNAR1", "IFNAR2",
  # Chemokines
  "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL11", "CCL13", "CCL17", "CXCL10",
  # Other downstream / reported genes
  "CYP1A1", "CYP1B1"
)

df <- read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Sandrine_Quad/Quadriceps.RFL-Quadriceps.FC/DEGs.xlsx")
df <- df %>% dplyr::filter(tolower(SYMBOL) %in% tolower(cGAS_STING_human), padj <= 0.05)
