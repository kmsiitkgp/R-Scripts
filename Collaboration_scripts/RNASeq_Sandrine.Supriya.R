# Initially, I analyzed Sandrine and Supriya samples together. However, 
# not many DEGs were identified when comparing groups within Supriya dataset
# like TQuadRFL vs TQuadFC. However, if I analyze Supriya dataset separately
# from Sandrine dataset, I see lot of DEGs. As Supriya dataset was generated
# at a different time point than Sandrine dataset, I decided to perform their
# analysis separately.

proj.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Sandrine_Quad"
gmt.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets"

proj.params <- list(
  proj                    = "Sandrine_Quad",               # Project name
  species                 = "Mus musculus",                # Species name (e.g. "Mus musculus", "Homo sapiens")
  counts.dir              = file.path(proj.dir, "counts"), # Directory containing count files
  contrast                = c("Quadriceps.RF-Quadriceps.FC",
                              "Quadriceps.RFL-Quadriceps.FC"),  # Vector of contrasts for DE analysis
  deseq2.design           = "Comparisons",                     # DESeq2 design formula or column name
  deseq2.lfc.cutoff       = 0,                                 # Log fold change cutoff for significance
  deseq2.padj.cutoff      = 0.1,                               # Adjusted p-value cutoff for significance
  deseq2.batch.correct    = FALSE,                             # Boolean, whether to apply batch correction
  
  heatmap.force.log       = FALSE,                       # Force log transform on heatmap data (default FALSE, auto detect)
  heatmap.col.ann         = c("Treatment"),              # Columns from metadata used as column annotation
  heatmap.row.ann         = NULL,                        # Columns from metadata_row for row annotation
  heatmap.col.gaps        = NULL,                        # Columns to define gaps in heatmap columns
  heatmap.row.gaps        = NULL,                        # Columns to define gaps in heatmap rows
  heatmap.col.cluster     = c("Treatment"),              # Clustering for columns ("all", "alphabetical", or metadata column)
  heatmap.row.cluster     = "all",                       # Clustering for rows ("all", "alphabetical", or metadata column)
  heatmap.palette         = "rdbu",                      # Color palette for heatmap ("rdbu" or "vrds")
  heatmap.ann.palette     = "discrete",                  # Annotation palette type ("discrete" or "sequential")
  heatmap.border.color    = NA,                          # Border color of heatmap cells, NA for no border
  heatmap.show.expr.legend= TRUE,                        # Show expression legend on heatmap (set FALSE if overlapping annotations)
  heatmap.title           = NA,                          # Title for heatmap (default NA = no title)
  heatmap.format          = "tiff",                      # Output format for heatmap ("tiff", "jpeg", etc.)
  
  volcano.lfc.cutoff      = 0.58,                         # Log fold change cutoff for volcano plot
  volcano.padj.cutoff     = 0.05,                         # Adjusted p-value cutoff for volcano plot
  volcano.color           = "vrds",                       # Color palette for volcano plot ("vrds", etc.)
  volcano.label.genes     = NULL                          # Optional vector of genes to label on volcano plot
)

# Remove bad samples [CASE by CASE BASIS]
meta_data <- read.xlsx(file.path(proj.dir, paste0(proj.params$proj, ".Metadata.xlsx")))
read_data <- read.xlsx(file.path(proj.dir, paste0(proj.params$proj, ".raw.counts.xlsx")))
read_data <- read_data %>% 
  dplyr::select(everything(), -c("SBQuadFc2", "SBQuadFc4"))


# # Some samples need to be removed in Sandrine's dataset based on PCA
# # NOTE: Although TQuadFc1 and TQuadFc2 are similar to SBQuadFc2 and SBQuadFc4
# # based on pCA plot, we DO NOT remove them because we are comparing TQuadFc1 and
# # TquadFc2 with TQuadRFL1 and TQuadRFL2
# read_data <- read_data %>% 
#   dplyr::select(everything(), -c("SBQuadFc2", "SBQuadFc4"))
# 
# DEG.params  <- list(contrast    = c("Quadriceps.RFL-Quadriceps.FC",
#                                     "Heart.RFL-Heart.FC"),
#                     design      = "Comparisons",
#                     design.ref  = c("Treatment:FC"),
#                     lfc.cutoff  = 0,
#                     padj.cutoff = 0.1,
#                     deseq2.batch.correct = FALSE,
#                     proj        = "RNASeq_Supriya",
#                     species     = "Mus musculus")
# 
# heatmap.params <- list(anno.row       = NULL,        # NULL, c("Group")
#                        anno.column    = c("Source", "Tissue", "Treatment"),
#                        row.split      = NA,     
#                        col.split      = NA, #c("Treatment"),
#                        row.cluster    = c("all"),           # c("alphabetical", "group", "all")
#                        col.cluster    = c("all"),  # c("alphabetical", "group", "all")
#                        discrete_panel = FALSE, 
#                        log.transform  = TRUE,
#                        scale          = TRUE,
#                        border_color   = "white",
#                        bar_width      = NA,              # NA , 5
#                        bar_height     = NA,              # NA , 5
#                        width          = 5,              # NA
#                        height         = 5,              # NA 
#                        matrix_color   = "rdbu",          # c("vrds", "rdbu")
#                        expr_legend    = TRUE,  # set FALSE if it overlaps with annotation legends
#                        file_format    = "tiff")
# 
# 
# data_path1 <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Sandrine/"
# data_path2 <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/RNASeq_Supriya/"
# output_path <- data_path1
# 
# file_suffix <- "Supriya.vs.Sandrine.RFL.vs.FC"
# df1 <- read.xlsx(paste0(data_path1, "DEGs.Quadriceps.RFL-Quadriceps.FC_.xlsx"))
# df2 <- read.xlsx(paste0(data_path2, "DEGs.Quadriceps.RFL-Quadriceps.FC_.xlsx"))
# 
# file_suffix <- "Supriya.vs.Sandrine.RFL.vs.FC"
# df1 <- read.xlsx(paste0(data_path1, "DEGs.Quadriceps.RFL-Quadriceps.FC_.xlsx"))
# df2 <- read.xlsx(paste0(data_path2, "DEGs.Quadriceps.RFL-Quadriceps.FC_.xlsx"))