gmt_dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets"

# ---- LOAD PACKAGES ----

pkgs <- c(
  "BiocManager", "remotes", "AnnotationHub", "ensembldb", "org.Hs.eg.db",
  "org.Mm.eg.db", "fgsea", "clusterProfiler", "DESeq2", "sva", "GSVA", 
  "RcisTarget", "glmGamPoi", "Seurat", "harmony", "hdf5r", "scCustomize", 
  "reticulate", "ashr", "infercnv", "UCell", "scDblFinder", "DropletUtils", 
  "CellChat", "SeuratWrappers", "presto", "DoubletFinder", "SeuratData", 
  "oligo", "oligoData", "illuminaHumanv4.db", "hgu133plus2.db", "GEOquery", 
  "affy", "lumi", "openxlsx", "dplyr", "tibble", "stringr", "purrr", "ggplot2",
  "ggplotify", "ggrepel", "ggpubr", "ggfortify", "ggridges", "ggbeeswarm",
  "pheatmap", "VennDiagram", "survival", "survminer", "UpSetR", "umap", 
  "plot3D", "cowplot", "viridis", "RColorBrewer", "colorspace", 
  "enrichplot", "ComplexHeatmap", "NanoStringNCTools", "GeomxTools", 
  "GeoMxWorkflows", "networkD3", "httr", "decoupleR", "OmnipathR", "SeuratDisk",
  "clustree", "crayon"
)

for (pkg in pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    message(paste("Loaded", pkg))
  } else {
    message(paste("Package", pkg, "is not installed — skipping"))
  }
}

# NOTE: survminer handles %++% while dplyr handles %>%

# ---- Custom Palette and ggplot Theme ----

custom_theme <- ggplot2::theme(
  plot.title    = element_text(family = "sans", face = "bold",  colour = "black", size = 15, hjust = 0.5),
  legend.title  = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0,   vjust = 1, angle = 0),
  axis.title.x  = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0.5, vjust = 0, angle = 0),
  axis.title.y  = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0.5, vjust = 1, angle = 90),
  legend.text   = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5),
  axis.text.x   = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5, vjust = 0.5, angle = 45),
  axis.text.y   = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5, vjust = 0.5, angle = 0)
)

custom_palette <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173", "#5254a3", "#6b6ecf", "#333333", "#cedb9c", "#8ca252",
  "#a55194", "#e5e56f", "#66a61e", "#e6ab02", "#a6761d", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ffff33",
  "#f781bf", "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3", "#1e90ff",
  "#ff4500", "#32cd32", "#ff0000", "#8a2be2", "#a0522d", "#ff66cc", "#9c9ede", "#adff2f", "#00ced1", "#ffd700",
  "#6699cc", "#cc6644", "#66aa66", "#cc6666", "#9966cc", "#996633", "#cc99cc", "#99cc44", "#66cccc", "#cccc66",
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
  "#ffffcc", "#e7ba52", "#ce6dbd", "#d6616b", "#b5cf6b", "#dbdb5c", "#e7cb94", "#ad494a", "#bd9e39", "#de9ed6",
  "#e7969c", "#cedb9c", "#33a02c", "#b2df8a", "#fdbf6f", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7",
  "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5",
  "#ffed6f", "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#999999"
)

tab_palettes <- c(
  # tab20
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#ffff33",
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#ffffcc",
  # tab20 (saturated)
  "#1e90ff", "#ff4500", "#32cd32", "#ff0000", "#8a2be2", "#a0522d", "#ff33ff", "#000000", "#adff2f", "#00ced1", "#ffd700",
  "#63b3ff", "#ff8052", "#66e066", "#ff6666", "#b199e8", "#cd8866", "#ff99ff", "#999999", "#d4ff7f", "#66f0f0", "#ffee66",
  # tab20 (muted)
  "#6699cc", "#cc6644", "#66aa66", "#cc6666", "#9966cc", "#996633", "#cc66cc", "#666666", "#99cc44", "#66cccc","#cccc66",
  "#99bbee", "#ddbbaa", "#aacc99", "#ee9999", "#bbaadd", "#ccaa88", "#ee99ee", "#bbbbbb", "#ddffaa", "#66dddd","#ffffaa",
  # tab20c (modified)
  "#393b79", "#e7ba52", "#637939", "#843c39", "#6b6ecf", "#8c6d31", "#ce6dbd", "#d6616b", "#b5cf6b", "#7b4173", "#dbdb5c",
  "#5254a3", "#e7cb94", "#8ca252", "#ad494a", "#9c9ede", "#bd9e39", "#de9ed6", "#e7969c", "#cedb9c", "#a55194", "#e5e56f"
)

# ---- Logging Related Functions ----

# Suppress warnings and messages
quiet_msg <- function(expr) {
  # create a temp file and open a connection
  tmp <- tempfile()
  con <- file(tmp, open = "wt")
  
  # sink both output and message streams to the same connection
  sink(con, type = "output")
  sink(con, type = "message")
  
  result <- NULL
  tryCatch({
    result <- expr
  }, finally = {
    # restore normal output in the correct order
    sink(type = "message")
    sink(type = "output")
    close(con)
  })
  
  return(result)
}

# Log info messages (green)
log_info <- function(sample, step, message) {
  prefix <- green(formatC("[INFO]", width = 7, flag = " "))
  message(glue::glue("{prefix} [{sample} | {toupper(step)}] {message}"))
}

# Log warning messages (yellow)
log_warn <- function(sample, step, message) {
  prefix <- yellow(formatC("[WARN]", width = 7, flag = " "))
  message(glue::glue("{prefix} [{sample} | {toupper(step)}] {message}"))
}

# Log error messages (red)
log_error <- function(sample, step, message) {
  prefix <- red(formatC("[ERROR]", width = 7, flag = " "))
  message(glue::glue("{prefix} [{sample} | {toupper(step)}] {message}"))
  stop("Workflow Stopped.", call. = FALSE)
}

# Optional: header for sample processing
log_sample_header <- function(sample) {
  cat(blue$bold(glue::glue("\n--- Processing Sample: {sample} ---\n\n")))
  
}

# Optional: section header
log_section <- function(section_name) {
  cat(magenta$bold(glue::glue("\n[{toupper(section_name)}]\n")))
}

# ---- BULK RNA SEQ ANALYSIS RELATED FUNCTIONS ----

# # Copy and paste in wrapper for "User Override Project Directories & Parameters"
# proj <- ""
# species <- ""
# contrasts <- c("Treatment1-Control1",
#                "Treatment2-Control2")
# 
# parent_dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data"
# gmt_dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets"
# 
# # DESeq2 overrides
# deseq2.override <- list(
#   contrasts     = contrasts
#   #design        = "Comparisons",            # DESeq2 design formula or column name
#   #lfc.cutoff    = 0,                        # Log fold change cutoff for significance
#   #padj.cutoff   = 0.1,                      # Adjusted p-value cutoff for significance
#   #batch.correct = FALSE                     # Boolean, whether to apply batch correction
# )
# 
# # Heatmap overrides
# heatmap.override <- list(
#   #force.log        = TRUE,                  # Force log transformation
#   col.ann          = NULL,                  # Column annotation
#   #row.ann          = NULL,                  # Row annotation
#   col.gaps         = NULL,                  # Column gaps
#   #row.gaps         = NULL,                  # Row gaps
#   col.cluster      = "all",                 # Column clustering
#   #row.cluster      = "all",                 # Row clustering
#   #palette         = "rdbu",                # Heatmap palette
#   #ann.palette     = "discrete",            # Annotation palette
#   #border.color    = NA,                    # Cell border color
#   #show.expr.legend = TRUE,                  # Show expression legend
#   #title           = "",                    # Heatmap title
#   #format           = "tiff"                 # Output file format
# )
# 
# # Volcano plot overrides
# volcano.override <- list(
#   #lfc.cutoff  = 0.58,                         # Minimum log2 fold-change to highlight
#   #padj.cutoff = 0.05,                      # Adjusted p-value cutoff
#   #color       = "vrds",                    # Color palette
#   label.genes = c()                         # Genes to label on the plot
# )
# # Setup project
# proj.params <- setup_project(
#   proj = proj,
#   species = species,
#   contrasts = contrasts,
#   parent_dir = parent_dir,
#   gmt_dir = gmt_dir,
#   deseq2.override = deseq2.override,
#   heatmap.override = heatmap.override,
#   volcano.override = volcano.override
# )

# Default Project Directories & Parameters
setup_project <- function(proj, species, contrasts,
                          parent_dir, gmt_dir, scripts_dir = NULL,
                          deseq2.override  = list(), heatmap.override = list(),
                          volcano.override = list(), pathway.override = list()) {
  
  # ---- 🛠️ Global Environment Configuration ----
  
  options(future.globals.maxSize = 1e15)            # Increase future globals size for parallelization
  options(Seurat.object.assay.version = "v5")       # Ensure Seurat uses v5 assay format
  set.seed(1234)                                    # Set seed for reproducibility
  
  # ---- 📊 Default Parameter Definition ----
  
  # Default DESeq2 Parameters
  default.deseq2 <- list(
    contrasts     = c("Treatment-Control"),   # Vector of contrasts for DE analysis
    design        = "Comparisons",            # DESeq2 design formula or column name
    lfc.cutoff    = 0,                        # Log fold change cutoff for significance
    padj.cutoff   = 0.1,                      # Adjusted p-value cutoff for significance
    batch.correct = FALSE                     # Boolean, whether to apply batch correction
  )
  
  # Default Heatmap Parameters
  default.heatmap <- list(
    force.log        = FALSE,       # Force log transform (default FALSE i.e. auto detect)
    col.ann          = NULL,        # Single/multiple columns from metadata_col for column annotation
    row.ann          = NULL,        # Single/multiple columns from metadata_row for row annotation
    col.gaps         = NULL,        # Single column from metadata_col to define column gaps in heatmap
    row.gaps         = NULL,        # Single column from metadata_row to define row gaps in heatmap
    col.cluster      = "all",       # Single column from metadata_col, "all", "alphabetical" for clustering columns
    row.cluster      = "all",       # Single column from metadata_row, "all", "alphabetical" for clustering columns
    palette          = "rdbu",      # Color palette for heatmap matrix ("rdbu" or "vrds")
    ann.palette      = "discrete",  # Color palette for heatmap annotation ("discrete" or "sequential")
    border.color     = NA,          # Color of heatmap cell borders (default NA i.e. no border)
    show.expr.legend = TRUE,        # Show expression legend (set FALSE if annotations overlap)
    title            = NA           # Title for heatmap (default NA i.e. no title)
  )
  
  # Default Volcano Parameters
  default.volcano <- list(
    lfc.cutoff   = 0.58,            # Log fold change cutoff for volcano plot
    padj.cutoff  = 0.05,            # Adjusted p-value cutoff for volcano plot
    color        = "vrds",          # Color palette for volcano plot ("vrds", etc.)
    label.genes  = NULL             # Optional vector of genes to label on volcano plot
  )
  
  # Apply Overrides
  deseq2.params  <- modifyList(default.deseq2,  deseq2.override)
  heatmap.params <- modifyList(default.heatmap, heatmap.override)
  volcano.params <- modifyList(default.volcano, volcano.override)
  
  # ---- 📂 Directory Structure Setup ----
  
  proj_dir      <- file.path(parent_dir, proj)
  
  # Bulk RNA-seq Directories
  counts_dir    <- file.path(proj_dir, "counts")               # Directory containing count files
  contrast_dir  <- file.path(proj_dir, contrasts)              # Directory to store results for each contrast
  deseq2_dir    <- file.path(contrast_dir, "DEG_Analysis")     # Directory to store DESeq2 results
  pathway_dir   <- file.path(contrast_dir, "Pathway_Analysis") # Directory to store Pathway analysis results
  tf_dir        <- file.path(contrast_dir, "TF_Analysis")      # Directory to store TF analysis results
  
  # Single-Cell RNA-seq Directories
  diagnostics_dir        <- file.path(proj_dir, "04.Diagnostics")
  demux_dir              <- file.path(proj_dir, "05.Demux")
  raw_matrix_dir         <- file.path(proj_dir, "06.Matrix", "raw_feature_bc_matrix")
  filt_matrix_dir        <- file.path(proj_dir, "06.Matrix", "filt_feature_bc_matrix")
  hto_matrix_dir         <- file.path(proj_dir, "06.Matrix", "HTO_bc_matrix")
  seurat_dir             <- file.path(proj_dir, "Seurat")
  sc_deseq2_dir          <- file.path(proj_dir, "DESeq2")
  pyscenic_dir           <- file.path(proj_dir, "pySCENIC")
  scvelo_dir             <- file.path(proj_dir, "scVelo")
  velocyto_dir           <- file.path(proj_dir, "velocyto")
  cellphonedb_dir        <- file.path(proj_dir, "CellphoneDB")
  cellchat_dir           <- file.path(proj_dir, "CellChat")
  
  # ---- 🧬 Cell Marker and Gene Set Extraction ----
 
  markers <- list()
  s_genes <- c()
  g2m_genes <- c()
  
  if (!is.null(scripts_dir)){ 
    metafile               <- file.path(scripts_dir, "scRNASeq_Metadata", paste0(proj, "_Metadata.xlsx"))
    cell.cycle.marker.file <- file.path(scripts_dir, "Cell_Cycle_Markers.xlsx")
    cell.type.marker.file  <- file.path(scripts_dir, "Cell_Type_Markers.xlsx")
    
    # Extract Cell type Genes
    if (file.exists(cell.type.marker.file)) {
      markers <- openxlsx::read.xlsx(cell.type.marker.file)
    }
    
    # Extract Cell Cycle Genes (Human and Mouse)
    if (file.exists(cell.cycle.marker.file)) {
      cc_markers <- openxlsx::read.xlsx(cell.cycle.marker.file)
      
      # Combine Human and Mouse genes for each phase
      s_genes <- c(cc_markers$Human_Gene[cc_markers$Phase == "S"],
                   cc_markers$Mouse_Gene[cc_markers$Phase == "S"])
      
      g2m_genes <- c(cc_markers$Human_Gene[cc_markers$Phase == "G2/M"],
                     cc_markers$Mouse_Gene[cc_markers$Phase == "G2/M"])
      
      # Clean up potential NA/empty values resulting from unlist/subsetting
      s_genes <- s_genes[!is.na(s_genes) & s_genes != ""]
      g2m_genes <- g2m_genes[!is.na(g2m_genes) & g2m_genes != ""]
    }
  }
  
  # ---- 📦 Final Project Parameters Construction ----
  
  proj.params <- c(
    list(
      # Project info
      proj        = proj,
      species     = species,
      
      # Bulk RNA-seq directories
      proj_dir    = normalizePath(proj_dir,     mustWork = FALSE),
      counts_dir  = normalizePath(counts_dir,   mustWork = FALSE),
      gmt_dir     = normalizePath(gmt_dir,      mustWork = FALSE),
      contrast_dir= normalizePath(contrast_dir, mustWork = FALSE),
      deseq2_dir  = normalizePath(deseq2_dir,   mustWork = FALSE),
      pathway_dir = normalizePath(pathway_dir,  mustWork = FALSE),
      tf_dir      = normalizePath(tf_dir,       mustWork = FALSE),
      
      # Bulk RNA-seq parameters
      contrasts   = contrasts,
      deseq2      = deseq2.params,
      heatmap     = heatmap.params,
      volcano     = volcano.params,
      
      # Single cell RNA-seq directories
      diagnostics_dir  = normalizePath(diagnostics_dir, mustWork = FALSE),
      demux_dir        = normalizePath(demux_dir,       mustWork = FALSE),
      raw_matrix_dir   = normalizePath(raw_matrix_dir,  mustWork = FALSE), 
      filt_matrix_dir  = normalizePath(filt_matrix_dir, mustWork = FALSE), 
      hto_matrix_dir   = normalizePath(hto_matrix_dir,  mustWork = FALSE), 
      seurat_dir       = normalizePath(seurat_dir,      mustWork = FALSE), 
      sc_deseq2_dir    = normalizePath(sc_deseq2_dir,   mustWork = FALSE), 
      pyscenic_dir     = normalizePath(pyscenic_dir,    mustWork = FALSE), 
      scvelo_dir       = normalizePath(scvelo_dir,      mustWork = FALSE), 
      velocyto_dir     = normalizePath(velocyto_dir,    mustWork = FALSE), 
      cellphonedb_dir  = normalizePath(cellphonedb_dir, mustWork = FALSE), 
      cellchat_dir     = normalizePath(cellchat_dir,    mustWork = FALSE),
      
      # Single cell data
      metafile = normalizePath(metafile,       mustWork = FALSE),
      markers = markers,
      cell_cycle = list(S = s_genes, G2M = g2m_genes)
    )
  )
  
  # ---- ⚠️ Final Checks and Return ----
  
  if (proj == "") warning("⚠️ Project name is empty")
  if (species == "") warning("⚠️ Species is not set")
  if (length(contrasts) == 0) warning("⚠️ No contrasts specified")
  
  return(proj.params)
}


# Helper functions #1
save_xlsx <- function(data, file, sheet_name, row_names) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = sheet_name)
  openxlsx::writeData(wb, sheet = sheet_name, x = data, rowNames = row_names)
  if(row_names){
    openxlsx::writeData(wb, sheet = sheet_name, x = "SYMBOL", startCol = 1, startRow = 1)
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

# Helper functions #2
filter_samples_by_contrast <- function(metadata, contrast) {
  contrast <- base::gsub(pattern = "\\(|\\)", replacement = "", x = contrast)
  metadata %>%
    dplyr::filter(Comparisons %in% stringr::str_split(contrast, "-")[[1]]) %>%
    dplyr::pull(Sample_ID) %>%
    as.character()
}

merge_counts <- function(proj.params) {
  
  # Check required proj.params attributes ----
  required_attrs <- c("counts_dir", "proj_dir", "proj")
  for (attr in required_attrs) {
    if (is.null(proj.params[[attr]])) {
      stop("⚠️ Missing required proj.params attribute: ", attr)
    } else {
      message("✔ Found proj.params attribute: ", attr)
    }
  }
  
  set.seed(1234)
  
  # Get count files
  count_files <- list.files(path = proj.params$counts_dir, pattern = "\\.txt$|ReadsPerGene\\.out\\.tab$", full.names = TRUE)
  if (length( count_files) == 0) {
    stop("No count files found in the directory.")
  }
  
  # ---- Initialize ----
  all_counts <- list()
  gene_lists <- list()
  sample_ids <- character()
  
  # ---- Define special counters generated by HTSeq and STAR outputs ----
  special_counters <- c("__no_feature", "__ambiguous", "__too_low_aQual", 
                        "__not_aligned", "__alignment_not_unique", "__assignment_counts",
                        "N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
  
  # ---- Parse Files ----
  for (count_file in count_files) {
    
    # Read count file
    df <- tryCatch({
      read.table(file = count_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    }, error = function(e) {
      stop("Error reading file: ", count_file, " — ", e$message)
    })
    
    if (ncol(df) < 4) stop("count file does not have expected 4 columns: ", count_files[i])
    
    # Remove special counters
    df <- df %>% dplyr::filter(!(.[[1]] %in% special_counters))
    
    sample_id <- gsub("\\..*$|ReadsPerGene\\.out\\.tab", "", basename(count_file))
    gene_ids <- df[[1]]
    strand_sums <- colSums(df[2:4], na.rm = TRUE)
    
    # Determine strandedness
    if (abs((strand_sums[1]/strand_sums[2]) - (strand_sums[1]/strand_sums[3])) < 2) {
      message("Detected unstranded library for: ", sample_id)
      counts <- df[[2]]
    } else if (strand_sums[2] > 3 * strand_sums[3]) {
      message("Detected positively stranded library for: ", sample_id)
      counts <- df[[3]]
    } else if (strand_sums[3] > 3 * strand_sums[2]) {
      message("Detected negatively stranded library for: ", sample_id)
      counts <- df[[4]]
    } else {
      stop("Could not determine strandedness for: ", f)
    }
    
    all_counts[[sample_id]] <- counts
    gene_lists[[sample_id]] <- gene_ids
    sample_ids <- c(sample_ids, sample_id)
  }
  
  # ---- Check Gene Consistency ----
  ref_genes <- gene_lists[[1]]
  for (i in seq_along(gene_lists)) {
    if (!identical(ref_genes, gene_lists[[i]])) {
      stop("Gene mismatch detected in file: ", names(gene_lists)[i])
    }
  }
  
  # ---- Build Count Matrix ----
  count_matrix <- do.call(cbind, all_counts)
  colnames(count_matrix) <- sample_ids
  count_matrix <- data.frame(SYMBOL = ref_genes, count_matrix, stringsAsFactors = FALSE)
  
  # ---- Filter Rows and Columns with All Zeros ----
  count_matrix <- count_matrix[rowSums(count_matrix[,-1]) > 0, , drop = FALSE]
  count_matrix <- count_matrix[, c(TRUE, colSums(count_matrix[,-1]) > 0), drop = FALSE]
  
  # ---- Export ----
  filename <- file.path(proj.params$proj_dir, paste0(proj.params$proj, "_Raw_counts.xlsx"))
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "Raw_counts")
  openxlsx::writeData(wb = wb, sheet = "Raw_counts", x = count_matrix)
  openxlsx::saveWorkbook(wb = wb, file = filename, overwrite = TRUE)
  message("✅ Saved counts to: ", file.path(proj.params$proj_dir, proj.params$proj))
  
  return(count_matrix)
}

prepare_deseq2_input <- function(meta_data, read_data, project_params) {
  
  set.seed(1234)
  
  # ---- Input checks ----
  if (!"Sample_ID" %in% colnames(meta_data)) {
    stop("`meta_data` must contain a 'Sample_ID' column.")
  }
  
  if (!"SYMBOL" %in% colnames(read_data)) {
    stop("`read_data` must contain a 'SYMBOL' column.")
  }
  if (anyDuplicated(read_data$SYMBOL)) {
    stop("The 'SYMBOL' column in `read_data` must not contain duplicate values.")
  }
  
  if (is.null(project_params$deseq2$design)) {
    stop("project_params$deseq2$design is missing. Please define the DESeq2 design formula.")
  }
  
  # ---- Initial summary ----
  message("Samples before filtering: ", nrow(meta_data))
  message("Genes before filtering: ", nrow(read_data))
  
  # ---- Clean and align meta_data ----
  meta_data <- meta_data %>%
    dplyr::filter(!is.na(Sample_ID)) %>%
    dplyr::mutate(Sample_ID = make.names(Sample_ID, unique = TRUE)) %>%
    dplyr::filter(Sample_ID %in% make.names(colnames(read_data)))
  rownames(meta_data) <- make.names(meta_data$Sample_ID)
  
  # Add Batch column if missing
  if (!"Batch" %in% colnames(meta_data)) {
    meta_data$Batch <- 1
    warning("No 'Batch' column found. Assigning all samples to Batch 1.")
  }
  
  # ---- Clean and align read_data ----
  colnames(read_data) <- make.names(colnames(read_data))
  valid_samples <- intersect(colnames(read_data), rownames(meta_data))
  
  read_data <- read_data %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "SYMBOL") %>%
    dplyr::select(all_of(valid_samples)) %>%
    replace(is.na(.), 0)
  
  read_data <- read_data[rowSums(read_data) != 0, ]
  
  # ---- Remove zero-count samples ----
  zero_samples <- which(colSums(read_data) == 0)
  if (length(zero_samples) > 0) {
    read_data <- read_data[, -zero_samples, drop = FALSE]
    meta_data <- meta_data[-zero_samples, , drop = FALSE]
  }
  
  # ---- Remove samples with NA in design variables ----
  design_vars <- stringr::str_split(string = project_params$deseq2$design, pattern = "[+*:]")[[1]]
  na_samples <- unique(unlist(lapply(design_vars, function(var) {
    if (var %in% colnames(meta_data)) {
      which(is.na(meta_data[[var]]))
    } else {
      warning(glue::glue("Variable '{var}' not found in meta_data."))
      integer(0)
    }
  })))
  
  if (length(na_samples) > 0) {
    read_data <- read_data[, -na_samples, drop = FALSE]
    meta_data <- meta_data[-na_samples, , drop = FALSE]
  }
  
  # ---- Remove sizeFactor column if present ----
  meta_data <- meta_data[, colnames(meta_data) != "sizeFactor", drop = FALSE]
  
  # ---- Convert meta_data columns to factors ----
  message("Structure of meta_data before conversion:")
  str(meta_data)
  meta_data[] <- lapply(meta_data, as.factor)
  message("Structure of meta_data after conversion:")
  str(meta_data)
  
  # ---- Reorder read_data to match meta_data ----
  read_data <- read_data[, rownames(meta_data), drop = FALSE]
  
  # ---- Sanity checks ----
  if (!is.data.frame(read_data) || !is.data.frame(meta_data)) {
    stop("`read_data` and `meta_data` must be data.frames")
  }
  if (!all(colnames(read_data) %in% rownames(meta_data))) {
    stop("Some samples in `read_data` are missing in `meta_data`")
  }
  if (!all(colnames(read_data) == rownames(meta_data))) {
    stop("Sample order mismatch between `read_data` and `meta_data`")
  }
  
  # ---- Return cleaned data ----
  return(invisible(list(meta_data = meta_data, read_data = read_data)))
}

plot_pca <- function(meta_data, read_data, output_dir, top_n_genes = 500, file_name = "PCA_Plots.pdf") {
  
  set.seed(1234)
  
  # ---- Input Checks ----
  if (!is.data.frame(meta_data)) {
    stop("`meta_data` must be a data.frame")
  }
  if (!is.data.frame(read_data)) {
    stop("`read_data` must be a data.frame")
  }
  if (!"Sample_ID" %in% colnames(meta_data)) {
    stop("`meta_data` must contain a 'Sample_ID' column.")
  }
  
  # Warn if output path does not exist
  if (!dir.exists(output_dir)) {
    warning("Output path does not exist. Attempting to create: ", output_dir)
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ---- DESeq2 Object Preparation ----
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = read_data,
    colData = meta_data,
    design = ~1
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  
  # ---- VST Transformation ----
  vsd <- DESeq2::vst(dds, blind = TRUE)
  
  # ---- PCA Data Preparation ----
  vst_mat <- SummarizedExperiment::assay(vsd) %>%
    as.data.frame() %>%
    dplyr::mutate(row_variance = matrixStats::rowVars(as.matrix(.))) %>%
    dplyr::slice_max(order_by = row_variance, n = top_n_genes) %>%
    dplyr::select(-row_variance)
  
  # ---- PCA Calculation ----
  pca <- stats::prcomp(x = t(vst_mat), center = TRUE, scale. = FALSE)
  
  # ---- Merge PCA Output with Metadata ----
  pca_df <- pca$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample_ID")
  
  df <- dplyr::inner_join(meta_data, pca_df, by=c("Sample_ID"="Sample_ID"))
  
  # ---- Variance explained ----
  percentVar <- round(100 * summary(pca)$importance[2, 1:2])
  comp_variables <- setdiff(colnames(meta_data), "Sample_ID")
  
  # ---- Save all plots to a single PDF ----
  output_file <- file.path(output_dir, file_name)
  message("Saving PCA plots to: ", output_file)
  
  grDevices::pdf(file = output_file, width = 8.5, height = 11)
  
  for (var in comp_variables) {
    
    n_unique <- length(unique(meta_data[[var]]))
    if (n_unique < 2 || n_unique == nrow(meta_data)) {
      warning("Variable ", var, " is either constant or has one unique value 
      per sample; skipping PCA plot.")
      next
    }
    
    # Define color palette
    pca_palette <- custom_palette[1:length(unique(meta_data[[var]]))]
    names(pca_palette) <- as.character(unique(meta_data[[var]]))
    
    # Create PCA plot
    p <- ggplot2::ggplot(data = df, mapping = aes(x = PC1, y = PC2, color = get(var))) +
      ggplot2::geom_point(size = 3, shape = 16) +
      ggrepel::geom_text_repel(ggplot2::aes(label = Sample_ID), show.legend = FALSE) +
      ggplot2::theme_light() +
      ggplot2::labs(
        color = var,
        x = paste0("PC1: ", percentVar[1], "% variance"),
        y = paste0("PC2: ", percentVar[2], "% variance")
      ) +
      custom_theme +
      ggplot2::scale_color_manual(values = pca_palette)
    
    print(p)  # print to PDF (one page per plot)
  }
  
  dev.off()
}

plot_ma <- function(dds, output_dir, file_name = "MA_Plot.pdf") {
  
  set.seed(1234)
  
  # ---- Input Checks ----
  if (!inherits(dds, "DESeqDataSet")) {
    stop("`dds` must be a DESeqDataSet object.")
  }
  
  if (!dir.exists(output_dir)) {
    warning("Output path does not exist. Creating: ", output_dir)
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ---- MA Plot ----
  output_file <- file.path(output_dir, file_name)
  message("Saving MA plot to: ", output_file)
  
  grDevices::pdf(file = output_file, width = 8.5, height = 11)
  
  DESeq2::plotMA(
    object = dds,
    alpha  = 0.1, # FDR threshold (blue dots for significant genes)
    main   = "MA Plot",
    xlab   = "Mean of Normalized Counts",
    MLE    = FALSE
  )
  
  grDevices::dev.off()
  message("MA plot completed successfully.")
}

plot_dispersion <- function(dds, output_dir, file_name = "Dispersion_Plot.pdf") {
  
  set.seed(1234)
  
  # ---- Input Checks ----
  if (!inherits(dds, "DESeqDataSet")) {
    stop("`dds` must be a DESeqDataSet object.")
  }
  
  if (!dir.exists(output_dir)) {
    warning("Output path does not exist. Creating: ", output_dir)
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ---- Dispersion Plot ----
  output_file <- file.path(output_dir, file_name)
  message("Saving dispersion plot to: ", output_file)
  
  grDevices::pdf(file = output_file, width = 8.5, height = 11)
  
  DESeq2::plotDispEsts(
    object  = dds,
    genecol = "black",
    fitcol  = "red",
    finalcol= "dodgerblue",
    legend  = TRUE,
    xlab    = "Mean of Normalized Counts",
    ylab    = "Dispersion",
    log     = "xy",
    cex     = 0.45
  )
  
  grDevices::dev.off()
  
  message("Dispersion plot completed successfully.")
  # Expected results: Higher the mean, lower the dispersion
}

plot_volcano <- function(DEGs_df, proj.params, contrast = "Target-Reference", 
                         output_dir, top_n = 5, file_prefix = "Volcano_Plot") {
  
  set.seed(1234)
  
  # ---- Input Checks ----
  required_cols <- c("log2FoldChange", "padj", "SYMBOL")
  if (!all(required_cols %in% colnames(DEGs_df))) {
    stop("DEGs_df must contain columns: log2FoldChange, padj, SYMBOL")
  }
  
  if (!dir.exists(output_dir)) {
    warning("Output path does not exist. Creating: ", output_dir)
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ---- Check required proj.params attributes ----
  required_attrs <- c("lfc.cutoff", "padj.cutoff")
  if (is.null(proj.params$volcano)) {
    stop("⚠️ proj.params must contain a 'volcano' list")
  }
  
  for (attr in required_attrs) {
    if (is.null(proj.params$volcano[[attr]])) {
      stop("⚠️ Missing required proj.params$volcano attribute: ", attr)
    } else {
      val <- proj.params$volcano[[attr]]
      message("✔ Found proj.params$volcano$", attr, " = ", paste(val, collapse = ", "))
    }
  }
  
  
  # ---- Volcano Parameters ----
  padj_cutoff <- proj.params$volcano$padj.cutoff
  lfc_cutoff  <- proj.params$volcano$lfc.cutoff
  label_genes <- if(!is.null(proj.params$volcano$label.genes)) proj.params$volcano$label.genes else NULL
  
  target <- stringr::str_split(string = contrast, pattern = "-")[[1]][1]
  reference <- stringr::str_split(string = contrast, pattern = "-")[[1]][2]
  
  # ---- Format DEGs ----
  DEGs_df <- DEGs_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(
      Direction = dplyr::case_when(
        padj < padj_cutoff & log2FoldChange > lfc_cutoff  ~ paste0("Up in ", target),
        padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ paste0("Up in ", reference),
        TRUE ~ "Not Significant"
      ),
      padj = dplyr::case_when(padj == 0 ~ min(padj[padj > 0], na.rm = TRUE), 
                              TRUE ~ padj),
      Significance = dplyr::case_when(
        abs(log2FoldChange) >= lfc_cutoff & padj <= 0.001 ~ "FDR < 0.001",
        abs(log2FoldChange) >= lfc_cutoff & padj <= 0.01  ~ "FDR < 0.01",
        abs(log2FoldChange) >= lfc_cutoff & padj <= 0.05  ~ "FDR < 0.05",
        TRUE ~ "Not Significant"
      ),
      Relevance = abs(log2FoldChange) * -log10(padj)
    )
  
  # ---- Color Palettes ----
  volcano_palette <- c(
    viridis::viridis(100)[50],  # Up in reference
    viridis::viridis(100)[1],   # Up in target
    viridis::viridis(100)[100]  # Not Significant
  )
  names(volcano_palette) <- c(
    paste0("Up in ", reference),
    paste0("Up in ", target),
    "Not Significant"
  )
  
  alpha_palette <- c("FDR < 0.001" = 1, "FDR < 0.01" = 0.8, 
                     "FDR < 0.05" = 0.6, "Not Significant" = 0.4)
  
  # ---- Axis limits and breaks ----
  x_vals <- DEGs_df$log2FoldChange
  y_vals <- -log10(DEGs_df$padj)
  
  # Keep only finite values
  x_vals <- x_vals[is.finite(x_vals)]
  y_vals <- y_vals[is.finite(y_vals)]
  
  # Round to nearest integer for nice axis limits
  x_min <- floor(min(x_vals, na.rm = TRUE))
  x_max <- ceiling(max(x_vals, na.rm = TRUE))
  y_min <- 0  # Start y-axis at 0 for volcano
  y_max <- ceiling(max(y_vals, na.rm = TRUE))
  
  # Set x-axis breaks in reasonable bins (approx 5 units each)
  x_bin <- max(abs(floor(x_min / 5)), abs(ceiling(x_max / 5)))
  x_breaks <- seq(from = -max(x_max, abs(x_min)), to = max(x_max, abs(x_min)), by = x_bin)
  x_breaks <- x_breaks[!x_breaks <= x_min-x_bin]
  x_breaks <- x_breaks[!x_breaks >= x_max+x_bin] 
  
  # Set y-axis breaks (dynamic, based on magnitude)
  y_bin <- if (y_max > 100) 100 else if (y_max > 10) 10 else 1
  y_breaks <- seq(from = y_min, to = ceiling(y_max/y_bin)*y_bin, by = y_bin)
  
  # ---- Build Base Plot ----
  p <- ggplot2::ggplot(data = DEGs_df, 
                       mapping = aes(x = log2FoldChange, y = -log10(padj),
                                     color = Direction, alpha = Significance,
                                     size = Relevance)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.05, height = 0.05)) +
    ggplot2::theme_classic() +
    custom_theme +
    ggplot2::labs(x = expression("log"[2]*"FC"),
                  y = expression("-log"[10]*"padj"),
                  fill = "Direction",
                  title = contrast) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), 
                        linetype = "dotted", color = "black", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), 
                        linetype = "dotted", color = "black", linewidth = 0.5) +
    ggplot2::scale_color_manual(values = volcano_palette) +
    ggplot2::scale_alpha_manual(values = alpha_palette) +
    ggplot2::scale_size_continuous(range = c(0, 3)) +
    ggplot2::coord_cartesian(xlim = c(min(x_breaks), max(x_breaks)),
                             ylim = c(min(y_breaks), max(y_breaks))) +
    ggplot2::scale_x_continuous(breaks = x_breaks, 
                                labels = function(x) { base::ifelse(x %% 1 == 0, as.integer(x), format(x, digits = 2)) }) +
    ggplot2::scale_y_continuous(breaks = y_breaks) +
    ggplot2::guides(size = "none",
                    shape = guide_legend(override.aes = list(size = 3)),
                    fill = guide_colourbar(theme = theme(legend.key.width = unit(0.75, "lines"),
                                                         legend.key.height = unit(10, "lines"),
                                                         legend.ticks = element_blank(),
                                                         legend.frame = element_rect(colour = "Black", linewidth = 1)))) 
  
  # ---- Identify Top 5 Up and Down regulated Genes ----
  predicted_gene_pattern <- "^(Gm[0-9]+|ENSMUSG[0-9]+|ENSG[0-9]+|LOC[0-9]+|C[0-9]+orf[0-9]+|RP[0-9]+-)|Rik$"
  
  # Filter out predicted/placeholder genes and non-significant genes
  filtered_DEGs <- DEGs_df %>%
    dplyr::filter(!stringr::str_detect(string = SYMBOL, pattern = predicted_gene_pattern)) %>%
    dplyr::filter(padj < padj_cutoff)
  
  # Top 5 up-regulated (log2FC > lfc_cutoff)
  top_up <- filtered_DEGs %>%
    dplyr::filter(log2FoldChange > lfc_cutoff) %>%
    #dplyr::slice_min(padj, n = 5, with_ties = FALSE) %>%
    dplyr::slice_max(Relevance, n = 5, with_ties = FALSE) %>%
    dplyr::pull(SYMBOL)
  
  # Top 5 down-regulated (log2FC < -lfc_cutoff)
  top_down <- filtered_DEGs %>%
    dplyr::filter(log2FoldChange < -lfc_cutoff) %>%
    #dplyr::slice_min(padj, n = 5, with_ties = FALSE) %>%
    dplyr::slice_max(Relevance, n = 5, with_ties = FALSE) %>%
    dplyr::pull(SYMBOL)
  
  # Combine
  top_genes <- unique(c(top_up, top_down))
  
  # ---- Label Top Genes ----
  if (!is.null(label_genes)) {
    genes_to_label <- intersect(label_genes, DEGs_df$SYMBOL)
  } else {
    genes_to_label <- top_genes
  }
  
  q <- p + ggrepel::geom_text_repel(
    data = DEGs_df %>% dplyr::filter(SYMBOL %in% genes_to_label),
    aes(label = SYMBOL),
    direction = "both",
    box.padding = 0.8,                # ↓ smaller padding around label
    point.padding = 0.1,              # minimal space between point and line start
    max.overlaps = nrow(DEGs_df),
    show.legend = FALSE,
    min.segment.length = 0,           # Only draw segments longer than this
    segment.curvature = -0.5,         # Negative = curve upward, positive = downward
    segment.ncp = 50,                 # More control points = smoother curves
    segment.angle = 20,               # Affects entry/exit angles
    segment.size = 0.5,               # Optional: line thickness
    size = 4,                         # text size in mm (1 mm = 2.83 points)
    position = ggbeeswarm::position_quasirandom(width = 0.1, varwidth = TRUE)
  )
  
  # ---- Save Plots ----
  ggplot2::ggsave(file.path(output_dir, paste0(file_prefix, "_", contrast, ".pdf")),
                  plot = p, width = 7, height = 7, device = "pdf")
  
  ggplot2::ggsave(file.path(output_dir, paste0(file_prefix, "_top_", contrast, ".pdf")),
                  plot = q, width = 7, height = 7, device = "pdf")
  
  message("Volcano plots saved successfully.")
  
  return(invisible(p))
}

# metadata has column Sample_ID and columns defined in proj.params$heatmap.col.ann
# metadata_row has column SYMBOL and columns defined in proj.params$heatmap.row.ann
# metadata_row <- NULL if no row annotations needed
# disp_genes either empty vector c() or vector of genes
# norm_counts is a matrix with gene symbols as rownames
plot_heatmap <- function(norm_counts, proj.params, metadata_col = NULL, metadata_row = NULL, disp_genes = c()) {
  
  set.seed(1234)
  
  # ---- Check required proj.params attributes ----
  if (is.null(proj.params$heatmap)) {
    stop("proj.params must contain a 'heatmap' list")
  }
  required_attrs <- c("force.log","col.ann","ann.palette","palette",
                      "col.cluster","row.cluster",
                      "title","border.color","show.expr.legend")
  
  for (attr in required_attrs){
    if(is.null(proj.params$heatmap[[attr]])){
      stop("Missing proj.params$heatmap$", attr)
    } else {
      message("✔ Found proj.params$heatmap$", attr, " = ", paste(proj.params$heatmap[[attr]], collapse=", "))
    }
  }
  
  # ---- Input check ----
  if (nrow(norm_counts) < 2){
    message("Input data frame has less than 2 genes. Skipping plotting.")
    return(NULL)
  }
  
  # ---- Prepare matrix ----
  mat <- norm_counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("SYMBOL") %>%
    base::replace(is.na(.), 0) %>%
    # This section retains gene copy with highest expression in case of duplicates
    # dplyr::mutate(n = rowSums(.[, -1])) %>%
    # dplyr::group_by(SYMBOL) %>%
    # dplyr::slice_max(n) %>%
    # dplyr::ungroup() %>%
    # dplyr::filter(n != 0) %>%
    # dplyr::select(-n) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("SYMBOL") %>%
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
  
  # Remove genes that have 0 across all samples
  mat <- mat[rowSums(mat) != 0, ]
  
  rownames(mat) <- make.names(rownames(mat), unique = TRUE)
  colnames(mat) <- make.names(colnames(mat), unique = TRUE)
  
  # Log transform if necessary
  quantiles <- stats::quantile(x = as.vector(as.matrix(mat)), probs = c(0, 0.01, 0.99, 1), na.rm = TRUE)
  huge_range <- (quantiles[4] - quantiles[1]) > 100   # Range of values greater than 100
  only_pos <- quantiles[1] >= 0                       # Min value greater than 0
  if ((huge_range & only_pos) | proj.params$heatmap$force.log){
    mat <- log2(1 + mat)
  }
  
  # Scale every feature across samples 
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[is.na(mat_scaled)] <- 0
  
  # ---- Annotations ----
  # Define column annotations 
  if(is.null(metadata_col)) {
    col_annotation <- NULL
  } else{
    col_annotation <- metadata_col %>%
      dplyr::select(Sample_ID, all_of(proj.params$heatmap$col.ann)) %>%
      dplyr::mutate(Sample_ID = make.names(Sample_ID)) %>%
      dplyr::filter(Sample_ID %in% colnames(mat)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("Sample_ID") %>%
      as.data.frame() %>%
      mutate(across(where(is.factor), as.character))
  }
  
  # Define row annotations
  if(is.null(metadata_row)) {
    row_annotation <- NULL
  } else{
    row_annotation <- metadata_row %>%
      dplyr::select(SYMBOL, all_of(proj.params$heatmap$row.ann)) %>%
      dplyr::mutate(SYMBOL = make.names(SYMBOL, unique = TRUE)) %>%
      dplyr::filter(SYMBOL %in% rownames(mat)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("SYMBOL") %>%
      as.data.frame() %>%
      mutate(across(where(is.factor), as.character))
  }
  
  # ---- Annotation Palettes ----
  # Define Color Palette for Annotation 
  # This is an example of how ann_colors should be specified
  ann_colors <- list(CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
                     GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E"))
  
  ann_colors <- list()
  base_colors <- c("#E08214", "#762A83", "#C51B7D", "#7FBC41", "#35978F", "#BF812D", "#542788",
                   "#D6604D", "#4393C3", "#878787", "#E41A1C", "#F781BF", "#4DAF4A", "#FFFFBF",
                   "#377EB8", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#999999", "#66C2A5",
                   "#FC8D62", "#000000", "#9E0142", "#1A1A1A", "#74006F", "#FFC606", "#F6D2E0",
                   "#C8E7F5")
  base_colors <- custom_palette
  
  col_list <- base::lapply(X = as.list(col_annotation), FUN = function(x) { as.character(x) %>% unique})
  row_list <- base::lapply(X = as.list(row_annotation), FUN = function(x) { as.character(x) %>% unique})
  ann_list <- c(row_list, col_list)
  
  color_index <- 1
  for (i in seq_along(ann_list)) {  # Iterate through each annotation variable (Eg: CellType) 
    levels <- sort(ann_list[[i]])   # Get levels within each annotation variable (Eg: CT1, CT2)
    n_levels <- length(levels)      # Get number of levels within each annotation variable
    
    palette_colors <- if (proj.params$heatmap$ann.palette == "discrete" | n_levels == 1){
      base_colors[color_index:(color_index + n_levels - 1)]
    } else{
      alphas <- seq(1 / n_levels, 1, length.out = n_levels)
      base::sapply(X = alphas, 
                   FUN = function(x) { colorspace::adjust_transparency(col = base_colors[color_index], alpha = x) })
    }
    
    names(palette_colors) <- levels                   # Name each color with levels
    ann_colors <- c(ann_colors, list(palette_colors)) # Append named color palette
    names(ann_colors)[i] <- names(ann_list)[i]        # Name the color palette with corresponding annotation variable name
    color_index <- color_index + n_levels             # Move to next color
  }
  
  # ---- Heatmap Palette ----
  # Define Color Palette for Heatmap 
  valid_palettes <- c("vrds", "rdbu")
  
  if (!proj.params$heatmap$palette %in% valid_palettes) {
    stop("Invalid heatmap palette. proj.params$heatmap$palette must be either 'vrds' or 'rdbu'.")
  }
  
  heatmap_palette <- switch(proj.params$heatmap$palette,
                            vrds = viridis::viridis(100),
                            rdbu = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100))
  
  # ---- Breaks ----
  # Define Color Breaks 
  n_breaks <- 100
  heatmap_palette <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(n_breaks)
  
  # Handle min and max thresholds with soft clamping
  mat_min <- min(mat_scaled, na.rm = TRUE)
  mat_max <- max(mat_scaled, na.rm = TRUE)
  mat_min <- dplyr::case_when(mat_min >= 0 ~ 0,
                              mat_min <= -3 ~ -3,
                              TRUE ~ mat_min)
  mat_max <- dplyr::case_when(mat_max <= 0 ~ 0,
                              mat_max >= 3 ~ 3,
                              TRUE ~ mat_max)
  
  if (mat_max == 0){
    breaks <- seq(from = floor(mat_min), to = 0, length.out = n_breaks)
  } else if (mat_min == 0){
    breaks <- seq(from = 0, to = ceiling(mat_max), length.out = n_breaks)
  } else{
    breaks <- c(seq(from = floor(mat_min),   to = 0,        length.out = n_breaks / 2),
                seq(from = mat_max / n_breaks, to = ceiling(mat_max), length.out = n_breaks / 2))
  }
  
  # ---- Gaps ----
  # Define gaps in heatmap 
  gaps_col <- if (!gtools::invalid(proj.params$heatmap$col.gaps) & proj.params$heatmap$col.cluster %in% colnames(col_annotation)) {
    if (all(proj.params$heatmap$col.gaps %in% colnames(col_annotation))) {
      col_annotation %>%
        dplyr::count(get(proj.params$heatmap$col.gaps)) %>%
        dplyr::mutate(n = cumsum(n)) %>%
        dplyr::pull(n) %>%
        .[. < ncol(mat_scaled)]
    } else{
      message("Gaps not introduced as gap identifiers missing in column annotation")
    }
  }
  
  gaps_row <- if (!gtools::invalid(proj.params$heatmap$row.gaps) & proj.params$heatmap$row.cluster %in% colnames(row_annotation)) {
    if (all(proj.params$heatmap$row.gaps %in% colnames(row_annotation))) {
      row_annotation %>%
        dplyr::count(get(proj.params$heatmap$row.gaps)) %>%
        dplyr::mutate(n = cumsum(n)) %>%
        dplyr::pull(n) %>%
        .[. < nrow(mat_scaled)]
    } else{
      message("Gaps not introduced as gap identifiers missing in row annotation")
    }
  }
  
  # ---- Clustering ----
  # Determine Ordering 
  if (proj.params$heatmap$col.cluster == "all"){
    colclust <- hclust(dist(t(mat_scaled)))
    col_order <- colnames(mat_scaled)[colclust$order]
  }
  if (proj.params$heatmap$col.cluster == "alphabetical"){
    col_order <- sort(colnames(mat_scaled))
  }
  if (proj.params$heatmap$col.cluster %in% colnames(col_annotation)){
    
    # NOTE: While calculating gaps_col, we use count(). It sorts alphabetically.
    # So, WE MUST sort col_elements to match gaps_col
    col_order <- c()
    col_elements <- col_annotation %>% 
      dplyr::pull(proj.params$heatmap$col.cluster) %>%
      unique() %>% sort()
    
    for (g in col_elements){
      
      samples <- rownames(col_annotation)[col_annotation %>% dplyr::pull(proj.params$heatmap$col.cluster) == g]
      if (length(samples) == 0) next
      temp_mat <- mat_scaled[, samples]
      
      if (length(samples) > 1){
        colclust <- hclust(dist(t(temp_mat)))
        col_order <- c(col_order, colnames(temp_mat)[colclust$order])
      } else if(length(samples) == 1){
        col_order <- c(col_order, samples)
      }
    }
  }
  
  if (proj.params$heatmap$row.cluster == "all"){
    rowclust <- hclust(dist(mat_scaled))
    row_order <- rownames(mat_scaled)[rowclust$order]
  }
  if (proj.params$heatmap$row.cluster == "alphabetical"){
    row_order <- sort(rownames(mat_scaled))
  }
  if (proj.params$heatmap$row.cluster %in% colnames(row_annotation)){
    
    # NOTE: While calculating gaps_row, we use count(). It sorts alphabetically.
    # So, WE MUST sort row_elements to match gaps_row
    row_order <- c()
    row_elements <- row_annotation %>% 
      dplyr::pull(proj.params$heatmap$row.cluster) %>%
      unique() %>% sort()
    
    for (g in row_elements){
      
      genes <- rownames(row_annotation)[row_annotation %>% dplyr::pull(proj.params$heatmap$row.cluster) == g]
      if (length(genes) == 0) next
      temp_mat <- mat_scaled[rownames(mat_scaled) %in% genes,]
      
      if (length(genes) > 1){
        rowclust <- hclust(dist(temp_mat))
        row_order <- c(row_order, rownames(temp_mat)[rowclust$order])
      } else if(length(genes) == 1){
        row_order <- c(row_order, genes)
      }
    }
  }
  
  # ---- Prepare Matrix for Plotting ---- 
  reordered <- mat_scaled[row_order, col_order]
  
  # ---- Set font sizes ----
  fontsize <- 10
  fontsize.row <- fontsize*1
  fontsize.col <- fontsize*1
  fontsize.number <- fontsize*0.8
  angle.col <- 45       # column label angle
  
  # ---- Set cell width, height (in points) dynamically ----
  
  # NOTE: A page is 8.5inch wide x 11 inch high
  # Giving 2 inch width for figure legend, 2 inch height for col annotations,
  # 0.5 inch for margins, actual plot area ~ 6inch wide x 8 inch high
  # For text to be readable, fontsize >= 10. We keep cell.width and cell.height
  # 5 point size larger than fontsize so it looks pretty
  
  optimal_col.width <- fontsize+5
  optimal_row.height <- fontsize+5
  max_cols <- 6 * 72 / optimal_col.width
  max_rows <- 8 * 72 / optimal_row.height
  
  if (ncol(reordered) <= max_cols){
    cell.width <- optimal_col.width
  } else{
    cell.width <- 6 * 72 / ncol(reordered)
  }
  
  if (nrow(reordered) <= max_rows){
    cell.height <- optimal_row.height
  } else{
    cell.height <- 8 * 72 / nrow(reordered)
  }
 
  # ---- Truncate long row and column labels ----
  main.title <- stringr::str_wrap(string = proj.params$heatmap$title, width = 20)
  if (ncol(reordered) <= max_cols){
    labels.col <- stringr::str_trunc(string = colnames(reordered), width = 15, side = "right", ellipsis = "…")
  } else{
    labels.col <- rep(x="", times = ncol(reordered))
  }
  if (nrow(reordered) <= max_rows & length(disp_genes) > 0){
    labels.row <- stringr::str_trunc(string = rownames(reordered), width = 15, side = "right", ellipsis = "…")
    labels.row <- dplyr::if_else(rownames(reordered) %in% make.names(disp_genes), rownames(reordered), " ")
  } else {
    labels.row <- rep(x="", times = nrow(reordered))
  }
  
  # ---- Heatmap plotting ----
  ph <- pheatmap::pheatmap(mat               = reordered,
                           color             = heatmap_palette,
                           breaks            = breaks,
                           annotation_row    = row_annotation,
                           annotation_col    = col_annotation,
                           annotation_colors = ann_colors,
                           gaps_row          = gaps_row,
                           gaps_col          = gaps_col,
                           
                           cellwidth         = cell.width,     
                           cellheight        = cell.height,  
                           show_rownames     = cell.height >= fontsize,
                           show_colnames     = cell.width >= fontsize,
                           labels_row        = labels.row,
                           labels_col        = labels.col,
                           angle_col         = angle.col,        # column label angle
                           fontsize          = fontsize,         # points; 72 points = 1 inch
                           fontsize_row      = fontsize.row,     # points
                           fontsize_col      = fontsize.col,     # points
                           fontsize_number   = fontsize.number,  # points
                           silent            = TRUE, 
                           
                           main              = main.title,
                           border_color      = proj.params$heatmap$border.color,
                           legend            = proj.params$heatmap$show.expr.legend,
                           
                           scale                    = "none",
                           cluster_rows             = FALSE,
                           cluster_cols             = FALSE,
                           clustering_distance_rows = "correlation", #"euclidean",
                           clustering_distance_cols = "correlation", #"euclidean",
                           clustering_method        = "average",     #complete",
                           annotation_legend        = TRUE,
                           annotation_names_row     = FALSE,
                           annotation_names_col     = FALSE,
                           width                    = NA,               # inches
                           height                   = NA,               # inches
                           filename                 = NA)
  
  # ---- Return matrix for Excel ----
  ph_mat <- if(ncol(reordered) > nrow(reordered)) t(reordered) else reordered
  
  return(invisible(list(ph = ph, mat = ph_mat)))
}

run_deseq2 <- function(meta_data, read_data, proj.params, n = 1) {
  
  set.seed(1234)
  
  # ---- Input Checks ----
  if (!is.list(proj.params) || !"deseq2" %in% names(proj.params)) {
    stop("⚠️ proj.params must contain a 'deseq2' list")
  }
  
  required_attrs <- c("design", "lfc.cutoff", "padj.cutoff", "contrasts")
  for (attr in required_attrs) {
    if (is.null(proj.params$deseq2[[attr]])) {
      stop("⚠️ Missing required proj.params$deseq2 attribute: ", attr)
    } else {
      val <- proj.params$deseq2[[attr]]
      message("✔ Found proj.params$deseq2$", attr, " = ", paste(val, collapse = ", "))
    }
  }
  
  contrast <- proj.params$deseq2$contrasts[n]
  contrast_dir <- proj.params$contrast_dir[n]
  
  if (!dir.exists(contrast_dir)) {
    warning("Contrast directory does not exist. Creating: ", contrast_dir)
    dir.create(contrast_dir, recursive = TRUE)
  }
  
  # ---- DESeq2 Object Preparation ----
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = read_data,
                                        colData   = meta_data,
                                        design    = ~1)
  design(dds) <- as.formula(paste0("~", proj.params$deseq2$design))
  
  # ---- Size Factor Estimation ----
  if (all(rowSums(read_data == 0) > 0)) {
    message("Cannot compute log geometric means as every gene contains at least
            one zero. Using 'poscounts' for size factor estimation.")
    dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
  }
  
  # ---- Pre-filter Lowly Expressed Genes ----
  # NOTE: This improves sizefactor estimation
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  # ---- DESeq2 Fit: Parametric and Local ----
  dds_para <- DESeq2::DESeq(object = dds, test = "Wald", fitType = "parametric",
                            betaPrior = FALSE, minReplicatesForReplace = 7)
  dds_local <- DESeq2::DESeq(object = dds, test = "Wald", fitType = "local",
                             betaPrior = FALSE, minReplicatesForReplace = 7)
  
  # ---- Select Best Fit Based on Residuals ----
  residual_para  <- mcols(dds_para)$dispGeneEst - mcols(dds_para)$dispFit
  residual_local <- mcols(dds_local)$dispGeneEst - mcols(dds_local)$dispFit
  dds <- if (median(residual_para^2, na.rm = TRUE) <= median(residual_local^2, na.rm = TRUE)) {
    dds_para
  } else {
    dds_local
  }
  
  # ---- Prepare Contrast Vector ----
  mod_mat <- model.matrix(design(dds), colData(dds))
  design_factors <- stringr::str_split(string = proj.params$deseq2$design, 
                                       pattern = "[+*:]")[[1]] %>% unique()
  
  # Create a replicate of meta_data with all possible groups that could be compared based on the design
  df <- colData(dds) %>%
    as.data.frame() %>%
    tidyr::unite(col = "Groups", all_of(design_factors), sep = ".")
  
  # Define all possible groups that could be compared based on the design
  groups <- unique(df$Groups)
  
  # Get all possible coefficient vectors
  group_coef_list <- lapply(groups, function(i) colMeans(as.matrix(mod_mat[df$Groups == i, , drop = FALSE])))
  names(group_coef_list) <- groups
  
  target <- stringr::str_split(string = contrast, pattern = "-")[[1]][1]
  ref    <- stringr::str_split(string = contrast, pattern = "-")[[1]][2]
  
  # Below line works for simple contrasts like (A - B). Fails if contrast is (A-B)-(C-D)
  # contrast_vec <- group_coef_list[[target]] - group_coef_list[[ref]] 
  
  replace_symbols <- function(node) {
    if (is.symbol(node)) {
      nm <- as.character(node)
      if (nm %in% names(group_coef_list)) {
        return(group_coef_list[[nm]])
      } else {
        # it's an operator like "-" or "+" → return unchanged
        return(node)
      }
    } else if (is.call(node)) {
      return(as.call(lapply(node, replace_symbols)))
    } else {
      return(node)
    }
  }
  
  parsed <- base::parse(text = contrast)[[1]]
  expr_sub <- replace_symbols(parsed)
  contrast_vec <- base::eval(expr_sub)
  
  # # [ALTERNATIVE METHOD] Get all possible coefficient vectors
  # for (i in groups) {
  #   val <- colMeans(as.matrix(mod_mat[df$Groups == i, , drop =FALSE]))
  #   assign(x = i, value = val)
  # }
  # # Define the contrast vector
  # contrast_vec <- base::eval(expr = base::parse(text = contrast))
  
  # ---- DESeq2 Results with LFC Threshold & Shrinkage ----
  res <- DESeq2::results(object = dds, 
                         contrast = contrast_vec,
                         lfcThreshold = proj.params$deseq2$lfc.cutoff,
                         altHypothesis = "greaterAbs",
                         cooksCutoff = TRUE,
                         independentFiltering = TRUE,
                         alpha = proj.params$deseq2$padj.cutoff,
                         pAdjustMethod = "BH")
  
  # Safe LFC shrinkage
  set.seed(1234)
  res <- DESeq2::lfcShrink(dds = dds, res = res, type = "ashr")
  summary(res)
  
  # ---- Differential Expressed Genes (DEGs) ----
  DEGs_df <- res %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    add_annotation() %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarize(across(.cols = where(is.numeric), .fns = mean, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(padj = ifelse(padj == 0, min(padj[padj > 0], na.rm = TRUE), padj))
  
  save_xlsx(DEGs_df, file.path(contrast_dir, "DEGs.xlsx"), "DEGs", row_names = FALSE)
  
  # ---- VST Counts (Non-blind) ----
  vsd <- DESeq2::vst(dds, blind = FALSE)
  vst_counts <- SummarizedExperiment::assay(vsd) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    add_annotation() %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarize(across(.cols = where(is.numeric), .fns = mean, na.rm = TRUE), .groups = "drop") %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("SYMBOL") %>%
    dplyr::select(-dplyr::starts_with("ENSEMBL"), -dplyr::starts_with("ENTREZ")) %>%
    as.matrix()
  
  save_xlsx(vst_counts, file.path(contrast_dir, "VST_counts.xlsx"), "VST_Nonblind", row_names = TRUE)
  
  # ---- Return Results ----
  invisible(list(degs = DEGs_df, vst = vst_counts, dds = dds))
}

# DEGs_df with column SYMBOL, padj, log2FoldChange
# k <- # overlapping genes between input and pathway
# n <- # overlapping genes between input and collection
# K <- # genes in pathway
# N <- # genes in collection
pathway_analysis <- function(DEGs_df, proj.params, output_dir) {
  
  # Check required proj.params attributes
  required_attrs <- c("gmt_dir", "species")
  for (attr in required_attrs) {
    if (is.null(proj.params[[attr]])) {
      stop("⚠️ Missing required proj.params attribute: ", attr)
    } else {
      val <- proj.params[[attr]]
      message("✔ Found proj.params attribute: ", attr, " = ", val)
    }
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  set.seed(1234)
  
  gmt_files <- list.files(file.path(proj.params$gmt_dir, proj.params$species), full.names = TRUE)
  
  # Initialize result dataframes 
  fgsea_df <- data.frame()
  gsea_df <- data.frame()
  ora_df_up <- data.frame()
  ora_df_down <- data.frame()
  concise_fgsea_df <- data.frame()
  
  # Define input genes for GSEA (Ranked list of all genes) 
  # IMPORTANT: Rank genes from high LFC to low lFC, so +NES ~ up-regulated
  ranked_df <- DEGs_df %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    dplyr::filter(!is.na(padj), !is.na(SYMBOL)) %>%
    dplyr::mutate(log2FoldChange = as.numeric(log2FoldChange),
                  padj = as.numeric(padj)) %>%
    dplyr::arrange(desc(log2FoldChange))  
  
  ranked_list <- ranked_df$log2FoldChange
  names(ranked_list) <- ranked_df$SYMBOL
  
  # Define input and universe genes for ORA (significant genes only) 
  sig_genes_up <- DEGs_df %>%
    dplyr::filter(padj <= 0.05, log2FoldChange > 0) %>%
    dplyr::pull(SYMBOL)
  
  sig_genes_down <- DEGs_df %>%
    dplyr::filter(padj <= 0.05, log2FoldChange < 0) %>%
    dplyr::pull(SYMBOL)
  
  universe_genes <- unique(DEGs_df$SYMBOL)
  
  # Iterate over GMT files and perform enrichment analyses
  for (gmt_file in gmt_files) {
    
    # Extract gene set name
    # gmt_name <- gsub(pattern = "^.*/|.v[0-9].*$", replacement = "", x = gmt_file)
    gmt_name <- gsub(pattern = "^.*/|", replacement = "", x = gmt_file)
    
    # Format gene sets for fgsea and keep only genes present in ranked_list
    gmt <- fgsea::gmtPathways(gmt_file)
    gmt <- lapply(X = gmt, FUN = function(x){x[x %in% names(ranked_list)]})
    
    # Format gene sets for clusterProfiler and keep only genes present in ranked_list
    pathway_gene_df <- data.frame(pathways = base::rep(x = names(gmt), times = base::unname(obj = lengths(gmt))),
                                  genes = unlist(gmt, use.names = FALSE))
    
    # Run fgseaMultilevel (GSEA)
    fgsea_res <- fgsea::fgseaMultilevel(pathways = gmt,
                                        stats = ranked_list,
                                        scoreType = dplyr::case_when(min(ranked_list) > 0 ~ "pos",
                                                                     max(ranked_list) < 0 ~ "neg",
                                                                     TRUE ~ "std"),
                                        sampleSize = 101,
                                        minSize = 1,
                                        maxSize = 500, # recommended 500 genes max
                                        eps = 1e-50,
                                        nproc = 0,
                                        gseaParam = 1,
                                        BPPARAM = NULL,
                                        nPermSimple = 10000)
    
    # Identify overlapping pathways and collapse into major pathways
    concise_fgsea_res <- fgsea::collapsePathways(fgseaRes = fgsea_res,
                                                 pathways = gmt,
                                                 stats = ranked_list)
    concise_fgsea_res <- fgsea_res %>%
      dplyr::filter(pathway %in% concise_fgsea_res$mainPathways)
    
    # Run clusterProfiler GSEA
    gsea_res <- clusterProfiler::GSEA(geneList = ranked_list,
                                      exponent = 1,
                                      minGSSize = 10,
                                      maxGSSize = 500,
                                      eps = 1e-10,
                                      pvalueCutoff = 0.05,
                                      pAdjustMethod = "BH",
                                      TERM2GENE = pathway_gene_df,
                                      TERM2NAME = NA,
                                      verbose = FALSE,
                                      seed = FALSE,
                                      by = "fgsea")
    
    # Run clusterProfiler ORA (enricher)
    # NOTE: Avoid using clusterProfiler::enrichGO() as it doesnt use proper 
    # background in universe parameter and includes GO terms outside the 
    # intended gene set collection.
    ora_res_up <- clusterProfiler::enricher(gene = sig_genes_up,
                                            pvalueCutoff = 0.05,
                                            pAdjustMethod = "BH",
                                            universe = universe_genes,
                                            minGSSize = 10,
                                            maxGSSize = 500,
                                            qvalueCutoff = 0.2,
                                            TERM2GENE = pathway_gene_df,
                                            TERM2NAME = NA)
    
    ora_res_down <- clusterProfiler::enricher(gene = sig_genes_down,
                                              pvalueCutoff = 0.05,
                                              pAdjustMethod = "BH",
                                              universe = universe_genes,
                                              minGSSize = 10,
                                              maxGSSize = 500,
                                              qvalueCutoff = 0.2,
                                              TERM2GENE = pathway_gene_df,
                                              TERM2NAME = NA)
    
    # Bind significant results to respective dataframes
    fgsea_df <- dplyr::bind_rows(fgsea_df, fgsea_res)
    concise_fgsea_df <- dplyr::bind_rows(concise_fgsea_df, concise_fgsea_res)
    
    if (!is.null(gsea_res)) {
      gsea_df <- dplyr::bind_rows(gsea_df, gsea_res@result)
    }
    
    if (!is.null(ora_res_up)) {
      ora_df_up <- dplyr::bind_rows(ora_df_up, ora_res_up@result)
    }
    
    if (!is.null(ora_res_down)) {
      ora_df_down <- dplyr::bind_rows(ora_df_down, ora_res_down@result)
    }
  } 
  
  # Format the results
  # Rename columns consistently across different methods
  lookup <- c(pathway = "ID",
              geneID = "leadingEdge", geneID = "core_enrichment",
              K = "size", K = "setSize",
              padj = "p.adjust", 
              pval = "pvalue")
  
  ora_df <- dplyr::bind_rows(ora_df_up %>% dplyr::mutate(Direction = "Upregulated"), 
                             ora_df_down %>% dplyr::mutate(Direction = "Downregulated"))
  if (nrow(ora_df) > 0){
    ora_df <- ora_df %>%
      tidyr::separate(col = GeneRatio, into = c("k", "n")) %>%
      tidyr::separate(col = BgRatio, into = c("K", "N")) %>%
      dplyr::mutate_at(c("k", "n", "K", "N"), as.numeric) %>%
      dplyr::mutate(GeneRatio = k / n,
                    BackgroundRatio = K / N,
                    EnrichmentRatio = GeneRatio / BackgroundRatio,
                    combined_score = GeneRatio * -log10(p.adjust),
                    NES = NA_integer_) %>%
      dplyr::rename(any_of(lookup))
  }
  
  for (i in c("fgsea_df", "gsea_df")){
    df <- get(i) %>% 
      dplyr::mutate(Direction = dplyr::case_when(NES > 0 ~ "Upregulated",
                                                 NES < 0 ~ "Downregulated",
                                                 TRUE ~ "No change")) %>%
      dplyr::rename(any_of(lookup))
    assign(x = i, value = df)
  }
  
  for (i in c("fgsea_df", "gsea_df", "ora_df")){
    
    df <- get(i)
    if (nrow(df) > 0){
      df <- df %>%
        dplyr::filter(padj <= 0.05) %>%
        tibble::remove_rownames() %>%
        tidyr::separate(col = pathway, into = c("Collection", "Description"), sep = "_", extra = "merge") %>%
        dplyr::mutate(Description = base::gsub(pattern = "_", replacement = " ", x= Description),
                      geneID = base::sapply(X = geneID, FUN = paste, collapse = "/")) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(leading_edge_size = length(unlist(stringr::str_split(geneID, "/")))) %>%
        dplyr::ungroup() %>%
        as.data.frame() %>%
        dplyr::select(Collection, Description, leading_edge_size, K, pval, padj, NES, Direction, everything(), -geneID, geneID)
      
      max_len <-  max(df$leading_edge_size, na.rm = TRUE)
      sep_char <- ifelse(grepl(",", df$geneID[1], fixed = TRUE), ",", "/")
      
      # across() allows you to apply a function(s) to multiple columns at once. 
      # map_chr() iterates over a list/vector and applies a function to each 
      # element (in this context, each cell within the selected columns),
      # returning a character vector.
      if (is.finite(max_len) & max_len > 0){
        df <- df %>% 
          tidyr::separate(col = geneID, into = paste0("gene", 1:max_len), sep = sep_char, remove = TRUE, fill = "right") %>%
          dplyr::mutate(across(.cols = starts_with("gene", ignore.case = FALSE), 
                               .fns = function (col) { col %>% 
                                   purrr::map_chr(.f = function (x){gsub('c\\(', "", x) }) %>%
                                   purrr::map_chr(.f = function (x){gsub('\\)', "", x) }) %>%
                                   purrr::map_chr(.f = function (x){gsub('"', "", x) }) %>%
                                   trimws() })) %>%
          dplyr::select(Collection, Description, leading_edge_size, K, padj, NES, Direction, everything())
      }
      assign(x = i, value = df)
    }
  }
  
  # Summarize results from fgsea, gsea and ora 
  consensus_df <- dplyr::bind_rows(fgsea_df %>% dplyr::mutate(method = "FGSEA"), 
                                   gsea_df %>% dplyr::mutate(method = "GSEA"), 
                                   ora_df %>% dplyr::mutate(method = "ORA")) %>%
    dplyr::add_count(Collection, Description, Direction, name = "n_methods") %>%
    dplyr::filter(n_methods > 1) %>%
    dplyr::mutate(Consensus = Direction) %>%
    dplyr::arrange(Collection, Description, desc(NES)) %>%
    dplyr::select(n_methods, method, Consensus, Collection, Description, 
                  leading_edge_size, K, padj, NES, Direction, everything(), 
                  -starts_with("gene",ignore.case = FALSE),
                  starts_with("gene", ignore.case = FALSE)) 
  
  gsea_res_list <- list(consensus = consensus_df,
                        fgsea = fgsea_df, 
                        gsea = gsea_df, 
                        ora = ora_df)
  # Save results
  wb <- openxlsx::createWorkbook()
  for (i in seq_along(gsea_res_list)) {
    openxlsx::addWorksheet(wb, sheetName = names(gsea_res_list)[i])
    openxlsx::writeData(wb, sheet = names(gsea_res_list)[i], x = gsea_res_list[[i]], rowNames = FALSE)
  }
  openxlsx::saveWorkbook(wb, file.path(output_dir, "Pathway_results.xlsx"), overwrite = TRUE)
  
  # Return results
  return(invisible(list(fgsea = fgsea_df, 
                        gsea = gsea_df, 
                        ora = ora_df,
                        consensus = consensus_df)))
}

### decoupleR Analysis [RECOMMENDED as it can run multiple algorithms & give consensus]
# input can be vst_counts or DEGs with t-statistic (-log10padj*log2FC)
tf_analysis <- function(input, species = "Homo sapiens", top = 500) {
  
  set.seed(1234)
  
  # --- Input checks ---
  if (!is.matrix(input)) stop("`input` must be a matrix.")
  if (!species %in% c("Homo sapiens", "Mus musculus")) {
    stop("`species` must be either 'Homo sapiens' or 'Mus musculus'.")
  }
  
  # --- Map species to decoupleR format ---
  organism <- dplyr::case_when(species == "Homo sapiens" ~ "human",
                               species == "Mus musculus" ~ "mouse",
                               TRUE ~ "rat")
  
  # Load network models
  # get_collectri returns mor column without edge weights (-1 or +1)
  # get_dorothea returns mor column with edge weights (-1 through +1)
  # viper works better with edge weights but dorothea's edge weights are solely based on confidence
  # get_collectri has only high confidence entries. So, get_collecttri is better.
  progeny_pathway_net <- decoupleR::get_progeny(organism = organism, top = top)
  collectri_tf_net <- decoupleR::get_collectri(organism = organism, split_complexes = FALSE)
  #dorothea_tf_net <- decoupleR::get_dorothea(organism = organism, levels = c('A', 'B', 'C'))
  
  # Run decoupleR
  # stats <- c("aucell", fgsea", "gsva", "mdt", "mlm", "ora", "udt", "ulm", "viper", "wmean", "wsum")
  stats <- c("ulm", "mlm", "viper") # 3 best methods to get consensus
  pathway_df <- decoupleR::decouple(mat = input, network = progeny_pathway_net,
                                    statistics = stats, minsize = 5)
  
  tf_df <- decoupleR::decouple(mat = input, network = collectri_tf_net,
                               statistics = stats, minsize = 5)
  
  # tf_df <- decoupleR::decouple(mat = input, network = dorothea_tf_net,
  #                              statistics = "viper", minsize = 5)
  
  # Remove insignificant entries
  pathway_sig <- pathway_df %>%
    dplyr::group_by(statistic) %>%
    dplyr::mutate(padj = stats::p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup() %>% 
    dplyr::filter(padj <= 0.05)
  
  tf_sig <- tf_df %>%
    dplyr::group_by(statistic) %>%
    dplyr::mutate(padj = stats::p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup() %>% 
    dplyr::filter(padj <= 0.05)
  
  # Return results 
  result <- list(all_pathways = pathway_df,
                 all_tfs = tf_df,
                 sig_pathways = pathway_sig,
                 sig_tf = tf_sig)
  
  return(invisible(result))
} 

# No "_" in Description column so that str_wrap works
# Collection column needed
# (combined_score, GeneRatio, k) OR (NES, leading_edge_size) columns needed
plot_pathway <- function(pathway_df, vst_counts, meta_data, samples, method, output_dir){
  
  set.seed(1234)
  
  if (nrow(pathway_df) == 0){
    message("Input data frame is empty. Skipping plotting.")
    return(NULL)
  }
  
  if (is.null(method) || method == "") {
    stop("⚠️ 'method' must be defined before creating output directories.")
  }
  
  # ---- Create output directories ----
  bar_output_dir <- file.path(output_dir, method, "Bar_Plots")
  dot_output_dir <- file.path(output_dir, method, "Dot_Plots")
  heatmap_output_dir <- file.path(output_dir, method, "Heatmaps")
  dirs <- c(output_dir,
            bar_output_dir,
            dot_output_dir,
            heatmap_output_dir)
  
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
      message("📂 Created directory: ", d)
    }
  }
  
  # ----   ----
  dot_plot_list <- list()
  bar_plot_list <- list()
  plot_colors <- c("Upregulated" = "#E69F00", "Downregulated" = "#56B4E9")
  
  # str_wrap() wraps only between words, and it defines "words" based on spaces (" ").
  # If all words are connected by "_", it wont split.
  pathway_df <- pathway_df %>%
    dplyr::mutate(Description = gsub(pattern = "_", replacement = " ", x = Description),
                  Description = stringr::str_wrap(string = Description, width = 30))
  
  # Decide if there are multiple collections or single collection.
  # If multiple, plot each collection separately,
  collections <- unique(pathway_df$Collection)
  n_collections <- length(unique(pathway_df$Collection))
  
  for (i in seq_len(n_collections)){
    
    # Save heatmap for pathways in each collection
    heatmap_plot_list <- list()
    
    plot_df <- pathway_df %>% dplyr::filter(Collection == collections[i])
    
    # Choose x axis column, size column and labels
    if (method == "ORA"){
      x_col <- sym("GeneRatio")
      x_label <- "GeneRatio"
      size_col <- sym("leading_edge_size")
      color_col <- sym("Direction")
      score_var <- dplyr::case_when("GeneRatio" %in% colnames(pathway_df) ~ "GeneRatio",
                                    "combined_score" %in% colnames(pathway_df) ~ "combined_score", 
                                    TRUE ~ NA_character_)
      x_limits <- c(0, NA)  # start at 0, auto end
      
      pathway_df <- pathway_df %>%
        dplyr::filter(!is.na(.data[[score_var]])) %>%
        dplyr::arrange(dplyr::desc(.data[[score_var]]))
      
    } else if (method == "GSEA"){
      x_col <- sym("NES")
      x_label <- "Normalized Enrichment Score (NES)"
      size_col <- sym("leading_edge_size")
      color_col <- sym("Direction")
      score_var <- "NES"
      x_min <- ifelse(floor(min(plot_df$NES, na.rm = TRUE)) > 0, 0, floor(min(plot_df$NES, na.rm = TRUE)))
      x_limits <- c(x_min, NA)
      
    }
    
    # Pad with empty rows if fewer than 15 pathways 
    n_missing <- 20 - nrow(plot_df)
    if(n_missing > 0){
      empty_df <- matrix(data = "", 
                         nrow = n_missing,
                         ncol = ncol(plot_df)) %>%
        as.data.frame() %>%
        dplyr::mutate(Description = paste0("", seq_len(n_missing)))
      plot_df <- dplyr::bind_rows(plot_df, empty_df)
    }
    
    max_label_len <- max(nchar(plot_df$Description), na.rm = TRUE)
    y_text_size <- dplyr::case_when(max_label_len > 50 ~ 6,
                                    max_label_len > 35 ~ 7,
                                    max_label_len > 25 ~ 8,
                                    TRUE ~ 10)
    
    # ---- Plot bar plot ---- 
    p1 <- ggplot2::ggplot(data = plot_df,
                          aes(x = !!x_col,
                              y = reorder(Description, !!x_col),
                              fill = !!color_col,
                              alpha = -log10(padj))) +
      ggplot2::geom_col(width = 0.75, na.rm = TRUE) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = x_label,
                    y = "",
                    title = paste("Top", collections[i], "Pathways"),
                    fill = "Direction") +
      custom_theme +
      theme(axis.text.y = element_text(size = y_text_size)) +
      coord_cartesian(clip = "off") +
      scale_x_continuous(limits = x_limits, expand = expansion(mult = c(0, 0.05))) +
      scale_alpha_continuous(range = c(0.5, 1)) +
      scale_fill_manual(values = plot_colors) +
      guides(fill = guide_legend(override.aes = list(shape = 22, size = 6)),
             color = guide_legend(override.aes = list(shape = 22, size = 6)),
             alpha = guide_legend(override.aes = list(shape = 22, size = 6))) +
      ggplot2::geom_text(aes(label =  !!size_col), x = 0, hjust = -0.1, size = 3, show.legend = FALSE)
    
    bar_plot_list <- c(bar_plot_list, list(p1))
    
    # ---- Plot dot plot ---- 
    vals <- c(min(plot_df[[size_col]], na.rm = TRUE), max(plot_df[[size_col]], na.rm = TRUE))
    breaks <- as.vector(floor(quantile(vals) / 10) * 10)
    
    p2 <- ggplot2::ggplot(data = plot_df,
                          aes(x = !!x_col,
                              y = reorder(Description, !!x_col),
                              fill = !!color_col,
                              alpha = -log10(padj),
                              color = !!color_col,
                              size = !!size_col)) +
      ggplot2::geom_point() +
      ggplot2::theme_classic() +
      ggplot2::labs(x = x_label ,
                    y = "",
                    title = paste("Top", collections[i], "Pathways"),
                    color = "Direction",
                    size = "Counts") +
      custom_theme +
      theme(axis.text.y = element_text(size = y_text_size)) +
      coord_cartesian(clip = "off") + 
      scale_x_continuous(limits = x_limits, expand = expansion(mult = c(0, 0.05))) +
      scale_alpha_continuous(range = c(0.5, 1)) +
      scale_color_manual(values = plot_colors) +
      scale_fill_manual(values = plot_colors) + # need for coloring the legend
      guides(fill = guide_legend(override.aes = list(shape = 22, size = 6)),
             color = guide_legend(override.aes = list(shape = 22, size = 6)),
             alpha = guide_legend(override.aes = list(shape = 15, size = 6))) +
      ggplot2::scale_size(breaks = breaks) 
    
    dot_plot_list <- c(dot_plot_list, list(p2))
    
    # ---- Plot heatmap ----
    
    pathways <- pathway_df %>% 
      dplyr::filter(Collection == collections[i]) %>%
      dplyr::pull(Description) %>%
      unique()
    
    for (pathway in pathways) {
      plot_genes <- pathway_df %>% 
        dplyr::filter(Collection == collections[i], Description == pathway) %>%
        dplyr::select(dplyr::starts_with("gene", ignore.case = FALSE)) %>%
        unlist(use.names = FALSE) %>%
        trimws() %>%
        na.omit() %>% 
        unique()
      
      if (length(plot_genes) >= 2){
        norm_counts <- vst_counts[plot_genes, samples, drop = FALSE]
        metadata_col <- meta_data %>% dplyr::filter(Sample_ID %in% samples)
        metadata_row <- NULL
        disp_genes <- c()
        proj.params$heatmap$title <- stringr::str_wrap(string = pathway, width = 30)
        
        ph <- plot_heatmap(norm_counts, proj.params, metadata_col, metadata_row, disp_genes)
        heatmap_plot_list[[length(heatmap_plot_list) + 1]] <- ph$ph$gtable
      }
    }
    
    if (length(heatmap_plot_list) > 0) {
      pdf(file.path(heatmap_output_dir, paste0("Heatmaps_", collections[i], "_", method, ".pdf")),
          width = 8.5, height = 11)
      for (ht in heatmap_plot_list) {
        grid::grid.newpage()
        grid::grid.draw(ht)
      }
      dev.off()
    }
  }
  
  bar_plots <- cowplot::plot_grid(plotlist = bar_plot_list, ncol = 3, nrow = 3)
  ggplot2::ggsave(filename = file.path(bar_output_dir, paste0("Bar_plot_pathways_", method, ".tiff")),
                  plot = bar_plots,
                  device = "jpeg",
                  width = 3*7,
                  height = 3*7,
                  units = "in",
                  dpi = 300,
                  bg = "white")
  
  dot_plots <- cowplot::plot_grid(plotlist = dot_plot_list, ncol = 3, nrow = 3)
  ggplot2::ggsave(filename = file.path(dot_output_dir, paste0("Dot_plot_pathways_", method, ".tiff")),
                  plot = dot_plots,
                  device = "jpeg",
                  width = 3*7,
                  height = 3*7,
                  units = "in",
                  dpi = 300,
                  bg = "white")
}

plot_tf <- function(input, contrast, meta_data, samples, output_dir, n_tfs = 20){
  
  set.seed(1234)
  
  # ---- Create output directories ----
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("📂 Created directory: ", output_dir)
  }
  
  # NOTE: wsum returns wsum, norm_wsum and corr_wsum.
  # wsum (DONT USE): Biased toward larger gene sets (more genes → bigger sum)
  # norm_wsum (DONT USE): Adjusts for pathway length so small and large gene sets are comparable
  # corr_sum (DONT USE): corrects for high correlation as it can make enrichment appear stronger
  
  target <- stringr::str_split(string = contrast, pattern = "-")[[1]][1]
  ref <- stringr::str_split(string = contrast, pattern = "-")[[1]][2]
  stats <- c("consensus", "ulm", "mlm", "viper") #unique(input$all_tfs$statistic)
  bar_plot_list <- list()
  heatmap_plot_list <- list()
  
  for (stat in stats) {
    
    if (length(unique(input$all_tfs[["condition"]])) == 1){
      top_tf <- input$all_tfs %>%
        dplyr::mutate(Direction = dplyr::case_when(score < 0 ~ "Downregulated",
                                                   score > 0 ~ "Upregulated",
                                                   TRUE ~ "No change")) %>%
        dplyr::group_by(statistic, Direction) %>%
        dplyr::slice_max(order_by = abs(score), n = n_tfs, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::filter(statistic == stat)
      
      p <- ggplot(data = top_tf, aes(x = reorder(source, score), y = score, fill = score)) +
        geom_col(width = 0.75, na.rm = TRUE) +
        scale_fill_gradient2(low = "darkblue", high = "indianred", mid = "whitesmoke", midpoint = 0) +
        labs(x = "", y = "Score", title = paste0("Top TFs (", stat, " method)"), fill = "Score") +
        custom_theme +
        coord_cartesian(clip = "off") +
        geom_text(label = paste0("Activated in ", target),
                  x = top_tf$source[which.max(top_tf$score)],
                  y = ceiling(max(top_tf$score)) + 1, hjust = 1, color = "indianred", fontface = "bold") +
        geom_text(label = paste0("Activated in ", ref),
                  x = top_tf$source[which.min(top_tf$score)],
                  y = ceiling(max(top_tf$score)) + 1, hjust = 0, color = "darkblue", fontface = "bold")
      
      bar_plot_list <- c(bar_plot_list, list(p))
    }
    
    else{
      top_tf <- input$all_tfs %>%
        dplyr::filter(statistic == stat) %>%
        dplyr::group_by(source) %>%
        dplyr::summarise(std = sd(score, na.rm = TRUE), .groups = "drop") %>%
        dplyr::slice_max(order_by = abs(std), n = n_tfs, with_ties = FALSE) %>%
        dplyr::pull(source)
      
      tf_mat <- input$all_tfs %>%
        dplyr::filter(statistic == stat) %>%
        tidyr::pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
        tibble::column_to_rownames("condition") %>%
        dplyr::select(all_of(top_tf)) %>%
        as.matrix() %>%
        t()
      
      norm_counts <- tf_mat[top_tf, samples, drop = FALSE]
      metadata_col <- meta_data %>% dplyr::filter(Sample_ID %in% samples)
      metadata_row <- NULL
      disp_genes <- c()
      proj.params$heatmap$title <- paste0("Top TFs (", stat, ") method")
      
      ph <- plot_heatmap(norm_counts, proj.params, metadata_col, metadata_row, disp_genes)
      heatmap_plot_list[[length(heatmap_plot_list) + 1]] <- ph$ph$gtable
    }
  }
  
  # Save combined heatmaps
  if (length(heatmap_plot_list) > 0) {
    pdf(file.path(output_dir, "Heatmap_TFs.pdf"),
        width = 10, height = 8)
    for (ht in heatmap_plot_list) {
      grid::grid.newpage()
      grid::grid.draw(ht)
    }
    dev.off()
  }
  
  # Save combined barplots
  if (length(bar_plot_list) > 0) {
    bar_plots <- cowplot::plot_grid(plotlist = bar_plot_list, ncol = 1, nrow = length(bar_plot_list))
    ggplot2::ggsave(filename = file.path(output_dir, "Bar_plot_TFs.tiff"),
                    plot = bar_plots,
                    device = "jpeg",
                    width = 10,
                    height = 4 * 3,
                    units = "in",
                    dpi = 300,
                    bg = "white")
  }
  
  return(invisible(list(barplots = bar_plot_list, heatmaps = heatmap_plot_list)))
}

main_analysis <- function(meta_data, read_data, proj.params, trial) {
  
  # Compile raw counts if read_data not available
  if (is.null(read_data)) {
    if (dir.exists(proj.params$counts_dir)) {
      read_data <- merge_counts(proj.params)
    } else {
      stop("The specified count directory does not exist: ", proj.params$counts_dir)
    }
  }
  
  # Prepare input
  deseq2_inputs <- prepare_deseq2_input(meta_data, read_data, proj.params)
  contrasts <- proj.params$deseq2$contrast
  
  # ---- Visualization: PCA plot ----
  meta_data <- deseq2_inputs$meta_data
  read_data <- deseq2_inputs$read_data
  output_dir <- proj.params$proj_dir
  plot_pca(meta_data, read_data, output_dir)
  
  if(!trial){
    
    # Loop through contrasts
    for (n in seq_along(contrasts)) {
      
      # ---- Differential Expression (DESeq2) ----
      meta_data <- deseq2_inputs$meta_data
      read_data <- deseq2_inputs$read_data
      deseq2_results <- run_deseq2(meta_data, read_data, proj.params, n)
      
      # ---- Visualization: MA Plot ----
      dds <- deseq2_results$dds
      output_dir <- proj.params$deseq2_dir[n]
      plot_ma(dds, output_dir)
      
      # ---- Visualization: Volcano Plot ----
      DEGs_df <- deseq2_results$degs
      contrast <- proj.params$deseq2$contrasts[n]
      output_dir <- proj.params$deseq2_dir[n]
      plot_volcano(DEGs_df, proj.params, contrast, output_dir)
      
      
      # ---- Visualization: Heatmap ----
      contrast <- proj.params$deseq2$contrasts[n]
      DEGs_df <- deseq2_results$degs
      vst_counts <- deseq2_results$vst
      samples <- filter_samples_by_contrast(meta_data, contrast)
      samples <- intersect(colnames(vst_counts), samples)
      sig_genes <- DEGs_df %>% dplyr::filter(padj <= 0.05) %>% dplyr::pull(SYMBOL)
      
      if (length(sig_genes) != 0){
        
        norm_counts <- vst_counts[sig_genes, samples, drop = FALSE]
        metadata_col <- meta_data %>% dplyr::filter(Sample_ID %in% samples)
        metadata_row <- NULL
        disp_genes <- c()
        output_dir <- proj.params$deseq2_dir[n]
        
        ph <- plot_heatmap(norm_counts, proj.params, metadata_col, metadata_row, disp_genes)
        
        jpeg(file.path(output_dir, "Heatmap.jpeg"), width = 5, height = 7, units = "in", res = 300)
        gridExtra::grid.arrange(grobs = list(ph$ph$gtable), ncol = 1)
        dev.off()
        
        save_xlsx(ph$mat, file.path(output_dir, "Heatmap_Matrix.xlsx"), "Heatmap_matrix", row_names = TRUE)
      } 
      
      # ---- Pathway Analysis ----
      
      # Perform pathway analysis
      output_dir <- proj.params$pathway_dir[n]
      gsea_res_list <- pathway_analysis(DEGs_df, proj.params, output_dir)
      
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
      
      output_dir <- proj.params$pathway_dir[n]
      samples <- filter_samples_by_contrast(meta_data, contrast)
      samples <- intersect(colnames(vst_counts), samples)
      plot_pathway(top_pathways_gsea, vst_counts, meta_data, samples, "GSEA", output_dir)
      plot_pathway(top_pathways_ora, vst_counts, meta_data, samples, "ORA", output_dir)
      
      # ---- Transcription Factor (TF) Analysis ----
      
      # Perform TF analysis on DEGs using t-statistics
      t_stats_mat <- DEGs_df %>% 
        as.data.frame() %>%
        dplyr::mutate(t = -log10(padj) * log2FoldChange) %>%
        dplyr::filter(!is.na(t)) %>%
        dplyr::select(SYMBOL, t) %>%
        tibble::column_to_rownames("SYMBOL") %>%
        as.matrix()
      
      tf_res_degs <- tf_analysis(t_stats_mat, proj.params$species)
      
      # Perform TF analysis on vst counts 
      samples <- filter_samples_by_contrast(meta_data, contrast)
      samples <- intersect(colnames(vst_counts), samples)
      norm_counts_sub <- vst_counts[, samples]
      
      tf_res_counts <- tf_analysis(norm_counts_sub, proj.params$species)
      
      # Plot Bar, Heatmap for top 20 Up and Down TFs
      contrast <- contrasts[n]
      output_dir <- proj.params$tf_dir[n]
      samples <- filter_samples_by_contrast(meta_data, contrast)
      samples <- intersect(colnames(vst_counts), samples)
      plot_tf(tf_res_degs, contrast, meta_data, samples, output_dir, n_tfs = 20)
      plot_tf(tf_res_counts, contrast, meta_data, samples, output_dir, n_tfs = 20)
    }
    
    # ---- Visualization: Dispersion Estimates ----
    dds <- deseq2_results$dds
    output_dir <- proj.params$proj_dir
    plot_dispersion(dds, output_dir)
    
    # ---- VST Counts (Blind) ----
    # NOTE: vst counts are affected by design ONLY when blind=FALSE
    dds <- deseq2_results$dds
    output_dir <- proj.params$proj_dir
    vsd_blind <- DESeq2::vst(dds, blind = TRUE)
    vst_counts_blind <- SummarizedExperiment::assay(vsd_blind) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("ID") %>%
      add_annotation() %>%
      dplyr::filter(!is.na(SYMBOL)) %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::summarize(across(.cols = where(is.numeric), .fns = mean, na.rm = TRUE), .groups = "drop") %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("SYMBOL") %>%
      dplyr::select(-dplyr::starts_with("ENSEMBL"), -dplyr::starts_with("ENTREZ")) %>%
      as.matrix()
    
    save_xlsx(vst_counts_blind, file.path(output_dir, "VST_counts_blind.xlsx"), "VST_blind", row_names = TRUE)
    
    # ---- Perform LRT Test ----
    dds <- deseq2_results$dds
    dds_LRT <- DESeq2::DESeq(dds, test = "LRT", reduced = ~1)
    res_LRT <- DESeq2::results(dds_LRT) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("ID") %>%
      add_annotation()
    
    # ---- Visualization: Heatmap ----
    output_dir <- proj.params$proj_dir
    samples <- colnames(vst_counts_blind)
    sig_genes <- res_LRT %>% dplyr::filter(padj <= 0.05) %>% dplyr::pull(SYMBOL)
    norm_counts <- vst_counts_blind[sig_genes, samples, drop = FALSE]
    metadata_col <- meta_data %>% dplyr::filter(Sample_ID %in% samples)
    metadata_row <- NULL
    disp_genes <- c()
    proj.params$heatmap.title <- ""
    if (length(sig_genes) != 0){
      
      ph <- plot_heatmap(norm_counts, proj.params, metadata_col, metadata_row, disp_genes)
      
      jpeg(file.path(output_dir, "Heatmap.jpeg"), width = 5, height = 7, units = "in", res = 300)
      gridExtra::grid.arrange(grobs = list(ph$ph$gtable), ncol = 1)
      dev.off()
      
      save_xlsx(ph$mat, file.path(output_dir, "Heatmap_Matrix.xlsx"), "Heatmap_matrix", row_names = TRUE)
    }
  }
}

get_annotations <- function() {
  
  set.seed(1234)
  
  # Initialize 
  species_list <- c("Homo sapiens", "Mus musculus")
  annotations_list <- list()
  
  for (species in species_list) {
    
    # Connect to AnnotationHub and Fetch Ensembl DB 
    ah <- AnnotationHub::AnnotationHub()
    ah_db <- AnnotationHub::query(x = ah, 
                                  pattern = c(species, "EnsDb"), 
                                  ignore.case = TRUE)
    
    # Acquire the latest annotation files
    latest <- ah_db %>%
      mcols() %>%
      rownames() %>%
      tail(n=1)
    
    # Download the appropriate Ensembldb database
    edb <- ah[[latest]]
    
    # Extract ENSEMBL Annotations 
    ensembl <- ensembldb::genes(x = edb, 
                                return.type = "data.frame") %>%
      dplyr::rename(ENSEMBL_ID    = gene_id,
                    ENSEMBL_SYMBOL  = gene_name,
                    ENSEMBL_BIOTYPE  = gene_biotype,
                    START       = gene_seq_start,
                    END        = gene_seq_end,
                    CHR        = seq_name,
                    STRAND      = seq_strand,
                    DESCRIPTION    = description,
                    ENSEMBL_TRANSCRIPT = canonical_transcript) %>%
      dplyr::mutate(ENSEMBL_SYMBOL = dplyr::case_when(nchar(ENSEMBL_SYMBOL) == 0 ~ NA,
                                                      TRUE ~ ENSEMBL_SYMBOL)) %>%
      dplyr::select(ENSEMBL_ID, ENSEMBL_TRANSCRIPT, ENSEMBL_SYMBOL, ENSEMBL_BIOTYPE,
                    START, END, CHR, STRAND, DESCRIPTION)
    
    # Extract ENTREZ Annotations 
    org_db <- if (species == "Homo sapiens") org.Hs.eg.db else org.Mm.eg.db
    entrez <- AnnotationDbi::select(x = org_db,
                                    keys = AnnotationDbi::keys(org_db),
                                    columns = c("ENSEMBL", "SYMBOL", "GENETYPE")) %>%
      dplyr::rename(ENTREZ_ID    = ENTREZID,
                    ENSEMBL_ID   = ENSEMBL,
                    ENTREZ_SYMBOL  = SYMBOL,
                    ENTREZ_BIOTYPE = GENETYPE)
    
    # Merge Annotations 
    annotations <- dplyr::full_join(ensembl, entrez, by=c("ENSEMBL_ID"="ENSEMBL_ID")) %>%
      dplyr::select(ENSEMBL_ID, ENSEMBL_TRANSCRIPT, ENTREZ_ID, ENSEMBL_SYMBOL, 
                    ENTREZ_SYMBOL, ENSEMBL_BIOTYPE, ENTREZ_BIOTYPE,START, END,
                    CHR, STRAND, DESCRIPTION)
    
    # Store Output 
    annotations_list[[species]] <- annotations
  }
  
  # Return: list(human, mouse) 
  return(invisible(annotations_list))
}

add_annotation <- function(normalized_counts) {
  
  set.seed(1234)
  
  # Retrieve Annotations 
  annotations <- get_annotations()
  
  # Flatten all annotation data into one named list 
  named_lists <- list(ensembl_id_human   = annotations$`Homo sapiens`$ENSEMBL_ID,
                      entrez_id_human   = annotations$`Homo sapiens`$ENTREZ_ID,
                      ensembl_symbol_human = annotations$`Homo sapiens`$ENSEMBL_SYMBOL,
                      entrez_symbol_human = annotations$`Homo sapiens`$ENTREZ_SYMBOL,
                      ensembl_id_mouse   = annotations$`Mus musculus`$ENSEMBL_ID,
                      entrez_id_mouse   = annotations$`Mus musculus`$ENTREZ_ID,
                      ensembl_symbol_mouse = annotations$`Mus musculus`$ENSEMBL_SYMBOL,
                      entrez_symbol_mouse = annotations$`Mus musculus`$ENTREZ_SYMBOL)
  
  # Compute intersection counts 
  overlap_counts <- sapply(X = named_lists, 
                           FUN = function(x) {length(intersect(x, normalized_counts$ID))})
  
  # Identify the best matching column 
  best_match <- names(which.max(overlap_counts))
  message("Best match: ", best_match)
  message("Overlap counts:\n", paste(names(overlap_counts), overlap_counts, sep = " \t: ", collapse = "\n"))
  
  # Select Species-Specific Annotation 
  if (grepl(pattern = "human", x = best_match)) {
    ann <- annotations$`Homo sapiens`
  } else if (grepl(pattern = "mouse", x = best_match)) {
    ann <- annotations$`Mus musculus`
  } else {
    stop("Unable to determine organism (human/mouse) from ID column.")
  }
  
  # Select ID and SYMBOL Columns 
  if (grepl(pattern = "ensembl", x = best_match)) {
    id_col <- "ENSEMBL_ID"
    symbol_col <- "ENSEMBL_SYMBOL"
  } else if (grepl(pattern ="entrez", x = best_match)) {
    id_col <- "ENTREZ_ID"
    symbol_col <- "ENTREZ_SYMBOL"
  } else {
    stop("Unable to determine ID type (Ensembl/Entrez) from ID column.")
  }
  
  # Finalize Columns 
  keep_ids <- c(id_col, symbol_col)
  drop_ids <- setdiff(colnames(ann), keep_ids)
  
  # Join Annotations and Create SYMBOL Column 
  normalized_counts <- ann %>%
    dplyr::right_join(normalized_counts, by = stats::setNames("ID", id_col), multiple = "all") %>%
    dplyr::mutate(SYMBOL = dplyr::case_when(is.na(ENSEMBL_SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                            is.na(ENSEMBL_SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                            TRUE ~ ENSEMBL_SYMBOL)) %>%
    dplyr::select(SYMBOL, all_of(keep_ids), dplyr::everything(), -all_of(drop_ids)) %>%
    dplyr::distinct_at(id_col, .keep_all = TRUE)
  
  # Return norm counts with annotation added 
  return(normalized_counts)
}

norm_counts_DESeq2 <- function(meta_data, read_data, proj.params) {
  
  set.seed(1234)
  
  # Create DESeq2 Dataset 
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = read_data,
                                        colData = meta_data,
                                        design = ~1)
  dds <- DESeq2::estimateSizeFactors(object = dds)
  
  # Normalized Counts 
  normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
  
  # VST-Transformed Counts 
  vsd <- DESeq2::vst(dds, blind = TRUE)
  vst_counts <- SummarizedExperiment::assay(vsd)
  
  # Batch Correction (if applicable) 
  if ("Batch" %in% colnames(meta_data) && length(unique(meta_data$Batch)) > 1) {
    normalized_counts_batch <- limma::removeBatchEffect(x = log2(normalized_counts + 1),
                                                        batch = dds$Batch)
  } else {
    normalized_counts_batch <- normalized_counts
    message("No batch correction performed")
  }
  
  # Convert to Data Frames 
  normalized_counts_df <- normalized_counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    add_annotation()
  
  normalized_counts_batch_df <- normalized_counts_batch %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    add_annotation()
  
  vst_counts_df <- vst_counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    add_annotation()
  
  # Write to Excel 
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "VST_counts(blind=TRUE)")
  openxlsx::writeData(wb, sheet = "VST_counts(blind=TRUE)", x = vst_counts_df, rowNames = FALSE)
  
  openxlsx::addWorksheet(wb, sheetName = "Norm_counts")
  openxlsx::writeData(wb, sheet = "Norm_counts", x = normalized_counts_df, rowNames = FALSE)
  
  openxlsx::addWorksheet(wb, sheetName = "Norm_counts_batch_corrected")
  openxlsx::writeData(wb, sheet = "Norm_counts_batch_corrected", x = normalized_counts_batch_df, rowNames = FALSE)
  
  openxlsx::saveWorkbook(wb, file = file.path(proj.params$dir, proj.params$proj, "Normalized.counts.DESeq2.xlsx"), overwrite = TRUE)
  
  
  return(invisible(NULL))
}

add_major_pathway <- function(df){
  
  set.seed(1234)
  
  df <- df %>% 
    dplyr::mutate(MAJOR_PATHWAY = dplyr::case_when(
      grepl(pattern = "hallmark", x = base::tolower(Description)) ~ "HALLMARK",
      grepl(pattern = "interferon|interleukin|cytokine|chemokine|immune|toll|antigen|leukocyte|lymphocyte|macrophage", x = base::tolower(Description)) ~ "IMMUNE RELATED",
      grepl(pattern = "metaboli|purine|carbohydrate", x = base::tolower(Description)) ~ "METABOLISM",
      grepl(pattern = "translation", x = base::tolower(Description)) ~ "PROTEIN REGULATION",
      grepl(pattern = "transcription", x = base::tolower(Description)) ~ "GENE REGULATION",
      grepl(pattern = "mitotic|cell cycle", x = base::tolower(Description)) ~ "CELL CYCLE",
      grepl(pattern = "muscle", x = base::tolower(Description)) ~ "MUSCLE",
      grepl(pattern = "cardiac", x = base::tolower(Description)) ~ "HEART",
      grepl(pattern = "angiogenesis|blood_vessel", x = base::tolower(Description)) ~ "ANGIOGENESIS",
      grepl(pattern = "actin_", x = base::tolower(Description)) ~ "CYTOSKELETAN",
      grepl(pattern = "glycosyl", x = base::tolower(Description)) ~ "GLYCOSYLATION",
      grepl(pattern = "dna_", x = base::tolower(Description)) ~ "DNA DAMAGE/REPAIR",
      grepl(pattern = "rna_", x = base::tolower(Description)) ~ "RNA REGULATION",
      grepl(pattern = "ubiquitin|proteasome", x = base::tolower(Description)) ~ "PROTEIN DEGRATION",
      grepl(pattern = "transport", x = base::tolower(Description)) ~ "LOCALIZATION",
      grepl(pattern = "phagy|apopto", x = base::tolower(Description)) ~ "CELL DEATH",
      grepl(pattern = "ribosom", x = base::tolower(Description)) ~ "RIBOSOME",
      grepl(pattern = "gtpase", x = base::tolower(Description)) ~ "GTPASE",
      TRUE ~ "UNCLASSIFIED")) %>%
    dplyr::select(MAJOR_PATHWAY, everything())
  
  return(invisible(df))
}

### PROGENy analysis [NOT RECOMMENDED as decoupleR is better]
progeny_analysis <- function(norm_counts, assay = "RNA", species = "Homo sapiens") {
  
  set.seed(1234)
  
  # --- Error checking ---
  if (!is.matrix(norm_counts)) stop("`norm_counts` must be a numeric matrix (genes x samples).")
  if (!species %in% c("Homo sapiens", "Mus musculus")) stop("`species` must be 'Homo sapiens' or 'Mus musculus'.")
  
  organism <- dplyr::if_else(species == "Homo sapiens", "Human", "Mouse")
  
  # --- Raw PROGENy scores (no permutation) ---
  progeny_scores <- progeny(expr = norm_counts,
                            scale = FALSE,  # Refer section below to understand how it works
                            organism = organism,
                            top = 500,
                            perm = 1,
                            verbose = FALSE,
                            z_scores = FALSE,
                            get_nulldist = FALSE,
                            assay_name = assay,
                            return_assay = FALSE)
  
  # --- Significance scores from permutations ---
  progeny_sig_scores <- progeny(expr = norm_counts,
                                scale = FALSE,  # Refer section below to understand how it works
                                organism = organism,
                                top = 100,
                                perm = 10000,  # Returns list where [[1]] has p values
                                verbose = FALSE,
                                z_scores = FALSE,
                                get_nulldist = FALSE, # TRUE swaps rows and columns
                                assay_name = assay,
                                return_assay = FALSE)
  
  # --- Empirical two-sided p-values ---
  progeny_pvals <- (1-abs(progeny_sig_scores))/2
  
  # --- Adjusted p-values (FDR) ---
  pval_vec <- as.vector(progeny_pvals)
  padj_vec <- stats::p.adjust(pval_vec, method = "BH")
  
  progeny_padj <- matrix(data = padj_vec, 
                         nrow = nrow(progeny_pvals), 
                         ncol = ncol(progeny_pvals),
                         dimnames = dimnames(progeny_pvals))
  
  # --- Return as named list ---
  progeny_list <- list(scores = progeny_scores,
                       padj = progeny_padj,
                       significance = progeny_sig_scores,
                       pval = progeny_pvals)
  
  return(invisible(progeny_list))
  
  # NOTE: Below section is for understanding how scale parameter works
  if (FALSE) {
    
    # Sample expr data
    human_input <- as.matrix(read.csv(system.file("extdata", "human_input.csv", 
                                                  package = "progeny"), 
                                      row.names = 1))
    mouse_input <- as.matrix(read.csv(system.file("extdata", "mouse_input.csv", 
                                                  package = "progeny"), 
                                      row.names = 1))
    
    # Expected pathway scores with scale=TRUE and top = 10
    human_def_expected <- read.csv(system.file("extdata", "human_def_expected.csv", 
                                               package = "progeny"),
                                   row.names = 1,
                                   check.names = FALSE)
    mouse_def_expected <- read.csv(system.file("extdata", "mouse_def_expected.csv", 
                                               package = "progeny"),
                                   row.names = 1, 
                                   check.names = FALSE)
    
    # With scaling
    scores_scaled <- progeny(human_input, scale = TRUE, organism = "Human", top = 10)
    scores_scaled_subset <- progeny(human_input[,1:5], scale = TRUE, organism = "Human", top = 10)
    
    # Without scaling
    scores_raw <- progeny(human_input, scale = FALSE, organism = "Human", top = 10)
    scores_raw_subset <- progeny(human_input[,1:5], scale = FALSE, organism = "Human", top = 10)
    scores_raw_perm <- progeny(human_input, scale = FALSE, organism = "Human", top = 10, perm=10)
    
    # Check mean and sd
    apply(X=scores_scaled, MARGIN=2, FUN=function(x){mean(x) %>% round(.,2)})
    apply(X=scores_scaled, MARGIN=2, FUN=function(x){sd(x) %>% round(.,2)})
    apply(X=scores_raw, MARGIN=2, FUN=function(x){mean(x) %>% round(.,2)})
    apply(X=scores_raw, MARGIN=2, FUN=function(x){sd(x) %>% round(.,2)})
    
    # Scale the raw scores
    scores_raw_scaled <- scale(x=scores_raw, center=TRUE, scale=TRUE)
    
    # Verify
    glimpse(scores_raw_scaled)
    glimpse(scores_scaled)
    identical(scores_raw_scaled %>% as.data.frame(), scores_scaled %>% as.data.frame())
    
    # [1] This shows scale parameter ONLY scales the calculated pathway scores.
    # and doesnt scale the input gene expression data.
    # [2] raw scores DONT CHANGE based on number of samples but scaled scores
    # vary if number of samples change.
    # RECOMMENDATION: Set scale=FALSE and scale the scores later when plotting
    
  }
}

# ---- 🧬 SINGLE CELL & SPATIAL ANALYSIS RELATED FUNCTIONS ----

validate_inputs <- function(step, sample = NULL, gene.column = NULL, bin = NULL, 
                            assay = NULL, n_pcs = NULL, pN = NULL,
                            reference_samples = NULL,
                            features = NULL, s_genes = NULL, g2m_genes = NULL, 
                            resolution = NULL, reduction = NULL, filename = NULL,
                            seurat_object = NULL, seurat_list = NULL, 
                            raw_metadata = NULL,  
                            matrix_dir = NULL, output_dir = NULL,
                            metafile = NULL, metadata_cols = NULL) {
  
  caller <- as.character(sys.call(-1)[[1]])  # returns the calling function object as a string
  step <- paste0(caller, " : validate_inputs")
  
  # ---- 1️⃣ Basic Parameter Checks ----
  
  if (!is.null(sample) && (!is.character(sample) || length(sample) != 1 || nchar(sample) == 0)) {
    log_error(sample, step, "❌ INPUT ERROR: 'sample' must be a single, non-empty character string.")
  }
  
  if (!is.null(gene.column) && (!is.numeric(gene.column) || length(gene.column) != 1 || !(gene.column %in% c(1,2)))) {
    log_error(sample, step, "❌ INPUT ERROR: 'gene.column' must be 1 (Ensembl IDs) or 2 (Gene symbols).")
  }
  
  if (!is.null(bin) && (!is.numeric(bin) || length(bin) != 1 || !(bin %in% c(2,8,16)))) {
    log_error(sample, step, "❌ INPUT ERROR: Invalid 'bin' size provided. Must be one of 2, 8, or 16.")
  }
  
  if (!is.null(assay) && (!is.character(assay) || length(assay) != 1 || nchar(assay) == 0)) {
    log_error(sample, step, "❌ INPUT ERROR: 'assay' must be a single, non-empty character string.")
  }
  
  if (!is.null(n_pcs) && (!is.numeric(n_pcs) || length(n_pcs) != 1 || n_pcs <= 0 || n_pcs %% 1 != 0 || n_pcs > 50)) {
    log_error(sample, step, "❌ INPUT ERROR: 'n_pcs' must be a positive integer less than 50.")
  }
  
  if (!is.null(pN) && (!is.numeric(pN) || length(pN) != 1 || pN <= 0 || pN > 1)) {
    log_error(sample, step, "❌ INPUT ERROR: 'pN' must be a numeric value between 0 and 1 (exclusive of 0).")
  }
  
  if (!is.null(reduction) && (!is.character(reduction) || length(reduction) != 1 || nchar(reduction) == 0)) {
    log_error(sample, step, "❌ INPUT ERROR: 'reduction' must be a single, non-empty character string.")
  }
  
  if (!is.null(reference_samples) && (!is.character(reference_samples) || length(reference_samples) == 0)) {
    log_error(sample, step, "❌ INPUT ERROR: 'reference_samples' must be a non-empty character vector.")
  }
  
  # ---- 2️⃣ Feature Checks ----
  
  # Helper Function for Validation
  validate_gene_set <- function(gene_set, set_name, available_features) {
    
    if (is.null(gene_set)) return(invisible(NULL))
    
    existing_features <- base::intersect(gene_set, available_features)
    missing_features <- base::setdiff(gene_set, existing_features)
    
    if (length(existing_features) == 0) {
      log_error("", step, glue::glue("❌ INPUT ERROR: None of the '{length(gene_set)}' provided features were found in the Seurat object."))
    }
  }

  # Check features, s_genes, g2m_genes
  if (!is.null(seurat_object)) {
    
    # Collect ALL available features once (genes from default assay + metadata columns)
    available_features <- c(SeuratObject::Features(seurat_object),
                            colnames(seurat_object@meta.data))
    
    # Run validation for each set using the helper function
    if (!is.null(features)) {
      features <- validate_gene_set(features, "features", available_features)
    }
    
    if (!is.null(s_genes)) {
      s_genes <- validate_gene_set(s_genes, "s_genes", available_features)
    }
    
    if (!is.null(g2m_genes)) {
      g2m_genes <- validate_gene_set(g2m_genes, "g2m_genes", available_features)
    }
  }
  
  # ---- 3️⃣ 'resolution' Check ----
  
  if (!is.null(resolution)) {
    allowed_res <- c(0.2,0.4,0.6,0.8,1.0,1.2)
    resolution_num <- as.numeric(resolution)
    if (is.na(resolution_num) || !(resolution_num %in% allowed_res)) {
      log_error(sample, step, glue::glue("❌ INPUT ERROR: Specified resolution '{resolution}' is NOT one of the **valid** resolution : '{paste(allowed_res, collapse = ", ")}'."))
    }
  }
  
  # ---- 4️⃣ 'filename' Check ----
  
  if (!is.null(filename)) {
    if (!is.character(filename) || length(filename) != 1 || nchar(filename) == 0) {
      log_error(sample, step, "❌ 'filename' must be a single, non-empty string.")
    }
    if (grepl("[<>:\"/\\\\|?*]", filename)) {
      log_error(sample, step, "❌ 'filename' contains illegal characters for file paths.")
    }
  }
  
  # ---- 5️⃣ Seurat Object and List Integrity Checks ----
  
  # Check Seurat object class (only if provided)
  if (!is.null(seurat_object)) {
    if (!inherits(seurat_object, "Seurat")) {
      log_error(sample, step, "❌ INPUT ERROR: 'seurat_object' must be a Seurat object.")
    }
    
    # Check for Non-Empty Seurat Object
    if (ncol(seurat_object) == 0) {
      log_error(sample, step, "❌ DATA ERROR: The input Seurat object is empty (0 cells).")
    }
    
    if (!is.null(assay)) {
      # First, check if the string is structurally sound (Done in Section 1)
      # Second, check if the assay name exists in the object:
      if (!(assay %in% names(seurat_object@assays))) {
        log_error(sample, step, 
                  glue::glue("❌ ASSAY ERROR: Specified assay '{assay}' is NOT one of the **valid assays** : '{paste(names(seurat_object@assays), collapse = ', ')}'."))
      }
    }
  }
  
  # Check 'seurat_list' structure and naming (only if provided)
  if (!is.null(seurat_list)) {
    if (!is.list(seurat_list) || length(seurat_list) == 0) {
      log_error(sample, step, "❌ INPUT ERROR: 'seurat_list' must be a non-empty list of Seurat objects.")
    }
    if (is.null(names(seurat_list)) || any(nchar(names(seurat_list)) == 0)) {
      log_error(sample, step, "❌ INPUT ERROR: 'seurat_list' must be a *named* list with non-empty names.")
    }
    # Check all elements are Seurat objects
    if (!all(vapply(seurat_list, inherits, logical(1), what = "Seurat"))) {
      log_error(sample, step, "❌ INPUT ERROR: All elements in 'seurat_list' must be Seurat objects.")
    }
  }
  
  # ---- 6️⃣ 'raw_metadata' Check ----
  
  if (!is.null(raw_metadata)) {
    if (!is.data.frame(raw_metadata)) {
      log_error(sample, step, "❌ INPUT ERROR: 'raw_metadata' must be a data.frame.")
    }
    if (nrow(raw_metadata) == 0) {
      log_error(sample, step, "❌ DATA ERROR: 'raw_metadata' must be a non-empty data.frame (0 rows detected).")
    }
  }
  
  # ---- 7️⃣ 'matrix_dir' Check ----
  
  if (!is.null(matrix_dir) && !is.null(sample)) {
    
    data_dir <- base::file.path(matrix_dir, sample)
    
    # Check for existence of the specific data directory
    if (!dir.exists(data_dir)) {
      log_error(sample, step, glue::glue("❌ PATH ERROR: The data directory '{data_dir}' does not exist!"))
    }
    
    # Check for required file structure (HDF5 OR 3-file structure - for CellRanger)
    h5_files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE, ignore.case = TRUE)
    expected_files_mtx <- c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")
    standard_files_found <- all(sapply(file.path(data_dir, expected_files_mtx), file.exists))
    
    if (length(h5_files) == 0 && !standard_files_found) {
      log_error(sample, step, glue::glue("❌ FILE STRUCTURE ERROR: Directory '{data_dir}' does not contain a valid 10X matrix."))
    } else if (length(h5_files) > 1) {
      log_warn(sample, step, glue::glue("⚠️ FILE STRUCTURE WARNING: Multiple HDF5 files found. Only using the first one : '{h5_files[1]}'."))
    }
    
    # # Check for REQUIRED Space Ranger files (for load_spaceranger)
    # required_spatial_files <- c("filtered_feature_bc_matrix.h5", "tissue_positions_list.csv", "tissue_lowres_image.png")
    # required_files_exist <- all(sapply(paste0(data_dir, "/", required_spatial_files), file.exists))
    # 
    # if (!required_files_exist) {
    #   stop("❌ FILE STRUCTURE ERROR: The directory '", data_dir, "' is missing one or more of the required Space Ranger files: ", paste(required_spatial_files, collapse = ", "), ".")
    # }
  }
 
  # ---- 8️⃣ 'output_dir' Check ----
  
  if (!is.null(output_dir)) {
    
    # Check 1: Directory must exist
    if (!dir.exists(output_dir)) {
      log_error(sample, step, glue::glue("❌ PATH ERROR: Output directory '{output_dir}' does not exist."))
    }
    
    # Check 2: Directory must be writable
    # file.access(path, mode=2) checks for write permission. Returns 0 on success.
    if (file.access(output_dir, 2) != 0) {
      log_error(sample, step, glue::glue("❌ PATH ERROR: Output directory '{output_dir}' is NOT writable."))
    }
  }
  
  # ---- 9️⃣ 'metafile' Check ---- 
  
  if (!is.null(metafile)) {
    
    # Check 1: File must exist
    if (!file.exists(metafile)) {
      log_error(sample, step, glue::glue("❌ FILE ERROR: External metadata file not found at :  '{metafile}'."))
    }
    
    # Check 2: File must be in xlsx format
    if (!grepl("\\.xlsx$", metafile, ignore.case = TRUE)) {
      log_error(sample, step, "❌ FILE ERROR: 'metafile' must be an .xlsx file.")
    }
  }
  
  # ---- 🔟 'metadata_cols' Check ----
  
  if (!is.null(seurat_object)) {
    
    required_cols <- unlist(metadata_cols)
    required_cols <- required_cols[!sapply(required_cols, is.null)]
    
    if (length(required_cols) > 0) {
      missing_cols <- base::setdiff(required_cols, colnames(seurat_object@meta.data))
      
      if (length(missing_cols) > 0) {
        log_error(sample, step, glue::glue("❌ INPUT ERROR: Required column(s) '{paste(missing_cols, collapse = ", ")}' NOT found in metadata.")) 
      }
    }
  }

}

load_cellranger <- function(sample, matrix_dir, gene.column = 2){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(sample = sample, matrix_dir = matrix_dir, gene.column = gene.column)
  
  # ---- 📥 Read the Feature-Barcode Matrix ----
  
  # gene.column = 1 → Ensembl IDs
  # gene.column = 2 → Gene symbols (required for mito/ribo/heme ratios)
  
  data_dir <- base::file.path(matrix_dir, sample)
  
  # Look for HDF5 files in the directory
  h5_files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE, ignore.case = TRUE)
  
  counts <- NULL
  
  # If HDF5 files are found, try to read them
  if (length(h5_files) > 0) {
    log_info(sample = sample, 
            step = "load_cellranger",
            message = glue::glue("Reading HDF5 matrix from '{h5_files[1]}'."))
    
    counts <- tryCatch({
      # Read the gene expression matrix from the HDF5 file
      Seurat::Read10X_h5(filename = h5_files[1],
                         use.names = TRUE,
                         unique.features = TRUE)
    }, error = function(e) {
      log_error(sample = sample, 
                step = "load_cellranger",
                message = glue::glue("Failed to read HDF5 matrix from '{h5_files[1]}'."))
    })
    
    # Error check: Ensure the read returns a matrix, not a list
    if (is.list(counts)) {
      log_error(sample = sample, 
              step = "load_cellranger",
              message = "HDF5 file returned a list. Check the file content.")
    }
    
  } else {
    # If no HDF5 files, fall back to the standard 10X directory
    log_info(sample = sample, 
             step = "load_cellranger",
             message = "No HDF5 file found. Reading standard 10X matrix from directory.")
    
    counts <- tryCatch({
      Seurat::Read10X(data.dir = data_dir,
                      gene.column = gene.column, # Parameterized (Default: 2 for gene symbols)
                      cell.column = 1,
                      unique.features = TRUE,
                      strip.suffix = FALSE
      )
    }, error = function(e) {
      log_error(sample = sample, 
                step = "load_cellranger",
                message = glue::glue("Failed to read 10X matrix from '{data_dir}'."))
    })
  }
  
  # ---- 🏗️ Create Seurat Object ----
  
  sample_seurat <- SeuratObject::CreateSeuratObject(counts = counts,
                                                    project = sample,
                                                    assay = "RNA",
                                                    names.field = 1,
                                                    names.delim = "_",
                                                    meta.data = NULL,
                                                    min.cells = 0,
                                                    min.features = 0)
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = sample, 
           step = "load_cellranger",
           message = glue::glue("Successfully created Seurat object for sample : '{sample}'."))
  return(invisible(sample_seurat))
}

load_spaceranger <- function(sample, bin, matrix_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(sample = sample, bin = bin, matrix_dir = matrix_dir)
  
  # ---- 📥 Load Seurat Object with Filtered Matrix ----
  
  data_dir <- base::file.path(matrix_dir, sample)
  
  # Load10X_Spatial automatically finds the necessary files (images, coordinates, matrix)
  sample_seurat <- tryCatch({
    Seurat::Load10X_Spatial(data.dir = data_dir,
                            filename = "filtered_feature_bc_matrix.h5",
                            assay = "Spatial",
                            slice = sample,
                            bin.size = bin,
                            filter.matrix = TRUE,
                            to.upper = FALSE,
                            image = NULL)
  }, error = function(e) {
    pad <- base::strrep(" ", nchar(sample) + 30)
    log_error(sample = sample, 
              step = "load_spaceranger",
              message = glue::glue("Failed to load binned data from '{data_dir}' for bin size '{bin}'.\n",
                                   "{pad}Please verify the Space Ranger output structure."))
  })
  
  # Add sample identifier to Seurat object
  sample_seurat[["orig.ident"]] <- sample
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = sample, 
           step = "load_spaceranger",
           message = glue::glue("Successfully created Seurat object for sample : '{sample}' ('{bin}' bin)."))
  return(invisible(sample_seurat))
}

classify_dropletutils <- function(sample_seurat){ 
  
  # Set seed for reproducible stochastic processes (emptyDrops)
  set.seed(100)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = sample_seurat)
  
  sample <- as.character(unique(sample_seurat@meta.data$orig.ident))
  
  # ---- 🧪 Iterative emptyDrops Run ----
  
  # Convert Seurat Object to SingleCellExperiment (MUST contains all barcodes)
  sce <- Seurat::as.SingleCellExperiment(x = sample_seurat)
  
  # Initialize Parameters for Iterative emptyDrops
  niters <- 10000
  
  # Wrapper to safely run emptyDrops
  df <- tryCatch({
    
    # Run emptyDrops to classify droplets as empty or non-empty
    e.out <- DropletUtils::emptyDrops(m = SingleCellExperiment::counts(sce),
                                      niters = niters)
    
    # Convert emptyDrops Output to Data Frame for easier filtering
    df <- as.data.frame(e.out)
    
    # Loop to handle limited droplets if needed
    # NOTE: If FDR > 0.05 AND Limited == TRUE, the FDR of these droplets can be 
    # reduced with more iterations.
    n_improve <- nrow(df %>% filter(Limited == TRUE, FDR > 0.05))
    while (n_improve > 0) {
      niters <- niters + 10000
      e.out <- DropletUtils::emptyDrops(m = SingleCellExperiment::counts(sce),
                                        niters = niters)
      df <- as.data.frame(e.out)
      n_improve <- nrow(df %>% filter(Limited == TRUE, FDR > 0.05))
      pad <- base::strrep(" ", nchar(sample) + 35)
      log_info(sample = sample, 
               step = "classify_dropletutils",
               message = glue::glue("emptyDrops check: '{n_improve}' droplets need more iterations.\n",
                                    "{pad}Current niters = '{niters}'."))
    }
    
    df
  }, error = function(e) {
    
    # Catch the specific ambient profile error
    if (grepl("no counts available to estimate the ambient profile", e$message)) {
      log_warn(sample = sample, 
               step = "classify_dropletutils",
               message = paste("emptyDrops failed: no counts available for ambient profile.\n",
                               base::strrep(" ", 35), "Setting all barcodes as non-empty (FDR < 0.05)."))
      
      # Create dummy dataframe with FDR < 0.05
      data.frame(Total = colSums(SingleCellExperiment::counts(sce)),
                 LogProb = NA,
                 PValue = NA,
                 FDR = 0,          # Mark all as significant / non-empty
                 Limited = FALSE,
                 row.names = colnames(sce))
    } else {
      # Re-throw other errors
      log_error(sample = sample, 
               step = "classify_dropletutils",
               message = glue::glue("Caught error ({class(e)[1]}): {e$message}")) 
    }
  })
  
  # Identify true (non-empty) cells with FDR ≤ 0.05
  # NOTE: Filtering out NA ensures only tested barcodes are considered.
  true_barcodes <- df %>%
    dplyr::filter(FDR <= 0.05, !is.na(FDR)) %>% 
    rownames()
  
  # ---- 🔍 Build Classification Vector ----
  
  # All barcodes present in the Seurat object — includes empty droplets
  all_barcodes <- colnames(sample_seurat)
  
  # Vectorized classification
  classification_vector <- dplyr::case_when(all_barcodes %in% true_barcodes ~ "Non-Empty Droplet",
                                            TRUE ~ "Empty Droplet")
  
  # Assign names so Seurat maps the metadata correctly by barcode
  names(classification_vector) <- all_barcodes
  
  # ---- 📥 Add Metadata to Seurat Object ----
  
  # Add the classification vector to the original Seurat object
  sample_seurat <- SeuratObject::AddMetaData(object   = sample_seurat,
                                             metadata = classification_vector,
                                             col.name = "DropletUtils")
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = sample, 
           step = "classify_dropletutils",
           message =  glue::glue("DropletUtils classification complete for sample : '{sample}'."))
  return(invisible(sample_seurat))
}

classify_cellranger <- function(sample_seurat, filt_matrix_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = sample_seurat, matrix_dir = filt_matrix_dir)
  
  # ---- 📥 Read Filtered Matrix for Barcodes ----
  
  # Get sample name from sample_seurat
  sample <- as.character(unique(sample_seurat@meta.data$orig.ident))
  
  # Read the filtered matrix for the sample
  sample_seurat_filt <- load_cellranger(sample, matrix_dir = filt_matrix_dir)
  
  # ---- 🔍 Build Classification Vector ----
  
  # All barcodes present in the Seurat object — includes empty droplets
  all_barcodes <- colnames(sample_seurat)
  
  # Barcodes identified as 'cells' by Cell Ranger's internal algorithm
  filtered_barcodes <- colnames(sample_seurat_filt)
  
  # Vectorized classification
  classification_vector <- dplyr::case_when(all_barcodes %in% filtered_barcodes ~ "Non-Empty Droplet",
                                            TRUE ~ "Empty Droplet")
  
  # Assign names so Seurat maps the metadata correctly by barcode
  names(classification_vector) <- all_barcodes
  
  # ---- 📥 Add Metadata to Seurat Object ----
  
  # Add the classification vector to the original Seurat object
  sample_seurat <- SeuratObject::AddMetaData(object   = sample_seurat,
                                             metadata = classification_vector,
                                             col.name = "CellRanger")

  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = sample, 
           step = "classify_cellranger",
           message =  glue::glue("CellRanger empty droplets identified for sample : '{sample}'."))
  return(invisible(sample_seurat))
}

calc_qc_metrics <- function(sample_seurat, assay){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = sample_seurat, assay = assay)
  
  sample <- as.character(unique(sample_seurat@meta.data$orig.ident))
  
  # ---- 📥 Populate missing optional metadata columns ----
  
  # NOTE: Use "ND" instead of NA for safer subsetting
  optional_cols <- c("DropletUtils", "CellRanger")
  for (col in optional_cols) {
    if (!col %in% colnames(sample_seurat@meta.data)) {
      sample_seurat <- Seurat::AddMetaData(object   = sample_seurat,
                                           metadata = "ND",
                                           col.name = col)
    }
  }
  
  # ---- 🔍 Calculate Percentage Feature Sets ----
  
  # Compute mitochondrial percent
  sample_seurat <- Seurat::PercentageFeatureSet(object = sample_seurat,
                                                pattern = "^[Mm][Tt]-",
                                                col.name = "MitoPercent",
                                                assay = assay)
  
  # Compute ribosomal percent
  sample_seurat <- Seurat::PercentageFeatureSet(object = sample_seurat,
                                                pattern = "^[Rr][Pp][SsLl]",
                                                col.name = "RiboPercent",
                                                assay = assay)
  
  # Compute hemoglobin percent
  sample_seurat <- Seurat::PercentageFeatureSet(object = sample_seurat,
                                                pattern = "^[Hh][Bb][AaBb]-",
                                                col.name = "HemePercent",
                                                assay = assay)
  
  # ---- 📥 Pull Metadata and Compute QC Metrics ----
  
  # (i)   Cell       : unique identifiers corresponding to each cell i.e. barcodes
  # (ii)  Sample     : sample names
  # (iii) nUMIs      : number of transcripts per cell
  # (iv)  nGenes     : number of genes per cell
  # (v)   nHTO_UMIs  : number of HTO reads per cell
  # (vi)  nHTOs      : number of HTO types per cell
  # (vii) MitoRatio   : MitoPercent / 100
  # (viii) RiboRatio  : RiboPercent / 100
  # (ix)  HemeRatio   : HemePercent / 100
  # (x)   Novelty     : log ratio of genes per UMI
  
  sample_metadata <- sample_seurat@meta.data %>%
    tibble::rownames_to_column(var = "barcode") %>%
    dplyr::mutate(Cell = paste0(orig.ident, "_", barcode),
                  Sample = orig.ident,
                  nUMIs = .data[[paste0("nCount_", assay)]],
                  nGenes = .data[[paste0("nFeature_", assay)]],
                  MitoRatio = MitoPercent / 100,
                  RiboRatio = RiboPercent / 100,
                  HemeRatio = HemePercent / 100,
                  Novelty = log10(nGenes) / log10(nUMIs))
  
  # ---- 📥 Handle HTO Metadata (if present) ----
  
  if ("nCount_HTO" %in% colnames(sample_metadata)){
    sample_metadata <- sample_metadata %>%
      dplyr::rename(nHTO_UMIs = nCount_HTO,
                    nHTOs = nFeature_HTO)
  }
  
  # ---- 📥 Add Spatial Coordinates (if available) ----
  
  if (length(sample_seurat@images) > 0){
    
    # Collect spatial coordinates into a single data frame
    df_coords <- data.frame()
    
    for (image_name in names(sample_seurat@images)){
      df <- data.frame(barcode = sample_seurat@images[[image_name]]@boundaries$centroids@cells,
                       X = sample_seurat@images[[image_name]]@boundaries$centroids@coords[,1],
                       Y = sample_seurat@images[[image_name]]@boundaries$centroids@coords[,2])
      df_coords <- dplyr::bind_rows(df_coords, df)
    }
    
    # Join coordinates based on barcode
    sample_metadata <- sample_metadata %>%
      dplyr::left_join(df_coords, by=c("barcode"="barcode"))
  }
  
  # ---- 🔍 Define QC Cutoffs and Classify Cells ----
  
  is_spatial <- length(names(sample_seurat@images)) > 0
  
  # Define cutoffs based on spatial or single-cell conditions
  if (is_spatial) {
    gene_cutoff     <- 25
    umi_cutoff      <- 50
    mito_cutoff     <- 0.2
    novelty_cutoff  <- 0.8
    ribo_cutoff     <- 0.05
  } else {
    gene_cutoff     <- 250
    umi_cutoff      <- 500
    mito_cutoff     <- 0.2
    novelty_cutoff  <- 0.8
    ribo_cutoff     <- 0.05
  }
  
  # Log the applied cutoffs
  pad <- base::strrep(" ", nchar(sample) + 46)
  log_info(sample = sample, 
           step = "calc_qc_metrics",
           message =  glue::glue("Cutoffs applied: nGenes ≥ {gene_cutoff}\n",
                                 "{pad}nUMIs ≥ {umi_cutoff}\n",
                                 "{pad}MitoRatio ≤ {mito_cutoff}\n",
                                 "{pad}Novelty ≥ {novelty_cutoff}"))
  
  # Create Quality column, drop 'barcode' column and restore it as rownames
  sample_metadata <- sample_metadata %>%
    dplyr::mutate(passes_qc = nGenes  >= gene_cutoff & 
                    nUMIs   >= umi_cutoff &
                    MitoRatio <= mito_cutoff &
                    Novelty   >= novelty_cutoff,
                  Quality = dplyr::case_when(
                    DropletUtils == "Empty Droplet" & CellRanger == "Empty Droplet" ~ "Empty Droplet",
                    !passes_qc                                                      ~ "Low Quality",
                    TRUE                                                            ~ "High Quality")) %>%
    tibble::column_to_rownames(var = "barcode")
  
  # ---- 📥 Add Metadata to Seurat Object ----
  
  # Add the classification vector to the original Seurat object
  sample_seurat <- SeuratObject::AddMetaData(object   = sample_seurat,
                                             metadata = sample_metadata)
                                             
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  
  log_info(sample = sample, 
           step = "calc_qc_metrics",
           message =  glue::glue("Cell-level QC metrics calculated for sample : '{sample}'."))
  return(invisible(sample_seurat))
}

classify_doubletfinder <- function(sample_seurat, n_pcs = NULL, pN = 0.25){
  
  # Set seed for reproducible stochastic processes (PCA, UMAP, Clustering, DoubletFinder)
  set.seed(100)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = sample_seurat, n_pcs = n_pcs, pN = pN)
  
  sample <- as.character(unique(sample_seurat@meta.data$Sample))
  
  # Ensure 'Quality' column is present in metadata
  required_cols <- c("Quality")
  if (!all(required_cols %in% colnames(sample_seurat@meta.data))) {
    log_error(sample = sample, 
             step = "classify_doubletfinder",
             message = "Missing 'Quality' column in metadata. Please run calc_qc_metrics().")
  }
  
  # ---- 🔍 Filter Out Empty Droplets and Low Quality Cells ----
  
  subset_seurat <- subset(x = sample_seurat,
                          subset = (Quality == "High Quality"))
  
  # Ensure enough cells remain
  if (ncol(subset_seurat) < 50) {
    log_warn(sample = sample, 
              step = "classify_doubletfinder",
              message = "Fewer than 50 cells remain after filtering empty droplets. Skipping DoubletFinder.")
    return(invisible(sample_seurat))
  }

  # ---- 🔍 Preprocessing Required for DoubletFinder ----
  
  subset_seurat <- subset_seurat %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    Seurat::FindVariableFeatures(verbose = FALSE) %>%
    Seurat::ScaleData(verbose = FALSE) %>%
    Seurat::RunPCA(npcs = 50, verbose = FALSE)  # Run enough PCs for the sweep/clustering
  
  # ---- 🔍 Find Optimal pK ----
  
  # Artificial doublets are introduced into the real dataset at varying proportions.
  # Data is then preprocessed and the proportion of artificial nearest neighbors 
  # (pANN) is computed for multiple combinations of pK and pN.
  
  # The optimal pK corresponds to the value yielding the maximum bimodality coefficient 
  # (BCmvn), indicating the strongest separation between real and artificial doublets.
  
  # Find significant PCs
  if (is.null(n_pcs)) {
    stdev_pc <- subset_seurat@reductions$pca@stdev
    percent_stdev_pc <- (stdev_pc / sum(stdev_pc)) * 100
    cumulative_stdev_pc <- cumsum(percent_stdev_pc)
    
    pc1 <- which(cumulative_stdev_pc > 90 & percent_stdev_pc < 5)[1]
    pc2 <- sort(which((percent_stdev_pc[1:(length(percent_stdev_pc)-1)] - 
                         percent_stdev_pc[2:length(percent_stdev_pc)]) > 0.1),
                decreasing = TRUE)[1] + 1
    n_pcs <- min(pc1, pc2)
  }
  
  # Run UMAP, FindNeighbors, and FindClusters is only needed for the sweep
  subset_seurat <- subset_seurat %>%
    Seurat::RunUMAP(dims = 1:n_pcs, verbose = FALSE) %>%
    Seurat::FindNeighbors(dims = 1:n_pcs, verbose = FALSE) %>%
    Seurat::FindClusters(resolution = 0.1, verbose = FALSE)
  
  # Run paramSweep
  sweep.res <- quiet_msg(DoubletFinder::paramSweep(subset_seurat, PCs = 1:n_pcs, sct = FALSE))
  sweep.stats <- quiet_msg(DoubletFinder::summarizeSweep(sweep.res, GT = FALSE))
  bcmvn <- quiet_msg(DoubletFinder::find.pK(sweep.stats))
  
  # Find optimal pK using base R functions (safer than dplyr)
  optimal_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  optimal_pK <- as.numeric(as.character(optimal_pK))  # Ensure numeric
  log_info(sample = sample, 
           step = "classify_doubletfinder",
           message = glue::glue("Optimal pK found : '{optimal_pK}'."))
  
  # ---- 🔍 Calculate Adjusted nExp (Expected Doublets) ----
  
  # https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled
  # From the 10x Genomics documentation, the multiplet rate is ~8*10^-6 per cell
  n_cells <- nrow(subset_seurat@meta.data)
  multiplet_rate_per_cell <- 8e-6               # Increment in doublet fraction per additional recovered cell
  multiplet_fraction <- 8e-6 * n_cells          # Fraction of cells expected to be doublets
  n_exp <- round(multiplet_fraction * n_cells)  # Number of doublets
  
  # Adjust for homotypic doublets
  homotypic.prop <- DoubletFinder::modelHomotypic(subset_seurat@meta.data$seurat_clusters)
  n_exp_adj <- round(n_exp * (1 - homotypic.prop))
  log_info(sample = sample, 
           step = "classify_doubletfinder",
           message = glue::glue("Expected doublets (nExp_adj) : '{n_exp_adj}'."))
  
  # ---- 🔍 Run DoubletFinder ----
  
  subset_seurat <- quiet_msg(DoubletFinder::doubletFinder(seu = subset_seurat,
                                                PCs = 1:n_pcs,
                                                pN = pN,  #0.25 default
                                                pK = optimal_pK,
                                                nExp = n_exp_adj))
  
  # Extract DoubletFinder classification column
  df_col <- grep("^DF.classifications", colnames(subset_seurat@meta.data), value = TRUE)
  
  if (length(df_col) != 1) {
    log_error(sample = sample, 
              step = "classify_doubletfinder",
              message = paste("More than 1 DoubletFinder classification columns present.",
                              "Cannot uniquely identify the DoubletFinder classification column.",
                              "Please check if DoubletFinder was already run.",
                              sep = "\n"))
  }
  
  # ---- 🔍 Build Classification Vector ----
  
  # All barcodes present in the Seurat object — includes empty droplets
  all_barcodes <- colnames(sample_seurat)
  
  # Barcodes identified as 'low quality' as per cutoffs
  quality_calls <- sample_seurat@meta.data[["Quality"]]
  low_quality_barcodes <- rownames(sample_seurat@meta.data)[quality_calls == "Low Quality"]

  # Barcodes identified as 'singlets' and 'doublets' by DoubletFinder's algorithm
  df_calls <- subset_seurat@meta.data[[df_col]]
  singlet_barcodes <- rownames(subset_seurat@meta.data)[df_calls == "Singlet"]
  doublet_barcodes <- rownames(subset_seurat@meta.data)[df_calls == "Doublet"]
  
  # Vectorized classification
  classification_vector <- dplyr::case_when(all_barcodes %in% low_quality_barcodes ~ "Low Quality",
                                            all_barcodes %in% doublet_barcodes ~ "Doublet",
                                            all_barcodes %in% singlet_barcodes ~ "Singlet",
                                            TRUE ~ "Empty Droplet")
  
  # Assign names so Seurat maps the metadata correctly by barcode
  names(classification_vector) <- all_barcodes
  
  # ---- 📥 Add Metadata to Seurat Object ----
  
  # Add the classification vector to the original Seurat object
  sample_seurat <- Seurat::AddMetaData(object = sample_seurat,
                                       metadata = classification_vector,
                                       col.name = "DoubletFinder")
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = sample, 
           step = "classify_doubletfinder",
           message =  glue::glue("DoubletFinder doublets identified for sample : '{sample}'."))
  return(invisible(sample_seurat))
}

classify_scdblfinder <- function(sample_seurat){
  
  # Set seed for reproducible stochastic processes 
  set.seed(100)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = sample_seurat)
  
  sample <- as.character(unique(sample_seurat@meta.data$Sample))
 
  # Ensure 'Quality' column is present in metadata
  required_cols <- c("Quality")
  if (!all(required_cols %in% colnames(sample_seurat@meta.data))) {
    stop("Missing 'Quality' classification columns. Please run calc_qc_metrics first.")
  }
  
  # ---- 🔍 Filter Out Empty Droplets and Low Quality Cells ----
  
  subset_seurat <- subset(x = sample_seurat,
                          subset = (Quality == "High Quality"))
  
  # Ensure enough cells remain
  if (ncol(subset_seurat) < 50) {
    log_warn(sample = sample, 
             step = "classify_scdblfinder",
             message = "Fewer than 50 cells remain after filtering empty droplets. Skipping scDblFinder.")
    return(invisible(sample_seurat))
  }
    
  # ---- 🔍 Run scDblFinder ----
  
  # Convert Seurat Object to SingleCellExperiment
  sce <- Seurat::as.SingleCellExperiment(x = subset_seurat)
  
  sce <- quiet_msg(scDblFinder::scDblFinder(sce = sce,
                                  clusters = NULL,
                                  samples = NULL,
                                  dbr = NULL))
  
  # ---- 🔍 Build Classification Vector ----
  
  # All barcodes present in the Seurat object — includes empty droplets
  all_barcodes <- colnames(sample_seurat)
  
  # Barcodes identified as 'low quality' as per cutoffs
  quality_calls <- sample_seurat@meta.data[["Quality"]]
  low_quality_barcodes <- rownames(sample_seurat@meta.data)[quality_calls == "Low Quality"]
  
  # Barcodes identified as 'singlets' and 'doublets' by scDblFinder's algorithm
  df_calls <- base::tolower(as.character(sce$scDblFinder.class))   # scDblFinder classification column
  singlet_barcodes <- colnames(sce)[df_calls == "singlet"]
  doublet_barcodes <- colnames(sce)[df_calls == "doublet"]
  
  # Vectorized classification
  classification_vector <- dplyr::case_when(all_barcodes %in% low_quality_barcodes ~ "Low Quality",
                                            all_barcodes %in% doublet_barcodes ~ "Doublet",
                                            all_barcodes %in% singlet_barcodes ~ "Singlet",
                                            TRUE ~ "Empty Droplet")
  
  # Assign names so Seurat maps the metadata correctly by barcode
  names(classification_vector) <- all_barcodes
  
  # ---- 📥 Add Metadata to Seurat Object ----
  
  # Add the classification vector to the original Seurat object
  sample_seurat <- Seurat::AddMetaData(object = sample_seurat,
                                       metadata = classification_vector,
                                       col.name = "scDblFinder")
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = sample, 
           step = "classify_scdblfinder",
           message = glue::glue("scDblFinder doublets identified for sample : '{sample}'."))
  return(invisible(sample_seurat))
}

filter_singlets <- function(sample_seurat){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = sample_seurat)
  
  sample <- as.character(unique(sample_seurat@meta.data$Sample))
  
  # ---- 📥 Populate missing optional metadata columns ----
  
  # NOTE: Use "ND" instead of NA for safer subsetting
  optional_cols <- c("DropletUtils", "CellRanger", "HTO_Final")
  for (col in optional_cols) {
    if (!col %in% colnames(sample_seurat@meta.data)) {
      sample_seurat <- Seurat::AddMetaData(object   = sample_seurat,
                                           metadata = "ND",
                                           col.name = col)
    }
  }
  
  # ---- 🧪 Ensure QC columns exist ----
  
  required_cols <- c("Quality", "DoubletFinder", "scDblFinder")
  if (!all(required_cols %in% colnames(sample_seurat@meta.data))) {
    log_error(sample = sample, 
             step = "filter_singlets",
             message = paste("Missing required columns in metadata: Quality, DoubletFinder, scDblFinder.",
                             "Please run calc_qc_metrics(), classify_doubletfinder() and classify_scdblfinder()",
                             sep = "\n"))
  }
  
  # Columns to retain (if available)
  keep_cols <- c("barcode", "Cell", "Sample", "QC", "nUMIs", "nGenes", 
                 "MitoRatio", "RiboRatio", "HemeRatio", "Novelty", "Quality",
                 "DropletUtils", "CellRanger", "DoubletFinder", "scDblFinder",
                 "nHTO_UMIs", "nHTOs", "HTO_Final", "X", "Y")
  
  # ---- 🔍 Subset to retain only high-quality singlets ----
  
  # QC CONFIDENCE STRATEGY (Final Classification)
  # This strategy defines the final 'QC' status based on a prioritized hierarchy 
  # designed to balance specificity (for retention) and sensitivity (for removal).
  #
  # 1. Empty Droplets: Require BOTH DropletUtils AND CellRanger to agree ('&') 
  #    to maximize **Specificity**. High confidence ensures true empty droplets 
  #    are removed while avoiding accidental exclusion of low-UMI real cells.
  #
  # 2. Low Quality: Flagged by failure of initial quality metrics (e.g., MitoRatio, nUMIs).
  #    These cells are prioritized for removal over doublet calls.
  #
  # 3. Doublets: Remove if EITHER DoubletFinder OR scDblFinder calls doublet ('|') 
  #    to maximize **Sensitivity**. This ensures all suspicious cells are excluded, 
  #    preventing doublets from contaminating downstream analyses.
  #
  # 4. Singlets: Any cell passing the Empty Droplet, Low Quality, and Doublet criteria.
  
  # Assign final classification (QC)
  metadata <- sample_seurat@meta.data %>%
    dplyr::mutate(QC = dplyr::case_when(DropletUtils == "Empty Droplet" & CellRanger == "Empty Droplet" ~ "Empty Droplet",
                                        Quality == "Low Quality" ~ "Low Quality",
                                        DoubletFinder == "Doublet" | scDblFinder == "Doublet" ~ "Doublet",
                                        TRUE ~ "Singlet")) 
  
 
  # Select available columns
  available_cols <- base::intersect(keep_cols, colnames(metadata))
  metadata <- metadata %>%
    dplyr::select(all_of(available_cols)) 
  
  # Assign the cleaned metadata back to the Seurat object
  sample_seurat@meta.data <- metadata
  
  # Subset high quality singlets
  sample_seurat <- base::subset(x = sample_seurat, QC == "Singlet")
  
  # ---- 🪵 Log Output and Return Seurat Object, raw_metadata ----
  
  log_info(sample = sample,
           step = "filter_singlets",
           message = glue::glue("Retained high-quality singlets for sample : '{sample}'."))
  return(list(sample_seurat = sample_seurat,
              raw_metadata  = metadata))
}

merge_filtered <- function(seurat_list, assay, metafile, output_dir){
  
  # Set seed for reproducible stochastic processes
  set.seed(100)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_list = seurat_list, assay = assay, 
                  metafile = metafile, output_dir = output_dir)
  
  # ---- 🤝 Merge Seurat Objects ----
 
  # Create sample ids from sample names
  samples <- names(seurat_list)
  sample_ids <- gsub(pattern = ".Spatial.*", replacement = "", x = samples)
  
  # NOTE: Since multiple samples can have same barcodes, add sample name to 
  # barcodes to keep track of cell identities (i.e., barcodes) coming from each 
  # sample after merging
  merged_seurat <- tryCatch({
    merge(x = seurat_list[[1]],       # get(paste0(samples[1])
          y = seurat_list[-1],        # lapply(paste0(samples[2:length(samples)]), get)
          add.cell.ids = sample_ids,  # add unique prefix to barcodes
          merge.data = FALSE)
  }, error = function(e){
    log_error(sample = "",
             step = "merge_filtered",
             message = glue::glue("Merge failed : '{e$message}'."))
  }) 
  
  # ---- 🗑️ Remove HTO assay (if present) ----
  
  if ("HTO" %in% Seurat::Assays(merged_seurat)){
    merged_seurat[["HTO"]] <- NULL
    log_info(sample = "",
              step = "merge_filtered",
              message = "Removed HTO assay prior to integration.")
  }
  
  # ---- 📥 Load Extra Metadata ----
  
  extra_metadata <- tryCatch({
    openxlsx::read.xlsx(xlsxFile = metafile) %>%
      dplyr::select(-dplyr::any_of("Comments"))
  }, error = function(e){
    log_warn(sample = "",
             step = "merge_filtered",
             message = glue::glue("Failed to read metadata file : '{e$message}'."))
    return(data.frame())  # Return empty dataframe to skip merge
  })
  
  # ---- 🔗 Join Extra Metadata ----
  
  if (nrow(extra_metadata) > 0){
    
    meta_data <- merged_seurat@meta.data %>%
      tibble::rownames_to_column(var = "temp_barcode_id")
    
    if (assay == "RNA"){
      # NOTE: For droplet-based RNA, join using Sample + HTO if available
      meta_data <- meta_data %>%
        dplyr::mutate(Unique_ID = dplyr::case_when(HTO_Final != "ND" ~ paste0(Sample, "_", HTO_Final),
                                                   TRUE ~ Sample)) %>%
        dplyr::left_join(extra_metadata, by = c("Unique_ID" = "Unique_ID"))
    } else {
      # NOTE: For Spatial assays, join by Slide/Sample
      meta_data <- meta_data %>%
        dplyr::left_join(extra_metadata, by = c("Sample" = "Slide"))
      
      log_info(sample = "",
               step = "merge_filtered",
               message = "Extra metadata joined for Spatial analysis. No cell filtering performed.")
    }
    
    # ---- 🧹 Clean Metadata ----
    
    # Remove columns that are entirely NA
    meta_data <- meta_data %>%
      dplyr::select(where(~ !all(is.na(.))),         # keep columns that are not entirely NA
                    -dplyr::any_of(c("Initial filename", "Comments")))  # remove unwanted columns
      
    # Re-assign metadata with original barcodes as rownames
    merged_seurat@meta.data <- meta_data %>%
      tibble::column_to_rownames(var = "temp_barcode_id")
    
    log_info(sample = "",
             step = "merge_filtered",
             message = "Extra metadata added successfully.")
  } else {
    log_info(sample = "",
             step = "merge_filtered",
             message = "No external metadata was loaded or joined.")
  }
  
  # ---- 💾 Save Merged Seurat Object ----
  
  # Determine file name based on assay
  filename <- if (assay == "RNA"){
    "filtered_seurat.rds"
  } else{   # For Spatial.008um and Spatial.016um assays
    paste0("filtered_seurat.", assay, ".rds")
  }
  
  base::saveRDS(object = merged_seurat, file = file.path(output_dir, filename))
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = "",
           step = "merge_filtered",
           message = glue::glue("Filtered Seurat object saved to : '{filename}'."))
  return(invisible(merged_seurat))
}

plot_qc <- function(raw_metadata, output_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(raw_metadata = raw_metadata, output_dir = output_dir)
  
  # Ensure required columns exist in metadata
  required_cols <- c("Sample", "QC", "nUMIs", "nGenes", "MitoRatio", "RiboRatio", "Novelty")
  missing_cols <- setdiff(required_cols, colnames(raw_metadata))
  if (length(missing_cols) > 0) {
    log_error(sample = "",
             step = "plot_qc",
             message = glue::glue("Missing required metadata columns : '{paste(missing_cols, collapse = ", ")}'."))
  }
  
  # ---- 🛠️ Helper Functions ----
  
  # Define classification levels and color palette
  qc_levels <- c("Empty Droplet", "Doublet", "Singlet", "Low Quality")
  fill_colors <- c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
                   "Doublet" = "#1F77B4", "Low Quality" = "#D62728")
  
  # Visualize cell counts per Sample
  cell_qc <- function(meta){
    
    # Count cells by Sample and QC category and fill missing levels
    df <- meta %>%
      dplyr::count(Sample, QC) %>%
      data.frame() %>%
      dplyr::mutate(QC = factor(QC, levels = qc_levels)) %>% 
      tidyr::complete(Sample, QC, fill = list(n = 1))
    
    # Numeric positions for vertical lines between samples
    sample_levels <- levels(factor(df$Sample))
    vline_positions <- seq(1.5, length(sample_levels) - 0.5, by = 1)
    
    # NOTE: position = "dodge" for grouped bars; "stack" for stacked bar
    # stat = "identity" because y is precomputed; "count" if y axis determined based on X axis frequency
    p <- ggplot(data = df, aes(x = Sample, y = n, fill = QC)) +
      geom_bar(stat = "identity", position = position_dodge(0.9)) +
      theme_classic() +
      custom_theme +
      labs(x = "Sample", y = "Cell Counts", title = "Number of Cells") +
      coord_cartesian(ylim = c(1,1e7), clip = "off", expand = FALSE) +
      scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000, 1000000)) +
      scale_fill_manual(values = fill_colors) +
      #geom_text(stat ="count", aes(label = after_stat(count)), y = 0, hjust = 0, angle = 90)
      geom_text(mapping = aes(label = n, ymin = 0.1, ymax = 1),
                position = position_dodge(width = 0.9),
                y = 0.1, hjust = 0, angle = 90) +
      geom_vline(xintercept = vline_positions, linetype = "dotted", color = "gray50")
    
    return(p)
  }
  
  # Visualize nUMIs, nGenes, MitoRatio, RiboRatio, Novelty per sample
  violin_qc <- function(meta, yvar, ylab, title, cutoff = NULL, ylog = TRUE, ylim = NULL){
    
    # Fill missing levels
    df <- meta %>% 
      dplyr::mutate(QC = factor(QC, levels = qc_levels))
    
    # Numeric positions for vertical lines between samples
    sample_levels <- levels(factor(df$Sample))
    vline_positions <- seq(1.5, length(sample_levels) - 0.5, by = 1)
    
    p <- ggplot(data = df, aes(x = Sample, y = .data[[yvar]], fill = QC)) +
      geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) +
      geom_boxplot(position = position_dodge(0.9), width = 0.05, outlier.size = 0.5) +
      theme_classic() + 
      custom_theme +
      labs(x = "Sample", y = ylab, title = title) +
      scale_fill_manual(values = fill_colors) +
      geom_vline(xintercept = vline_positions, linetype = "dotted", color = "gray50")
    
    if (!is.null(cutoff)) p <- p + geom_hline(yintercept = cutoff, linetype = 2)
    if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim, clip = "off")
    if (ylog) p <- p + scale_y_log10()
    
    return(p)
  }
  
  # Visualize number of genes/cell, number of UMIs/cell & MitoRatio together.
  # Bottom left quadrant : Poor quality cells with low genes & UMIs per cell 
  # Top right quadrant   : Good quality cells with high genes & UMIs per cell
  # Bottom right quadrant: Cells with low genes but high UMIs per cell. These 
  # could be dying cells or population of low complexity cells (i.e erythrocytes)
  gene_umi_mito_qc <- function(meta){
    
    umi_cutoff <- 500
    gene_cutoff <- 250
    
    df <- meta %>%
      dplyr::mutate(QC = factor(QC, levels = qc_levels))
    
    ggplot(data = df, aes(x = nUMIs, y = nGenes, color = MitoRatio)) +
      geom_point(alpha = 0.5, size = 0.25) +
      theme_classic() +
      custom_theme +
      labs(x = "Number of UMIs", y = "Number of Genes", title = "UMIs vs Genes (Colored by MitoRatio)") +
      coord_cartesian(xlim = c(1,1e6), ylim = c(1,20000), clip = "off") +
      scale_x_log10(breaks = c(1,10,100,1000,10000,100000,1e6)) + 
      scale_y_log10(breaks = c(1,10,100,1000,10000,100000)) + 
      scale_color_viridis(option = "D", limits = c(0,1)) +
      #facet_wrap(.~Sample, nrow = 4) +   #split the plot by X-axis label
      facet_wrap(Sample ~ QC, ncol = 4, drop = FALSE) +
      geom_vline(xintercept = umi_cutoff, linetype = "dashed") +
      geom_hline(yintercept = gene_cutoff, linetype = "dashed") +
      stat_smooth(method = lm, color = "yellow", se = FALSE)
  }
  
  # ---- 💾 Generate and Save Plots ----

  # Plot generators
  plot_list <- list(
    Cell_Counts                      = cell_qc,
    UMI_Distribution                 = function(x) violin_qc(x, "nUMIs", "Number of UMIs", "UMI Distribution", 500, TRUE, c(1, 1e6)),
    Gene_Distribution                = function(x) violin_qc(x, "nGenes", "Number of Genes", "Gene Distribution", 250, TRUE, c(1, 30000)),
    MitoRatio_Distribution           = function(x) violin_qc(x, "MitoRatio", "MitoRatio", "MitoRatio Distribution", 0.2, FALSE, c(0, 1)), 
    RiboRatio_Distribution           = function(x) violin_qc(x, "RiboRatio", "RiboRatio", "RiboRatio Distribution", 0.05, FALSE, c(0, 1)),
    Novelty_Score_Distribution       = function(x) violin_qc(x, "Novelty", "Novelty", "Novelty Score Distribution", 0.8, FALSE, c(0, 1)),
    Genes_UMI_MitoRatio_Distribution = gene_umi_mito_qc)
  
  # Generate plots
  for (plot_name in names(plot_list)) {
    p <- plot_list[[plot_name]](raw_metadata)    #  p <- get(funcs[i])(raw_metadata)
    
    # Find number of samples
    n_samples <- length(unique(raw_metadata$Sample))
    
    # Set plot width based on number of samples, e.g., 0.8 inch per sample
    max_pdf_width <- 30
    max_pdf_height <- 30
    
    # Default: 0.8 inch per sample
    desired_width_per_sample <- 0.8
    pdf_width <- min(n_samples * desired_width_per_sample + 2, max_pdf_width)  # 2 inch for margin/legend
    pdf_height <- 8

    # Special case for Genes_UMI_MitoRatio_Distribution (4 panels per sample)
    if (plot_name == "Genes_UMI_MitoRatio_Distribution"){
      desired_width_per_panel <- 2
      desired_height_per_sample <- 2
      pdf_width <- min(4 * desired_width_per_panel + 2, max_pdf_width)
      pdf_height <- min(n_samples * desired_height_per_sample + 2, max_pdf_height)  # 2 inch for axis etc
    }
    
    # Optional: enforce minimum dimensions
    pdf_width  <- max(pdf_width, 6)
    pdf_height <- max(pdf_height, 5)
    
    # Save plot
    ggplot2::ggsave(filename = paste0("QC_", plot_name, ".pdf"),
                    plot = p,
                    device = "pdf",
                    path = output_dir,
                    width =  pdf_width,
                    height = pdf_height,
                    dpi = 600,
                    units = "in")
  }
  
  # ---- 🪵 Log Output ----
  
  log_info(sample = "",
           step = "plot_qc",
           message = glue::glue("QC plots generated and saved to : '{output_dir}'."))
}

run_sctransform <- function(filtered_seurat, assay, s_genes, g2m_genes){
  
  # Set seed for reproducible stochastic processes (PCA, UMAP)
  set.seed(1234)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = filtered_seurat, assay = assay,
                  s_genes = s_genes, g2m_genes = g2m_genes)
  
  # ---- 🔍 Check for Cell Cycle Genes ----
  
  # Step 1: Remove duplicates ignoring case (currently mix of human and mouse)
  s_genes_unique <- s_genes[!duplicated(toupper(s_genes))]
  g2m_genes_unique <- g2m_genes[!duplicated(toupper(g2m_genes))]
  
  # Step 2: keep only genes present in the Seurat object
  #present_features <- rownames(filtered_seurat@assays[[assay]])
  present_features <- SeuratObject::Features(filtered_seurat, assay = assay)
  s_genes <- intersect(s_genes, present_features)
  g2m_genes <- intersect(g2m_genes, present_features)
  
  if (length(s_genes) == 0) {
    pad <- base::strrep(" ", 29)
    log_error(sample = "",
              step = "run_sctransform",
              message = glue::glue("Zero S-phase genes were found in '{assay}' assay\n",
                                   "{pad}Cannot perform cell cycle scoring."))
                              
  } else {
    log_info(sample = "",
             step = "run_sctransform",
             message = glue::glue("{length(s_genes)} of {length(s_genes_unique)} S-phase genes found and will be used for scoring."))
  }
  
  if (length(g2m_genes) == 0) {
    pad <- base::strrep(" ", 29)
    log_error(sample = "",
              step = "run_sctransform",
              message = glue::glue("Zero G2M-phase genes were found in '{assay}' assay\n",
                                   "{pad}Cannot perform cell cycle scoring."))
  } else {
    log_info(sample = "",
             step = "run_sctransform",
             message = glue::glue("{length(g2m_genes)} of {length(g2m_genes_unique)} G2M-phase genes found and will be used for scoring."))
  }

  # ---- SCTransfrom Workflow ----
  
  # 1️⃣ Normalize Data using LogNormalize (for cell cycle scoring)
  filtered_seurat <- Seurat::NormalizeData(object = filtered_seurat,
                                           assay = assay,
                                           normalization.method = "LogNormalize",
                                           scale.factor = 10000,
                                           margin = 1,
                                           verbose = FALSE)
  
  # 2️⃣ Join Layers (for complete data view before cell cycle scoring)
  # NOTE: CellCycleScoring uses a single data layer (log-normalized counts), 
  # but currently, data layers for each sample may be stored separately. Join them first.
  filtered_seurat@assays[[assay]] <- SeuratObject::JoinLayers(object = filtered_seurat@assays[[assay]])
  
  # 3️⃣ Score Cell Cycle
  filtered_seurat <- Seurat::CellCycleScoring(object = filtered_seurat,
                                              s.features = s_genes,
                                              g2m.features = g2m_genes,
                                              ctrl = NULL)
  
  # 4️⃣ Calculate CC.Score to regress out cell cycle differences
  # https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette
  filtered_seurat$CC.Score <- filtered_seurat$G2M.Score - filtered_seurat$S.Score
 
  
  # 5️⃣ Split Object by 'Sample' (for batch-aware SCTransform)
  # NOTE: All cells within the same batch MUST be analyzed together
  filtered_seurat@assays[[assay]] <- base::split(x = filtered_seurat@assays[[assay]],
                                                 f = filtered_seurat@meta.data[["Sample"]])
  
  # 6️⃣ Run SCTransform (with CC.Score and MitoRatio regression)
  # https://github.com/satijalab/seurat/issues/7342
  sct_seurat <- quiet_msg(Seurat::SCTransform(object = filtered_seurat,
                                    assay = assay,
                                    new.assay.name = "SCT",
                                    do.correct.umi = TRUE,
                                    ncells = 5000,
                                    variable.features.n = 3000,
                                    vars.to.regress = c("CC.Score", "MitoRatio"),
                                    do.scale = FALSE,
                                    do.center = TRUE,
                                    vst.flavor = "v2",
                                    return.only.var.genes = TRUE,
                                    verbose = FALSE))
  
  # 7️⃣ Prepare SCT Assay for Differential Expression
  sct_seurat <- Seurat::PrepSCTFindMarkers(object = sct_seurat,
                                           assay = "SCT",
                                           verbose = FALSE)
  
  # 8️⃣ Filter Out Unwanted Variable Features (Ribo, Mito, etc.)
  # NOTE: PCA, UMAP, and clustering should not be influenced by these genes.
  var_f <- Seurat::VariableFeatures(sct_seurat, assay = "SCT")
  var_f <- var_f[!grepl(pattern = "^[Rr][Pp][SsLl]|R[Ii][Kk]$|^[Mm][Tt]-|^G[Mm][0-9.]+$", x = var_f)]
  Seurat::VariableFeatures(sct_seurat, assay = "SCT") <- var_f
  #sct_seurat@assays[["SCT"]]@var.features <- var_f
  cat("\nFinal number of variable features:", length(var_f), "\n")
  
  # 9️⃣ Scale and Run PCA on Original Assay (for export compatibility)
  # NOTE: Populates 'scale.data' slot and generates PCA reduction on log-normalized data.
  # This is REQUIRED for seamless export to AnnData/Scanpy/scVI environments.
  sct_seurat <- Seurat::ScaleData(object = sct_seurat,
                                  assay = assay,
                                  features = Seurat::VariableFeatures(sct_seurat, assay = "SCT"), # Use SCT features
                                  verbose = FALSE)
  
  # Use a consistent naming scheme for dimensional reduction keys. See NOTE on 🔟.
  sct_seurat <- Seurat::RunPCA(object = sct_seurat,
                               assay = assay,
                               features = Seurat::VariableFeatures(sct_seurat, assay = "SCT"),
                               reduction.name = paste0(tolower(assay), "_pca"),
                               reduction.key = paste0(toupper(assay), "pca_"),
                               verbose = FALSE)
  
  # 🔟 Run PCA on SCT Assay
  # NOTE: IntegrateLayers() sets reduction.key as "integcca_" when reduction.name = "integ_cca".
  # For a consistent naming, we set reduction.key as "sctpca_" when reduction.name = "sct_pca".
  sct_seurat <- Seurat::RunPCA(object = sct_seurat,
                               assay = "SCT",
                               features = Seurat::VariableFeatures(sct_seurat, assay = "SCT"),
                               reduction.name = "sct_pca",
                               reduction.key = "sctpca_", 
                               verbose = FALSE)
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = "",
           step = "run_sctransform",
           message = "SCTransform workflow completed successfully.")
  return(invisible(sct_seurat))
}

integrate_sct_data <- function(sct_seurat, assay, reference_samples=NULL){
  
  # Set seed for reproducible stochastic processes (Harmony, PCA)
  set.seed(1234)
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = sct_seurat, assay = assay,
                  reference_samples = reference_samples)
  
  # ---- 🔍 Ensure 'sct_pca' reduction exists ----
  if (!"sct_pca" %in% names(sct_seurat@reductions)) {
    
    log_error(sample = "",
             step = "integrate_sct_data",
             message = "Reduction 'sct_pca' is missing. Please run PCA on SCT assay before integration.")
  }
  
  # ---- ⚖️ Compute optimal k.weight for integration ----
  
  # NOTE: k.weight is typically half the number of cells in the smallest sample, capped at 100
  kweight <- min(sct_seurat@meta.data %>% 
                   dplyr::count(Sample) %>% 
                   dplyr::pull(n) %>% 
                   min()/2, 100) 
  log_info(sample = "",
           step = "integrate_sct_data",
           message = glue::glue("Dataset k.weight : '{kweight}'."))
  
  # ---- 📏 Determine maximum number of dimensions for integration ----
  
  # NOTE: Seurat::JointPCAIntegration() gives errors if fewer than 30 dims are provided
  max_dims <- min(30, ncol(Seurat::Embeddings(sct_seurat, reduction = "sct_pca")))
  #max_dims <- min(30, ncol(sct_seurat@reductions$sct_pca@cell.embeddings))
  
  # ---- 🔁 Loop over Integration Methods ----
  
  # Methods supported:
  #   - CCA      : Canonical Correlation Analysis, aligns datasets based on correlated components
  #   - RPCA     : Reciprocal PCA, preserves sample-specific variance while integrating
  #   - Harmony  : Batch correction in PCA space, fast and scalable
  #   - JointPCA : Joint PCA integration, combines PCA embeddings from multiple datasets
  
  # Initialize Seurat object for integration
  integrated_seurat <- sct_seurat
  Seurat::DefaultAssay(integrated_seurat) <- "SCT"
  
  integration_methods <- c("CCA", "RPCA", "Harmony", "JointPCA")
  
  for (method in integration_methods){
    
    # Set the name of the reduction after integration
    reduction_name <- paste0("integ_", base::tolower(method))
    
    # Perform integration using IntegrateLayers
    integrated_seurat <- Seurat::IntegrateLayers(object = integrated_seurat,
                                                 method = paste0(method, "Integration"),
                                                 normalization.method = "SCT",
                                                 orig.reduction = "sct_pca", 
                                                 new.reduction = reduction_name,
                                                 reference = reference_samples,
                                                 k.weight = kweight,    # relevant for RPCA
                                                 dims = 1:max_dims,
                                                 verbose = FALSE)
  }
  
  # ---- Optional integration for scVI and FastMNN ----
  
  # NOTE: scVI needs raw counts. Vignette also uses it on RNA assay
  # NOTE: We use variable features of SCT assay for integration.
  # NOTE: We use pca reduction from RNA assay (derived using variable features of SCT assay)
  # FastMNN throws error "Error in checkBatchConsistency(batches, cells.in.columns = TRUE)"
  
  # for (method in c("scVI", "FastMNN")) {
  #   DefaultAssay(integrated_seurat) <- assay
  #   integrated_seurat <- Seurat::IntegrateLayers(object = integrated_seurat,
  #                                                method = paste0(r, "Integration"),
  #                                                normalization.method = "LogNormalize",
  #                                                orig.reduction = paste0(assay, ".pca"),
  #                                                features = integrated_seurat@assays$SCT@var.features,
  #                                                new.reduction = base::tolower(method),
  #                                                reference = reference_samples,
  #                                                k.weight = kweight,                                    # for RPCA
  #                                                conda_env = "/hpc/home/kailasamms/miniconda3/envs/R",  # for scVI
  #                                                verbose = FALSE)
  # }
  
  # ---- 🔗 Merge layers after integration (for RNA/Spatial assay) ----
  
  # NOTE: Only applicable for RNA or Spatial assays, not for SCT assay itself
  integrated_seurat@assays[[assay]] <- SeuratObject::JoinLayers(integrated_seurat@assays[[assay]])
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = "",
           step = "integrate_sct_data",
           message = "Integration completed successfully.")
  return(invisible(integrated_seurat))
}

cluster_integrated_data <- function(integrated_seurat, assay){
  
  # Set seed for reproducible stochastic processes (Leiden clustering, UMAP)
  set.seed(1234)
  
  # ---- ⚙️ Validate Input Parameters ----
  validate_inputs(seurat_object = integrated_seurat, assay = assay)
  
  # ---- 🔬 Find Nearest Neighbors (for every cell) ----
  
  # Integration methods used during integration
  integration_methods <- c("CCA", "RPCA", "Harmony", "JointPCA")
  
  # Define clustering resolutions to explore
  resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)

  # Loop over integration methods to calculate the nearest neighbors
  for (method in integration_methods) {
    
    # Define name of reduction, NN and SNN graphs after integration
    reduction_name <- paste0("integ_", base::tolower(method))
    nn_name <- paste0(base::tolower(method), ".nn")   # kNN graph
    snn_name <- paste0(base::tolower(method), ".snn") # Shared nearest neighbor (SNN) graph
    
    # Determine max dimensions available for this reduction
    max_dims <- min(40, ncol(SeuratObject::Embeddings(integrated_seurat, reduction = reduction_name)))
    #max_dims <- min(40, ncol(integrated_seurat@reductions[[reduction_name]]@cell.embeddings))
    
    # Ensure SCT assay is active
    DefaultAssay(integrated_seurat) <- "SCT"
    
    # Compute nearest neighbors
    integrated_seurat <- Seurat::FindNeighbors(object = integrated_seurat,
                                               reduction = reduction_name,
                                               dims = 1:max_dims,
                                               k.param = 30,   # Default k for kNN
                                               graph.name = c(nn_name, snn_name),
                                               verbose = FALSE)
  }
  
  # ---- 🧩 Find Clusters using the SNN graph ----
  
  # Loop over integration methods and clustering resolutions
  for (method in integration_methods) {
    for (res in resolutions) {
      
      snn_name <- paste0(base::tolower(method), ".snn")    # SNN graph used for clustering
      cluster_name <- paste0(base::tolower(method), res)   # Cluster label name
      
      # Leiden clustering
      integrated_seurat <- Seurat::FindClusters(object = integrated_seurat,
                                                resolution = res,
                                                graph.name = snn_name,
                                                cluster.name = cluster_name,
                                                modularity.fxn = 1,
                                                algorithm = 4,     # 4 = Leiden algorithm (recommended)
                                                verbose = FALSE)
    }
  }
  
  # ---- 🗺️ Run UMAP for dimensionality reduction (visualization) ----
  
  for (method in integration_methods) {
    
    # Define the name of the reduction after integration
    reduction_name <- paste0("integ_", base::tolower(method))
    
    # Determine max dimensions available for this reduction
    max_dims <- min(40, ncol(SeuratObject::Embeddings(integrated_seurat, reduction = reduction_name)))
    #max_dims <- min(40, ncol(integrated_seurat@reductions[[reduction_name]]@cell.embeddings))
    
    integrated_seurat <- Seurat::RunUMAP(object = integrated_seurat,
                                         dims = 1:max_dims,
                                         n.neighbors = 30L,   # Number of neighbors for UMAP graph
                                         reduction = reduction_name,
                                         reduction.name = paste0("umap_", base::tolower(method)),
                                         verbose = FALSE)
  }
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = "",
           step = "cluster_integrated_data",
           message = "Clustering completed successfully.")
  return(invisible(integrated_seurat))
}

remove_sparse_clusters <- function(integrated_seurat, assay, output_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  validate_inputs(seurat_object = integrated_seurat, assay = assay,
                  output_dir = output_dir)
  
  # ---- 🔎 Identify all cluster columns dynamically ----
 
  # Define all possible starting prefixes based on your integration methods
  prefixes <- c("cca", "rpca", "harmony", "jointpca")
  pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
  
  # Select cluster columns that contain numeric cluster IDs
  cluster_cols <- colnames(integrated_seurat@meta.data) %>%
    stringr::str_subset(pattern) %>%
    stringr::str_subset("[0-9]")
  
  if (length(cluster_cols) == 0) {
    log_error(sample = "",
             step = "remove_sparse_clusters",
             message = glue::glue("No clustering columns found with prefix : '{pattern}'."))
  }
  
  # ---- 🗑️ Find and Remove Sparse Cells ----
  
  sparse_cells <- character(0)
  
  for (col in cluster_cols) {
    
    # Count number of cells per cluster in this column 
    # NOTE: Use table() for fast frequency counting
    cluster_counts <- table(integrated_seurat@meta.data[[col]])
    
    # Identify clusters with very few cells (<= 5)
    sparse_clusters_ids <- names(cluster_counts[cluster_counts <= 5])
    
    # If sparse clusters exist, identify the cells belonging to them
    if (length(sparse_clusters_ids) > 0) {
      
      # Get barcodes of cells in sparse clusters
      cells_in_sparse_clusters <- rownames(integrated_seurat@meta.data)[integrated_seurat@meta.data[[col]] %in% sparse_clusters_ids]
      
      # Append these cells to the master sparse cell list
      sparse_cells <- c(sparse_cells, cells_in_sparse_clusters)
      
      # Optional progress message
      log_info(sample = "",
               step = "remove_sparse_clusters",
               message = glue::glue("Found '{length(cells_in_sparse_clusters)}' cells to remove from column : '{col}'."))
    }
  }
  
  # Identify barcodes in non-sparse clusters
  unique_sparse_cells <- unique(sparse_cells)
  all_cells <- SeuratObject::Cells(x = integrated_seurat)
  keep_cells <- base::setdiff(all_cells, unique_sparse_cells)
  
  # Subset Seurat object safely using barcode name (i.e. rownames)
  if (length(unique_sparse_cells) > 0) {
    integrated_seurat <- subset(x = integrated_seurat,
                                cells = keep_cells)
  }
  
  log_info(sample = "",
           step = "remove_sparse_clusters",
           message = glue::glue("Total cells removed from sparse clusters : '{length(unique_sparse_cells)}'."))

  # ---- 💾 Save integrated Object ----
  
  # Determine file name based on assay
  filename <- if (assay == "RNA") {
    "integrated_seurat.rds"
  } else {
    paste0("integrated_seurat.", assay, ".rds")
  }
  
  base::saveRDS(integrated_seurat, file = file.path(output_dir, filename))
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = "",
           step = "remove_sparse_clusters",
           message = "Successfully saved integrated seurat object after removing sparse cells.")
  return(invisible(integrated_seurat))
}
  
calc_optimal_resolution <- function(integrated_seurat, reduction, output_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = integrated_seurat, reduction = reduction,
                  output_dir = output_dir)
  
  # ---- 🔎 Identify all cluster columns dynamically ----
  
  # Define all possible starting prefixes based on your integration method
  prefixes <- base::tolower(reduction)
  pattern <- stringr::str_remove(prefixes, "^(umap_|integrated_|sct_|pca_)*")
  
  # Select cluster columns that contain numeric cluster IDs
  cluster_cols <- colnames(integrated_seurat@meta.data) %>%
    stringr::str_subset(pattern) %>%
    stringr::str_subset("[0-9]")
  
  if (length(cluster_cols) == 0) {
    log_error(sample = "",
             step = "calc_optimal_resolution",
             message = glue::glue("No clustering columns found with prefix : '{pattern}'."))
  }
  
  # ---- 🌳 Visualize Clustering Stability (Clustree) ----
  
  # Sort cluster columns by numeric resolution (REQUIRED for SC3 stability scoring)
  cluster_cols <- cluster_cols[order(as.numeric(stringr::str_extract(cluster_cols, "[0-9.]+")))]
  
  # Run clustree
  clustree_plot <- clustree::clustree(x = integrated_seurat,
                                      prefix = pattern)
  
  # Save clustree plot
  ggplot2::ggsave(filename = paste0("Clustree_", pattern, ".pdf"),
                  plot = clustree_plot,
                  path = output_dir,
                  device = "pdf",
                  width = 8,
                  height = 11,
                  dpi = 300,
                  units = "in")
  
  # ---- 📈 Quantify Clustering Stability (SC3) ----
  
  # Convert factor/character cluster columns to integer IDs
  cluster_df <- integrated_seurat@meta.data[, cluster_cols, drop = FALSE]
  cluster_df[] <- lapply(cluster_df, function(x) as.integer(as.factor(x)))
  cluster_matrix <- as.matrix(cluster_df)
  
  # Calculate SC3 stability
  # NOTE: stability_mat uses resolution index like 1, 2, 3 instead of actual
  # resolution names like "harmony0.4" etc...
  stability_mat <- clustree:::calc_sc3_stability(cluster_matrix)
  
  # Compute mean stability per resolution and map resolution names
  stability_df <- stability_mat %>%
    as.data.frame() %>%
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric)) %>%
    dplyr::group_by(resolution) %>%
    dplyr::summarise(mean_stability = mean(stability, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(resolution_name = colnames(cluster_matrix)[resolution])

  # Pick the resolution with highest mean stability
  optimal_res <- stability_df %>%
    dplyr::slice_max(mean_stability) %>%
    dplyr::pull(resolution_name)
  
  optimal_res_numeric <- as.numeric(stringr::str_extract(optimal_res, "[0-9.]+"))
  
  # Visualize stability across resolutions
  stability_plot <- ggplot(stability_df, aes(x = resolution_name, y = mean_stability, group = 1)) +
    geom_line(color = "blue") +
    geom_point(size = 3, color = "red") +
    theme_classic() +
    labs(title = "Mean Cluster Stability per Resolution",
         x = "Resolution",
         y = "Mean SC3 Stability")
  
  # Save stability scores plot
  ggplot2::ggsave(filename = paste0("Clustree_Stability_scores_", pattern, ".pdf"),
                  plot = stability_plot,
                  path = output_dir,
                  device = "pdf",
                  width = 8,
                  height = 11,
                  dpi = 300,
                  units = "in")
  
  # ---- 🪵 Log Output and Return Seurat Object ----
  
  log_info(sample = "",
           step = "calc_optimal_resolution",
           message = glue::glue("Optimal resolution identified : '{optimal_res_numeric}'."))
  return(invisible(optimal_res_numeric))
}

identify_markers <- function(integrated_seurat, resolution, reduction, output_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = integrated_seurat, resolution = resolution,
                  reduction = reduction, output_dir = output_dir)
  
  # ---- 🧪 Detect Appropriate Assay for FindAllMarkers ----
  
  all_assays <- names(integrated_seurat@assays)
  
  if ("SCT" %in% all_assays) {
    active_assay <- "SCT"
  } else if ("RNA" %in% all_assays) {
    active_assay <- "RNA"
  } else if (length(all_assays) > 0) {
    active_assay <- all_assays[1]
    log_info(sample = "",
              step = "identify_markers",
              message = glue::glue("No SCT/RNA assay found. Using assay : '{active_assay}'."))
  } else {
    log_error(sample = "",
              step = "identify_markers",
              message = "No assays found in Seurat object for DE analysis.")
  }
  
  # Set active assay
  DefaultAssay(integrated_seurat) <- active_assay
  
  # ---- 🏷️ Set Idents Robustly ----
  
  # Define all possible starting prefixes based on your integration method
  prefixes <- base::tolower(reduction)
  pattern <- stringr::str_remove(prefixes, "^(umap_|integrated_|sct_|pca_)*")
  
  # Construct the full cluster column name (e.g., "harmony0.8")
  cluster_cols <- colnames(integrated_seurat@meta.data) %>%
    stringr::str_subset(pattern) %>%
    stringr::str_subset(as.character(resolution))
  
  # Check if identifier column exists
  if (length(cluster_cols) == 0) {
    log_error(sample = "",
              step = "identify_markers",
              message = glue::glue("Cluster column not found in metadata for method '{pattern}' and resolution '{resolution}'."))
  } else if (length(cluster_cols) > 1) {
    log_error(sample = "",
              step = "identify_markers",
              message = glue::glue("More than 1 matching column found in metadata : '{paste(cluster_cols, collapse = ", ")}'"))
  }
  
  # Set active ident
  Idents(object = integrated_seurat) <- cluster_cols
  
  # ---- 🔍 Find ALL Markers ----
  
  all_markers <- quiet_msg(Seurat::FindAllMarkers(object = integrated_seurat,
                                        assay = active_assay,
                                        logfc.threshold = 0.25, 
                                        test.use = "wilcox",
                                        slot = "data",
                                        min.pct = 0.1,
                                        min.diff.pct = -Inf,
                                        only.pos = TRUE))
  
  if (nrow(all_markers) == 0) {
    
    log_warn(sample = "",
             step = "identify_markers",
             message = glue::glue("No markers found for resolution '{resolution}' using method '{reduction}'."))
    return(invisible(FALSE))
  }
  
  # ---- 🧹 Annotation and Filtering ----
  
  # Get annotations from ENSEMBL
  annotations_list <- get_annotations()
  
  # Compute number of matching genes for each annotation list
  matches <- sapply(annotations_list, function(x) length(intersect(x$SYMBOL, all_markers$gene)))
  
  # Pick the annotation list with the most matches
  annotations <- annotations_list[[which.max(matches)]]
  
  sig_markers <- all_markers %>%
    # Correct division by zero warning before calculating ratio
    dplyr::mutate(pct.1 = dplyr::if_else(pct.1 == 0, 0.001, pct.1),
                  pct.2 = dplyr::if_else(pct.2 == 0, 0.001, pct.2),
                  ratio = pct.1 / pct.2) %>%
    dplyr::filter(p_val_adj <= 0.05) %>%
    dplyr::left_join(y = annotations, by = c("gene" = "ENSEMBL_SYMBOL")) %>%
    dplyr::relocate(cluster, gene, CHR, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, ratio, DESCRIPTION) %>%
    dplyr::distinct(cluster, gene, avg_log2FC, pct.1, pct.2, .keep_all = TRUE)
  
  if (nrow(sig_markers) == 0) {
    log_warn(sample = "",
             step = "identify_markers",
             message = glue::glue("No markers (p_val_adj < 0.05) identified at resolution '{cluster_cols}'. Skipping saving."))
    return(invisible(FALSE))
  }
  
  # Find top 30 markers for each major cluster
  top_markers <- sig_markers %>%
    dplyr::filter(avg_log2FC >= 0.58, p_val_adj <= 0.05) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    dplyr::slice_head(n = 30) %>%
    ungroup()
  
  # ---- 📊 Create Ordered Marker Matrix for Heatmap ----
  
  mat <- sig_markers %>%
    tidyr::pivot_wider(id_cols = cluster,
                       names_from = gene,
                       values_from = avg_log2FC,
                       values_fill = 0) %>%
    tibble::column_to_rownames("cluster") %>%
    scale()   # column wise scaling so each gene has mean 0, stdev = 1
  
  # Cluster rows and columns
  row_dist <- stats::dist(x = mat, method = "euclidean")      # distance between rows based on columns
  row_clust <- hclust(row_dist)                               # clustering based on distance calculated
  col_dist <- stats::dist(x = t(mat), method = "euclidean")
  col_clust <- hclust(col_dist)
  
  # Reorder matrix
  row_order <- rownames(mat[row_clust$order,])
  col_order <- colnames(mat[,col_clust$order])
  mat <- mat[row_order, col_order]
  
  # Truncate negative values and round
  mat[mat < 0] <- 0
  mat <- mat %>%
    t() %>%               # transpose so genes are rows
    as.data.frame() %>%   
    dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 2)))
  
  # ---- 💾 Save Markers to Excel ----
  
  file_name <- file.path(output_dir, paste0("Markers.All.", cluster_cols, ".xlsx"))
  wb <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(wb, "All_Markers")
  openxlsx::writeData(wb, "All_Markers", sig_markers)
  
  openxlsx::addWorksheet(wb, "Top_Markers")
  openxlsx::writeData(wb, "Top_Markers", top_markers)
  
  openxlsx::addWorksheet(wb, "Matrix")
  openxlsx::writeData(wb, "Matrix", mat, rowNames = TRUE)
  
  openxlsx::saveWorkbook(wb, file = file_name, overwrite = TRUE)
  
  # ---- 🪵 Log Output ----
  
  log_info(sample = "",
           step = "identify_markers",
           message = glue::glue("Marker analysis complete. Results saved to : '{file_name}'."))
}

plot_umap <- function(integrated_seurat, reduction, color.col, filename, output_dir, split.col = NULL){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = integrated_seurat, reduction = reduction, 
                  metadata_cols = c(color.col, split.col), 
                  filename = filename, output_dir = output_dir)
  
  # ---- 🔄 Check and Update Reduction Name ----
  
  if (!(reduction %in% names(integrated_seurat@reductions))) {
    message("Reduction '", reduction, "' is NOT present in the Seurat object.")
    
    # Check alternative reduction name
    alt_reduction <- paste0("umap_", tolower(reduction))
    
    if (alt_reduction %in% names(integrated_seurat@reductions)) {
      log_info(sample = "",
               step = "plot_umap",
               message = glue::glue("Using alternative reduction : '{alt_reduction}'."))
      reduction <- alt_reduction
    } else {
      log_error(sample = "",
                step = "plot_umap",
                message = glue::glue("Alternative reduction '{alt_reduction}' is NOT present."))
    }
  }
  
  # ---- 👥 Determine Groups for Plotting ----
  
  if (is.null(split.col)) {
    # No splitting → single panel
    groups <- "All"
  } else if (length(split.col) == 1) {
    # Single-column split → one panel per unique value
    groups <- unique(integrated_seurat@meta.data[[split.col]])
  } else {
    # Multi-column split → one panel per column name (no subsetting)
    groups <- split.col
  }
  
  # ---- 🖼️ Generate Plots for each group ----
  
  all_plots <- list()
  for (g in groups) {
    
    if (is.null(split.col)) {
      
      # No split → use full object
      subset_obj <- integrated_seurat
      
      # Determine levels of color.col
      all_levels <- sort(unique(integrated_seurat@meta.data[[color.col]]))
      
    } else if (length(split.col) == 1) {
      
      # Single-column split → subset by value
      subset_obj <- integrated_seurat[, integrated_seurat@meta.data[[split.col]] == g]
      
      # Determine levels of color.col
      all_levels <- sort(unique(integrated_seurat@meta.data[[color.col]]))
      
    } else{
      
      # multi-column split → use full object
      subset_obj <- integrated_seurat
      
      # Determine levels of color.col
      color.col <- g
      all_levels <- sort(unique(integrated_seurat@meta.data[[color.col]]))
      
    }
    
    # Fix factor levels for consistent coloring across splits
    subset_obj@meta.data[[color.col]] <- factor(subset_obj@meta.data[[color.col]], levels = all_levels)
    
    # Create UMAP plot
    p <- quiet_msg(Seurat::DimPlot(object     = subset_obj,
                         reduction  = reduction,
                         cols       = custom_palette,
                         label      = FALSE,
                         group.by   = color.col,
                         pt.size    = 0.2,
                         label.size = 5,
                         repel      = FALSE,
                         raster     = FALSE) +
      ggplot2::labs(color = "CLUSTERS", x = "UMAP_1", y = "UMAP_2", title = g) +
      custom_theme +
      ggplot2::coord_fixed(ratio = 1))  # 1 unit on x-axis = 1 unit on y-axis, ensuring a square aspect
    
    all_plots[[g]] <- p
    
  }
  
  # ---- 🌐 Combine Plots Using cowplot ----
  
  n_plots <- length(all_plots)
  ncol_plots <- ceiling(sqrt(n_plots))
  nrow_plots <- ceiling(n_plots / ncol_plots)
  
  # Combine all plots
  combined_plot <- cowplot::plot_grid(plotlist = all_plots, 
                                      ncol = ncol_plots, 
                                      nrow = nrow_plots)
  
  # ---- 💾 Save Combined Plot ----
  
  file_name <- file.path(output_dir, paste0(filename, ".pdf"))
  ggplot2::ggsave(filename  = file_name,
                  plot      = combined_plot,
                  width     = ncol_plots * 10, # extra 2 inch for legend
                  height    = nrow_plots * 8, 
                  units     = "in",
                  limitsize = FALSE,
                  bg        = "white")
  
  # ---- 🪵 Log Output ----
  
  log_info(sample = "",
           step = "plot_umap",
           message = glue::glue("UMAP plot saved successfully to : '{file_name}'."))
  
}

plot_features <- function(integrated_seurat, features, reduction, filename, output_dir, split.col=NULL){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = integrated_seurat, reduction = reduction, 
                  features = features, metadata_cols = c(split.col),
                  filename = filename, output_dir = output_dir)
  
  # ---- 🔄 Check and Update Reduction Name ----
  
  if (!(reduction %in% names(integrated_seurat@reductions))) {
    message("Reduction '", reduction, "' is NOT present in the Seurat object.")
    
    # Check alternative reduction name
    alt_reduction <- paste0("umap_", tolower(reduction))
    
    if (alt_reduction %in% names(integrated_seurat@reductions)) {
      log_info(sample = "",
               step = "plot_features",
               message = glue::glue("using alternative reduction : '{alt_reduction}'."))
      reduction <- alt_reduction
    } else {
      log_error(sample = "",
                step = "plot_features",
                message = glue::glue("Alternative reduction '{alt_reduction}' is NOT present."))
    }
  }
  
  # ---- 👥 Determine Groups for Plotting ----
  
  if (is.null(split.col)) {
    groups <- "All"
  } else {
    groups <- unique(integrated_seurat@meta.data[[split.col]])
  }
  
  # ---- 🖼️ Generate Plots for each group ----
  
  all_plots <- list()
  for (feature in features) {
    
    # Determine feature type
    is_gene <- feature %in% rownames(SeuratObject::GetAssayData(integrated_seurat, assay = "RNA", slot = "data"))
    is_metadata <- feature %in% colnames(integrated_seurat@meta.data)
    
    if (!is_gene & !is_metadata) next  # Skip missing feature
    
    for (g in groups) {
      
      # Subset object if splitting
      if (g == "All") {
        subset_obj <- integrated_seurat
        group_label <- ""
      } else {
        subset_obj <- integrated_seurat[, integrated_seurat@meta.data[[split.col]] == g]
        group_label <- g
      }
      
      # Determine expression range
      if (is_gene) {
        expr_data <- SeuratObject::GetAssayData(subset_obj, assay = "RNA", slot = "data")
        max_expr <- max(expr_data[feature, ], na.rm = TRUE)
        min_expr <- 0
      } else if (is_metadata) {
        max_expr <- max(subset_obj@meta.data[[feature]])
        min_expr <- min(subset_obj@meta.data[[feature]])
      }
      
      # Color scale
      cols <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
      
      # Create FeaturePlot
      p <- quiet_msg(Seurat::FeaturePlot(object = subset_obj,
                               features = feature,
                               slot = "data",
                               reduction = reduction,
                               #cols = c("grey", viridis(n = 10, option = "C", direction = -1)),
                               #cols =  c("#440154FF", viridis(n = 10, option = "C", direction = 1)),
                               pt.size = 0.2,
                               label = FALSE,
                               order = TRUE,
                               raster = FALSE,
                               combine = TRUE) +
        ggplot2::labs(title = paste0(feature, group_label), x = "UMAP_1", y = "UMAP_2") +
        custom_theme +
        {
          if (min_expr < 0) {
            # Case 1: Diverging Scale (Expression includes negative values)
            # Use all colors, centered at 0
            ggplot2::scale_colour_gradientn(
              colours = cols,
              values = scales::rescale(c(min_expr, 0, max_expr)),
              limits = c(min_expr, max_expr))
          } else {
            # Case 2: Sequential Scale (All expression is >= 0)
            # Start the scale at 0 (using the center color: cols[6])
            ggplot2::scale_colour_gradientn(
              colours = cols[6:11],
              values = scales::rescale(c(0, max_expr)),
              limits = c(0, max_expr)) # Start limits at 0 to match scale
          }
        })
      
      # Store plot
      all_plots[[paste0(feature, group_label)]] <- p
    }
  }
  
  # ---- 🌐 Combine Plots Using cowplot ----
  
  n_plots <- length(all_plots)
  ncol_plots <- ceiling(sqrt(n_plots))
  nrow_plots <- ceiling(n_plots / ncol_plots)
  
  # Restrict unsupported combination
  if (ncol_plots > 10 && nrow_plots > 10) {
    
    log_error(sample = "",
              step = "plot_features",
              message = "Image size too large. More than 100 plots cannot be viewed in a single figure")
  }
  
  # Combine all plots
  combined_plot <- cowplot::plot_grid(plotlist = all_plots, 
                                      ncol = ncol_plots, 
                                      nrow = nrow_plots)
  
  # ---- 💾 Save Combined Plot ----
  
  file_name <- file.path(output_dir, paste0(filename, ".pdf"))
  ggplot2::ggsave(filename  = file_name,
                  plot      = combined_plot,
                  width     = ncol_plots * 6,
                  height    = nrow_plots * 6,
                  units     = "in",
                  limitsize = FALSE,
                  bg        = "white")
  
  # ---- 🪵 Log Output ----
  
  log_info(sample = "",
           step = "plot_features",
           message = glue::glue("Feature plots saved to : '{file_name}'."))
}

plot_metrics_post_integration <- function(integrated_seurat, output_dir){
  
  # ---- ⚙️ Validate Input Parameters ----
  
  validate_inputs(seurat_object = integrated_seurat, output_dir = output_dir)
  
  # ---- 1️⃣ Feature Plots for QC Metrics ----
  
  qc_features <- c("nUMIs", "nGenes", "S.Score", "G2M.Score", "CC.Score", "MitoRatio")
  plot_features(integrated_seurat, features = qc_features, reduction = "umap_harmony",
                filename = "UMAP.Numerical.Metrics", output_dir = output_dir, split.col = NULL)
  
  # ---- 2️⃣ UMAP Plots for Cluster Visualization ----
  
  plot_umap(integrated_seurat, reduction = "sct_pca",       color.col = "harmony0.8", filename = "Pre.Integ.PCA",  output_dir, split.col = "Sample")
  plot_umap(integrated_seurat, reduction = "integ_harmony", color.col = "harmony0.8", filename = "Post.Integ.PCA", output_dir, split.col = "Sample")
  plot_umap(integrated_seurat, reduction = "umap_harmony",  color.col = "harmony0.8", filename = "UMAP.Sample",    output_dir, split.col = "Sample")
  plot_umap(integrated_seurat, reduction = "umap_harmony",  color.col = "harmony0.8", filename = "UMAP.Phase",     output_dir, split.col = "Phase")
  plot_umap(integrated_seurat, reduction = "umap_cca",      color.col = NULL,         filename = "UMAP.CCA",       output_dir, split.col = paste0("cca",      seq(0.2, 1.2, by = 0.2)))
  plot_umap(integrated_seurat, reduction = "umap_rpca",     color.col = NULL,         filename = "UMAP.RPCA",      output_dir, split.col = paste0("rpca",     seq(0.2, 1.2, by = 0.2)))
  plot_umap(integrated_seurat, reduction = "umap_jointpca", color.col = NULL,         filename = "UMAP.JointPCA",  output_dir, split.col = paste0("jointpca", seq(0.2, 1.2, by = 0.2)))
  plot_umap(integrated_seurat, reduction = "umap_harmony",  color.col = NULL,         filename = "UMAP.Harmony",   output_dir, split.col = paste0("harmony",  seq(0.2, 1.2, by = 0.2)))
  plot_umap(integrated_seurat, reduction = "umap_harmony",  color.col = NULL,         filename = "UMAP.Cell.QC",   output_dir, split.col = c("DropletUtils", "CellRanger", "Quality",
                                                                                                                                             "DoubletFinder", "scDblFinder", "QC"))
  # ---- 🪵 Log Output ----
  
  log_info(sample = "",
          step = "plot_metrics_post_integration",
          message = "All post-integration QC and UMAP plots generated successfully.")
}




cluster_by_umap_dist <- function(integrated_seurat,  reduction, resolution, 
                                 percentile = 0.25, connectedness_thresh = 1) {
  
  # Check original reduction
  if (!(reduction %in% names(integrated_seurat@reductions))) {
    message("Reduction '", reduction, "' is NOT present in the Seurat object.")
    
    # Check alternative reduction name
    alt_reduction <- paste0("umap_", tolower(reduction))
    
    if (alt_reduction %in% names(integrated_seurat@reductions)) {
      message("Alternative reduction '", alt_reduction, "' is present. Using this instead.")
      reduction <- alt_reduction
    } else {
      stop("Alternative reduction '", alt_reduction, "' is NOT present.")
    }
  }
  
  # ---- Set Idents Robustly ----
  
  # Extract method name from reduction (e.g., "integrated_harmony" -> "harmony").
  # This assumes the clustering columns were named as 'harmony0.8', 'cca0.6', etc.
  method_name <- base::tolower(reduction)
  method_name <- stringr::str_remove(method_name, "^(umap_|integrated_|sct_|pca_)*")
  
  # Construct the full cluster column name (e.g., "harmony0.8")
  idents <- paste0(method_name, resolution)
  
  # Check if identifier column exists
  if (!idents %in% colnames(integrated_seurat@meta.data)) {
    stop(paste0("Cluster identifier column '", idents, "' not found in metadata. Ensure clustering was run and named correctly."))
  }
  
  # 1️⃣ Extract UMAP coordinates
  umap <- SeuratObject::Embeddings(object = integrated_seurat, reduction = reduction)
  
  # 2️⃣ Extract clusters: if resolution provided, use that metadata column
  clusters <- integrated_seurat@meta.data[[idents]]
  
  # 3️⃣ Build df
  df <- data.frame(UMAP_1 = umap[,1],
                   UMAP_2 = umap[,2],
                   cluster = clusters)
  
  # 4️⃣ Compute cluster centroids
  centroids <- df %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(across(1:2, mean))
  
  # 5️⃣ Compute distance matrix
  dist_matrix <- as.matrix(dist(centroids[, c("UMAP_1", "UMAP_2")]))
  rownames(dist_matrix) <- centroids$cluster
  colnames(dist_matrix) <- centroids$cluster
  
  # 6️⃣ Compute connectedness
  connectedness <- base::colSums(dist_matrix < stats::quantile(dist_matrix, percentile))
  
  # 7️⃣ Identify distinct clusters first
  distinct_clusters <- names(connectedness[connectedness == connectedness_thresh])
  
  # Identify big clusters
  big_clusters <- names(connectedness[connectedness > connectedness_thresh])
  
  centroids <- df %>%
    dplyr::filter(cluster %in% big_clusters) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(across(1:2, mean))
  
  
  # 8️⃣ Assign dist_cluster
  integrated_seurat$dist_cluster <- ifelse(
    clusters %in% main_cluster,
    "BigCluster",
    paste0("Outlier_", clusters)
  )
  
  # 9️⃣ Optional plot
  if(plot){
    DimPlot(integrated_seurat, group.by = "dist_cluster", label = TRUE)
  }
  
  return(integrated_seurat)
  
  
  
  
  
  
  
}

calc_module_scores <- function(integrated_seurat, assay, proj.params, output_dir){
  
  set.seed(1234)
  
  # Read marker file
  marker_df <- read.xlsx(xlsxFile = proj.params$cell.type.marker.file)
  
  # Set default assay
  DefaultAssay(integrated_seurat) <- assay
  
  # Iterate through each celltype and plot its module scores
  for (i in 1:ncol(marker_df)){
    
    features <- marker_df[[i]] %>% na.omit()
    features <- rownames(integrated_seurat@assays[[assay]]$data)[tolower(rownames(integrated_seurat@assays[[assay]]$data)) %in% 
                                                                   tolower(features)]
    features <- list(sort(features))
    
    # Calculate module scores
    if (lengths(features) > 1){
      integrated_seurat <- Seurat::AddModuleScore(object=integrated_seurat,
                                                  features=features,
                                                  assay=assay,
                                                  slot="data",
                                                  name=make.names(colnames(marker_df)[i]))
      
      # names(features) <- make.names(colnames(marker_df)[i])
      # integrated_seurat <- UCell::AddModuleScore_UCell(obj=integrated_seurat,
      #                                                  features=features,
      #                                                  assay="RNA",
      #                                                  slot="data",
      #                                                  name="_UCell")
    }
  }
  return(integrated_seurat)
}

plot_module_scores <- function(integrated_seurat, assay, resolution, reduction, proj.params, output_dir){
  
  set.seed(1234)
  
  # Read marker file
  marker_df <- read.xlsx(xlsxFile = proj.params$cell.type.marker.file)
  
  modules <- intersect(colnames(integrated_seurat@meta.data), paste0(make.names(colnames(marker_df)),1))
  
  # Change active.ident
  idents <- paste0("cluster.", resolution, ".", base::tolower(reduction))
  Idents(object = integrated_seurat) <- idents
  
  # Set default assay
  DefaultAssay(integrated_seurat) <- assay
  
  module_score_seurat <- function(module){
    
    Seurat::FeaturePlot(object = integrated_seurat,
                        slot = "data",
                        features = module,
                        reduction = paste0("umap.", base::tolower(reduction)),
                        #cols= c("grey", viridis(n = 10, option="C", direction = 1)),
                        pt.size = 0.2,
                        order = TRUE,
                        label = TRUE,
                        raster = FALSE,
                        combine = TRUE) +  
      #BUG: if raster = TRUE, order = TRUE is ignored. So, set raster = FALSE
      scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name="RdBu"))[5:11])
    # Dont set limits like in plot_features as module scores can be negative too
  }
  
  combined_plot <- purrr::map(.x = modules, 
                              .f = module_score_seurat) %>%         #get(funcs[j]))
    cowplot::plot_grid(plotlist=.,
                       align="hv",
                       axis="tblr",
                       nrow = NULL,
                       ncol = dplyr::if_else(ncol(marker_df) > 10, 5, ceiling(ncol(marker_df)/3)),
                       rel_widths = 1,
                       rel_heights = 1,
                       greedy = TRUE,
                       byrow = TRUE)
  
  ggplot2::ggsave(filename = "Module_plot.jpg",
                  plot = combined_plot,
                  device = "jpeg",
                  path = output_dir,
                  width = 8.5*4,
                  height = 11*2,
                  units = c("in"),
                  dpi = 300,
                  limitsize = FALSE,
                  bg = "white")
}

annotate_sc_sp <- function(integrated_seurat, assay, resolution, reduction, output_dir){
  
  set.seed(1234)
  
  # Get columns containing module scores fro each cell type
  score_cols <- integrated_seurat@meta.data %>% 
    dplyr::select(ends_with("1")) %>% 
    colnames()
  
  # Extract numeric score matrix
  score_mat <- as.matrix(integrated_seurat@meta.data[, score_cols])
  
  # Get index of max per row
  max_idx <- max.col(score_mat, ties.method = "first")
  
  # Map to column names
  integrated_seurat@meta.data$Est.Cell.Type <- colnames(score_mat)[max_idx]
  
  # Identify col to be used for annotation
  cluster.col <- paste0("cluster.", resolution, ".", base::tolower(reduction))
  
  # ---- Cluster-level summary ----
  cluster_summary <- integrated_seurat@meta.data %>%
    dplyr::count(cluster = .data[[cluster.col]], Est.Cell.Type, name = "Cell.Count") %>%
    dplyr::arrange(cluster, desc(Cell.Count)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(
      Total.Cells = sum(Cell.Count),
      Primary.Cell.Count = Cell.Count[1],
      Secondary.Cell.Count = if_else(length(Cell.Count) > 1, Cell.Count[2], NA_integer_),
      Primary.Cell.Percent  = Primary.Cell.Count / Total.Cells * 100,
      Secondary.Cell.Percent  = Secondary.Cell.Count / Total.Cells * 100,
      Primary.Cell.Type = Est.Cell.Type[1],
      Secondary.Cell.Type = if_else(length(Cell.Count) > 1, Est.Cell.Type[2], NA_character_),
      Difference.Percent = Primary.Cell.Percent - Secondary.Cell.Percent,
      .groups = "drop") %>%
    as.data.frame() %>%
    dplyr::mutate(Confidence = case_when(Primary.Cell.Percent >= 80 ~ "HIGH",
                                         (Primary.Cell.Percent >= 60 & Secondary.Cell.Percent < 30 )~ "MEDIUM",
                                         TRUE ~ "LOW"))
  
  cluster_info <- cluster_summary %>% 
    dplyr::select(cluster, Primary.Cell.Type, Secondary.Cell.Type, Primary.Cell.Percent, Secondary.Cell.Percent, Confidence)
  
  integrated_seurat@meta.data <- integrated_seurat@meta.data %>% 
    dplyr::left_join(cluster_info, by = setNames("cluster", cluster.col))
  rownames(integrated_seurat@meta.data) <- integrated_seurat@meta.data$Cell
  
  # Create .rds object for integrated seurat object
  if (assay %in% c("RNA", "SCT")){
    saveRDS(integrated_seurat, file=file.path(output_dir, "integrated_seurat.ann.rds"))
  } else{
    # For Spatial.008um and Spatial.016um assays
    saveRDS(integrated_seurat, file=file.path(output_dir, paste0("integrated_seurat.", assay, ".ann.rds")))
  }
  
  return(integrated_seurat)
}

### Annotate based on clusters variable defined by user
# clusters <- list("Hepatocytes"         = c(),
#                  "Pancreatic.Acinar"   = c(),
#                  "Pancreatic.Islet"    = c(),
#                  "B.Plasma"            = c(),
#                  "T.NK"                = c(),
#                  "Fibroblasts"         = c(),
#                  "Macrophages"         = c(),
#                  "Dendritic"           = c(),
#                  "Endothelial"         = c(),
#                  "Lymph.Endothelial"   = c(),
#                  "Myocytes"            = c(),
#                  "CAFs"                = c(),
#                  "Epithelial"          = c(),
#                  "Neurons"             = c(),
#                  "Epithelial.I"        = c(),
#                  "Epithelial.II"       = c(),
#                  "Epithelial.III"      = c(),
#                  "Epithelial.IV"       = c(),
#                  "Unclassified"        = c())

annotate_manual_sc_sp <- function(integrated_seurat, assay, clusters, resolution, reduction, output_dir){
  
  set.seed(1234)
  
  idents <- paste0("cluster.", resolution, ".", base::tolower(reduction))
  
  # Make sure you have assigned all clusters to one of the cell types
  # NOTE: "integrated_snn_res.1.4" etc are stored as factors. 
  # So, use as.character() and then as.numeric() to get accurate cluster values
  list_1 <- integrated_seurat@meta.data %>% 
    dplyr::count(get(idents)) %>% 
    dplyr::select(identity(1)) %>% 
    unlist(use.names=FALSE) %>% 
    as.character() %>% 
    as.numeric() %>% 
    sort()
  
  list_2 <- clusters %>% 
    unlist(., use.names=FALSE) %>% 
    sort()
  
  # Proceed with annotation ONLY if all clusters have been annotated
  if (identical(list_1, list_2)){
    print("All Clusters have been annotated")
    
    # Extract metadata from Seurat object, assign appropriate resolution to
    # seurat_clusters column and add Cell.Type, Cell.Subtype columns
    data <- integrated_seurat@meta.data %>% 
      dplyr::mutate(seurat_clusters=get(idents),
                    Cell.Type=NA, Cell.Subtype=NA)
    
    # Assign cell type based on cluster numbers within seurat_clusters column
    for (j in 1:nrow(data)){
      for (i in 1:length(clusters)){
        if (as.numeric(as.character(data$seurat_clusters[j])) %in% clusters[[i]]){
          data$Cell.Type[j] <- names(clusters[i])
        }
      }
    }
    
    # Check summary of cell counts
    print(data %>% dplyr::count(Cell.Type) %>% dplyr::arrange(n))
    cat("\n")
    
    # Import the metadata into Seurat object
    integrated_seurat@meta.data <- data
    
    # Create .rds object for integrated seurat object
    if (assay == "RNA"){
      saveRDS(integrated_seurat, file=file.path(output_dir, "integrated_seurat.ann.rds"))
    } else{
      # For Spatial.008um and Spatial.016um assays
      saveRDS(integrated_seurat, file=file.path(output_dir, paste0("integrated_seurat.", assay, ".ann.rds")))
    }
    
  } else {
    cat("\nYou missed annotating these clusters:\t", setdiff(list_1, list_2))
    cat("\nThese clusters are not present in data:\t", setdiff(list_2, list_1))
    cat("\nThese clusters have duplicate annotation:\t", list_2[duplicated(list_2)])
  }
  
  return(integrated_seurat)
} 

# Use to check how good your markers are able to calc module score for each cell
# and how accurately they align with known cell annotations
estimate_cell_type <- function(integrated_seurat){
  
  set.seed(1234)
  
  # Get columns containing module scores fro each cell type
  score_cols <- integrated_seurat@meta.data %>% 
    dplyr::select(ends_with("1")) %>% 
    colnames()
  
  # Extract numeric score matrix
  score_mat <- as.matrix(integrated_seurat@meta.data[, score_cols])
  
  # Get index of max per row
  max_idx <- max.col(score_mat, ties.method = "first")
  
  # Map to column names
  integrated_seurat@meta.data$Est.Cell.Type <- colnames(score_mat)[max_idx]
  
  cleaned_score_cols <- gsub("1$", "", score_cols)
  cell.types <- unique(integrated_seurat@meta.data$Cell.Type)
  
  # Key:value dataframe
  map_df <- list("Epithelial1" = c("Epithelial.I", "Epithelial.II"),
                 "Macrophages1" = c("Macrophages"),
                 "Myofibroblasts1" = c("Myocytes.Myofibroblasts"),
                 "Endothelial1" = c("Endothelial"),
                 "Lymphatic.Endothelial1" = c("Endothelial"),
                 "T1" = c("T.NK"), 
                 "NK1" = c("T.NK"),
                 "Fibroblasts1" = c("Fibroblasts"),
                 "B1" = c("B.Plasma"),
                 "Plasma1" = c("B.Plasma"),
                 "Granulocytes1" = c("Granulocytes", "Mast"),
                 "Erythrocytes1" = c("Erythrocytes"),
                 "Neurons1" = c("Neurons"),
                 "Dendritic1" = c("Dendritic"),
                 "Mast1" = c("Mast"))
  
  
  integrated_seurat@meta.data <- integrated_seurat@meta.data %>%
    dplyr::mutate(Confidence = ifelse(mapply(function(est, true){
      true %in% map_df[[est]]}, Est.Cell.Type, Cell.Type), "HIGH", "LOW"))
  
  
  return(integrated_seurat) 
}

plot_dot_plot <- function(integrated_seurat, idents, features, filename, output_dir, gene_y_axis = FALSE, split.col = NULL){
  
  set.seed(1234)
  
  # Re-order the active ident alphabetically
  Idents(integrated_seurat) <- idents
  Idents(integrated_seurat) <- base::factor(x = integrated_seurat@active.ident, 
                                            levels= sort(levels(integrated_seurat@active.ident)))
  
  # Determine grouping variable
  if (is.null(split.col)) {
    groups <- "All"
  } else {
    groups <- unique(integrated_seurat@meta.data[[split.col]])
  }
  
  # Loop over groups
  plot_list <- list()
  for (g in groups) {
    if (g == "All") {
      subset_obj <- integrated_seurat
      plot_title <- "All"
    } else {
      subset_obj <- integrated_seurat[, integrated_seurat@meta.data[[split.col]] == g]
      plot_title <- g
    }
    
    p <- Seurat::DotPlot(object = subset_obj,
                         assay = "RNA",
                         features = features,
                         dot.min = 0,
                         dot.scale = 4,
                         scale = TRUE,         # if plotting 1+ gene, always SCALE
                         scale.by = "size",    # scale percent cells expressed i.e. size(area) of dots
                         scale.min = 0,        # set this to 0, so pct.size = 0 has smallest size
                         scale.max = NA) +     # if you set this to 100, points become too small
      ggplot2::geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.25) +
      ggplot2::labs(fill = "Avg. Expression",         # rename color legend
                    size  = "% Cells Expressed",
                    title = plot_title) +    # rename size legend   
      ggplot2::scale_colour_distiller(type = "div", palette = "RdYlGn", direction = -1) +
      ggplot2::guides(size = ggplot2::guide_legend(
        override.aes = list(shape = 21, colour = "black", fill = "white", stroke = 0.75))) +
      custom_theme +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    if (gene_y_axis){
      p <- p + ggplot2::coord_flip()
    }
    
    plot_list[[plot_title]] <- p
  }
  
  # Per panel dimensions
  n_y <- length(unique(integrated_seurat@active.ident))  
  height_per_element <- 0.5
  width_per_element <- 0.2 
  if (gene_y_axis){
    height_panel <- length(features) * width_per_element + 2
    width_panel <- n_y * height_per_element + 6
  } else{
    height_panel <- n_y * height_per_element + 2              # extra for x-axis labels
    width_panel <- length(features) * width_per_element + 6 #extra for y-axis labels, legend
  }
  
  # Automatically calculate ncol and nrow for cowplot
  n_plots <- length(plot_list)
  ncol_plots <- ceiling(sqrt(n_plots))
  nrow_plots <- ceiling(n_plots / ncol_plots)
  
  # Total height
  total_height <- height_panel * nrow_plots
  total_width  <- width_panel * ncol_plots
  
  # Combine all plots
  combined_plot <- cowplot::plot_grid(plotlist = plot_list, 
                                      ncol = ncol_plots, 
                                      nrow = nrow_plots)
  
  # Save combined plot
  ggsave(filename = file.path(output_dir, paste0(filename,".pdf")),
         plot = combined_plot,
         width = total_width,
         height = total_height,
         limitsize = FALSE,
         units = "in",
         bg        = "white")
}

plot_dot_custom <- function(integrated_seurat, assay, proj.params, ident.1, ident.2, features, filename, output_dir, CellType.col = "Cell.Type", split.col=NULL){
  
  set.seed(1234)
  
  # Keep only cell types present in your dataset and marker file
  marker_df <- openxlsx::read.xlsx(file.path(proj.params$cell.type.marker.file))
  celltypes_in_seurat <- unique(integrated_seurat[[CellType.col]])
  celltypes_to_keep <- intersect(colnames(marker_df), celltypes_in_seurat)
  
  # Determine groups based on split.col
  groups <- if (!is.null(split.col) && split.col %in% colnames(integrated_seurat@meta.data)) {
    as.character(unique(integrated_seurat@meta.data[[split.col]]))
  } else {
    "All"
  }
  
  # Determine features present in data set
  features <- features[base::tolower(features) %in% base::tolower(SeuratObject::Features(integrated_seurat))] %>% 
    na.omit() %>%
    as.vector()
  
  if(length(features) == 0) stop("No matching features found in the assay.")
  
  # Loop over groups
  plot_list <- list()
  for (g in groups) {
    
    # Subset object
    if (g == "All") {
      subset_obj <- integrated_seurat
      plot_title <- ""
    } else {
      subset_obj <- integrated_seurat[, integrated_seurat@meta.data[[split.col]] == g]
      plot_title <- g
    }
    
    # Determine levels based on ident.1
    levels <- if (!is.null(ident.1) && ident.1 %in% colnames(integrated_seurat@meta.data)) {
      subset_obj@meta.data[[ident.1]] %>% unique() %>% as.character()
    } else {
      stop(ident.1, " missing in metadata")
    }
    
    # Determine levels.extra based on ident.2
    levels.extra <- if (!is.null(ident.2) && ident.2 %in% colnames(integrated_seurat@meta.data)) {
      subset_obj@meta.data[[ident.2]] %>% unique() %>% as.character()
    } else {
      levels.extra <- NA # set NA so loop runs once later below
    }
    
    # Limit to one feature if ident.2 is used
    if (length(features) > 1 & !is.null(ident.2)){
      warning("You cannot plot more than 1 feature using variables ", ident.1,  " and ", ident.2, ".\n")
      warning("Proceeding to generate dot plot using only one feature ", features[[1]], ".")
      features <- features[1]
    }
    
    # Prepare dot plot data
    pct.exp <- list()
    avg.exp <- list()
    features.plot <- list()
    id <- list()
    id.extra <- list()
    
    for (le in levels.extra){
      for (l in levels){
        for (f in features){
          
          # Cells to retain
          cells <- subset_obj@meta.data %>% 
            dplyr::filter(.data[[ident.1]] == l) 
          if (!is.null(ident.2)){
            cells <- cells %>% 
              dplyr::filter(.data[[ident.2]] == le) 
          }
          cells <- rownames(cells)
          
          # Subset features and expression data
          expr_data <- Seurat::GetAssayData(object = subset_obj, assay = assay, layer = "data")
          subset_data <- expr_data[f, cells, drop = FALSE]
          pct <- sum(subset_data > 0)/length(subset_data) * 100
          avg <- mean(expm1(subset_data), na.rm = TRUE)   # back-transform from log1p
          
          pct.exp <- c(pct.exp, pct)
          avg.exp <- c(avg.exp, avg)
          features.plot <- c(features.plot, f)
          id <- c(id, l)
          id.extra <- c(id.extra, le)
        }
      }
    }
    
    # Combine and scale data
    dotplot_data <- data.frame(avg.exp = unlist(avg.exp),
                               pct.exp = unlist(pct.exp),
                               features.plot = unlist(features.plot),
                               id = unlist(id),
                               id.extra = unlist(id.extra)) %>%
      dplyr::group_by(features.plot, id.extra) %>%
      dplyr::mutate(avg.exp.scaled = scale(avg.exp)) %>%
      dplyr::mutate(avg.exp.scaled = dplyr::case_when(avg.exp.scaled > 2.5 ~ 2.5,
                                                      avg.exp.scaled < -2.5 ~ -2.5,
                                                      TRUE ~ avg.exp.scaled)) %>%
      as.data.frame() %>%
      dplyr::mutate(features.plot = factor(features.plot, levels = features))
    
    if (ident.1 == "Cell.Type"){
      dotplot_data <- dotplot_data %>%
        filter(id %in% celltypes_to_keep)
    }
    
    # Determine X and Y variables
    if (is.null(ident.2) && length(features) > length(levels) & length(features) > length(levels.extra)){
      y_var <- "features.plot"
      x_var <- "id"
    } else if (!is.null(ident.2) && length(levels.extra) > length(levels)){
      y_var <- "id.extra"
      x_var <- "id"
    } else if (!is.null(ident.2) && length(levels.extra) < length(levels)){
      y_var <- "id"
      x_var <- "id.extra"
    } else{
      y_var <- "id"
      x_var <- "features.plot"
    }
    
    # Calculate panel dimensions
    height_per_element <- 0.5
    width_per_element <- 0.5
    height_panel <- length(unique(dotplot_data[[y_var]])) * height_per_element + 2
    width_panel <- length(unique(dotplot_data[[x_var]])) * width_per_element + 6
    x_label_size <- 72 * width_per_element * 0.5
    y_label_size <- 72 * height_per_element * 0.5  # 0.5 scaling so it doesnt get too large
   
    # Create ggplot
    p <- ggplot(data = dotplot_data, 
                mapping = aes(x = !!sym(x_var), y = !!sym(y_var), 
                              size = pct.exp, fill = avg.exp.scaled)) +
      geom_point(shape = 21, colour = "black", stroke = 0.25) +
      labs(x = "", y = "", title = plot_title,
           fill = "Avg. Expression",
           size = "% Cells Expressed") +
      scale_fill_distiller(type = "div", palette = "RdYlGn", direction = -1) +
      guides(size = guide_legend(override.aes = list(shape = 21, colour = "black", fill = "white", stroke = 0.75))) +
      scale_size_area(max_size = 15) + #min_size = 1) +
      custom_theme +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = x_label_size ),
            axis.text.y = element_text(size = y_label_size),
            legend.text  = element_text(size = 0.8 * y_label_size),   # slightly smaller than y-axis labels
            legend.title = element_text(size = y_label_size))
    
    plot_list[[plot_title]] <- p
  }
  
  # Determine layout
  n_plots <- length(plot_list)
  ncol_plots <- ceiling(sqrt(n_plots))
  nrow_plots <- ceiling(n_plots / ncol_plots)
  
  # Total dimensions
  total_height <- height_panel * nrow_plots
  total_width  <- width_panel * ncol_plots
  
  # Combine plots
  combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = ncol_plots, nrow = nrow_plots)  
  
  # Save plot
  filename <- gsub("[/\\?<>\\:*|\"]", "_", filename)
  ggsave(filename = file.path(output_dir, paste0(filename,".pdf")),
         plot = combined_plot,
         width = total_width,
         height = total_height,
         limitsize = FALSE,
         units = "in",
         bg = "white")
}

find_top_features <- function(integrated_seurat, assay, proj.params){
  
  set.seed(1234)
  
  marker_df <- openxlsx::read.xlsx(file.path(proj.params$cell.type.marker.file))
  expr_data <- Seurat::GetAssayData(object = integrated_seurat, assay = assay, layer = "data")
  
  top_features <- lapply(marker_df, FUN = function(f){
    
    # Subset features and expression data
    f <- f[base::tolower(f) %in% base::tolower(SeuratObject::Features(integrated_seurat))] %>% 
      na.omit() %>%
      as.vector()
    
    if(length(f) == 0) return(NULL)
    
    expr_sub <- expr_data[f, , drop = FALSE]
    
    # Find top 3 expressing markers
    rowMeans(expr_sub, na.rm = TRUE) %>% 
      sort(., decreasing = TRUE) %>% 
      head(3) %>% 
      names()
  })
  
  # Keep only cell types present in your dataset
  celltypes_in_seurat <- unique(integrated_seurat$Cell.Type)
  top_features <- top_features[intersect(names(top_features), celltypes_in_seurat)]
  
  return(unlist(top_features, use.names = FALSE))
}


# Input is seurat object of a single slide with columns X, Y, Sample, Group
plot_spatial_map <- function(plot.seurat, x1, y1, x2, y2, suffix, output_dir){
  
  set.seed(1234)
  
  Sample <- unique(plot.seurat@meta.data$Sample)
  
  # Get an idea of co-ordinates
  xmin <- (min(plot.seurat@meta.data$X)%/%100-1)*100
  xmax <- (max(plot.seurat@meta.data$X)%/%100+1)*100
  ymin <- (min(plot.seurat@meta.data$Y)%/%100-1)*100
  ymax <- (max(plot.seurat@meta.data$Y)%/%100+1)*100
  
  if (Sample == "TMA1-A1"){
    # Flip on X axis to get correct orientation for TMA1-A1
    p1 <- ggplot(plot.seurat@meta.data, mapping=aes(x=X, y=Y, color = Sample)) +
      geom_point() +
      scale_x_reverse(limits =c(xmax, xmin), breaks=seq(xmax, xmin, -100), position="top") +
      scale_y_continuous(limits =c(ymin, ymax), breaks=seq(ymin, ymax, 100), position="left") +
      geom_vline(xintercept = x1) +
      geom_hline(yintercept = y1) +
      scale_color_manual(values=custom_palette)
    
    p2 <- ggplot(plot.seurat@meta.data, mapping=aes(x=X, y=Y, color = Treatment)) +
      geom_point() +
      scale_x_reverse(limits =c(xmax, xmin), breaks=seq(xmax, xmin, -100), position="top") +
      scale_y_continuous(limits =c(ymin, ymax), breaks=seq(ymin, ymax, 100), position="left") +
      geom_vline(xintercept = x1) +
      geom_hline(yintercept = y1) +
      scale_color_manual(values=custom_palette)
  } else if (Sample == "TMA1-D1"){
    # Flip on Y axis to get correct orientation for TMA1-D1
    p1 <- ggplot(plot.seurat@meta.data, mapping=aes(x=X, y=Y, color = Sample)) +
      geom_point() +
      scale_x_continuous(limits =c(xmin, xmax), breaks=seq(xmin, xmax, 100), position="top") +
      scale_y_reverse(limits =c(ymax, ymin), breaks=seq(ymax, ymin, -100), position="left") +
      geom_vline(xintercept = x2) +
      geom_hline(yintercept = y2) +
      scale_color_manual(values=custom_palette)
    
    p2 <- ggplot(plot.seurat@meta.data, mapping=aes(x=X, y=Y, color = Treatment)) +
      geom_point() +
      scale_x_continuous(limits =c(xmin, xmax), breaks=seq(xmin, xmax, 100), position="top") +
      scale_y_reverse(limits =c(ymax, ymin), breaks=seq(ymax, ymin, -100), position="left") +
      geom_vline(xintercept = x2) +
      geom_hline(yintercept = y2) +
      scale_color_manual(values=custom_palette)
  }
  
  # Save the plot
  ggplot2::ggsave(filename=paste0("Sample.Map.", Sample, ".", suffix, ".tiff"),
                  plot=p1+p2,
                  device="jpeg",
                  path=output_dir,
                  width=25,
                  height=8.5,
                  units=c("in"),
                  dpi=600,
                  limitsize=TRUE,
                  bg="white")
}  

tabulate_frequency <- function(integrated_seurat, split.cols, output_dir){
  
  set.seed(1234)
  
  # Create a new workbook
  wb <- createWorkbook()
  
  for (split.col in split.cols){
    if(!is.null(integrated_seurat@meta.data[[split.col]] %>% unique())){
      
      counts <- integrated_seurat@meta.data %>% 
        dplyr::count(.data[[split.col]], Cell.Type) %>% 
        dplyr::group_by(.data[[split.col]]) %>% 
        dplyr::mutate(Percent = round(100 * n / sum(n), 2)) %>%
        tidyr::pivot_wider(names_from = .data[[split.col]], values_from = c(n, Percent)) %>%
        as.data.frame()
      
      
      # Make sheet name safe
      sheet_name <- substr(gsub("[/\\?<>\\:*|\"]", "_", split.col), 1, 31)  # Excel sheet name limit = 31
      addWorksheet(wb, sheetName = sheet_name)
      writeData(wb, sheet = sheet_name, counts)
    }
  }
  
  # Save workbook
  saveWorkbook(wb, file = file.path(output_dir, "Population_Frequencies.xlsx"), overwrite = TRUE)
}

prep_pseudobulk <- function(integrated_seurat, comparison.col, assay = "RNA"){
  
  set.seed(1234)
  
  # Input checks
  if (!inherits(integrated_seurat, "Seurat")) {
    stop("'integrated_seurat' must be a Seurat object.")
  }
  if (!comparison.col %in% colnames(integrated_seurat@meta.data)) {
    stop(paste0("'", comparison.col, "' not found in metadata columns."))
  }
  if (!"Sample" %in% colnames(integrated_seurat@meta.data)) {
    stop("'Sample' column is required in metadata.")
  }
  if (!assay %in% names(integrated_seurat@assays)) {
    stop(paste0("Assay '", assay, "' not found in Seurat object."))
  }
  
  # Initialize storage
  counts_matrix <- integrated_seurat@assays[[assay]]$counts
  meta_data_full <- data.frame()
  read_data_full <- data.frame(SYMBOL = rownames(counts_matrix))
  
  # Extract unique groups from comparison.col that will be compared in DE analysis
  groups <- unique(integrated_seurat@meta.data[[comparison.col]])
  groups <- groups[!is.na(groups)]  # remove NA groups if present
  
  # Loop through each group
  for (g in groups) {
    subset.seurat <- integrated_seurat[, integrated_seurat@meta.data[[comparison.col]] == g]
    
    # Generate metadata
    meta_data <- subset.seurat@meta.data %>%
      tibble::rownames_to_column("barcodes") %>%
      dplyr::filter(!is.na(Sample)) %>%
      dplyr::mutate(Sample_ID = paste0(Sample, "_", g),
                    Comparisons = .data[[comparison.col]]) %>%
      dplyr::add_count(Sample_ID) %>%
      dplyr::filter(n >= 100)
    
    ### Generate read data
    # read data will have "the reads of all cells belonging to a single sample" 
    # merged together in each column. 
    
    # First, create a list of samples
    samples <-  meta_data[["Sample_ID"]] %>%
      unique()
    
    # Second, initialize read counts matrix
    read_data <- matrix(0,
                        nrow = nrow(counts_matrix),
                        ncol = length(samples),
                        dimnames = list(rownames(counts_matrix), samples))
    
    # Third, add row-wise, the counts of each gene for each sample
    for(i in samples){
      
      # Create a list of cells for each sample
      cells_subset <- meta_data %>% 
        dplyr::filter(Sample_ID == i) %>% 
        dplyr::pull(barcodes)
      
      # Subset counts
      subset_counts <- counts_matrix[, cells_subset, drop=FALSE]
      
      # Sum counts
      read_data[, i] <- Matrix::rowSums(subset_counts)
    }
    
    # Fourth, format metadata and readdata
    cols_to_keep <- intersect(c("Sample_ID", "Comparisons", "Patient", "Condition", "Treatment", "Disease", 
                                "Tissue", "Strain", "Cell_line", "Sex", "n", "barcodes"),
                              colnames(meta_data))
    meta_data <- meta_data %>%
      dplyr::distinct(Sample_ID, .keep_all = TRUE) %>%
      dplyr::select(all_of(cols_to_keep))
    
    read_data <- as.data.frame(read_data) %>% 
      tibble::rownames_to_column("SYMBOL")
    
    # Finally, merge results
    meta_data_full <- dplyr::bind_rows(meta_data_full, meta_data)
    read_data_full <- dplyr::left_join(read_data_full, read_data, by=c("SYMBOL"="SYMBOL"))
  }
  
  return(list(meta_data = meta_data_full, read_data = read_data_full))
}

find_degs_seurat <- function(integrated_seurat, comparison.col, output_dir, celltype.col = "Cell.Type", assay = "RNA"){
  
  set.seed(1234)
  
  # Setup log file
  log_file <- file.path(output_dir, "find_degs_seurat.log")
  log_con <- file(log_file, open = "wt")
  log_msg <- function(...) {
    msg <- paste0(...)
    message(msg)                 # print to console
    writeLines(msg, log_con)     # also save to log file
  }
  on.exit(close(log_con), add = TRUE)  # close file safely when function ends
  
  # Initialize storage
  deg_MAST_list <- list()
  deg_Wilcox_list <- list()
  deg_DESeq2_list <- list()
  counts_list <- list()
  vst_counts_list <- list()
  skipped <- c()
  
  # Define cell types and comparisons
  celltype_levels <- unique(integrated_seurat@meta.data[[celltype.col]])
  comparison_levels <- unique(integrated_seurat@meta.data[[comparison.col]])
  comparisons <- utils::combn(x = comparison_levels, m = 2, simplify = FALSE)
  
  for(celltype in celltype_levels){
    for(comparison in comparisons){
      
      target <- comparison[1]
      reference <- comparison[2]
      
      # Subset Seurat object
      subset_obj <- integrated_seurat[, integrated_seurat@meta.data[[celltype.col]] == celltype]
      subset_obj <- subset_obj[, subset_obj@meta.data[[comparison.col]] %in% c(target, reference)]
      
      # Count cells
      n_cells_ref <- sum(subset_obj@meta.data[[comparison.col]] == reference)
      n_cells_target <- sum(subset_obj@meta.data[[comparison.col]] == target)
      
      # Count samples with >= 100 cells
      n_samples_ref <- subset_obj@meta.data %>%
        dplyr::filter(.data[[comparison.col]] == reference) %>%
        dplyr::count(Sample) %>% dplyr::filter(n >= 100) %>% nrow()
      n_samples_target <- subset_obj@meta.data %>%
        dplyr::filter(.data[[comparison.col]] == target) %>%
        dplyr::count(Sample) %>% dplyr::filter(n >= 100) %>% nrow()
      
      # Skip if too few cells
      if(n_cells_ref < 100 | n_cells_target < 100){
        log_msg("Skipping subset ", celltype, " | ",
                target, " (", n_cells_target, " cells) vs ",
                reference, " (", n_cells_ref, " cells) due to <100 cells")
        skipped <- c(skipped, paste(celltype, target, reference, sep="_"))
        next
      }
      
      # ---- MAST (always run) ----
      deg_obj <- subset_obj
      log_msg("Running MAST for ", celltype, " | ",
              target, " (", n_cells_target, " cells) vs ",
              reference, " (", n_cells_ref, " cells)")
      
      mast_degs <- FindMarkers(
        object = deg_obj,
        ident.1 = target,
        ident.2 = reference,
        group.by = comparison.col,
        assay = assay,
        test.use = "MAST",
        slot = "data",         # MAST must be run on log transformed data
        min.cells.group = 2,
        min.pct = 0.1
      )
      mast_degs <- mast_degs %>%
        tibble::rownames_to_column("SYMBOL") %>%
        dplyr::mutate(Comparison = paste0(celltype, ".", target, ".vs.", reference))
      deg_MAST_list[[paste0(celltype, ".", target, ".vs.", reference)]] <- mast_degs
      
      # ---- Wilcoxon (always run) ----
      deg_obj <- subset_obj
      log_msg("Running Wilcoxon for ", celltype, " | ",
              target, " (", n_cells_target, " cells) vs ",
              reference, " (", n_cells_ref, " cells)")
      wilcox_degs <- FindMarkers(
        object = deg_obj,
        ident.1 = target,
        ident.2 = reference,
        group.by = comparison.col,
        assay = assay,
        test.use = "wilcox",
        slot = "data",
        min.cells.group = 2,
        min.pct = 0.1
      )
      wilcox_degs <- wilcox_degs %>%
        tibble::rownames_to_column("SYMBOL") %>%
        dplyr::mutate(Comparison = paste0(celltype, ".", target, ".vs.", reference))
      deg_Wilcox_list[[paste0(celltype, ".", target, ".vs.", reference)]] <- wilcox_degs
      
      
      # ---- Counts (always run) ----
      deg_obj <- Seurat::AggregateExpression(
        object = subset_obj,
        group.by = c("Sample", comparison.col),
        assays = assay,
        slot = "counts",
        return.seurat = TRUE
      )
      
      raw_counts <- as.data.frame(as.matrix(SeuratObject::GetAssayData(deg_obj, assay = assay, slot = "counts"))) 
      counts_list[[paste0(celltype, ".", target, ".vs.", reference)]] <- raw_counts %>%
        tibble::rownames_to_column(var = "SYMBOL")
      
      # VST normalized counts
      meta_data <- deg_obj@meta.data[colnames(raw_counts), , drop = FALSE]
      rownames(meta_data) <- colnames(raw_counts)
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts,
                                            colData = meta_data,
                                            design = ~1)
      vst_counts <- tryCatch(
        {
          vst <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
          as.data.frame(SummarizedExperiment::assay(vst)) %>%
            rownames_to_column(var = "SYMBOL")
        },
        error = function(e) {
          message("VST failed due to zeros in counts. Returning empty data frame instead.")
          data.frame(SYMBOL = character(0))  # ensures at least SYMBOL column exists
        }
      )
      vst_counts_list[[paste0(celltype, ".", target, ".vs.", reference)]] <- vst_counts
      
      # ---- DESeq2 pseudobulk ----
      if(n_samples_ref >= 2 & n_samples_target >= 2){
        log_msg("Running DESeq2 for ", celltype, " | ",
                target, " (", n_cells_target, " cells) vs ",
                reference, " (", n_cells_ref, " cells)")
        
        deseq2_degs <- FindMarkers(
          object = deg_obj,
          ident.1 = target,
          ident.2 = reference,
          group.by = comparison.col,
          assay = assay,
          test.use = "DESeq2",
          slot = "counts",      # DESeq2 must be run on raw counts
          min.cells.group = 2   # DESeq2 uses this argument to determine number of groups
        )
        deseq2_degs <- deseq2_degs %>%
          tibble::rownames_to_column("SYMBOL") %>%
          dplyr::mutate(Comparison = paste0(celltype, ".", target, ".vs.", reference))
        deg_DESeq2_list[[paste0(celltype, ".", target, ".vs.", reference)]] <- deseq2_degs
      }
    }
  }
  
  # ---- Save to separate Excel files ----
  save_list_to_xlsx <- function(lst, file_path) {
    wb <- openxlsx::createWorkbook()
    
    # Create index sheet first
    index_df <- data.frame(
      Sheet = paste0("Sheet", seq_along(lst)),
      Comparison = names(lst),
      stringsAsFactors = FALSE
    )
    openxlsx::addWorksheet(wb, "Index")
    openxlsx::writeData(wb, "Index", index_df)
    
    # Add each comparison as Sheet1, Sheet2, ...
    for (i in seq_along(lst)) {
      sheet_nm <- paste0("Sheet", i)
      openxlsx::addWorksheet(wb, sheetName = sheet_nm)
      openxlsx::writeData(wb, sheet = sheet_nm, lst[[i]])
    }
    openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
  }
  
  all_DESeq2 <- dplyr::bind_rows(deg_DESeq2_list)
  all_MAST <- dplyr::bind_rows(deg_MAST_list)
  all_Wilcox <- dplyr::bind_rows(deg_Wilcox_list)
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "DEGs_DESeq2")
  openxlsx::writeData(wb, "DEGs_DESeq2", all_DESeq2)
  openxlsx::addWorksheet(wb, "DEGs_MAST")
  openxlsx::writeData(wb, "DEGs_MAST", all_MAST)
  openxlsx::addWorksheet(wb, "DEGs_Wilcox")
  openxlsx::writeData(wb, "DEGs_Wilcox", all_Wilcox)
  openxlsx::saveWorkbook(wb, file.path(output_dir, "Seurat_DEGs.xlsx"), overwrite = TRUE)
  
  
  if (length(counts_list) > 0) save_list_to_xlsx(counts_list, file.path(output_dir, "Seurat_Raw_Counts.xlsx"))
  if (length(vst_counts_list) > 0) save_list_to_xlsx(vst_counts_list, file.path(output_dir, "Seurat_VST_Counts.xlsx"))
  
  # ---- Report skipped ----
  if (length(skipped) > 0) {
    log_msg("Skipped comparisons (low cell counts):")
    log_msg(paste(" -", skipped, collapse = "\n"))
  }
  
  # ---- Return ----
  return(list(
    DESeq2 = deg_DESeq2_list,
    MAST   = deg_MAST_list,
    Wilcox = deg_Wilcox_list
  ))
}

# We need to run SCT on each assay. But for cell cycle scoring, NA is present
# in barcodes from other assays like 002um and 016um. So, it throws error.
# Therefore, create a separate seurat object for each bin size so that each
# seurat object has only one assay equivalent to "RNA" assay of single cell data

h5ad_to_seurat_batch <- function(path_to_h5ad) {
  
  # Find all .h5ad files (path must be a directory)
  h5ad_files <- list.files(path = path_to_h5ad,
                           pattern = "\\.h5ad$",
                           full.names = TRUE)
  
  if(length(h5ad_files) == 0) {
    stop("No .h5ad files found in the specified path.")
  }
  
  seurat_list <- list()
  
  for(h5ad_file in h5ad_files) {
    
    message("Processing file: ", h5ad_file)
    
    # Prepare .h5seurat filename
    h5seurat_file <- sub("\\.h5ad$", ".h5seurat", h5ad_file)
    
    # Convert .h5ad to .h5seurat
    if(!file.exists(h5seurat_file)) {
      message("Converting .h5ad to .h5seurat...")
      # IMPORTANT CORRECTION: 'dest' must be the full path to the output file
      # or simply the filename if the working directory is correct, but passing 
      # the full path to the source is safer.
      SeuratDisk::Convert(source = h5ad_file,
                          dest = "h5seurat", 
                          assay = "RNA",
                          overwrite = FALSE)
    } else {
      message(".h5seurat file already exists: ", h5seurat_file)
    }
    
    # Load Seurat object
    message("Loading Seurat object...")
    # If you get error, "Error: Missing required datasets 'levels' and 'values'",
    # set meta.data = FALSE, misc = FALSE within SeuratDisk::LoadH5Seurat()
    seurat_obj <- SeuratDisk::LoadH5Seurat(file = h5seurat_file, 
                                           meta.data = FALSE, 
                                           misc = FALSE)
    
    # Extract and Assign Metadata (adata.obs to seurat_obj@meta.data)
    message("Extracting and assigning metadata from AnnData (.obs) to Seurat (@meta.data)...")
    adata <- anndata::read_h5ad(h5ad_file)
    obs_df <- as.data.frame(adata$obs)
    # adata.obs.to_excel("/hpc/home/kailasamms/scratch/scRNASeq_PanCancer/Metadata.xlsx", index=True)
    # seurat_obj@meta.data <- read.xlsx("/hpc/home/kailasamms/scratch/scRNASeq_PanCancer/Metadata.xlsx") %>%
    #   tibble::column_to_rownames("X1")
    
    # Ensure rownames match Seurat object cells
    if (!all(rownames(obs_df) %in% colnames(seurat_obj))) {
      warning("Cell names in metadata do not match Seurat object. Check rownames!")
    }
    seurat_obj@meta.data <- obs_df
    
    # Save Seurat object as RDS
    rds_file <- base::gsub(pattern = "\\.h5ad$", replacement = ".rds", x = h5ad_file)
    message("Saving Seurat object to RDS: ", rds_file)
    saveRDS(object = seurat_obj, file = rds_file)
    
    seurat_list[[basename(rds_file)]] <- seurat_obj
  }
  
  message("All .h5ad files processed. Returning list of Seurat objects.")
  return(seurat_list)
}

# ---- SURVIVAL RELATED FUNCTIONS ----

# NOTE: When plotting KM curves for individual genes, non-transformed data,
# log-transformed or median centered log-transformed data give identical results
# for all cutoff methods except thirds. RECOMMEND using vst counts from DESeq2 or
# log transformed counts

# When plotting KM curves for signature scores, RECOMMEND using vst counts from 
# DESeq2 or log transformed counts. advanced_z() automatically does the 
# necessary median centering before signature score calculation.

# Display pre-defined survival analysis scenarios
show_survival_scenarios <- function() {
  scenarios <- tibble::tribble(
    ~Scenario, ~Description, ~stratify_var, ~substratify_var, ~facet_var, ~multiple_cutoff, ~Curves_per_plot, ~Plots_facets, ~Curve_labels,
    "(i)", "Survival based on gene A", "Gene A", "–", "–", "–", 2, 1, "HIGH vs LOW",
    "(ii)", "Survival based on gene A + Sex", "Gene A", "Sex", "–", "FALSE / TRUE", 4, 1, "HIGH/LOW × Male/Female",
    "(iii)", "Survival based on gene A, faceted by Sex", "Gene A", "–", "Sex", "FALSE / TRUE", 2, 2, "HIGH vs LOW",
    "(iv)", "Survival based on gene A + Smoking, faceted by Sex", "Gene A", "Smoking", "Sex", "FALSE / TRUE", 4, 2, "HIGH/LOW × Yes/No",
    "(v)", "Survival based on Sex", "Sex", "–", "–", "–", 2, 1, "Male vs Female",
    "(vi)", "Survival based on Sex + Race", "Sex", "Race", "–", "–", 4, 1, "Male/Female × White/Black",
    "(vii)", "Survival based on Sex, faceted by Race", "Sex", "–", "Race", "–", 2, 2, "Male vs Female",
    "(viii)", "Survival based on Sex + Smoking, faceted by Race", "Sex", "Smoking", "Race", "FALSE / TRUE", 4, 2, "Male/Female × Yes/No"
  )
  return(scenarios)
}

# Calculate multi-gene signature scores
# Described in Levine et al https://doi.org/10.1186/gb-2006-7-10-r93
advanced_z <- function(gene_set, expr_matrix) {
  
  ix <- toupper(rownames(expr_matrix)) %in% toupper(gene_set)
  cat("Genes found:", sum(ix), "\n")
  
  # --------------------------------------------------------------------------------------------------------------------------------------------|
  # Method                 | What it does                                 | Goal                                  | Use Case
  # ---------------------- | ------------------------------------------------------------------------------------ | ----------------------------|
  # Row-Wise Centering     | Subtracts the median expression of a gene    | To look at how a gene's expression in | Differential Expression,    |
  # (MARGIN=1)             | across all samples from that gene's          | a specific sample deviates from its   | Clustering genes, Heatmaps  |
  #                        | expression in each sample.                   | typical expression across all samples.| (to compare gene behavior). |
  #
  # Column-Wise Centering  | Subtracts the median expression of all genes | To look at how a gene's expression    | Signature Scoring,          |
  # (MARGIN=2)             | in a sample from each gene's expression in   | average deviates from the overall     | Single-Sample Z-Score       |
  #                        | that sample.                                 | transcriptome of that specific sample.| (ssGSEA-like).              |
  # ---------------------------------------------------------------------------------------------------------------------------------------------
  
  # Median-center each sample across genes
  expr_matrix_centered <- base::sweep(x = expr_matrix, 
                                      MARGIN = 2, 
                                      STATS = apply(expr_matrix, 2, median, na.rm = TRUE), 
                                      FUN = "-")
  
  if (sum(ix) > 1) {
    avg_gene_set <- base::apply(X = expr_matrix[ix, , drop = FALSE], MARGIN = 2, FUN = mean, na.rm = TRUE)
    avg_all      <- base::apply(X = expr_matrix,                     MARGIN = 2, FUN = mean, na.rm = TRUE)
    sd_all       <- base::apply(X = expr_matrix,                     MARGIN = 2, FUN = sd,   na.rm = TRUE)
    
    z <- (avg_gene_set - avg_all) * sqrt(sum(ix)) / sd_all
  } else {
    stop("Cannot use sum(ix) genes for signature score calculation")
  }
  return(z)
}

# Calculate cutoffs (expression based survival) for each facet
calc_cutoffs <- function(cutoff_df, stratify_var, survival_params){
  
  expr_values <- cutoff_df[[stratify_var]]
  
  # Compute quantiles
  qs <- stats::quantile(expr_values,
                        probs = c(0, 0.25, 0.33, 0.5, 0.66, 0.75),
                        na.rm = TRUE)
  iqr <- stats::IQR(expr_values, na.rm = TRUE)
  
  # Determine cutoffs
  cutoffs <- switch(
    survival_params$cutoff_method,
    "median"   = list(lower = qs["50%"], upper = qs["50%"], middle = NA),
    "tertile"  = list(lower = qs["33%"], upper = qs["66%"], middle = NA),
    "quartile" = list(lower = qs["25%"], upper = qs["75%"], middle = qs["50%"]),
    "thirds"   = {
      lower <- max(qs["0%"] - 1.5*iqr, 0)
      upper <- qs["75%"] + 1.5*iqr
      list(lower = lower + (upper-lower)/3, upper = lower + (upper-lower)*2/3, middle = NA)
    },
    "optimal" = tryCatch({
      res <- survminer::surv_cutpoint(data = cutoff_df,
                                      time = survival_params$time_col,
                                      event = survival_params$status_col,
                                      variables = stratify_var)
      list(lower = res$cutpoint$cutpoint, upper = res$cutpoint$cutpoint, middle = NA)
    }, error = function(e) {
      list(lower = NA, upper = NA, middle = NA)
    }))
  
  # Categorize expression into bins
  model_col <- paste0("model_", stratify_var)
  cutoff_df <- cutoff_df %>%
    dplyr::filter(!is.na(.data[[stratify_var]])) %>%
    dplyr::mutate(!!model_col := dplyr::case_when(.data[[stratify_var]] > cutoffs$upper ~ "HIGH",
                                           .data[[stratify_var]] <= cutoffs$lower ~ "LOW",
                                           (!is.na(cutoffs$middle) & .data[[stratify_var]] > cutoffs$middle) ~ "MED_HIGH",
                                           (!is.na(cutoffs$middle) & .data[[stratify_var]] <= cutoffs$middle) ~ "MED_LOW",
                                           TRUE ~ "MID"),
                  !!model_col := factor(.data[[model_col]], 
                                        levels = c("LOW", "MED_LOW", "MID", "MED_HIGH", "HIGH"))) %>%
    dplyr::select(Sample_ID, all_of(model_col))
  
  # Optionally plot only HIGH and LOW
  if (!survival_params$show_all_bins) {
    cutoff_df <- cutoff_df %>%
      dplyr::filter(.data[[model_col]] %in% c("HIGH", "LOW"))
  }
  return(cutoff_df)
}

# Calculate cox model stats (HR, CI, p vals) for each facet
calc_cox_stats <- function(facet_df, stratify_var, surv_formula, survival_params){
  
  # ---- Cox model ----
  
  # Ensure model column is a factor
  model_col <- paste0("model_", stratify_var)
  facet_df[[model_col]] <- factor(facet_df[[model_col]])
  
  # Fit Cox model
  cox_model <- survival::coxph(formula = surv_formula, data = facet_df)
  
  # Cox model coefficients
  cox_coef_df <- summary(cox_model)$coefficients
  cox_ci_df <- as.data.frame(confint(cox_model))
  baseline <- levels(facet_df[[model_col]])[1]  # the factor baseline
  cox_df <- data.frame(Gene = stratify_var,
                       Target = gsub(paste0("^model_", stratify_var), "", rownames(cox_coef_df)), # remove "model" prefix if present
                       Reference = baseline,
                       HR = exp(cox_coef_df[, "coef"]),
                       CI_lower = exp(cox_ci_df[, 1]),
                       CI_upper = exp(cox_ci_df[, 2]),
                       pval = cox_coef_df[, "Pr(>|z|)"],
                       stringsAsFactors = FALSE) %>%
    dplyr::mutate(contrast = paste0(Target, " / ", Reference)) %>%
    dplyr::select(Gene, contrast, HR, pval, CI_lower, CI_upper, Target, Reference) %>%
    tibble::remove_rownames() 
  
  # ---- Optional: emmeans for pairwise contrasts ----
  # Step 1: estimated marginal means on log-hazard scale
  emm <- emmeans::emmeans(object = cox_model, specs = as.formula(paste0("~", model_col)))
  
  # Step 2: pairwise contrasts on HR scale
  pairwise <- emmeans::contrast(object = emm, method = "pairwise", type = "response")
  
  # Step 3: add confidence intervals
  pairwise_ci <- stats::confint(pairwise)
  
  # Step 4: combine into clean data.frame
  pairwise_df <- dplyr::left_join(x = as.data.frame(pairwise) %>%
                                    dplyr::select(contrast, ratio, p.value),
                                  y = as.data.frame(pairwise_ci) %>%
                                    dplyr::select(contrast, asymp.LCL, asymp.UCL),
                                  by = c("contrast"="contrast")) %>%
    dplyr::rename(HR = ratio,
                  CI_lower = asymp.LCL,
                  CI_upper = asymp.UCL,
                  pval = p.value) %>%
    dplyr::mutate(Target = sub(" / .*", "", contrast),
                  Reference = sub(".* / ", "", contrast),
                  Gene = stratify_var)
  
  # Step 5: Calculate HR, CI for reverse contrasts
  reversed <- pairwise_df %>%
    dplyr::mutate(Reference = sub(" / .*", "", contrast),
                  Target = sub(".* / ", "", contrast),
                  contrast = paste(Target, "/", Reference),
                  # reciprocal HR
                  HR = 1 / HR,
                  # swap CI bounds
                  CI_lower = 1 / CI_upper,
                  CI_upper = 1 / CI_lower)
  
  # Step 6: Merge all stats
  emmeans_df <- dplyr::bind_rows(pairwise_df, reversed)
  
  # Step 7: Compute non-parametric p-values
  emmeans_df <- sapply(X = emmeans_df$contrast, FUN = calc_pvals, 
                       stratify_var =  stratify_var,
                       facet_df = facet_df, survival_params = survival_params) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("contrast") %>%
    dplyr::right_join(emmeans_df, by=c("contrast" = "contrast")) %>%
    dplyr::select(Gene, contrast, Target, Reference, HR, pval, CI_lower, CI_upper, everything())
  
  cox_df <- sapply(X = cox_df$contrast, FUN = calc_pvals,
                   stratify_var =  stratify_var,
                   facet_df = facet_df, survival_params = survival_params) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("contrast") %>%
    dplyr::right_join(cox_df, by=c("contrast" = "contrast")) %>%
    dplyr::select(Gene, contrast, Target, Reference, HR, pval, CI_lower, CI_upper, everything())
  
  # ---- Return both ----
  list(cox_model_df = cox_df,
       emmeans_df = emmeans_df)
}

# Calculate all 7 non-parametric p-values for each contrast in cox_df or emmeans_df
calc_pvals <- function(contrast, facet_df, stratify_var, survival_params) {
  
  model_col <- paste0("model_", stratify_var)
  g1 <- sub(".* / ", "", contrast)
  g2 <- sub(" / .*", "", contrast)
  df_pair <- subset(facet_df, facet_df[[model_col]] %in% c(g1, g2))
  surv_obj <- survival::Surv(
    time   = df_pair[[survival_params$time_col]],
    event  = df_pair[[survival_params$status_col]],
    type   = "right",
    origin = 0
  )
  surv_form <- as.formula(paste("surv_obj ~", model_col))
  fit <- survminer::surv_fit(formula = surv_form , data = df_pair)
  
  test_methods <- c("survdiff", "1", "n", "sqrtN", "S1", "S2", "FH_p=1_q=1")
  sapply(X = test_methods, FUN = function(test.method) {
    res <- tryCatch(
      survminer::surv_pvalue(fit = fit,
                             method = test.method,
                             test.for.trend = FALSE,
                             combine = FALSE)[[2]],
      error = function(e) NA
    )})
}

# Plot survival for each facet
plot_facets <- function(facet_df, stratify_var, surv_curve, cox_df, surv_type, facet_group, survival_params){
  
  if (survival_params$plot_curve) {
    
    model_col <- paste0("model_", stratify_var)
    
    # Legend labels
    # IMPORTANT: ggsurvplot() labels the groups in alphabetical order. So,
    # when we want to use custom labels, initialize them in alphabetical order.
    # Eg: c("High", "Low") instead of c("Low, "High")
    legend_label <- sort(unique(facet_df[[model_col]]))
    
    # Legend title
    if (surv_type %in% c("single_gene", "multi_gene", "signature")){
      legend_title <- paste0(c("Expression", survival_params$substratify_var, survival_params$facet_var), collapse = "_")
    } else if (surv_type == "meta"){
      legend_title <- paste0(c(survival_params$stratify_var, survival_params$substratify_var, survival_params$facet_var), collapse = "_")
    }
    
    # X-axis breaks
    # We want a maximum of 10 timepoint intervals that are multiples of 12
    max_time <- max(facet_df[[survival_params$time_col]], na.rm = TRUE)
    n <- max(floor(max_time / 10 / 12) * 12, 1)
    breaks <- if (max_time %/% n <= 10) n else n + 12
    
    # Plot KM curve
    # NOTE: Use survminer::ggsurvplot() instead of base::plot()
    # ggsurvplot() produces a list of ggplot objects: survival curve and risk table
    surv_plot <- survminer::ggsurvplot(
      fit = surv_curve,
      pval = FALSE,
      palette = survival_params$color_palette,
      linetype = "solid",
      size = 1.5,                  # thickness of line
      
      # Format the legend
      legend = "top",
      legend.title = legend_title,
      legend.labs = legend_label,
      
      # Format the axes
      break.time.by = breaks,       # break X axis in time intervals of 12 months
      xlab = "Time",
      ylab = "Survival Probability",
      title =  paste(na.omit(c(stratify_var, survival_params$substratify_var, 
                               facet_group, survival_params$cutoff_method)), 
                     collapse = "."),
      
      # Format confidence intervals
      conf.int = survival_params$conf_interval,
      #conf.int.fill = ?,               # color to fill confidence interval
      conf.int.style = "ribbon",        # confidence interval style
      conf.int.alpha = 0.3,             # confidence fill color transparency
      
      # Format the risk table
      risk.table = survival_params$plot_risk_table,
      risk.table.title = "Number at risk",
      risk.table.y.text.col = TRUE,     # color of risk table text annotations
      risk.table.pos = "out",           # draw risk table outside survival plot
      
      # Format the censor points
      censor = TRUE,
      censor.shape = '|',
      censor.size = 4
    )
    
    # Adjust x-axis limits for survival plot and risk table
    surv_plot$table <- surv_plot$table +
      coord_cartesian(x = c(0, ceiling(max_time / breaks) * breaks), clip = "off")
    
    surv_plot$plot <- surv_plot$plot +
      coord_cartesian(x = c(0, ceiling(max_time / breaks) * breaks), clip = "off")
    
    # Create a text annotation showing p-value, HR with CI, and method
    # [Accurate ONLY if 2 curves are present. If more, refer excel file and manually add to image]
    if (dplyr::n_distinct(facet_df[[model_col]]) == 2){
      method_plot <- "log-rank"
      p_plot <- formatC(cox_df$pval[1], format = "e", digits = 1)
      hr_plot <- round(cox_df$HR[1], 1)
      ci_lower_plot <- round(cox_df$CI_lower[1], 1)
      ci_upper_plot <- round(cox_df$CI_upper[1], 1)
      
      survplot_stats_grob <- grobTree(textGrob(
        label = paste0("p = ", p_plot,
                       "\nHR = ", hr_plot, " [", ci_lower_plot, ", ", ci_upper_plot, "]",
                       "\nMethod = ", method_plot),
        x = 0.50, y = 0.90, hjust = 0,
        gp = grid::gpar(fontfamily = "Times", fontface = "bold", col = "black", fontsize = 10)
      ))
      
      # Merge text annotation with survival plot
      surv_plot$plot <- surv_plot$plot %++%
        ggplot2::annotation_custom(survplot_stats_grob)
    }
    
    # NOTE: Using cowplot() and then ggsave() works nicely as compared to 
    # saving directly using ggsave()
    p <- cowplot::plot_grid(plotlist = surv_plot,
                       align = "hv",
                       axis = "tblr",
                       nrow = 2,
                       ncol = 1,
                       rel_widths = 1,
                       rel_heights = c(1, 0.45),
                       labels = NULL,
                       label_size = 14,
                       label_fontface = "bold")

    return(p)
  }
}

survival_analysis <- function(meta_data, expr_data = NULL, survival_params) {
  
  # ---- Input checks & parameter extraction ----
  
  # Create output directory if missing
  if (!dir.exists(survival_params$output_dir)) {
    dir.create(survival_params$output_dir, recursive = TRUE)
  }
  
  rownames(expr_data) <- make.names(rownames(expr_data))
  colnames(expr_data) <- make.names(colnames(expr_data))
  meta_data$Sample_ID <- make.names(meta_data$Sample_ID)
  survival_params$stratify_var <- make.names(survival_params$stratify_var)
  
  # Extract parameters
  stratify_vars   <- survival_params$stratify_var
  substratify_var <- survival_params$substratify_var
  facet_var       <- survival_params$facet_var
  time_col        <- survival_params$time_col
  status_col      <- survival_params$status_col
  cutoff_method   <- survival_params$cutoff_method
  multiple_cutoff <- survival_params$multiple_cutoff
  show_all_bins   <- survival_params$show_all_bins
  sig_score       <- survival_params$sig_score
  
  # Must provide stratify_var
  if (is.null(stratify_vars) || length(stratify_vars) == 0) {
    stop("Must provide a non-empty stratify_var")
  }
  
  # substratify_var, if provided, must exist in meta_data
  if (!is.null(substratify_var) && !substratify_var %in% colnames(meta_data)) {
    stop("substratify_var not found in meta_data columns")
  }
  
  # facet_var, if provided, must exist in meta_data
  if (!is.null(facet_var) && !facet_var %in% colnames(meta_data)) {
    stop("facet_var not found in meta_data columns")
  }
  
  # Must define substratify_var if multiple_cutoff = TRUE
  if (isTRUE(multiple_cutoff) && is.null(substratify_var)) {
    stop("multiple_cutoff = TRUE requires substratify_var")
  }
  
  # Check stratify_var against both metadata and expression data and determine
  # type of survival analysis
  missing_genes <- setdiff(stratify_vars, rownames(expr_data))
  valid_genes   <- intersect(stratify_vars, rownames(expr_data))
  in_meta       <- all(stratify_vars %in% colnames(meta_data))
  
  if (!in_meta) {
    # Expression-based stratification
    if (length(missing_genes) > 0) {
      warning("Some requested genes not found in expr_data: ", paste(missing_genes, collapse = ", "))
    }
    
    if (length(valid_genes) == 0) {
      stop("stratify_var must match either a column in meta_data OR genes in expr_data.")
      
    } else if (length(valid_genes) > 1 && isTRUE(sig_score)) {
      message("Proceeding with signature expression-based survival analysis (", length(valid_genes), " valid genes).")
      surv_type <- "signature"
      
    } else if (length(valid_genes) > 1 && !isTRUE(sig_score)) {
      message("Proceeding with expression-based survival analysis for ", length(valid_genes), " gene(s).")
      surv_type <- "multi_gene"
      
    } else if (length(valid_genes) == 1) {
      message("Proceeding with expression-based survival analysis for single gene.")
      surv_type <- "single_gene"
    }
    
  } else {
    # Metadata-based stratification
    if (length(valid_genes) == 0) {
      message("Proceeding with metadata-based survival analysis.")
      surv_type <- "meta"
    } else {
      stop("stratify_var matches BOTH a column in meta_data and gene(s) in expr_data — ambiguous.")
    }
  }
  
  # Define facet groups
  if (!is.null(facet_var) && facet_var %in% colnames(surv_df)) {
    facet_groups  <- unique(surv_df[[facet_var]])
  } else {
    facet_groups  <- NA_character_ # placeholder for whole dataset
  }
  
  # ---- Format metadata ----
  
  meta_data <- meta_data %>% 
    dplyr::mutate(Sample_ID = make.names(names = Sample_ID), 
                  !!time_col := as.numeric(.data[[time_col]])) %>%
    dplyr::filter(.data[[time_col]] > 0 & !is.na(.data[[time_col]])) %>%
    dplyr::distinct(Sample_ID, .keep_all = TRUE)
  
  # ---- Prepare expression data ----
  if (surv_type == "signature") {
    # Signature score survival (whole dataset, even if substratify_var defined)
    sig_scores <- advanced_z(gene_set = stratify_vars, expr_matrix = expr_data)
    expr_df    <- as.data.frame(sig_scores, check.names = FALSE) %>%
      dplyr::rename(sig_score = identity(1)) %>%
      tibble::rownames_to_column(var = "Sample_ID")
    
    # Update local variable and survival_params
    stratify_vars <- "sig_score"
    survival_params$stratify_var <- stratify_vars
    
  } else if (surv_type %in% c("single_gene", "multi_gene")) {
    # Single/Multi-gene survival (expression-based)
    expr_df <- expr_data[stratify_vars, , drop = FALSE] %>%
      t() %>%
      as.data.frame(check.names = FALSE) %>%
      tibble::rownames_to_column(var = "Sample_ID")
    
  } else if (surv_type == "meta") {
    # Metadata-based survival
    expr_df <- meta_data %>%
      dplyr::select(Sample_ID, dplyr::all_of(stratify_vars))
    
  } else {
    stop("Invalid stratify_var: must be gene(s) in expr_data or a column in meta_data.")
  }
  
  colnames(expr_df) <- make.names(colnames(expr_df))
  
  # ---- Merge expression data with metadata ----
  
  keep_cols <- unique(c(time_col, status_col, stratify_vars, substratify_var, facet_var))
  
  surv_df <- expr_df %>%
    dplyr::inner_join(meta_data, by = c("Sample_ID"="Sample_ID")) %>%
    dplyr::select(Sample_ID, dplyr::all_of(keep_cols))
  
  if (nrow(surv_df) == 0) stop("No overlapping Sample_IDs between expr_data and meta_data.")
  
  
  # ---- Define model column for metadata based survival ----
  
  if (surv_type == "meta") {
    model_col <- paste0("model_", stratify_vars)
    surv_df <- surv_df %>%
      dplyr::mutate(!!model_col := paste(.data[[stratify_vars]],
                                  if (!is.null(substratify_var)) .data[[substratify_var]] else NULL,
                                  if (!is.null(facet_var)) .data[[facet_var]] else NULL,
                                  sep = "_"))
  }

  # ---- Define model column for expression based survival ----
  
  if (surv_type %in% c("signature", "single_gene", "multi_gene")) {
    
    full_merged_df <- data.frame(Sample_ID = surv_df$Sample_ID)
    # Loop over each gene
    for (stratify_var in stratify_vars){
      
      model_col <- paste0("model_", stratify_var)
      merged_df <- tibble::tibble()
      # Loop over each facet group
      for (facet_group in facet_groups){
        
        # Subset by facet (or use full data if NA)
        if (is.na(facet_group)){
          facet_df <- surv_df
        } else{
          facet_df <- surv_df %>%
            dplyr::filter(.data[[facet_var]] == facet_group)
        }
        
        # Determine cutoff groups
        if (!is.null(substratify_var) && substratify_var %in% colnames(surv_df) && isTRUE(multiple_cutoff)) {
          cutoff_groups <- unique(facet_df[[substratify_var]])
        } else {
          cutoff_groups <- NA_character_
        }
        
        # Calculate cutoffs for each cutoff_group
        for (cutoff_group in cutoff_groups){
          
          # Subset by cutoff_group (or use facet data if NA)
          if (is.na(cutoff_group)) {
            cutoff_df <- facet_df
          } else {
            cutoff_df <- facet_df %>% 
              dplyr::filter(.data[[substratify_var]] == cutoff_group)
          }
          
          # Calculate cutoffs for this subset
          df <- calc_cutoffs(cutoff_df = cutoff_df,
                             stratify_var = stratify_var,
                             survival_params = survival_params)
          
          merged_df <- dplyr::bind_rows(merged_df, df)
        }
      }
      
      # Append substratify_var to model if defined
      if (!is.null(substratify_var)) {
        merged_df <- merged_df %>%
          dplyr::mutate(!!model_col := paste0(model_col, "_", .data[[substratify_var]]))
      }
      
      # Add model columns for every gene
      full_merged_df <- dplyr::left_join(x = full_merged_df, 
                                         y = merged_df,
                                         by = c("Sample_ID" = "Sample_ID"))
    }
    
    # Merge classifications into surv_df
    surv_df <- surv_df %>%
      dplyr::left_join(full_merged_df, by=c("Sample_ID"="Sample_ID"))
  }
  
  
  # ---- Main loop: facets, compute stats, save, plot ----
  
  # Initialize empty data frames to store all stats
  all_cox_df     <- tibble::tibble()
  all_emmeans_df <- tibble::tibble()
  
  # Open PDF for all KM plots
  pdf(file = file.path(survival_params$output_dir, "KM_curves.pdf"),
      width = 7, height = 7)  # open PDF device
  
  # Loop over each gene
  for (stratify_var in stratify_vars){
    
    model_col <- paste0("model_", stratify_var)
    
    # Plot a survival curve for each facet
    for (facet_group in facet_groups){
      
      # Subset by facet (or use full data if NA)
      if (is.na(facet_group)){
        facet_df <- surv_df
      } else{
        facet_df <- surv_df %>%
          dplyr::filter(.data[[facet_var]] == facet_group)
      }
      
      # Check if each facet has least 2 groups in model column
      if (dplyr::n_distinct(facet_df[[model_col]]) < 2) {
        message("Skipping facet (single group): ", facet_group)
        next
      } else if (base::setequal(x = unique(facet_df[[model_col]]), y = c("HIGH", "LOW"))) {
        survival_params$color_palette <- c(HIGH = "#d73027", LOW = "#0c2c84")
      }
      
      # Create a survival object (Alive = 0, Dead = 1)
      surv_object <- survival::Surv(time   = facet_df[[time_col]],
                                    event  = facet_df[[status_col]],
                                    type   = "right",
                                    origin = 0)
      
      # Create a formula for survival analysis
      surv_formula <-  as.formula(paste("surv_object ~", model_col))
      
      # Create a fit for Kaplan-Meier curve
      # NOTE: survival::survfit() gives error in ggsurvplot()
      surv_curve <- survminer::surv_fit(formula   = surv_formula,
                                        data      = facet_df)
      
      # Compute Cox & emmeans stats for each facet
      stats_list <- calc_cox_stats(facet_df, stratify_var, surv_formula, survival_params)
      
      # Add identifiers for traceability
      stats_list$cox_model_df$Facet <- facet_group
      stats_list$emmeans_df$Facet   <- facet_group
      
      # Accumulate stats
      all_cox_df     <- dplyr::bind_rows(all_cox_df, stats_list$cox_model_df)
      all_emmeans_df <- dplyr::bind_rows(all_emmeans_df, stats_list$emmeans_df)
      
      # Plot survival curve for each facet
      cox_df <- stats_list$cox_model_df
      p <- plot_facets(facet_df, stratify_var, surv_curve, cox_df, surv_type, facet_group, survival_params)
      
      # Print the plot to the current PDF page
      print(p)
    }
  }
  
  dev.off()  # close PDF device
  
  # Save combined stats to a single Excel workbook
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "cox_stats")
  openxlsx::writeData(wb, sheet = "cox_stats", x = all_cox_df)
  openxlsx::addWorksheet(wb, sheetName = "emmeans_stats")
  openxlsx::writeData(wb, sheet = "emmeans_stats", x = all_emmeans_df)
  openxlsx::addWorksheet(wb, sheetName = "surv_df")
  openxlsx::writeData(wb, sheet = "surv_df", x = surv_df)
  openxlsx::saveWorkbook(wb, file = file.path(survival_params$output_dir, "Survival_Stats.xlsx"),
                         overwrite = TRUE)
}

survival_params <- list(
  
  # Stratification (Expression + Metadata-based survival)
  stratify_var     = NULL,          # one or more genes or metadata columns
  sig_score        = FALSE,         # TRUE = combine genes into one signature score
  substratify_var  = NULL,          # optional metadata column for sub-stratification
  facet_var        = NULL,          # optional faceting variable
  
  # Cutoff settings (ONLY for Expression-based survival)
  cutoff_method    = "optimal",      # mean, median, quartile, tertile, optimal, thirds
  # median   : splits samples into 2 bins (below 50%, above 50%)
  # tertile  : splits samples into 3 bins (below 33%, 33%-67%, above 67%)
  # quartile : splits samples into 4 bins (below 25%, 25%-50%, 50%-75%, above 75%)
  # optimal  : splits samples into 2 bins (above & below optimum cutoff)
  # thirds   : splits samples into 3 bins (bottom 33%, middle33%, top 33% based on expression range)
  show_all_bins    = FALSE,          # TRUE = plot all bins (LOW, HIGH, MID/MED_HIGH+MED_LOW)
  multiple_cutoff  = FALSE,          # TRUE = compute cutoffs separately for substratify_var
  
  # Plot settings
  sig_score        = FALSE,          # TRUE = combine genes into one signature score
  conf_interval    = FALSE,          # TRUE = show confidence interval in survival curve
  plot_curve       = TRUE,           # TRUE = plot the survival curve
  plot_risk_table  = TRUE,           # TRUE = plot the risk table below the curve
  color_palette    = custom_palette, # vector of colors for groups c("#d73027","#0c2c84")
  
  # Survival data columns
  time_col         = "Time",         # metadata column containing Time values
  status_col       = "Status",       # metadata column containing Status values
  
  # Output
  prefix           = "",
  output_dir      = "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop"
)

# Run to understand how to define parameters for the survival function
show_survival_scenarios()

# Run your survival analysis
#survival_analysis(meta_data, expr_data, survival_params)

#******************************************************************************#
#                       SURVIVAL CURVE RELATED FUNCTIONS                       #
#******************************************************************************#

# Read this paper for survival analysis
# https://doi.org/10.1093/jncimonographs/lgu024

# NOTE:  Output of prep_expr_df is df
#log_norm_counts is matrix with SYMBOLS as rownames
prep_expr_df <- function(log_norm_counts, meta_data, plot_genes, survival_params){
  
  # Merge expression data with survival data
  if (survival_params$gene_sig_score == TRUE){
    
    # Calculate gene signature score
    expr_df <- as.data.frame(advanced_Z(plot_genes, log_norm_counts))
    
    expr_df <- expr_df %>%
      data.frame() %>%
      dplyr::rename(combined.exp = identity(1)) %>%
      tibble::rownames_to_column("Sample_ID") %>%
      dplyr::inner_join(meta_data, by=c("Sample_ID"="Sample_ID")) %>%
      dplyr::select(Sample_ID, combined.exp, Time, Status)
  } else {
    expr_df <- log_norm_counts %>%
      t() %>%
      data.frame() %>%
      tibble::rownames_to_column("Sample_ID") %>%
      dplyr::inner_join(meta_data, by=c("Sample_ID"="Sample_ID")) %>%
      dplyr::select(Sample_ID, all_of(plot_genes), Time, Status)
  }
  
  # Add split_by column to expr_df to define groups in order to calculate multiple_cutoff
  if (!is.na(survival_params$split_by)){
    expr_df <- expr_df %>% 
      dplyr::left_join(meta_data %>% dplyr::select(Sample_ID, survival_params$split_by),
                       by=c("Sample_ID"="Sample_ID"))
  }
  
  return(expr_df)
}

# NOTE:  Output of calc_cutoffs is list(df,ls)
# If plotting by Sex, make sure to create column "model" based on which lines will be split
# calc_cutoffs <- function(df, gene, group, survival_params){
#   
#   # Identify upper & lower cutoffs based on stratify_criteria
#   #*************************Split samples by median**************************#
#   if(survival_params$stratify_criteria == "m"){
#     quartiles <- stats::quantile(x = df[[gene]],
#                                  probs = c(0, 0.25, 0.50, 0.75, 1),
#                                  na.rm=TRUE)
#     
#     cutoff_lower_end <- quartiles[[3]]
#     cutoff_upper_end <- quartiles[[3]]
#     cutoff_middle <- "NA"
#   }
#   
#   #****************Split samples into top and bottom tertiles****************#
#   else if(survival_params$stratify_criteria == "t"){
#     tertiles <- stats::quantile(x = df[[gene]],
#                                 probs = c(0, 0.33, 0.66, 1),
#                                 na.rm=TRUE)
#     
#     cutoff_lower_end <- tertiles[[2]]
#     cutoff_upper_end <- tertiles[[3]]
#     cutoff_middle <- "NA"
#   }
#   
#   #***************Split samples into top and bottom quartiles****************#
#   else if(survival_params$stratify_criteria == "q"){
#     quartiles <- stats::quantile(x = df[[gene]], 
#                                  probs = c(0, 0.25, 0.50, 0.75, 1),
#                                  na.rm=TRUE)
#     
#     cutoff_lower_end <- quartiles[[2]]
#     cutoff_upper_end <- quartiles[[4]]
#     cutoff_middle <- quartiles[[3]]
#   }
#   
#   #*********************Split expression range by thirds*********************#
#   else if(survival_params$stratify_criteria == "th"){
#     quartiles <- stats::quantile(x = df[[gene]], 
#                                  probs = c(0, 0.25, 0.50, 0.75, 1),
#                                  na.rm=TRUE)
#     iqr <- stats::IQR(x = df[[gene]],
#                       na.rm=TRUE)
#     
#     # Normal range of expression values lie between cutoff_lower & cutoff_upper
#     cutoff_upper <- quartiles[[4]]+1.5*iqr
#     cutoff_lower <- dplyr::if_else(quartiles[[1]]-1.5*iqr > 0, quartiles[[1]]-1.5*iqr, 0)
#     
#     # Based on normal range of expression, identify onethird & twothird cutoff
#     cutoff_lower_end <- cutoff_lower + (cutoff_upper-cutoff_lower)/3
#     cutoff_upper_end <- cutoff_lower + (cutoff_upper-cutoff_lower)*2/3
#     cutoff_middle <- "NA"
#   }
#   
#   #***************Split expression range using optimal cutoff****************#
#   else if(survival_params$stratify_criteria == "o"){
#     
#     # Sometimes quartiles will look like: 
#     # 0%       25%      50%      75%     100% 
#     # 0.000000 0.000000 0.000000 0.000000 3.495493 
#     # In such cases, surv_cutpoint() will fail. So, we add extra if() here.
#     quartiles <- stats::quantile(x = df[[gene]], 
#                                  probs = c(0, 0.25, 0.50, 0.75, 1),
#                                  na.rm=TRUE)
#     
#     if (quartiles[[4]] > quartiles[[2]]){
#       res.cut <- survminer::surv_cutpoint(data = df,
#                                           time = "Time",
#                                           event = "Status",
#                                           variables = gene)
#       
#       cutoff_lower_end <- res.cut$cutpoint$cutpoint
#       cutoff_upper_end <- res.cut$cutpoint$cutpoint
#       cutoff_middle <- "NA"
#     } else{
#       #cat("Surv cutpoint unable to detect optimum cutoff")
#       cutoff_lower_end <- "NA"
#       cutoff_upper_end <- "NA"
#       cutoff_middle <- "NA"
#     }
#   }
#   
#   # Categorize the sample based on above cutoffs
#   if (survival_params$plot_all_bins == TRUE & survival_params$stratify_criteria == "q"){
#     df <- df %>% 
#       dplyr::mutate(Expression = dplyr::case_when(get(gene) > cutoff_upper_end ~ "HIGH",
#                                                   get(gene) <= cutoff_lower_end ~ "LOW",
#                                                   get(gene) <= cutoff_middle ~ "MED_LOW",
#                                                   TRUE ~ "MED_HIGH"))
#     
#   } else if (survival_params$plot_all_bins == TRUE) {
#     df <- df %>% 
#       dplyr::mutate(Expression = dplyr::case_when(get(gene) > cutoff_upper_end ~ "HIGH",
#                                                   get(gene) <= cutoff_lower_end ~ "LOW",
#                                                   TRUE ~ "MID"))
#     
#   } else if (survival_params$stratify_criteria == "none") {
#     #When plotting by Sex, Treatment response, we dont use expression data.
#     df <- df %>% 
#       dplyr::mutate(Expression = model)
#     cutoff_lower_end <- NA
#     cutoff_upper_end <- NA
#     cutoff_middle <- NA
#     
#   } else {
#     df <- df %>% 
#       dplyr::mutate(Expression = dplyr::case_when(get(gene) > cutoff_upper_end ~ "HIGH",
#                                                   get(gene) <= cutoff_lower_end ~ "LOW",
#                                                   TRUE ~ "MID")) %>%
#       dplyr::filter(Expression != "MID")
#   }
#   
#   # # Print the cutoffs
#   # cat("\nGene         :", gene)
#   # cat("\nGroup        :", group)
#   # cat("\nLower cutoff :", cutoff_lower_end)
#   # cat("\nUpper cutoff :", cutoff_upper_end)
#   # cat("\nMiddle cutoff:", cutoff_middle)
#   
#   # Create a list to store cutoff values
#   ls <- list("group" = c(), 
#              "gene" = c(), 
#              "lower" = c(), 
#              "upper" = c(), 
#              "middle" = c())
#   
#   ls$group <- c(group)
#   ls$gene <- c(gene)
#   ls$lower <- c(cutoff_lower_end)
#   ls$upper <- c(cutoff_upper_end)
#   ls$middle <- c(cutoff_middle)
#   
#   # Return the df and the cutoffs
#   return(list(df, ls))
# }

# NOTE:  Output of calc_surv_stats is list. 
# It also generate survival plot with risk table
calc_surv_stats <- function(df, group, prefix, output_dir){
  
  # If all samples belong to one group (like LOW or HIGH or males or female),
  # then quit the function as comparison cannot be done
  if (nrow(df %>% dplyr::count(model)) > 1){
    
    # Create a survival object where Alive = 0, Dead = 1
    surv_object <- survival::Surv(time = df$Time,
                                  event = df$Status,
                                  type = "right",
                                  origin = 0)
    
    # Create a formula for plotting survival curve
    surv_formula <- surv_object ~ model
    
    # Create a fit for survival curve.
    # NOTE: survival::survfit() gives error in ggsurvplot(). Use survminer::surv_fit()
    surv_curve <- survminer::surv_fit(formula = surv_formula,
                                      data = df,
                                      type = "kaplan-meier",
                                      group.by = NULL,
                                      match.fd = FALSE)
    
    # # Check summary of the survival curve with time duration of our interest
    # cat("\nRange of survival (months):", range(df$Time, na.rm=TRUE), "\n")
    # base::summary(surv_curve, times = base::seq(from = floor(range(df$Time, na.rm=TRUE)[[1]]),
    #                                                 to = ceiling(range(df$Time, na.rm=TRUE)[[2]]),
    #                                                 by = 3))
    
    # Create a Cox model for the survival curve and calculate stats
    cox_model <- survival::coxph(formula = surv_formula,
                                 data = df)
    #print(summary(cox_model))
    cat("\n")
    
    # Calculate HR, 95% CI for HR, p-val
    # NOTE: Variable mentioned in summary(cox_model) is numerator in h1(t)/h0(t).
    # The reference variable h0(t) will not be mentioned in co-efficients.
    # Make sure this is not the reference level i.e. low expression. If this is
    # the reference, then we need to reverse the HR ratio, legend labels
    #print(names(cox_model$coefficients))  
    
    # If samples belong to more than 2 groups (like LOW, MID, HIGH), then we 
    # cannot have survival stats. So, we set them to 0.
    if (nrow(df %>% dplyr::count(model)) == 2){
      # Store HR and CI
      if (stringr::str_detect(names(cox_model$coefficients), survival_params$reference)){
        HR <- round(exp(-cox_model$coefficients[[1]]), 2)
        CI <- round(exp(-confint(cox_model)), 2)
        CI_1 <- CI[1]
        CI[1] <- CI[2]
        CI[2] <- CI_1
      } else {
        HR <- round(exp(cox_model$coefficients[[1]]),2)
        CI <- round(exp(confint(cox_model)),2)
      }
      
      # Store pvalues by different methods
      pvals <- c()
      for (test.method in c("survdiff", "1", "n", "sqrtN", "S1","S2", "FH_p=1_q=1")){
        
        # Some of the methods give error. So, we catch them and skip
        if (sum(str_detect(string = class(tryCatch(survminer::surv_pvalue(fit = surv_curve,
                                                                          method = test.method,
                                                                          test.for.trend = FALSE,
                                                                          combine = FALSE), error = function(e) e)),
                           pattern = "error")) == 0){
          p_val <- survminer::surv_pvalue(fit = surv_curve,
                                          method = test.method,
                                          test.for.trend = FALSE,
                                          combine = FALSE)
          pvals <- c(pvals, p_val[[2]])
        } else{
          pvals <- c(pvals, 1)
          print(test.method)
        } 
      } 
    } else {
      HR <- 0
      CI <- c(0, 0)
      pvals <- c(0,0,0,0,0,0,0) # 7 pvalues for each method
    }
    
    # Plot the survival curve using survminer::ggsurvplot() instead of base::plot()
    # ggsurvplot() produces a list of ggplot objects: a survival curve and a risk table
    # Saving it using cowplot() first and then using ggsave() works nicely as
    # compared to saving directly using ggsave()
    if(survival_params$plot_curve == TRUE){
      
      # Plot the survival curve
      legend_label <- df %>% 
        dplyr::count(model) %>% 
        dplyr::select(model) %>% 
        unlist(.,use.names=FALSE)
      
      # We identify proper breaks based on max duration of the dataset
      # We want a maximum of 10 timepoint intervals that are multiples of 12
      max_time <- max(df$Time,na.rm=TRUE)
      n <- floor(max_time/10/12)*12
      if(max_time %/% n <= 10){
        breaks <- n
      } else{
        breaks <- n+12
      }
      
      surv_plot <- survminer::ggsurvplot(fit = surv_curve,
                                         pval = FALSE,
                                         palette = survival_params$color_palette,
                                         linetype = "solid",
                                         size = 1.5,                       # thickness of line
                                         
                                         # Format the legend
                                         legend  = "top",                  # position of legend
                                         legend.title = survival_params$legend_title,
                                         legend.labs = survival_params$legend_label,
                                         
                                         # Format the axes
                                         break.time.by = breaks,           # break X axis in time intervals of 12 months
                                         xlab = "Time (Months)",           # customize X axis label
                                         ylab = "Survival Probability",    # customize Y axis label
                                         title = dplyr::if_else(gene == "combined.exp", "", gene),
                                         
                                         # Format confidence intervals
                                         conf.int = survival_params$conf_interval,
                                         #conf.int.fill = ?,               # color to fill confidence interval
                                         conf.int.style = "ribbon",        # confidence interval style
                                         conf.int.alpha = 0.3,             # confidence fill color transparency
                                         
                                         # Format the risk table
                                         risk.table = survival_params$plot_risk_table,
                                         risk.table.title = "Number at risk",
                                         risk.table.y.text.col = TRUE,     # color of risk table text annotations
                                         risk.table.pos = "out",           # draw risk table outside survival plot
                                         
                                         # Format the censor points
                                         censor = TRUE,
                                         censor.shape = '|',
                                         censor.size = 5)
      
      surv_plot$table <- surv_plot$table + 
        coord_cartesian(x=c(0,ceiling(max_time/breaks)*breaks), clip = "off")
      
      surv_plot$plot <- surv_plot$plot + 
        coord_cartesian(x=c(0,ceiling(max_time/breaks)*breaks), clip = "off")
      
      # Plot p and HR value
      method_plot <- "log-rank"
      p_plot <- pvals[1]  
      
      survplot_stats_grob <- grobTree(textGrob(label = paste0("p = ", formatC(p_plot, format = "e", digits = 1),
                                                "\nHR = ", round(HR,1), " [", round(CI[1],1), ", ", round(CI[2],1), "]",
                                                "\nMethod = ", method_plot),
                                 x = 0.50,
                                 y = 0.90,
                                 hjust = 0,
                                 gp = grid::gpar(fontfamily="Times", fontface="bold", col="black", fontsize=10)))
      
      # Add p values and HR values to plot
      surv_plot$plot <- surv_plot$plot %++%
        ggplot2::annotation_custom(survplot_stats_grob)
      
      cowplot::plot_grid(plotlist = surv_plot,
                         align = "hv",
                         axis = "tblr",
                         nrow = 2,  
                         ncol = 1, 
                         rel_widths = 1,
                         rel_heights = c(1,0.45),
                         labels = NULL,
                         label_size = 14,
                         label_fontfamily = NULL,
                         label_fontface = "bold",
                         label_colour = NULL)
      
      f_name <- paste0(prefix, "_", group, "_", survival_params$stratify_criteria, ".pdf")
      f_name <- gsub("/", "-", x=f_name)
      
      # Save the plot
      ggplot2::ggsave(filename = f_name,
                      plot = last_plot(),
                      device = "pdf",
                      path = output_dir,
                      width = 7,
                      height = 7,
                      units = c("in"),
                      dpi = 300,
                      limitsize = TRUE,
                      bg = NULL)
    }
  } else {
    HR <- 0
    CI <- c(0, 0)
    pvals <- c(0,0,0,0,0,0,0) # 7 pvalues for each method
  }
  
  # Create a list to store survival stats
  ls <- list("group" = c(), 
             "HR" = c(), 
             "CI_lower" = c(), 
             "CI_upper" = c(), 
             "pvalue" =c(), 
             "logrank" = c(), 
             "reg_logrank.late" = c(), 
             "Gehan_Breslow.early" = c(),
             "Tarone_Ware.early" = c(), 
             "Peto_Peto.early" = c(),  
             "modified_Peto_Peto" = c(), 
             "Fleming_Harrington" = c())
  
  ls$group               <- c(group)
  ls$HR                  <- c(HR)
  ls$CI_lower            <- c(CI[1])
  ls$CI_upper            <- c(CI[2])
  ls$logrank             <- c(pvals[1])
  ls$reg_logrank.late    <- c(pvals[2])
  ls$Gehan_Breslow.early <- c(pvals[3])
  ls$Tarone_Ware.early   <- c(pvals[4])
  ls$Peto_Peto.early     <- c(pvals[5])
  ls$modified_Peto_Peto  <- c(pvals[6])
  ls$Fleming_Harrington  <-c(pvals[7])
  
  return(ls)
}

# NOTE: Output of plot_survival is list(df,ls)
plot_survival <- function(expr_df, gene, survival_params, prefix, output_dir){
  
  # Create an empty dataframe to store expr_df and classification from calculate_cutoffs()
  survival_data <- data.frame(model = " ")
  
  # Create a list to store results of calculate_cutoffs, surv_plot et
  stats <- list("gene" = c(), 
                "group" = c(),
                "lower_cutoff" = c(),
                "middle_cutoff" = c(),
                "upper_cutoff" = c(),
                "HR" = c(), 
                "CI_lower" = c(), 
                "CI_upper" = c(), 
                "logrank" = c(), 
                "reg_logrank.late" = c(), 
                "Gehan_Breslow.early" = c(),
                "Tarone_Ware.early" = c(), 
                "Peto_Peto.early" = c(),  
                "modified_Peto_Peto" = c(), 
                "Fleming_Harrington" = c())
  
  # Create a list of groups for multiple_cutoff calculation 
  if (survival_params$multiple_cutoff == TRUE & !is.na(survival_params$split_by)){
    cutoff_groups <- expr_df %>% 
      dplyr::add_count(get(survival_params$split_by)) %>%
      dplyr::filter(n>2) %>%
      dplyr::select(all_of(survival_params$split_by)) %>% 
      unlist(use.names=FALSE) %>% 
      unique()
  } else if (survival_params$multiple_cutoff == TRUE & is.na(survival_params$split_by)){
    print("Please define survival_params$split_by to calculate multiple cutoffs")
  } else {
    cutoff_groups <- c(NA)
  }
  
  # STEP 1: Calculate cutoffs
  # If cutoffs need to be calculated for each group, subset the expr_df and pass
  # it to calculate_cutoffs(). Else, pass entire expr_df to calculate_cutoffs()
  for (group in cutoff_groups){
    
    # Subset the expr_df for each group to calculate cutoffs
    if (survival_params$multiple_cutoff == TRUE & !is.na(survival_params$split_by)){
      df <- expr_df %>% dplyr::filter(get(survival_params$split_by) == group)
    } else if (survival_params$multiple_cutoff == TRUE & is.na(survival_params$split_by)){
      cat("\n 'split_by' variable is undefined but 'multiple_cutoff' is set to TRUE")
    } else{
      df <- expr_df
    }
    
    # Calculate cutoffs for each group
    mat <- calc_cutoffs(df, gene, group, survival_params)
    
    ##### Save the data from output of calculate_cutoffs()
    survival_data       <- dplyr::bind_rows(survival_data, mat[[1]])
    stats$gene          <- c(stats$gene,          mat[[2]]$gene)
    #stats$group         <- c(stats$group,         mat[[2]]$group)
    stats$lower_cutoff  <- c(stats$lower_cutoff,  mat[[2]]$lower)
    stats$middle_cutoff <- c(stats$middle_cutoff, mat[[2]]$middle)
    stats$upper_cutoff  <- c(stats$upper_cutoff,  mat[[2]]$upper)
  }
  
  # Populate the model variable by concatenating "Expression" and "split_by"
  if (!is.na(survival_params$split_by)){
    survival_data <- survival_data %>%
      dplyr::mutate(model = paste0(Expression, "_", get(survival_params$split_by))) %>%
      dplyr::filter(!is.na(Sample_ID))
  } else {
    survival_data <- survival_data %>%
      dplyr::mutate(model = Expression) %>%
      dplyr::filter(!is.na(Sample_ID))
  }
  
  # Create a list of groups for plotting survival curves 
  if (!is.na(survival_params$split_by)){
    plot_groups <- expr_df %>% 
      dplyr::add_count(get(survival_params$split_by)) %>%
      dplyr::filter(n>2) %>%
      dplyr::select(all_of(survival_params$split_by)) %>% 
      unlist(use.names=FALSE) %>% 
      unique()
  } else {
    plot_groups <- c(NA)
  }
  
  # STEP 2: Calculate survival stats
  # If each group has to be plotted in separate plots, subset the survival_data
  # and pass it to calc_surv_stats(). Else, pass entire survival_data to 
  # calc_surv_stats().
  for (group in plot_groups){
    
    # Subset the survival_data for each group to generate separate plots
    if (survival_params$split_plot == TRUE & !is.na(survival_params$split_by)){
      df <- survival_data %>% dplyr::filter(get(survival_params$split_by) == group)
    } else if (survival_params$split_plot == TRUE & is.na(survival_params$split_by)){
      cat("\n 'split_by' variable is undefined but 'split_plot' is set to TRUE")
    } else{
      df <- survival_data
    }
    
    # Calculate survival stats for each group
    cox_stats <- calc_surv_stats(df, group, prefix, output_dir)
    
    ##### Save the data from output of calc_surv_stats()
    stats$group               <- c(stats$group,               cox_stats$group)
    stats$HR                  <- c(stats$HR,                  cox_stats$HR)
    stats$CI_lower            <- c(stats$CI_lower,            cox_stats$CI_lower)
    stats$CI_upper            <- c(stats$CI_upper,            cox_stats$CI_upper)
    stats$logrank             <- c(stats$logrank,             cox_stats$logrank)
    stats$reg_logrank.late    <- c(stats$reg_logrank.late,    cox_stats$reg_logrank.late)
    stats$Gehan_Breslow.early <- c(stats$Gehan_Breslow.early, cox_stats$Gehan_Breslow.early)
    stats$Tarone_Ware.early   <- c(stats$Tarone_Ware.early,   cox_stats$Tarone_Ware.early)
    stats$Peto_Peto.early     <- c(stats$Peto_Peto.early,     cox_stats$Peto_Peto.early)
    stats$modified_Peto_Peto  <- c(stats$modified_Peto_Peto,  cox_stats$modified_Peto_Peto)
    stats$Fleming_Harrington  <- c(stats$Fleming_Harrington,  cox_stats$Fleming_Harrington)
  }
  return(list(survival_data, stats))
}

# We need to calculate a combined expression value of all genes in the gene
# signature for each sample. Next, we use surv_cutpoint() to classify samples 
# into high and low groups and plot survival curves.
# There are several approaches to calculate the combined expression value. While
# normal z-scaling is logical, https://doi.org/10.1186/gb-2006-7-10-r93 
# recommends a more advanced z-scaling which is implemented below. Refer
# "Measuring gene set expression" in "Materials & stratify_criteria" section.

# ---- CRISPR FUNCTIONS ----

prep_crispr_guides <- function(data, output_dir){
  
  # Ensure required columns are present
  required_cols <- c("ID", "Seq", "Gene")
  
  if (!all(required_cols %in% colnames(data))) {
    stop("Input data must contain columns: ID, Seq, Gene")
  }
  
  # Define a simple reverse complement function
  reverse_complement <- function(seq) {
    seq_chars <- unlist(strsplit(seq, split = ""))
    comp <- c(A = "T", T = "A", G = "C", C = "G")
    rc <- rev(comp[seq_chars])
    paste(rc, collapse = "")
  }
  
  # Compute reverse complements
  data <- data %>%
    dplyr::mutate(Rev_Seq = vapply(Seq, reverse_complement, character(1)),
                  Mageck_format = paste(ID, Seq, Gene, sep = ","),
                  Mageck_format_rev = paste(ID, Rev_Seq, Gene, sep = ",")
    )
  
  # Write to Excel
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Guides")
  openxlsx::writeData(wb, sheet = "Guides", x = data, rowNames = FALSE)
  openxlsx::saveWorkbook(wb, 
                         file = file.path(output_dir, "CRISPR_Guides.xlsx"), 
                         overwrite = TRUE)
  
  message("CRISPR guide data saved to ", file.path(output_dir, "CRISPR_Guides.xlsx"))
}

# Dataframe as input
# Column "Gene" MUST be present
# Column Gene MUST have control sgRNAs labelled as "none" and/or "safe"
calc_t_score <- function(data){
  
  # Create a dataframe of control sgRNAs
  data_control <- data %>%
    dplyr::filter(Gene %in% c("none", "safe"))
  
  median_ctrl <- median(data_control$LFC, na.rm=TRUE)
  sd_ctrl <- sd(data_control$LFC, na.rm=TRUE)
  
  # Normalize to control sgRNAs
  data <- data %>%
    dplyr::mutate(pZ = (LFC-median_ctrl)/sd_ctrl)
  
  data_control <- data_control %>%
    dplyr::mutate(pZ = (LFC-median_ctrl)/sd_ctrl)
  
  U_ctrl <- median(data_control$pZ)
  Var_ctrl <- var(data_control$pZ)
  N_ctrl <- mean((data %>% dplyr::count(Gene))$n)
  # Nctrl is the average number of sgRNAs per gene in a given screen
  
  data <- data %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(U_gene = median(pZ),
                  Var_gene = var(pZ),
                  N_gene = n(),
                  U_ctrl = U_ctrl,
                  Var_ctrl = Var_ctrl,
                  N_ctrl = N_ctrl,
                  S_gene = (Var_gene*(N_gene-1)) + (Var_ctrl*(N_ctrl-1)),
                  t_score = (U_gene - U_ctrl)/sqrt(S_gene/N_gene + S_gene/N_ctrl),
                  Abs_t_score = abs(t_score)) %>%
    dplyr::select(Gene, U_gene, Var_gene, N_gene, U_ctrl, Var_ctrl, N_ctrl, S_gene, t_score, Abs_t_score) %>%
    dplyr::distinct_at("Gene", .keep_all = TRUE)
  
  return(data)
}

# data is a dataframe output of calc_t_score
plot_t_score <- function(data, disp_genes, suffix, save_path){
  
  y_cutoff <- sort(data$Abs_t_score, decreasing = TRUE)[100]
  xmin <- floor(min(data$U_gene))
  xmax <- ceiling(max(data$U_gene))
  ymin <- 0
  ymax <- max(data$Abs_t_score)
  
  color_breaks <- c(-20,0,20)
  p <- ggplot2::ggplot(data = data,
                       aes(x = U_gene, 
                           y = Abs_t_score,
                           size = Abs_t_score,
                           #color = pz,
                           fill = U_gene)) +
    # Plot dot plot
    ggplot2::geom_point(col="black", 
                        shape=21,
                        stroke=0.5,
                        position=position_jitter(h=0.01,w=0.01)) +
    # Define the theme of plot
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "U_gene") +
    coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax), clip = "off") +
    #scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
    #scale_y_continuous(breaks = seq(0, 5, by = 1)) +
    ggplot2::guides(size = "none",
                    fill = guide_colourbar(theme = theme(legend.key.width  = unit(0.75, "lines"),
                                                         legend.key.height = unit(10, "lines"),
                                                         legend.ticks = element_blank(),
                                                         legend.frame = element_rect(colour = "Black",
                                                                                     linewidth = 0.5)))) +
    # Define the color of the dots
    ggplot2::scale_fill_viridis_c(option="turbo", limits =c(-5,3))
  
  if (length(disp_genes) > 0){
    p <- p + ggrepel::geom_text_repel(data = data %>% dplyr::filter(Gene %in% disp_genes),
                                      mapping = aes(label = Gene),
                                      size = 5,
                                      show.legend = FALSE,
                                      direction = "both",   #"y"
                                      box.padding = 2.5,      # increases line length somehow
                                      point.padding = 0.1,  # distance around point = dist between line and point
                                      max.overlaps = nrow(data),
                                      position = position_quasirandom())
    
    
  }
  #geom_hline(yintercept= y_cutoff, linetype ="dotted")
  
  # scale_fill_gradientn(colors=c("#007ba7", "Black","#FFFF00"), 
  #                      limits=c(-20, 20), 
  #                      values=c(0, scales::rescale(color_breaks, from = range(color_breaks)), 1))
  #scale_fill_gradient2(low="#007ba7", mid="Black", high="Yellow", midpoint = 0, limits=c(-5, 2))
  #scale_fill_continuous_diverging(palette = "Tofino")
  
  ggplot2::ggsave(filename = paste0(suffix, ".jpg"),
                  plot = p,
                  device = "tiff",
                  path = save_path,
                  width = 7,
                  height = 7,
                  units = c("in"),
                  dpi = 300,
                  limitsize = TRUE,
                  bg = NULL)
}

# ---- VENN DIAGRAM ----

# Use this website for 3 set data if this script's venn diagram isnt pretty
# http://www.ehbio.com/test/venn/#/

plot_venn <- function(data, path, suffix){
  
  # Input validation
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (!dir.exists(path)) stop("`path` directory does not exist.")
  if (!is.character(suffix) || length(suffix) != 1) stop("`suffix` must be a single string.")
  
  plot_title <- suffix
  ncol <- ncol(data)
  
  # Clean column names (replace _ and . with space)
  colnames(data) <- stringr::str_replace_all(colnames(data), c("_" = " ", "\\." = " "))
  
  # Set cat.pos, cat.dist, cex and palette based on number of columns
  if (ncol == 4){
    pos <- c(330, 15, 330, 15)
    dist <- c(0.27, 0.25, 0.15, 0.13)
    cex = 2
    palette1 <- c("#C8E7F5", "#00008C", "#F6D2E0", "#E75480")        
  } else if (ncol == 3){
    pos <- c(0, 0, 180)
    dist <- c(0.1, 0.1, 0.1)
    cex = 2
    palette1 <- c("#C8E7F5", "#F6D2E0", "#db6d00")     
  } else if (ncol == 2){
    pos <- c(0, 0)
    dist <- c(0.05, 0.05)
    cex = 2.75
    palette1 <- c("#C8E7F5", "#db6d00")                             
  } else if (ncol == 1){
    pos <- c(0)
    dist <- c(0.1)
    cex = 2.75
    palette1 <- c("#F6D2E0")                                         
  } else {
    stop("`data` must have between 1 and 4 columns.")
  }
  
  # Create a dataframe to store the wrapped column names
  annotation <- data.frame(Labels = stringr::str_wrap(colnames(data), width = 10))
  
  # Convert the data frame to a named list (removing NAs)
  genes <- base::vector(mode = "list", length = ncol(data))
  names(genes) <- annotation$Labels
  
  for (i in 1:ncol(data)){
    
    # remove NA values and create a list of genes for each label
    genes[[i]] <- data[!is.na(data[i]),i]
  }
  
  # Plot the venn diagram
  VennDiagram::venn.diagram(x = genes,
                            main = plot_title, 
                            category.names = annotation$Labels,
                            filename = file.path(path, paste0("Venn_Diagram_", suffix, ".tiff")),
                            output = TRUE,
                            scaled = FALSE,
                            imagetype = "tiff",
                            height = 11, 
                            width = 11,
                            units = "in",
                            resolution = 600,
                            compression = "lzw",
                            margin = 0.3,    #amount of white space around Venn Diagram in grid units
                            
                            # Formatting the shapes of venn diagram
                            lwd = 1.5,                 #thickness of line
                            lty = 1,                   #type of line
                            col = "black",             #color of line
                            
                            # Formatting numbers inside venn diagram
                            cex = cex,                 #font size (2 or 2.75)
                            fontface = "bold",         #font style
                            fontfamily = "sans",       #font type
                            
                            # Formatting title of venn diagram
                            main.cex = 2,              #font size
                            main.fontface = "bold",    #font style
                            main.fontfamily = "sans",  #font type
                            main.col = "black",        #font color
                            
                            # Formatting category of venn diagram
                            cat.cex = 2,               #font size
                            cat.fontface = "bold",     #font style
                            cat.fontfamily = "sans",   #font type
                            cat.col = palette1,  #"black",
                            
                            # Formatting colors of venn diagram
                            fill = palette1,
                            alpha = rep(0.5, ncol), #0.5=50% transparency, 1=0% transparency
                            #cat.default.pos = "outer",    
                            
                            cat.pos = pos,    
                            cat.dist = dist, 
                            disable.logging = TRUE,
                            ext.text = TRUE)
  
  #******************************************************************************#
  #                          SAVE THE OVERLAPPING GENES                          #
  #******************************************************************************#
  
  # Save the list of overlapping genes. NOTE: You need to manually figure out
  # which genes belong to which overlap based on number of genes overlapping
  
  # Calculate overlaps
  overlap <- VennDiagram::calculate.overlap(x = genes)
  
  # Identify maximum number of genes present in any overlap
  max_len = max(lengths(overlap))
  
  # Create an dataframe of size length(overlap), max with NAs
  results = data.frame(matrix("", nrow = max_len, ncol = length(overlap)))
  rownames(results) <- paste0("Gene#", seq(max_len))
  colnames(results) <- paste0("Intersection#", seq(length(overlap)))
  
  # Populate the dataframe with gene names
  for (i in 1:length(overlap)){
    if (length(overlap[[i]]) > 0){
      for (j in 1:length(overlap[[i]])){
        results[j,i] <- overlap[[i]][j]
        #results[[j,i]] <- overlap[[i]][j]
      }
    }
  }
  
  # Save results to Excel
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Output")
  openxlsx::writeData(wb, sheet = "Output", x = results, rowNames = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "Input")
  openxlsx::writeData(wb, sheet = "Input", x = data, rowNames = FALSE)
  openxlsx::saveWorkbook(wb, file = file.path(path, paste0("Overlap_", suffix, ".xlsx")),
                         overwrite = TRUE,  returnValue = FALSE)
  
  invisible(NULL)
}

# ---- UPSET PLOT ----

# NOTE: The names of list will be on Y axis of bottom graph
# The intersection size based on Y values of top graph
# All possible intersections will be displayed in bottom graph
plot_upset <- function(listInput = NULL, selected_sets = NULL,
                       min_intersection_size = NULL, filename = "Upset_Plot",
                       output_dir = getwd()) {
  
  # Default example list if none provided
  if (is.null(listInput)) {
    listInput <- list(
      DMPK = c(1,2,3,4,5,6,7,8,9,11,16,17),
      MAPK1 = c(1,2,3,4,5,6,7,8,9,16,17,18),
      RAF1 = c(1,2,3,4,5,6,7,8,9,16,17,18),
      ROCK1 = c(1,2,3,4,5,6,7,8,9,17,18),
      PIM1 = c(1,2,3,4,5,6,7,8,9,17),
      DYRK2 = c(1,2,3,4,5,6,7,8,9,17),
      STK26 = c(1,2,3,4,5,6,7,8,9,17),
      MAP2K2 = c(1,2,3,4,5,6,7,8,9),
      LIMK1 = c(1,2,3,4,5,6,7,8,9),
      MYO3A = c(1,2,3,4,5,6,7,8,9),
      TSSK1B = c(1,2,3,4,5,6,7,8,9),
      HIPK4 = c(1,2,3,4,5,6,7,8,9),
      DCLK1 = c(1,2,3,4,5,6,7,8,9),
      NEK3 = c(1,2,3,4,5,6,7,8,9),
      CDKL1 = c(1,2,3,4,5,6,7,8,9),
      TRPM7 = c(1,2,3,4,5,6,7,8,10,11,12,13,14,15),
      MAPK7 = c(1,2,3,4,5,6,7,8,17),
      IKBKE = c(1,2,4,5,6,7,8,9,17),
      PRKAG3 = c(1,2,4,5,6,7,8,9),
      TAF1L = c(1,2,4,5),
      JAK2 = c(2,4,5,6,7,8,9,16,17,18)
    )
    
    selected_sets <- c("DCLK1", "MAPK1", "CDKL1", "ROCK1", "MAPK7", 
                       "HIPK4", "MAP2K2", "DYRK2")
    min_intersection_size <- 5
  }
 
  # Convert list to UpSetR input format
  upset_data <- UpSetR::fromList(listInput)
  
  # Create UpSet plot
  p <- UpSetR::upset(data = upset_data,
                     empty.intersections = "on",
                     cutoff = min_intersection_size,      # minimum intersection size to show
                     mb.ratio = c(0.5, 0.5),              # ratio between main bar and sets bar
                     sets = selected_sets,                # sets to display (order matters)
                     order.by = "freq",                   # order intersections by frequency
                     main.bar.color = "#1f78b4",          # nicer color for bars
                     sets.bar.color = "#33a02c",
                     #nintersects = 5,                    # number of groups on X axis
                     #nsets = 21,                         # number of groups on Y axis
                     text.scale = c(2, 2, 1.5, 1.5, 2, 1.5) # scale axis and text sizes
  )
  
  # Convert to ggplot object for saving
  ggplot_obj <- ggplotify::as.ggplot(p)
  
  # Save plot
  ggsave(filename = file.path(output_dir, paste0(filename, ".jpg")),
         plot = ggplot_obj,
         height = 11,
         width = 11
  )
  
  message("UpSet plot saved to: ", file.path(output_dir, paste0(filename, ".jpg")))
  
  return(ggplot_obj)
}

# ---- DEPRECATED FUNCTIONS ----

# DEPRECATED (used during Seurat v3)
v3_sctransform_singlecell <- function(filtered_seurat){
  
  # Seurat v5 stores counts of each sample in separate layers. Merge them.
  filtered_seurat@assays$RNA <- SeuratObject::JoinLayers(filtered_seurat@assays$RNA)
  
  # Split each sample into a seurat object to get a list of seurat object
  split.seurat <- Seurat::SplitObject(object = filtered_seurat,
                                      split.by = "Sample")
  
  # Remove samples with less than 50 cells so that RunPCA() doesnt give error
  split.seurat <- split.seurat[names(split.seurat)[sapply(split.seurat, ncol) > 50]]
  
  for (i in 1:length(split.seurat)){
    
    # Normalize the data
    split.seurat[[i]] <- Seurat::NormalizeData(object = split.seurat[[i]],
                                               assay = "RNA",
                                               normalization.method = "LogNormalize",
                                               scale.factor = 10000,
                                               margin = 1,
                                               verbose = FALSE)
    
    # Perform cell cycle scoring
    split.seurat[[i]]  <- Seurat::CellCycleScoring(object = split.seurat[[i]],
                                                   s.features = intersect(s_genes,rownames(split.seurat[[i]]@assays$RNA@features)),
                                                   g2m.features = intersect(g2m_genes, rownames(split.seurat[[i]]@assays$RNA@features)),
                                                   ctrl = NULL)
    
    split.seurat[[i]]$CC.Score <- split.seurat[[i]]$G2M.Score-split.seurat[[i]]$S.Score
    
    # SCTransform() is better than FindVariableFeatures() & ScaleData()
    # split.seurat[[i]] <- Seurat::FindVariableFeatures(object = split.seurat[[i]],
    #                                                   assay = "RNA",
    #                                                   selection.method = "vst",
    #                                                   nfeatures = 2000)
    # split.seurat[[i]] <- Seurat::ScaleData(object = split.seurat[[i]],
    #                                        features = NULL,
    #                                        assay = "RNA",
    #                                        vars.to.regress = NULL)
    
    # Perform scaling & variable feature identification usign SCTransform()
    split.seurat[[i]] <- Seurat::SCTransform(object =  split.seurat[[i]],
                                             assay = "RNA",
                                             new.assay.name = "SCT",
                                             do.correct.umi = TRUE,
                                             ncells = 5000,
                                             variable.features.n = 3000,
                                             vars.to.regress = c("CC.Score","MitoRatio"),
                                             do.scale = FALSE,
                                             do.center = TRUE,
                                             vst.flavor = "v2",
                                             return.only.var.genes = TRUE)
    
    
    # Remove ribosomal, Riken, predicted and mitochondrial genes from
    # VariableFeatures so that PCA, UMAP and hence clustering are not affected
    var_f <- split.seurat[[i]]@assays$SCT@var.features
    var_f <- var_f[!grepl(pattern = "^[Rr][Pp][SsLl]|R[Ii][Kk]$|^[Mm][Tt]-|^G[Mm][0-9.]+$", 
                          x=var_f)]
    
    split.seurat[[i]]@assays$SCT@var.features <- var_f
    cat("\nFinal number of variable features:", length(var_f), "\n")
    
    # Perform dimensional reduction using PCA on SCT assay variable features
    split.seurat[[i]] <- Seurat::RunPCA(object = split.seurat[[i]],
                                        assay = "SCT",
                                        features = NULL,
                                        ndims.print = 1,
                                        nfeatures.print = 1,
                                        reduction.name = "pca",
                                        reduction.key = "PC_",
                                        verbose = FALSE)
    
    # Perform dimensional reduction using UMAP on PCA dimensions
    split.seurat[[i]] <- Seurat::RunUMAP(object = split.seurat[[i]],
                                         dims = 1:40,
                                         reduction = "pca",
                                         reduction.name = "umap",
                                         reduction.key = "UMAP_",
                                         verbose = FALSE)
    
  }
  
  return(split.seurat)
}

# DEPRECATED (used during Seurat v3)
v3_integrate_singlecell <- function(sct, ref_samples){
  
  # NOTE: In v3, Harmony integration was not possible as FindIntegrationAnchors()
  # doesnt support reduction="harmony". Moreover, integration using rpca, cca, 
  # jpca couldnt be stored in same object since FindIntegrationAnchors() output
  # varies for each method
  
  #***STEP 7B: SELECT 3000 MOST VARIABLE GENES TO USE FOR INTEGRATING THE DATA***#
  integ_features <- Seurat::SelectIntegrationFeatures(object.list=split.seurat,
                                                      nfeatures=3000,
                                                      assay=NULL) #c("SCT", "SCT"),
  
  
  #******************STEP 7C: FIND RESIDUALS FOR MISSING GENES*******************#
  split.seurat <- Seurat::PrepSCTIntegration(object.list=split.seurat,
                                             assay="SCT",
                                             anchor.features=integ_features)
  
  #******STEP 7D: FIND COMMON ANCHORS BETWEEN SAMPLES TO INTEGRATE THE DATA******#
  integ_anchors <- Seurat::FindIntegrationAnchors(object.list=split.seurat,
                                                  reference=ref_samples,
                                                  anchor.features=integ_features,
                                                  scale=TRUE,
                                                  normalization.method="SCT",
                                                  sct.clip.range=NULL,
                                                  reduction="rpca", #"cca", "jpca", "rlsi"
                                                  l2.norm=TRUE,
                                                  dims=1:30)
  
  #******STEP 7E: FIND OPTIMUM k.weight FOR USE IN Seurat::IntegrateData()*****#
  # Find minimum anchors between 2 datasets
  kweight1 <- as.data.frame(integ_anchors@anchors) %>%
    dplyr::group_by(dataset1, dataset2) %>%
    distinct_at("cell1", .keep_all=TRUE) %>%
    dplyr::summarize(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::select(n) %>%
    unlist(use.names=FALSE) %>%
    min()
  
  # Find half of number of cells in sample with least cell count
  kweight2 <- filtered_seurat@meta.data %>%
    dplyr::count(Sample) %>%
    dplyr::filter(n >=50) %>%
    dplyr::select(n) %>%
    unlist(use.names=FALSE) %>%
    min()
  
  kweight2 <- floor(kweight2/2)
  
  kweight <- base::min(kweight1, kweight2, 100)
  dplyr::if_else(kweight >= 100, 100, kweight)
  cat("\n", celltype, "\tkweight1:", kweight1, "\tkweight2:", kweight2, "\tkweight:",  kweight, "\n")
  
  # NOTE: Integration will not fail anymore. If it fails, identify the 2
  # datasets that are involved in the error and use kweight=number of anchors
  # for these 2 datasets.
  cat("\nNumber of unique anchors between datasets\n")
  print(as.data.frame(integ_anchors@anchors) %>%
          dplyr::group_by(dataset1, dataset2) %>%
          distinct_at("cell1", .keep_all=TRUE) %>%
          dplyr::summarize(n=n()), n=1000)
  
  #************************STEP 7F: INTEGRATE THE DATA*************************#
  # NOTE: weight.reduction=NULL means new PCA will be calculated & used to
  # calculate anchor weights
  integrated_seurat.rpca <- Seurat::IntegrateData(anchorset=integ_anchors.rpca,
                                                  new.assay.name="integrated",
                                                  normalization.method="SCT",
                                                  features=NULL,
                                                  features.to.integrate=NULL,
                                                  dims=1:30,
                                                  k.weight=kweight, #default is 100
                                                  weight.reduction=NULL,
                                                  sd.weight=1)
  
  #**STEP 7G: RUN PCA USING 3000 INTEGRATION FEATURES & UMAP USING FIRST 40 PCs**#
  integrated_seurat <- Seurat::RunPCA(object=integrated_seurat,
                                      assay="integrated",
                                      features=NULL)
  
  integrated_seurat <- Seurat::RunUMAP(object=integrated_seurat,
                                       dims=1:40,
                                       reduction="pca")
  return(integrated_seurat)
}

# DEPRECATED (used during Seurat v3)
### Generate whitelist for CITESeq
# Input is filtered seurat object
# Output is a list of csv files - one per batch containing valid barcodes
v3_generate_whitelist <- function(filtered_seurat, output_dir){
  
  # Extract barcodes and split by "_"
  bc <- filtered_seurat@meta.data$Cell
  
  # Adjust this based on how your samples are named
  # NOTE: There will be multiple samples within each batch
  barcodes <- data.frame(stringr::str_split_fixed(bc, "_", 2)) %>%
    dplyr::rename(Batch = identity(1), Barcodes = identity(2)) %>%
    dplyr::mutate(Barcodes = stringr::str_replace(Barcodes, "-1", ""),
                  Batch = gsub(pattern="-GEX.*", replacement="", x=Batch))
  
  # Remove duplicate barcodes within each batch
  barcodes <- barcodes %>% 
    dplyr::group_by(Batch) %>% 
    dplyr::distinct_at("Barcodes", .keep_all=TRUE) %>% 
    as.data.frame()
  
  # Check how many barcodes are present in each batch
  barcodes %>% dplyr::group_by(Batch) %>% dplyr::count()
  
  # Save barcodes from each batch to individual csv files
  for (i in unique(barcodes$Batch)){
    whitelist <- barcodes %>%
      dplyr::filter(Batch == i) %>%
      dplyr::select(Barcodes)
    
    write.table(x = whitelist,
                file = paste0(scripts_path, proj, "_", i, "_whitelist.csv"),
                row.names = FALSE,
                col.names = FALSE)
  }
}

# DEPRECATED (used during Seurat v3)
# Input is path to folder containing h5ad
# Output is a raw seurat object
v3_read_h5ad <- function(input_path){
  
  # Load h5ad (useful if analyzing collaborator data in h5ad format)
  SeuratDisk::Convert(source = paste0(input_path, proj, ".h5ad"),
                      dest = "h5seurat",
                      assay="RNA",
                      overwrite = FALSE)
  
  raw_seurat <- SeuratDisk::LoadH5Seurat(file = paste0(input_path, proj, ".h5seurat"))
  
  return(raw_seurat)
}

### Import data from output of citeseq
# Input is path to demux_results folder
# Ouput is seurat object of each sample
v3_read_citeseq <- function(input_path){
  
  # Create a list of samples that have been demultiplexed already
  files <- list.files(path = paste0(input_path, "singlets/"),
                      full.names = FALSE)
  samples <- gsub(pattern="\\..*", replacement="", x=files)
  
  # Loop through each of the individual object in demux directory & import data
  for (i in 1:length(files)){
    
    # Read the seurat object containing demultiplexed singlets
    sample.seurat <- readRDS(file = paste0(demux_results, "singlets/", files[i]))
    
    # Assign the seurat object to its corresponding variable
    assign(samples[i], sample.seurat)
  }
}

SpatialFeaturePlotBlend <- function(cells_obj, column_1, column_2, combine=TRUE){
  
  # Convert decimal number to hexadecimal. Pad with 0s if only a single
  # character following conversion.
  as_hex <- function(num) {
    hex_str <- as.character(as.hexmode(num))
    if (nchar(hex_str) == 1) {
      hex_str <- paste0("0", hex_str)
    }
    
    return(hex_str)
  }
  
  metadata_to_hexadecimal <- function(in_dat) {
    apply(in_dat, 2,
          function(x) {
            # Make minimum 0
            x - min(x)
          }) %>%
      apply(2,
            function(x) {
              # Constrain to range [0, 255]
              round(255 * (x / max(x)))
            }) %>%
      apply(1,
            function(x) {
              # Convert to hexadecimal codes
              toupper(paste0("#", as_hex(x[1]), as_hex(x[2]), "00"))
            })
  }
  
  blend_plot_theme <- theme(legend.position="none",
                            plot.title=element_text(hjust=0.5))
  
  plot_list <- lapply(c(column_1, column_2),
                      function(column) {
                        max_color <- if_else(column == column_1,
                                             "#FF0000", "#00FF00")
                        SpatialFeaturePlot(cells_obj, column) +
                          scale_fill_gradient(low="#000000",
                                              high=max_color) +
                          ggtitle(column) +
                          blend_plot_theme
                      })
  
  dat <- FetchData(cells_obj, c(column_1, column_2))
  colors <- as.matrix(dat) %>% metadata_to_hexadecimal()
  
  new_md_column <- paste0(column_1, "_vs_", column_2)
  cells_obj[[new_md_column]] <- colors
  names(colors) <- as.character(colors)
  
  plot_list[[3]] <- SpatialDimPlot(cells_obj, new_md_column, cols=colors) +
    ggtitle(paste0(column_1, "_", column_2)) +
    blend_plot_theme
  
  side_length <- 100
  legend_grid <- expand.grid(seq(from=min(dat[, column_1]),
                                 to=max(dat[, column_1]),
                                 length.out=side_length),
                             seq(from=min(dat[, column_2]),
                                 to=max(dat[, column_2]),
                                 length.out=side_length))
  colnames(legend_grid) <- c(column_1, column_2)
  legend_colors <- metadata_to_hexadecimal(legend_grid)
  legend_grid$color <- legend_colors
  names(legend_colors) <- legend_colors
  
  legend <- ggplot(legend_grid,
                   aes(x=.data[[column_1]], y=.data[[column_2]],
                       color=color)) +
    geom_point(shape=15, size=1.9) +
    scale_color_manual(values=legend_colors) +
    coord_cartesian(expand=FALSE) +
    theme(legend.position="none", aspect.ratio=1,
          panel.background=element_blank())
  
  plot_list[[4]] <- wrap_plots(ggplot() + theme_void(), legend,
                               ggplot() + theme_void(), ncol=1,
                               heights=c(0.2, 0.6, 0.2))
  
  if (combine == FALSE) {
    return(plot_list)
  } else {
    p <- wrap_plots(plot_list, nrow=1,
                    widths=c(0.28, 0.28, 0.28, 0.16))
    return(p)
  }
}

# ---- NOTES ----

#******************************************************************************#
#                         TO SUBSET DATA OR NOT BEFORE DESEQ2                          
#******************************************************************************#

# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
# https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html

# NOTE: If any one group has high within-group variability in PCA plot, we 
# SHOULD exclude those samples by sub-setting before creating dds object and 
# calculating dispersion estimates. Else, use full dataset for modelling.

# NOTE: 
# One factor   : design ~ Condition
# Two factors  : design ~ Condition + Treatment + Condition:Treatment (OR) design ~ Condition*Treatment
# Three factors: design ~ Cell.Line + Condition + Treatment + Condition:Treatment +.. (OR) design ~ Cell.Line*Condition*Treatment
# Using * will include all possible interaction terms and is RECOMMENDED
# The results are the same irrespective of the order i.e 
# Cell.Line*Condition*Treatment gives same results as Cell.Line*Treatment*Condition

# NOTE: If you have Cell.Line, Treatment, Condition variables and want to 
# include them in design, ideally design MUST be ~ Cell.Line*Treatment*Condition.
# However, if one cell line doesnt have a specific treatment, then DESeq2 will
# throw "Full model martix is less than full rank" error. So, ALWAYS create a 
# new column "Comparisons" in meta_data and use design ~ Comparisons ALWAYS.
# Populate the "Comparisons" column in meta_data by concatenating as below:
# paste0(Cell.Line, Treatment, Condition, sep=".") 

# NOTE: contrast terms CANNOT start with numbers. So, if you have cell line name
# "22RV1" rename it as "RV1_22" etc

# NOTE: You MUST have atleast 2 replicates per group. Else, DESeq2 will throw
# "Error in base::colMeans(x, na.rm = na.rm, dims = dims, ...) : 'x' must be an
# array of at least two dimensions"

# for (n in 1:length(DEG.params$contrast)){
#   
#   if (DEG.params$deseq2.batch.correct == TRUE){
#     
#     # Perform DESeq2() using sva modelled surrogate variables
#     # Create DESeq2 object with surrogate variables in design
#     sva.formula_string <- formula_string
#     sva_dds <- batch_correct_sva(meta_data, read_data, sva.formula_string)
#     sva_dds <- run_deseq2(sva_dds, meta_data, DEG.params, n, "sva", degs_dir)
#     
#     # Perform DESeq2() using combat corrected counts
#     # Create DESeq2 object with appropiate variables in design
#     # Get combat corrected raw reads
#     read_data_combat <- batch_correct_combat(meta_data, read_data, combat.formula_string)
#     combat.formula_string <- formula_string
#     combat_dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data_combat,
#                                                  colData=meta_data, 
#                                                  design=~1)
#     design(combat_dds) <- as.formula(combat.formula_string)
#     dds <- run_deseq2(combat_dds, meta_data, DEG.params, n, "combat", degs_dir)
#   }
# }

# You have 3 groups A, B & C but want to perform DEG between group A and B only.
# Should you exclude group C samples and perform DEG analysis??
# NOTE: Subsetting data affects (i) sizeFactors (ii) Dispersion estimates and
# hence the final DEGs

# The links below explain when data needs to be subset.
# https://support.bioconductor.org/p/108471/
# https://support.bioconductor.org/p/81038/
# https://support.bioconductor.org/p/9156036/
# https://support.bioconductor.org/p/69341/

# The authors of DESeq2 and EdgeR recommend NOT to subset meta_data & 
# read_data but rather use contrasts to get DEGs between groups
# You should subset your data ONLY when one of the groups has large variance 
# (like group C) in link below
# https://master.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#if-i-have-multiple-groups-should-i-run-all-together-or-split-into-pairs-of-groups

# Dispersion of gene X in group A ~ Variance of gene X in group A/mean of gene X group A
# Dispersion of gene X in group B ~ Variance of gene X in group B/mean of gene X group B
# Dispersion of gene X in group C ~ Variance of gene X in group C/mean of gene X group C
# Dispersion of gene X for experiment ~ (Dispersion of gene X in groups A, B, C)
# Samples within Group A and B have low variance in overall expression profile 
# (clustered closely in PCA plot) but samples in group C have different
# expression profile. So, the final dispersion estimates of most genes will be 
# affected a lot by group C. In such a scenario, if we are only doing DEG 
# comparison between groups A & B, then we can subset only group A and B samples
# before doing DESeq2. 
# SO, CHECK PCA PLOT FOR EVERY EXPERIMENT TO DETERMINE IF DATA NEED TO BE SUBSET.

# meta_data <- meta %>% dplyr::filter(get(Variable) %in% c(Comparisons$Target[n], Comparisons$Reference[n]))
# corrected_read_data <- read %>% dplyr::select(rownames(meta_data))
# NOTE: sizefactors MUST be calculated ONLY using samples being compared for 
# differential expression. So, make sure read_data and meta_data ONLY have
# samples being compared for differential expression.

# dds <- DESeq2::DESeqDataSetFromMatrix(countData=corrected_read_data, colData=meta_data, design=~ Tissue + Age)
# dds <- DESeq(dds, test="LRT", reduced=~ Tissue)
# res <- DESeq2::results(object=dds)

# Study effect of Treatment ignoring effect of Sex
#dds <- DESeqDataSetFromMatrix(countData=corrected_read_data, colData=meta_data, design=~ Treatment)

# Study effect of Treatment after removing effect of Sex (assumes equal effect of Sex on different treatments)
#dds <- DESeqDataSetFromMatrix(countData=corrected_read_data, colData=meta_data, design=~ Sex + Treatment)

# Study effect of Sex on treatment?
# Study effect of Treatment after removing effect of Sex (assumes different effect of Sex on different treatments)
#dds <- DESeqDataSetFromMatrix(countData=corrected_read_data, colData=meta_data, design=~ Sex + Treatment + Sex:Treatment)

# #*******************************DIAGNOSTIC TESTS DESEQ2 *******************************#

# YET TO IMPLEMENT

# meta_data MUST be xlsx file, MUST have "Sample_ID" column with sample names
# without any duplication
# Example: If samples were isolated on different days & prepared using different
# kits, "Batch" column must have values like "1_Ribo", "1_Poly-A", "2_Ribo", etc
# NOTE: Make sure there are no white spaces in the Target and Reference columns
# in excel file. R will change any white space (" ") to dot ("."). So, replace 
# white space (" ") with underscore ("_") before importing into R.

# NOTE: The normalized counts you get from sva_dds and DESeq2 dds are NOT batch 
# corrected and they will be identical. Both DESeq2() and svaseq() ONLY model
# for batch factors/surrogate variables in the design formula, they DO NOT 
# modify counts like combatseq. The normalized counts you get from combatseq
# are batch corrected. 

# for (n in 1:length(Comparisons$Variable)){
#   
#   # This generates a new column "id" that has info on samples being comparared
#   meta_data_comp <- meta_data %>%
#     dplyr::mutate(id=get(Comparisons$Variable[n]))
#   
#   #combat_corrected_read_data <- combatseq_batch(read_data, meta_data_comp)
#   #sva_dds <- svaseq_batch(read_data, meta_data_comp)
#   
#   # Perform DESeq2() using in-built batch modelling
#   approach <- "DESeq2_modelled"
#   if (length(unique(meta_data_comp$Batch)) > 1){
#     dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
#                                           colData=meta_data_comp, 
#                                           design=~ Batch+id)
#   } else {
#     dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
#                                           colData=meta_data_comp, 
#                                           design=~ id)
#   }
#   dds <- run_deseq2(dds, meta_data_comp, annotations, Comparisons, n, approach, prefix, results_path)
#   deseq2_norm_counts(dds, annotations, approach, suffix) # batch corrected if you more than 1 batch
#   plot_qc(dds, meta_data_comp, approach, suffix)
#   
#   # Perform DESeq2() using combat corrected read counts
#   if (!identical(read_data, combat_corrected_read_data)){
#     approach <- "combat_corrected"
#     combat_dds <- DESeq2::DESeqDataSetFromMatrix(countData=combat_corrected_read_data,
#                                                  colData=meta_data_comp, 
#                                                  design=~ id)
#     combat_dds <- run_deseq2(combat_dds, meta_data_comp, annotations, Comparisons, n, lfc.cutoff, padj.cutoff, suffix)
#     combatseq_norm_counts(combat_dds, annotations, approach, suffix)  #combat batch corrected
#     plot_qc(combat_dds, meta_data_comp, approach, suffix)
#   }
#   
#   # Perform DESeq2() using sva modelled surrogate variables SV1 and SV2
#   approach <- "sva_modelled"
#   sva_dds <- run_deseq2(sva_dds, meta_data_comp, annotations, Comparisons, n, lfc.cutoff, padj.cutoff, approach, suffix)
#   # calc_norm_counts(sva_dds, annotations, approach, suffix)   # uncorrected
#   svaseq_norm_counts(sva_dds, annotations, approach, suffix)   # sva seq batch corrected
#   plot_qc(sva_dds, meta_data_comp, approach, suffix)
# }


# # (i) To view counts of specific gene across samples
# plotCounts(dds, gene=which.min(res$padj), intgroup=Variable)           # gene with lowest padj
# plotCounts(dds, gene=which.min(res$log2FoldChange), intgroup=Variable) # gene with lowest log2FC
# plotCounts(dds, gene=which.max(res$log2FoldChange), intgroup=Variable) # gene with highest log2FC
# plotCounts(dds, gene="ENSMUSG00000030598", intgroup=Variable)          # a specific gene
# 

# # To identify the genes interactively, run the 2 lines below. 
# # Then click on multiple dots and click Finish. A list of genes will be displayed 
# # idx <- identify(res$baseMean, res$log2FoldChange)
# # rownames(res)[idx]

# # (iv) Hierarchical clustering of samples using rld or vst
# rld_mat <- DESeq2::assay(rld)     # extract the matrix from a DESeq2 object
# rld_cor <- cor(x=rld_mat,       # compute pairwise correlation values
#                y=NULL,
#                use="everything",
#                method="pearson") 
# 
# # Check the output of cor(), make note of the rownames and colnames
# head(rld_cor)
# 
# vst_mat <- assay(vst)  
# vst_cor <- cor(vst_mat) 
# 
# # Check the output of cor(), make note of the rownames and colnames
# head(vst_cor)    
# 
# for (process in c("rld", "vst")){
#   
#   pheatmap::pheatmap(mat=get(paste0(process, "_cor")),
#                      color=colorRampPalette(rev(brewer.pal(n=11, name ="RdYlBu")))(100),
#                      breaks=NA, 
#                      border_color="white", #"grey60",
#                      cellwidth=NA, 
#                      cellheight=NA, 
#                      scale="none",   
#                      cluster_rows=TRUE,   #cluster the rows
#                      cluster_cols=TRUE,   #cluster the columns
#                      clustering_distance_rows="euclidean",
#                      clustering_distance_cols="euclidean",
#                      clustering_method="complete",
#                      legend=TRUE, 
#                      legend_breaks=NA,
#                      legend_labels=NA, 
#                      #annotation_row=,  
#                      #annotation_col=, 
#                      annotation_colors=dplyr::if_else(nrow(col_annotation)+nrow(row_annotation) > 0, ann_colors, NA),
#                      annotation_legend=TRUE,
#                      annotation_names_row=TRUE,
#                      annotation_names_col=TRUE,
#                      show_rownames=dplyr::if_else(nrow(mat)<80, TRUE, FALSE, missing=NULL), 
#                      show_colnames=dplyr::if_else(ncol(mat)<30, TRUE, FALSE, missing=NULL),
#                      fontsize=8, 
#                      fontsize_row=8, 
#                      fontsize_col=8,
#                      angle_col=c("270", "0", "45", "90", "315"),
#                      fontsize_number=0.8*fontsize, 
#                      #labels_row=display_row,
#                      #labels_col=display_col,
#                      filename=paste0(diagnostics_path, "Diagnostic_Correlation_Heatmap_using_", process, "_", celltype, ".pdf"))
# }
# 
# # (v) Checking if mean < variance (for NB model) or Mean=Variance (for Poisson model). 
# # Each point is a gene denoted by (x,y) where 
# # x=mean_count of gene across all samples & y=variance_count of gene across all samples
# 
# mean_counts <- apply(read_data[,], 1, mean)
# variance_counts <- apply(read_data[,], 1, var)
# df <- data.frame(mean_counts, variance_counts)
# ggplot(df) +
#   geom_point(aes(x=mean_counts, y=variance_counts)) +
#   geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
#   scale_y_log10() +
#   scale_x_log10()
# 
# ggplot2::ggsave(filename=paste0("Diagnostic_Scatter_plot_", celltype, ".pdf"),
#                 plot=last_plot(),
#                 device="pdf",
#                 path=diagnostics_path,
#                 scale=1,
#                 #width=8.5,
#                 #height=11,
#                 units=c("in"),
#                 dpi=300,
#                 limitsize=TRUE,
#                 bg="white")
# 
# # # (vii) Plot p-value histogram
# # hist(res$pvalue, col="lightblue", breaks=20)

# # # (ix) Plot a histogram for one of the samples to see how the counts are distributed. Adjust xlim if needed
# # ggplot(read_data, aes(x=results.S10_R1_trimmed.fastq.gz.csv)) +
# #   geom_histogram(stat="bin", bins=200) +
# #   xlim(-5,1000) +
# #   xlab("Raw expression counts") +
# #   ylab("Number of genes") +
# #   scale_color_brewer(palette="Dark2")
# # 
# # # (x) Extract counts, size factors, etc.. from dds
# # dds <- estimateSizeFactors(dds)         #redundant if you already ran dds <- DESeq(dds)
# # counts <- counts(dds)                   #counts[i,j]=raw_count of gene i in sample j
# # sizefactors <- sizeFactors(dds)         #sizefactors[j]=median (counts[,j]/geometric_mean(counts[i,]))
# # colSums(counts(dds))                    #Total number of raw counts per sample
# # colSums(counts(dds, normalized=TRUE))   #Total number of normalized counts per sample

# ---- CODE TESTING ---- 
# 
# DEG.params  <- list(Variable    = c("Comparisons"),
#                     Target      = c("ARCaPM.4Gy.NDRG1_mut"),
#                     Reference   = c("ARCaPM.0Gy.NDRG1_mut"),
#                     contrast    = c("ARCaPM.NDRG1_mut.4Gy-ARCaPM.NDRG1_mut.0Gy"),
#                     lfc.cutoff  = 0,
#                     padj.cutoff = 0.1,
#                     design      = "Cell.Line*Treatment*Condition", #"Comparisons", #"0+Condition:Treatment"
#                     design.ref  = c("Condition:WT", "Treatment:0Gy", "Comparisons:ARCaPM.0Gy.WT"),
#                     deseq2.batch.correct = FALSE,
#                     proj        = "RNASeq_Manish_22RV1_ARCaPM",
#                     species     = "Homo sapiens")
# 
# meta_data <- meta_data %>% filter(Cell.Line != "22RV1")
# n <- 1
# 
# 
# dds.new <- run_deseq2(dds, meta_data, DEG.params, n, approach, data_dir)
# design(dds) <- as.formula(~Comparisons)
# dds.old <- run_deseq2_old(dds, meta_data, DEG.params, n, approach, data_dir)
# 
# 
# # The following 3 approaches give identical results
# 
# # Method 1 => Similar way of defining contrasts like method2. Easy to compare 
# # samples but difference of difference not possible
# design <- "0+Condition:Treatment" 
# dds.test <- dds
# contrast1 <- c("Condition", "NDRG1_mut.Treatment4Gy", "NDRG1_mut.Treatment0Gy")
# contrast2 <- c("Condition", "NDRG1_mut.Treatment4Gy", "WT.Treatment0Gy")
# contrast3 <- c("Condition", "NDRG1_mut.Treatment4Gy", "WT.Treatment4Gy")
# contrast4 <- c("Condition", "NDRG1_mut.Treatment0Gy", "WT.Treatment4Gy")
# contrast5 <- c("Condition", "NDRG1_mut.Treatment0Gy", "WT.Treatment0Gy")
# contrast6 <- c("Condition", "WT.Treatment4Gy", "WT.Treatment0Gy")
# res1 <- DESeq2::results(dds, contrast = contrast1)
# res2 <- DESeq2::results(dds, contrast = contrast2)
# res3 <- DESeq2::results(dds, contrast = contrast3)
# res4 <- DESeq2::results(dds, contrast = contrast4)
# res5 <- DESeq2::results(dds, contrast = contrast5)
# res6 <- DESeq2::results(dds, contrast = contrast6)
# 
# # Method 2= > combine COndition and Treatment to new column Comparisons. 
# # Similar way of defining contrasts like method1. Easy to compare 
# # samples but difference of difference not possible
# design <- "Comparisons"
# contrast1A <- c("Comparisons", "ARCaPM.4Gy.NDRG1_mut", "ARCaPM.0Gy.NDRG1_mut")
# contrast2A <- c("Comparisons", "ARCaPM.4Gy.NDRG1_mut", "ARCaPM.0Gy.WT")
# contrast3A <- c("Comparisons", "ARCaPM.4Gy.NDRG1_mut", "ARCaPM.4Gy.WT")
# contrast4A <- c("Comparisons", "ARCaPM.0Gy.NDRG1_mut", "ARCaPM.4Gy.WT")
# contrast5A <- c("Comparisons", "ARCaPM.0Gy.NDRG1_mut", "ARCaPM.0Gy.WT")
# contrast6A <- c("Comparisons", "ARCaPM.4Gy.WT", "ARCaPM.0Gy.WT")
# res1A <- DESeq2::results(dds, contrast = contrast1A)
# res2A <- DESeq2::results(dds, contrast = contrast2A)
# res3A <- DESeq2::results(dds, contrast = contrast3A)
# res4A <- DESeq2::results(dds, contrast = contrast4A)
# res5A <- DESeq2::results(dds, contrast = contrast5A)
# res6A <- DESeq2::results(dds, contrast = contrast6A)
# 
# # Method 3 => conventional design, contrast needs to calculated for each sample
# # difference of difference is easy
# design <- "Condition+Treatment+Condition:Treatment"
# dds.standard <- dds
# mod_mat <- model.matrix(design(dds), colData(dds))
# NDRG1_0Gy <- colMeans(mod_mat[dds$Condition == "NDRG1_mut" & dds$Treatment == "0Gy",])
# NDRG1_4Gy <- colMeans(mod_mat[dds$Condition == "NDRG1_mut" & dds$Treatment == "4Gy",])
# WT_0Gy <- colMeans(mod_mat[dds$Condition == "WT" & dds$Treatment == "0Gy",])
# WT_4Gy <- colMeans(mod_mat[dds$Condition == "WT" & dds$Treatment == "4Gy",])
# 
# contrast1B <- NDRG1_4Gy - NDRG1_0Gy 
# contrast2B <- NDRG1_4Gy - WT_0Gy
# contrast3B <- NDRG1_4Gy - WT_4Gy
# contrast4B <- NDRG1_0Gy - WT_4Gy 
# contrast5B <- NDRG1_0Gy - WT_0Gy
# contrast6B <- WT_4Gy - WT_0Gy
# contrast7B <- (NDRG1_4Gy - NDRG1_0Gy) - (WT_4Gy - WT_0Gy)
# contrast8B <- 
#   
#   res1B. <- DESeq2::results(dds, contrast = contrast1B)
# res2B. <- DESeq2::results(dds, contrast = contrast2B)
# res3B. <- DESeq2::results(dds, contrast = contrast3B)
# res4B. <- DESeq2::results(dds, contrast = contrast4B)
# res5B. <- DESeq2::results(dds, contrast = contrast5B)
# res6B. <- DESeq2::results(dds, contrast = contrast6B)
# res7B. <- DESeq2::results(dds, contrast = contrast7B)
# 
# 
# # Check for differences in results
# df1 <- data.frame(res1) %>% dplyr::mutate_all(~(round(.,2)))
# df1A <- data.frame(res1A) %>% dplyr::mutate_all(~(round(.,2)))
# df1B <- data.frame(res1B) %>% dplyr::mutate_all(~(round(.,2)))
# 
# df2 <- data.frame(res2) %>% dplyr::mutate_all(~(round(.,2)))
# df2A <- data.frame(res2A) %>% dplyr::mutate_all(~(round(.,2)))
# df2B <- data.frame(res2B) %>% dplyr::mutate_all(~(round(.,2)))
# 
# df3 <- data.frame(res3) %>% dplyr::mutate_all(~(round(.,2)))
# df3A <- data.frame(res3A) %>% dplyr::mutate_all(~(round(.,2)))
# df3B <- data.frame(res3B) %>% dplyr::mutate_all(~(round(.,2)))
# 
# df4 <- data.frame(res4) %>% dplyr::mutate_all(~(round(.,2)))
# df4A <- data.frame(res4A) %>% dplyr::mutate_all(~(round(.,2)))
# df4B <- data.frame(res4B) %>% dplyr::mutate_all(~(round(.,2)))
# 
# df5 <- data.frame(res5) %>% dplyr::mutate_all(~(round(.,2)))
# df5A <- data.frame(res5A) %>% dplyr::mutate_all(~(round(.,2)))
# df5B <- data.frame(res5B) %>% dplyr::mutate_all(~(round(.,2)))
# 
# df6 <- data.frame(res6) %>% dplyr::mutate_all(~(round(.,2)))
# df6A <- data.frame(res6A) %>% dplyr::mutate_all(~(round(.,2)))
# df6B <- data.frame(res6B) %>% dplyr::mutate_all(~(round(.,2)))
# 
# # Display only rows and columns that are different
# dplyr::full_join(df1[rowSums(df1-df1A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df1A[rowSums(df1-df1A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df2[rowSums(df2-df2A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df2A[rowSums(df2-df2A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df3[rowSums(df3-df3A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df3A[rowSums(df3-df3A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df4[rowSums(df4-df4A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df4A[rowSums(df4-df4A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df5[rowSums(df5-df5A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df5A[rowSums(df5-df5A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df6[rowSums(df6-df6A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df6A[rowSums(df6-df6A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df1[rowSums(df1-df1B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df1B[rowSums(df1-df1B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df2[rowSums(df2-df2B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df2B[rowSums(df2-df2B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df3[rowSums(df3-df3B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df3B[rowSums(df3-df3B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df4[rowSums(df4-df4B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df4B[rowSums(df4-df4B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df5[rowSums(df5-df5B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df5B[rowSums(df5-df5B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# dplyr::full_join(df6[rowSums(df6-df6B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
#                  df6B[rowSums(df6-df6B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
#                  by=c("SYMBOL"="SYMBOL")) %>%
#   dplyr::transmute(SYMBOL              = SYMBOL,
#                    baseMean.diff       = baseMean.x-baseMean.y,
#                    log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
#                    lfcSE.diff          = lfcSE.x-lfcSE.y,
#                    stat.diff           = stat.x-stat.y,
#                    pvalue.diff         = pvalue.x-pvalue.y,
#                    padj.diff           = padj.x-padj.y)
# 
# 
# 
# res1 <- DESeq2::lfcShrink(dds=dds.test, res=res1, type="ashr")
# res1B <- DESeq2::lfcShrink(dds=dds.standard, res=res1B, type="ashr")
# res2 <- DESeq2::lfcShrink(dds=dds.test, res=res2, type="ashr")
# res2B <- DESeq2::lfcShrink(dds=dds.standard, res=res2B, type="ashr")
# res3 <- DESeq2::lfcShrink(dds=dds.test, res=res3, type="ashr")
# res3B <- DESeq2::lfcShrink(dds=dds.standard, res=res3B, type="ashr")
# res4 <- DESeq2::lfcShrink(dds=dds.test, res=res4, type="ashr")
# res4B <- DESeq2::lfcShrink(dds=dds.standard, res=res4B, type="ashr")
# res5 <- DESeq2::lfcShrink(dds=dds.test, res=res5, type="ashr")
# res5B <- DESeq2::lfcShrink(dds=dds.standard, res=res5B, type="ashr")
# res6 <- DESeq2::lfcShrink(dds=dds.test, res=res6, type="ashr")
# res6B <- DESeq2::lfcShrink(dds=dds.standard, res=res6B, type="ashr")
# res7B <- DESeq2::lfcShrink(dds=dds.standard, res=res7B, type="ashr")
# 
# 
# summary(res1)
# summary(res1B)
# summary(res1B.)
# summary(res2)
# summary(res2B)
# summary(res2B.)
# summary(res3)
# summary(res3B)
# summary(res3B.)
# summary(res4)
# summary(res4B)
# summary(res4B.)
# summary(res5)
# summary(res5B)
# summary(res5B.)
# summary(res6)
# summary(res6B)
# summary(res6B.)
# summary(res7B)
# summary(res7B.)
# res1B. <- DESeq2::lfcShrink(dds=dds, res=res1B., type="ashr")
# res2B. <- DESeq2::lfcShrink(dds=dds, res=res2B., type="ashr")
# res3B. <- DESeq2::lfcShrink(dds=dds, res=res3B., type="ashr")
# res4B. <- DESeq2::lfcShrink(dds=dds, res=res4B., type="ashr")
# res5B. <- DESeq2::lfcShrink(dds=dds, res=res5B., type="ashr")
# res6B. <- DESeq2::lfcShrink(dds=dds, res=res6B., type="ashr")
# res7B. <- DESeq2::lfcShrink(dds=dds, res=res7B., type="ashr")
# 
# setdiff(rownames(data.frame(res1B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
#         rownames(data.frame(res1B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res1B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
#         rownames(data.frame(res1B.) %>% filter(log2FoldChange < -0.58, padj < 0.1)))
# 
# setdiff(rownames(data.frame(res2B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
#         rownames(data.frame(res2B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res2B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
#         rownames(data.frame(res2B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 
# 
# 
# setdiff(rownames(data.frame(res3B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
#         rownames(data.frame(res3B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res3B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
#         rownames(data.frame(res3B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 
# 
# 
# setdiff(rownames(data.frame(res4B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
#         rownames(data.frame(res4B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res4B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
#         rownames(data.frame(res4B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res5B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
#         rownames(data.frame(res5B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res5B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
#         rownames(data.frame(res5B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 
# 
# 
# setdiff(rownames(data.frame(res6B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
#         rownames(data.frame(res6B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res6B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
#         rownames(data.frame(res6B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res7B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
#         rownames(data.frame(res7B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 
# 
# setdiff(rownames(data.frame(res7B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
#         rownames(data.frame(res7B.) %>% filter(log2FoldChange < -0.58, padj < 0.1)))