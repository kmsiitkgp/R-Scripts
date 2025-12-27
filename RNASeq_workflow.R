#!/usr/bin/env Rscript

# Read and store variables from command line interface (CLI)
cli <- base::commandArgs(trailingOnly = TRUE) 
args <- base::strsplit(x = cli, split = "=", fixed = TRUE)

for (e in args){
  argname <- e[1]
  argval <- e[2]
  assign(argname, argval)
}

# ---- PROJECT SET UP ----

if (.Platform$OS.type == "windows") {
  parent_dir  <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data"
  gmt_dir     <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GSEA_genesets"
  scripts_dir <- NULL
  script_file <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/Custom_Functions.R"
} else {  # Linux/macOS (e.g., HPC)
  parent_dir  <- "/hpc/home/kailasamms/scratch"
  gmt_dir     <- "/hpc/home/kailasamms/projects/GSEA_genesets"
  scripts_dir <- "/hpc/home/kailasamms/projects/scRNASeq"
  script_file <- "/hpc/home/kailasamms/projects/scRNASeq/Custom_Functions.R"
}

if (file.exists(script_file)) {
  source(script_file)
} else{
  stop(paste("Custom_Functions.R not found at:", script_file))
}

# Define these 4 variables in project specific script
# contrasts <- c()
# deseq2.override <- list()
# heatmap.override <- list()
# volcano.override <- list()

proj.params <- setup_project(proj             = proj,
                             species          = species,  #"Mus musculus", "Homo sapiens"
                             contrasts        = contrasts,
                             parent_dir       = parent_dir,
                             gmt_dir          = gmt_dir,
                             scripts_dir      = scripts_dir,
                             deseq2.override  = deseq2.override,
                             heatmap.override = heatmap.override,
                             volcano.override = volcano.override)

# ---- BULK RNA SEQ WORKFLOW ----

# Import metadata
metadata <- read.xlsx(file.path(proj.params$proj_dir, paste0(proj, "_Metadata.xlsx")))

# Compile raw counts if available
raw_counts_mat <- merge_counts(counts_dir = proj.params$counts_dir,
                               filename   = proj.params$proj, 
                               output_dir = proj.params$proj_dir)

# Import raw counts
raw_counts_mat <- read.xlsx(file.path(proj.params$proj_dir, paste0(proj, "_Raw_counts.xlsx"))) %>%
  tibble::column_to_rownames("SYMBOL")

# Reformat raw_counts_mat and metadata for DESeq2
deseq2_data <- prepare_deseq2_input(expr_mat = raw_counts_mat,
                                    metadata = metadata,
                                    design   = proj.params$deseq2$design)

# Visualization : Plot PCA
plot_pca(expr_mat    = deseq2_data$expr_mat, 
         metadata    = deseq2_data$metadata, 
         top_n_genes = 500,
         perform_vst = TRUE, 
         skip_plot   = FALSE,
         filename    = proj.params$proj,
         output_dir  = proj.params$proj_dir)

# Remove bad samples based on PCA [CASE by CASE BASIS]
# remove_samples <- c("SBQuadFc2", "SBQuadFc4")
remove_samples <- NULL
raw_counts_mat <- raw_counts_mat[, !(colnames(raw_counts_mat) %in% remove_samples), drop = FALSE]

# Run DESeq2 on each contrast
contrasts <- proj.params$deseq2$contrast

for (contrast in contrasts) {
  
  # ---- Differential Expression (DESeq2) ----
  
  output_dir  <- file.path(proj.params$proj_dir, contrast)
  deseq2_results <- run_deseq2(expr_mat    = deseq2_data$expr_mat, 
                               metadata    = deseq2_data$metadata,
                               design      = proj.params$deseq2$design, 
                               contrast    = contrast,
                               output_dir  = output_dir,
                               lfc.cutoff  = 0, 
                               padj.cutoff = 0.1)
  
  # ---- Visualization : MA Plot ----
  
  output_dir <- file.path(proj.params$proj_dir, contrast, "DEG_Analysis")
  plot_ma(dds        = deseq2_results$dds, 
          output_dir = output_dir)
  
  # ---- Visualization : Volcano Plot ----
  
  output_dir <- file.path(proj.params$proj_dir, contrast, "DEG_Analysis")
  plot_volcano(DEGs_df    = deseq2_results$degs, 
               contrast   = contrast, 
               output_dir = output_dir)


}
