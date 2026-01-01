#!/usr/bin/env Rscript

# Read and store variables from command line interface (CLI)
cli <- base::commandArgs(trailingOnly = TRUE) 
args <- base::strsplit(x = cli, split = "=", fixed = TRUE)

for (e in args){
  argname <- e[1]
  argval <- e[2]
  assign(argname, argval)
}

# ==============================================================================
# 🧬 BULK RNA-SEQ WORKFLOW
# ==============================================================================

# ---- ⚙️️ Project Setup ---- 

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

# ---- 🧾 Load Metadata ----

metadata_xlsx <- file.path(proj.params$proj_dir, paste0(proj, "_Metadata.xlsx"))
metadata <- openxlsx::read.xlsx(xlsxFile = metadata_xlsx)

# ---- 📥 Load raw counts ----

proj_name       <- proj.params$proj
proj_dir        <- proj.params$proj_dir
counts_dir      <- proj.params$counts_dir
raw_counts_xlsx <- file.path(proj.params$proj_dir, paste0(proj, "_Raw_counts.xlsx"))

if (dir.exists(counts_dir) && length(list.files(counts_dir)) > 0){
  
  # Compile raw counts from  txt files
  raw_counts_mat <- merge_counts(counts_dir = counts_dir,
                                 filename   = proj_name, 
                                 output_dir = proj_dir)
} else if (file.exists(raw_counts_xlsx)){
  
  # Read raw counts from excel file
  raw_counts_df <- openxlsx::read.xlsx(xlsxFile = raw_counts_xlsx)
  
  if (!"SYMBOL" %in% colnames(raw_counts_df)) {
    log_warn(sample = "",
             step   = "RNASeq_Workflow.R",
             msg    = "1st column of xlsx file should be named 'SYMBOL'")
  } else {
    raw_counts_mat <- raw_counts_df %>% 
      tibble::column_to_rownames("SYMBOL")
  }

} else {
  
  log_error(sample = "",
            step = "RNASeq_Workflow.R",
            msg = glue::glue("Please provide raw counts in excel file '{raw_counts_xlsx}' or as txt files in '{counts_dir}'"))
}

# ---- 🛠️ Prepare DESeq2 Input ----

design      <- proj.params$deseq2$design
deseq2_data <- prepare_deseq2_input(expr_mat = raw_counts_mat,
                                    metadata = metadata,
                                    design   = design)

# ---- 🔍 QC : Initial PCA & Batch Inspection ----

proj_name <- proj.params$proj
proj_dir  <- proj.params$proj_dir

plot_pca(expr_mat    = deseq2_data$expr_mat, 
         metadata    = deseq2_data$metadata,
         filename    = proj_name,
         output_dir  = proj_dir,
         top_n_genes = 500,
         perform_vst = TRUE, 
         skip_plot   = FALSE)

# ---- 🧹 Sample Filtering : Outlier Removal ----

remove_samples <- NULL
# remove_samples <- c("SBQuadFc2", "SBQuadFc4")

raw_counts_mat <- raw_counts_mat[, !(colnames(raw_counts_mat) %in% remove_samples), drop = FALSE]


# ==============================================================================
# 🔁 DIFFERENTIAL EXPRESSION PER CONTRAST
# ==============================================================================

contrasts <- proj.params$deseq2$contrasts
design    <- proj.params$deseq2$design
proj_dir  <- proj.params$proj_dir

for (contrast in contrasts) {
  
  log_info(sample = contrast,
           step   = "RNASeq_Workflow.R",
           msg    = glue::glue("Running DESeq2 for contrast: '{contrast}'"))
  
  contrast_dir <- file.path(proj_dir, contrast)
  
  # ---- 🧬 Differential Expression (DESeq2) ----
  deseq2_results <- run_deseq2(expr_mat    = deseq2_data$expr_mat, 
                               metadata    = deseq2_data$metadata,
                               design      = design, 
                               contrast    = contrast,
                               output_dir  = contrast_dir,
                               lfc_cutoff  = 0, 
                               padj_cutoff = 0.1)
  
  # ---- 🏹 Visualization : MA Plot ----
  plot_ma(dds        = deseq2_results$dds,
          filename   = paste(contrast, proj_name, sep = "_"), 
          output_dir = contrast_dir)
  
  # ---- 🌋 Visualization : Volcano Plot ----
  plot_volcano(res_df      = deseq2_results$degs, 
               filename    = paste(contrast, proj_name, sep = "_"), 
               output_dir  = contrast_dir, 
               contrast    = contrast,
               label_genes = proj.params$volcano$label_genes,
               top_n       = proj.params$volcano$top_n,
               lfc_cutoff  = proj.params$volcano$lfc_cutoff, 
               padj_cutoff = proj.params$volcano$padj_cutoff)
  
  # ---- 🔥 Visualization : Heatmap ----
  
  # Identify samples for current contrast
  samples <- metadata %>%
    dplyr::filter(Comparisons %in% all.vars(expr = as.formula(paste0("~", contrast)))) %>%
    dplyr::pull(Sample_ID) %>%
    as.character() %>%
    base::intersect(colnames(deseq2_results$vst))
  
  # Identify significant DEGs
  sig_genes <- deseq2_results$degs %>% 
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::pull(SYMBOL) %>%
    base::intersect(rownames(deseq2_results$vst))
  
  # Plot heatmap
  ph <- plot_heatmap(expr_mat            = deseq2_results$vst[sig_genes, samples], 
                     label_genes         = NULL,
                     filename            = paste(contrast, proj_name, sep = "_"),
                     output_dir          = contrast_dir,
                     metadata_col        = metadata, 
                     metadata_row        = NULL,
                     col_annotations     = proj.params$heatmap$col_annotations,
                     row_annotations     = proj.params$heatmap$row_annotations,
                     col_gap_by          = proj.params$heatmap$col_gap_by,
                     row_gap_by          = proj.params$heatmap$row_gap_by,
                     col_cluster_by      = proj.params$heatmap$col_cluster_by,
                     row_cluster_by      = proj.params$heatmap$row_cluster_by,
                     plot_title          = proj.params$heatmap$plot_title,
                     heatmap_palette     = proj.params$heatmap$heatmap_palette,
                     annotation_palette  = proj.params$heatmap$annotation_palette,
                     border_color        = proj.params$heatmap$border_color,
                     force_log           = proj.params$heatmap$force_log,
                     show_expr_legend    = proj.params$heatmap$show_expr_legend,
                     save_plot           = proj.params$heatmap$save_plot,
                     save_matrix         = proj.params$heatmap$save_matrix)
  
  # ---- 🛤️ Pathway Analysis (GSEA & ORA) ----
  
  species <- proj.params$species
  gmt_dir <- proj.params$gmt_dir
  pathway_dir <- file.path(proj_dir, contrast, "Pathway_Analysis")
  pathway_results <- analyze_pathway(res_df     = deseq2_results$degs,
                                     species    = species, 
                                     gmt_dir    = gmt_dir,
                                     output_dir = pathway_dir,
                                     minsize    = 15, 
                                     maxsize    = 500)
  
  # ---- 🌿 Visualization : Pathway Enrichment ----
  
  # Identify top 10 Up & top 10 Down GSEA pathways for each collection
  top_gsea <- pathway_results$consensus %>%
    dplyr::filter(method != "ORA") %>%
    # Deduplicate by method priority
    dplyr::group_by(Collection, Consensus, Description) %>%
    dplyr::slice_min(order_by = match(method, c("FGSEA", "GSEA", "ORA")), n = 1) %>%
    dplyr::ungroup() %>%
    # Rank based on abs(NES) for each direction
    dplyr::group_by(Collection, Consensus) %>%
    dplyr::slice_max(order_by = abs(NES), n = 10, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # Identify top 10 Up & top 10 Down ORA pathways for each collection
  top_ora <- pathway_results$consensus %>%
    dplyr::filter(method == "ORA") %>%
    # Rank based on padj for each direction
    dplyr::group_by(Collection, Consensus) %>%
    dplyr::slice_min(order_by = padj, n = 10, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # Identify samples for current contrast
  samples <- metadata %>%
    dplyr::filter(Comparisons %in% all.vars(expr = as.formula(paste0("~", contrast)))) %>%
    dplyr::pull(Sample_ID) %>%
    as.character() %>%
    base::intersect(colnames(deseq2_results$vst))
  
  pathway_dir <- file.path(proj_dir, contrast, "Pathway_Analysis")
  plot_pathway(pathway_df = top_gsea, 
               method     = "GSEA",
               expr_mat   = deseq2_results$vst[, samples], 
               metadata   = metadata,
               output_dir = pathway_dir)
  
  plot_pathway(pathway_df = top_ora, 
               method     = "ORA",
               expr_mat   = deseq2_results$vst[, samples], 
               metadata   = metadata,
               output_dir = pathway_dir)
 

  # ---- 📡 Regulatory Network Analysis ----
  
  # Identify samples for current contrast
  samples <- metadata %>%
    dplyr::filter(Comparisons %in% all.vars(expr = as.formula(paste0("~", contrast)))) %>%
    dplyr::pull(Sample_ID) %>%
    as.character() %>%
    base::intersect(colnames(deseq2_results$vst))
  
  species <- proj.params$species
  tf_activity_samples <- analyze_tf(expr_mat = deseq2_results$vst[, samples], 
                                    res_df   = NULL, 
                                    species  = species, 
                                    stats    = c("ulm", "mlm", "viper"), 
                                    minsize  = 5, 
                                    top_n    = 500) 
  
  tf_activity_degs <- analyze_tf(expr_mat = NULL, 
                                 res_df   = deseq2_results$degs, 
                                 species  = species, 
                                 stats    = c("ulm", "mlm", "viper"), 
                                 minsize  = 5, 
                                 top_n    = 500)
  
  # ---- 🌿 Visualization : Regulatory Networks ----
   
  # Identify samples for current contrast
  samples <- metadata %>%
    dplyr::filter(Comparisons %in% all.vars(expr = as.formula(paste0("~", contrast)))) %>%
    dplyr::pull(Sample_ID) %>%
    as.character() %>%
    base::intersect(colnames(deseq2_results$vst))
  
  tf_dir <- file.path(proj_dir, contrast, "Regulatory_Network_Analysis")
  
  # Plot Bar, Heatmap for top 20 Up and Down TFs
  contrast <- contrasts[n]
  output_dir <- proj.params$tf_dir[n]
 
  plot_tf(tf_df = tf_activity_samples, 
          contrast = contrast, 
          metadata = metadata, 
          samples, 
          output_dir = tf_dir,
          n_tfs      = 20)
  plot_tf(tf_res_counts, contrast, metadata, samples, output_dir, n_tfs = 20)
  
}

# ==============================================================================
# 📉 GLOBAL QUALITY & EXPLORATION
# ==============================================================================

# ---- 🧬 QC : Dispersion Estimates ----
plot_dispersion(dds        = deseq2_results$dds,
                filename   = proj_name,
                output_dir = proj_dir)

# ---- 🔥 Visualization : Heatmap (Global) ----

# Extract VST Counts (blind)
# NOTE: vst counts are affected by design ONLY when blind=FALSE
vsd_blind <- DESeq2::vst(object = deseq2_results$dds, 
                         blind = TRUE)
vst_counts_blind <- SummarizedExperiment::assay(vsd_blind) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  add_annotation() %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarize(across(.cols = where(is.numeric), .fns = mean, na.rm = TRUE), .groups = "drop") %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("SYMBOL") %>%
  dplyr::select(-dplyr::matches("ENSEMBL|ENTREZ|GENEID")) %>%
  as.matrix()

# Perform LRT
# 'reduced = ~1' tests if the 'design(dds)' significantly improves the model over no groups at all
dds <- deseq2_results$dds
dds_LRT <- DESeq2::DESeq(dds, test = "LRT", reduced = ~1)
res_LRT <- DESeq2::results(dds_LRT) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  add_annotation()

# Identify significant DEGs
sig_genes <- res_LRT %>% 
  dplyr::filter(padj <= 0.05) %>% 
  dplyr::pull(SYMBOL) %>%
  base::intersect(rownames(vst_counts_blind))

# Plot heatmap
ph <- plot_heatmap(expr_mat            = vst_counts_blind[sig_genes, , drop = FALSE], 
                   label_genes         = NULL,
                   filename            = proj_name,
                   output_dir          = proj_dir,
                   metadata_col        = metadata, 
                   metadata_row        = NULL,
                   col_annotations     = proj.params$heatmap$col_annotations,
                   row_annotations     = proj.params$heatmap$row_annotations,
                   col_gap_by          = proj.params$heatmap$col_gap_by,
                   row_gap_by          = proj.params$heatmap$row_gap_by,
                   col_cluster_by      = proj.params$heatmap$col_cluster_by,
                   row_cluster_by      = proj.params$heatmap$row_cluster_by,
                   plot_title          = proj.params$heatmap$plot_title,
                   heatmap_palette     = proj.params$heatmap$heatmap_palette,
                   annotation_palette  = proj.params$heatmap$annotation_palette,
                   border_color        = proj.params$heatmap$border_color,
                   force_log           = proj.params$heatmap$force_log,
                   show_expr_legend    = proj.params$heatmap$show_expr_legend,
                   save_plot           = proj.params$heatmap$save_plot,
                   save_matrix         = proj.params$heatmap$save_matrix)

# ---- Save Normalized Counts (Global) ----

# NOTE: Normalized counts are primarily used for individual gene visualization 
# (box/violin plots) and "ground truth" reporting. NEVER use them for any other
# analysis.
# For PCA, Heatmaps, Sample-level TF activity), always use VST/rlog. 
# For differential expression, GSEA, and ORA, use the statistical results (res_df). 

norm_counts_DESeq2 <- function(metadata, read_data, proj.params)
  
  
  
# ---- THE END ----