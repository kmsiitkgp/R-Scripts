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

# ---- Helper functions ----
save_xlsx <- function(data, file, sheet_name, row_names) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = sheet_name)
  openxlsx::writeData(wb, sheet = sheet_name, x = data, rowNames = row_names)
  if(row_names){
    openxlsx::writeData(wb, sheet = sheet_name, x = "SYMBOL", startCol = 1, startRow = 1)
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

get_samples <- function(metadata, contrast) {
  metadata %>%
    dplyr::filter(Comparisons %in% stringr::str_split(contrast, "-")[[1]]) %>%
    dplyr::pull(Sample.ID) %>%
    as.character()
}

trim_heatmap_whitespace <- function(gtable) {
  gtable$widths <- grid::unit.pmax(gtable$widths, grid::unit(0.1, "cm"))
  gtable$heights <- grid::unit.pmax(gtable$heights, grid::unit(0.1, "cm"))
  gtable
}

# ---- Pre-analysis ----

meta_data <- read.xlsx(file.path(proj.dir, paste0(proj.params$proj, ".Metadata.xlsx")))
read_data <- NULL

# Compile raw counts if not provided
if (is.null(read_data)) {
  if (dir.exists(proj.params$counts.dir)) {
    read_data <- merge_counts(proj.dir, proj.params)
  } else {
    stop("The specified count directory does not exist: ", proj.params$counts.dir)
  }
}

# Format metadata and readdata prior to DESeq2 analysis 
meta_data <- prep_metadata(meta_data, read_data)
read_data <- prep_readdata(meta_data, read_data)
l <- check_data(meta_data, read_data, proj.params)
meta_data <- l[[1]]
read_data <- l[[2]]

# PCA plot for entire dataset
plotPCA_DESeq2(meta_data, read_data, proj.dir)

# ---- Main analysis ----

# Remove bad samples [CASE by CASE BASIS]
meta_data <- read.xlsx(file.path(proj.dir, paste0(proj.params$proj, ".Metadata.xlsx")))
read_data <- read.xlsx(file.path(proj.dir, paste0(proj.params$proj, ".raw.counts.xlsx")))
# read_data <- read_data %>%
#   dplyr::select(everything(), -c("SBQuadFc2", "SBQuadFc4"))

analyze_data <- function(proj.dir, meta_data, read_data, proj.params) {
  
  # Prepare input
  meta_data <- prep_metadata(meta_data, read_data)
  read_data <- prep_readdata(meta_data, read_data)
  l <- check_data(meta_data, read_data, proj.params)
  meta_data <- l[[1]]
  read_data <- l[[2]]
  
  # Loop through contrasts
  for (n in seq_along(proj.params$contrast)) {
    
    contrast <- proj.params$contrast[n]
    save.dir <- file.path(proj.dir, contrast)
    if (!dir.exists(save.dir)) dir.create(save.dir, recursive = TRUE)
    
    # ---- Create DESeq2 object ----
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = read_data,
                                          colData = meta_data, 
                                          design = ~1)
    design(dds) <- as.formula(paste0("~", proj.params$deseq2.design))
    
    # ---- Run DESeq2 ----
    dds.list <- run_deseq2(dds, contrast, proj.params)
    
    # ---- Annotate DEGs ----
    DEGs_df <- dds.list$res %>% 
      as.data.frame() %>%
      tibble::rownames_to_column("ID") %>% 
      add_annotation() %>%
      dplyr::filter(!is.na(SYMBOL)) %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::summarize(across(.cols = where(is.numeric), .fns = mean), .groups = "drop")
    
    save_xlsx(DEGs_df, file.path(save.dir, "DEGs.xlsx"), "DEGs", row_names = FALSE)
    
    # ---- VST (non-blind) counts ----
    vsd <- DESeq2::vst(dds.list$dds, blind = FALSE)
    vst_counts <- SummarizedExperiment::assay(vsd) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("ID") %>%
      add_annotation() %>%
      dplyr::filter(!is.na(SYMBOL)) %>%
      dplyr::group_by(SYMBOL) %>%
      dplyr::summarize(across(.cols = where(is.numeric), mean), .groups = "drop") %>%
      tibble::column_to_rownames("SYMBOL") %>%
      dplyr::select(-dplyr::starts_with("ENSEMBL"), -dplyr::starts_with("ENTREZ")) %>%
      as.matrix()
    
    save_xlsx(vst_counts, file.path(save.dir, "VST_counts.xlsx"), "VST_Nonblind", row_names = TRUE)
    
    # ---- Plots: MA, Volcano, Heatmap ----
    plotMA_DESeq2(dds.list$dds, contrast, proj.dir)
    plot_volcano(DEGs_df, contrast, proj.dir)
    
    samples <- get_samples(meta_data, contrast)
    sig_genes <- DEGs_df %>% dplyr::filter(padj <= 0.05) %>% dplyr::pull(SYMBOL)
    norm_counts <- vst_counts[sig_genes, samples, drop = FALSE]
    
    metadata_col <- meta_data %>% dplyr::filter(Sample.ID %in% samples)
    ph <- plot_heatmap(norm_counts, proj.params, metadata_col)
    
    jpeg(file.path(save.dir, "Heatmap.jpeg"), width = 5, height = 7, units = "in", res = 300)
    gridExtra::grid.arrange(grobs = list(ph$ph$gtable), ncol = 1)
    dev.off()
    
    save_xlsx(ph$mat, file.path(save.dir, "Heatmap_Matrix.xlsx"), "Heatmap_matrix", row_names = TRUE)
    
    # ---- Pathway analysis ----
    gmt_files <- list.files(file.path(gmt.dir, proj.params$species), full.names = TRUE)
    gsea_res_list <- pathway_analysis(DEGs_df, gmt_files)
    
    wb <- openxlsx::createWorkbook()
    for (i in seq_along(gsea_res_list)) {
      openxlsx::addWorksheet(wb, sheetName = names(gsea_res_list)[i])
      openxlsx::writeData(wb, sheet = names(gsea_res_list)[i], x = gsea_res_list[[i]], rowNames = FALSE)
    }
    openxlsx::saveWorkbook(wb, file.path(save.dir, "Pathway_results.xlsx"), overwrite = TRUE)
    
    # ---- Pathway plots ----
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
    
    plot_pathways(top_pathways_gsea, "GSEA", save.dir)
    plot_pathways(top_pathways_ora, "ORA", save.dir)
    
    # ---- Pathway heatmaps ----
    top_pathways <- dplyr::bind_rows(top_pathways_gsea, top_pathways_ora) %>%
      dplyr::mutate(Pathway_name = stringr::str_wrap(paste0(Collection, ":", Description), width = 30))
    
    collections <- unique(top_pathways$Collection)
    
    for (col_name in collections) {
      top_pathways_collection <- top_pathways %>% dplyr::filter(Collection == col_name)
      plot_names <- unique(top_pathways_collection$Pathway_name)
      all_plots <- list()
      
      for (plot_name in plot_names) {
        plot_genes <- top_pathways_collection %>% 
          dplyr::filter(Pathway_name == plot_name) %>%
          dplyr::select(dplyr::starts_with("gene", ignore.case = FALSE)) %>%
          unlist(use.names = FALSE) %>%
          trimws() %>%
          na.omit() %>% 
          unique()
        
        norm_counts <- vst_counts[plot_genes, samples, drop = FALSE]
        metadata_col <- meta_data %>% dplyr::filter(Sample.ID %in% samples)
        proj.params$heatmap.title <- stringr::str_wrap(string = plot_name, width = 30)
        ph <- plot_heatmap(norm_counts, proj.params, metadata_col)
        all_plots[[length(all_plots) + 1]] <- trim_heatmap_whitespace(ph$ph$gtable)
      }
      
      jpeg(file.path(save.dir, paste0("Heatmap_", col_name, ".jpeg")),
           width = 7 * length(all_plots), 
           height = 10, units = "in", res = 300)
      gridExtra::grid.arrange(grobs = all_plots, ncol = length(all_plots))
      dev.off()
    }
    
    # ---- Perform TF analysis on DEGs using t-statistics ----
    t_stats_mat <- DEGs_df %>% 
      as.data.frame() %>%
      dplyr::mutate(t = -log10(padj) * log2FoldChange) %>%
      dplyr::filter(!is.na(t)) %>%
      dplyr::select(SYMBOL, t) %>%
      tibble::column_to_rownames("SYMBOL") %>%
      as.matrix()
    
    tf_res_degs <- tf_analysis(t_stats_mat, proj.params$species)
    
    # ---- Perform TF analysis on vst counts ----
    samples <- get_samples(meta_data, contrast[n])
    norm_counts_sub <- vst_counts[, samples]
    
    tf_res_counts <- tf_analysis(norm_counts_sub, proj.params$species)
    
    # ---- Plot barplots of top TFs from DEGs ----
    # NOTE: wsum returns wsum, norm_wsum and corr_wsum.
    # wsum (DONT USE): Biased toward larger gene sets (more genes → bigger sum)
    # norm_wsum (USE): Adjusts for pathway length so small and large gene sets are comparable
    # corr_sum (USE): corrects for high correlation as it can make enrichment appear stronger
    stats_degs <- unique(tf_res_degs$all_tfs$statistic)
    bar_plots <- list()
    
    # Extract treatment and control from current contrast string
    contrast_split <- stringr::str_split(string = contrast[n], pattern = "-")[[1]]
    treatment <- contrast_split[1]
    control <- contrast_split[2]
    
    for (stat in stats_degs) {
      n_tfs <- 20
      top_tf <- tf_res_degs$all_tfs %>%
        dplyr::mutate(Direction = dplyr::case_when(score < 0 ~ "Downregulated",
                                                   score > 0 ~ "Upregulated",
                                                   TRUE ~ "No change")) %>%
        dplyr::group_by(statistic, Direction) %>%
        dplyr::slice_max(order_by = abs(score), n = n_tfs, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::filter(statistic == stat)
      
      p1 <- ggplot(data = top_tf, aes(x = reorder(source, score), y = score, fill = score)) +
        geom_col(width = 0.75, na.rm = TRUE) +
        scale_fill_gradient2(low = "darkblue", high = "indianred", mid = "whitesmoke", midpoint = 0) +
        labs(x = "", y = "Score", title = paste0("Top TFs (", stat, " method)"), fill = "Score") +
        custom_theme +
        coord_cartesian(clip = "off") +
        geom_text(label = paste0("Activated in ", treatment),
                  x = top_tf$source[which.max(top_tf$score)],
                  y = ceiling(max(top_tf$score)) + 1, hjust = 1, color = "indianred", fontface = "bold") +
        geom_text(label = paste0("Activated in ", control),
                  x = top_tf$source[which.min(top_tf$score)],
                  y = ceiling(max(top_tf$score)) + 1, hjust = 0, color = "darkblue", fontface = "bold")
      
      bar_plots <- c(bar_plots, list(p1))
    }
    
    # Save combined barplots
    bar_plots_combined <- cowplot::plot_grid(plotlist = bar_plots, ncol = 1, nrow = 6)
    ggsave(filename = file.path(proj.dir, contrast[n], "Bar_plot_TFs.jpeg"),
           plot = bar_plots_combined,
           device = "jpeg",
           width = 10,
           height = 6 * 3,
           units = "in",
           dpi = 300,
           bg = "white")
    
    # ---- Plot heatmaps of top variable TFs across samples ----
    heatmap.params <- list(heatmap.force.log = FALSE,
                           heatmap.col.ann   = c("Treatment"),
                           heatmap.row.ann   = NULL,
                           heatmap.col.gaps  = NULL,  
                           heatmap.row.gaps  = NULL,  
                           heatmap.col.cluster  = c("Treatment"),
                           heatmap.row.cluster  = "all",
                           heatmap.palette      = "rdbu",
                           heatmap.ann.palette  = "discrete",
                           heatmap.border.color = NA,
                           heatmap.show.expr.legend = TRUE,
                           heatmap.title        = NA,
                           heatmap.format       = "tiff")
    
    stats_counts <- unique(tf_res_counts$all_tfs$statistic)
    heatmap_plots <- list()
    
    for (stat in stats_counts) {
      n_tfs <- 20
      top_tf <- tf_res_counts$all_tfs %>%
        dplyr::filter(statistic == stat) %>%
        dplyr::group_by(source) %>%
        dplyr::summarise(std = sd(score)) %>%
        dplyr::slice_max(order_by = abs(std), n = n_tfs, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::pull(source)
      
      plot_df <- tf_res_counts$all_tfs %>%
        dplyr::filter(statistic == stat) %>%
        tidyr::pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
        tibble::column_to_rownames("condition") %>%
        dplyr::select(all_of(top_tf)) %>%
        as.matrix() %>%
        t()
      
      samples <- get_samples(meta_data, contrast[n])
      metadata_col <- meta_data %>% dplyr::filter(Sample.ID %in% samples)
      heatmap.params$heatmap.title <- paste0("Top TFs (", stat, ") method")
      ph <- plot_heatmap(plot_df, heatmap.params, metadata_col)
      
      # Remove whitespace around heatmap
      ph$ph$gtable$widths <- grid::unit.pmax(ph$ph$gtable$widths, grid::unit(0.1, "cm"))
      ph$ph$gtable$heights <- grid::unit.pmax(ph$ph$gtable$heights, grid::unit(0.1, "cm"))
      
      heatmap_plots <- c(heatmap_plots, list(ph$ph$gtable))
    }
    
    # Save combined heatmaps ----
    n_plots <- length(heatmap_plots)
    jpeg(filename = file.path(proj.dir, contrast[n], "Heatmap_TFs.jpeg"),
         width = 7 * n_plots,
         height = 10,
         units = "in",
         res = 300)
    gridExtra::grid.arrange(grobs = heatmap_plots, ncol = n_plots, nrow = 1)
    dev.off()
  }
  
  # ---- Plot Dispersion Estimates for Entire Dataset ----
  plotDispEst_DESeq2(dds.list$dds, proj.dir)
  
  # ---- Generate and Annotate Blinded VST Counts ----
  # NOTE: vst counts are affected by design ONLY when blind=FALSE
  vsd_blind <- DESeq2::vst(dds.list$dds, blind = TRUE)
  vst_counts_blind <- SummarizedExperiment::assay(vsd_blind) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    add_annotation() %>%
    dplyr::select(SYMBOL, everything(), -dplyr::any_of(c("ENSEMBL_SYMBOL", "ENSEMBL_ID", "ENTREZ_SYMBOL", "ENTREZ_ID"))) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarize(across(.cols = where(is.numeric), .fns = mean)) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("SYMBOL") %>%
    as.matrix()
  
  # ---- Save Blinded VST Counts to Excel Using Helper ----
  save_xlsx(vst_counts_blind,
            file = file.path(proj.dir, "VST_counts_blind.xlsx"),
            sheet_name = "VST_blind",
            row_names = TRUE)
  
  # ---- Perform LRT Test for Entire Dataset ----
  dds_LRT <- DESeq2::DESeq(dds.list$dds, test = "LRT", reduced = ~1)
  res_LRT <- DESeq2::results(dds_LRT) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    add_annotation()
  
  # ---- Identify Significant Genes for Heatmap ----
  sig_genes <- res_LRT %>%
    dplyr::filter(padj <= 0.05) %>%
    dplyr::pull(SYMBOL)
  
  norm_counts <- vst_counts_blind[rownames(vst_counts_blind) %in% sig_genes, ]
  metadata_col <- meta_data
  proj.params$heatmap.title <- ""
  
  # ---- Plot and Save Heatmap Using Helper ----
  ph <- plot_heatmap(norm_counts, proj.params, metadata_col)
  
  jpeg(filename = file.path(proj.dir, "Heatmap.jpeg"),
       width = 14, height = 10, units = "in", res = 300)
  gridExtra::grid.arrange(grobs = list(ph$ph$gtable), ncol = 1, nrow = 1)
  dev.off()
}

# ---- Custom Theme ----
custom_theme <- ggplot2::theme(
  plot.title   = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0.5),
  legend.title = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0,   vjust = 1, angle = 0),
  axis.title.x = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0.5, vjust = 0, angle = 0),
  axis.title.y = element_text(family = "sans", face = "bold",  colour = "black", size = 12, hjust = 0.5, vjust = 1, angle = 90),
  legend.text  = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5),
  axis.text.x  = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 0.5, vjust = 0.5, angle = 45),
  axis.text.y  = element_text(family = "sans", face = "plain", colour = "black", size = 10, hjust = 1,   vjust = 0.5, angle = 0))

# ---- Custom Color Palette ----
scanpy_default_102 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                        "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
                        "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", 
                        "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173", 
                        "#3182bd", "#e6550d", "#31a354", "#756bb1", "#636363", 
                        "#9e9ac8", "#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4",
                        "#fed9a6", "#ffffcc", "#e5d8bd", "#fddaec", "#f2f2f2",
                        "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", 
                        "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", 
                        "#ffff99", "#b15928", "#8dd3c7", "#ffffb3", "#bebada",
                        "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
                        "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "#a65628",
                        "#ffff99", "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
                        "#66a61e", "#e6ab02", "#a6761d", "#666666", "#e41a1c", 
                        "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", 
                        "#a65628", "#f781bf", "#999999", "#66c2a5", "#fc8d62",
                        "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494",
                        "#b3b3b3", "#1b9e77", "#d95f02", "#7570b3", "#e7298a", 
                        "#66a61e", "#e6ab02", "#a6761d", "#666666")

# ---- Key Functions ----
merge_counts <- function(proj.dir, proj.params) {
  
  # Get count files
  count_files <- list.files(path = proj.params$counts.dir, pattern = "\\.txt$|ReadsPerGene\\.out\\.tab$", full.names = TRUE)
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
  filename <- file.path(proj.dir, paste0(proj.params$proj, ".raw.counts.xlsx"))
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb = wb, sheetName = "Raw_counts")
  openxlsx::writeData(wb = wb, sheet = "Raw_counts", x = count_matrix)
  openxlsx::saveWorkbook(wb = wb, file = filename, overwrite = TRUE)
  message("✅ Saved counts to: ", file.path(proj.dir, proj.params$proj))  
  
  return(count_matrix)
}

prep_metadata <- function(meta_data, read_data) { 
  
  # Input checks 
  if (!"Sample.ID" %in% colnames(meta_data)) {
    stop("`meta_data` must contain a 'Sample.ID' column.")
  }
  
  # Filter, clean and align with read_data 
  meta_data <- meta_data %>%
    dplyr::filter(!is.na(Sample.ID)) %>%
    dplyr::mutate(Sample.ID = make.names(Sample.ID, unique = TRUE)) %>%
    dplyr::filter(Sample.ID %in% make.names(colnames(read_data))) %>%
    dplyr::mutate(row_names = Sample.ID) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "row_names")
  
  # Add Batch column if missing 
  if (!"Batch" %in% colnames(meta_data)) {
    meta_data$Batch <- 1
    warning("No 'Batch' column found. Assigning all samples to Batch 1.")
  }
  
  # Return cleaned metadata 
  return(invisible(meta_data))
}

prep_readdata <- function(meta_data, read_data) {
  
  # Input checks 
  if (!"SYMBOL" %in% colnames(read_data)) {
    stop("Input read_data must contain a 'SYMBOL' column.")
  }
  if (anyDuplicated(read_data$SYMBOL)) {
    stop("The 'SYMBOL' column in read_data must not contain duplicate values.")
  }
  
  # Filter and clean data 
  colnames(read_data) <- make.names(colnames(read_data))
  valid_samples <- intersect(colnames(read_data), make.names(rownames(meta_data)))
  
  read_data <- read_data %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "SYMBOL") %>%
    dplyr::select(all_of(valid_samples)) %>%
    base::replace(is.na(.), 0)
  
  # Return cleaned read data 
  return(invisible(read_data))
  
}

check_data <- function(meta_data, read_data, proj.params) {
  
  # Remove samples with zero counts 
  zero_cols <- which(colSums(read_data) == 0)
  if (length(zero_cols) > 0) {
    read_data <- read_data[, -zero_cols, drop = FALSE]
    meta_data <- meta_data[-zero_cols, , drop = FALSE]
  }
  
  # Remove samples with NA in any DEG variable 
  NA_samples <- integer()
  DEG_variables <- stringr::str_split(string = proj.params$deseq2.design, pattern = "[+*:]")[[1]]
  
  for (var in DEG_variables) {
    if (var %in% colnames(meta_data)) {
      NA_samples <- c(NA_samples, which(is.na(meta_data[[var]])))
    } else {
      warning(sprintf("Variable '%s' not found in metadata.", var))
    }
  }
  NA_samples <- unique(NA_samples)
  if (length(NA_samples) > 0) {
    read_data <- read_data[, -NA_samples, drop = FALSE]
    meta_data <- meta_data[-NA_samples, , drop = FALSE]
  }
  
  # Remove sizeFactor column if present 
  meta_data <- meta_data[, setdiff(colnames(meta_data), "sizeFactor"), drop = FALSE]
  
  # Convert all metadata columns to factors 
  message("Structure of meta_data before conversion:")
  str(meta_data)
  for (i in base::seq_len(ncol(meta_data))) {
    meta_data[[i]] <- as.factor(meta_data[[i]])
  }
  message("Structure of meta_data after conversion:")
  str(meta_data)
  
  # Rearrange read_data columns to match meta_data 
  read_data <- read_data[, rownames(meta_data), drop = FALSE]
  
  # Sanity checks 
  if (!is.data.frame(read_data)) {
    stop("`read_data` must be a data.frame")
  }
  if (!is.data.frame(meta_data)) {
    stop("`meta_data` must be a data.frame")
  }
  if (!all(colnames(read_data) %in% rownames(meta_data))) {
    stop("Some samples in `read_data` are missing in `meta_data`")
  }
  if (!all(colnames(read_data) == rownames(meta_data))) {
    stop("Sample order mismatch between `read_data` and `meta_data`")
  }
  
  # Return DESeq2 ready metadata and read data 
  return(invisible(list(meta_data = meta_data, read_data = read_data)))
}

plotPCA_DESeq2 <- function(meta_data, read_data, proj.dir) {
  
  # Input Checks 
  stopifnot(is.data.frame(meta_data))
  stopifnot(is.data.frame(read_data))
  if (!"Sample.ID" %in% colnames(meta_data)) {
    stop("meta_data must contain a 'Sample.ID' column.")
  }
  
  # DESeq2 Object Preparation 
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = read_data,
                                        colData = meta_data,
                                        design = ~1)
  dds <- DESeq2::estimateSizeFactors(dds)
  
  # VST Transformation 
  vsd <- DESeq2::vst(dds, blind = TRUE)
  
  # PCA Data Preparation 
  vst_mat <- SummarizedExperiment::assay(vsd) %>%
    as.data.frame() %>%
    dplyr::mutate(row_variance = matrixStats::rowVars(as.matrix(.))) %>%
    dplyr::slice_max(order_by = row_variance, n = 500) %>%
    dplyr::select(-row_variance)
  
  # PCA Calculation 
  pca <- stats::prcomp(x = t(vst_mat), center = TRUE, scale. = FALSE)
  
  # Merge PCA Output with Metadata 
  pca_df <- pca$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample.ID")
  
  df <- dplyr::inner_join(meta_data, pca_df, by=c("Sample.ID"="Sample.ID"))
  
  # PCA Plot 
  percentVar <- round(100 * summary(pca)$importance[2, 1:2])
  
  comp_variables <- setdiff(colnames(meta_data), "Sample.ID")
  
  for (var in comp_variables) {
    
    # Define color palette
    pca_palette <- scanpy_default_102[1 : length(unique(meta_data[[var]]))]
    names(pca_palette) <- as.character(unique(meta_data[[var]]))
    
    p <- ggplot2::ggplot(data = df, mapping = aes(x = PC1, y = PC2, color = get(var))) +
      ggplot2::geom_point(size = 3, shape = 16) +
      ggrepel::geom_text_repel(ggplot2::aes(label = Sample.ID), show.legend = FALSE) +
      ggplot2::theme_light() +
      ggplot2::labs(color = var,
                    x = paste0("PC1: ", percentVar[1], "% variance"),
                    y = paste0("PC2: ", percentVar[2], "% variance")) +
      custom_theme + 
      scale_color_manual(values = pca_palette)
    
    # Save Plot 
    ggplot2::ggsave(filename = paste0("PCA_Plot_", var, ".tiff"),
                    plot = p,
                    device = "tiff",
                    path = file.path(proj.dir),
                    width = 7,
                    height = 7,
                    units = "in",
                    dpi = 300,
                    limitsize = TRUE,
                    bg = "white")
  }
}

get_annotations <- function() {
  
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

run_deseq2 <- function(dds, contrast, proj.params) {
  
  # Pre-filter to remove lowly expressed genes to improve sizefactor estimation
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # Run DESeq2 with Parametric and Local Fit 
  dds_para <- DESeq2::DESeq(object = dds, 
                            test = "Wald", 
                            fitType = "parametric", 
                            betaPrior = FALSE, 
                            minReplicatesForReplace = 7)
  
  dds_local <- DESeq2::DESeq(object = dds,
                             test = "Wald",
                             fitType = "local",
                             betaPrior = FALSE,
                             minReplicatesForReplace = 7)
  
  # Choose Best Fit Type Based on Residuals 
  residual_para <- mcols(dds_para)$dispGeneEst - mcols(dds_para)$dispFit
  residual_local <- mcols(dds_local)$dispGeneEst - mcols(dds_local)$dispFit
  dds <- if (median(residual_para^2, na.rm = TRUE) <= median(residual_local^2, na.rm = TRUE)) {
    dds_para
  } else {
    dds_local
  }
  
  # Prepare Contrast Vector 
  mod_mat <- model.matrix(design(dds), colData(dds))
  design.factors <- stringr::str_split(string = proj.params$deseq2.design, 
                                       pattern = "[+*:]")[[1]] %>% unique()
  
  df <- colData(dds) %>% 
    as.data.frame() %>% 
    tidyr::unite(col = "Groups", design.factors, sep = ".")
  
  groups <- unique(df$Groups)
  for (i in groups) {
    val <- colMeans(mod_mat[df$Groups == i, ])
    assign(x = i, value = val)
  }
  
  contrast_vec <- base::eval(expr = base::parse(text = contrast))
  
  # Run DESeq2 Results with LFC Threshold and Shrinkage 
  res <- DESeq2::results(object = dds,
                         contrast = contrast_vec,
                         lfcThreshold = proj.params$deseq2.lfc.cutoff,
                         altHypothesis = "greaterAbs",
                         cooksCutoff = TRUE,
                         independentFiltering = TRUE,
                         alpha = proj.params$deseq2.padj.cutoff,
                         pAdjustMethod = "BH")
  
  set.seed(1234)
  res <- DESeq2::lfcShrink(dds = dds, res = res, type = "ashr")
  summary(res)
  
  # Return dds and results 
  return(invisible(list(dds = dds, res = res)))
}

plotMA_DESeq2 <- function(dds, contrast = "", proj.dir) {
  
  # Input Check 
  if (!inherits(dds, "DESeqDataSet")) {
    stop("dds must be a DESeqDataSet object.")
  }
  save.dir <- file.path(proj.dir, contrast)
  if (!dir.exists(save.dir)) {
    dir.create(save.dir, recursive = TRUE)
  }
  
  # Plot MA 
  grDevices::tiff(filename = file.path(save.dir, "MA_Plot.tiff"),
                  width = 8.5,
                  height = 11,
                  units = "in",
                  bg = "white",
                  res = 600)
  
  DESeq2::plotMA( object = dds,
                  alpha = 0.1, # FDR threshold (blue dots for significant genes)
                  main = "MA Plot",
                  xlab = "Mean of Normalized Counts",
                  MLE = FALSE)
  
  grDevices::dev.off()
}

plotDispEst_DESeq2 <- function(dds, proj.dir) {
  
  # Input Check 
  if (!inherits(dds, "DESeqDataSet")) {
    stop("`dds` must be a DESeqDataSet object.")
  }
  
  # Create JPEG Output 
  grDevices::jpeg(filename = file.path(proj.dir, "Dispersion_Plot.tiff"), 
                  width = 8.5,
                  height = 11,
                  units = "in",
                  quality = 75,
                  bg = "white",
                  res = 300)
  
  # Plot Dispersion Estimates 
  # Expected results: Higher the mean, lower the dispersion
  DESeq2::plotDispEsts(object = dds,
                       genecol = "black",
                       fitcol = "red",
                       finalcol = "dodgerblue",
                       legend = TRUE,
                       xlab = "Mean of Normalized Counts",
                       ylab = "Dispersion",
                       log = "xy",
                       cex = 0.45)
  
  # Close Device 
  grDevices::dev.off()
}

# DEGs_df with column SYMBOL, padj, log2FoldChange
# k <- # overlapping genes between input and pathway
# n <- # overlapping genes between input and collection
# K <- # genes in pathway
# N <- # genes in collection
pathway_analysis <- function(DEGs_df, gmt_files) {
  
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
                                      verbose = TRUE,
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
  ora_df <- dplyr::bind_rows(ora_df_up %>% dplyr::mutate(Direction = "Upregulated"), 
                             ora_df_down %>% dplyr::mutate(Direction = "Downregulated")) %>%
    tidyr::separate(col = GeneRatio, into = c("k", "n")) %>%
    tidyr::separate(col = BgRatio, into = c("K", "N")) %>%
    dplyr::mutate_at(c("k", "n", "K", "N"), as.numeric) %>%
    dplyr::mutate(GeneRatio = k / n,
                  BackgroundRatio = K / N,
                  EnrichmentRatio = GeneRatio / BackgroundRatio,
                  combined_score = GeneRatio * -log10(p.adjust),
                  NES = NA_integer_)
  
  for (i in c("fgsea_df", "gsea_df")){
    df <- get(i) %>% 
      dplyr::mutate(Direction = dplyr::case_when(NES > 0 ~ "Upregulated",
                                                 NES < 0 ~ "Downregulated",
                                                 TRUE ~ "No change"))
    assign(x = i, value = df)
  }
  
  # Rename columns consistently across different methods
  lookup <- c(pathway = "ID",
              geneID = "leadingEdge", geneID = "core_enrichment",
              K = "size", K = "setSize",
              padj = "p.adjust", 
              pval = "pvalue")
  
  for (i in c("fgsea_df", "gsea_df", "ora_df")){
    
    df <- get(i) %>%
      dplyr::rename(any_of(lookup)) %>%
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
    
    max_len <- max(df$leading_edge_size)
    sep_char <- ifelse(grepl(",", df$geneID[1], fixed = TRUE), ",", "/")
    
    # NOTE: across() allows you to apply a function(s) to multiple columns at 
    # once. map_chr() iterates over a list/vector and applies a function to 
    # each element (in this context, each cell within the selected columns),
    # returning a character vector.
    df <- df %>% 
      tidyr::separate(col = geneID, into = paste0("gene", 1:max_len), sep = sep_char, remove = TRUE, fill = "right") %>%
      dplyr::mutate(across(.cols = starts_with("gene", ignore.case = FALSE), 
                           .fns = function (col) { col %>% 
                               purrr::map_chr(.f = function (x){gsub('c\\(', "", x) }) %>%
                               purrr::map_chr(.f = function (x){gsub('\\)', "", x) }) %>%
                               purrr::map_chr(.f = function (x){gsub('"', "", x) }) %>%
                               trimws() })) %>%
      dplyr::select(Collection, Description, leading_edge_size, K, padj, NES, Direction, everything())
    
    assign(x = i, value = df)
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
  
  # Return results
  return(invisible(list(fgsea = fgsea_df, 
                        gsea = gsea_df, 
                        ora = ora_df,
                        consensus = consensus_df)))
}

add_major_pathway <- function(df){
  
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

### decoupleR Analysis [RECOMMENDED as it can run multiple algorithms & give consensus]
# input can be vst_counts or DEGs with t-statistic (-log10padj*log2FC)
tf_analysis <- function(input, species = "Homo sapiens", top = 500) {
  
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
  pathway_net <- decoupleR::get_progeny(organism = organism, top = top)
  tf_net <- decoupleR::get_collectri(organism = organism, split_complexes = FALSE)
  
  # Run decoupleR consensus 
  # stats <- c("aucell", fgsea", "gsva", "mdt", "mlm", "ora", "udt", "ulm", "viper", "wmean", "wsum")
  stats <- c("mlm", "ulm", "wsum")
  pathway_df <- decoupleR::decouple(mat = input, network = pathway_net,
                                    statistics = stats, minsize = 5)
  
  tf_df <- decoupleR::decouple(mat = input, network = tf_net,
                               statistics = stats, minsize = 5)
  
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

plot_volcano <- function(DEGs_df, contrast = "", proj.dir){
  
  # Input checks 
  required_cols <- c("log2FoldChange", "padj", "SYMBOL")
  if (!all(required_cols %in% colnames(DEGs_df))) {
    stop("DEGs_df must contain columns: log2FoldChange, padj, SYMBOL")
  }
  save.dir <- file.path(proj.dir, contrast)
  if (!dir.exists(save.dir)) {
    dir.create(save.dir, recursive = TRUE)
  }
  
  padj.cutoff <- proj.params$volcano.padj.cutoff
  lfc.cutoff <- proj.params$volcano.lfc.cutoff
  target <- stringr::str_split(string = contrast, pattern = "-")[[1]][1]
  reference <- stringr::str_split(string = contrast, pattern = "-")[[1]][2]
  
  # Format data 
  DEGs_df <- DEGs_df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(Direction = dplyr::case_when(padj < padj.cutoff & log2FoldChange > lfc.cutoff ~ paste0("Up in ", target),
                                               padj < padj.cutoff & log2FoldChange < -lfc.cutoff ~ paste0("Up in ", reference),
                                               TRUE ~ "Not Significant"),
                  padj = dplyr::case_when(padj == 0 ~ sort(unique(padj))[2], 
                                          TRUE ~ padj),
                  Significance = dplyr::case_when(abs(log2FoldChange) >= lfc.cutoff & padj <= 0.001 ~ "FDR < 0.001",
                                                  abs(log2FoldChange) >= lfc.cutoff & padj <= 0.01 ~ "FDR < 0.01",
                                                  abs(log2FoldChange) >= lfc.cutoff & padj <= 0.05 ~ "FDR < 0.05",
                                                  TRUE ~ "Not Significant"),
                  Relevance = abs(log2FoldChange) * -log10(padj))
  
  # Axis limits and breaks
  x_vals <- DEGs_df$log2FoldChange
  y_vals <- -log10(DEGs_df$padj)
  x_vals <- x_vals[is.finite(x_vals)]
  y_vals <- y_vals[is.finite(y_vals)]
  x_lims <- stats::quantile(x = x_vals, probs = c(0, 1), na.rm = TRUE) %>% unname()
  y_lims <- stats::quantile(x = y_vals, probs = c(0, 1), na.rm = TRUE) %>% unname()
  x_min <- floor(x_lims[1])
  x_max <- ceiling(x_lims[2])
  y_min <- floor(y_lims[1])
  y_max <- ceiling(y_lims[2])
  
  bin <- max(abs(floor(x_min / 5)), abs(ceiling(x_max / 5)))
  x_breaks <- seq(from = -max(x_max, abs(x_min)), to = max(x_max, abs(x_min)), by = bin)
  x_breaks <- x_breaks[!x_breaks <= x_min-bin]
  x_breaks <- x_breaks[!x_breaks >= x_max+bin] 
  
  y_breaks <- seq(from = y_min, to = ceiling(y_max/10)*10, by = dplyr::case_when(y_max%/%100 > 0 ~ 100,
                                                                                 y_max%/%10 > 0 ~ 10,
                                                                                 TRUE ~ 1))
  
  # Color palettes 
  volcano_palette <- c(viridis_pal()(100)[100], viridis_pal()(100)[50], viridis_pal()(100)[1])
  names(volcano_palette) <- c("Not Significant",
                              paste0("Up in ", reference),
                              paste0("Up in ", target))
  
  alpha_palette <- c("FDR < 0.001" = 1, 
                     "FDR < 0.01" = 0.8,
                     "FDR < 0.05" = 0.6,
                     "Not Significant" = 0.4)
  
  # Build ggplot object 
  p <- ggplot2::ggplot(data = DEGs_df, 
                       mapping = aes(x = log2FoldChange, 
                                     y = -log10(padj), 
                                     color = Direction, 
                                     alpha = Significance,
                                     size = Relevance)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.05, height = 0.05)) +
    ggplot2::theme_classic() +
    custom_theme +
    ggplot2::labs(x = expression("log"[2]*"FC"),
                  y = expression("-log"[10]*"padj"),
                  fill = "Direction",
                  title = contrast) +
    ggplot2::geom_vline(xintercept = c(-lfc.cutoff, lfc.cutoff),
                        color = "black", linetype = "dotted", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(padj.cutoff),
                        color = "black", linetype = "dotted", linewidth = 0.5) +
    ggplot2::coord_cartesian(xlim = c(min(x_breaks), max(x_breaks)),
                             ylim = c(min(y_breaks), max(y_breaks))) +
    ggplot2::scale_x_continuous(breaks = x_breaks, 
                                labels = function(x) { base::ifelse(x %% 1 == 0, as.integer(x), format(x, digits = 2)) }) +
    ggplot2::scale_y_continuous(breaks = y_breaks) +
    ggplot2::scale_color_manual(values = volcano_palette) + 
    ggplot2::scale_alpha_manual(values = alpha_palette) +
    ggplot2::scale_size_continuous(range = c(0, 3)) +
    ggplot2::guides(size = "none",
                    shape = guide_legend(override.aes = list(size = 3)),
                    fill = guide_colourbar(theme = theme(legend.key.width = unit(0.75, "lines"),
                                                         legend.key.height = unit(10, "lines"),
                                                         legend.ticks = element_blank(),
                                                         legend.frame = element_rect(colour = "Black", linewidth = 1)))) 
  
  # Identify top and bottom 5 genes 
  predicted_gene_pattern <- "^(Gm[0-9]+|ENSMUSG[0-9]+|ENSG[0-9]+|LOC[0-9]+|C[0-9]+orf[0-9]+|RP[0-9]+-)|Rik$"
  proper_DEGs_df <- DEGs_df %>%
    dplyr::filter(!(stringr::str_detect(string = SYMBOL, pattern = predicted_gene_pattern))) %>%
    dplyr::filter(padj < padj.cutoff)
  top_genes <- c(proper_DEGs_df %>%
                   dplyr::filter(log2FoldChange < -lfc.cutoff) %>%
                   #dplyr::slice_min(padj, n = 5, with_ties = FALSE) %>%
                   dplyr::slice_max(Relevance, n = 5, with_ties = FALSE) %>%
                   dplyr::pull(SYMBOL),
                 DEGs_df %>%
                   dplyr::filter(log2FoldChange > lfc.cutoff) %>%
                   #dplyr::slice_min(padj, n = 5, with_ties = FALSE) %>%
                   dplyr::slice_max(Relevance, n = 5, with_ties = FALSE) %>%
                   dplyr::pull(SYMBOL))
  
  # Add labels for custom or top genes 
  if (!is.null(proj.params$volcano.label.genes)) {
    disp_genes <- proj.params$volcano.label.genes
  } else {
    disp_genes <- top_genes
  }
  
  q <- p + ggrepel::geom_text_repel(data = DEGs_df %>% 
                                      dplyr::filter(SYMBOL %in% base::intersect(disp_genes, DEGs_df$SYMBOL)),
                                    aes(label = SYMBOL),
                                    direction = "both",
                                    box.padding = 0.8,         # ↓ smaller padding around label
                                    point.padding = 0.1,       # minimal space between point and line start
                                    max.overlaps = nrow(DEGs_df),
                                    show.legend = FALSE,
                                    min.segment.length = 0,    # Only draw segments longer than this
                                    segment.curvature = -0.5,  # Negative = curve upward, positive = downward
                                    segment.ncp = 50,          # More control points = smoother curves
                                    segment.angle = 20,        # Affects entry/exit angles
                                    segment.size = 0.5,        # Optional: line thickness
                                    size = 4,                  # text size in mm (1 mm = 2.83 points)
                                    position = ggbeeswarm::position_quasirandom(width = 0.1, varwidth = TRUE))
  
  # Save plot 
  ggplot2::ggsave(filename = file.path(proj.dir, contrast, "Volcano_Plot.tiff"),
                  plot = p,
                  device = "tiff",
                  width = 7,
                  height = 7,
                  units = "in",
                  dpi = 300,
                  limitsize = TRUE,
                  bg = "white")
  
  ggplot2::ggsave(filename = file.path(proj.dir, contrast, "Volcano_Plot_top.tiff"),
                  plot = q,
                  device = "tiff",
                  width = 7,
                  height = 7,
                  units = "in",
                  dpi = 300,
                  limitsize = TRUE,
                  bg = "white")
  return(invisible(p))
}

# metadata has column Sample.ID and columns defined in proj.params$heatmap.col.ann
# metadata_row has column SYMBOL and columns defined in proj.params$heatmap.row.ann
# metadata_row <- data.frame(SYMBOL = c("")) if no row annotations needed
# disp_genes either empty vector c() or vector of genes
plot_heatmap <- function(norm_counts, proj.params, metadata_col = NULL, metadata_row = NULL, disp_genes = c()) {
  
  # Prepare Matrix: Filter, Normalize, Transform 
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
  
  rownames(mat) <- make.names(rownames(mat), unique = TRUE)
  colnames(mat) <- make.names(colnames(mat), unique = TRUE)
  
  # log transform if needed 
  quantiles <- stats::quantile(x = as.vector(as.matrix(mat)), probs = c(0., 0.01, 0.99, 1.0), na.rm=T)
  huge_range <- quantiles[4] - quantiles[1] > 100     # Range of values greater than 100
  only_pos_values <- quantiles[1] >= 0          # Min value greater than 0
  if((huge_range & only_pos_values) | proj.params$heatmap.force.log){
    mat <- log2(1 + mat)
  }
  
  # Scale every feature across samples 
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[is.na(mat_scaled)] <- 0
  
  # Define column annotations 
  if(is.null(metadata_col)) {
    col_annotation <- NULL
  } else{
    col_annotation <- metadata_col %>%
      dplyr::select(Sample.ID, all_of(proj.params$heatmap.col.ann)) %>%
      dplyr::mutate(Sample.ID = make.names(Sample.ID)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("Sample.ID") %>%
      as.data.frame() %>%
      mutate(across(where(is.factor), as.character))
  }
  
  # Define row annotations
  if(is.null(metadata_row)) {
    row_annotation <- NULL
  } else{
    row_annotation <- metadata_row %>%
      dplyr::select(SYMBOL, all_of(proj.params$heatmap.row.ann)) %>%
      dplyr::mutate(SYMBOL = make.names(SYMBOL, unique = TRUE)) %>%
      tibble::remove_rownames() %>%
      tibble::column_to_rownames("SYMBOL") %>%
      as.data.frame() %>%
      mutate(across(where(is.factor), as.character))
  }
  
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
  base_colors <- scanpy_default_102
  
  col_list <- base::lapply(X = as.list(col_annotation), FUN = function(x) { as.character(x) %>% unique})
  row_list <- base::lapply(X = as.list(row_annotation), FUN = function(x) { as.character(x) %>% unique})
  ann_list <- c(col_list, row_list)
  
  color_index <- 1
  for (i in seq_along(ann_list)) {  # Iterate through each annotation variable (Eg: CellType) 
    levels <- sort(ann_list[[i]])  # Get levels within each annotation variable (Eg: CT1, CT2)
    n_levels <- length(levels)    # Get number of levels within each annotation variable
    
    palette_colors <- if (proj.params$heatmap.ann.palette == "discrete" | n_levels == 1){
      base_colors[color_index:(color_index + n_levels - 1)]
    } else{
      alphas <- seq(1 / n_levels, 1, length.out = n_levels)
      base::sapply(X = alphas, 
                   FUN = function(x) { colorspace::adjust_transparency(col = base_colors[color_index], alpha = x) })
    }
    
    names(palette_colors) <- levels          # Name each color with levels
    ann_colors <- c(ann_colors, list(palette_colors)) # Append named color palette
    names(ann_colors)[i] <- names(ann_list)[i]     # Name the color palette with corresponding annotation variable name
    color_index <- color_index + n_levels       # Move to next color
  }
  
  # Define Color Palette for Heatmap 
  valid_palettes <- c("vrds", "rdbu")
  
  if (!proj.params$heatmap.palette %in% valid_palettes) {
    stop("Invalid heatmap palette. proj.params$heatmap.palette must be either 'vrds' or 'rdbu'.")
  }
  
  heatmap_palette <- switch(proj.params$heatmap.palette,
                            vrds = viridis::viridis(100),
                            rdbu = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100))
  
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
  
  # Define gaps in heatmap 
  gaps_col <- if (!gtools::invalid(proj.params$heatmap.col.gaps) & proj.params$heatmap.col.cluster %in% colnames(col_annotation)) {
    if (all(proj.params$heatmap.col.gaps %in% colnames(col_annotation))) {
      col_annotation %>%
        dplyr::count(get(proj.params$heatmap.col.gaps)) %>%
        dplyr::mutate(n = cumsum(n)) %>%
        dplyr::pull(n) %>%
        .[. < ncol(mat_scaled)]
    } else{
      message("Gaps not introduced as gap identifiers missing in column annotation")
    }
  }
  
  gaps_row <- if (!gtools::invalid(proj.params$heatmap.row.gaps) & proj.params$heatmap.row.cluster %in% colnames(row_annotation)) {
    if (all(proj.params$heatmap.row.gaps %in% colnames(row_annotation))) {
      row_annotation %>%
        dplyr::count(get(proj.params$heatmap.row.gaps)) %>%
        dplyr::mutate(n = cumsum(n)) %>%
        dplyr::pull(n) %>%
        .[. < nrow(mat_scaled)]
    } else{
      message("Gaps not introduced as gap identifiers missing in row annotation")
    }
  }
  
  # Determine Ordering 
  if (proj.params$heatmap.col.cluster == "all"){
    colclust <- hclust(dist(t(mat_scaled)))
    col_order <- colnames(mat_scaled)[colclust$order]
  }
  if (proj.params$heatmap.col.cluster == "alphabetical"){
    col_order <- sort(colnames(mat_scaled))
  }
  if (proj.params$heatmap.col.cluster %in% colnames(col_annotation)){
    
    # NOTE: While calculating gaps_col, we use count(). It sorts alphabetically.
    # So, WE MUST sort col_elements to match gaps_col
    col_order <- c()
    col_elements <- col_annotation %>% 
      dplyr::pull(proj.params$heatmap.col.cluster) %>%
      unique() %>% sort()
    
    for (g in col_elements){
      
      samples <- rownames(col_annotation)[col_annotation %>% dplyr::pull(proj.params$heatmap.col.cluster) == g]
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
  
  if (proj.params$heatmap.row.cluster == "all"){
    rowclust <- hclust(dist(mat_scaled))
    row_order <- rownames(mat_scaled)[rowclust$order]
  }
  if (proj.params$heatmap.row.cluster == "alphabetical"){
    row_order <- sort(rownames(mat_scaled))
  }
  if (proj.params$heatmap.row.cluster %in% colnames(row_annotation)){
    
    # NOTE: While calculating gaps_row, we use count(). It sorts alphabetically.
    # So, WE MUST sort row_elements to match gaps_row
    row_order <- c()
    row_elements <- row_annotation %>% 
      dplyr::pull(proj.params$heatmap.row.cluster) %>%
      unique() %>% sort()
    
    for (g in row_elements){
      
      genes <- rownames(row_annotation)[row_annotation %>% dplyr::pull(proj.params$heatmap.row.cluster) == g]
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
  
  # Prepare Matrix for Plotting 
  reordered <- mat_scaled[row_order, col_order]
  
  # Final formatting
  # Set font sizes
  fontsize <- 10
  fontsize.row <- fontsize*1
  fontsize.col <- fontsize*1
  fontsize.number <- fontsize*0.8
  
  # Set cell width, height (in points) dynamically
  cell.width <- dplyr::if_else(ncol(reordered) <= 60, fontsize+2, NA)
  cell.height <- dplyr::if_else(nrow(reordered) <= 48, fontsize+2, NA)
  
  # Truncate long row and column labels
  main.title <- stringr::str_wrap(string = proj.params$heatmap.title, width = 20)
  labels.col <- stringr::str_trunc(string = colnames(reordered), width = 15, side = "right", ellipsis = "…")
  labels.row <- stringr::str_trunc(string = rownames(reordered), width = 15, side = "right", ellipsis = "…")
  if(length(disp_genes) > 0){
    labels.row <- dplyr::if_else(rownames(reordered) %in% make.names(disp_genes), rownames(reordered), " ")
  }
  
  # Set column label angel
  angle.col <- 45
  
  # Plot Heatmap 
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
                           show_rownames     = !is.na(cell.height),
                           show_colnames     = !is.na(cell.width),
                           labels_row        = labels.row,
                           labels_col        = labels.col,
                           angle_col         = angle.col,
                           fontsize          = fontsize,         # points; 72 points = 1 inch
                           fontsize_row      = fontsize.row,     # points
                           fontsize_col      = fontsize.col,     # points
                           fontsize_number   = fontsize.number,  # points
                           silent            = TRUE, 
                           
                           main              = main.title,
                           border_color      = proj.params$heatmap.border.color,
                           legend            = proj.params$heatmap.show.expr.legend,
                           
                           scale                    = "none",
                           cluster_rows             = FALSE,
                           cluster_cols             = FALSE,
                           clustering_distance_rows = "euclidean",
                           clustering_distance_cols = "euclidean",
                           clustering_method        = "complete",
                           annotation_legend        = TRUE,
                           annotation_names_row     = FALSE,
                           annotation_names_col     = FALSE,
                           width                    = NA,               # inches
                           height                   = NA,               # inches
                           filename                 = NA)
  
  # Save Matrix to Excel 
  write_data <- if (ncol(reordered) > nrow(reordered)) t(reordered) else reordered
   
  return(invisible(list(ph = ph, mat = write_data)))
}

# No "_" in Description column so that str_wrap works
# Collection column needed
# (combined_score, GeneRatio, k) OR (NES, leading_edge_size) columns needed
plot_pathways <- function(df, method, save.dir){
  
  dot_plot_list <- list()
  bar_plots <- list()
  plot_colors <- c("Upregulated" = "#E69F00", "Downregulated" = "#56B4E9")
  
  # str_wrap() wraps only between words, and it defines "words" based on spaces (" ").
  # If all words are connected by "_", it wont split.
  df <- df %>%
    dplyr::mutate(Description = gsub(pattern = "_", replacement = " ", x = Description),
                  Description = stringr::str_wrap(string = Description, width = 30))
  
  # Decide if there are multiple collections or single collection.
  # If multiple, plot each collection separately,
  collections <- unique(df$Collection)
  n_collections <- length(unique(df$Collection))
  
  for (i in seq_len(n_collections)){
    
    plot_df <- df %>% dplyr::filter(Collection %in% collections[i])
    
    # Choose x axis column, size column and labels
    if (method == "ORA"){
      x_col <- sym("GeneRatio")
      x_label <- "GeneRatio"
      size_col <- sym("leading_edge_size")
      color_col <- sym("Direction")
      score_var <- dplyr::case_when("GeneRatio" %in% colnames(df) ~ "GeneRatio",
                                    "combined_score" %in% colnames(df) ~ "combined_score", 
                                    TRUE ~ NA_character_)
      x_limits <- c(0, NA)  # start at 0, auto end
      
      df <- df %>%
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
    
    # Plot bar plot 
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
    
    bar_plots <- c(bar_plots, list(p1))
    
    # Plot dot plot 
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
  }
  
  bar_plots <- cowplot::plot_grid(plotlist = bar_plots, ncol = 3, nrow = 3)
  ggplot2::ggsave(filename = paste0("Bar_plot_pathways_", method, ".tiff"),
                  plot = bar_plots,
                  device = "jpeg",
                  path = file.path(save.dir),
                  width = 3*7,
                  height = 3*7,
                  units = "in",
                  dpi = 300,
                  bg = "white")
  
  dot_plots <- cowplot::plot_grid(plotlist = dot_plot_list, ncol = 3, nrow = 3)
  ggplot2::ggsave(filename = paste0("Dot_plot_pathways_", method, ".tiff"),
                  plot = dot_plots,
                  device = "jpeg",
                  path = file.path(save.dir),
                  width = 3*7,
                  height = 3*7,
                  units = "in",
                  dpi = 300,
                  bg = "white")
}