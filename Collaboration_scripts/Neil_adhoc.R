
### Josh data ####
proj.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Neil_adhoc/"
gmt.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets"
save.dir <- proj.dir

proj.params <- list(
  proj                    = "Neil_adhoc",               # Project name
  species                 = "Homo sapiens"              # Species name (e.g. "Mus musculus", "Homo sapiens")
)

data <- read.xlsx(file.path(proj.dir, "DEGs_LMvsOthers_PO1.xlsx"), sheet= "DEGs_LMvsOthers_PO1ALL")
  
DEGs_df <- data %>%
  dplyr::rename(SYMBOL = gene, padj = adjp_forVolPlot, log2FoldChange = avg_log2FC)
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

# ---- Perform TF analysis on DEGs using t-statistics ----
t_stats_mat <- DEGs_df %>% 
  as.data.frame() %>%
  dplyr::mutate(t = -log10(padj) * log2FoldChange) %>%
  dplyr::filter(!is.na(t)) %>%
  dplyr::select(SYMBOL, t) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  as.matrix()

tf_res_degs <- tf_analysis(t_stats_mat, proj.params$species)

# ---- Plot barplots of top TFs from DEGs ----
# NOTE: wsum returns wsum, norm_wsum and corr_wsum.
# wsum (DONT USE): Biased toward larger gene sets (more genes → bigger sum)
# norm_wsum (USE): Adjusts for pathway length so small and large gene sets are comparable
# corr_sum (USE): corrects for high correlation as it can make enrichment appear stronger
stats_degs <- unique(tf_res_degs$all_tfs$statistic)
bar_plots <- list()

# Extract treatment and control from current contrast string
contrast_split <- stringr::str_split(string = "LM-Rest", pattern = "-")[[1]]
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
ggsave(filename = file.path(save.dir, "Bar_plot_TFs.jpeg"),
       plot = bar_plots_combined,
       device = "jpeg",
       width = 10,
       height = 6 * 3,
       units = "in",
       dpi = 300,
       bg = "white")




### Spatial data ####
# ---- VOlcano plot macrophages MAST test ----
save.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Neil_Spatial"
DEGs_df <- read.xlsx(file.path(save.dir, "deseq2/Seurat_DEGs.xlsx"),
                     sheet = "DEGs_MAST") %>%
  dplyr::filter(Comparison == "Macrophages.Metastatic.Tumor.Pancreas.vs.Normal.Liver") %>%
  dplyr::mutate(log2FoldChange = avg_log2FC,
                padj = p_val_adj)

volcano.params <- list(volcano.padj.cutoff = 0.05,
                       volcano.lfc.cutoff = 0.58)
plot_volcano(DEGs_df, volcano.params, contrast = "Met.Macrophage-Liver.Macrophage", save.dir)

# ---- VOlcano plot macrophages Wilcoxon test ----
DEGs_df <- read.xlsx(file.path(save.dir, "deseq2/Seurat_DEGs.xlsx"),
                     sheet = "DEGs_Wilcox") %>%
  dplyr::filter(Comparison == "Macrophages.Metastatic.Tumor.Pancreas.vs.Normal.Liver") %>%
  dplyr::mutate(log2FoldChange = avg_log2FC,
                padj = p_val_adj)

volcano.params <- list(volcano.padj.cutoff = 0.05,
                       volcano.lfc.cutoff = 0.58)
plot_volcano(DEGs_df, volcano.params, contrast = "Met.Macrophage-Liver.Macrophage", save.dir)


# ---- Heatmap macrophages VST counts ----
disp_genes <- c("SPP1", "C1QA", "C1QB", "C1QC", "MMP19", "SDC2")
metadata_col <- NULL
metadata_row <- NULL
norm_counts <- read.xlsx(file.path(save.dir, "deseq2/Seurat_VST_Counts.xlsx"),
                         sheet = "Sheet26") %>%
  dplyr::filter(SYMBOL %in% disp_genes) %>% 
  tibble::column_to_rownames("SYMBOL")
heatmap.params <- list(heatmap.force.log       = FALSE,                       # Force log transform on heatmap data (default FALSE, auto detect)
                       heatmap.col.ann         = NULL,                        # Columns from metadata used as column annotation
                       heatmap.row.ann         = NULL,                        # Columns from metadata_row for row annotation
                       heatmap.col.gaps        = NULL,                        # Columns to define gaps in heatmap columns
                       heatmap.row.gaps        = NULL,                        # Columns to define gaps in heatmap rows
                       heatmap.col.cluster     = "all",                       # Clustering for columns ("all", "alphabetical", or metadata column)
                       heatmap.row.cluster     = "all",                       # Clustering for rows ("all", "alphabetical", or metadata column)
                       heatmap.palette         = "rdbu",                      # Color palette for heatmap ("rdbu" or "vrds")
                       heatmap.ann.palette     = "discrete",                  # Annotation palette type ("discrete" or "sequential")
                       heatmap.border.color    = NA,                          # Border color of heatmap cells, NA for no border
                       heatmap.show.expr.legend= TRUE,                        # Show expression legend on heatmap (set FALSE if overlapping annotations)
                       heatmap.title           = NA,                          # Title for heatmap (default NA = no title)
                       heatmap.format          = "tiff"
)
ph <- plot_heatmap(norm_counts, heatmap.params, metadata_col, metadata_row, disp_genes)
# Remove whitespace around heatmap
ph$ph$gtable$widths <- grid::unit.pmax(ph$ph$gtable$widths, grid::unit(0.1, "cm"))
ph$ph$gtable$heights <- grid::unit.pmax(ph$ph$gtable$heights, grid::unit(0.1, "cm"))
heatmap_plots <- list()
if (!is.null(ph) && !is.null(ph$ph) && !is.null(ph$ph$gtable)) {
  heatmap_plots <- c(heatmap_plots, list(ph$ph$gtable))
}
n_plots <- length(heatmap_plots)
jpeg(filename = file.path(save.dir, "Heatmap.jpeg"),
     width = 7 * n_plots,
     height = 10,
     units = "in",
     res = 300)
gridExtra::grid.arrange(grobs = heatmap_plots, ncol = n_plots, nrow = 1)
dev.off()

save_xlsx(ph$mat, file.path(save.dir, "Heatmap_Matrix.xlsx"), "Heatmap_matrix", row_names = TRUE)

# ---- Heatmap all Epithelial VST counts ----

assay <- "Spatial.016um"
ctoi <- c("Epithelial", "Met.Panc", "Primary.Panc", "Pancreatic.Acinar", "Pancreatic.Islet", "Hepatocytes")
coi <- integrated.seurat@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::filter(Cell.Type %in% ctoi) %>%
  
  # Count cells per group (Treatment + Cell.Type + Patient)
  dplyr::add_count(Treatment, Cell.Type, Patient, name = "cells_per_group") %>%
  dplyr::filter(cells_per_group >= 100) %>%
  
  # Count number of patients per group (Treatment + Cell.Type)
  dplyr::group_by(Treatment, Cell.Type) %>%
  dplyr::mutate(replicates_per_group = n_distinct(Patient)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(replicates_per_group >= 2) %>%
  
  # Keep all barcodes passing the filters
  dplyr::pull(barcode)

# Subset to keep only cells of interest
subset_obj <- integrated.seurat[, rownames(integrated.seurat@meta.data) %in% coi]
subset_obj@meta.data$Comparisons <- subset_obj@meta.data %>%
  dplyr::mutate(Comparisons = paste(Treatment, Cell.Type, Patient, sep="_")) %>%
  dplyr::pull(Comparisons)

# Aggregate expression
deg_obj <- Seurat::AggregateExpression(
  object = subset_obj,
  group.by = c("Treatment", "Cell.Type", "Patient", "Comparisons"),
  assays = assay,
  slot = "counts",
  return.seurat = TRUE
)

# Create dds
raw_counts <- as.data.frame(as.matrix(GetAssayData(deg_obj, assay = assay, slot = "counts")))
meta_data <- deg_obj@meta.data[colnames(raw_counts), , drop = FALSE]
rownames(meta_data) <- colnames(raw_counts)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts,
                                      colData = meta_data,
                                      design = ~Cell.Type)

# --- Generate and Annotate Blinded VST Counts ---
# NOTE: vst counts are affected by design ONLY when blind=FALSE
vsd_blind <- DESeq2::vst(dds, blind = TRUE)
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

# --- Perform LRT Test for Entire Dataset ---
dds_LRT <- DESeq2::DESeq(dds, test = "LRT", reduced = ~1)
res_LRT <- DESeq2::results(dds_LRT) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  add_annotation()

# --- Identify Significant Genes for Heatmap ---
sig_genes <- res_LRT %>%
  dplyr::filter(padj <= 0.05) %>%
  dplyr::pull(SYMBOL)

norm_counts <- vst_counts_blind[rownames(vst_counts_blind) %in% sig_genes, ]
metadata_col <- meta_data %>%
  dplyr::mutate(Sample.ID = rownames(.))
proj.params$heatmap.title <- ""

# --- Plot and Save Heatmap Using Helper ---
heatmap.params <- list(heatmap.force.log       = FALSE,                       # Force log transform on heatmap data (default FALSE, auto detect)
                       heatmap.col.ann         = c("Treatment", "Cell.Type", "Patient"), 
                       heatmap.row.ann         = NULL,                        # Columns from metadata_row for row annotation
                       heatmap.col.gaps        = c("Cell.Type"),
                       heatmap.row.gaps        = NULL,                        # Columns to define gaps in heatmap rows
                       heatmap.col.cluster     = c("Cell.Type"),              # Clustering for columns ("all", "alphabetical", or metadata column)
                       heatmap.row.cluster     = "all",                       # Clustering for rows ("all", "alphabetical", or metadata column)
                       heatmap.palette         = "rdbu",                      # Color palette for heatmap ("rdbu" or "vrds")
                       heatmap.ann.palette     = "discrete",                  # Annotation palette type ("discrete" or "sequential")
                       heatmap.border.color    = NA,                          # Border color of heatmap cells, NA for no border
                       heatmap.show.expr.legend= TRUE,                        # Show expression legend on heatmap (set FALSE if overlapping annotations)
                       heatmap.title           = NA,                          # Title for heatmap (default NA = no title)
                       heatmap.format          = "tiff"
)
ph <- plot_heatmap(norm_counts, heatmap.params, metadata_col)

jpeg(filename = file.path(proj.params$deseq2_dir, "Heatmap_CellType.jpeg"),
     width = 14, height = 10, units = "in", res = 300)
gridExtra::grid.arrange(grobs = list(ph$ph$gtable), ncol = 1, nrow = 1)
dev.off()

save_xlsx(ph$mat, file.path(proj.params$deseq2_dir, "Heatmap_Matrix_CellType.xlsx"), "Heatmap_matrix", row_names = TRUE)






# Heatmap FGFR
path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Neil_Spatial/Metastatic.Tumor.Pancreas.Met.Panc-Primary.Tumor.Pancreas.Primary.Panc"
norm_counts <- read.xlsx(file.path(path, "VST_counts.xlsx")) %>%
  dplyr::select(SYMBOL, contains("Hepatocytes")) %>%
  dplyr::filter(SYMBOL %in% c("FGFR1", "FGFR2", "FGFR3", "FGFR4")) %>%
  tibble::column_to_rownames("SYMBOL")

metadata_col <- data.frame(Type = c(rep("Normal", 4), rep("Metastatic", 3)),
                           Sample_ID = colnames(norm_counts))


heatmap.params <- list(force.log       = FALSE,                       # Force log transform on heatmap data (default FALSE, auto detect)
                       col.ann         = c("Type"), 
                       row.ann         = NULL,                        # Columns from metadata_row for row annotation
                       col.gaps        = "Type",
                       row.gaps        = NULL,                        # Columns to define gaps in heatmap rows
                       col.cluster     = "Type",                       # Clustering for columns ("all", "alphabetical", or metadata column)
                       row.cluster     = "alphabetical",                       # Clustering for rows ("all", "alphabetical", or metadata column)
                       palette         = "rdbu",                      # Color palette for heatmap ("rdbu" or "vrds")
                       ann.palette     = "discrete",                  # Annotation palette type ("discrete" or "sequential")
                       border.color    = NA,                          # Border color of heatmap cells, NA for no border
                       show.expr.legend= TRUE,                        # Show expression legend on heatmap (set FALSE if overlapping annotations)
                       title           = NA,                          # Title for heatmap (default NA = no title)
                       format          = "tiff"
)
proj.params <- c(list(heatmap = heatmap.params))
ph <- plot_heatmap(norm_counts, proj.params, metadata_col)

jpeg(filename = file.path(path, "Heatmap_FGFR.jpeg"),
     width = 14, height = 10, units = "in", res = 300)
gridExtra::grid.arrange(grobs = list(ph$ph$gtable), ncol = 1, nrow = 1)
dev.off()

save_xlsx(ph$mat, file.path(path, "Heatmap_FGFR.xlsx"), "Heatmap_matrix", row_names = TRUE)






















