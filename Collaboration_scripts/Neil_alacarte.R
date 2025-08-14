proj.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Neil_alacarte/"
gmt.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets"
save.dir <- proj.dir

proj.params <- list(
  proj                    = "Neil_alacarte",               # Project name
  species                 = "Homo sapiens"               # Species name (e.g. "Mus musculus", "Homo sapiens")
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
