# Define comparisons and design
Comparisons <- list(Variable             = c("Treatment",           "Condition",            "Condition"),
                    Target               = c("Growth.hormone",      "Line1_Growth.hormone", "Line2_Growth.hormone"),
                    Reference            = c("Control",             "Line1_Control",        "Line2_Control"),
                    lfc.cutoff           = 0,
                    padj.cutoff          = 0.1,
                    design               = c("Cell.Line+Treatment", "Condition",            "Condition"),
                    deseq2.batch.correct = FALSE, 
                    proj                 = "RNASeq_Vera",
                    species              = "Homo sapiens")

Comparisons <- list(Variable    = c("Condition"),
                    Target      = c("Y_Negative"),
                    Reference   = c("Y_Positive"),
                    lfc.cutoff  = 0,
                    padj.cutoff = 0.1,
                    design      = c("Condition"),
                    deseq2.batch.correct = FALSE,
                    proj        = "RNASeq_Hany",
                    species     = "Mus musculus")

Comparisons <- list(Variable    = c("Comparisons", "Comparisons"),
                    Target      = c("Quad.RF",     "Quad.RFL"),
                    Reference   = c("Quad.Fc",     "Quad.Fc"),
                    lfc.cutoff  = 0,
                    padj.cutoff = 0.1,
                    design      = c("Comparisons", "Comparisons"),
                    deseq2.batch.correct = FALSE,
                    proj        = "RNASeq_Sandrine.Supriya",
                    species     = "Mus musculus")



#******************************************************************************#
#               SINGLE REPLICATE DIFFERENTIAL EXPRESSION ANALYSIS              #
#******************************************************************************#

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Dr.Vera/"
species <- "Homo sapiens"
annotations <- get_annotations(species)

# Read the normalized counts
norm_count <- read.xlsx(paste0(data_path, "Norm_counts_DESeq2.xlsx"))

# Calculate log2FC
deg_file <- norm_count %>%
  tibble::column_to_rownames("ENSEMBL_ID") %>%
  dplyr::select(everything(), -c(SYMBOL, ENSEMBL_SYMBOL)) %>%
  dplyr::mutate(across(.cols=everything(), ~ log2(1+.))) %>%
  dplyr::mutate(Line1_FC1 = Line1_GH1-Line1_Control1,
                Line1_FC2 = Line1_GH2-Line1_Control2,
                Line1_FC3 = Line1_GH3-Line1_Control3,
                Line2_FC1 = Line2_GH1-Line2_Control1,
                Line2_FC2 = Line2_GH2-Line2_Control2,
                Line2_FC3 = Line2_GH3-Line2_Control3) %>%
  dplyr::mutate(log2FoldChange_L1 = (Line1_FC1+Line1_FC2+Line1_FC3)/3,
                log2FoldChange_L2 = (Line2_FC1+Line2_FC2+Line2_FC3)/3)

# Calculate pvalues
subset <- deg_file[,13:15]
t1 <- c()
for (i in 1:nrow(subset)){
  t <- stats::t.test(x = subset[i,], alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE)
  t1 <- c(t1,t$p.value)
}

subset <- deg_file[,16:18]
t2 <- c()
for (i in 1:nrow(subset)){
  t <- stats::t.test(x = subset[i,], alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE)
  t2 <- c(t2,t$p.value)
}

stats_df <- data.frame(ID = rownames(subset), pval_L1 = t1, pval_L2 = t2)
stats_df$padj_L1 <- stats::p.adjust(p = stats_df$pval_L1, method = "fdr", n = length(stats_df$pval_L1))
stats_df$padj_L2 <- stats::p.adjust(p = stats_df$pval_L2, method = "fdr", n = length(stats_df$pval_L2))

deg_file <- deg_file %>% 
  tibble::rownames_to_column("ID") %>%
  dplyr::left_join(stats_df, by=c("ID"="ID"))

# Get gene symbols
deg_file <- add_annotation(deg_file, annotations)

# Remove unwanted genes
deg_file <- deg_file %>%
  dplyr::filter(!(SYMBOL %in% c("Y_RNA", "Metazoa_SRP")))

# Identify genes that have consistent trend in expression upon GH treatment
deg_file <- deg_file %>%
  dplyr::mutate(Status = dplyr::case_when(Line1_FC1>0 & Line1_FC2>0 & Line1_FC3>0 & Line2_FC1>0 & Line2_FC2>0 & Line2_FC3>0 ~ "Up",
                                          Line1_FC1<0 & Line1_FC2<0 & Line1_FC3<0 & Line2_FC1<0 & Line2_FC2<0 & Line2_FC3<0 ~ "Down",
                                          TRUE ~ "Ambiguous"))

# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName="DEGs")
openxlsx::writeData(wb, sheet="DEGs", x=deg_file, rowNames=FALSE)
openxlsx::saveWorkbook(wb, file=paste0(data_path, "Similar.Trend.Genes.xlsx"), 
                       overwrite=TRUE)

#*************************OPTIONAL*********************************************#

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

proj <- "RNASeq_Vera"
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Dr.Vera/"
species <- "Homo sapiens"
annotations <- get_annotations(species)
file_suffix <- ""

# List pathways to plot
plot_pathways <- c("REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
                   "REACTOME_COLLAGEN_DEGRADATION",
                   "REACTOME_COLLAGEN_CHAIN_TRIMERIZATION",
                   "GOBP_CELL_SUBSTRATE_ADHESION",
                   "GOBP_CELL_MATRIX_ADHESION",
                   "GOBP_COLLAGEN_FIBRIL_ORGANIZATION",
                   "GOBP_CELL_SUBSTRATE_JUNCTION_ORGANIZATION",
                   "GOBP_NEGATIVE_REGULATION_OF_CELL_SUBSTRATE_ADHESION",
                   "GOCC_COLLAGEN_TRIMER",
                   "GOCC_COMPLEX_OF_COLLAGEN_TRIMERS",
                   "GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT_CONFERRING_TENSILE_STRENGTH")

plot_pathways <- c("REACTOME_DNA_DOUBLE_STRAND_BREAK_RESPONSE",
                   "GOBP_NEGATIVE_REGULATION_OF_DNA_REPAIR",
                   "GOBP_REGULATION_OF_DNA_REPAIR")

# Read ORA results file
ora_results <- read.xlsx(paste0(data_path, "Pathway_Analysis_Results_Similar.Trend.xlsx")) %>%
  dplyr::filter(Description %in% plot_pathways)

#**********************Plot enrichment plot of ORA pathways********************#

plot_ora(ora_results, file_suffix, data_path)

#**********************Plot heatmap of genes in pathways**********************#

# Get genes in these pathways
ora_df <- convert_ora_genelist_to_df(ora_results)
ora_list <- ora_df %>% unlist(use.names = FALSE) %>% unique()
ora_list <- ora_list[!is.na(ora_list)]

# Read normalized counts
norm_counts1 <- read.xlsx(paste0(data_path, "Norm_counts_DESeq2.xlsx")) %>%
  dplyr::select(SYMBOL, starts_with(c("Line1"))) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::filter(!(SYMBOL %in% c("Y_RNA", "Metazoa_SRP")))

# Read DEGs
df1 <- read.xlsx(paste0(data_path, "Similar.Trend.Genes.xlsx")) %>%
  dplyr::filter(pval_L1 < 0.05 | pval_L2 < 0.05)

# Heatmap parameters
plot_genes <- ora_list
disp_genes <- ora_list
file_suffix <- ""
output_path <- data_path
heatmap_params <- list(anno.row       = c("Group"),        # NULL, c("Group")
                       anno.column    = c("Treatment", "Cell.Line"),
                       row.split      = c("Group"),        # c(), NULL, NA     
                       col.split      = c("Treatment"),    # c(), NULL, NA
                       row.cluster    = c("group"),        # c("alphabetical", "group", "all")
                       col.cluster    = c("alphabetical"), # c("alphabetical", "group", "all")
                       log.transform  = TRUE,
                       scale          = TRUE,
                       border_color   = NA, #"white",
                       bar_width      = 10,              # NA , 5
                       bar_height     = 10,              # NA , 5
                       width          = 5,              # NA
                       height         = 5,              # NA 
                       matrix_color   = "rdbu",         # c("vrds", "rdbu")
                       expr_legend    = TRUE,           # set FALSE if it overlaps with annotation legends
                       file_format    = "tiff")

metadata_column <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, proj, "_Metadata.xlsx")) %>%
  dplyr::distinct_at("Sample.ID", .keep_all=TRUE)

metadata_row1 <- df1 %>% 
  dplyr::mutate(Group = Status) %>%
  dplyr::select(SYMBOL, Group) %>%
  dplyr::filter(Group != "Ambiguous",  SYMBOL %in% plot_genes) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE)

plot_heatmap(norm_counts1, metadata_column, metadata_row1, heatmap_params,
             plot_genes, disp_genes, "ORA.Pathway.Genes.L1", output_path)

#**********************Plot volcano of genes in pathways**********************#

# # Read DEGs
# df <- read.xlsx(paste0(data_path, "Similar.Trend.Genes.xlsx")) %>%
#   dplyr::filter(SYMBOL %in% ora_list)
# 
# # Define volcano plot parameters
# disp_genes <- c()
# file_suffix <- ""
# output_path <- data_path
# volcano_params <- list(Target          = "Treatment",
#                        Reference       = "Control", 
#                        lfc.cutoff      = 0.58,
#                        padj.cutoff     = 0.05)
# 
# # Input dataframe MUST have SYMBOL, padj, log2FoldChange columns
# L1 <- df %>% 
#   dplyr::select(SYMBOL, pval_L1, log2FoldChange_L1) %>%
#   dplyr::rename(SYMBOL = SYMBOL, padj = pval_L1, log2FoldChange = log2FoldChange_L1) %>%
#   dplyr::select(SYMBOL, padj, log2FoldChange)
# 
# plot_volcano(L1, volcano_params, disp_genes, "ORA.L1", data_path)
# 
# L2 <- df %>% 
#   dplyr::select(ENSEMBL_SYMBOL, pval_L2, log2FoldChange_L2) %>%
#   dplyr::rename(SYMBOL = ENSEMBL_SYMBOL, padj = pval_L2, log2FoldChange = log2FoldChange_L2) %>%
#   dplyr::select(SYMBOL, padj, log2FoldChange)
# 
# plot_volcano(L2, volcano_params, disp_genes, "ORA.L2", data_path)

#******************************************************************************#








# # Plot enrichment plots for significant pathways (padj < 0.05)
# sig_pathways <- gsea_results %>% 
#   dplyr::filter(padj < 0.05) %>% 
#   dplyr::slice_max(abs_NES, n = 12) %>%
#   dplyr::select(pathway) %>% 
#   unlist(use.names=FALSE)

# if (length(sig_pathways) > 0){
#   for (p in sig_pathways){
#     
#     # Create Enrichment plots for all significant pathways
#     fgsea::plotEnrichment(pathway = gsea_results %>%
#                             dplyr::filter(pathway == p) %>%
#                             dplyr::select(leadingEdge) %>%
#                             unlist(., use.names = FALSE),
#                           stats = DEGs_list,
#                           gseaParam = 1,
#                           ticksSize = 0.2)
#     
#     ggplot2::ggsave(filename = paste0("GSEA_NES_Plot_", p, ".tiff"),
#                     plot = last_plot(),
#                     device = "jpeg",
#                     path = data_path,
#                     scale = 1,
#                     width = 6,
#                     height = 7,
#                     units = c("in"),
#                     dpi = 600,
#                     limitsize = TRUE,
#                     bg = NULL)
#     
#     # # Find geneSetID corresponding to pathway in gseaobject from clusterprofiler
#     # enrichplot::gseaplot2(x = gsea,
#     #                       geneSetID = 1,
#     #                       title = p,
#     #                       color = "green",
#     #                       base_size = 11,
#     #                       rel_heights = c(1.5, 0.5, 1),
#     #                       subplots = 1:3,
#     #                       pvalue_table = FALSE,
#     #                       ES_geom = "line")
#   }
#   
#   # Plot bar plots for top 12 significant pathways by NES (padj < 0.05)
#   # It is difficult to control the width of bars in ggplot. Since, we plot
#   # top 12 pathways, we insert dummy entries to make the data frame have 12
#   # pathways "if" the data frame has less than 12 pathways
#   if (length(sig_pathways) < 12){
#     gsea_results <- gsea_results %>% 
#       dplyr::filter(padj < 0.05)
#     nrows <- nrow(gsea_results)
#     dummy_pathway   <- paste0("None.", seq(1:(12-nrows)))
#     dummy_NES       <- rep(0, times=12-nrows)
#     dummy_Direction <- rep("Downregulated", times=12-nrows)
#     dummy_df <- data.frame(pathway = dummy_pathway, NES = dummy_NES, Direction = dummy_Direction)
#     gsea_results <- dplyr::bind_rows(gsea_results, dummy_df)  
#   } else {
#     gsea_results <- gsea_results %>% 
#       dplyr::filter(padj < 0.05) %>% 
#       dplyr::slice_max(abs_NES, n = 12)
#   }
#   
#   # Modify pathway names to make the plot pretty
#   gsea_results_pretty <- gsea_results %>% 
#     data.frame() %>%
#     dplyr::mutate(pathway = gsub("HALLMARK_|SA_|SIG_|NABA_|GOBP_|GOMF_", "", pathway),
#                   pathway = gsub("_", " ", pathway),
#                   pathway = gsub("ENDOPLASMIC RETICULUM", "ER", pathway),
#                   #pathway = stringr::str_trunc(pathway, 45, "right"),
#                   #pathway = stringr::str_to_title(pathway),
#                   #length = stringr::str_length(pathway),
#                   pathway = stringr::str_wrap(pathway, width = 22))
#   
#   ggplot2::ggplot(data = gsea_results_pretty,
#                   aes(x = NES, y = reorder(pathway, NES), fill = Direction)) +
#     # fill = direction means direction will be arranged in alphabetical order.
#     # So, if you had labeled direction as "Upregulated" and "Downregulated",
#     # then first color in scale_fill_manual() will be assigned to
#     # "downregulated" and it will be labeled as "Activated in Males" in the
#     # plot. So, be careful.
#     ggplot2::geom_col(width = 0.75) +
#     ggplot2::theme_classic() +
#     ggplot2::labs(x = "Normalized Enrichment Score(NES)",
#                   y = "",
#                   title = "GSEA",
#                   fill = "") +
#     ggplot2::coord_cartesian(xlim = c(floor(-max(abs(gsea_results$NES), na.rm=TRUE)), ceiling(max(abs(gsea_results$NES), na.rm=TRUE)))) +
#     ggplot2::theme(#aspect.ratio = 2,
#       plot.title =   element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
#       plot.caption = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0),
#       axis.title.x = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
#       axis.title.y = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
#       axis.text.x =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
#       axis.text.y =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 1),
#       legend.title = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
#       legend.text =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
#       #legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
#       legend.position = "bottom",
#       legend.justification = "left",
#       legend.direction = "horizontal",
#       legend.key.height= unit(0.5, 'cm'),
#       legend.key.width= unit(1.25, 'cm')) +
#     ggplot2::scale_fill_manual(labels=c("Upregulated", "Downregulated"),
#                                values = c(RColorBrewer::brewer.pal(11, "RdYlBu")[c(1)],
#                                           RColorBrewer::brewer.pal(11, "RdYlBu")[c(11)]))
#   
#   # Save the plot
#   ggplot2::ggsave(filename = paste0("GSEA_", gmt_name, "_", file_suffix, ".tiff"),
#                   plot = last_plot(),
#                   device = "jpeg",
#                   path = data_path,
#                   scale = 1,
#                   width = 6,
#                   height = 7,
#                   units = c("in"),
#                   dpi = 600,
#                   limitsize = TRUE,
#                   bg = NULL)
# }