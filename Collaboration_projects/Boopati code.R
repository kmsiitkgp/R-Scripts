source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

#*************** Metabolite analysis

# USe metaboanalyst to normalize by sum, log10 transform, perform pareto-scaling,
# calculate fold change, pvalue, FDR as described in Nature paper

#************************GSEA on bulk RNASeq

parent_path <- "C:/Users/kailasamms/Box/Boopati_paper/Transcriptomics/"
results_path <- parent_path

species <- "Mus musculus"
gmt_files <- list.files("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets/Mouse/", full.names = TRUE)
annotations <- get_annotations(species)

for (gmt_file in gmt_files){
  
  #DEGs_df <- read.xlsx(paste0(parent_path, "CRISPR_LOY/Results_id_Y_Negative_vs_Y_Positive_DESeq2_modelled__DEGs.xlsx"))
  #DEGs_df <- read.xlsx(paste0(parent_path, "Natural_LOY/Results_id_Y_Negative_vs_Y_Positive_DESeq2_modelled__DEGs.xlsx"))
  DEGs_df <- read.xlsx(paste0(parent_path, "DDR2KO/Results_id_KO_vs_Control_DESeq2_modelled__DEGs.xlsx"))
  
  fgsea(DEGs_df, gmt_file, annotations, results_path)
}

#*********************Volcano plots of metabolites

parent_path <- "C:/Users/kailasamms/Box/Boopati_paper/Metabolomics/"
results_path <- parent_path










# data <- read.xlsx(paste0(parent_path, "Data_metabolomics_DDR2KO_vs_SCR.xlsx")) %>%
#   tibble::column_to_rownames("compound") %>%
#   dplyr::select(everything(), -Pathway)
# 
# meta_data <- data.frame("Sample" = colnames(data), 
#                         "Condition" = c("Y+","Y+", "Y+", "Y-", "Y-", "Y-"))
# 
# 
# # Plot box plot before normalization
# b_data <- log(1+data, base=2) %>%
#   tibble::rownames_to_column("compound") %>%
#   tidyr::pivot_longer(cols = !compound, names_to="Sample", values_to="Intensity")
# boxplot(b_data$Intensity ~ b_data$Sample)
# 
# 
# # Perform quantile normalization
# quant_norm <- TRUE
# n_data <- quantile_norm(data, meta_data, quant_norm)
# 
# # Plot box plot before normalization
# b_data <- log(1+n_data, base=2) %>%
#   tibble::rownames_to_column("compound") %>%
#   tidyr::pivot_longer(cols = !compound, names_to="Sample", values_to="Intensity")
# boxplot(b_data$Intensity ~ b_data$Sample)
# 
# # Calculate padj and log2FC
# Target <- "Y-"
# Reference <- "Y+"
# n_data <- log(1+n_data, base=2)
# log2_transformed_already <- FALSE
# f_data <- calc_stats(n_data, meta_data, Target, Reference,log2_transformed_already)