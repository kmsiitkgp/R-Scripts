

# data <- openxlsx::read.xlsx("C:/Users/KailasammS/Desktop/Prince_Mass_Spec.xlsx")
# data <- data %>% dplyr::mutate(Genes = make.unique(Genes, sep = "."))
# data <- data %>% dplyr::select(Genes, colnames(data)[grepl(pattern="KO|pos|neg", x=colnames(data))])
# 
# data <- data %>% 
#   dplyr::filter(!is.na(Genes)) %>% 
#   tibble::column_to_rownames("Genes")
# 
# data <- log(1+data, base=2)
# 
# data1 <- data[,1:6]
# data2 <- data[,7:12]
# 
# metadata1 <- data.frame("Sample" = colnames(data1), 
#                         "Condition" = c("Ypos","Ypos", "Ypos", "Yneg", "Yneg", "Yneg"))
# metadata2 <- data.frame("Sample" = colnames(data2), 
#                         "Condition" = c("Ypos","Ypos", "Ypos", "Yneg", "Yneg", "Yneg"))
# 
# 
# results_path <- "C:/Users/KailasammS/Desktop/"
# Target <- "Yneg"
# Reference <- "Ypos"
# metadata <- metadata1
# file_suffix <- "CRISPR"
# a <- calc_stats(data1, metadata, file_suffix)
# a <- a %>% dplyr::mutate(Gene = gsub(pattern="\\..*$", replacement="", x=Gene)) %>% 
#   dplyr::arrange(Gene) %>%
#   dplyr::group_by(Gene) %>%
#   dplyr::mutate(group_n = n()) %>%
#   dplyr::mutate(sig=case_when(padj <= 0.05 ~ 1, TRUE ~ 0)) %>%
#   dplyr::mutate(sig_n = sum(sig)) %>%
#   dplyr::select(everything(), -sig)
# # Save the results
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb, sheetName = "1")
# openxlsx::writeData(wb, sheet = "1", x = a, rowNames = FALSE)
# openxlsx::saveWorkbook(wb, 
#                        file = paste0(results_path, "Volcano_Results_", file_suffix, ".xlsx"), 
#                        overwrite = TRUE)
# 
# file_suffix <- "Natural"
# a <- calc_stats(data2, metadata2, file_suffix)
# a <- a %>% dplyr::mutate(Gene = gsub(pattern="\\..*$", replacement="", x=Gene)) %>% 
#   dplyr::arrange(Gene) %>%
#   dplyr::group_by(Gene) %>%
#   dplyr::mutate(group_n = n()) %>%
#   dplyr::mutate(sig=case_when(padj <= 0.05 ~ 1, TRUE ~ 0)) %>%
#   dplyr::mutate(sig_n = sum(sig)) %>%
#   dplyr::select(everything(), -sig)
# # Save the results
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb, sheetName = "1")
# openxlsx::writeData(wb, sheet = "1", x = a, rowNames = FALSE)
# openxlsx::saveWorkbook(wb, 
#                        file = paste0(results_path, "Volcano_Results_", file_suffix, ".xlsx"), 
#                        overwrite = TRUE)
# 
# 
# # Read Expr data from Depmap
# r <- read.csv(paste0(parent_path, "OmicsExpressionGenesExpectedCountProfile.csv")) %>% 
#   t() %>% 
#   as.data.frame() 
# 
# colnames(r) <- r[1,]
# r <- r[-1,]
# r <- r %>%
#   dplyr::mutate(across(.cols=everything(), .fns=as.numeric)) %>%
#   dplyr::mutate(across(.cols=everything(), .fns=round))
# 
# # Save expr data
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb, sheetName = "expr")
# openxlsx::writeData(wb, sheet = "expr", x = r,rowNames = TRUE)
# openxlsx::saveWorkbook(wb,
#                        file = paste0(results_path, "Readdata.xlsx"),
#                        overwrite = TRUE)
# 
# # Read metadata from Depmap
# m <- read.csv(paste0(parent_path, "Model.csv"))
# p <- read.csv(paste0(parent_path, "OmicsProfiles.csv")) %>%
#   dplyr::filter(Datatype == "rna")
# m <- m %>% dplyr::left_join(p, by=c("ModelID"="ModelID")) %>% 
#   dplyr::filter(ProfileID %in% colnames(r))
# 

# 
# #******************************************************************************#
# # Classify male patients into Y+ and Y-
# # Get Y genes from ensembl
# y_ensembl <- annotations %>% 
#   dplyr::filter(CHR == "Y") %>%
#   dplyr::filter(nchar(SYMBOL) > 0) %>%
#   # dplyr::filter(!grepl(pattern = "predicted|pseudogene|RIKEN|rRNA", x = DESCRIPTION)) %>%
#   # dplyr::arrange(SYMBOL) %>%
#   dplyr::select(SYMBOL) %>%
#   unlist(., use.names = FALSE)
# 
# # Find Y genes present in read data
# y_present <- y_ensembl[sapply(X=lapply(X=y_ensembl, FUN=grepl, x = rownames(r)),FUN=any)]
# 
# # keep expr data only for Y genes
# keep_r <- c()
# for (i in 1:length(y_ensembl)){
#   keep_r <- c(keep_r, rownames(r)[grepl(pattern=y_ensembl[i], x=rownames(r))])
# }
# y_expr <- r[keep_r,]
# y_expr <- y_expr[rowSums(y_expr) != 0,]
# 
# 
# # Find expr_data for y_genes of female cell lines
# y_expr_female <- y_expr[, m[m$Sex == "Female",]$ProfileID]
# # Calculate average expr of each Y gene in female
# y_expr_female[y_expr_female == 0] <- NA
# y_expr_average <- rowMeans(y_expr_female, na.rm=TRUE)
# #******************************************************************************#
# 
# 
# # Add YPos , Yneg classification by Mythreye
# 
# m_myhtree <- read.xlsx(paste0(parent_path, "!Metadata_Mythree.xlsx")) %>%
#   dplyr::select(ModelID, Condition)
# 
# meta_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Metadata.xlsx")) %>% 
#   dplyr::left_join(m_myhtree,by=c("ModelID"="ModelID")) %>%
#   dplyr::filter(!is.na(Condition)) %>%
  dplyr::distinct(pick(ModelID), .keep_all = TRUE)


#******************************************************************************#

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
  
# parent directory : directory where input files, results, etc are stored
parent_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Prince/"
results_path <- parent_path

#*****************Run GSEA on YKO vs scr RNASeq of MB49 cells******************#


#******Calculate normalized counts for YKO vs scr Proteomics of MB49 cells*****#
# 
# # Import raw counts
# raw_counts <- read.xlsx(paste0(parent_path, "Results_id_YKO_centromere_vs_scr_KO_centromere_pg_matrix.xlsx"))
# 
# # Remove rows that have multiple genes or missing genes
# raw_counts <- raw_counts[!grepl(pattern=";", x=raw_counts$Genes),] %>%
#   dplyr::filter(!is.na(Genes)) %>%
#   dplyr::arrange(Genes)
# rownames(raw_counts) <- seq(1:nrow(raw_counts))
# 
# # Create metadata
# # metadata MUST have "Sample" and "Condition" columns used by quantile_norm()
# metadata <- data.frame("Sample" = colnames(raw_counts)[-1], 
#                        "Condition" = c("SCR","SCR", "SCR", "YKO", "YKO", "YKO"))
# 
# # Notice that some genes have duplicates. This is because they correspond to 
# # different peptides detected by the mass spec algorithm. We keep the row that
# # contains the highest intensity for most samples and discard rest of duplicates
# dup_genes <- raw_counts %>% 
#   dplyr::count(Genes) %>% 
#   dplyr::filter(n>1) %>%
#   dplyr::select(Genes) %>%
#   unlist(use.names = FALSE)
# 
# # Create empty dataframe to store max intensities of duplicated genes
# df <- data.frame(matrix(NA, nrow = 1, ncol = ncol(raw_counts)))
# colnames(df) <- colnames(raw_counts)
# 
# # Iterate through each duplicated gene and find which row has max intensities
# # across multiple samples 
# for (gene in dup_genes){
#   
#   test <- raw_counts %>% 
#     dplyr::filter(Genes %in% gene) %>%
#     base::replace(is.na(.), 0)
#   
#   row_id <- c()
#   for (i in 2:ncol(test)){
#     row_id <- c(row_id, which.max(test[,i]))
#   }
#   
#   row_id <- as.data.frame(table(row_id)) %>% 
#     dplyr::slice_max(Freq) %>%  # if there are 2 max values both are selected
#     dplyr::select(row_id) %>%
#     unlist(use.names=FALSE) %>% 
#     as.character() %>%
#     as.numeric()
#   
#   df1 <- test[row_id,]
#   if (all(colnames(df) == colnames(df1))){
#     df <- dplyr::bind_rows(df,df1)
#   }
# }
# 
# # Remove the dummy column with NAs
# df <- df[-1,]
# 
# # If some genes still have 2 rows, average them
# df <- df %>% 
#   dplyr::group_by(Genes) %>% 
#   summarize(across(.cols=everything(), .fns=mean)) %>%
#   dplyr::ungroup()
# 
# # Remove duplicated genes from raw_counts    
# raw_counts <- raw_counts %>% 
#   dplyr::filter(!Genes %in% dup_genes)
# 
# # Add max values we calcualted to raw_counts
# if (all(colnames(df) == colnames(raw_counts))){
#   raw_counts <- dplyr::bind_rows(df, raw_counts)
# } 
# 
# # Replace all 0 with NA again
# raw_counts <- raw_counts %>% 
#   dplyr::mutate(across(.cols=where(is.numeric), ~na_if(., 0))) %>%
#   tibble::column_to_rownames("Genes")
#  
# # Perform quantile normalization before imputation
# quant_norm <- TRUE
# quant_norm_counts <- quantile_norm(raw_counts, metadata, quant_norm)
# 
# # Impute missing values in biological replicates using mean
# imputed_counts <- impute_with_mean(quant_norm_counts)
# 
# # Calculate padj and log2FoldChange
# Target <- "YKO" 
# Reference <- "SCR" 
# results <- calc_stats(imputed_counts, metadata, Target, Reference)
# 
# # Save the results
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb, sheetName = "results")
# openxlsx::writeData(wb, sheet = "results", x = results, rowNames = FALSE)
# openxlsx::addWorksheet(wb, sheetName = "raw_intensites")
# openxlsx::writeData(wb, sheet = "raw_intensites", x = raw_counts, rowNames=TRUE)
# openxlsx::addWorksheet(wb, sheetName = "imputed_intensities")
# openxlsx::writeData(wb, sheet = "imputed_intensities", x = imputed_counts, rowNames = FALSE)
# openxlsx::addWorksheet(wb, sheetName = "quant_norm_intensities")
# openxlsx::writeData(wb, sheet = "quant_norm_intensities", x = quant_norm_counts, rowNames = FALSE)
# openxlsx::saveWorkbook(wb, file = paste0(results_path, "pg_matrix.xlsx"),
#                        overwrite = TRUE)

###########Venn diagram of DEGs and DEPs

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
parent_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Prince/Venn_DEG_DEP/"
results_path <- parent_path

data_up <- read.xlsx(paste0(parent_path, "Venn_Input.xlsx"), sheet = "Sheet1")
data_down <- read.xlsx(paste0(parent_path, "Venn_Input.xlsx"), sheet = "Sheet2")

suffix_up <- "up"
suffix_down <- "down"
plot_venn(data_up, parent_path, suffix_up)
plot_venn(data_down, parent_path, suffix_down)

#***************Run GSEA on Transcriptomics and Proteomics of MB49 cells****************#

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
parent_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Prince/"
results_path <- parent_path

species <- "Mus musculus"
gmt_files <- list.files("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets/Mouse/", full.names = TRUE)
annotations <- get_annotations(species)

for (gmt_file in gmt_files){
  
  # Import DEGs dataframe with "SYMBOL", "padj" and "log2FoldChange" columns
  #DEGs_df <- read.xlsx(paste0(parent_path, "Natural_LOY/Natural_LOY_DEPs.xlsx"))
  DEGs_df <- read.xlsx(paste0(parent_path, "CRISPR_LOY/CRISPR_LOY_DEPs.xlsx"))
  #DEGs_df <- read.xlsx(paste0(parent_path, "CRISPR_LOY/Results_id_Y_Negative_vs_Y_Positive_DESeq2_modelled__DEGs.xlsx"))
  #DEGs_df <- read.xlsx(paste0(parent_path, "Natural_LOY/Results_id_Y_Negative_vs_Y_Positive_DESeq2_modelled__DEGs.xlsx"))
  
  fgsea(DEGs_df, gmt_file, annotations, results_path)
}

#***********GSEA RNA and protein overlap

files <- list.files(parent_path)
RNA_files <- files[grepl(pattern="_RNASeq",x=files)]
proteomics_files <- files[grepl(pattern="_Proteomics",x=files)]

rna_df <- data.frame(matrix(data=NA, nrow=1, ncol=10))
colnames(rna_df) <- c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge", "abs_NES", "Direction")

for (r in RNA_files){
  rna <- read.xlsx(paste0(parent_path, r))
  rna_df <- dplyr::bind_rows(rna_df, rna)
}

protein_df <- data.frame(matrix(data=NA, nrow=1, ncol=10))
colnames(protein_df) <- c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge", "abs_NES", "Direction")

for (r in proteomics_files){
  protein <- read.xlsx(paste0(parent_path, r))
  protein_df <- dplyr::bind_rows(protein_df, protein)
}

rna_df <- rna_df %>% dplyr::filter(padj <= 0.05)
#protein_df <- protein_df %>% dplyr::filter(padj <= 0.05)

rna_up <- rna_df %>% dplyr::filter(NES > 0) %>% dplyr::select(pathway) %>% unlist(use.names=FALSE)
rna_down <- rna_df %>% dplyr::filter(NES < 0) %>% dplyr::select(pathway) %>% unlist(use.names=FALSE)
protein_up <- protein_df %>% dplyr::filter(NES > 0) %>% dplyr::select(pathway) %>% unlist(use.names=FALSE)
protein_down <- protein_df %>% dplyr::filter(NES < 0) %>% dplyr::select(pathway) %>% unlist(use.names=FALSE)

list <- list(rna_up, rna_down, protein_up, protein_down)
names(list) <- c("rna_up", "rna_down", "protein_up", "protein_down")

# Function to merge list of differing lengths to a dataframe
merge_list_to_df <- function(list){
  
  max_len <- max(lengths(list))
  
  for (i in 1:length(list)){
    list[[i]] <- c(list[[i]], rep(x=NA, times=max_len-length(list[[i]])))
  }
  
  return(data.frame(list))
  
}

data <- merge_list_to_df(list)
#Run the script plot_venn.R

# #**********Gene and protein overlap*******************
# 
# rna_df <- read.xlsx(paste0(parent_path, "Results_id_YKO_centromere_vs_scr_KO_centromere_DESeq2_modelled_DEGs.xlsx"))
# protein_df <- read.xlsx(paste0(parent_path, "pg_matrix.xlsx"))
# 
# rna_df <- rna_df %>% dplyr::filter(padj <= 0.05)
# #protein_df <- protein_df %>% dplyr::filter(padj <= 0.05)
# 
# rna_up <- rna_df %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::select(SYMBOL) %>% unlist(use.names=FALSE)
# rna_down <- rna_df %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::select(SYMBOL) %>% unlist(use.names=FALSE)
# protein_up <- protein_df %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::select(SYMBOL) %>% unlist(use.names=FALSE)
# protein_down <- protein_df %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::select(SYMBOL) %>% unlist(use.names=FALSE)
# 
# list <- list(rna_up, rna_down, protein_up, protein_down)
# names(list) <- c("rna_up", "rna_down", "protein_up", "protein_down")
# 
# # Function to merge list of differing lengths to a dataframe
# merge_list_to_df <- function(list){
#   
#   max_len <- max(lengths(list))
#   
#   for (i in 1:length(list)){
#     list[[i]] <- c(list[[i]], rep(x=NA, times=max_len-length(list[[i]])))
#   }
#   
#   return(data.frame(list))
#   
# }
# 
# data <- merge_list_to_df(list)
# #Run the script plot_venn.R
# 
# 
# #******Find top pathways containing these common up and common down genes****#
# 
# common_up <- intersect(rna_up,protein_up)
# common_down <- intersect(rna_down,protein_down)
# 
# rna_df <- rna_df %>% dplyr::select(pathway, padj, NES, leadingEdge) 
# n_up <- c()
# n_down <- c()
# 
# for (i in 1:nrow(rna_df)){
#   g <- unlist(stringr::str_split(string=rna_df$leadingEdge[i], pattern = ","))
#   g <- gsub(pattern=" ", replacement ="", x=g)
#   up <- length(intersect(g, common_up))
#   down <- length(intersect(g, common_down))
#   n_up <- c(n_up, up)
#   n_down <- c(n_down, down)
# }
# 
# rna_df$n_up <- n_up
# rna_df$n_down <- n_down
# rna_df <- rna_df %>% 
#   dplyr::mutate(n_total = n_up+n_down) %>% 
#   dplyr::arrange(desc(n_total))
# 
# protein_df <- protein_df %>% dplyr::select(pathway, padj, NES, leadingEdge)
# n_up <- c()
# n_down <- c()
# 
# for (i in 1:nrow(protein_df)){
#   g <- unlist(stringr::str_split(string=protein_df$leadingEdge[i], pattern = ","))
#   g <- gsub(pattern=" ", replacement ="", x=g)
#   up <- length(intersect(g, common_up))
#   down <- length(intersect(g, common_down))
#   n_up <- c(n_up, up)
#   n_down <- c(n_down, down)
# }
# 
# protein_df$n_up <- n_up
# protein_df$n_down <- n_down
# protein_df <- protein_df %>% 
#   dplyr::mutate(n_total = n_up+n_down) %>% 
#   dplyr::arrange(desc(n_total))
# 
# # Save the results
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb, sheetName = "RNA_pathways")
# openxlsx::writeData(wb, sheet = "RNA_pathways", x = rna_df, rowNames = FALSE)
# openxlsx::addWorksheet(wb, sheetName = "Protein_pathways")
# openxlsx::writeData(wb, sheet = "Protein_pathways", x = protein_df, rowNames=FALSE)
# openxlsx::saveWorkbook(wb, file = paste0(results_path, "Pathways_with_common_genes.xlsx"),
#                        overwrite = TRUE)

#***************************Gene wise correlation

# Get normalized protein counts
protein_norm <- read.xlsx(paste0(parent_path, "pg_matrix.xlsx"))
protein_norm <- protein_norm[,c(1:7)]

# Get normalized rna counts
rna_norm <- read.xlsx(paste0(parent_path, "Normalized_Counts_DESeq2_modelled.xlsx"))
rna_norm <- rna_norm[,-c(2:5, 9:11)]
colnames(rna_norm) <- colnames(protein_norm)
rna_norm <- rna_norm %>% 
  dplyr::group_by(SYMBOL) %>% 
  dplyr::summarise(across(.cols=everything(), ~max(.))) %>%
  dplyr::ungroup()

common_genes <- intersect(protein_norm$SYMBOL, rna_norm$SYMBOL)

protein_norm <- protein_norm %>% 
  dplyr::filter(SYMBOL %in% common_genes) %>% 
  dplyr::mutate(SYMBOL = paste0(SYMBOL, "_p")) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  t() %>%
  data.frame()

rna_norm <- rna_norm %>% 
  dplyr::filter(SYMBOL %in% common_genes) %>% 
  dplyr::mutate(SYMBOL = paste0(SYMBOL, "_r")) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  t() %>%
  data.frame()

if (all(rownames(protein_norm) == rownames(rna_norm))){
  data_norm <- dplyr::bind_cols(protein_norm, rna_norm)
}
data_norm <- data_norm[, sort(colnames(data_norm))]

# Find Spearman correlation between protein and RNA expression for each gene
res_p <- c()
res_cor <- c()
# Genes like H2-D1_p will be converted t H2.D1_p in column names of my_data
for (g in make.names(common_genes)){
  
  my_data <- data_norm %>% 
    dplyr::select(c(paste0(g, "_p"), paste0(g, "_r")))
  
  res <- cor.test(my_data[,1], my_data[,2], method = "spearman")
  
  res_p <- c(res_p, res$p.value)
  res_cor <- c(res_cor, res$estimate[[1]])
}

cor_df <- data.frame(common_genes, res_cor, res_p) %>% 
  dplyr::filter(res_p<= 0.05)
# cor_df$padj <- stats::p.adjust(p=cor_df$res_p, method="fdr", n=length(cor_df$res_p))
# cor_df <- cor_df %>%
#   dplyr::filter(padj <= 0.05)
 

# Median of correlation co-efficients  
median(cor_df$res_cor)

# % of gene-protein showing positive correlation
nrow(cor_df %>% dplyr::filter(res_cor >0))*100/nrow(cor_df)

ggplot(data = cor_df, aes(x=res_cor)) +
  geom_histogram(aes(y=2.95*..count../sum(..count..)), position="identity", alpha=0.5, fill = "orange", color = "black") +
  geom_density(alpha=0.6, aes(y=..density..)) +
  theme_classic()



