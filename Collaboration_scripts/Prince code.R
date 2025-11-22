

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

files <- list.files(paste0(parent_path,"CRISPR_LOY"), full.names = TRUE)
files <- list.files(paste0(parent_path,"Natural_LOY"), full.names = TRUE)
RNA_files <- files[grepl(pattern="_RNA",x=files)]
proteomics_files <- files[grepl(pattern="_Protein",x=files)]

rna_df <- data.frame(matrix(data=NA, nrow=1, ncol=10))
colnames(rna_df) <- c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge", "abs_NES", "Direction")

for (r in RNA_files){
  rna <- read.xlsx(r)
  rna_df <- dplyr::bind_rows(rna_df, rna)
}

protein_df <- data.frame(matrix(data=NA, nrow=1, ncol=10))
colnames(protein_df) <- c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge", "abs_NES", "Direction")

for (r in proteomics_files){
  protein <- read.xlsx(r)
  protein_df <- dplyr::bind_rows(protein_df, protein)
}

rna_df <- rna_df %>% dplyr::filter(padj <= 0.05)
protein_df <- protein_df %>% dplyr::filter(pval <= 0.05)

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
data_crispr <- data
data_natural <- data

# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "cripr_pathways")
openxlsx::writeData(wb, sheet = "cripr_pathways", x = data_crispr, rowNames = TRUE)
openxlsx::addWorksheet(wb, sheetName = "natural_pathways")
openxlsx::writeData(wb, sheet = "natural_pathways", x = data_natural, rowNames=TRUE)
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Pathways.xlsx"), overwrite = TRUE)

###########Venn diagram of DEGs and DEPs

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
parent_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Prince/Venn_DPG_DPE/"
results_path <- parent_path

data_up <- read.xlsx(paste0(parent_path, "Venn_Input.xlsx"), sheet = "Sheet1")
data_down <- read.xlsx(paste0(parent_path, "Venn_Input.xlsx"), sheet = "Sheet2")

suffix_up <- "up"
suffix_down <- "down"
plot_venn(data_up, parent_path, suffix_up)
plot_venn(data_down, parent_path, suffix_down)

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


#**********DTKPA2, DTKPB1 method comparisons

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

parent_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/DEG_results/"
results_path <- parent_path
libraries <- c("DTKPA2", "DTKPB1")

# Save to excel file
wb <- openxlsx::createWorkbook()

for (l in libraries){
  
  # Create empty dataframe to store hits
  df <- data.frame(Genes = c(rep(x="Gene", times=1000)))
  
  files <- list.files(paste0(parent_path, l), full.names = TRUE)
  for (f in files){
    
    print(f)
    
    # Create column name from filename
    sample <- gsub(parent_path, "", f)
    sample <- gsub(l, "", sample)
    sample <- gsub(".gene_summary.txt", "", sample)
    sample <- gsub("\\/", "", sample)
    
    # Read file, filter and select gene names
    data <- read.table(f, header=TRUE)
    data <- data %>%
      dplyr::filter(neg.goodsgrna >=4, neg.p.value <= 0.05) %>%
      dplyr::select(id) %>%
      dplyr::rename(!!rlang::sym(sample) := id)
    
    # Increase number of rows to 1000 so all columns have same length
    data[(nrow(data)+1):1000,1] <- NA
    
    # Join hits
    df <- dplyr::bind_cols(df, data)
  }
  
  openxlsx::addWorksheet(wb, sheetName = l)
  openxlsx::writeData(wb, sheet = l, x = df, rowNames = FALSE)
}
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Method comparison.xlsx"),
                       overwrite = TRUE)

# I checked the hits for all 9 approaches
# (i) Total - median, alphamedian, second best
# (ii) Median - median, alphamedian, second best
# (iii) Control - median, alphamedian, second best

# Filter genes with pvalue <=0.05, neg.good.sgrna >=4
# There is no difference in hits between median, alphamedian and secondbest other than logFC
# There is 10%-20% difference in hits between total vs median vs control
# Stick to median normalization. I notice with control or total normalization, there
# are additional genes that show up as hits but I dont see difference in normalized counts
# Median normalization and alphamedian log2FC seems best

#*************scRNASeq Simon

proj <- "scRNASeq_Simon"
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")

reduc <- "harmony"
integrated_seurat <- base::readRDS(paste0(seurat_results, "Simon_integrated_seurat_snn.rds"))
integrated_seurat <- subset(x = integrated_seurat,
                            Sex == "Male")

# # Plot UMAP of celltypes and subtypes
# Seurat::DimPlot(object=integrated_seurat,
#                 reduction="umap",
#                 cols=my_palette,
#                 label=FALSE,
#                 group.by="celltype", #"subtype", 
#                 split.by=NULL,#split,
#                 shape.by=NULL,
#                 pt.size=0.2,
#                 label.size=5,
#                 repel=FALSE,
#                 raster=FALSE) +
#   ggplot2::labs(fill="CLUSTERS",
#                 x="UMAP_1",
#                 y="UMAP_2") +
#   my_theme

y_genes <- c("DDX3Y", "EIF1AY", "HSFY2", "KDM5D", "UTY", "NLGN4Y", 
             "PCDH11Y", "RPS4Y1", "TBL1Y", "TMSB4Y", "USP9Y", "ZFY", 
             "DAZ1", "DAZ2", "DAZ3", "DAZ4", "PRY2", "RBMY1A1")
y_genes <- rownames(integrated_seurat@assays$RNA$data)[tolower(rownames(integrated_seurat@assays$RNA$data)) %in% 
                                                         tolower(y_genes)]
integrated_seurat <- Seurat::AddModuleScore(object=integrated_seurat,
                                            features=list(y_genes),
                                            assay="RNA",
                                            slot="data",
                                            name="Yscore")

df <- integrated_seurat@assays$RNA$data
df <- df[y_genes,]
df <- df[rowSums(df) != 0,]
df[df>0] <- 1
ypos_cells <- colnames(df)[colSums(df)!=0]
yneg_cells <- colnames(df)[colSums(df)==0]

metadata <- integrated_seurat@meta.data %>%
  dplyr::mutate(Ystatus = dplyr::case_when(barcode %in% ypos_cells ~ "Ypos",
                                           barcode %in% yneg_cells ~ "Yneg"),
                Ystatus_score = dplyr::case_when(Yscore1 > 0 ~ "Ypos",
                                                 TRUE ~ "Yneg"))

integrated_seurat@meta.data <- metadata

# Ystatus classifies cell as Ypos and Yneg based on normalized counts
# Ystatus_score classifies cells as Ypos and Yneg based on the modulescore of 18 Y genes
# metadata %>% dplyr::count(Ystatus, Ystatus_score) clearly shows all Yneg cells
# based on Ystatus are also Y neg based on Ystatus_score, but some cells expressing
# low levels of Y genes are classified as Yneg by Ystatus_score but Ypos by 
# Ystatus.

# Plot UMAPs of individual genes
features <- c("AARS1", "ABAT", "ACP4","APRT", "ASRGL1", "ATM", "ATOX1", 
              "ATP5F1C", "BDKRB2", "BHMT2", "BMPR2", "CCR5", "CDH5", "CDK5R2", 
              "CDKL5","CHERP", "CHKB", "CHUK", "CKM", "CLPS", "CNKSR3", "COMMD8",
              "CTSH", "CUTA", "DAPK3", "DARS1", "DARS2", "DCAF11", "DESI1", 
              "DUS22", "EGLN3", "F5", "FBP2", "FGFR2", "G6PC2", "GABRA1", "GAPDH",
              "GLRB", "GNGT1", "GP1BA", "GRIA1", "HASPIN", "GSK3B", "HLA-A", 
              "HLA-B", "HLA-C", "HINT1", "HPRT1", "HS3ST3A1", "INPP5K", "KLK6",
              "LIPA", "LSM2", "MAN1B1", "MAP3K19", "MMUT", "NAGS", "NDUFA10",
              "NEK11", "OGFOD1", "PBK", "PDHA1", "PHKA1", "PLEKHA4", "PPP1R14C",
              "PPP1R9B", "PPP4R1", "PRKAR2A", "PRKD2", "PTK2B", "RDH14", 
              "RETSAT", "RIPK1", "SDHA", "SLC19A3", "SLC25A4", "SPHKAP",
              "TNFSF13B", "TNK1", "TRIB3","TSPO", "POLR1F", "YIPF3")
features <- rownames(integrated_seurat@assays$RNA$data)[tolower(rownames(integrated_seurat@assays$RNA$data)) %in% 
                                                          tolower(features)]
feature_plot <- function(i){
  
  split <- "Condition" #NULL
  split <- NULL
  Seurat::FeaturePlot(object=integrated_seurat,
                      slot="data",
                      features=i,
                      split.by=split,
                      #cols= c("grey", viridis(n=10, option="C", direction=-1)),
                      pt.size=0.4,
                      order=TRUE,
                      min.cutoff='q10',
                      reduction="umap",
                      label=TRUE,
                      combine=TRUE,
                      raster=FALSE) +
    scale_colour_gradientn(colours=rev(brewer.pal(n=11, name="RdBu"))[5:11])
}
for (n in 1:ceiling((length(features)/10))){
  j <- 10*n-9
  k <- 10*n
  vec <- features[j:k]
  purrr::map(.x=vec[!is.na(vec)], 
             .f=feature_plot) %>%
    cowplot::plot_grid(plotlist=.,
                       align="hv",
                       axis="tblr",
                       nrow=2,
                       ncol=5,
                       rel_widths=1,
                       rel_heights=1,
                       greedy=TRUE,
                       byrow=TRUE)
  
  ggplot2::ggsave(filename=paste0("Feature_plot_",n, ".jpg"),
                  plot=last_plot(),
                  device="jpeg",
                  #path=diagnostics_path,
                  width=8.5*4,
                  height=11*2,
                  units=c("in"),
                  dpi=300,
                  limitsize=FALSE,
                  bg="white")
}

# Plot UMAP of module score
integrated_seurat <- Seurat::AddModuleScore(object=integrated_seurat,
                                            features=list(features),
                                            assay="RNA",
                                            slot="data",
                                            name="LOY_sig")

Seurat::FeaturePlot(object=subset(integrated_seurat, Ystatus == "Yneg"),
                    slot="data",
                    features="LOY_sig1",
                    #cols= c("grey", viridis(n=10, option="C", direction=1)),
                    pt.size=0.4,
                    order=TRUE,
                    min.cutoff='q10',
                    reduction="umap",
                    label=FALSE,
                    combine=TRUE,
                    raster=FALSE) +  
  #BUG: if raster=TRUE, order=TRUE is ignored. So, set raster=FALSE
  scale_colour_gradientn(colours=rev(brewer.pal(n=11, name="RdBu"))[5:11])

ggsave("Yneg.tiff")

Seurat::FeaturePlot(object=subset(integrated_seurat, Ystatus == "Ypos"),
                    slot="data",
                    features="LOY_sig1",
                    #cols= c("grey", viridis(n=10, option="C", direction=1)),
                    pt.size=0.4,
                    order=TRUE,
                    min.cutoff='q10',
                    reduction="umap",
                    label=FALSE,
                    combine=TRUE,
                    raster=FALSE) +  
  #BUG: if raster=TRUE, order=TRUE is ignored. So, set raster=FALSE
  scale_colour_gradientn(colours=rev(brewer.pal(n=11, name="RdBu"))[5:11])

ggsave("Ypos.tiff")

#*************scRNASeq Simon

flop <- function(){
  proj <- "scRNASeq_Jinfen"
  source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")
  
  reduc <- "harmony"
  integrated_seurat <- base::readRDS(paste0(seurat_results, "integrated_seurat_snn.rds"))
  integrated_seurat1 <- subset(x = integrated_seurat,
                               subset = cluster.0.4.harmony %in% c(1,3,4,5,7,9,15))
  integrated_seurat1 <- subset(x = integrated_seurat,
                               subset = Sex == "Male" & cluster.0.4.harmony %in% c(1,3,4,5,7,9,15))
  y_genes <- c("Ddx3y", "Eif2s3y", "Kdm5d","Uty","Rbmy","Sly","Sry","Uba1y",
               "Usp9y","Zfy1","Zfy2","H2al2b","H2al2c","Orly","Rbm31y","Srsy",
               "Ssty1","Ssty2")
  # integrated_seurat <- base::readRDS(paste0(seurat_results, "Simon_integrated_seurat_snn.rds"))
  # integrated_seurat1 <- subset(x = integrated_seurat,
  #                             subset = Sex == "Male" & subtype %in% c("CD8T", "Naive_Tcell","Treg"))
  # y_genes <- read.xlsx(paste0(seurat_results,"Y Genes & Paralogs_207genes.xlsx"),
  #                      sheet="Sheet1")
  # y_genes <- y_genes$SYMBOL
  remove_y_genes <- c("PCDH11Y","TTTY14", "Pcdh11y" )
  y_genes <- setdiff(y_genes, remove_y_genes)
  features <- intersect(y_genes, rownames(integrated_seurat@assays$RNA$data))
  
  seurat_obj <- integrated_seurat1
  
  # Extract expression inf, keep only Y genes that have expression
  df <- seurat_obj@assays$RNA$data
  df <- df[features,]
  df <- df[rowSums(df) != 0,]
  
  # Convert count matrix to binary format [1=Expressed, 0=Not expressed]
  #df <- as.data.frame(df) %>% replace(.> 0, 1)
  df[df>0] <- 1
  
  yneg_cells <- colnames(df)[colSums(df) <= 0]  # 97% of female cells have colSums <2 fop Simon data
  ypos_cells <- colnames(df)[colSums(df) > 0]
  
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    dplyr::mutate(Y_status = dplyr::case_when(Cell %in% yneg_cells ~ "Yneg",
                                              Cell %in% ypos_cells ~ "Ypos"))
  # seurat_obj@meta.data <- seurat_obj@meta.data %>% 
  #   dplyr::mutate(Y_status = dplyr::case_when(barcode %in% yneg_cells ~ "Yneg",
  #                                             barcode %in% ypos_cells ~ "Ypos"))
  
  # Plot 2 genes
  Seurat::FeaturePlot(object=seurat_obj,
                      slot="data",
                      features="Gsk3b", #"GSK3B",
                      #cols= c("grey", viridis(n=10, option="C", direction=1)),
                      pt.size=0.4,
                      order=TRUE,
                      split.by="Y_status",
                      min.cutoff='q10',
                      reduction="umap",
                      label=FALSE,
                      combine=TRUE,raster=FALSE)
  #                    +  
  # #BUG: if raster=TRUE, order=TRUE is ignored. So, set raster=FALSE
  # scale_colour_gradientn(colours=rev(brewer.pal(n=11, name="RdBu"))[5:11])
  
  ggsave("GSK3B.tiff")
  
  Seurat::FeaturePlot(object=seurat_obj,
                      slot="data",
                      features="Pdcd1", #"PDCD1"
                      #cols= c("grey", viridis(n=10, option="C", direction=1)),
                      pt.size=0.4,
                      order=TRUE,
                      split.by="Y_status",
                      min.cutoff='q10',
                      reduction="umap",
                      label=FALSE,
                      combine=TRUE,raster=FALSE)
  #                    +  
  # #BUG: if raster=TRUE, order=TRUE is ignored. So, set raster=FALSE
  # scale_colour_gradientn(colours=rev(brewer.pal(n=11, name="RdBu"))[5:11])
  
  ggsave("PDCD1.tiff")
}

############ t score

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Prince/"

wb <- openxlsx::createWorkbook()
for (file in c("T14vsT0_SCR.DTKPA2.mageck.median.alphamedian.sgrna_summary.txt", "T14vsT0_SCR.DTKPB1.mageck.median.alphamedian.sgrna_summary.txt",
               "T14vsT0_YKO.DTKPA2.mageck.median.alphamedian.sgrna_summary.txt", "T14vsT0_YKO.DTKPB1.mageck.median.alphamedian.sgrna_summary.txt",
               "D14vsD0_2D.DTKP.mageck.median.alphamedian.sgrna_summary.txt", "D14vsD0_3D.DTKP.mageck.median.alphamedian.sgrna_summary.txt")){
  
  data <- read.table(paste0(data_path, file), header = TRUE) %>%
    dplyr::select(sgrna, Gene, control_mean, treat_mean, LFC)
  
  data <- calc_t_score(data)
  
  suffix <- gsub(pattern="mageck.median.alphamedian.sgrna_summary.txt", replacement="", x=file)
  openxlsx::addWorksheet(wb, sheetName = suffix)
  openxlsx::writeData(wb, sheet = suffix, x = data, rowNames = FALSE) 
}
openxlsx::saveWorkbook(wb, file = paste0(data_path, "t_score.xlsx"), overwrite = TRUE)

##### Volcano plot & Hit filtering

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Prince/"

# Find genes with negative efffect size in YKO and more negative effect size in YKO relative to Yscr
yscr_dtkpa2 <- read.xlsx(paste0(data_path, "t_score.xlsx"),
                         sheet="T14vsT0_SCR.DTKPA2.")
yscr_dtkpb1 <- read.xlsx(paste0(data_path, "t_score.xlsx"),
                         sheet="T14vsT0_SCR.DTKPB1.")
yko_dtkpa2 <- read.xlsx(paste0(data_path, "t_score.xlsx"),
                        sheet="T14vsT0_YKO.DTKPA2.")
yko_dtkpb1 <- read.xlsx(paste0(data_path, "t_score.xlsx"),
                        sheet="T14vsT0_YKO.DTKPB1.")

disp_genes <- c("Gsk3b", "Dcaf11", "Gngt1", "Pkn1", "Grik3", "Akr1b3", "Gad2", 
                "Mmp10", "Ptprg", "Uhmk1")
disp_genes <- c()

yscr <- dplyr::bind_rows(yscr_dtkpa2, yscr_dtkpb1)
plot_t_score(yscr, data_path, disp_genes, "yscr")
colnames(yscr) <- paste0(colnames(yscr), "_YSCR")
colnames(yscr)[1] <- "Gene"

yko <- dplyr::bind_rows(yko_dtkpa2, yko_dtkpb1)
plot_t_score(yko, data_path, disp_genes, "yko")
colnames(yko) <- paste0(colnames(yko), "_YKO")
colnames(yko)[1] <- "Gene"

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "hits_yscr")
openxlsx::writeData(wb, sheet = "hits_yscr", x = yscr, rowNames = FALSE) 
openxlsx::addWorksheet(wb, sheetName = "hits_yko")
openxlsx::writeData(wb, sheet = "hits_yko", x = yko, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(data_path, "Hits.xlsx"), overwrite = TRUE)

###### S curve

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Prince/"

# Read the t scores
dtkpa2_scr <- read.xlsx(paste0(data_path, "t_score.xlsx"),
                        sheet = "T14vsT0_SCR.DTKPA2.")
dtkpb1_scr <- read.xlsx(paste0(data_path, "t_score.xlsx"),
                        sheet = "T14vsT0_SCR.DTKPB1.")
tscore_scr <- dplyr::bind_rows(dtkpa2_scr, dtkpb1_scr)

dtkpa2_yko <- read.xlsx(paste0(data_path, "t_score.xlsx"),
                        sheet = "T14vsT0_YKO.DTKPA2.")
dtkpb1_yko <- read.xlsx(paste0(data_path, "t_score.xlsx"),
                        sheet = "T14vsT0_YKO.DTKPB1.")
tscore_yko <- dplyr::bind_rows(dtkpa2_yko, dtkpb1_yko)

# Read the logFC and pval
dtkpa2_scr <- read.table(paste0(data_path, "T14vsT0_SCR.DTKPA2.mageck.median.alphamedian.gene_summary.txt"), header = TRUE)
dtkpb1_scr <- read.table(paste0(data_path, "T14vsT0_SCR.DTKPB1.mageck.median.alphamedian.gene_summary.txt"), header = TRUE)
lfc_scr <- dplyr::bind_rows(dtkpa2_scr, dtkpb1_scr)

dtkpa2_yko <- read.table(paste0(data_path, "T14vsT0_YKO.DTKPA2.mageck.median.alphamedian.gene_summary.txt"), header = TRUE)
dtkpb1_yko <- read.table(paste0(data_path, "T14vsT0_YKO.DTKPB1.mageck.median.alphamedian.gene_summary.txt"), header = TRUE)
lfc_yko <- dplyr::bind_rows(dtkpa2_yko, dtkpb1_yko)

rm(dtkpa2_scr, dtkpb1_scr, dtkpa2_yko, dtkpb1_yko)

# Merge tscore and lfc
scr <- tscore_scr %>% dplyr::full_join(lfc_scr, by=c("Gene"="id"))
yko <- tscore_yko %>% dplyr::full_join(lfc_yko, by=c("Gene"="id"))

disp_genes <- c("Gsk3b", "Dcaf11", "Gngt1", "Pkn1", "Grik3", "Akr1b3", "Gad2", 
                "Mmp10", "Ptprg", "Uhmk1")

for (f in c("scr", "yko")){
  
  data <- get(f)
  
  for (c in c("U_gene", "t_score", "neg.lfc")){
    
    plot_data <- data %>% 
      dplyr::distinct_at("Gene", .keep_all = TRUE)
    
    # Set limits of x and y axis
    high <- max(plot_data %>% dplyr::select(all_of(c)) %>% unlist(use.names = FALSE), na.rm=TRUE)
    low <- min(plot_data %>% dplyr::select(all_of(c)) %>% unlist(use.names = FALSE), na.rm=TRUE)
    ylims <- c(-max(abs(low), abs(high)), max(abs(low), abs(high)))
    xlims <- c(-10, nrow(plot_data)+10)
    
    ggplot(data = plot_data, aes(x=reorder(x=Gene, X=get(c)),y=get(c))) +  #size = abs(get(c))
      geom_point() +
      geom_point(data = plot_data %>% dplyr::filter(Gene %in% disp_genes), color="red") +
      theme_classic() +
      ggplot2::labs(x = "Genes", y = c) +
      my_theme +
      coord_cartesian(xlim = xlims, ylim = ylims, clip = "off") +
      ggplot2::theme(axis.text.x= element_blank(),
                     axis.ticks.x = element_blank()) +
      geom_text_repel(data = plot_data %>% dplyr::filter(Gene %in% disp_genes),
                      mapping = aes(label = Gene),
                      size = 5,
                      force = 0.5,
                      point.size = 1,
                      angle = 0,
                      #vjust = 0,
                      #hjust = 0,
                      #direction = "y",
                      box.padding = 1,  # increases line length somehow
                      point.padding = 0.1,
                      max.overlaps = Inf,
                      xlim = c(NA, NA),
                      ylim = c(-Inf,NA),
                      min.segment.length = 0.2,
                      #min.segment.length = dplyr::if_else(draw_line == "TRUE", 0.2, Inf),
                      #arrow = arrow(length = unit(0.015, "npc")),
                      position = position_quasirandom())
    
    ggsave(paste0(data_path, f, "_", c, ".tiff"),
           width = 10,
           height = 8)
  }
}


#### Depmap of Yhigh  Ylow and gene effect score of 419 hits
#### TCGA BLCA of Yhigh Ylow  and heatmap of 419 hits

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Prince/"

# Read the expression data for DepMap cell lines & classify into Yhigh, Ylow
exp_df <- read.table(paste0(data_path, "DepMap.OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv"), 
                     header=TRUE, sep=",", quote="", skip=0, fill=TRUE)
colnames(exp_df) <- gsub(pattern = "\\.[0-9]*", replacement = "", colnames(exp_df))
exp_df <- exp_df %>%
  tibble::column_to_rownames("X") %>%
  t()

#...Run this for classification of YhighYlow
y_genes <- c("DDX3Y", "UTY", "KDM5D", "USP9Y", "ZFY", "RPS4Y1", "TMSB4Y", "EIF1AY", "NLGN4Y")

# Generate a list of gene sets
gs <- list(y_genes)
names(gs) <- "Y.sig"

# Calculate GSVA scores
gsvaPar <- GSVA::gsvaParam(exprData = as.matrix(exp_df), 
                           geneSets = as.list(gs))
gsva.scores <- gsva(gsvaPar, 
                    verbose=TRUE)

gsva.scores <- t(gsva.scores) %>% 
  data.frame() %>%
  dplyr::mutate(gsva.class = dplyr::case_when(Y.sig > 0 ~ "YHigh",
                                              TRUE ~ "Ylow")) %>%
  tibble::rownames_to_column("ModelID") %>%
  dplyr::rename(gsva.score = Y.sig)

# Calculate SSGSEA scores
ssgseaPar <- GSVA::ssgseaParam(exprData = as.matrix(exp_df),
                               geneSets = as.list(gs))
ssgsea.scores <- GSVA::gsva(ssgseaPar, 
                            verbose=TRUE)

ssgsea.scores <- t(ssgsea.scores) %>% 
  data.frame() %>%
  dplyr::mutate(ssgsea.class = dplyr::case_when(Y.sig > 0 ~ "YHigh",
                                                TRUE ~ "Ylow")) %>%
  tibble::rownames_to_column("ModelID") %>%
  dplyr::rename(ssgsea.score= Y.sig)
#...


# Read the human orthologs of the 419 hits
# yscr_genes <- read.xlsx(paste0(data_path, "Hits.xlsx"), sheet= "hits_yscr") %>%
#   dplyr::filter(U_gene_YSCR > 0)
# yko_genes <- read.xlsx(paste0(data_path, "Hits.xlsx"), sheet= "hits_yko") %>%
#   dplyr::filter(U_gene_YKO < 0)
# hits <- intersect(yko_genes$Gene,yscr_genes$Gene)
hits <- read.xlsx(paste0(data_path, "Human_orthologs.xlsx")) %>%
  dplyr::mutate(Human = make.names(Human)) %>%
  dplyr::select(Human) %>%
  unlist(use.names = FALSE) %>%
  unique() 

# Read the gene effect scores
main_df <- read.table(paste0(data_path, "DepMap.CRISPRGeneEffect.csv"), 
                      header=TRUE, sep=",", quote="", skip=0, fill=TRUE)
colnames(main_df) <- gsub(pattern = "\\.[0-9]*", replacement = "", colnames(main_df))
# Keep only gene effect score of hit genes
main_df <- main_df %>% 
  dplyr::select(X, intersect(hits, colnames(main_df)))

# Read the clinical data
metadata_df <- read.xlsx(paste0(data_path, "DepMap.Model.xlsx")) %>%
  dplyr::select(ModelID, StrippedCellLineName, DepmapModelType, Sex)

meta_data <- openxlsx::read.xlsx(xlsxFile = paste0(data_path, "TCGA.BLCA.Metadata.xlsx"))

# Merge Y status and clinical data with gene effect score
main_df <- main_df %>%
  dplyr::left_join(gsva.scores, by=c("X"="ModelID")) %>%
  dplyr::left_join(ssgsea.scores, by=c("X"="ModelID")) %>%
  dplyr::left_join(metadata_df, by=c("X"="ModelID")) %>%
  dplyr::rename(ModelID = X, CellLine = StrippedCellLineName, ProjectID = DepmapModelType) %>%
  dplyr::select(ModelID, CellLine, ProjectID, Sex, gsva.class, gsva.score, ssgsea.class, ssgsea.score, everything())

# Save the results as xlsx file
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName="Gene Effect Score")
openxlsx::writeData(wb, sheet="Gene Effect Score", x=main_df)
openxlsx::saveWorkbook(wb, file=paste0(data_path, "Gene_effect_scores.xlsx"),
                       overwrite=TRUE)

#### Heatmap of 419 hits in TCGA

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Prince/"

# Read the expression data for TCGA BLCA patients & classify into Yhigh, Ylow
exp_df <- openxlsx::read.xlsx(xlsxFile = paste0(data_path, "TCGA.BLCA.Normalized.counts.xlsx"))
colnames(exp_df)[1] <- "SYMBOL"
exp_df <- exp_df %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL, unique=TRUE)) %>%
  tibble::column_to_rownames("SYMBOL")

#...Run this for classification of YhighYlow
y_genes <- c("DDX3Y", "UTY", "KDM5D", "USP9Y", "ZFY", "RPS4Y1", "TMSB4Y", "EIF1AY", "NLGN4Y")

# Generate a list of gene sets
gs <- list(y_genes)
names(gs) <- "Y.sig"

# Calculate GSVA scores
gsvaPar <- GSVA::gsvaParam(exprData = as.matrix(exp_df), 
                           geneSets = as.list(gs))
gsva.scores <- gsva(gsvaPar, 
                    verbose=TRUE)

gsva.scores <- t(gsva.scores) %>% 
  data.frame() %>%
  dplyr::mutate(gsva.class = dplyr::case_when(Y.sig > 0 ~ "YHigh",
                                              TRUE ~ "Ylow")) %>%
  tibble::rownames_to_column("ModelID") %>%
  dplyr::rename(gsva.score = Y.sig)

# Calculate SSGSEA scores
ssgseaPar <- GSVA::ssgseaParam(exprData = as.matrix(exp_df),
                               geneSets = as.list(gs))
ssgsea.scores <- GSVA::gsva(ssgseaPar, 
                            verbose=TRUE)

ssgsea.scores <- t(ssgsea.scores) %>% 
  data.frame() %>%
  dplyr::mutate(ssgsea.class = dplyr::case_when(Y.sig > 0 ~ "YHigh",
                                                TRUE ~ "Ylow")) %>%
  tibble::rownames_to_column("ModelID") %>%
  dplyr::rename(ssgsea.score= Y.sig)

# Calculate z scores
normalized_counts <- exp_df
normalized_counts <- normalized_counts[!rowSums(normalized_counts, na.rm=TRUE) == 0,]
normalized_counts <- log(1+normalized_counts, base=2)
t <- base::apply(X=normalized_counts, MARGIN=1, FUN=median, na.rm=TRUE)
normalized_counts <- base::sweep(x=normalized_counts, MARGIN=1, FUN="-", STATS=t)
z.scores <- as.data.frame(advanced_Z(y_genes, normalized_counts)) %>%
  data.frame() %>%
  dplyr::rename(z.score= identity(1)) %>%
  dplyr::mutate(z.class = dplyr::case_when(z.score > 0 ~ "YHigh",
                                           TRUE ~ "Ylow")) %>%
  tibble::rownames_to_column("ModelID")

#...

# Read the human orthologs of the 419 hits
# yscr_genes <- read.xlsx(paste0(data_path, "Hits.xlsx"), sheet= "hits_yscr") %>%
#   dplyr::filter(U_gene_YSCR > 0)
# yko_genes <- read.xlsx(paste0(data_path, "Hits.xlsx"), sheet= "hits_yko") %>%
#   dplyr::filter(U_gene_YKO < 0)
# hits <- intersect(yko_genes$Gene,yscr_genes$Gene)
hits <- read.xlsx(paste0(data_path, "Human_orthologs.xlsx"), sheet="New") %>%
  dplyr::mutate(Human = make.names(Human)) %>%
  dplyr::select(Human) %>%
  unlist(use.names = FALSE) %>%
  unique() 

# Read the clinical data
#metadata_df <- openxlsx::read.xlsx(xlsxFile = paste0(data_path, "TCGA.BLCA.Metadata.xlsx")) 
metadata_df <- openxlsx::read.xlsx(xlsxFile = paste0(data_path, "TCGA_Metadata_Yiling.xlsx")) %>%
  dplyr::filter(CancerType == "BLCA") %>%
  dplyr::mutate(Sample_ID = Sample)

# Read the expr data
main_df <- openxlsx::read.xlsx(xlsxFile = paste0(data_path, "TCGA.BLCA.Normalized.counts.xlsx"))
colnames(main_df)[1] <- "SYMBOL"
main_df <-  main_df %>% 
  dplyr::filter(SYMBOL %in% hits) %>%
  dplyr::select(SYMBOL, intersect(metadata_df$Sample_ID, colnames(main_df)))

# Merge Y status and clinical data with gene effect score
# metadata_df <- metadata_df %>%
#   dplyr::left_join(gsva.scores, by=c("Sample_ID"="ModelID")) %>%
#   dplyr::left_join(ssgsea.scores, by=c("Sample_ID"="ModelID")) %>%
#   dplyr::left_join(z.scores, by=c("Sample_ID"="ModelID")) %>%
#   dplyr::select(Sample_ID, Sex, gsva.class, gsva.score, ssgsea.class, ssgsea.score, z.class, z.score, everything())

# Heatmap parameters
plot_genes <- main_df$SYMBOL
disp_genes <- plot_genes
file_suffix <- ""
output_path <- data_path
heatmap_params <- list(anno.row       = NULL,        # NULL, c("Group")
                       anno.column    = c("Y_status"), #gsva.class"),
                       row.split      = NA,     
                       col.split      = c("Y_status"),
                       row.cluster    = c("all"),           # c("alphabetical", "group", "all")
                       col.cluster    = c("group"),  # c("alphabetical", "group", "all")
                       log.transform  = TRUE,
                       scale          = TRUE,
                       border_color   = "white",
                       bar_width      = NA,              # NA , 5
                       bar_height     = NA,              # NA , 5
                       width          = NA, #5,              # NA
                       height         = NA, #5,              # NA 
                       matrix_color   = "rdbu",          # c("vrds", "rdbu")
                       expr_legend    = TRUE,  # set FALSE if it overlaps with annotation legends
                       file_format    = "tiff")

metadata_column <- metadata_df %>% dplyr::filter(Sample_ID %in% colnames(main_df))
metadata_row <- data.frame(SYMBOL = "")
norm_counts <- main_df %>% dplyr::select(SYMBOL, all_of(metadata_column$Sample_ID))
plot_heatmap(norm_counts, metadata_column, metadata_row, heatmap_params,
             plot_genes, disp_genes, file_suffix, output_path)

##### KM plots of 419 hits 

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Prince/"

# Classify TCGA BLCA into Yhigh and Ylow
# Human Y gene signature
plot_genes <- c("DDX3Y", "EIF1AY", "HSFY2", "KDM5D", "UTY", "NLGN4Y", 
                "PCDH11Y", "RPS4Y1", "TBL1Y", "TMSB4Y", "USP9Y", "ZFY", 
                "DAZ1", "DAZ2", "DAZ3", "DAZ4", "PRY2", "RBMY1A1")

# Import read_data
read_data <- openxlsx::read.xlsx(xlsxFile = paste0(data_path, "TCGA.BLCA.Normalized.counts.xlsx"))
#read_data <- read_data %>% dplyr::select(SYMBOL, starts_with("TCGA"))
colnames(read_data)[1] <- "SYMBOL"

# Import meta_data and subset meta_data if needed
meta_data <- openxlsx::read.xlsx(xlsxFile = paste0(data_path, "TCGA.BLCA.Metadata.xlsx")) 

# Reformat metadata 
meta_data <- meta_data %>% 
  dplyr::mutate(Sample_ID = make.names(Sample_ID)) %>%
  dplyr::mutate(Time = as.numeric(Time)) %>%
  dplyr::filter(Time > 0 & !is.na(Time)) %>%
  dplyr::distinct_at("Sample_ID", .keep_all = TRUE)

meta_data <- meta_data %>%
  dplyr::filter(Sex == "Male")

# Reformat read data
normalized_counts <- read_data %>%
  dplyr::mutate(SYMBOL = make.names(names = SYMBOL, unique = TRUE)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  tibble::column_to_rownames(var = "SYMBOL") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
colnames(normalized_counts) <- base::make.names(names = colnames(normalized_counts))
normalized_counts <- normalized_counts[,intersect(make.names(meta_data$Sample_ID), colnames(normalized_counts))]
normalized_counts <- normalized_counts[!rowSums(normalized_counts, na.rm=TRUE) == 0,]

normalized_counts <- log(1+normalized_counts, base=2)
t <- base::apply(X=normalized_counts, MARGIN=1, FUN=median, na.rm=TRUE)
normalized_counts <- base::sweep(x=normalized_counts, MARGIN=1, FUN="-", STATS=t)

# Calculate score for Y genes
expr_df <- as.data.frame(advanced_Z(plot_genes, normalized_counts))

# Merge expression data with survival data
expr_df <- expr_df %>%
  data.frame() %>%
  dplyr::rename(combined.exp = identity(1)) %>%
  tibble::rownames_to_column("Sample_ID") %>%
  dplyr::inner_join(meta_data, by=c("Sample_ID"="Sample_ID")) %>%
  dplyr::select(Sample_ID, combined.exp, Time, Status)

# Standard parameters for survival
subset_group <- NA
subset_value <- NA
plot_by <- "Expression"
split_by <- subset_group
combine_plot <- FALSE 
multiple_cutoff <- TRUE
stratify_criteria <- "o"
reference <- "T1"
reference <- dplyr::if_else(plot_by == "Expression", "LOW", reference)
confidence_interval <- FALSE
legend_title <- "Expression"
plot_risk_table <- TRUE
plot_curve <- TRUE
all_quartiles <- FALSE
gene_signature <- TRUE
color_palette <- c("#d73027","#0c2c84")
color_palette <- c(color_palette, 
                   colorspace::adjust_transparency(col = color_palette, alpha = 0.3), 
                   colorspace::adjust_transparency(col = color_palette, alpha = 0.6))
variable_x <- "Time"
variable_y <- "Status" 

# Case specific parameters
gene <- "combined.exp"
prefix <- ""

summary <- wrangle_data(expr_df, gene, stratify_criteria, prefix, data_path)

surv_df <- summary[[1]] %>%
  dplyr::mutate(Ystatus = dplyr::case_when(model == "HIGH" ~ "Yhigh",
                                           TRUE ~ "Ylow")) %>%
  dplyr::select(Sample_ID, Ystatus, everything(), -model)

# Plot individual KMplots for 419 hits
yscr_genes <- read.xlsx(paste0(data_path, "Hits.xlsx"), sheet= "hits_yscr") %>%
  dplyr::filter(U_gene_YSCR > 0)
yko_genes <- read.xlsx(paste0(data_path, "Hits.xlsx"), sheet= "hits_yko") %>%
  dplyr::filter(U_gene_YKO < 0)
plot_genes <- intersect(yko_genes$Gene,yscr_genes$Gene)

# Save as excel and get human homologs
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "mouse hits")
openxlsx::writeData(wb, sheet = "mouse hits", x = data.frame(plot_genes))
openxlsx::saveWorkbook(wb, file = paste0(data_path, "Get_human_homologs.xlsx"), overwrite = TRUE)

# Import human orthologs
plot_genes <- read.xlsx(paste0(data_path, "Human_orthologs.xlsx")) %>%
  dplyr::mutate(Human = make.names(Human)) %>%
  dplyr::select(Human) %>%
  unlist(use.names = FALSE) %>%
  unique()

# Keep only genes that have expression data
plot_genes <- intersect(plot_genes, rownames(normalized_counts))

for (status in c("Yhigh", "Ylow")){
  
  # Create a list to store survminer cutoffs, coxph stats, etc..
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
                "Tarone_Ware.early"  = c(),
                "Peto_Peto.early" = c(),
                "modified_Peto_Peto"  = c(),
                "Fleming_Harrington" = c())
  
  # Subset Yhigh or Ylow patients
  meta_data <- surv_df %>%
    dplyr::filter(Ystatus == status)
  
  # Merge expression data with survival data
  expr_df <- normalized_counts %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::inner_join(meta_data, by=c("Sample_ID"="Sample_ID")) %>%
    dplyr::select(Sample_ID, plot_genes, Time, Status)
  
  classification_df <- expr_df %>% 
    dplyr::select(Sample_ID) %>%
    dplyr::mutate(Dummy_col = 0)
  
  # Plot survival curves
  for (gene in plot_genes) {
    plot_curve <- FALSE
    prefix <- status
    summary <- wrangle_data(expr_df, gene, stratify_criteria, prefix, data_path)
    
    class_df <- summary[[1]] %>%
      dplyr::select(Sample_ID, all_of(gene), model) %>%
      dplyr::rename(!!paste0(gene, "_model") := model)
    
    classification_df <- classification_df %>% 
      dplyr::left_join(class_df, by=("Sample_ID"="Sample_ID"))
    
    stats$gene          <- c(stats$gene, summary[[2]]$gene)
    stats$group         <- c(stats$group, summary[[2]]$group)
    stats$lower_cutoff  <- c(stats$lower_cutoff, summary[[2]]$lower)
    stats$middle_cutoff <- c(stats$middle_cutoff, summary[[2]]$middle)
    stats$upper_cutoff  <- c(stats$upper_cutoff, summary[[2]]$upper)
    stats$HR            <- c(stats$HR, summary[[2]]$HR )
    stats$CI_lower      <- c(stats$CI_lower, summary[[2]]$CI_lower)
    stats$CI_upper      <- c(stats$CI_upper, summary[[2]]$CI_upper)
    stats$pvalue        <- c(stats$pvalue, summary[[2]]$pvalue)
    stats$logrank       <- c(stats$logrank, summary[[2]]$logrank) 
    stats$reg_logrank.late    <- c(stats$reg_logrank.late, summary[[2]]$reg_logrank.late)
    stats$Gehan_Breslow.early <- c(stats$Gehan_Breslow.early, summary[[2]]$Gehan_Breslow.early)
    stats$Tarone_Ware.early   <- c(stats$Tarone_Ware.early, summary[[2]]$Tarone_Ware.early)
    stats$Peto_Peto.early     <- c(stats$Peto_Peto.early, summary[[2]]$Peto_Peto.early)
    stats$modified_Peto_Peto  <- c(stats$modified_Peto_Peto, summary[[2]]$modified_Peto_Peto)
    stats$Fleming_Harrington  <- c(stats$Fleming_Harrington, summary[[2]]$Fleming_Harrington)
  }
  
  stats_df <- data.frame(stats)
  val_df <- normalized_counts[intersect(plot_genes, rownames(normalized_counts)),] %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column("SYMBOL")
  
  # Save the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Summary")
  openxlsx::writeData(wb, sheet = "Summary", x = stats_df)
  openxlsx::addWorksheet(wb, sheetName = "Classification")
  openxlsx::writeData(wb, sheet = "Classification", x = classification_df)
  openxlsx::addWorksheet(wb, sheetName = "Norm_counts")
  openxlsx::writeData(wb, sheet = "Norm_counts", x = val_df)
  openxlsx::saveWorkbook(wb, file = paste0(data_path, "Individual_stats_", status, ".xlsx"), overwrite = TRUE)
}

#******************ORA analysis

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Prince/"

# Read the gene set files
species <- "Mouse"
gmt_dir <- paste0("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets/", species, "/")
gmt_files <- list.files(gmt_dir, full.names = TRUE)

# Define universe genes
universe_genes <- read.xlsx(paste0(data_path, "Hits.xlsx")) %>% 
  dplyr::select(Gene) %>%
  unlist(use.names=FALSE)

input_genes_258 <- read.xlsx(paste0(data_path, "Prince_Hits.xlsx")) %>% dplyr::select(Genes_258) %>% unlist(use.names = FALSE) %>% unique()
input_genes_421 <- read.xlsx(paste0(data_path, "Prince_Hits.xlsx")) %>% dplyr::select(Genes_421) %>% unlist(use.names = FALSE) %>% unique()

ora_result_258 <- ora(input_genes_258, universe_genes, gmt_files)
ora_result_421 <- ora(input_genes_258, universe_genes, gmt_files)

# Create workbook to store results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName="ora_result_258")
openxlsx::writeData(wb, sheet="ora_result_258", x=ora_result_258, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName="ora_result_421")
openxlsx::writeData(wb, sheet="ora_result_421", x=ora_result_421, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(data_path, "Pathway_Analysis_Results.xlsx"),
                       overwrite = TRUE)


for (p in c("BIOCARTA", "REACTOME", "WP", "GOBP", "GOCC", "GOMF", "HALLMARK")){
  
  df1 <- ora_result_258 %>% dplyr::filter(grepl(pattern=p, x=Description)) %>% dplyr::slice_max(order_by=k.K, n=25)
  df2 <- ora_result_421 %>% dplyr::filter(grepl(pattern=p, x=Description)) %>% dplyr::slice_max(order_by=k.K, n=25)
  
  plot_ora(df1, paste0("ora_258_",p), data_path)
  plot_ora(df2, paste0("ora_421_",p), data_path)
  
}


#######Heatmap pancancer single cell

parent.dir  <- "/hpc/home/kailasamms/scratch"
gmt.dir     <- "/hpc/home/kailasamms/projects/GSEA_genesets"
scripts.dir <- "/hpc/home/kailasamms/projects/scRNASeq"

# Run the Custom_Functions.R script
path1 <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/Custom_Functions.R"
path2 <- "/hpc/home/kailasamms/projects/scRNASeq/Custom_Functions.R"
if (file.exists(path1)) {
  source(path1)
} else if (file.exists(path2)) {
  source(path2)
}

proj <- "scRNASeq_PanCancer"
species <- "Homo sapiens"
contrasts <- c("Treatment-Reference")

# DESeq2 overrides
deseq2.override <- list(
  contrasts     = contrasts,
  #design        = "Comparisons",            # DESeq2 design formula or column name
  #lfc.cutoff    = 0,                        # Log fold change cutoff for significance
  #padj.cutoff   = 0.1,                      # Adjusted p-value cutoff for significance
  batch.correct = FALSE                     # Boolean, whether to apply batch correction
)

# Heatmap overrides
heatmap.override <- list(
  #force.log        = TRUE,                  # Force log transformation
  col.ann          = c("Cancertype_TCGA", "LOY", "Age"),                  # Column annotation
  #row.ann          = NULL,                  # Row annotation
  #col.gaps         = NULL,                  # Column gaps
  #row.gaps         = NULL,                  # Row gaps
  col.cluster      = c("LOY"),                 # Column clustering
  #row.cluster      = "all",                 # Row clustering
  #palette         = "rdbu",                # Heatmap palette
  #ann.palette     = "discrete",            # Annotation palette
  #border.color    = NA,                    # Cell border color
  #show.expr.legend = TRUE,                  # Show expression legend
  #title           = "",                    # Heatmap title
  format           = "tiff"                 # Output file format
)

# Volcano plot overrides
volcano.override <- list(
  #lfc.cutoff  = 0.58,                         # Minimum log2 fold-change to highlight
  #padj.cutoff = 0.05,                      # Adjusted p-value cutoff
  #color       = "vrds",                    # Color palette
  #label.genes = c()                         # Genes to label on the plot
)

proj.params <- setup_project(
  proj               = proj,
  species = species,
  contrasts = contrasts,
  parent.dir = parent.dir,
  gmt.dir = gmt.dir,
  scripts.dir      = scripts.dir,
  deseq2.override  = deseq2.override,
  heatmap.override = heatmap.override,
  volcano.override = volcano.override
)


save.dir <- "/hpc/home/kailasamms/scratch/scRNASeq_PanCancer/"

samples <- integ.ann@meta.data %>% 
  dplyr::filter(Sex == "Male") %>%
  dplyr::pull(SampleID) %>%
  unique()

# assay <- "RNA"
# read_data <- data.frame(matrix(NA, nrow=nrow(degs_seurat@assays[[assay]]@data), ncol=length(samples)))
# rownames(read_data) <- rownames(degs_seurat@assays[[assay]]@data)
# colnames(read_data) <- samples 
# 
# for(i in samples){
#   
#   # Create a list of cells for each sample
#   cells_subset <- degs_seurat@meta.data %>% dplyr::filter(SampleID == i) %>% rownames()
#   
#   # Use data.frame to convert "." in sparse matrix to "0"
#   subset <- data.frame(degs_seurat@assays[[assay]]@data[,cells_subset])
#   read_data[,i]  <- rowSums(subset)
#   
#   print(i)
# }
# read_data <- read_data[rowSums(read_data) != 0,]
plot_genes <- intersect(plot_genes, rownames(integ.ann@assays[["RNA"]]@data))
# norm_counts <- read_data[plot_genes, samples]
# metadata_col <- seurat_obj@meta.data %>% 
#   dplyr::filter(Sex == "Male") %>%
#   dplyr::rename(Sample.ID = SampleID) %>%
#   dplyr::distinct_at("Sample.ID", .keep_all = TRUE)
# ph <- plot_heatmap(norm_counts, proj.params, metadata_col)
# jpeg(file.path(save.dir, "Heatmap.jpeg"), width = 5, height = 7, units = "in", res = 300)
# gridExtra::grid.arrange(grobs = list(ph$ph$gtable), ncol = 1)
# dev.off()

integ.ann <- readRDS(file.path(proj.params$seurat_dir, "scRNAseq_pancan_normalized_after_scanvi.rds"))
integ.ann <- subset(integ.ann,
                    Sex == "Male")

integ.ann.coi <- subset(integ.ann, Cancertype_TCGA %in% coi)
epi_seurat <- subset(integ.ann, Celltype1 == "Epithelium")
epi_seurat_coi <- subset(integ.ann, Celltype1 == "Epithelium" & Cancertype_TCGA %in% coi)
epi_seurat_blca <- subset(integ.ann, Celltype1 == "Epithelium" & Cancertype_TCGA %in% "BLCA")
TNK_seurat_blca <- subset(integ.ann, Celltype1 == "T_NK_cell" & Cancertype_TCGA %in% "BLCA")

coi <- c("HNSC", "LIHC", "CHOL", "CRC", "LC", "STAD", "BLCA", "PDAC")
plot_genes <- c("AKR1B1", "CAMK1D", "CAMKK1", "CDKN1B", "CHEK2", "CHRNA7", 
                "CLK1", "CP", "CRYZ", "DBI", "DCAF11", "DGKZ", "DUSP19", 
                "EHHADH", "GAD2", "GNGT1", "GRIA1", "GRID1", "GRIK3", "GSK3B",
                "IYD", "KCND3", "LATS1", "MAP3K6", "MKNK2", "MMP10", "MTMR1", 
                "MTNR1B", "MTR", "PDP1", "PFKFB3", "PKN1", "PPM1A", "PPP1R15B",
                "PRKACB", "PRKD1", "PTPRG", "RARRES1", "RNASE4", "RPS6KA6", 
                "STK31", "TAOK1", "TESK2", "UHMK1")

# --- Dot plot 44 hits for 12 cancer type separately LOY vs WT ---
assay <- "RNA"
ident.1 <- "LOY"
ident.2 <- NULL
features <- plot_genes
output_path <- proj.params$seurat_dir
CellType.col <- "Celltype1"
split.col <- NULL
for (i in as.character(unique(epi_seurat@meta.data[["Cancertype_TCGA"]]))){
  subset_obj <- epi_seurat[, epi_seurat@meta.data[["Cancertype_TCGA"]] == i]
  filename <- paste0("Dot_plot_44hits_LOY_", i)
  plot_dot_custom(subset_obj, assay, proj.params, ident.1, ident.2, features, filename, output_path, CellType.col, split.col)
}

# --- Dot plot 44 hits for 12 cancer type separately by celltype ---
assay <- "RNA"
ident.1 <- "Celltype1"
ident.2 <- NULL
features <- plot_genes
output_path <- proj.params$seurat_dir
CellType.col <- "Celltype1"
split.col <- NULL
for (i in as.character(unique(integ.ann@meta.data[["Cancertype_TCGA"]]))){
  subset_obj <- integ.ann[, integ.ann@meta.data[["Cancertype_TCGA"]] == i]
  filename <- paste0("Dot_plot_44hits_celltype_", i)
  plot_dot_custom(subset_obj, assay, proj.params, ident.1, ident.2, features, filename, output_path, CellType.col, split.col)
}

#--- Dot plot 44 hits for 12/8 cancer types combined LOY vs WT ---
assay <- "RNA"
ident.1 <- "LOY"
ident.2 <- NULL
features <- plot_genes
output_path <- proj.params$seurat_dir
CellType.col <- "Celltype1"
split.col <- NULL

filename <- paste0("Dot_plot_44hits_LOY_Combined_8_cancers")
plot_dot_custom(epi_seurat_coi, assay, proj.params, ident.1, ident.2, features, filename, output_path, CellType.col, split.col)

filename <- paste0("Dot_plot_44hits_LOY_Combined_12_cancers")
plot_dot_custom(epi_seurat, assay, proj.params, ident.1, ident.2, features, filename, output_path, CellType.col, split.col)

# --- Dot plot GSK3B for each cancer type  LOY vs WT ---
assay <- "RNA"
ident.1 <- "LOY"
ident.2 <- "Cancertype_TCGA"
features <- "GSK3B"
output_path <- proj.params$seurat_dir
CellType.col <- "Celltype1"
split.col <- NULL

filename <- "Dot_plot_GSK3B_LOY_8_cancers"
plot_dot_custom(epi_seurat_coi, assay, proj.params, ident.1, ident.2, features, filename, output_path, CellType.col, split.col)

filename <- "Dot_plot_GSK3B_LOY_12_cancers"
plot_dot_custom(epi_seurat, assay, proj.params, ident.1, ident.2, features, filename, output_path, CellType.col, split.col)

# --- Feature Plots 44 hits for 12/8 cancer types combined LOY vs WT ---
reduction <- "umap"
output_path <- proj.params$seurat_dir
split.col <- "LOY"

features <- plot_genes
filename <- "Feature_plot_44hits_LOY_Combined_12_cancers"
features <- "GSK3B"
filename <- "Feature_plot_GSK3B_LOY_Combined_12_cancers"
plot_features(integ.ann, features, reduction, filename, output_path, split.col)

features <- plot_genes
filename <- "Feature_plot_44hits_LOY_Combined_8_cancers"
features <- "GSK3B"
filename <- "Feature_plot_GSK3B_LOY_Combined_8_cancers"
plot_features(integ.ann.coi, features, reduction, filename, output_path, split.col)

features <- "GSK3B"
filename <- "Feature_plot_GSK3B_LOY_Combined_8_cancers_Epi_only"
plot_features(epi_seurat_coi, features, reduction, filename, output_path, split.col)

features <- c("GSK3B", "PDCD1", "CD274")
filename <- "Feature_plot_GSK3B,PDCD1,CD274_LOY_BLCA_Epi_only"
plot_features(epi_seurat_blca, features, reduction, filename, output_path, split.col)

features <- c("GSK3B", "PDCD1", "CD274")
filename <- "Feature_plot_GSK3B,PDCD1,CD274_LOY_BLCA_TNK_only"
plot_features(TNK_seurat_blca, features, reduction, filename, output_path, split.col)

# --- UMAP cell types 8 cancer types ---
reduction <- "umap"
color.col <- "Celltype1"
filename <- "UMAP_CellType"
output_path <- proj.params$seurat_dir
plot_umap(integ.ann.coi, reduction, color.col, filename, output_path)

# --- UMAP LOY 8 cancer types ---
reduction <- "umap"
color.col <- "LOY"
filename <- "UMAP_LOY"
output_path <- proj.params$seurat_dir
plot_umap(integ.ann.coi, reduction, color.col, filename, output_path)

# --- UMAP TCAGCancerType 8 cancer types ---
reduction <- "umap"
color.col <- "Cancertype_TCGA"
filename <- "UMAP_CancerType"
output_path <- proj.params$seurat_dir
plot_umap(integ.ann.coi, reduction, color.col, filename, output_path)

# --- UMAP cell types 12 cancer types ---
reduction <- "umap"
color.col <- "Celltype1"
filename <- "UMAP_CellType_12_cancer"
output_path <- proj.params$seurat_dir
split.col<- "Cancertype_TCGA"
plot_umap(integ.ann.coi, reduction, color.col, filename, output_path, split.col)

# --- Survival curves
proj <- "TCGA_BLCA"
meta_data <- openxlsx::read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/OneDrive_1_9-19-2025/TCGA_BLCA_DNA_Metadata.xlsx")
read_data <- openxlsx::read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/OneDrive_1_9-19-2025/TCGA_BLCA_DNA_Normalized.xlsx")
output_path <-  "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/OneDrive_1_9-19-2025/WT/"
#output_path <-  "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/OneDrive_1_9-19-2025/LOY/"

survival_params <- list(plot_by             = c(NA), 
                        split_by            = c(NA),
                        split_plot          = FALSE,
                        multiple_cutoff     = FALSE,
                        stratify_criteria   = c("o"),
                        reference           = c("LOW"),
                        conf_interval       = FALSE,
                        plot_curve          = TRUE,
                        plot_risk_table     = TRUE,
                        legend_title        = "Expression",
                        legend_label        = c("High", "Low"),
                        color_palette       = c("#d73027","#0c2c84"),
                        plot_all_bins       = FALSE,
                        plot_all_quartiles  = FALSE,
                        gene_sig_score      = FALSE)

# Reformat metadata 
meta_data <- meta_data %>% 
  dplyr::mutate(Sample_ID = make.names(names = Sample_ID)) %>%
  dplyr::mutate(Time = as.numeric(Time)) %>%
  dplyr::filter(Time > 0 & !is.na(Time)) %>%
  dplyr::distinct_at("Sample_ID", .keep_all = TRUE) %>%
  dplyr::filter(Sex == "Male", Y_status == "WT")
#dplyr::filter(Sex == "Male", Y_status == "LOY")

# Reformat read data
norm_counts <- read_data %>%
  dplyr::mutate(SYMBOL = make.names(names = SYMBOL, unique = TRUE)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  tibble::column_to_rownames(var = "SYMBOL") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
colnames(norm_counts) <- base::make.names(names = colnames(norm_counts))
norm_counts <- norm_counts[,intersect(make.names(meta_data$Sample_ID), colnames(norm_counts))]
norm_counts <- norm_counts[!rowSums(norm_counts, na.rm=TRUE) == 0,]

#log_norm_counts <- log(1+norm_counts, base=2)
#t <- base::apply(X=log_norm_counts, MARGIN=1, FUN=median, na.rm=TRUE)
#log_norm_counts <- base::sweep(x=log_norm_counts, MARGIN=1, FUN="-", STATS=t)

plot_genes <- c("AKR1B1", "CAMK1D", "CAMKK1", "CDKN1B", "CHEK2", "CHRNA7", 
                "CLK1", "CP", "CRYZ", "DBI", "DCAF11", "DGKZ", "DUSP19", 
                "EHHADH", "GAD2", "GNGT1", "GRIA1", "GRID1", "GRIK3", "GSK3B",
                "IYD", "KCND3", "LATS1", "MAP3K6", "MKNK2", "MMP10", "MTMR1", 
                "MTNR1B", "MTR", "PDP1", "PFKFB3", "PKN1", "PPM1A", "PPP1R15B",
                "PRKACB", "PRKD1", "PTPRG", "RARRES1", "RNASE4", "RPS6KA6", 
                "STK31", "TAOK1", "TESK2", "UHMK1")

# Generate expr_df
expr_df <- prep_expr_df(norm_counts, meta_data, plot_genes, survival_params)

# Create a list to store survminer cutoffs, coxph stats, etc..
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
              "Tarone_Ware.early"  = c(),
              "Peto_Peto.early" = c(),
              "modified_Peto_Peto"  = c(),
              "Fleming_Harrington" = c())

# Create a dataframe to classification info
classification_df <- expr_df %>% 
  dplyr::select(Sample_ID) %>%
  dplyr::mutate(Dummy_col = 0)

# Plot survival curves
if (survival_params$gene_sig_score != TRUE){
  
  for (gene in plot_genes) {
    
    prefix <- paste0(proj, "_", gene)
    summary <- plot_survival(expr_df, gene, survival_params, prefix, output_path)
    
    # Get sample classification info for the individual gene
    class_df <- summary[[1]] %>%
      dplyr::select(Sample_ID, all_of(gene), model) %>%
      dplyr::rename(!!paste0(gene, "_model") := model)
    
    # Merge classification info to parent dataframe
    classification_df <- classification_df %>% 
      dplyr::left_join(class_df, by=("Sample_ID"="Sample_ID"))
    
    # Store stats for individual gene
    stats$gene          <- c(stats$gene, summary[[2]]$gene)
    stats$group         <- c(stats$group, summary[[2]]$group)
    stats$lower_cutoff  <- c(stats$lower_cutoff, summary[[2]]$lower)
    stats$middle_cutoff <- c(stats$middle_cutoff, summary[[2]]$middle)
    stats$upper_cutoff  <- c(stats$upper_cutoff, summary[[2]]$upper)
    stats$HR            <- c(stats$HR, summary[[2]]$HR )
    stats$CI_lower      <- c(stats$CI_lower, summary[[2]]$CI_lower)
    stats$CI_upper      <- c(stats$CI_upper, summary[[2]]$CI_upper)
    stats$pvalue        <- c(stats$pvalue, summary[[2]]$pvalue)
    stats$logrank       <- c(stats$logrank, summary[[2]]$logrank) 
    stats$reg_logrank.late    <- c(stats$reg_logrank.late, summary[[2]]$reg_logrank.late)
    stats$Gehan_Breslow.early <- c(stats$Gehan_Breslow.early, summary[[2]]$Gehan_Breslow.early)
    stats$Tarone_Ware.early   <- c(stats$Tarone_Ware.early, summary[[2]]$Tarone_Ware.early)
    stats$Peto_Peto.early     <- c(stats$Peto_Peto.early, summary[[2]]$Peto_Peto.early)
    stats$modified_Peto_Peto  <- c(stats$modified_Peto_Peto, summary[[2]]$modified_Peto_Peto)
    stats$Fleming_Harrington  <- c(stats$Fleming_Harrington, summary[[2]]$Fleming_Harrington)
  }
  
  # Merge all stats into a dataframe
  stats_df <- data.frame(stats)
  
  # Create a dataframe of normalized counts of genes plotted
  norm_df <- norm_counts[intersect(plot_genes, rownames(norm_counts)),] %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column("SYMBOL")
  
  # Save the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Summary")
  openxlsx::writeData(wb, sheet = "Summary", x = stats_df)
  openxlsx::addWorksheet(wb, sheetName = "Classification")
  openxlsx::writeData(wb, sheet = "Classification", x = classification_df)
  openxlsx::addWorksheet(wb, sheetName = "Norm_counts")
  openxlsx::writeData(wb, sheet = "Norm_counts", x = norm_df)
  openxlsx::saveWorkbook(wb, file = paste0(output_path, proj, "_Individual_stats.xlsx"), overwrite = TRUE)
}