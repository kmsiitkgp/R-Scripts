

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

source("C:/Users/kailasamms/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
  
# parent directory : directory where input files, results, etc are stored
parent_path <- "C:/Users/kailasamms/Desktop/Collaboration projects/Prince/"
results_path <- parent_path

#*****************Run GSEA on YKO vs scr RNASeq of MB49 cells******************#

# Import DEGs dataframe with "SYMBOL", "padj" and "log2FoldChange" columns
DEGs_df <- read.xlsx(paste0(parent_path, "Results_id_YKO_centromere_vs_scr_KO_centromere_DESeq2_modelled_DEGs.xlsx"))

# List all gmt files
gmt_files <- list.files("C:/Users/kailasamms/Documents/GitHub/R-Scripts/GSEA_genesets/Mouse", full.names = TRUE)

species <- "Mus musculus"
annotations <- get_annotation(species)

# Loop over each gmt file and perform GSEA analysis
for (f in gmt_files){
  fgsea(DEGs_df, f, annotations)
}

#******Calculate normalized counts for YKO vs scr Proteomics of MB49 cells*****#

# Import raw counts
raw_counts <- read.xlsx(paste0(parent_path, "Results_id_YKO_centromere_vs_scr_KO_centromere_pg_matrix.xlsx"))

# Remove rows that have multiple genes or missing genes
raw_counts <- raw_counts[!grepl(pattern=";", x=raw_counts$Genes),] %>%
  dplyr::filter(!is.na(Genes)) %>%
  dplyr::arrange(Genes)
rownames(raw_counts) <- seq(1:nrow(raw_counts))

# Create metadata
metadata <- data.frame("Sample" = colnames(raw_counts)[-1], 
                       "Condition" = c("SCR","SCR", "SCR", "YKO", "YKO", "YKO"))

# Notice that some genes have duplicates. This is because they correspond to 
# different peptides detected by the mass spec algorithm. We keep the row that
# contains the highest intensity for most samples and discard rest of duplicates
dup_genes <- raw_counts %>% 
  dplyr::count(Genes) %>% 
  dplyr::filter(n>1) %>%
  dplyr::select(Genes) %>%
  unlist(use.names = FALSE)

# Create empty dataframe to store max intensities of duplicated genes
df <- data.frame(matrix(NA, nrow = 1, ncol = ncol(raw_counts)))
colnames(df) <- colnames(raw_counts)

# Iterate through each duplicated gene and find which row has max intensities
# across multiple samples 
for (gene in dup_genes){
  
  test <- raw_counts %>% 
    dplyr::filter(Genes %in% gene) %>%
    base::replace(is.na(.), 0)
  
  row_id <- c()
  for (i in 2:ncol(test)){
    row_id <- c(row_id, which.max(test[,i]))
  }
  
  row_id <- as.data.frame(table(row_id)) %>% 
    dplyr::slice_max(Freq) %>%  # if there are 2 max values both are selected
    dplyr::select(row_id) %>%
    unlist(use.names=FALSE) %>% 
    as.character() %>%
    as.numeric()
  
  df1 <- test[row_id,]
  if (all(colnames(df) == colnames(df1))){
    df <- dplyr::bind_rows(df,df1)
  }
}

# Remove the dummy column with NAs
df <- df[-1,]

# If some genes still have 2 rows, average them
df <- df %>% 
  dplyr::group_by(Genes) %>% 
  summarize(across(.cols=everything(), .fns=mean)) %>%
  dplyr::ungroup()

# Remove duplicated genes from raw_counts    
raw_counts <- raw_counts %>% 
  dplyr::filter(!Genes %in% dup_genes)

# Add max values we calcualted to raw_counts
if (all(colnames(df) == colnames(raw_counts))){
  raw_counts <- dplyr::bind_rows(df, raw_counts)
} 

# Replace all 0 with NA again
raw_counts <- raw_counts %>% 
  dplyr::mutate(across(.cols=where(is.numeric), ~na_if(., 0))) %>%
  tibble::column_to_rownames("Genes")
 
# Impute missing values in biological replicates
imputed_counts <- impute_with_mean(raw_counts)

# Perform quantile normalization
quant_norm <- TRUE
quant_norm_counts <- quantile_norm(imputed_counts, quant_norm)

# Calculate padj and log2FoldChange
Target <- "YKO" 
Reference <- "SCR" 
results <- calc_stats(quant_norm_counts, metadata, Target, Reference)

# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "results")
openxlsx::writeData(wb, sheet = "results", x = results, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName = "raw_intensites")
openxlsx::writeData(wb, sheet = "raw_intensites", x = raw_counts, rowNames=TRUE)
openxlsx::addWorksheet(wb, sheetName = "imputed_intensities")
openxlsx::writeData(wb, sheet = "imputed_intensities", x = imputed_counts, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName = "quant_norm_intensities")
openxlsx::writeData(wb, sheet = "quant_norm_intensities", x = quant_norm_counts, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(results_path, "pg_matrix.xlsx"),
                       overwrite = TRUE)

#***************Run GSEA on YKO vs scr Proteomics of MB49 cells****************#

# Import DEGs dataframe with "SYMBOL", "padj" and "log2FoldChange" columns
DEGs_df <- read.xlsx(paste0(parent_path, "pg_matrix.xlsx"))

# List all gmt files
gmt_files <- list.files("C:/Users/kailasamms/Documents/GitHub/R-Scripts/GSEA_genesets/Mouse", full.names = TRUE)

species <- "Mus musculus"
annotations <- get_annotations(species)

# Loop over each gmt file and perform GSEA analysis
for (f in gmt_files){
  fgsea(DEGs_df, f, annotations)
}