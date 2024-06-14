library(openxlsx)

annotations <- get_annotations(species)

x_GOI <- annotations %>% 
  dplyr::filter(CHR=="X" ) %>%
  dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) ~ ENSEMBL_ID, !is.na(SYMBOL) ~ SYMBOL)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL))

y_GOI <- annotations %>% 
  dplyr::filter(CHR=="Y") %>%
  dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) ~ ENSEMBL_ID, !is.na(SYMBOL) ~ SYMBOL)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL))

GOI <- annotations %>% 
  #dplyr::filter(CHR=="X" | CHR == "Y") %>%
  dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) ~ ENSEMBL_ID, !is.na(SYMBOL) ~ SYMBOL)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL))

#TCGA
tcga <- read.xlsx("/hpc/home/kailasamms/scratch/biomarker/TCGA_Normalized_Counts.xlsx") %>%
  dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) ~ ENSEMBL_ID, !is.na(SYMBOL) ~ SYMBOL)) %>%
  dplyr::select(everything(), -c("ENSEMBL_ID", "ENTREZ_ID", "SYMBOL_ENTREZ")) %>%
  # If there are duplicated genes, keep only data for highest expressing copy
  dplyr::mutate(n = rowSums(.[,-1])) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::slice_max(n) %>%
  dplyr::ungroup() %>%
  # Duplicated genes with 0 expression in all samples still remain, remove them
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::select(everything(), -n) %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL)) %>%
  dplyr::filter(SYMBOL %in% GOI$SYMBOL) %>%
  data.frame()

tcga_metadata <- read.xlsx("/hpc/home/kailasamms/scratch/biomarker/TCGA_Metadata.xlsx")

# Nat Comms
celltype <- "Epithelial"
seurat_results <- "/hpc/home/kailasamms/scratch/scRNASeq_Simon_old/results_seurat/"
integrated_seurat <- readRDS(paste0(seurat_results, "integrated_seurat_snn",
                                    dplyr::if_else(is.null(celltype), ".rds", paste0("_", celltype, ".rds"))))

male_integrated_seurat <- subset(x=integrated_seurat,
                                 Sex == "Male")
female_integrated_seurat <- subset(x=integrated_seurat,
                                   Sex == "Female")

male_nat_comm <- male_integrated_seurat@assays$RNA@data %>%
  data.frame() %>%
  tibble::rownames_to_column("SYMBOL") %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL)) %>%
  dplyr::filter(SYMBOL %in% GOI$SYMBOL)

female_nat_comm <- female_integrated_seurat@assays$RNA@data %>%
  data.frame() %>%
  tibble::rownames_to_column("SYMBOL") %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL)) %>%
  dplyr::filter(SYMBOL %in% GOI$SYMBOL)

male_nat_comm_metadata <- male_integrated_seurat@meta.data %>%
  dplyr::rename(Sample1 = Sample) %>%
  tibble::rownames_to_column("Sample")

female_nat_comm_metadata <- female_integrated_seurat@meta.data %>%
  dplyr::rename(Sample1 = Sample) %>%
  tibble::rownames_to_column("Sample")

# DepMap
model <- read.csv("/hpc/home/kailasamms/scratch/biomarker/Model.csv")
model <- model %>% dplyr::filter(DepmapModelType=="BLCA")

profile <- read.csv("/hpc/home/kailasamms/scratch/biomarker/OmicsProfiles.csv")
profile <- dplyr::inner_join(profile, model, by=c("ModelID"="ModelID")) %>%
  dplyr::select(StrippedCellLineName, ModelID, ProfileID, Sex, PrimaryOrMetastasis)

ccle <- read.csv("/hpc/home/kailasamms/scratch/biomarker/OmicsExpressionAllGenesTPMLogp1Profile.csv")
colnames(ccle) <- gsub("..ENSG.*$","",colnames(ccle))
colnames(ccle)[1] <- "ProfileID"
ccle <- ccle %>% dplyr::select(ProfileID, intersect(GOI$SYMBOL, make.names(colnames(ccle))))
ccle <- dplyr::inner_join(profile, ccle, by=c("ProfileID"="ProfileID")) 

ccle_metadata <- ccle[,1:5] %>%
  dplyr::rename(Sample = StrippedCellLineName)

ccle <- ccle[,-c(2,3,4,5)] %>% 
  tibble::column_to_rownames("StrippedCellLineName") %>%
  t() %>%
  data.frame() %>%
  tibble::rownames_to_column("SYMBOL") %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL))

# Find good Y genes
for (file_suffix in c("tcga", "ccle", "nat_comm")){
  
  results_path <- "/hpc/home/kailasamms/"
  # Save the clustered scores in xlsx
  wb <- openxlsx::createWorkbook()
  
  if(file_suffix != "nat_comm"){
    
    metadata_column <- get(paste0(file_suffix, "_metadata"))
    normalized_counts <- get(file_suffix)
    colnames(normalized_counts)[1] <- "SYMBOL"
    colnames(normalized_counts) <- make.names(colnames(normalized_counts))
    
    for (sex in c("Male", "Female")){
      
      samples <- metadata_column %>%
        dplyr::filter(Sex == sex)
      
      avg_exp <- rowMeans(as.matrix(normalized_counts[,make.names(samples$Sample)])) %>%
        data.frame()
      colnames(avg_exp)[1] <- "avg_exp"
      rownames(avg_exp) <- normalized_counts$SYMBOL
      avg_exp <- avg_exp %>% 
        tibble::rownames_to_column("SYMBOL")
      
      # median_exp <- rowMedians(x = as.matrix(normalized_counts[,make.names(samples$Sample)])) %>%
      #   data.frame() 
      # colnames(median_exp)[1] <- "median_exp"
      # rownames(median_exp) <- normalized_counts$SYMBOL
      # median_exp <- median_exp %>% 
      #   tibble::rownames_to_column("SYMBOL")
      
      quantile_exp <- rowQuantiles(x = as.matrix(normalized_counts[,make.names(samples$Sample)])) %>%
        data.frame()
      colnames(quantile_exp) <- c("Q0", "Q25", "Q50", "Q75", "Q100")
      rownames(quantile_exp) <- normalized_counts$SYMBOL
      quantile_exp <- quantile_exp %>% 
        tibble::rownames_to_column("SYMBOL")
      
      percent_exp <- normalized_counts %>% 
        tibble::column_to_rownames("SYMBOL") %>% 
        t() %>% 
        scale() %>%
        t() %>%
        data.frame() %>%
        tibble::rownames_to_column("SYMBOL") %>%
        dplyr::mutate_if(is.numeric, ~ dplyr::if_else(. > 0, 1, 0)) %>%
        dplyr::select(SYMBOL, make.names(samples$Sample)) %>%
        dplyr::mutate(percent_exp = rowMeans(select_if(., is.numeric),na.rm=TRUE)*100) %>%
        dplyr::select(SYMBOL, percent_exp)
      
      avg <- normalized_counts %>% 
        dplyr::select(SYMBOL, make.names(samples$Sample)) %>%
        dplyr::left_join(GOI %>% dplyr::select(SYMBOL, CHR, BIOTYPE), by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::left_join(avg_exp,   by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::left_join(percent_exp, by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::left_join(quantile_exp,   by=c("SYMBOL"="SYMBOL")) %>%
        #dplyr::left_join(median_exp,  by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::select(SYMBOL, CHR, BIOTYPE, avg_exp, percent_exp, Q0, Q25, Q50, Q75, Q100) #,everything())
      
      openxlsx::addWorksheet(wb, sheetName = sex)
      openxlsx::writeData(wb, sheet = sex, x = avg, rowNames = FALSE)
    }
    openxlsx::saveWorkbook(wb, file = paste0(results_path, file_suffix, "_Summary.xlsx"), 
                           overwrite = TRUE)
  }
  
  if (file_suffix == "nat_comm"){ 
    
    for (sex in c("Male", "Female")){
      
      normalized_counts <- get(paste0(tolower(sex), "_", file_suffix))
      
      avg_exp <- rowMeans(x = as.matrix(normalized_counts[,-1])) %>%
        data.frame()
      colnames(avg_exp)[1] <- "avg_exp"
      rownames(avg_exp) <- normalized_counts$SYMBOL
      avg_exp <- avg_exp %>% 
        tibble::rownames_to_column("SYMBOL")
      
      # median_exp <- rowMedians(x = as.matrix(normalized_counts[,-1])) %>%
      #   data.frame() 
      # colnames(median_exp)[1] <- "median_exp"
      # rownames(median_exp) <- normalized_counts$SYMBOL
      # median_exp <- median_exp %>% 
      #   tibble::rownames_to_column("SYMBOL")
      
      quantile_exp <- rowQuantiles(x = as.matrix(normalized_counts[,-1])) %>%
        data.frame()
      colnames(quantile_exp) <- c("Q0", "Q25", "Q50", "Q75", "Q100")
      rownames(quantile_exp) <- normalized_counts$SYMBOL
      quantile_exp <- quantile_exp %>% 
        tibble::rownames_to_column("SYMBOL")
      
      percent_exp <- normalized_counts %>% 
        tibble::column_to_rownames("SYMBOL") %>% 
        t() %>% 
        scale() %>%
        t()
      
      percent_exp[percent_exp > 0] <- 1
      percent_exp[percent_exp <= 0] <- 0
      
      percent_exp <- percent_exp %>%
        data.frame() %>%
        dplyr::mutate(percent_exp = rowMeans(., na.rm=TRUE)*100) %>%
        tibble::rownames_to_column("SYMBOL") %>%
        dplyr::select(SYMBOL, percent_exp)
      
      avg <- normalized_counts %>% 
        dplyr::left_join(GOI %>% dplyr::select(SYMBOL, CHR, BIOTYPE), by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::left_join(avg_exp,   by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::left_join(percent_exp, by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::left_join(quantile_exp,   by=c("SYMBOL"="SYMBOL")) %>%
        #dplyr::left_join(median_exp,  by=c("SYMBOL"="SYMBOL")) %>%
        dplyr::select(SYMBOL, CHR, BIOTYPE, avg_exp, percent_exp, Q0, Q25, Q50, Q75, Q100)
      
      openxlsx::addWorksheet(wb, sheetName = sex)
      openxlsx::writeData(wb, sheet = sex, x = avg, rowNames = FALSE)
    }
    openxlsx::saveWorkbook(wb, file = paste0(results_path, file_suffix, "_Summary.xlsx"), 
                           overwrite = TRUE)
  }
}

# Find good reference gene
m <- read.xlsx("/hpc/home/kailasamms/scratch/biomarker/Mutated,CNV_Genes.xlsx", sheet="Mutated_Genes")
c <- read.xlsx("/hpc/home/kailasamms/scratch/biomarker/Mutated,CNV_Genes.xlsx", sheet="CNA_Genes")

goi <- setdiff(annotations$SYMBOL, c$Gene)
goi <- setdiff(goi, m$Gene)
goi <- as.data.frame(goi)
colnames(goi) <- "SYMBOL"

GOI <- annotations %>% 
  dplyr::filter(SYMBOL %in% goi$SYMBOL) %>%
  dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) ~ ENSEMBL_ID, !is.na(SYMBOL) ~ SYMBOL)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::mutate(SYMBOL = make.names(SYMBOL))

# Run lines 23 through 94





wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Unmutated")
openxlsx::writeData(wb, sheet = "Unmutated", x = goi, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = "C:/Users/KailasammS/Desktop/Unmutated_Genes.xlsx", 
                       overwrite = TRUE)

#..good reference genes
gse <- "TCGA_BLCA"
parent_path <- "C:/Users/KailasammS/Box/Saravana@cedars/10. Ongoing Projects/BBN project/BLCA_Cohorts_correct/"

# Import read_data
read_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, gse, "_Normalized.xlsx"))
colnames(read_data)[1] <- "SYMBOL"
read_data <- read_data %>% 
  dplyr::mutate(SYMBOL = make.names(SYMBOL, unique=TRUE)) %>%
  tibble::column_to_rownames("SYMBOL")

mean <- rowMeans(read_data)
sd <- apply(X=read_data,MARGIN=1,FUN=sd)
df <- data.frame(mean,sd)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Unmutated")
openxlsx::writeData(wb, sheet = "Unmutated", x = df, rowNames = TRUE)
openxlsx::saveWorkbook(wb, file = "C:/Users/KailasammS/Desktop/Unmutated_Genes.xlsx", 
                       overwrite = TRUE)

ccle <- ccle %>%
  tibble::column_to_rownames("SYMBOL")
mean <- rowMeans(ccle)
sd <- apply(X=ccle,MARGIN=1,FUN=sd)
df <- data.frame(mean,sd)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Unmutated")
openxlsx::writeData(wb, sheet = "Unmutated", x = df, rowNames = TRUE)
openxlsx::saveWorkbook(wb, file = "/hpc/home/kailasamms/Unmutated_Genes.xlsx", 
                       overwrite = TRUE)


# # Remove genes with zero expression in all samples
# read_data <- read_data[rowSums(read_data) != 0,]
#
# # Remove genes not expressed in all samples
# binary1 <- read_data
# binary1[binary1 > 0] <- 1
# sum(rowSums(binary1) >= 0.9*ncol(binary1))
# g1 <- rownames(binary1[rowSums(binary1) >= 0.9*ncol(binary1),])
# 
# # Remove low expressed genes
# binary2 <- read_data
# q <- quantile(as.vector(as.matrix(binary2)))
# binary2[binary2 < 1000] <- 0
# binary2[binary2 > 1000] <- 1
# sum(rowSums(binary2) >= 0.9*ncol(binary2))
# g2 <- rownames(binary2[rowSums(binary2) >= 0.9*ncol(binary2),])
# 
# g <- intersect(g1,g2)
# t <- read_data[g,]
# 
# 
# 
# t <- apply(X=t,MARGIN=1,FUN=sd)
# sort(t)[1:100]
# 
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb, sheetName = "Unmutated")
# openxlsx::writeData(wb, sheet = "Unmutated", x = as.data.frame(t), rowNames = TRUE)
# openxlsx::saveWorkbook(wb, file = "C:/Users/KailasammS/Desktop/Unmutated_Genes.xlsx", 
#                        overwrite = TRUE)
# 



