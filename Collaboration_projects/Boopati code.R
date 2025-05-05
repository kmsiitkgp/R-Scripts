########################## HEATMAP FOR BIPUL GENESETS ##########################

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Boopati/"

# Get DESEq2 normalized counts
rna <- read.xlsx(paste0(data_path, "Normalized_Counts_DESeq2_CRISPR_LOY.xlsx"))
rna <- rna %>% 
  dplyr::select(SYMBOL, contains(c("SCR_", "Y_"))) %>%
  #dplyr::select(SYMBOL, contains(c("Yneg", "Ypos"))) %>%
  # If there are duplicated genes, keep only data for highest expressing copy
  dplyr::mutate(n = rowSums(.[,-1])) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::slice_max(n) %>%
  dplyr::ungroup() %>%
  # Duplicated genes with 0 expression in all samples still remain, remove them
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::select(everything(), -n)
colnames(rna) <- c("SYMBOL",
                   "MB49_Ypos1","MB49_Ypos2","MB49_Ypos3",
                   "MB49_Yneg1","MB49_Yneg2", "MB49_Yneg3")
# colnames(rna) <- c("SYMBOL", 
#                    "MB49_Yneg1","MB49_Yneg2", "MB49_Yneg3", "MB49_Yneg4","MB49_Yneg5", 
#                    "MB49_Ypos1","MB49_Ypos2","MB49_Ypos3","MB49_Ypos4","MB49_Ypos5") 

# Read Bipul gene sets
gene_sets <- read.xlsx(paste0(data_path, "Bipul_Gene_Sets.xlsx"))

# Format the matrix for heatmap
mat <- rna %>% 
  dplyr::filter(SYMBOL %in% gene_sets$Mouse) %>%
  tibble::column_to_rownames("SYMBOL")

mat <- log(1+mat, base = 2)
mat <- mat %>% t() %>% scale() %>% t()
mat[is.na(mat)] <- 0
mat <- mat[rowSums(mat) != 0,]

# Row annotation
row_ann <- gene_sets %>%
  dplyr::select(Mouse, GROUP) %>%
  dplyr::mutate(Mouse = make.names(Mouse), Dummy = 0) 

row_elements <- row_ann %>% 
  dplyr::distinct_at("GROUP", .keep_all = TRUE) %>%
  dplyr::select(GROUP) %>%
  unlist(use.names=FALSE)

# # Cluster within each pathway
# row_order <- c()
# for (g in row_elements){
#   temp_mat <- mat[rownames(row_ann)[which(row_ann$GROUP == g)],]
#   rowclust <- hclust(dist(temp_mat))
#   row_order <- c(row_order,rownames(temp_mat[rowclust$order,]))
# }

# # Re-order rows
# mat <- mat[row_order,]

# # Gaps
# element_names <- c()
# element_counts <- c()
# c <- 0
# anno_rows <- "GROUP"
# for (i in 1:nrow(row_ann)){
#   if (!(row_ann[,anno_rows][i] %in% element_names)){
#     element_names <- c(element_names, row_ann[,anno_rows][i])
#     element_counts <- c(element_counts, c)
#     c <- 1
#   } else{
#     c <- c+1
#   }
# }
# element_counts <- c(element_counts,c)
# gaps_row <- data.frame("Description" = element_names, n = element_counts[-1]) %>%
#   dplyr::mutate(n = cumsum(n)) %>%
#   dplyr::select(n) %>%
#   unlist(use.names = FALSE)

# Save 
wb <- openxlsx::createWorkbook()

for (r in row_elements){
  
  g <- gene_sets %>% dplyr::filter(GROUP == r ) %>% dplyr::select(Mouse) %>% unlist(use.names=FALSE)
  subset_mat <- mat[rownames(mat) %in% g,]
  
  rowclust <- hclust(dist(subset_mat))
  colclust <- hclust(dist(t(subset_mat)))
  subset_mat <- subset_mat[rowclust$order,colclust$order]
  
  # Breaks
  breaks <- c(seq(from = -1, to = 0, length.out = 50), seq(from = 1/100, to = 1, length.out = 50))
  breaks <- c(seq(from = min(mat), to = 0, length.out = 50), seq(from = max(mat)/100, to = 1, length.out = 50))
  vrds <- viridis_pal()(100)
  rdbu <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
  
  for (p in c("vrds", "rdbu")){
    pheatmap(subset_mat, 
             scale = "column",
             border_color = NA,
             cluster_rows = T,
             cluster_cols = T,
             annotation_row = row_ann %>% 
               dplyr::filter(GROUP == r) %>% 
               dplyr::distinct_at("Mouse", .keep_all=TRUE) %>% 
               tibble::column_to_rownames("Mouse") %>%
               dplyr::select(GROUP), 
             show_rownames = F, show_colnames = T, annotation_legend = TRUE,
             breaks = breaks, 
             color = get(p),
             gaps_row = gaps_row, 
             treeheight_row = 0,
             treeheight_col = 0,
             filename=paste0(data_path, "Heatmap_Bipul_gene_sets_", r, "_", p, ".tiff"))
  }
  
  openxlsx::addWorksheet(wb, sheetName = paste0(r, "_matrix"))
  openxlsx::writeData(wb, sheet = paste0(r, "_matrix"), x = subset_mat, rowNames = TRUE)
}
openxlsx::saveWorkbook(wb = wb,
                       file = paste0(data_path, "Heatmaps using CRISPR LOY RNA.xlsx"),
                       overwrite = TRUE)

######################### METABOLITE DATA PREPARATION ##########################

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Boopati/"

# Combine metabolite data from all assays and get proper metabolite names from Metaboanalyst

data1 <- read.xlsx(paste0(data_path, "Data_metabolomics_MB49_Yneg_vs_Ypos_Repeat1.xlsx"))
data2 <- read.xlsx(paste0(data_path, "Data_metabolomics_MB49Yneg_DDR2KO_vs_SCR_Repeat1.xlsx"))
data3 <- read.xlsx(paste0(data_path, "Data_metabolomics_MB49Yneg_DDR2KO_vs_SCR_Repeat2.xlsx"),
                   sheet="Global metabolomics")

data1 <- data1[,c(2:9,13:21)]
data2 <- data2[,c(2:9,13:18)]
data3 <- data3[,c(2:4,8:13)]

colnames(data1) <- data1[2,]
colnames(data2) <- data2[1,]
colnames(data3)[1:3] <- data3[2,c(1:3)]

data1 <- data1[-c(1:2),]
data2 <- data2[-1,]
data3 <- data3[-c(1:2),]

colnames(data1)[9:17] <- c("MB49_Parental1", "MB49_Parental2", "MB49_Parental3", 
                           "MB49_Ypos1",  "MB49_Ypos2", "MB49_Ypos3",
                           "MB49_Yneg1", "MB49_Yneg2", "MB49_Yneg3")
colnames(data2)[9:14] <- c("MB49_Yneg_scr1",  "MB49_Yneg_scr2", "MB49_Yneg_scr3",
                           "MB49_Yneg_DDR2KO1", "MB49_Yneg_DDR2KO2", "MB49_Yneg_DDR2KO3")
colnames(data3)[4:9] <- c("MB49_Yneg_scr1",  "MB49_Yneg_scr2", "MB49_Yneg_scr3",
                          "MB49_Yneg_DDR2KO1", "MB49_Yneg_DDR2KO2", "MB49_Yneg_DDR2KO3")

metabolites1 <- data1[,c(1:8)]
metabolites2 <- data2[,c(1:8)]
metabolites3 <- data3[,c(1:3)]

full_metabolites <- metabolites1 %>%
  dplyr::bind_rows(metabolites2) %>%
  dplyr::bind_rows(metabolites3) %>%
  dplyr::distinct_at("CmpdID", .keep_all = TRUE) %>%
  dplyr::filter(!is.na(compound))

# Save the metabolites
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Metabolites")
openxlsx::writeData(wb, sheet = "Metabolites", x = full_metabolites, rowNames = FALSE)
openxlsx::saveWorkbook(wb = wb,
                       file = paste0(data_path, "MB49_Metabolites_incomplete.xlsx"),
                       overwrite = TRUE)

# Using HMDBID, KEGGID, etc from above xlsx file, create a new file that 
# contains proper common name of metabolites along with their different ids
# This needs to be done manually
name_mapping <- read.xlsx(paste0(data_path, "Metaboanalyst_name_map.xlsx"))

# Recursively merge name_mapping with pathway info for the metabolites
name_mapping1 <- name_mapping %>%
  dplyr::left_join(full_metabolites %>% dplyr::select(PubChem, Pathway, compound) %>% dplyr::filter(!is.na(PubChem)) %>%
                     dplyr::mutate(PubChem = as.numeric(PubChem)),
                   by=c("PubChem"="PubChem"))

name_mapping2 <- name_mapping1 %>%
  dplyr::filter(is.na(Pathway)) %>%
  dplyr::select(everything(), -c(Pathway, compound)) %>%
  dplyr::left_join(full_metabolites %>% dplyr::select(CmpdID, Pathway, compound) %>% dplyr::filter(!is.na(CmpdID)),
                   by=c("KEGG"="CmpdID"))

name_mapping3 <- name_mapping2 %>%
  dplyr::filter(is.na(Pathway)) %>%
  dplyr::select(everything(), -c(Pathway, compound)) %>%
  dplyr::left_join(full_metabolites %>% dplyr::select(CmpdID, Pathway, compound) %>% dplyr::filter(!is.na(CmpdID)),
                   by=c("Query"="CmpdID"))

sum(is.na(name_mapping1$Pathway))
sum(is.na(name_mapping2$Pathway))
sum(is.na(name_mapping3$Pathway))

name_mapping1 <- name_mapping1 %>% dplyr::filter(!is.na(Pathway))
name_mapping2 <- name_mapping2 %>% dplyr::filter(!is.na(Pathway))
name_mapping3 <- name_mapping3 %>% dplyr::filter(!is.na(Pathway))

metadata <- dplyr::bind_rows(name_mapping1, name_mapping2, name_mapping3) %>%
  dplyr::rename(Initial_compound = compound, Compound = Match) %>%
  dplyr::select(Initial_compound, Compound, PubChem, HMDB, KEGG, ChEBI, METLIN, Pathway)

# Replace the metabolites names in data1, data2 and data3 with the correct info
# from metadata
data1 <- data1 %>% 
  dplyr::select(everything(), -c(CmpdID, PubChem, HMDB, KEGG, ChEBI, METLIN, Pathway)) %>%
  dplyr::left_join(metadata, by=c("compound" = "Initial_compound")) %>%
  dplyr::select(Compound, everything(), -compound)
data2 <- data2 %>% 
  dplyr::select(everything(), -c(CmpdID, PubChem, HMDB, KEGG, ChEBI, METLIN, Pathway)) %>%
  dplyr::left_join(metadata, by=c("compound" = "Initial_compound")) %>%
  dplyr::select(Compound, everything(), -compound)
data3 <- data3 %>% 
  dplyr::select(everything(), -c(CmpdID, Pathway)) %>%
  dplyr::left_join(metadata, by=c("compound" = "Initial_compound")) %>%
  dplyr::select(Compound, everything(), -compound)

# Get proper metabolite names from Metaboanalyst for CCLE data 
cell_lines <- read.xlsx(paste0(data_path, "41591_2019_404_MOESM2_ESM.xlsx"),
                        sheet="1-cell line annotations") %>%
  dplyr::filter(Gender =="male", Classifications == "urinary_tract")

ccle <- read.xlsx(paste0(data_path, "41591_2019_404_MOESM2_ESM.xlsx"),
                  sheet="1-raw data") %>%
  tidyr::separate_wider_delim(cols=X1, delim="_", names=c("Name"), too_many="drop", too_few="align_start") %>%
  dplyr::filter(Name %in% cell_lines$Name) %>%
  tibble::column_to_rownames("Name") %>%
  t() %>%
  data.frame()

# Save the incomplete metabolites
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Metabolites")
openxlsx::writeData(wb, sheet = "Metabolites", x = ccle, rowNames = TRUE)
openxlsx::saveWorkbook(wb = wb,
                       file = paste0(data_path, "CCLE_Metabolites_incomplete.xlsx"),
                       overwrite = TRUE)

# Get proper metabolite names from metaboanalyst
name_mapping  <- read.xlsx(paste0(data_path, "CCLE_Metaboanalyst_name_map.xlsx"))

# Replace metabolite names in CCLE with proper names from metaboanalyst
ccle <- ccle %>%
  tibble::rownames_to_column("Query") %>%
  dplyr::left_join(name_mapping, by=c("Query"="Query")) %>%
  dplyr::mutate(Compound = dplyr::case_when(is.na(Match) ~ Query,
                                            TRUE ~ Match)) %>%
  dplyr::select(Compound, everything(), -c(Query, Match, SMILES, Comment))

# Save the metabolites with proper name and ids
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Metabolites")
openxlsx::writeData(wb, sheet = "Metabolites", x = metadata, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName = "MB49_Yneg_vs_Ypos_Repeat1")
openxlsx::writeData(wb, sheet = "MB49_Yneg_vs_Ypos_Repeat1", x = data1, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName = "MB49Yneg_DDR2KO_vs_SCR_Repeat1")
openxlsx::writeData(wb, sheet = "MB49Yneg_DDR2KO_vs_SCR_Repeat1", x = data2, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName = "MB49Yneg_DDR2KO_vs_SCR_Repeat2")
openxlsx::writeData(wb, sheet = "MB49Yneg_DDR2KO_vs_SCR_Repeat2", x = data3, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName = "CCLE")
openxlsx::writeData(wb, sheet = "CCLE", x = ccle, rowNames = FALSE)
openxlsx::saveWorkbook(wb = wb,
                       file = paste0(data_path, "Metabolites_complete.xlsx"),
                       overwrite = TRUE)

###################### CCLE vs Y+Y- Metabolite Correlation #####################

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Boopati/"

# Get CCLE data with correct metabolite names
ccle <- read.xlsx(paste0(data_path, "Metabolites_complete.xlsx"),
                  sheet="CCLE") %>%
  dplyr::select(everything(), -c(PubChem, HMDB, KEGG, ChEBI, METLIN)) 
y <- read.xlsx(paste0(data_path, "Metabolites_complete.xlsx"),
               sheet="MB49_Yneg_vs_Ypos_Repeat1") %>%
  dplyr::select(everything(), -c(PubChem, HMDB, KEGG, ChEBI, METLIN, Pathway))

# Remove parental lines
y <- y [,-c(2:4)]

# Merge
data <- y %>%
  dplyr::inner_join(ccle, by=c("Compound"="Compound")) %>%
  tibble::column_to_rownames("Compound") %>%
  dplyr::mutate(across(.cols=everything(), .fns= as.numeric))

# Normalize metabolites by total
total <- colSums(data, na.rm=TRUE)
min_val <- min(total, na.rm=TRUE)
scale_factor <- total/min_val
norm_data <- data/scale_factor
norm_data <- na.omit(norm_data)
norm_data <- norm_data%>%
  dplyr::mutate(MB49_Ypos = rowMeans(norm_data[,1:3], na.rm=TRUE),
                MB49_Yneg = rowMeans(norm_data[,4:6], na.rm=TRUE))
norm_data <- norm_data[, -c(1:6)]
y_cor_mat <- round(cor(x = norm_data, method = "pearson"),2)

# Breaks
breaks <- c(seq(from = min(y_cor_mat), to = 0, length.out = 50), seq(from = max(y_cor_mat)/100, to = 1, length.out = 50))
vrds <- viridis_pal()(100)
rdbu <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

for (p in c("vrds", "rdbu")){
  select_ccle <- c("RT4", "X5637", "HT1197", "X647V", "UMUC1", "MB49_Ypos")
  pheatmap(y_cor_mat[select_ccle, select_ccle], 
           breaks = breaks, 
           color = get(p), 
           treeheight_row = 0,
           treeheight_col = 0,
           filename=paste0(data_path, "Fig 1_Ypos_", p, ".tiff"))
  select_ccle <- c("RT4", "X5637", "HT1197", "X647V", "UMUC1", "MB49_Yneg")
  pheatmap(y_cor_mat[select_ccle, select_ccle], 
           breaks = breaks, 
           color = get(p), 
           treeheight_row = 0,
           treeheight_col = 0,
           filename=paste0(data_path, "Fig 1_Yneg_", p, ".tiff"))
  pheatmap(y_cor_mat, 
           breaks = breaks, 
           color = get(p), 
           treeheight_row = 0,
           treeheight_col = 0,
           filename=paste0(data_path, "Fig 1_Yall_", p, ".tiff"))
}

# Save the correlation matrix
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "CCLE Correlation")
openxlsx::writeData(wb, sheet = "CCLE Correlation", x = y_cor_mat, rowNames = TRUE)
openxlsx::saveWorkbook(wb = wb,
                       file = paste0(data_path, "Fig1_Corr_Matrix.xlsx"),
                       overwrite = TRUE)

######################## Fig 2A: RNA vs  Metabolite Correlation ########################

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Boopati/"

# Get Y+Y- metabolite data
y_full <- read.xlsx(paste0(data_path, "Metabolites_complete.xlsx"),
                    sheet="MB49_Yneg_vs_Ypos_Repeat1")

y <- y_full %>%
  dplyr::select(everything(), -c(PubChem, HMDB, KEGG, ChEBI, METLIN, Pathway))

# Keep only 85 or 41 DEMs from Bipul
DEMs <- read.xlsx(paste0(data_path, "Bipul_DEMs.xlsx"))

# Remove parental lines
y <- y[,-c(2:4)]
y <- y %>%
  tibble::column_to_rownames("Compound") %>%
  dplyr::mutate(across(.cols=everything(), .fns= as.numeric))

# Normalize metabolics data
total <- colSums(y, na.rm=TRUE)
min_val <- min(y, na.rm=TRUE)
total <- total/min_val
y <- y/total
y <- na.omit(y)
y <- y %>%
  tibble::rownames_to_column("SYMBOL")

# Remove non DEMs
y <- y %>%
  dplyr::filter(SYMBOL %in% DEMs$Compound)

# Remove DEMs not in Amino acid, Glycolysis and Nucleotide pathways
y <- y %>%
  dplyr::filter(SYMBOL %in% (y_full %>%
                  dplyr::filter(Pathway %in% c("Amino acids", "Glycolysis", "Nucleotides")) %>%
                  dplyr::select(Compound) %>%
                  unlist(use.names=FALSE)))

# Get DESEq2 normalized counts
rna <- read.xlsx(paste0(data_path, "Normalized_Counts_DESeq2_CRISPR_LOY.xlsx"))

rna <- rna %>% 
  dplyr::select(SYMBOL, SCR_KO1, SCR_KO2, SCR_KO3, Y_KO1, Y_KO2, Y_KO3) %>%
  #dplyr::select(SYMBOL, SCR_KO1, SCR_KO2, SCR_KO3, DDR2_KO1, DDR2_KO2, DDR2_KO3) %>%
  # If there are duplicated genes, keep only data for highest expressing copy
  dplyr::mutate(n = rowSums(.[,-1])) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::slice_max(n) %>%
  dplyr::ungroup() %>%
  # Duplicated genes with 0 expression in all samples still remain, remove them
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  dplyr::select(everything(), -n)
colnames(rna) <- c("SYMBOL", "MB49_Ypos1","MB49_Ypos2","MB49_Ypos3","MB49_Yneg1","MB49_Yneg2", "MB49_Yneg3") 
#colnames(rna) <- c("SYMBOL", "MB49_Yneg_scr1", "MB49_Yneg_scr2", "MB49_Yneg_scr3", "MB49_Yneg_DDR2KO1", "MB49_Yneg_DDR2KO2", "MB49_Yneg_DDR2KO3") 

# Get list of DEGs
DEGs <- read.xlsx(paste0(data_path, "Results_id_CRISPR_Yneg_vs_CRISPR_YPos_DESeq2_modelled__DEGs.xlsx")) %>%
  dplyr::filter(padj < 0.05) 

# Filter out non-DEGs
rna <- rna %>%
  dplyr::filter(SYMBOL %in% make.names(DEGs$SYMBOL))

# Read Bipul gene sets
gene_sets <- read.xlsx(paste0(data_path, "Bipul_Gene_Sets.xlsx"))
  #dplyr::filter(GROUP %in% c("KEGG_Glycolysis", "Metabolism"))

# Filter out DEGs absent in gene sets
rna <- rna %>%
  dplyr::filter(SYMBOL %in% make.names(gene_sets$Mouse))

# Merge RNA and metabolite data
data <- dplyr::bind_rows(y, rna)

# Calculate correlation
y_all <- data %>%
  tibble::column_to_rownames("SYMBOL") %>%
  t() %>%
  data.frame()
y_all <- y_all[,colSums(y_all) != 0]
y_all_cor_mat <- psych::corr.test(x = y_all, method = "pearson", ci=FALSE)

# Extract ypos data
ypos <- data[,1:4] %>%
  tibble::column_to_rownames("SYMBOL") %>%
  t() %>%
  data.frame()
ypos <- ypos[,colSums(ypos) != 0]
#ypos_cor_mat <- round(cor(x = ypos, method = spearman),2)
ypos_cor_mat <- psych::corr.test(x = ypos, method = "pearson", ci=FALSE)

# Extract yneg data
yneg <- data[,c(1,5,6,7)] %>%
  tibble::column_to_rownames("SYMBOL") %>%
  t() %>%
  data.frame()
yneg <- yneg[,colSums(yneg) != 0]
#yneg_cor_mat <- round(cor(x = yneg, method = "spearman"),2)
yneg_cor_mat <- psych::corr.test(x = yneg, method = "pearson", ci=FALSE)

# Row annotation
row_ann <- y_full %>%
  dplyr::filter(Compound %in% DEMs$Compound) %>%
  dplyr::select(Compound, Pathway) %>%
  dplyr::mutate(Compound = make.names(Compound), Dummy = 0) %>%
  tibble::column_to_rownames("Compound") %>%
  dplyr::filter(Pathway %in% c("Amino acids", "Glycolysis", "Nucleotides"))

row_elements <- row_ann %>% 
  dplyr::count(Pathway) %>%
  dplyr::filter(n>1) %>%
  dplyr::select(Pathway) %>%
  unlist(use.names=FALSE)

# Column annotation
col_ann <- gene_sets %>%
  dplyr::filter(Mouse %in% DEGs$SYMBOL) %>%
  dplyr::select(Mouse, GROUP) %>%
  dplyr::mutate(Mouse = make.names(Mouse), Dummy = 0) %>%
  dplyr::distinct_at("Mouse", .keep_all = TRUE) %>%
  tibble::column_to_rownames("Mouse")

col_elements <- col_ann %>%
  dplyr::count(GROUP) %>%
  dplyr::filter(n>1) %>%
  dplyr::select(GROUP) %>%
  unlist(use.names=FALSE)

vrds <- viridis_pal()(100)
rdbu <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

# Save the correlation matrix
wb <- openxlsx::createWorkbook()
p <- "rdbu"
cluster <- "unclustered"  #unclustered
ypos_matrix <- data.frame(ypos_cor_mat$r)
yneg_matrix <- data.frame(yneg_cor_mat$r)

for (m in c("Amino acids", "Glycolysis", "Nucleotides")){
  for (d in c("yneg_matrix", "ypos_matrix")){
    
    dat <- get(d)
    
    if (m == "Amino acids"){
      pathway <- "KEGG_ARGININE_AND_PROLINE_METABOLISM"
    }else if(m == "Glycolysis"){
      pathway <- "KEGG_GLYCOLYSIS_GLUCONEOGENESIS"
    }else{
      pathway <- "KEGG_PURINE_METABOLISM"
    }
    
    dat <- dat[rownames(dat) %in% rownames(row_ann)[row_ann$Pathway == m],
               colnames(dat) %in% rownames(col_ann)[col_ann$GROUP == pathway]]
    
    # Breaks
    breaks <- c(seq(from = min(dat), to = 0, length.out = 50),
                seq(from = max(dat)/100, to = 1, length.out = 50))
    
    #mat <- dat[rownames(dat) %in% make.names(DEMs$Compound), colnames(dat) %in% rna$SYMBOL]
    mat <- dat
    
    if (d != "ypos_matrix" | cluster == "clustered"){
      # Cluster within each pathway
      row_order <- c()
      for (g in intersect(row_elements,m)){
        temp_mat <- mat[rownames(row_ann)[which(row_ann$Pathway == g)],]
        rowclust <- hclust(dist(temp_mat))
        row_order <- c(row_order,rownames(temp_mat[rowclust$order,]))
      }
      # # Add groups that have only 1gene/metabolite
      # row_order <- c(row_order, rownames(row_ann)[!(rownames(row_ann) %in% row_order)])
      
      # Cluster within each GROUP
      col_order <- c()
      for (g in intersect(col_elements,pathway)){
        temp_mat <- mat[,colnames(mat) %in% rownames(col_ann)[which(col_ann$GROUP == g)]]
        colclust <- hclust(dist(t(temp_mat)))
        col_order <- c(col_order,colnames(temp_mat[,colclust$order]))
      }
      mat <- mat[row_order, col_order]
      g <- ""
      #gaps_row <- (row_ann %>% dplyr::count(Pathway) %>% dplyr::mutate(n = cumsum(n)) %>% dplyr::select(n) %>% unlist(use.names = FALSE))
      #gaps_row <- gaps_row[gaps_row < nrow(mat)]
      #gaps_col <- (col_ann %>% dplyr::count(GROUP) %>% dplyr::mutate(n = cumsum(n)) %>% dplyr::select(n) %>% unlist(use.names = FALSE))
      #gaps_col <- gaps_col[gaps_col < ncol(mat)]
      pheatmap(mat, 
               cluster_rows = F,
               cluster_cols = F,
               border_color = NA,
               cellwidth = 5,
               cellheight = 5,
               #gaps_row = gaps_row, 
               #gaps_col = gaps_col,
               annotation_row = row_ann %>% dplyr::select(Pathway), 
               annotation_col = col_ann %>% dplyr::select(GROUP), 
               show_rownames = F, show_colnames = F, annotation_legend = TRUE,
               breaks = breaks, 
               color = get(p),
               annotation_names_row = FALSE, 
               annotation_names_col = FALSE,
               #treeheight_row = 0,
               #treeheight_col = 0,
               width = 7,
               height = 7,
               filename=paste0(data_path, "RNA-metabolite_", m, "_", d,"_", cluster, ".tiff"))
      
      openxlsx::addWorksheet(wb, sheetName = paste0(d, "_", m))
      openxlsx::writeData(wb, sheet = paste0(d, "_", m), x = mat, rowNames = TRUE)
    }else{
      mat <- mat[row_order, col_order]
      g <- ""
      #gaps_row <- (row_ann %>% dplyr::count(Pathway) %>% dplyr::mutate(n = cumsum(n)) %>% dplyr::select(n) %>% unlist(use.names = FALSE))
      #gaps_row <- gaps_row[gaps_row < nrow(mat)]
      #gaps_col <- (col_ann %>% dplyr::count(GROUP) %>% dplyr::mutate(n = cumsum(n)) %>% dplyr::select(n) %>% unlist(use.names = FALSE))
      #gaps_col <- gaps_col[gaps_col < ncol(mat)]
      pheatmap(mat, 
               cluster_rows = F,
               cluster_cols = F,
               border_color = NA,
               cellwidth = 5,
               cellheight = 5,
               #gaps_row = gaps_row, 
               #gaps_col = gaps_col,
               annotation_row = row_ann %>% dplyr::filter(Pathway == m) %>% dplyr::select(Pathway), 
               annotation_col = col_ann %>% dplyr::filter(GROUP == pathway) %>% dplyr::select(GROUP), 
               show_rownames = F, show_colnames = F, annotation_legend = TRUE,
               breaks = breaks, 
               color = get(p),
               annotation_names_row = FALSE, 
               annotation_names_col = FALSE,
               #treeheight_row = 0,
               #treeheight_col = 0,
               width = 7,
               height = 7,
               filename=paste0(data_path, "RNA-metabolite_", m, "_", d,"_", cluster, ".tiff"))
    }
  } 
}  
openxlsx::saveWorkbook(wb = wb,
                       file = paste0(data_path, "Fig2_Corr_Matrix.xlsx"),
                       overwrite = TRUE)



######## GSEA ANALYSIS

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Boopati/"
species <- "Mus musculus"

# Get annotations
annotations <- get_annotations(species)

# Get GMT files
gmt_files <- list.files("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets/Mouse/", full.names = TRUE)

# Import DEGs dataframe with "SYMBOL", "padj" and "log2FoldChange" columns
DEGs_df <- openxlsx::read.xlsx(xlsxFile = paste0(data_path, "Results_id_CRISPR_Yneg_vs_CRISPR_YPos_DESeq2_modelled__DEGs.xlsx"))

# Run fgsea and gsea on all collections. ora is usually run on hallmark collection
for (gmt_file in gmt_files){
  fgsea(DEGs_df, gmt_file, annotations, data_path)
}



# ################################# CIRCOS PLOT #################################
# 
# source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
# data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Boopati/"
# 
# #install.packages("Hmisc")
# #install.packages("corrplot")
# #install.packages("circlize")
# #install.packages("psych")
# library(Hmisc)
# library(corrplot)
# library(ggfortify)
# library(circlize)
# library(psych)
# 
# # Get Y+Y- metabolite data
# y_full <- read.xlsx(paste0(data_path, "Metabolites_complete.xlsx"),
#                     sheet="MB49_Yneg_vs_Ypos_Repeat1")
# 
# # Read Bipul gene sets
# gene_sets <- read.xlsx(paste0(data_path, "Bipul_Gene_Sets.xlsx"))
# 
# # Get correlation matrix from previous section
# mat <- read.xlsx(paste0(data_path, "Fig2_Corr_Matrix.xlsx"),
#                  sheet="Ypos_cor_matrix")
# p <- read.xlsx(paste0(data_path, "Fig2_Corr_Matrix.xlsx"),
#                sheet="Ypos_p_matrix") %>%
#   dplyr::rename(id = identity(1)) %>%
#   tibble::column_to_rownames("id")
# 
# p[upper.tri(p)] <- p[lower.tri(p)]
# 
# # Convert to long format
# mat_long <- mat %>%
#   data.frame() %>%
#   dplyr::rename(from = identity(1)) %>%
#   dplyr::filter(from %in% make.names(y_full$Compound)) %>%
#   dplyr::select(everything(), -c(all_of(make.names(y_full$Compound)))) %>%
#   tidyr::pivot_longer(names_to = "to", values_to = "cor", cols = c(everything(), -"from")) 
# 
# p_long <- p %>%
#   data.frame() %>%
#   tibble::rownames_to_column("from") %>%
#   dplyr::filter(from %in% make.names(y_full$Compound)) %>%
#   dplyr::select(everything(), -c(all_of(make.names(y_full$Compound)))) %>%
#   tidyr::pivot_longer(names_to = "to", values_to = "p", cols = c(everything(), -"from")) 
# 
# mat_p_long <- mat_long %>%
#   dplyr::left_join(p_long, by=c("from"="from", "to"="to")) %>%
#   dplyr::filter(p <=0.05)
# 
# # Row annotation
# row_ann <- y_full %>%
#   dplyr::select(Compound, Pathway) %>%
#   dplyr::mutate(Compound = make.names(Compound), Dummy = 0) %>%
#   tibble::column_to_rownames("Compound")
# 
# # Column annotation
# col_ann <- gene_sets %>%
#   dplyr::select(Mouse, GROUP) %>%
#   dplyr::mutate(Mouse = make.names(Mouse), Dummy = 0) %>%
#   dplyr::distinct_at("Mouse", .keep_all = TRUE) %>%
#   tibble::column_to_rownames("Mouse")
# 
# # Replace metabolites with pathways and genes with group
# data <- mat_p_long %>%
#   dplyr::left_join(row_ann %>% tibble::rownames_to_column("Compound") %>% dplyr::select(Compound, Pathway),
#                    by=c("from"="Compound")) %>%
#   dplyr::left_join(col_ann %>% tibble::rownames_to_column("Compound") %>% dplyr::select(Compound, GROUP),
#                    by=c("to"="Compound")) %>%
#   dplyr::select(Pathway, GROUP, cor)
# 
# # Plot chord diagram
# circos.clear()
# circos.par(gap.after=c(rep(1,length(rownames(data))-1),10,rep(1,length(colnames(data))-1),10))
# circlize::chordDiagram(x = data[,1:3], 
#                        #grid.col = mycolor,
#                        order = NULL,
#                        transparency = 0.25,  # 0.25 = 25% transparency
#                        directional = 1,      # 1  = direction is from 1st column to 2nd column; 
#                        # -1 = direction is from 2nd column to 1st column
#                        direction.type = c("arrows", "diffHeight"), 
#                        diffHeight  = -0.05,  # -ve value = link start site is closer to circle;
#                        # +ve value = link start site is farther from circle
#                        annotationTrack = "grid", 
#                        annotationTrackHeight = c(0.05, 0.1),
#                        link.arr.type = "big.arrow", 
#                        link.sort = TRUE, 
#                        link.largest.ontop = TRUE)
# 
# 
# 

# quant_norm <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(data)))
# rownames(quant_norm) <- rownames(data)
# colnames(quant_norm) <- colnames(data)
# quant_norm <- na.omit(quant_norm)
# quant_norm <- quant_norm %>%
#   dplyr::mutate(MB49_Ypos = rowMeans(quant_norm[,1:3], na.rm=TRUE),
#                 MB49_Yneg = rowMeans(quant_norm[,4:6], na.rm=TRUE))
# quant_norm <- quant_norm[, -c(1:6)]
# y_cor_mat <- round(cor(quant_norm),2)
# pheatmap(y_cor_mat)



# # PCA plot
# df <- data
# pca_res <- prcomp(t(na.omit(df)), scale = TRUE)
# autoplot(pca_res, label=TRUE)
# ggsave("1.jpg")

############################## Metabolite analysis #############################

# Use metaboanalyst to normalize by sum, log10 transform, perform pareto-scaling,
# calculate fold change, pvalue, FDR as described in Nature paper

# #************************GSEA on bulk RNASeq
# 
# parent_path <- "C:/Users/kailasamms/Box/Boopati_paper/Transcriptomics/"
# results_path <- parent_path
# 
# species <- "Mus musculus"
# gmt_files <- list.files("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets/Mouse/", full.names = TRUE)
# annotations <- get_annotations(species)
# 
# for (gmt_file in gmt_files){
#   
#   #DEGs_df <- read.xlsx(paste0(parent_path, "CRISPR_LOY/Results_id_Y_Negative_vs_Y_Positive_DESeq2_modelled__DEGs.xlsx"))
#   #DEGs_df <- read.xlsx(paste0(parent_path, "Natural_LOY/Results_id_Y_Negative_vs_Y_Positive_DESeq2_modelled__DEGs.xlsx"))
#   DEGs_df <- read.xlsx(paste0(parent_path, "DDR2KO/Results_id_KO_vs_Control_DESeq2_modelled__DEGs.xlsx"))
#   
#   fgsea(DEGs_df, gmt_file, annotations, results_path)
# }
# 
# #*********************Volcano plots of metabolites
# 
# parent_path <- "C:/Users/kailasamms/Box/Boopati_paper/Metabolomics/"
# results_path <- parent_path
# 
# # data <- read.xlsx(paste0(parent_path, "Data_metabolomics_DDR2KO_vs_SCR.xlsx")) %>%
# #   tibble::column_to_rownames("compound") %>%
# #   dplyr::select(everything(), -Pathway)
# # 
# # meta_data <- data.frame("Sample" = colnames(data), 
# #                         "Condition" = c("Y+","Y+", "Y+", "Y-", "Y-", "Y-"))
# # 
# # 
# # # Plot box plot before normalization
# # b_data <- log(1+data, base=2) %>%
# #   tibble::rownames_to_column("compound") %>%
# #   tidyr::pivot_longer(cols = !compound, names_to="Sample", values_to="Intensity")
# # boxplot(b_data$Intensity ~ b_data$Sample)
# # 
# # 
# # # Perform quantile normalization
# # quant_norm <- TRUE
# # n_data <- quantile_norm(data, meta_data, quant_norm)
# # 
# # # Plot box plot before normalization
# # b_data <- log(1+n_data, base=2) %>%
# #   tibble::rownames_to_column("compound") %>%
# #   tidyr::pivot_longer(cols = !compound, names_to="Sample", values_to="Intensity")
# # boxplot(b_data$Intensity ~ b_data$Sample)
# # 
# # # Calculate padj and log2FC
# # Target <- "Y-"
# # Reference <- "Y+"
# # n_data <- log(1+n_data, base=2)
# # log2_transformed_already <- FALSE
# # f_data <- calc_stats(n_data, meta_data, Target, Reference,log2_transformed_already)
# ypos_cor <- rcorr(as.matrix(ypos))
# ypos_cor_mat <- ypos_cor$r
# ypos_cor_matp <- ypos_cor$P

# corrplot(ypos_cor_mat, 
#          type = "upper", 
#          order = "hclust", 
#          tl.col = "black", 
#          tl.srt = 45)

# p <- flattenCorrMatrix(p) %>%
#   dplyr::rename(Ypos_cor = cor)
# n <- flattenCorrMatrix(n) %>%
#   dplyr::rename(Yneg_cor = cor)
# pn <- p %>%
#   dplyr::full_join(n, by=c("row"="row", "column"="column")) %>%
#   dplyr::mutate(label = paste0(row, "_", column)) %>%
#   dplyr::filter((Ypos_cor > 0.25 & Yneg_cor < -0.25) | (Ypos_cor < -0.25 & Yneg_cor > 0.25))
#   #               diff = Ypos_cor - Yneg_cor) %>%
#   # dplyr::filter(abs(diff) >= 0.5)
# 
# ggplot(pn, aes(x=Ypos_cor, y=Yneg_cor)) + 
#   geom_point(size = 0.01)
# 
# ggsave("1.jpg")