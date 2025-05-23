source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

species <- "Homo sapiens"
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Maryam/"
proj <- "RNASeq_Maryam"

annotations <- get_annotations(species)
method <- "STAR"

# # Consolidate raw counts from STAR into excel file
# # Download the counts from HPC and store RNASeq_Vera/counts/STAR_counts
# read_data <- compile_raw_counts(data_path, proj, method)

# If already consolidated in excel file, read_data MUST have "SYMBOL" column 
# with gene names/ensembl_id, entrez_id without any duplication
read_data <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, proj, "_raw.counts.xlsx"))

# meta_data MUST be xlsx file and have "Sample_ID" column with sample names
meta_data <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, proj, "_Metadata.xlsx")) %>%
  dplyr::distinct_at("Sample_ID", .keep_all=TRUE)

# Rename appropiate column in meta_data to "Sample_ID"
# meta_data <- meta_data %>% dplyr::rename(Sample_ID = GEO_Accession__exp_)

# Rename appropiate column in read_data to "SYMBOL"
# read_data <- read_data %>% dplyr::rename(SYMBOL = GeneID)
# colnames(read_data)[1] <- "SYMBOL"

# Define comparisons and design
Comparisons <- list(Variable    = c("Details"),
                    Target      = c("TERT KO Clones"),
                    Reference   = c("NTC Clones"),
                    lfc.cutoff  = 0,
                    padj.cutoff = 0.1,
                    design      = c("Details"),
                    deseq2.batch.correct = FALSE)

# in vitro and in vivo data are widely different. So, normalize them separately
meta_data <- meta_data %>% dplyr::filter(grepl(pattern="Clones", x=Details))

Comparisons <- list(Variable    = c("Details"),
                    Target      = c("TERT KO Tumor"),
                    Reference   = c("NTC Tumor"),
                    lfc.cutoff  = 0,
                    padj.cutoff = 0.1,
                    design      = c("Details"),
                    deseq2.batch.correct = FALSE)

# in vitro and in vivo data are widely different. So, normalize them separately
meta_data <- meta_data %>% dplyr::filter(grepl(pattern="Tumor", x=Details))

# Perform QC
meta_data <- prep_metadata(meta_data, read_data)
read_data <- prep_readdata(read_data, meta_data)
l <- check_data(read_data, meta_data)
meta_data <- l[[2]]
read_data <- l[[1]]

#******************************************************************************#
#          DATA QUALITY ASSESSMENT BY SAMPLE CLUSTERING & VISUALZATION         #
#******************************************************************************#

# Plot PCA
plotPCA_DESeq2(meta_data, read_data, Comparisons, data_path)
 
# meta_data1 <- meta_data %>% dplyr::filter(!grepl(pattern="Pool", x=Details))
# read_data1 <- read_data %>% dplyr::select(!contains("Pool"))
# meta_data1 <- prep_metadata(meta_data1, read_data1)
# read_data1 <- prep_readdata(read_data1, meta_data1)
# l <- check_data(read_data1, meta_data1)
# meta_data1 <- l[[2]]
# read_data1 <- l[[1]]
# plotPCA_DESeq2(meta_data1, read_data1, Comparisons, data_path)  # remove text label

#******************************************************************************#
#                         CALCULATE NORMALIZED COUNTS                          #
#******************************************************************************#

# Get normalized counts
norm <- norm_counts_DESeq2(meta_data, read_data, annotations, data_path)
norm_counts_sva(meta_data, read_data, sva.formula_string, annotations, data_path)

#******************************************************************************#
#                       DIFFERENTIAL EXPRESSION ANALYSIS                       #
#******************************************************************************#

# Perform differential analysis for each comparison
for (n in 1:length(Comparisons$Variable)){
  
  # Create design formula 
  formula_string <- paste0("~", Comparisons$design[n])
  combat.formula_string <- formula_string
  sva.formula_string <- formula_string
  deseq2.formula_string <- dplyr::case_when(Comparisons$deseq2.batch.correct == TRUE ~ paste0("~", paste("Batch", Comparisons$design[n], sep = "+")),
                                            TRUE ~ formula_string)
  
  # Get combat corrected raw reads
  read_data_combat <- batch_correct_combat(meta_data, read_data, combat.formula_string)
  
  # Perform DESeq2() using in-built batch modelling
  prefix <- proj
  approach <- "DESeq2"
  # Create DESeq2 object with appropiate variables in design
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=as.formula(deseq2.formula_string))
  design(dds) <- as.formula(deseq2.formula_string)
  dds <- run_deseq2(dds, meta_data, annotations, Comparisons, n, approach, prefix, data_path)
  
  suffix <- paste(Comparisons$Target[n], "vs", Comparisons$Reference[n], sep='_')
  suffix <- gsub(pattern="/", replacement="", x=suffix)
  plotMA_DESeq2(dds, suffix, data_path)
}

#******************************************************************************#
#                                 VOLCANO PLOT                                 #
#******************************************************************************#

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Maryam/"

# Read DEGs
df1 <- read.xlsx(paste0(data_path, "RNA seq. clones.TKOvsNTC.DEG.xlsx"))
df2 <- read.xlsx(paste0(data_path, "RNA seqTKO_TuvsNTC_Tu.DEG.xlsx"))
#df1 <- read.xlsx(paste0(data_path, "RNASeq_Maryam.DEGs.TERT KO Clones_vs_NTC Clones_DESeq2.xlsx"))
#df2 <- read.xlsx(paste0(data_path, "RNASeq_Maryam.DEGs.TERT KO Tumor_vs_NTC Tumor_DESeq2.xlsx"))

#df1 <- read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/FromBox/Saravana@cedars/05. Bioinformatics/RNASeq/Deepak/id_TERT KO Tumor_vs_NTC Tumor_DEGs.xlsx")
#df2 <- read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/FromBox/Saravana@cedars/05. Bioinformatics/RNASeq/Deepak/id_TERT KO Clones_vs_NTC Clones_DEGs.xlsx")

# Define volcano plot parameters
disp_genes <- intersect(df1 %>% dplyr::filter(padj < 0.05, log2FoldChange > 0.58) %>% dplyr::select(SYMBOL) %>% unlist(use.names=FALSE),
                        df2 %>% dplyr::filter(padj < 0.05, log2FoldChange > 0.58) %>% dplyr::select(SYMBOL) %>% unlist(use.names=FALSE))

file_suffix <- ""
output_path <- data_path
volcano_params <- list(Target           = "Tert KO",
                       Reference        = "Control", 
                       matrix_color     = "vrds",
                       color_disp_genes = TRUE,
                       label_disp_genes = FALSE,
                       label_top        = FALSE,
                       lfc.cutoff       = 0.58,
                       padj.cutoff      = 0.05)

# Input dataframe MUST have SYMBOL, padj, log2FoldChange columns
plot_volcano(df1, volcano_params, disp_genes, "in.vitro", data_path)

plot_volcano(df2, volcano_params, disp_genes, "in.vivo", data_path)

#******************************************************************************#
#                                    HEATMAP                                   #
#******************************************************************************#

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

proj <- "RNASeq_Maryam"
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Maryam/"

# Read DEGs
df1 <- read.xlsx(paste0(data_path, "RNA seq. clones.TKOvsNTC.DEG.xlsx"))
df2 <- read.xlsx(paste0(data_path, "RNA seqTKO_TuvsNTC_Tu.DEG.xlsx"))
#df1 <- read.xlsx(paste0(data_path, "RNASeq_Maryam.DEGs.TERT KO Clones_vs_NTC Clones_DESeq2.xlsx"))
#df2 <- read.xlsx(paste0(data_path, "RNASeq_Maryam.DEGs.TERT KO Tumor_vs_NTC Tumor_DESeq2.xlsx"))

# Heatmap parameters
disp_genes <- intersect(df1 %>% dplyr::filter(padj < 0.05, log2FoldChange > 0.58) %>% dplyr::select(SYMBOL) %>% unlist(use.names=FALSE),
                        df2 %>% dplyr::filter(padj < 0.05, log2FoldChange > 0.58) %>% dplyr::select(SYMBOL) %>% unlist(use.names=FALSE))
disp_genes <- c(disp_genes, "TERT")
file_suffix <- ""
output_path <- data_path
heatmap_params <- list(anno.row       = NULL,        # NULL, c("Group")
                       anno.column    = c("Condition"),
                       row.split      = NA,     
                       col.split      = NA, #c("Condition"),
                       row.cluster    = c("all"),    # c("alphabetical", "group", "all")
                       col.cluster    = c("group"),  # c("alphabetical", "group", "all")
                       log.transform  = TRUE,
                       scale          = TRUE,
                       border_color   = "white",
                       bar_width      = NA,              # NA , 5
                       bar_height     = NA,              # NA , 5
                       matrix_color   = "rdbu",          # c("vrds", "rdbu")
                       expr_legend    = TRUE,  # set FALSE if it overlaps with annotation legends
                       file_format    = "tiff")

metadata_row <- data.frame(SYMBOL = "", Gene = "")

# Plot for clones
plot_genes <- df1 %>% dplyr::filter(padj < 0.05) %>% dplyr::select(SYMBOL) %>% unlist(use.names=FALSE)

metadata_column <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, proj, "_Metadata.xlsx")) %>%
  dplyr::distinct_at("Sample_ID", .keep_all=TRUE) %>%
  dplyr::filter(grepl(pattern="Clones", x=Details))

norm_counts1 <- read.xlsx(paste0(data_path, "Norm_counts_TERT KO Clones_vs_NTC Clones_DESeq2.xlsx")) %>%
  dplyr::select(SYMBOL, metadata_column$Sample_ID) %>%
  dplyr::filter(SYMBOL %in% (df1 %>% dplyr::filter(padj < 0.05) %>% dplyr::select(SYMBOL) %>% unlist(use.names = FALSE))) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE)  %>%
  dplyr::filter(!(SYMBOL %in% c("Y_RNA", "Metazoa_SRP")))

plot_heatmap(norm_counts1, metadata_column, metadata_row, heatmap_params,
             plot_genes, disp_genes, "clones", output_path)

# Plot for tumors
plot_genes <- df2 %>% dplyr::filter(padj < 0.05) %>% dplyr::select(SYMBOL) %>% unlist(use.names=FALSE)

metadata_column <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, proj, "_Metadata.xlsx")) %>%
  dplyr::distinct_at("Sample_ID", .keep_all=TRUE) %>%
  dplyr::filter(grepl(pattern="Tumor", x=Details))

norm_counts2 <- read.xlsx(paste0(data_path, "Norm_counts_TERT KO Tumor_vs_NTC Tumor_DESeq2.xlsx")) %>%
  dplyr::select(SYMBOL, metadata_column$Sample_ID) %>%
  dplyr::filter(SYMBOL %in% (df2 %>% dplyr::filter(padj < 0.05) %>% dplyr::select(SYMBOL) %>% unlist(use.names = FALSE))) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE)  %>%
  dplyr::filter(!(SYMBOL %in% c("Y_RNA", "Metazoa_SRP")))

plot_heatmap(norm_counts2 , metadata_column, metadata_row, heatmap_params,
             plot_genes, disp_genes, "tumors", output_path)

#******************************************************************************#

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

proj <- "RNASeq_Maryam"
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Maryam/"

# Read DEGs
df1 <- read.xlsx(paste0(data_path, "RNA seq. clones.TKOvsNTC.DEG.xlsx"))
df2 <- read.xlsx(paste0(data_path, "RNA seqTKO_TuvsNTC_Tu.DEG.xlsx"))

disp_genes <- intersect(df1 %>% dplyr::filter(padj < 0.05, log2FoldChange > 0.58) %>% dplyr::select(SYMBOL) %>% unlist(use.names=FALSE),
                        df2 %>% dplyr::filter(padj < 0.05, log2FoldChange > 0.58) %>% dplyr::select(SYMBOL) %>% unlist(use.names=FALSE))
disp_genes <- c(disp_genes, "TERT")

df1 <- df1 %>% dplyr::filter(SYMBOL %in% disp_genes)
df2 <- df2 %>% dplyr::filter(SYMBOL %in% disp_genes)

# Create workbook to store results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName="Norm_counts_clones")
openxlsx::writeData(wb, sheet="Norm_counts_clones", x=df1, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName="Norm_counts_tumor")
openxlsx::writeData(wb, sheet="Norm_counts_tumor", x=df2, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(data_path, "Fig3B.xlsx"),
                       overwrite = TRUE)


