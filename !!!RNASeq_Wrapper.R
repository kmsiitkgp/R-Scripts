# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
# https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html

# NOTE: If any one group has high within-group variability in PCA plot, we 
# SHOULD exclude those samples by sub-setting before creating dds object and 
# calculating dispersion estimates. Else, use full dataset for modelling.

# NOTE: 
# One factor   : design ~ Condition
# Two factors  : design ~ Condition + Treatment + Condition:Treatment (OR) design ~ Condition*Treatment
# Three factors: design ~ Cell.Line + Condition + Treatment + Condition:Treatment +.. (OR) design ~ Cell.Line*Condition*Treatment
# Using * will include all possible interaction terms and is RECOMMENDED
# The results are the same irrespective of the order i.e 
# Cell.Line*Condition*Treatment gives same results as Cell.Line*Treatment*Condition

# NOTE: If you have Cell.Line, Treatment, Condition variables and want to 
# include them in design, ideally design MUST be ~ Cell.Line*Treatment*Condition.
# However, if one cell line doesnt have a specific treatment, then DESeq2 will
# throw "Full model martix is less than full rank" error. So, ALWAYS create a 
# new column "Comparisons" in meta_data and use design ~ Comparisons ALWAYS.
# Populate the "Comparisons" column in meta_data by concatenating as below:
# paste0(Cell.Line, Treatment, Condition, sep=".") 

# NOTE: contrast terms CANNOT start with numbers. So, if you have cell line name
# "22RV1" rename it as "RV1_22" etc

# NOTE: You MUST have atleast 2 replicates per group. Else, DESeq2 will throw
# "Error in base::colMeans(x, na.rm = na.rm, dims = dims, ...) : 'x' must be an
# array of at least two dimensions"

proj <- ""
source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/Custom.Functions.R.Analysis.R")

# DEG.params  <- list(contrast    = c("ARCaPM.NDRG1.4Gy-ARCaPM.NDRG1.0Gy",
#                                     "(RV1_22.NDRG1.4Gy-RV1_22.NDRG1.0Gy)-(RV1_22.WT.4Gy-RV1_22.WT.0Gy)"),
#                     design      = "Comparisons",
#                     design.ref  = c("Cell.Line:ARCaPM", "Condition:WT", "Treatment:0Gy"),
#                     lfc.cutoff  = 0,
#                     padj.cutoff = 0.1,
#                     deseq2.batch.correct = FALSE,
#                     proj        = "RNASeq_Manish_22RV1_ARCaPM",
#                     species     = "Homo sapiens")

method <- "STAR"
data_dir <- paste0("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/", DEG.params$proj, "/")
count_dir <- paste0(data_dir, "counts/", method, "_raw_counts/")
degs_dir <- paste0(data_dir, "DEGs/")
base::dir.create(path = degs_dir)

annotations <- get_annotations()

# Consolidate raw counts from STAR into excel file
read_data <- compile_raw_counts(count_dir, DEG.params$proj, method)

# If raw counts are already consolidated in xlsx, read raw.counts.xlsx
# NOTE: "SYMBOL" column with gene names/ensembl_id/entrez_id without any duplication
read_data <- openxlsx::read.xlsx(xlsxFile = paste0(data_dir, DEG.params$proj, ".raw.counts.xlsx"))

# meta_data MUST be xlsx file and have "Sample.ID" column with sample names
meta_data <- openxlsx::read.xlsx(xlsxFile=paste0(data_dir, DEG.params$proj, ".Metadata.xlsx")) %>%
  dplyr::distinct_at("Sample.ID", .keep_all=TRUE)

# Rename appropiate column in read_data to "SYMBOL"
# read_data <- read_data %>% dplyr::rename(SYMBOL = GeneID)
# colnames(read_data)[1] <- "SYMBOL"

# STOP...STOP...STOP...Subset data based on PCA plot if needed before proceeding

# NOTE: Remove samples that fail QC (eg:PCA plot) or samples that belong to 
# different cell line or in vivo samples from dataset while analyzing in vitro samples
#read_data <- read_data %>% dplyr::select(everything(), -c("SBQuadFc2", "SBQuadFc4"))

#norm_counts_sva(meta_data, read_data, sva.formula_string, annotations, data_dir)

for (n in 1:length(DEG.params$contrast)){

  if (DEG.params$deseq2.batch.correct == TRUE){
    
    # Perform DESeq2() using in-built batch modelling
    # Create DESeq2 object with appropriate variables in design
    deseq2.formula_string <- dplyr::case_when(DEG.params$deseq2.batch.correct == TRUE ~ paste0("~", paste("Batch", DEG.params$design[n], sep = "+")),
                                              TRUE ~ formula_string)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                          colData=meta_data, 
                                          design=~1)
    design(dds) <- as.formula(deseq2.formula_string)
    dds <- run_deseq2(dds, meta_data, DEG.params, n, "DESeq2", degs_dir)
    
    # Perform DESeq2() using sva modelled surrogate variables
    # Create DESeq2 object with surrogate variables in design
    sva.formula_string <- formula_string
    sva_dds <- batch_correct_sva(meta_data, read_data, sva.formula_string)
    sva_dds <- run_deseq2(sva_dds, meta_data, DEG.params, n, "sva", degs_dir)
    
    # Perform DESeq2() using combat corrected counts
    # Create DESeq2 object with appropiate variables in design
    # Get combat corrected raw reads
    read_data_combat <- batch_correct_combat(meta_data, read_data, combat.formula_string)
    combat.formula_string <- formula_string
    combat_dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data_combat,
                                                 colData=meta_data, 
                                                 design=~1)
    design(combat_dds) <- as.formula(combat.formula_string)
    dds <- run_deseq2(combat_dds, meta_data, DEG.params, n, "combat", degs_dir)
  }
}
plotMA_DESeq2(dds, data_dir)

#******************************************************************************#
#                                 VOLCANO PLOT                                 #
#******************************************************************************#

data_dir <- paste0("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/", DEG.params$proj, "/")
volcano_dir <- paste0(data_dir, "Volcano plots/")
base::dir.create(path = volcano_dir)

for (n in 1:length(DEG.params$contrast)){
  
  # Read DEGs
  approach <- ""
  file_suffix <- DEG.params$contrast[n]
  file_suffix <- gsub(pattern="/", replacement="", x=file_suffix)
  df <- read.xlsx(paste0(degs_dir, "DEGs.", file_suffix , "_", approach, ".xlsx"))
  
  # Define volcano plot parameters
  disp_genes <- c()
  output_path <- volcano_dir
  volcano_params <- list(Target           = gsub(pattern="-.*$",replacement="", x=file_suffix),
                         Reference        = gsub(pattern="^.*-",replacement="", x=file_suffix),
                         plot.title       = file_suffix,
                         matrix_color     = "vrds",
                         color_disp_genes = TRUE,
                         label_disp_genes = FALSE,
                         label_top        = FALSE,
                         lfc.cutoff       = 0.58,
                         padj.cutoff      = 0.05)
  
  # Input dataframe MUST have SYMBOL, padj, log2FoldChange columns
  volcano_df <- df %>% 
    dplyr::select(SYMBOL, padj, log2FoldChange)
  
  # Plot the volcano plot
  plot_volcano(volcano_df, volcano_params, disp_genes, file_suffix, output_path)
}

#******************************************************************************#
#                                    HEATMAP                                   #
#******************************************************************************#

data_dir <- paste0("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/", DEG.params$proj, "/")
heatmap_dir <- paste0(data_dir, "Heatmaps/")
base::dir.create(path = heatmap_dir)

# Read normalized counts
norm_counts <- read.xlsx(paste0(data_dir, "Normalized.counts.DESeq2.xlsx"))

# Get list of DEG files
DEG.files <- list.files(degs_dir)
DEG.files <- DEG.files[grepl(pattern="DEGs", x=DEG.files)]
DEG.files <- DEG.files[!grepl(pattern="combat|sva", x=DEG.files)]

# Combine all DEGs from all comparisons
DEGs.list <- c()
for (f in DEG.files){
  df <- read.xlsx(paste0(degs_dir, f)) %>%
    dplyr::filter(padj <= 0.05)
  
  DEGs.list <- c(DEGs.list, df$SYMBOL) %>% unique()
}

# Heatmap parameters
plot_genes <- DEGs.list
disp_genes <- c()
file_suffix <- ""
output_path <- heatmap_dir

# Create metadata for columns
# MUST have column Sample.ID
# MUST have columns listed in heatmap_params$anno.column
metadata_column <- openxlsx::read.xlsx(xlsxFile=paste0(data_dir,DEG.params$proj, ".Metadata.xlsx")) %>%
  dplyr::distinct_at("Sample.ID", .keep_all=TRUE)

# Create metadata for rows
# MUST have column SYMBOL
# MUST have columns listed in heatmap_params$anno.row
metadata_row <- data.frame(SYMBOL = "")
#metadata_row <- data.frame(SYMBOL = "", Gene = "")
#metadata_row <- data.frame(SYMBOL = norm_counts$SYMBOL, Gene = rep(x=c("A", "B"), times=c(700, 799)))

# Plot norm_counts
mat <- norm_counts %>% 
  dplyr::select(everything(), -intersect(colnames(.), c("ENTREZ_ID", "ENTREZ_SYMBOL", "ENSEMBL_ID", "ENSEMBL_SYMBOL"))) %>%
  dplyr::filter(SYMBOL %in% plot_genes) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE)  %>%
  dplyr::filter(!(SYMBOL %in% c("Y_RNA", "Metazoa_SRP")))

plot_heatmap(mat, metadata_column, metadata_row, heatmap.params,
             plot_genes, disp_genes, "Full", output_path)

# Plot heatmap for individual comparisons
for (n in 1:length(DEG.params$contrast)){
  
  # Read DEGs
  approach <- ""
  file_suffix <- DEG.params$contrast[n]
  file_suffix <- gsub(pattern="/", replacement="", x=file_suffix)
  df <- read.xlsx(paste0(degs_dir, "DEGs.", file_suffix , "_", approach, ".xlsx")) %>%
      dplyr::filter(padj <= 0.05)
    
  DEGs.list <- df$SYMBOL %>% unique()
  
  metadata_column <- openxlsx::read.xlsx(xlsxFile=paste0(data_dir, DEG.params$proj, ".Metadata.xlsx")) %>%
    dplyr::distinct_at("Sample.ID", .keep_all=TRUE) %>%
    dplyr::filter(Comparisons %in% c(gsub(pattern="-.*$",replacement="", x=file_suffix),
                                     gsub(pattern="^.*-",replacement="", x=file_suffix)))
  
  metadata_row <- data.frame(SYMBOL = "")
  
  # Heatmap parameters
  plot_genes <- DEGs.list
  disp_genes <- c()
  output_path <- heatmap_dir
  
  mat <- norm_counts %>%
    dplyr::select(everything(), -intersect(colnames(.), c("ENTREZ_ID", "ENTREZ_SYMBOL", "ENSEMBL_ID", "ENSEMBL_SYMBOL"))) %>%
    dplyr::select(SYMBOL, intersect(colnames(.), metadata_column$Sample.ID)) %>%
    dplyr::filter(SYMBOL %in% plot_genes) %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE)  %>%
    dplyr::filter(!(SYMBOL %in% c("Y_RNA", "Metazoa_SRP")))

  plot_heatmap(mat, metadata_column, metadata_row, heatmap.params,
               plot_genes, disp_genes, file_suffix, output_path)

}

#******************************************************************************#
#                               PATHWAY ANALYSIS                               #
#******************************************************************************#

data_dir <- paste0("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/", DEG.params$proj, "/")
pathways_dir <- paste0(data_dir, "Pathways/")
degs_dir <- paste0(data_dir, "DEGs/")
base::dir.create(path = pathways_dir)

# Read the gene set files
species <- DEG.params$species
gmt_dir <- paste0("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets/", species, "/")
gmt_files <- list.files(gmt_dir, full.names = TRUE)

for (n in 1:length(DEG.params$contrast)){
  
  # Read DEGs
  approach <- ""
  file_suffix <- DEG.params$contrast[n]
  file_suffix <- gsub(pattern="/", replacement="", x=file_suffix)
  DEGs_df <- read.xlsx(paste0(degs_dir, "DEGs.", file_suffix , "_", approach, ".xlsx"))
  
  # Input dataframe MUST have SYMBOL, padj, log2FoldChange columns
  DEGs_df <- DEGs_df %>%
    dplyr::select(SYMBOL, padj, log2FoldChange)
  
  # Perform pathway analysis
  pathway_results <- pathway_analysis(DEGs_df, gmt_files)
  
  # Save the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="FGSEA")
  openxlsx::writeData(wb, sheet="FGSEA", x=pathway_results[[1]], rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName="GSEA")
  openxlsx::writeData(wb, sheet="GSEA", x=pathway_results[[2]], rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName="ORA")
  openxlsx::writeData(wb, sheet="ORA", x=pathway_results[[3]], rowNames = FALSE)
  openxlsx::saveWorkbook(wb, file = paste0(pathways_dir, "Pathways.", file_suffix, ".xlsx"),
                         overwrite = TRUE)
}

#******************************************************************************#
#               MASTER REGULATOR & PATHWAY ACTIVITY ANALYSIS                   #
#******************************************************************************#

data_dir <- paste0("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/", DEG.params$proj, "/")
pathways_dir <- paste0(data_dir, "Pathways/")
species <- DEG.params$species

# Read normalized counts
norm_counts <- read.xlsx(paste0(data_dir, "Normalized.counts.DESeq2.xlsx"),
                         sheet = "vst")

get_pathway_tf_scores(norm_counts, species)

#******************************************************************************#
#                             MULTIPLE DE ANALYSIS                             #
#******************************************************************************#

# Sometimes, we want to create a heatmap of samples from multiple conditions.
# In order to see an proper pattern in heatmap, we need to perform DE analysis 
# for every possible pair of conditions, retain most variable genes and then 
# plot the heatmap 

# Perform DESEQ2 on all possible combinations


#******************************************************************************#
#                                 COMPARE ANY 2 DEG RESULTS                    #
#******************************************************************************#

compare_deg_results(df1, df2, file_suffix, output_path)




#*************************PATHWAY COMPARISON ANALYSIS*************************#

# Merge pathways from all differential expression analysis for comparison
# Read Pathway analysis files
Pathway.files <- list.files(pathways_dir)
Pathway.files <- Pathway.files[grepl(pattern="Pathway.", x=Pathway.files)]

file_suffix <- gsub(pattern="Pathway.Analysis|.xlsx",replacement="", x=Pathway.files[1])
ora.pathway.df <- read.xlsx(paste0(data_dir, Pathway.files[1]), sheet = "ORA") %>%
  dplyr::filter(p.adjust <= 0.05) %>%
  dplyr::select(Description, k.K) %>%
  dplyr::rename(!!rlang::sym(file_suffix):= k.K)

for (f in Pathway.files[-1]){
  
  file_suffix <- gsub(pattern="Pathway.Analysis|.xlsx",replacement="", x=f)
  
  df <- read.xlsx(paste0(data_dir, f)) %>%
    dplyr::filter(p.adjust <= 0.05) %>%
    dplyr::select(Description, k.K) %>%
    dplyr::rename(!!rlang::sym(file_suffix):= k.K)
  
  ora.pathway.df <- ora.pathway.df %>% 
    dplyr::full_join(df, by=c("Description"="Description"))
}

fgsea.pathway.df <- read.xlsx(paste0(data_dir, Pathway.files[1]), sheet = "FGSEA") %>%
  dplyr::filter(padj <= 0.05) %>%
  dplyr::select(pathway, NES) %>%
  dplyr::rename(!!rlang::sym(file_suffix):= NES)

for (f in Pathway.files[-1]){
  
  file_suffix <- gsub(pattern="Pathway.Analysis|.xlsx",replacement="", x=f)
  
  df <- read.xlsx(paste0(data_dir, f), sheet = "FGSEA") %>%
    dplyr::filter(padj <= 0.05) %>%
    dplyr::select(pathway, NES) %>%
    dplyr::rename(!!rlang::sym(file_suffix):= NES)
  
  fgsea.pathway.df <- fgsea.pathway.df %>% 
    dplyr::full_join(df, by=c("pathway"="pathway"))
}

# Create workbook to store results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName="ORA")
openxlsx::writeData(wb, sheet="ORA", x=ora.pathway.df, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName="FGSEA")
openxlsx::writeData(wb, sheet="FGSEA", x=fgsea.pathway.df, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(data_dir, "Pathway.Analysis.Comparison.xlsx"),
                       overwrite = TRUE)


































#***********************************TESTING************************************#

DEG.params  <- list(Variable    = c("Comparisons"),
                    Target      = c("ARCaPM.4Gy.NDRG1_mut"),
                    Reference   = c("ARCaPM.0Gy.NDRG1_mut"),
                    contrast    = c("ARCaPM.NDRG1_mut.4Gy-ARCaPM.NDRG1_mut.0Gy"),
                    lfc.cutoff  = 0,
                    padj.cutoff = 0.1,
                    design      = "Cell.Line*Treatment*Condition", #"Comparisons", #"0+Condition:Treatment"
                    design.ref  = c("Condition:WT", "Treatment:0Gy", "Comparisons:ARCaPM.0Gy.WT"),
                    deseq2.batch.correct = FALSE,
                    proj        = "RNASeq_Manish_22RV1_ARCaPM",
                    species     = "Homo sapiens")

meta_data <- meta_data %>% filter(Cell.Line != "22RV1")
n <- 1


dds.new <- run_deseq2(dds, meta_data, DEG.params, n, approach, data_dir)
design(dds) <- as.formula(~Comparisons)
dds.old <- run_deseq2_old(dds, meta_data, DEG.params, n, approach, data_dir)


# The following 3 approaches give identical results

# Method 1 => Similar way of defining contrasts like method2. Easy to compare 
# samples but difference of difference not possible
design <- "0+Condition:Treatment" 
dds.test <- dds
contrast1 <- c("Condition", "NDRG1_mut.Treatment4Gy", "NDRG1_mut.Treatment0Gy")
contrast2 <- c("Condition", "NDRG1_mut.Treatment4Gy", "WT.Treatment0Gy")
contrast3 <- c("Condition", "NDRG1_mut.Treatment4Gy", "WT.Treatment4Gy")
contrast4 <- c("Condition", "NDRG1_mut.Treatment0Gy", "WT.Treatment4Gy")
contrast5 <- c("Condition", "NDRG1_mut.Treatment0Gy", "WT.Treatment0Gy")
contrast6 <- c("Condition", "WT.Treatment4Gy", "WT.Treatment0Gy")
res1 <- DESeq2::results(dds, contrast = contrast1)
res2 <- DESeq2::results(dds, contrast = contrast2)
res3 <- DESeq2::results(dds, contrast = contrast3)
res4 <- DESeq2::results(dds, contrast = contrast4)
res5 <- DESeq2::results(dds, contrast = contrast5)
res6 <- DESeq2::results(dds, contrast = contrast6)

# Method 2= > combine COndition and Treatment to new column Comparisons. 
# Similar way of defining contrasts like method1. Easy to compare 
# samples but difference of difference not possible
design <- "Comparisons"
contrast1A <- c("Comparisons", "ARCaPM.4Gy.NDRG1_mut", "ARCaPM.0Gy.NDRG1_mut")
contrast2A <- c("Comparisons", "ARCaPM.4Gy.NDRG1_mut", "ARCaPM.0Gy.WT")
contrast3A <- c("Comparisons", "ARCaPM.4Gy.NDRG1_mut", "ARCaPM.4Gy.WT")
contrast4A <- c("Comparisons", "ARCaPM.0Gy.NDRG1_mut", "ARCaPM.4Gy.WT")
contrast5A <- c("Comparisons", "ARCaPM.0Gy.NDRG1_mut", "ARCaPM.0Gy.WT")
contrast6A <- c("Comparisons", "ARCaPM.4Gy.WT", "ARCaPM.0Gy.WT")
res1A <- DESeq2::results(dds, contrast = contrast1A)
res2A <- DESeq2::results(dds, contrast = contrast2A)
res3A <- DESeq2::results(dds, contrast = contrast3A)
res4A <- DESeq2::results(dds, contrast = contrast4A)
res5A <- DESeq2::results(dds, contrast = contrast5A)
res6A <- DESeq2::results(dds, contrast = contrast6A)

# Method 3 => conventional design, contrast needs to calculated for each sample
# difference of difference is easy
design <- "Condition+Treatment+Condition:Treatment"
dds.standard <- dds
mod_mat <- model.matrix(design(dds), colData(dds))
NDRG1_0Gy <- colMeans(mod_mat[dds$Condition == "NDRG1_mut" & dds$Treatment == "0Gy",])
NDRG1_4Gy <- colMeans(mod_mat[dds$Condition == "NDRG1_mut" & dds$Treatment == "4Gy",])
WT_0Gy <- colMeans(mod_mat[dds$Condition == "WT" & dds$Treatment == "0Gy",])
WT_4Gy <- colMeans(mod_mat[dds$Condition == "WT" & dds$Treatment == "4Gy",])

contrast1B <- NDRG1_4Gy - NDRG1_0Gy 
contrast2B <- NDRG1_4Gy - WT_0Gy
contrast3B <- NDRG1_4Gy - WT_4Gy
contrast4B <- NDRG1_0Gy - WT_4Gy 
contrast5B <- NDRG1_0Gy - WT_0Gy
contrast6B <- WT_4Gy - WT_0Gy
contrast7B <- (NDRG1_4Gy - NDRG1_0Gy) - (WT_4Gy - WT_0Gy)
contrast8B <- 
  
res1B. <- DESeq2::results(dds, contrast = contrast1B)
res2B. <- DESeq2::results(dds, contrast = contrast2B)
res3B. <- DESeq2::results(dds, contrast = contrast3B)
res4B. <- DESeq2::results(dds, contrast = contrast4B)
res5B. <- DESeq2::results(dds, contrast = contrast5B)
res6B. <- DESeq2::results(dds, contrast = contrast6B)
res7B. <- DESeq2::results(dds, contrast = contrast7B)


# Check for differences in results
df1 <- data.frame(res1) %>% dplyr::mutate_all(~(round(.,2)))
df1A <- data.frame(res1A) %>% dplyr::mutate_all(~(round(.,2)))
df1B <- data.frame(res1B) %>% dplyr::mutate_all(~(round(.,2)))

df2 <- data.frame(res2) %>% dplyr::mutate_all(~(round(.,2)))
df2A <- data.frame(res2A) %>% dplyr::mutate_all(~(round(.,2)))
df2B <- data.frame(res2B) %>% dplyr::mutate_all(~(round(.,2)))

df3 <- data.frame(res3) %>% dplyr::mutate_all(~(round(.,2)))
df3A <- data.frame(res3A) %>% dplyr::mutate_all(~(round(.,2)))
df3B <- data.frame(res3B) %>% dplyr::mutate_all(~(round(.,2)))

df4 <- data.frame(res4) %>% dplyr::mutate_all(~(round(.,2)))
df4A <- data.frame(res4A) %>% dplyr::mutate_all(~(round(.,2)))
df4B <- data.frame(res4B) %>% dplyr::mutate_all(~(round(.,2)))

df5 <- data.frame(res5) %>% dplyr::mutate_all(~(round(.,2)))
df5A <- data.frame(res5A) %>% dplyr::mutate_all(~(round(.,2)))
df5B <- data.frame(res5B) %>% dplyr::mutate_all(~(round(.,2)))

df6 <- data.frame(res6) %>% dplyr::mutate_all(~(round(.,2)))
df6A <- data.frame(res6A) %>% dplyr::mutate_all(~(round(.,2)))
df6B <- data.frame(res6B) %>% dplyr::mutate_all(~(round(.,2)))

# Display only rows and columns that are different
dplyr::full_join(df1[rowSums(df1-df1A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
                 df1A[rowSums(df1-df1A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
                 by=c("SYMBOL"="SYMBOL")) %>%
  dplyr::transmute(SYMBOL              = SYMBOL,
                   baseMean.diff       = baseMean.x-baseMean.y,
                   log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
                   lfcSE.diff          = lfcSE.x-lfcSE.y,
                   stat.diff           = stat.x-stat.y,
                   pvalue.diff         = pvalue.x-pvalue.y,
                   padj.diff           = padj.x-padj.y)

dplyr::full_join(df2[rowSums(df2-df2A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
                 df2A[rowSums(df2-df2A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
                 by=c("SYMBOL"="SYMBOL")) %>%
  dplyr::transmute(SYMBOL              = SYMBOL,
                   baseMean.diff       = baseMean.x-baseMean.y,
                   log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
                   lfcSE.diff          = lfcSE.x-lfcSE.y,
                   stat.diff           = stat.x-stat.y,
                   pvalue.diff         = pvalue.x-pvalue.y,
                   padj.diff           = padj.x-padj.y)

dplyr::full_join(df3[rowSums(df3-df3A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
                 df3A[rowSums(df3-df3A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
                 by=c("SYMBOL"="SYMBOL")) %>%
  dplyr::transmute(SYMBOL              = SYMBOL,
                   baseMean.diff       = baseMean.x-baseMean.y,
                   log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
                   lfcSE.diff          = lfcSE.x-lfcSE.y,
                   stat.diff           = stat.x-stat.y,
                   pvalue.diff         = pvalue.x-pvalue.y,
                   padj.diff           = padj.x-padj.y)

dplyr::full_join(df4[rowSums(df4-df4A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
                 df4A[rowSums(df4-df4A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
                 by=c("SYMBOL"="SYMBOL")) %>%
  dplyr::transmute(SYMBOL              = SYMBOL,
                   baseMean.diff       = baseMean.x-baseMean.y,
                   log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
                   lfcSE.diff          = lfcSE.x-lfcSE.y,
                   stat.diff           = stat.x-stat.y,
                   pvalue.diff         = pvalue.x-pvalue.y,
                   padj.diff           = padj.x-padj.y)

dplyr::full_join(df5[rowSums(df5-df5A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
                 df5A[rowSums(df5-df5A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
                 by=c("SYMBOL"="SYMBOL")) %>%
  dplyr::transmute(SYMBOL              = SYMBOL,
                   baseMean.diff       = baseMean.x-baseMean.y,
                   log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
                   lfcSE.diff          = lfcSE.x-lfcSE.y,
                   stat.diff           = stat.x-stat.y,
                   pvalue.diff         = pvalue.x-pvalue.y,
                   padj.diff           = padj.x-padj.y)

dplyr::full_join(df6[rowSums(df6-df6A, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
                 df6A[rowSums(df6-df6A, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
                 by=c("SYMBOL"="SYMBOL")) %>%
  dplyr::transmute(SYMBOL              = SYMBOL,
                   baseMean.diff       = baseMean.x-baseMean.y,
                   log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
                   lfcSE.diff          = lfcSE.x-lfcSE.y,
                   stat.diff           = stat.x-stat.y,
                   pvalue.diff         = pvalue.x-pvalue.y,
                   padj.diff           = padj.x-padj.y)

dplyr::full_join(df1[rowSums(df1-df1B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
                 df1B[rowSums(df1-df1B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
                 by=c("SYMBOL"="SYMBOL")) %>%
  dplyr::transmute(SYMBOL              = SYMBOL,
                   baseMean.diff       = baseMean.x-baseMean.y,
                   log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
                   lfcSE.diff          = lfcSE.x-lfcSE.y,
                   stat.diff           = stat.x-stat.y,
                   pvalue.diff         = pvalue.x-pvalue.y,
                   padj.diff           = padj.x-padj.y)

dplyr::full_join(df2[rowSums(df2-df2B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
                 df2B[rowSums(df2-df2B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
                 by=c("SYMBOL"="SYMBOL")) %>%
  dplyr::transmute(SYMBOL              = SYMBOL,
                   baseMean.diff       = baseMean.x-baseMean.y,
                   log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
                   lfcSE.diff          = lfcSE.x-lfcSE.y,
                   stat.diff           = stat.x-stat.y,
                   pvalue.diff         = pvalue.x-pvalue.y,
                   padj.diff           = padj.x-padj.y)

dplyr::full_join(df3[rowSums(df3-df3B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
                 df3B[rowSums(df3-df3B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
                 by=c("SYMBOL"="SYMBOL")) %>%
  dplyr::transmute(SYMBOL              = SYMBOL,
                   baseMean.diff       = baseMean.x-baseMean.y,
                   log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
                   lfcSE.diff          = lfcSE.x-lfcSE.y,
                   stat.diff           = stat.x-stat.y,
                   pvalue.diff         = pvalue.x-pvalue.y,
                   padj.diff           = padj.x-padj.y)

dplyr::full_join(df4[rowSums(df4-df4B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
                 df4B[rowSums(df4-df4B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
                 by=c("SYMBOL"="SYMBOL")) %>%
  dplyr::transmute(SYMBOL              = SYMBOL,
                   baseMean.diff       = baseMean.x-baseMean.y,
                   log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
                   lfcSE.diff          = lfcSE.x-lfcSE.y,
                   stat.diff           = stat.x-stat.y,
                   pvalue.diff         = pvalue.x-pvalue.y,
                   padj.diff           = padj.x-padj.y)

dplyr::full_join(df5[rowSums(df5-df5B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
                 df5B[rowSums(df5-df5B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
                 by=c("SYMBOL"="SYMBOL")) %>%
  dplyr::transmute(SYMBOL              = SYMBOL,
                   baseMean.diff       = baseMean.x-baseMean.y,
                   log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
                   lfcSE.diff          = lfcSE.x-lfcSE.y,
                   stat.diff           = stat.x-stat.y,
                   pvalue.diff         = pvalue.x-pvalue.y,
                   padj.diff           = padj.x-padj.y)

dplyr::full_join(df6[rowSums(df6-df6B, na.rm=TRUE) !=0, ] %>%tibble::rownames_to_column("SYMBOL"), 
                 df6B[rowSums(df6-df6B, na.rm=TRUE) !=0, ] %>% tibble::rownames_to_column("SYMBOL"),
                 by=c("SYMBOL"="SYMBOL")) %>%
  dplyr::transmute(SYMBOL              = SYMBOL,
                   baseMean.diff       = baseMean.x-baseMean.y,
                   log2FoldChange.diff = log2FoldChange.x-log2FoldChange.y,
                   lfcSE.diff          = lfcSE.x-lfcSE.y,
                   stat.diff           = stat.x-stat.y,
                   pvalue.diff         = pvalue.x-pvalue.y,
                   padj.diff           = padj.x-padj.y)



res1 <- DESeq2::lfcShrink(dds=dds.test, res=res1, type="ashr")
res1B <- DESeq2::lfcShrink(dds=dds.standard, res=res1B, type="ashr")
res2 <- DESeq2::lfcShrink(dds=dds.test, res=res2, type="ashr")
res2B <- DESeq2::lfcShrink(dds=dds.standard, res=res2B, type="ashr")
res3 <- DESeq2::lfcShrink(dds=dds.test, res=res3, type="ashr")
res3B <- DESeq2::lfcShrink(dds=dds.standard, res=res3B, type="ashr")
res4 <- DESeq2::lfcShrink(dds=dds.test, res=res4, type="ashr")
res4B <- DESeq2::lfcShrink(dds=dds.standard, res=res4B, type="ashr")
res5 <- DESeq2::lfcShrink(dds=dds.test, res=res5, type="ashr")
res5B <- DESeq2::lfcShrink(dds=dds.standard, res=res5B, type="ashr")
res6 <- DESeq2::lfcShrink(dds=dds.test, res=res6, type="ashr")
res6B <- DESeq2::lfcShrink(dds=dds.standard, res=res6B, type="ashr")
res7B <- DESeq2::lfcShrink(dds=dds.standard, res=res7B, type="ashr")


summary(res1)
summary(res1B)
summary(res1B.)
summary(res2)
summary(res2B)
summary(res2B.)
summary(res3)
summary(res3B)
summary(res3B.)
summary(res4)
summary(res4B)
summary(res4B.)
summary(res5)
summary(res5B)
summary(res5B.)
summary(res6)
summary(res6B)
summary(res6B.)
summary(res7B)
summary(res7B.)
res1B. <- DESeq2::lfcShrink(dds=dds, res=res1B., type="ashr")
res2B. <- DESeq2::lfcShrink(dds=dds, res=res2B., type="ashr")
res3B. <- DESeq2::lfcShrink(dds=dds, res=res3B., type="ashr")
res4B. <- DESeq2::lfcShrink(dds=dds, res=res4B., type="ashr")
res5B. <- DESeq2::lfcShrink(dds=dds, res=res5B., type="ashr")
res6B. <- DESeq2::lfcShrink(dds=dds, res=res6B., type="ashr")
res7B. <- DESeq2::lfcShrink(dds=dds, res=res7B., type="ashr")

setdiff(rownames(data.frame(res1B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
        rownames(data.frame(res1B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 

setdiff(rownames(data.frame(res1B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
        rownames(data.frame(res1B.) %>% filter(log2FoldChange < -0.58, padj < 0.1)))

setdiff(rownames(data.frame(res2B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
        rownames(data.frame(res2B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 

setdiff(rownames(data.frame(res2B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
        rownames(data.frame(res2B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 


setdiff(rownames(data.frame(res3B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
        rownames(data.frame(res3B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 

setdiff(rownames(data.frame(res3B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
        rownames(data.frame(res3B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 


setdiff(rownames(data.frame(res4B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
        rownames(data.frame(res4B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 

setdiff(rownames(data.frame(res4B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
        rownames(data.frame(res4B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 

setdiff(rownames(data.frame(res5B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
        rownames(data.frame(res5B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 

setdiff(rownames(data.frame(res5B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
        rownames(data.frame(res5B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 


setdiff(rownames(data.frame(res6B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
        rownames(data.frame(res6B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 

setdiff(rownames(data.frame(res6B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
        rownames(data.frame(res6B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 

setdiff(rownames(data.frame(res7B) %>% filter(log2FoldChange > 0.58, padj < 0.1)),
        rownames(data.frame(res7B.) %>% filter(log2FoldChange > 0.58, padj < 0.1))) 

setdiff(rownames(data.frame(res7B) %>% filter(log2FoldChange < -0.58, padj < 0.1)),
        rownames(data.frame(res7B.) %>% filter(log2FoldChange < -0.58, padj < 0.1))) 