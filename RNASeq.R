# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
# https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html

#******************************************************************************#
#                                  RUN DESEQ2                                  #
#******************************************************************************#

# This function runs DESeq2 and plots volcano plot, heatmap, PCA plot and
# correlation plot.
# meta_data MUST have "Sample" column with sample names without any duplication
# meta_data MUST have "Batch" column with batch info
# read_data MUST have "SYMBOL" column with gene names without any duplication
# Several global variables need to be defined. See RNASeq.R

analyze_DESeq2 <- function(meta_data, read_data, f_suffix){
  
  #****************************************************************************#
  #                           STEP 1: IMPORT META DATA                         #
  #****************************************************************************#
  
  # SUBSET METADATA (VERY VERY IMPORTANT, OFTEN OVERLOOKED ASPECT OF ANALYSIS)
  # DESeq2 calculates sizefactors using all the samples in meta_data. So, final
  # results are heavily dependent on sizefactors. For instance, if you have 
  # 5 replicates of Vehicle treated samples, 5 replicates of Enzalutamide 
  # treated samples and 5 replicates of R1881 treated samples, and you want to 
  # find DEGs between R1881 vs vehicle, remove the Enzalutamide treated samples 
  # from meta_data & read_data. Since, we have to apply batch correction on the
  # whole dataset, we keep all samples for now but we remove them before 
  # differential expression analysis.
  # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#if-i-have-multiple-groups-should-i-run-all-together-or-split-into-pairs-of-groups
  
  meta_data <- meta_data %>%
    dplyr::filter(!is.na(Sample)) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames(var="Sample") %>%
    dplyr::mutate(id = dplyr::if_else(get(Variable) %in% c(Comparisons$Target, Comparisons$Reference),
                                      get(Variable), "ignore")) %>%
    dplyr::mutate(id = dplyr::if_else(is.na(id), "ignore", id)) %>%
    dplyr::filter(id != "ignore")
    
  
  #****************************************************************************#
  #                          STEP 2: IMPORT READ DATA                          #                                           
  #****************************************************************************#
  
  #colnames(read_data) <- gsub("_.*", "", colnames(read_data))
  read_data <- read_data %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames(var="SYMBOL") %>%
    dplyr::select(intersect(colnames(.), rownames(meta_data))) %>%
    replace(is.na(.), 0)
  
  #****************************************************************************#
  #              STEP 3: PREPARE META DATA and READ DATA FOR DESEQ2            #
  #****************************************************************************#
  
  # DESeq2 automatically removes genes that have 0 counts. However, we have to
  # make sure there are no samples which have 0 counts for all genes. For eg,
  # sample MB1 has 0 counts for all genes in cluster 35. So, we need to manually
  # remove this sample before DESeq2 analysis. Else, the geometric mean will be 
  # zero for all genes and DESeq2 will halt execution.
  
  if (TRUE %in% (colSums(read_data)==0)){
    read_data <- read_data[, -c(which(colSums(read_data)==0))]
    meta_data <- meta_data[-c(which(colSums(read_data)==0)),]
  }
  
  # Rearrange read_data to match meta_data
  read_data <- read_data[,rownames(meta_data)]
  
  # Remove samples that have "NA" for variable being analyzed in metadata.
  NA_sample_list <- which(is.na(meta_data[Variable]))
  if (length(NA_sample_list)!=0){
    meta_data <- meta_data[-c(NA_sample_list),]
    read_data <- read_data[,-c(NA_sample_list)]
  }
  
  # Remove any column named "sizeFactor" in meta_data as it interferes with DESeq2()
  if ("sizeFactor" %in% colnames(meta_data)){
    meta_data <- meta_data %>% dplyr::select(!sizeFactor)
  }
  
  # Check if all columns of meta table are factors. Else, convert to factors.
  print("meta_data before factor conversion")
  str(meta_data)
  for (i in 1:ncol(meta_data)){
    meta_data[,i] <- as.factor(meta_data[,i])
  }
  print("meta_data after factor conversion")
  str(meta_data)
  
  # Check (i) read_data & meta_data are dataframes
  # (ii) column names of read_data (i.e.sample name) are PRESENT in row names of meta_data (i.e. sample name)
  # (iii) column names of read_data (i.e.sample name) are PRESENT in same order as row names of meta_data (i.e. sample name)
  if (class(read_data) == "data.frame" & class(meta_data) == "data.frame" & 
      base::all(colnames(read_data) %in% rownames(meta_data)) &
      base::all(colnames(read_data) == rownames(meta_data))){
    print("ALL OK")
  } else if (class(read_data) != "data.frame"){
    print("read_data is NOT a dataframe")
  } else if (class(meta_data) != "data.frame"){
    print("meta_data is NOT a dataframe")
  } else if (!(all(colnames(read_data) %in% rownames(meta_data)))){
    print("Missing samples")
  } else if (all(colnames(read_data) != rownames(meta_data))){
    print("Varying sample order")
  } else{}
  
  #****************************************************************************#
  #                      STEP 4: PERFORM BATCH CORRECTION                      #
  #****************************************************************************#
  
  if (length(unique(as.numeric(meta_data$Batch))) > 1){
    # Adjust for batch effects using Combat-Seq
    corrected_read_data <- sva::ComBat_seq(counts = as.matrix(read_data), 
                                           batch = as.numeric(meta_data$Batch), 
                                           group = as.numeric(as.factor(meta_data$id)))
  } else{
    corrected_read_data <- read_data
  }
  
  #****************************************************************************#
  #                     STEP 5: CALCULATE NORMALIZED COUNTS                    #
  #****************************************************************************#
  
  # NOTE: sizeFactors MUST be estimated first before normalization.
  # sizeFactors will differ based on samples included i.e. all samples vs subset
  # of samples. This will in turn affect normalized counts.
  # sizeFactors doesnt depend on the design formula. So, design ~1 or design ~id
  # will give same normalized values.
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = corrected_read_data, colData = meta_data, design = ~ 1)
  dds <- DESeq2::estimateSizeFactors(dds) 
  
  normalized_counts <- DESeq2::counts(dds, normalized=TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID")
  
  ensembl <- length(intersect(normalized_counts$ID, annotations$ENSEMBL_ID))
  entrez <- length(intersect(normalized_counts$ID, annotations$ENTREZ_ID))
  symbol <- length(intersect(normalized_counts$ID, annotations$SYMBOL))
  
  if (ensembl > entrez & ensembl > symbol){
    normalized_counts <- annotations %>%
      dplyr::right_join(normalized_counts, by = c("ENSEMBL_ID"="ID"), multiple = "all") %>%
      dplyr::select(SYMBOL, everything(), -c(CHR, DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
      dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                              is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                              TRUE ~ SYMBOL)) %>%
      dplyr::distinct_at("ENSEMBL_ID", .keep_all = TRUE)
  }
  
  if (entrez > ensembl & entrez > symbol){
    normalized_counts <- annotations %>%
      dplyr::right_join(normalized_counts, by = c("ENTREZ_ID"="ID"), multiple = "all") %>%
      dplyr::select(SYMBOL, everything(), -c(CHR, DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
      dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                              is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                              TRUE ~ SYMBOL)) %>%
      dplyr::distinct_at("ENTREZ_ID", .keep_all = TRUE)
  }
  
  if (symbol > ensembl & symbol > entrez){
    normalized_counts <- annotations %>%
      dplyr::right_join(normalized_counts, by = c("SYMBOL"="ID"), multiple = "all") %>%
      dplyr::select(SYMBOL, everything(), -c(CHR, DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
      dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                              is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                              TRUE ~ SYMBOL)) %>%
      dplyr::distinct_at("SYMBOL", .keep_all = TRUE)
  }
  
  # Save batch corrected normalized counts for entire dataset
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "normalized")
  openxlsx::writeData(wb, sheet = "normalized", x = normalized_counts)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "Normalized_Counts.xlsx"),
                         overwrite = TRUE)
  
  #****************************************************************************#
  #                          STEP 6: PREPARE QC PLOTS                          #
  #****************************************************************************#
  
  # Plot PCA
  vst_mat <- assay(DESeq2::vst(dds, blind=FALSE))
  pca <- prcomp(t(vst_mat))
  df <- cbind(meta_data, pca$x)      # Create data frame with metadata, PC1 & PC2 values for ggplot
  ggplot(df, aes(x=PC1, y=PC2, color = id)) + 
    geom_point() +
    geom_text_repel(label = rownames(df))
  ggplot2::ggsave(filename =  paste0("PCA_Plot_overall.pdf"),
                  plot = last_plot(),
                  device = "pdf",
                  path = results_path,
                  scale = 1,
                  #width = 8.5,
                  #height = 11,
                  units = c("in"),	 
                  dpi = 300,
                  limitsize = TRUE,
                  bg = "white")
  
  
  #****************************************************************************#
  #                           STEP 7: PERFORM DESEQ2                           #
  #****************************************************************************#
  
  meta <- meta_data
  read <- corrected_read_data
  
  # Perform DEG analysis for each comparison with the dataset
  for (n in 1:length(Comparisons$Target)){
    
    # Subset meta_data
    meta_data <- meta %>% dplyr::filter(get(Variable) %in% c(Comparisons$Target[n], Comparisons$Reference[n]))
    
    # Subset read_data
    corrected_read_data <- read %>% dplyr::select(rownames(meta_data))
    
    #**************************************************************************#
    #                          STEP 5: PERFORM DESEQ2                          #
    #**************************************************************************#
    
    # Create DESeq2Dataset object, estimate Size factors and view it
    # NOTE: sizefactors MUST be calculated ONLY using samples being compared for 
    # differential expression. So, make sure read_data and meta_data ONLY have
    # samples being compared for differential expression.
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = corrected_read_data, colData = meta_data, design = ~ id)
    # dds <- DESeq2::DESeqDataSetFromMatrix(countData = corrected_read_data, colData = meta_data, design = ~ Tissue + Age)
    # dds <- DESeq(dds, test="LRT", reduced = ~ Tissue)
    # res <- DESeq2::results(object = dds)
    
    # Study effect of Treatment ignoring effect of Sex
    #dds <- DESeqDataSetFromMatrix(countData = corrected_read_data, colData = meta_data, design = ~ Treatment)
    
    # Study effect of Treatment after removing effect of Sex (assumes equal effect of Sex on different treatments)
    #dds <- DESeqDataSetFromMatrix(countData = corrected_read_data, colData = meta_data, design = ~ Sex + Treatment)
    
    # Study effect of Sex on treatment?
    # Study effect of Treatment after removing effect of Sex (assumes different effect of Sex on different treatments)
    #dds <- DESeqDataSetFromMatrix(countData = corrected_read_data, colData = meta_data, design = ~ Sex + Treatment + Sex:Treatment)
    
    # Run analysis using both local and parameteric fit
    # betaPrior: default=FALSE, shrunken LFCs are obtained later using lfcShrink
    # fitType: "local" may be better sometimes
    dds_para <- DESeq2::DESeq(object = dds, 
                              test = "Wald",
                              fitType = "parametric",
                              betaPrior = FALSE,
                              minReplicatesForReplace = 7)
    dds_local <- DESeq2::DESeq(object = dds, 
                               test = "Wald",
                               fitType = "local",
                               betaPrior = FALSE,
                               minReplicatesForReplace = 7)
    
    # Determine whether to use "local" or "parametric" fit
    residual_para <- mcols(dds_para)$dispGeneEst - mcols(dds_para)$dispFit
    residual_local <- mcols(dds_local)$dispGeneEst - mcols(dds_local)$dispFit
    
    if (median(residual_para^2, na.rm=TRUE) <= median(residual_local^2, na.rm=TRUE)){
      fit <- "parametric"
      dds <- dds_para
    } else{
      fit <- "local"
      dds <- dds_local
    }
    
    # Define contrast and coeff
    DE_levels <- as.vector(unique(meta_data$id))
    contrast <- c("id", 
                  DE_levels[grepl(pattern = Comparisons$Target[n], x=DE_levels)],
                  DE_levels[grepl(pattern = Comparisons$Reference[n], x=DE_levels)])
    #contrast <- contrast[-4]
    
    coeff <- paste(contrast[1], contrast[2], "vs", contrast[3], sep='_')
    cat("Coeff is ", coeff, "\n")
    
    # Calculate results
    res <- DESeq2::results(object = dds,
                           contrast = contrast,
                           lfcThreshold = lfc.cutoff,
                           altHypothesis = "greaterAbs",
                           cooksCutoff = TRUE,
                           independentFiltering = TRUE,
                           alpha = padj.cutoff,
                           pAdjustMethod = "BH")
    
    # Perform lfcshrinkage to account for variability between replicates
    # For ashr, if res is provided, then coef and contrast are ignored.
    # lfcshrinkage will not change the number of DEGs and affects only logFC
    res <- DESeq2::lfcShrink(dds = dds, res = res, type = "ashr")
    
    # Summarize results
    summary(res)
    
    # Export results as dataframe and replace ensembl IDs with Gene names
    results <- res %>%
      as.data.frame() %>%
      tibble::rownames_to_column("ID") 
    
    ensembl <- length(intersect(normalized_counts$ID, annotations$ENSEMBL_ID))
    entrez <- length(intersect(normalized_counts$ID, annotations$ENTREZ_ID))
    symbol <- length(intersect(normalized_counts$ID, annotations$SYMBOL))
    
    if (ensembl > entrez & ensembl > symbol){
      results <- annotations %>%
        dplyr::right_join(results, by = c("ENSEMBL_ID"="ID"), multiple = "all") %>%
        dplyr::select(SYMBOL, CHR, everything(), -c(DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
        dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                                is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                                TRUE ~ SYMBOL)) %>%
        dplyr::distinct_at("ENSEMBL_ID", .keep_all = TRUE)
    }
    
    if (entrez > ensembl & entrez > symbol){
      results <- annotations %>%
        dplyr::right_join(results, by = c("ENTREZ_ID"="ID"), multiple = "all") %>%
        dplyr::select(SYMBOL, CHR, everything(), -c(DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
        dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                                is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                                TRUE ~ SYMBOL)) %>%
        dplyr::distinct_at("ENTREZ_ID", .keep_all = TRUE)
    }
    
    if (symbol > ensembl & symbol > entrez){
      results <- annotations %>%
        dplyr::right_join(results, by = c("SYMBOL"="ID"), multiple = "all") %>%
        dplyr::select(SYMBOL, CHR, everything(), -c(DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
        dplyr::mutate(SYMBOL = dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                                is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                                TRUE ~ SYMBOL)) %>%
        dplyr::distinct_at("SYMBOL", .keep_all = TRUE)
    }
    
    # Define any filename you want added to final file
    file_suffix <- paste0(Variable2_value, "_", coeff)
    file_suffix <- gsub(pattern="/", replacement="", x=file_suffix)
    
    # Save the results
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "DEGs")
    openxlsx::writeData(wb, sheet = "DEGs", x = results, rowNames = FALSE)
    openxlsx::saveWorkbook(wb, 
                           file = paste0(results_path, "Results_", file_suffix, "_DEGs2_", f_suffix, ".xlsx"), 
                           overwrite = TRUE)
    
    #**************************************************************************#
    #             STEP 6: CREATE HEATMAPS, VOLCANO PLOTS, PCA PLOTS            #
    #**************************************************************************#
    
    if (heatmap_plot == TRUE){
      
      # Refer https://www.biostars.org/p/339065/ where the authors of DESeq2 recommend
      # using vst or rld rather than log1p transformation while plotting DESeq2 
      # normalized counts in heatmap.
      # NOTE: We scale the vst values and replace all NAs with 0 before feeding it to 
      # pheatmap(). If we perform the scaling in pheatmap, the NAs generated during
      # scaling will give error.
      # NOTE: While use of log1p transformation is NOT recommended, I notice only 
      # minor changes in color
      
      # Define any filename you want added to final file
      file_suffix <- paste0(Variable2_value, "_", coeff)
      file_suffix <- gsub(pattern="/", replacement="", x=file_suffix)
      file_suffix <- paste0(file_suffix, "_", f_suffix)
      
      # (I) Read metadata
      metadata <- meta_data %>%
        tibble::rownames_to_column("Sample")
      
      # (II) Define genes to plot
      plot_genes <- results %>%
        dplyr::filter(padj < 0.05) %>% 
        dplyr::select(SYMBOL) %>%
        unlist(use.names = FALSE)
      
      # (III) Genes to display in heatmap
      disp_genes <- c()
      
      # (IV) Read expr data
      #normalized_counts <- normalized_counts
      normalized_counts <- assay(DESeq2::vst(dds, blind=FALSE)) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("ENSEMBL_ID") %>%
        dplyr::left_join(annotations, by = c("ENSEMBL_ID"="ENSEMBL_ID"), multiple = "all") %>%
        dplyr::filter(nchar(SYMBOL) > 0) %>%
        dplyr::select(SYMBOL, everything(), -c(ID, CHR, DESCRIPTION, START, END, STRAND))
      colnames(normalized_counts)[1] <- "SYMBOL"
      
      # Run the function
      plot_heatmap(normalized_counts, metadata, plot_genes, disp_genes, file_suffix)
    }
    
    if (volcano_plot ==TRUE){
      
      # Define any filename you want added to final file
      file_suffix <- paste0(Variable2_value, "_", coeff)
      file_suffix <- gsub(pattern="/", replacement="", x=file_suffix)
      file_suffix <- paste0(file_suffix, "_", f_suffix)
      
      # (I) Import expression data with log2FC and pval
      volcano_df <- results %>%
        dplyr::rename(log2FC = log2FoldChange, Gene = SYMBOL)
      
      # (II) Define or import metadata
      # NOTE: Metadata is a dataframe with "Sample" and "Condition" columns
      metadata <- meta_data %>%
        tibble::rownames_to_column("Sample")
      
      # (III) Define any genes you want to mark in volcano plot
      disp_genes <- c()
      
      # Make volcano plots
      plot_volcano(volcano_df, disp_genes, file_suffix)
    }
    
    if (cor_plot==TRUE){
      
      #Determine number of cuts in heatmap
      cuts <- 20
      
      # Use vst as well as rld
      vst <- DESeq2::vst(dds)
      sig_genes <- results %>% 
        dplyr::filter(abs(log2FoldChange) >= 0.58 & padj < 0.05 & padj > 0 & !is.na(padj)) %>%
        dplyr::select(SYMBOL) %>%
        unlist(use.names=FALSE)
      vst <- vst[rownames(vst) %in% sig_genes]
      head(assay(vst))
      sampleDists <- dist(assay(vst))
      sampleDistMatrix <- as.matrix(sampleDists)
      rownames(sampleDistMatrix) <- dplyr::left_join(data.frame("ENSEMBL_ID" = rownames(sampleDistMatrix)), annotations, by=c("ENSEMBL_ID"="ENSEMBL_ID")) %>% 
        dplyr::select(SYMBOL) %>% 
        unlist(use.names = FALSE)
      colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
      
      # cluster and re-order rows
      rowclust <- hclust(dist(sampleDistMatrix))
      
      # cluster and re-order columns
      colclust <- hclust(dist(t(sampleDistMatrix)))
      reordered <- sampleDistMatrix[rowclust$order,colclust$order]
      
      # Save batch corrected normalized counts for entire dataset
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, sheetName = "Correlation")
      openxlsx::writeData(wb, sheet = "Correlation", x = reordered , rowNames = TRUE)
      openxlsx::saveWorkbook(wb,
                             file = paste0(results_path, "Results_", Variable2_value, "_Correlation_Matrix.xlsx"),
                             overwrite = TRUE)
      
      pheatmap::pheatmap(mat = reordered,
                         scale = "none",
                         cellwidth = 3,
                         cellheight = 2,
                         cutree_rows = cuts,
                         cutree_cols = cuts,
                         cluster_rows = TRUE,   #cluster the rows
                         cluster_cols = TRUE,   #cluster the columns
                         fontsize = 8, 
                         fontsize_row = 8, 
                         fontsize_col = 8,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         angle_col = c("270", "0", "45", "90", "315"),
                         fontsize_number = 0.8*fontsize,
                         width = 40,
                         height = 40,
                         filename = paste0(results_path, "Correlation_vst.jpg"))
    }
  }
}

#******************************************************************************#
#                           DECLARE GLOBAL VARIABLES                           #
#******************************************************************************#

# Choose if data is from human or mice. We will adjust gene names accordingly.
species <- "Homo sapiens"
species <- "Mus musculus"

f_suffix <- ""

#****************************DEFINE DESEQ2 VARIABLES***************************#

# NOTE: Make sure there are no white spaces in Target and Reference columns
# in excel file. R will change any white space (" ") to dot ("."). So, replace 
# white space (" ") with underscore ("_") before importing into R.
# NOTE: Metadata MUST have a column named "Batch" before importing into R.
# Example: If samples were isolated on different days & prepared using different
# kits, "Batch" column must have values like "1_Ribo", "1_Poly-A", "2_Ribo", etc

# Define method used for calculating counts
method <- "HTSEQ"
method <- "STAR"

# Define DESeq2 thresholds
padj.cutoff <- 0.1
lfc.cutoff <- 0     # 0.58 for more strict analysis

# Indicate if you want to plot heatmap, volcano plot, PCA plot
heatmap_plot <- FALSE
volcano_plot <- FALSE
cor_plot  <- FALSE

# Folder name corresponding to project containing input files
proj <- "TCGA"
proj <- "IMVigor210"
proj <- "Boopati_MB49Sigma_DDR2KD3"
proj <- "Boopati_MB49Y-_DDR2KO"
proj <- "Mukta_NA13_MB49_CDH12OE"
proj <- "Mukta_GSE95097"
proj <- "Hany_Y"
proj <- "Hany_YKO"
proj <- "Hany_KO"
proj <- "Hany_OE"
proj <- "Hany_OE_Tumor"
proj <- "TRAMP_GSE79756"
proj <- "Hany_Y_Tumor"
proj <- "Prince_DepMap"
proj <- "Jinfen_CRISPR"
proj <- "PRJNA587619"
proj <- "GSE75192"
proj <- "EGAD00001007575"
proj <- "EGAD00001003977"

# parent directory : directory where input files, results, etc are stored
parent_path <- paste0("C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/", proj, "/")
results_path <- paste0("C:/Users/KailasammS/Box/Saravana@cedars/05. Bioinformatics/RNASeq/", proj, "/")

# Define Target and Reference base don project
if (proj == "TCGA"){
  Variable <- "Sex"
  Comparisons <- list(Target = c("Male"),
                      Reference = c("Female"))
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj == "IMVigor210"){
  Variable <- "binaryResponse"
  Comparisons <- list(Target = c("SD/PD"),
                      Reference = c("CR/PR"))
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj =="Mukta_GSE95097"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("8hr", "4hr"),
                      Reference = c("0hr", "0hr"))
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj == "Boopati_MB49Y-_DDR2KO"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Ddr2KO"),
                      Reference = c("Control"))
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj == "Boopati_MB49Sigma_DDR2KD3"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Ddr2KD"),
                      Reference = c("Control"))
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="Mukta_NA13_MB49_CDH12OE"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("OE"),
                      Reference = c("Control"))
  Variable2 <- "Cell_line"
  Variable2_value <- "NA13"  #MB49"
}
if (proj=="Hany_KO"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Kdm5d_KO", "Uty_KO"), 
                      Reference = c("scr_KO", "scr_KO"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="Hany_YKO"){
  Variable <- "Condition"
  Comparisons <- list(Target =    c("YKO_centromere"), 
                      Reference = c("scr_KO_centromere"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="Hany_Y"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Y_neg"),
                      Reference = c("Y_pos"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="Hany_Y_Tumor"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Y_pos_Tumor", "Hany_Tumor"),
                      Reference = c("Y_neg_Tumor", "Sigma_Tumor"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="Hany_OE"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Kdm5d_OE", "Uty_OE"),
                      Reference = c("scr_OE", "scr_OE"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="Hany_OE_Tumor"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Kdm5d_OE_Tumor", "Uty_OE_Tumor"),
                      Reference = c("scr_OE_Tumor", "scr_OE_Tumor"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="TRAMP_GSE79756"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("PAX8-NFE2L2 fusion"),
                      Reference = c("Empty vector"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="TCGA"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Male"),
                      Reference = c("Female"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="Prince_DepMap"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Ypos", "Female"),
                      Reference = c("Yneg", "Yneg"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="scRNASeq_BBN_Nude"){
  Variable <- "Sex"
  Comparisons <- list(Target = c("Male"),
                      Reference = c("Female"))    
  Variable2 <- NULL
  Variable2_value <- NULL
}
if (proj=="scRNASeq"){
  Variable <- "Ystatus"
  Comparisons <- list(Target = c("Yneg"),
                      Reference = c("Ypos"))    
  Variable2 <- "Condition"
  Variable2_value <- "Tumor"
}
if (proj=="Jinfen_CRISPR"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("Day7_2D", "Day14_2D", "Day7_3D", "Day14_3D"),
                      Reference = c("Day0_2D", "Day0_2D", "Day0_3D", "Day0_3D"))    
  Variable2 <- ""
  Variable2_value <- ""
}
if (proj=="PRJNA587619"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("BBN"),
                      Reference = c("Control"))    
  Variable2 <- ""
  Variable2_value <- ""
}
if (proj=="GSE75192"){
  Variable <- "Condition"
  Comparisons <- list(Target = c("09m_Brain", "15m_Brain", "24m_Brain", "30m_Brain", 
                                 "09m_Liver", "15m_Liver", "24m_Liver", "30m_Liver",
                                 "09m_Skin",  "15m_Skin",  "24m_Skin",  "30m_Skin",
                                 "09m_Whole Blood", "15m_Whole Blood", "24m_Whole Blood", "30m_Whole Blood"),
                      Reference = c("02m_Brain", "02m_Brain", "02m_Brain", "02m_Brain",
                                    "02m_Liver", "02m_Liver", "02m_Liver", "02m_Liver",
                                    "02m_Skin",  "02m_Skin",  "02m_Skin", "02m_Skin",
                                    "02m_Whole Blood", "02m_Whole Blood", "02m_Whole Blood","02m_Whole Blood"))    
  Variable2 <- ""
  Variable2_value <- ""
}
if (proj=="EGAD00001007575"){
  Variable <- "Sex"
  Comparisons <- list(Target = c("Male"),
                      Reference = c("Female"))    
  Variable2 <- ""
  Variable2_value <- ""
}
if (proj=="EGAD00001003977"){
  Variable <- "Sex"
  Comparisons <- list(Target = c("Male"),
                      Reference = c("Female"))    
  Variable2 <- ""
  Variable2_value <- ""
}

#***************************DEFINE HEATMAP VARIABLES***************************#

# Define if data is log transformed already
already_log <- FALSE

# Define if data is scaled already
already_scaled <- FALSE

# Define if you want to perform unsupervised row and column clustering
# NOTE: If set to FALSE, samples (columns) and genes(rows) will be arranged 
# in alphabetical order (default) in heatmap. If you want to arranged in 
# specific order, define below.
row_clustering <- TRUE    # Usually TRUE
col_clustering <- TRUE    # Usually FALSE

# Define if you want genes or samples to be arranged in alphabetical order
# NOTE: If set to FALSE, write the plot_genes in order you want in heatmap
# NOTE: If row_clustering==TRUE, then row_clustering_alphabetical is irrelevant
row_clustering_alphabetical <- TRUE
col_clustering_alphabetical <- FALSE

# List annotations you want on heatmap
# NOTE: anno_columns MUST match one of the column names in metadata
anno_columns <- c("Condition")

# Define colors for heatmap
my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(200)
#my_palette <- viridis(200)

# Dimensions of each box in heatmap
cell_width <- NA
cell_height <- NA

#*************************DEFINE VOLCANO PLOT VARIABLES************************#

# Define the target and reference groups
# NOTE: These values MUST be present in "Condition" column of metadata
Target <- Comparisons$Target
Reference <- Comparisons$Reference

# Define if you want to color by padj values or direction (up vs down)
color_by <- "Significance"
color_by <- "Direction"

# Define cutoffs for padj and log2FC
padj_cutoff <- 0.05
log2_cutoff <- 0.58

#******************************************************************************#
#                            STEP 1: GET ANNOTATIONS                           #
#******************************************************************************#

annotations <- get_annotations(species)

#******************************************************************************#
#                STEP 2: COMPILE RAW COUNTS FROM MULTIPLE FILES                #                                           
#******************************************************************************#

# Create a list of all txt files within folder that will be analyzed
count_folder <- paste0(parent_path, "count_results/")
files <- list.files(path=count_folder)

# Create an empty dataframe with 0's
read_data <- data.frame(0)

# Create the reads table 
for (i in 1:length(files)){
  
  # Read the txt file
  temp_file <- read.table(file = paste0(count_folder, files[i]), header = FALSE, sep = "\t")  
  
  # Remove top 4 rows in STAR count output:  
  # N_unmapped, N_multimapping, N_noFeature, N_ambiguous 
  if (method == "STAR"){
    temp_file <- temp_file[5:nrow(temp_file),]
  }
  
  # Remove last 5 rows in HTSeq count output: __no_feature,  __ambiguous, 
  # __too_low_aQual, __not_aligned, __alignment_not_unique
  if (method == "HTSEQ"){
    temp_file <- temp_file[1:(nrow(temp_file)-5),]
  }
  
  # The 1st Column will have Ensembl ids. 
  # For HTSeq, gene counts may be in 2nd or 3rd column. 
  # For STAR, gene counts are in 2nd (unstranded), 3rd (+), 4th (-) column. 
  # Append appropriate column.
  if (method == "HTSEQ" & sum(temp_file[2], na.rm = TRUE) == 0 & sum(temp_file[3], na.rm = TRUE) > 0){
    read_data <- bind_cols(read_data, temp_file[,3])
  } else if (method == "HTSEQ" &  sum(temp_file[2], na.rm = TRUE) > 0 & sum(temp_file[3], na.rm = TRUE) ==0){
    read_data <- bind_cols(read_data, temp_file[,2])                  
  } else if (method == "STAR" & abs((sum(temp_file[2])/sum(temp_file[3])) - (sum(temp_file[2])/sum(temp_file[4]))) < 2){
    print("Unstranded")
    read_data <- bind_cols(read_data, temp_file[,2])                  
  } else if (method == "STAR" & sum(temp_file[3]) > 3*sum(temp_file[4])){
    print("Pos stranded")
    read_data <- bind_cols(read_data, temp_file[,3])                  
  } else if (method == "STAR" & sum(temp_file[4]) > 3*sum(temp_file[3])){
    print("Neg stranded")
    read_data <- bind_cols(read_data, temp_file[,4])                  
  } else{
    print("Error: Gene counts NOT PRESENT in either column 2 or 3 of count file")
  }
  
  # Rename the column names to sample names
  colnames(read_data)[i+1] <- gsub(pattern="\\..*$", replacement="", x=files[i])
}

# Check if all count files have same order of genes in the rows so that files can be merged together
temp_file <- read.table(file = paste0(count_folder, files[1]), header = FALSE, sep = "\t")
gene_list <- temp_file[, 1]
for (i in 1:length(files)){
  
  temp_file <- read.table(file = paste0(count_folder, files[i]), header = FALSE, sep = "\t")
  genes <- temp_file[, 1]
  
  if (!identical(gene_list, genes)){ 
    print("Gene order is different between the count files")
  }
}

if (method == "STAR"){
  gene_list <- gene_list[5:length(gene_list)]
} else if(method == "HTSEQ"){
  gene_list <-gene_list[1:(length(gene_list)-5)]
}

# Add gene names to 1st column
read_data[, 1] <- gene_list
colnames(read_data)[1] <- "SYMBOL"
colnames(read_data) <- gsub(pattern = "ReadsPerGene", replacement = "", x = colnames(read_data))

# Remove genes with 0 counts in all samples
read_data <- read_data[rowSums(read_data[,-1]) != 0,]

# Remove samples with 0 counts in all genes
read_data <- read_data[,colSums(read_data[,-1]) != 0]

# Save the results as xlsx file
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Raw_counts")
openxlsx::writeData(wb, sheet = "Raw_counts", x = read_data)
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Read_data.xlsx"),
                       overwrite = TRUE)

#******************************************************************************#
#                           STEP 3: IMPORT READ DATA                           #                                           
#******************************************************************************#

# read_data MUST be xlsx file
# read_data MUST have "SYMBOL" column with gene names without any duplication
read_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Read_data.xlsx"))

#******************************************************************************#
#                           STEP 4: IMPORT META DATA                           #
#******************************************************************************#

# meta_data MUST be xlsx file
# meta_data MUST have "Sample" column with sample names without any duplication
# meta_data MUST have "Batch" column with batch info
meta_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Meta_data.xlsx"))

#******************************************************************************#
#                              STEP 5: RUN DESEQ2                              #                                           
#******************************************************************************#

analyze_DESeq2(meta_data, read_data, f_suffix)

# #*******************************DIAGNOSTIC TESTS*******************************#
# 
# celltype <- NULL
# 
# # (i) To view counts of specific gene across samples
# plotCounts(dds, gene=which.min(res$padj), intgroup=Variable)           # gene with lowest padj
# plotCounts(dds, gene=which.min(res$log2FoldChange), intgroup=Variable) # gene with lowest log2FC
# plotCounts(dds, gene=which.max(res$log2FoldChange), intgroup=Variable) # gene with highest log2FC
# plotCounts(dds, gene="ENSMUSG00000030598", intgroup=Variable)          # a specific gene
# 
# # (ii) Plot MA
# # The output of plotMA is not a ggplot2 object. So, we cant use ggsave()
# pdf(filename = paste0(diagnostics_path, "Diagnostic_MA_Plot_", celltype, ".pdf"),
#     width = 8.5,
#     height = 11,
#     units = "in",
#     quality = 75,
#     bg = "white",
#     res = 600)
# 
# # Blue points indicate genes with adj p value <0.1
# plotMA(object = res,
#        alpha = 0.1,
#        main = "",
#        xlab = "Mean of Normalized Counts",
#        #ylim = c(-5,5), 
#        MLE = FALSE) 
# 
# dev.off()
# # To identify the genes interactively, run the 2 lines below. 
# # Then click on multiple dots and click Finish. A list of genes will be displayed 
# # idx <- identify(res$baseMean, res$log2FoldChange)
# # rownames(res)[idx]
# 
# # (iii) Plot PCA using rld or vst
# # Variable is defined at beginning of DESeq2. So, we set intgroup=Variable.
# rld <- rlog(object = dds,
#             blind = TRUE,
#             # intercept = ,
#             # betaPriorVar = ,
#             fitType = "parametric")
# 
# vst <- varianceStabilizingTransformation(object = dds,
#                                          blind = TRUE,
#                                          fitType = "parametric")
# 
# for (process in c("rld", "vst")){
#   
#   plotPCA(object = get(process),
#           intgroup = Variable,
#           ntop = 500,
#           returnData = FALSE) +
#     geom_text(aes(label = name), nudge_x = 2, nudge_y = 2)
#   
#   ggplot2::ggsave(filename = paste0("Diagnostic_PCA_plot_using_", process, "_", celltype, ".pdf"),
#                   plot = last_plot(),
#                   device = "pdf",
#                   path = diagnostics_path,
#                   scale = 1,
#                   #width = 8.5,
#                   #height = 11,
#                   units = c("in"),
#                   dpi = 300,
#                   limitsize = TRUE,
#                   bg = "white")
# }
# 
# # (iv) Hierarchical clustering of samples using rld or vst
# rld_mat <- DESeq2::assay(rld)     # extract the matrix from a DESeq2 object
# rld_cor <- cor(x = rld_mat,       # compute pairwise correlation values
#                y = NULL,
#                use = "everything",
#                method = "pearson") 
# 
# # Check the output of cor(), make note of the rownames and colnames
# head(rld_cor)
# 
# vst_mat <- assay(vst)  
# vst_cor <- cor(vst_mat) 
# 
# # Check the output of cor(), make note of the rownames and colnames
# head(vst_cor)    
# 
# for (process in c("rld", "vst")){
#   
#   pheatmap::pheatmap(mat = get(paste0(process, "_cor")),
#                      color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100),
#                      breaks = NA, 
#                      border_color = "white", #"grey60",
#                      cellwidth = NA, 
#                      cellheight = NA, 
#                      scale = "none",   
#                      cluster_rows = TRUE,   #cluster the rows
#                      cluster_cols = TRUE,   #cluster the columns
#                      clustering_distance_rows = "euclidean",
#                      clustering_distance_cols = "euclidean",
#                      clustering_method = "complete",
#                      legend = TRUE, 
#                      legend_breaks = NA,
#                      legend_labels = NA, 
#                      #annotation_row = ,  
#                      #annotation_col = , 
#                      annotation_colors = dplyr::if_else(nrow(col_annotation)+nrow(row_annotation) > 0, ann_colors, NA),
#                      annotation_legend = TRUE,
#                      annotation_names_row = TRUE,
#                      annotation_names_col = TRUE,
#                      show_rownames = dplyr::if_else(nrow(mat)<80, TRUE, FALSE, missing = NULL), 
#                      show_colnames = dplyr::if_else(ncol(mat)<30, TRUE, FALSE, missing = NULL),
#                      fontsize = 8, 
#                      fontsize_row = 8, 
#                      fontsize_col = 8,
#                      angle_col = c("270", "0", "45", "90", "315"),
#                      fontsize_number = 0.8*fontsize, 
#                      #labels_row = display_row,
#                      #labels_col = display_col,
#                      filename = paste0(diagnostics_path, "Diagnostic_Correlation_Heatmap_using_", process, "_", celltype, ".pdf"))
# }
# 
# # (v) Checking if mean < variance (for NB model) or Mean = Variance (for Poisson model). 
# # Each point is a gene denoted by (x,y) where 
# # x = mean_count of gene across all samples & y = variance_count of gene across all samples
# 
# mean_counts <- apply(read_data[,], 1, mean)
# variance_counts <- apply(read_data[,], 1, var)
# df <- data.frame(mean_counts, variance_counts)
# ggplot(df) +
#   geom_point(aes(x=mean_counts, y=variance_counts)) +
#   geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
#   scale_y_log10() +
#   scale_x_log10()
# 
# ggplot2::ggsave(filename = paste0("Diagnostic_Scatter_plot_", celltype, ".pdf"),
#                 plot = last_plot(),
#                 device = "pdf",
#                 path = diagnostics_path,
#                 scale = 1,
#                 #width = 8.5,
#                 #height = 11,
#                 units = c("in"),
#                 dpi = 300,
#                 limitsize = TRUE,
#                 bg = "white")
# 
# # (vi) Plot dispersion estimates. Higher the mean, lower the dispersion
# # The output of plotDispEsts() is not a ggplot2 object. So, we cant use ggsave()
# 
# pdf(filename = paste0(diagnostics_path, "Diagnostic_Dispersion_Estimates_", celltype, ".pdf"),
#     width = 8.5,
#     height = 11,
#     units = "in",
#     quality = 75,
#     bg = "white",
#     res = 600)
# 
# plotDispEsts(object = dds,
#              #ymin,
#              genecol = "black",
#              fitcol = "red",
#              finalcol = "dodgerblue",
#              legend = TRUE,
#              xlab = "Mean of Normalized counts" ,
#              ylab = "Dispersion",
#              log = "xy",
#              cex = 0.45)
# 
# dev.off()
# 
# # # (vii) Plot p-value histogram
# # hist(res$pvalue, col="lightblue", breaks=20)
# # 

# # 
# # # (ix) Plot a histogram for one of the samples to see how the counts are distributed. Adjust xlim if needed
# # ggplot(read_data, aes(x = results.S10_R1_trimmed.fastq.gz.csv)) +
# #   geom_histogram(stat = "bin", bins = 200) +
# #   xlim(-5,1000) +
# #   xlab("Raw expression counts") +
# #   ylab("Number of genes") +
# #   scale_color_brewer(palette="Dark2")
# # 
# # # (x) Extract counts, size factors, etc.. from dds
# # dds <- estimateSizeFactors(dds)         #redundant if you already ran dds <- DESeq(dds)
# # counts <- counts(dds)                   #counts[i,j] = raw_count of gene i in sample j
# # sizefactors <- sizeFactors(dds)         #sizefactors[j] = median (counts[,j]/geometric_mean(counts[i,]))
# # colSums(counts(dds))                    #Total number of raw counts per sample
# # colSums(counts(dds, normalized=TRUE))   #Total number of normalized counts per sample
