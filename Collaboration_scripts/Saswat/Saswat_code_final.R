
# Check "access_TCGA.R" for details on how TCGA data for individual and 
# PanCancer were downloaded and normalized.

#******************************************************************************#
#                 (ICB DATASETS): REFORMAT TCGA CLINICAL DATA                 #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"

gse <- c("GSE78220", "GSE91061", "GSE96619", "GSE115821", "GSE131521")

for (proj in gse){
  
  meta_data <- read.xlsx(paste0(data_path, proj, ".Metadata.xlsx"))
  clinical_data <- read.xlsx(paste0(data_path, "ICB_Clinical_Data.xlsx"))
  
  # Combine relevant info from both dataframe
  clinical_data <- clinical_data %>% 
    dplyr::filter(series_id %in% proj) %>%
    dplyr::mutate(Time = as.numeric(OS)/30,
                  Status = OS_FLAG,
                  Sample = toupper(Sample))
  
  meta_data <- meta_data %>%
    dplyr::rename(Sample_ID = GEO_Accession__exp_) %>%
    dplyr::mutate(Description = toupper(Description),
                  Run = toupper(Run))
  
  if(proj == "GSE78220"){
    metadata <- meta_data %>% dplyr::inner_join(clinical_data,by=c("Description"="Sample"))
  } 
  else {
    metadata <- meta_data %>% dplyr::inner_join(clinical_data,by=c("Run"="Sample"))
  }
  
  # Print number of missing samples
  print(nrow(metadata)-nrow(meta_data))
  print(nrow(metadata)-nrow(clinical_data))
  
  # Save metadata
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Metadata")
  openxlsx::writeData(wb, sheet="Metadata", x=metadata)
  openxlsx::saveWorkbook(wb, file=paste0(output_path, make.names(proj), ".Metadata.xlsx"), 
                         overwrite=TRUE)
}

#******************************************************************************#
#             (ICB DATASETS): CALCULATE DESEQ2 NORMALIZED COUNTS               #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"

gse <- c("GSE78220", "GSE91061", "GSE96619", "GSE115821", "GSE131521", "IMvigor010", "IMvigor210")
species <- "Homo sapiens"
annotations <- get_annotations(species)

for (proj in gse){
  
  # Import metadata and read data
  meta_data <- read.xlsx(paste0(data_path, proj, ".Metadata.xlsx"))
  read_data <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, proj,".raw.counts.xlsx"))
  
  # Rename appropiate column in meta_data to "Sample_ID"
  # meta_data <- meta_data %>% dplyr::rename(Sample_ID = GEO_Accession__exp_)
  
  # Rename appropiate column in read_data to "SYMBOL"
  if (proj %in% c("GSE78220", "GSE91061", "GSE96619", "GSE115821", "GSE131521")){
    read_data <- read_data %>% dplyr::rename(SYMBOL = GeneID)
  }
  
  # Define comparisons
  Comparisons <- list(Variable =c(NA),
                      Target   =c(NA),
                      Reference=c(NA))
  
  # Perform QC
  meta_data <- prep_metadata(meta_data, read_data)
  read_data <- prep_readdata(read_data, meta_data)
  l <- check_data(read_data, meta_data)
  meta_data <- l[[2]]
  read_data <- l[[1]]
  
  # Prepare DESeq2 object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~ 1)
  
  # Calculate normalized counts using DESeq2()
  norm_counts <- deseq2_norm_counts(dds, meta_data, annotations)
  
  # Save batch corrected normalized counts
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Normalized counts")
  openxlsx::writeData(wb, sheet="Normalized counts", x=norm_counts)
  openxlsx::saveWorkbook(wb, file=paste0(output_path, make.names(proj), ".Normalized.counts.xlsx"), 
                         overwrite=TRUE)
}

#******************************************************************************#
#             CLASSIFY INTO PATIENTS INTO IMMUNE HOT vs IMMUNE COLD            #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/results/"

# List the 4 groups of TME from https://doi.org/10.1016/j.ccell.2021.04.014
angio <- c("Angiogenesis","Endothelium","Cancer.associated.fibroblasts",
           "Matrix", "Matrix.remodeling")
pro_tumor_immune_infiltrate <- c("Tumor.associated.Macrophages", "Macrophage.and.DC.traffic",
                                 "Myeloid.cells.traffic","Immune.Suppression.by.Myeloid.Cells",
                                 "Th2.signature","Treg.and.Th2.traffic", "Treg", "M1.signature",
                                 "Protumor.cytokines", "Neutrophil.signature", "Granulocyte.traffic")
anti_tumor_immune_infiltrate <- c("MHCII", "Antitumor.cytokines", "Co.activation.molecules",
                                  "B.cells", "NK.cells", "Checkpoint.molecules", "Effector.cells",
                                  "T.cells", "Th1.signature", "Effector.cell.traffic", "MHCI")
malignant_cell <- c("EMT.signature", "Tumor.proliferation.rate")

# Read signature gene sets for 29 pathways but keep only immune infiltrate 
# pathways as Immune Enriched and Immune Depleted are defined based on them.
# NOTE: HLA-A, HLA-B etc are present in sig dataframe but appear as HLA.A, HLA.B
# in normalized counts dataframe. So, fix this.
sig <- read.xlsx(paste0(file_path, "Saswat_Sig_list.xlsx"))
sig <- sig[-1, -1] %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames(colnames(.)[1]) %>% 
  t() %>%
  data.frame() %>% 
  dplyr::select(all_of(pro_tumor_immune_infiltrate), all_of(anti_tumor_immune_infiltrate)) %>%
  dplyr::mutate(across(.cols=everything(), 
                       ~ gsub(pattern = "-", replacement = ".", x=.)))

# Read pan cancer metadata
meta_data <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"))
proj <- make.names(unique(meta_data$Project_ID))
proj <- c(proj, "IMvigor010", "IMvigor210", "GSE78220", "GSE91061", "GSE96619", "GSE115821", "GSE131521", "TCGA.PanAtlas")

# Import expression data that has already been normalized across samples
for (f in proj){
  
  meta_data <- read.xlsx(paste0(data_path, f, ".Metadata.xlsx"))
  
  if (f == "TCGA.PanAtlas"){
    read_data <- read.table(paste0(data_path, f, ".Normalized.counts.tsv"), header=TRUE)
    read_data <- read_data[,-c(2:3)]
  } else if(f %in% c("IMvigor010", "IMvigor210")){
    read_data <- read.xlsx(paste0(data_path, f, ".Normalized.counts.xlsx"))
    read_data <- read_data[,-c(2:4)]
  } else{
    read_data <- read.xlsx(paste0(data_path, f, ".Normalized.counts.xlsx"))
    read_data <- read_data[,-c(2:3)]
  }
  colnames(read_data) <- make.names(colnames(read_data))
  
  # Reformat readdata
  norm_counts <- read_data %>%
    # Replace NA with 0
    base::replace(is.na(.), 0) %>%                               
    # If there are duplicated genes, keep only row for highest expressing copy
    dplyr::mutate(n = rowSums(.[,-1])) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(n) %>%
    dplyr::ungroup() %>%
    # Remove rows with 0 expression across all columns
    dplyr::filter(n != 0) %>%
    dplyr::select(everything(), -n) %>%
    # Move gene names to rownames
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("SYMBOL") %>%
    # Make sure all columns are numeric
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric)) %>%
    # Keep only samples in metadata
    dplyr::select(all_of(intersect(make.names(meta_data$Sample_ID), make.names(colnames(read_data)))))
  
  colnames(norm_counts) <- base::make.names(names = colnames(norm_counts))
  
  log_norm_counts <- log(1+norm_counts, base=2)
  t <- base::apply(X=log_norm_counts, MARGIN=1, FUN=median, na.rm=TRUE)
  log_norm_counts <- base::sweep(x=log_norm_counts, MARGIN=1, FUN="-", STATS=t)
  
  # For Pancancer, save median centered log norm counts as this step is time consuming
  if (f == "TCGA.PanAtlas"){
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "median centered log normalized")
    openxlsx::writeData(wb, sheet = "median centered log normalized", x = log_norm_counts, rowNames=TRUE)
    openxlsx::saveWorkbook(wb, file = paste0(output_path, "TCGA.PanAtlas.median_centered_log_norm_counts.xlsx"), overwrite = TRUE)
  }
  
  #***** Calculate scores for each patient for each of the 22 signatures *****#
  
  # https://support.bioconductor.org/p/9141192/
  # https://support.bioconductor.org/p/107466/
  # NOTE: DESeq2 normalized counts are NOT normalized to gene length. So, smaller
  # genes have fewer reads and hence smaller normalized counts. This is not an
  # issue if we are comparing same gene across different samples, however, in 
  # ssgsea and gsva, we rank genes within a sample based on expression and 
  # calculate enrichment scores for each sample. So, it COULD BE AN ISSUE.
  
  # Generate a list of gene sets
  gs <- list()
  for (i in 1:ncol(sig)){
    plot_genes <- list(sig[,i][!is.na(sig[,i])])
    names(plot_genes) <- colnames(sig)[i]
    gs <- c(gs, plot_genes)
  }
  
  # Calculate GSVA scores
  # NOTE: gsva works better than ssgsea for normalized counts
  # NOTE: ssgsea works better than gsva for TPM counts
  # This observation is based on Y score of TCGA BLCA (normalized counts) and 
  # DepMap cell lines (TPM counts). gsva classifies all Female as Ylow for TCGA
  # BLCA but not for DepMap. Similarly, ssgsea classifies all Female as Ylow for
  # DepMap but not for TCGA BLCA.
  gsvaPar <- GSVA::gsvaParam(exprData = as.matrix(log_norm_counts), 
                             geneSets = as.list(gs))
  gsva.scores <- gsva(gsvaPar, 
                      verbose=TRUE)
  
  # Calculate SSGSEA scores
  ssgseaPar <- GSVA::ssgseaParam(exprData = as.matrix(log_norm_counts), 
                                 geneSets = as.list(gs))
  ssgsea.scores <- GSVA::gsva(ssgseaPar, 
                              verbose=TRUE)
  
  # Calculate z-scores (based on Levine et al https://doi.org/10.1186/gb-2006-7-10-r93)
  # Create empty dataframe of samples as rownames and pathway as column names
  z.scores <- log_norm_counts %>% 
    t() %>%
    data.frame() %>%
    dplyr::select(colnames(.)[1]) %>%  # Keep a dummy column to keep dataframe structure intact
    tibble::rownames_to_column("Sample")
  
  for (i in 1:ncol(sig)){
    
    plot_genes <- sig[,i][!is.na(sig[,i])] 
    expr_df <- as.data.frame(advanced_Z(plot_genes, log_norm_counts))
    colnames(expr_df) <- colnames(sig)[i]
    expr_df <- expr_df %>%
      tibble::rownames_to_column("Sample")
    
    z.scores <- z.scores %>% 
      dplyr::left_join(expr_df, by=c("Sample"="Sample"))
  }
  
  z.scores <- z.scores[,-2]  # remove dummy column
  z.scores <- z.scores %>% 
    tibble::column_to_rownames("Sample") %>%
    t()
  
  rownames(z.scores) <- make.names(rownames(z.scores))
  colnames(z.scores) <- make.names(colnames(z.scores))
  
  #***** Classify Patients into Immune Enriched vs Immune Depleted groups *****#
  
  for (s in c("gsva.scores", "ssgsea.scores", "z.scores")){
    
    scores <- get(s)
    mat <- scores
    
    # Col annotation for heatmap
    col_ann <- scores %>%
      t() %>%
      data.frame() %>%
      dplyr::mutate(Mean = rowMeans(as.matrix(.)),
                    Column_Groups = dplyr::case_when(Mean > 0 ~ "Immune Enriched",
                                                     TRUE ~ "Immune Depleted")) %>%
      dplyr::select(Column_Groups, Mean)
    
    col_elements <- col_ann %>%
      dplyr::count(Column_Groups) %>%
      dplyr::filter(n>1) %>%
      dplyr::select(Column_Groups) %>%
      unlist(use.names=FALSE)
    
    # Row annotation for heatmap
    row_ann <- data.frame(Signature=c(pro_tumor_immune_infiltrate, anti_tumor_immune_infiltrate),
                          Row_Groups=c(rep(x="Pro-tumor Immune infiltrate",times=11),
                                       rep(x="Anti-tumor Immune infiltrate",times=11)),
                          Dummy_col = NA) %>%
      tibble::column_to_rownames("Signature")
    
    row_elements <- row_ann %>% 
      dplyr::count(Row_Groups) %>% 
      dplyr::filter(n>1) %>%
      dplyr::select(Row_Groups) %>%
      unlist(use.names=FALSE)
    
    # Arrange columns in increasing order of mean
    col_order <- rownames(col_ann %>% dplyr::arrange(Mean))
    
    # Within group row clustering
    row_order <- c()
    for (g in row_elements){
      temp_mat <- mat[rownames(row_ann)[which(row_ann$Row_Groups == g)],]
      rowclust <- hclust(dist(temp_mat))
      row_order <- c(row_order,rownames(temp_mat[rowclust$order,]))
    }
    
    # We cluster all rows
    row_order <- c()
    temp_mat <- mat
    rowclust <- hclust(dist(temp_mat))
    row_order <- c(row_order,rownames(temp_mat[rowclust$order,]))
    
    # Rearrange rows and columns in mat
    mat <- mat[row_order, col_order]
    
    # Define breaks
    if(max(mat) == 0){
      breaks <- c(seq(from = min(mat), to = 0, length.out = 100))
      my_palette <- my_palette[1:50]
    } else if (min(mat) == 0){
      breaks <- c(seq(from = 0, to = max(mat), length.out = 100))
      my_palette <- my_palette[50:100]
    } else if(min(mat) < -3 | max(mat) > 3){
      breaks <- c(seq(-1.5, 0, length.out = 50), seq(1.5/100, 1.5, length.out = 50))
    } else{
      breaks <- c(seq(from = min(mat), to = 0, length.out = 50), seq(from = max(mat)/100, to = max(mat), length.out = 50))
    }
    
    # Define colors for heatmap
    vrds <- viridis_pal()(100)
    rdbu <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
    
    # Define colors for row and col legend
    ann_colors = list(Column_Groups = c(`Immune Depleted` = "#CB181D", `Immune Enriched` = "#A6D854"),
                      Row_Groups = c(`Pro-tumor Immune infiltrate` = "white", `Anti-tumor Immune infiltrate` = "white"))
    
    # Define gaps
    gaps_row <- (row_ann %>% dplyr::count(Row_Groups) %>% dplyr::mutate(n = cumsum(n)) %>% dplyr::select(n) %>% unlist(use.names = FALSE))
    gaps_row <- gaps_row[gaps_row < nrow(mat)]
    gaps_col <- (col_ann %>% dplyr::count(Column_Groups) %>% dplyr::mutate(n = cumsum(n)) %>% dplyr::select(n) %>% unlist(use.names = FALSE))
    gaps_col <- gaps_col[gaps_col < ncol(mat)]
    
    # Perform row scaling if you want to compare across columns
    row_mat <- t(scale(t(mat)))
    
    # Perform col scaling
    # col_mat <- scale(mat)
    
    # Plot heatmap
    pheatmap(mat = row_mat, 
             scale = "none",
             # scale = "row",    # equivalent to mat = t(scale(t(mat))) & scale = "none"
             # scale = "column", # equivalent to mat = scale(mat)       & scale = "none"
             breaks = breaks,
             color = rdbu,
             border_color = NA, #"white",
             #cellwidth = 1,  # defined in points. 1 inch = 72 points
             #cellheight = 1, # defined in points. 1 inch = 72 points
             fontsize = 10,   # defined in points. 1 inch = 72 points
             annotation_colors = ann_colors,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             #gaps_row = gaps_row, 
             gaps_col = gaps_col,
             #annotation_row = row_ann %>% dplyr::select(Row_Groups), 
             annotation_col = col_ann %>% dplyr::select(Column_Groups),
             annotation_legend = TRUE,
             show_rownames = FALSE, 
             show_colnames = FALSE, 
             annotation_names_row = FALSE, 
             annotation_names_col = FALSE,
             #treeheight_row = 0,
             #treeheight_col = 0,
             # IMPORTANT: Give atleast 10 inches extra in width, so that there is sufficient margin space
             width = 20,
             #width = ceiling(ncol(mat)*5/72) + 10,  # defined in inches. 1 inch = 72 points
             #height = ceiling(nrow(mat)*5/72) + 2,
             filename=paste0(output_path, f, "_Heatmap_", s, ".tiff"))
  }
  #************ Merge patient classifications based on 3 approaches ***********#
  #*** Append NPEPPS expression as well as 22 scores based on 3 approaches ***#
  #******* Identify top 33 and bottom 33 enriched and depleted patients *******#
  
  # Get NPEPPS median centered log norm counts
  npepps_zscore <- log_norm_counts %>% 
    tibble::rownames_to_column("SYMBOL") %>%
    dplyr::filter(SYMBOL == "NPEPPS") %>% 
    tibble::column_to_rownames("SYMBOL") %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::mutate(Sample = make.names(Sample))
  
  # Get NPEPPS normalized counts
  npepps <- norm_counts %>% 
    tibble::rownames_to_column("SYMBOL") %>% 
    dplyr::filter(SYMBOL == "NPEPPS") %>% 
    tibble::column_to_rownames("SYMBOL") %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::mutate(Sample = make.names(Sample))
  
  # Merge classifications from all 3 approaches
  summary_df <- data.frame(Sample = rownames(t(z.scores)))
  for (s in c("gsva.scores", "ssgsea.scores", "z.scores")){
    
    scores <- get(s)
    col_ann <- scores %>%
      t() %>%
      data.frame() %>%
      dplyr::mutate(Mean = rowMeans(as.matrix(.)),
                    Class = dplyr::case_when(Mean > 0 ~ "Immune Enriched",
                                             TRUE ~ "Immune Depleted")) %>%
      tibble::rownames_to_column("Sample") %>%
      dplyr::select(Sample, Class, Mean, Effector.cells) %>%
      dplyr::mutate(Sample = make.names(Sample)) %>%
      dplyr::rename(!!rlang::sym(gsub(pattern=".scores",replacement=".mean",          x=rlang::sym(s))) := Mean,
                    !!rlang::sym(gsub(pattern=".scores",replacement=".classification",x=rlang::sym(s))) := Class,
                    !!rlang::sym(gsub(pattern=".scores",replacement=".CTL",           x=rlang::sym(s))) := Effector.cells)
    
    summary_df <- summary_df %>% 
      dplyr::left_join(col_ann, by=c("Sample"="Sample"))
  }
  
  # Merge NPEPPS info with classification info
  summary_df <- summary_df %>% 
    dplyr::mutate(Comments = dplyr::case_when(gsva.classification == ssgsea.classification & 
                                                gsva.classification == z.classification ~ "OK",
                                              TRUE ~ "AMBIGUOS CLASSIFICATION")) %>%
    dplyr::mutate(gsva.33 = dplyr::case_when(gsva.mean <= quantile(summary_df$gsva.mean, c(0.33))[[1]] |
                                               gsva.mean >= quantile(summary_df$gsva.mean, c(0.66))[[1]] ~ gsva.classification,
                                             TRUE ~ NA),
                  ssgsea.33 = dplyr::case_when(ssgsea.mean <= quantile(summary_df$ssgsea.mean, c(0.33))[[1]] |
                                                 ssgsea.mean >= quantile(summary_df$ssgsea.mean, c(0.66))[[1]] ~ ssgsea.classification,
                                               TRUE ~ NA),
                  z.33 = dplyr::case_when(z.mean <= quantile(summary_df$z.mean, c(0.33))[[1]] | 
                                            z.mean >= quantile(summary_df$z.scores, c(0.66))[[1]] ~ z.classification,
                                          TRUE ~ NA)) %>%
    dplyr::left_join(npepps_zscore, by=c("Sample"="Sample")) %>%
    dplyr::left_join(npepps, by=c("Sample"="Sample"))
  
  cat(nrow(summary_df %>% dplyr::filter(Comments != "OK")), "out of ", nrow(summary_df), "samples ambiguous in ", f, "\n")
  
  # Save heatmap details and NPEPPS expression
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "NPEPPS")
  openxlsx::writeData(wb, sheet = "NPEPPS", x = summary_df, rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName = "gsva.scores")
  openxlsx::writeData(wb, sheet = "gsva.scores", x = t(gsva.scores), rowNames = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "ssgsea.scores")
  openxlsx::writeData(wb, sheet = "ssgsea.scores", x = t(ssgsea.scores), rowNames = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "z.scores")
  openxlsx::writeData(wb, sheet = "z.scores", x = t(z.scores), rowNames = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "Heatmap_matrix")
  openxlsx::writeData(wb, sheet = "Heatmap_matrix", x = row_mat, rowNames = TRUE)
  openxlsx::saveWorkbook(wb, file = paste0(output_path, f, ".Scores.&.Heatmap.details.xlsx"), overwrite = TRUE)
  
}

#******************************************************************************#
#      PERFORM DIFFERENTIAL EXPRESSION IMMUNE COLD vs IMMUNE HOT PATIENTS      #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/results/"

species <- "Homo sapiens"
annotations <- get_annotations(species)

# Read pan cancer metadata
meta_data <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"))
proj <- make.names(unique(meta_data$Project_ID))
proj <- c(proj, "IMvigor010", "IMvigor210")

# Import expression data that has already been normalized across samples
for (f in proj){
  
  # Merge immune enriched-depleted classification with metadata
  meta_data <- read.xlsx(paste0(data_path, f, ".Metadata.xlsx"))
  classification_data <- read.xlsx(paste0(output_path, f, ".Scores.&.Heatmap.details.xlsx"), sheet = "NPEPPS")
  meta_data <- meta_data %>% 
    dplyr::left_join(classification_data, by=c("Sample_ID"="Sample"))
  
  # Import read_data
  if (f == "TCGA.PanAtlas"){
    read_data <- read.table(paste0(data_path, f, ".raw.counts.tsv"), header=TRUE)
    read_data <- read_data[,-c(2:3)]
  } else if(f %in% c("IMvigor010", "IMvigor210")){
    read_data <- read.xlsx(paste0(data_path, f, ".raw.counts.xlsx"))
    read_data <- read_data[,-c(2:4)]
  } else{
    read_data <- read.xlsx(paste0(data_path, f, ".raw.counts.xlsx"))
    read_data <- read_data[,-c(2)]
  }
  
  # Rename appropiate column in meta_data to "Sample_ID"
  # meta_data <- meta_data %>% dplyr::rename(Sample_ID = GEO_Accession__exp_)
  
  # Rename appropiate column in read_data to "SYMBOL"
  #read_data <- read_data %>% dplyr::rename(SYMBOL = GeneID)
  colnames(read_data)[1] <- "SYMBOL"
  
  # Define comparisons
  Comparisons <- list(Variable    = c("gsva.classification", "ssgsea.classification", "z.classification"),
                      Target      = c("Immune Depleted", "Immune Depleted", "Immune Depleted"),
                      Reference   = c("Immune Enriched", "Immune Enriched", "Immune Enriched"),
                      lfc.cutoff  = 0,
                      padj.cutoff = 0.1)
  
  # Perform QC
  meta_data <- prep_metadata(meta_data, read_data)
  read_data <- prep_readdata(read_data, meta_data)
  l <- check_data(read_data, meta_data)
  meta_data <- l[[2]]
  read_data <- l[[1]]
  
  # Perform analysis for every comparison
  for (n in 1:length(Comparisons$Variable)){
    
    # This generates a new column "id" that has info on samples being comparared
    meta_data_comp <- meta_data %>%
      dplyr::mutate(id=get(Comparisons$Variable[n]))
    
    # Perform DESeq2() using in-built batch modelling
    approach <- "DESeq2"
    if (length(unique(meta_data_comp$Batch)) > 1){
      dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                            colData=meta_data_comp, 
                                            design=~ Batch+id)
    } else {
      dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                            colData=meta_data_comp, 
                                            design=~ id)
    }
    
    prefix <- make.names(f)
    dds <- run_deseq2(dds, meta_data_comp, annotations, Comparisons, n, approach, prefix, output_path)
  }
}

#******************************************************************************#
#                  IDENTIFY FREQUENCY OF DEGS ACROSS 33 CANCERS                #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/results/"

# Read pan cancer metadata
meta_data <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"))
proj <- make.names(unique(meta_data$Project_ID))
proj <- c(proj)

# Create workbook to store results
wb <- openxlsx::createWorkbook()

for (s in c("gsva", "ssgsea", "z")){
  
  # Create a list of all DEGs from all excel files
  deg_genes <- c()
  for (f in proj){
    degs <- read.xlsx(paste0(output_path, f, ".DEGs.", s, ".classification_Immune Depleted_vs_Immune Enriched_DESeq2_.xlsx")) %>%
      dplyr::filter(padj <=0.05)
    deg_genes <- c(deg_genes, degs$SYMBOL)
  }
  
  # Create an empty dataframe to store LFC of each DEG for each cancer
  all_degs <- unique(deg_genes)
  df <- data.frame(SYMBOL = all_degs, 
                   DUMMY = 0)
  
  for (f in proj){
    degs <- read.xlsx(paste0(output_path, f, ".DEGs.", s, ".classification_Immune Depleted_vs_Immune Enriched_DESeq2_.xlsx")) %>%
      dplyr::filter(padj <=0.05) %>%
      dplyr::select(SYMBOL, log2FoldChange) %>%
      dplyr::rename(!!rlang::sym(f) := "log2FoldChange") %>%
      dplyr::distinct_at("SYMBOL", .keep_all = TRUE)
    
    df <- df %>% 
      dplyr::left_join(degs, by=c("SYMBOL"="SYMBOL"))
  }  
  
  # Dataframe has genes as rownames and expression for each cancer on columns
  df <- df %>%
    dplyr::select(c(everything(), -c("DUMMY"))) %>%
    tibble::column_to_rownames("SYMBOL")
  
  # Create 3 dataframes
  # total: any expression is marked as 1, 
  # up   : upregulated genes are marked as 1
  # down : downregulated genes are marked as 1
  df_total <- df
  df_up <- df
  df_down <- df
  
  df_total[!is.na(df)] <- 1
  
  df_up[df > 0] <- 1
  df_up[df <= 0] <- 0
  
  df_down[df >= 0] <- 0
  df_down[df < 0] <- 1
  
  df <- df %>% 
    dplyr::mutate(n_TOTAL = rowSums(df_total, na.rm=TRUE), 
                  n_UP = rowSums(df_up, na.rm=TRUE), 
                  n_DOWN = rowSums(df_down, na.rm=TRUE)) %>%
    dplyr::select(n_TOTAL, n_UP, n_DOWN, everything()) %>%
    dplyr::arrange(desc(n_UP))
  
  # Save
  openxlsx::addWorksheet(wb, sheetName = paste0("Summary_", s))
  openxlsx::writeData(wb, sheet = paste0("Summary_", s), x = df, rowNames = TRUE)
}
openxlsx::saveWorkbook(wb, file = paste0(output_path, "zSummary_of_DEGs.xlsx"), overwrite = TRUE)

# gsva classification "NPEPPS" is UP in 13 Cancer and DOWN in 1 cancer. 
# So, choose gsva for further analysis.

#******************************************************************************#
#                  2D UMAP OF PATHWAY SCORES ACROSS 33 CANCERS                 #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/results/"

# Read pan cancer metadata
meta_data <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"))
proj <- make.names(unique(meta_data$Project_ID))
proj <- c(proj)

# Create workbook to store UMAP coordinates
wb <- openxlsx::createWorkbook()

# Choose classification method
s <- "gsva"

# Create an empty dataframe with Pathway names as columns
mat <- read.xlsx(paste0(output_path, proj[1], ".Scores.&.Heatmap.details.xlsx"),
                 sheet = paste0(s, ".scores"))
colnames(mat)[1] <- "Sample_ID"
mat <- mat[1,]

# Merge scores from all 33 cancers
for (f in proj){
  data <- read.xlsx(paste0(output_path, f, ".Scores.&.Heatmap.details.xlsx"),
                    sheet = paste0(s, ".scores"))
  colnames(data)[1] <- "Sample_ID"
  mat <- dplyr::bind_rows(mat, data)
}

mat <- mat %>% 
  dplyr::distinct_at("Sample_ID", .keep_all = TRUE) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("Sample_ID")

# PCA is performed on columns. So, have the variables i.e. genes on columns.
# "Error cannot rescale a constant/zero column to unit variance" is due to 
# genes that have 0 expression across all samples
# Also, NA must be replaced with 0
# prcomp() is better than princomp()

mat[is.na(mat)] <- 0
cat(s, ":", sum(rowSums(mat) == 0), "\n")
mat <- mat[rowSums(mat) != 0,]

# pr <- prcomp(x = mat, center = TRUE, scale = TRUE)
# umap <- umap(d = as.matrix(pr$x),
#              n_components = 3,
#              config = umap.defaults,
#              method = c("naive"),
#              preserve.seed = TRUE)

umap <- umap::umap(d = mat,
                   n_components = 3, # umap cordinates differ based on n_components
                   config = umap.defaults,
                   method = c("naive"),
                   preserve.seed = TRUE)

umap_val <- umap$layout %>%
  data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::left_join(meta_data %>% 
                     dplyr::select(Sample_ID, Project_ID) %>% 
                     dplyr::mutate(Sample_ID = make.names(Sample_ID)), 
                   by=c("Sample"="Sample_ID"))
colnames(umap_val) <- c("Sample", "UMAP1", "UMAP2", "UMAP3", "Project_ID")

ggplot(data = umap_val, aes(x=UMAP1, y=UMAP2, color=Project_ID, fill=Project_ID)) +
  geom_point(size=1) + 
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme_bw() +
  guides(color = guide_legend(override.aes=list(shape = 22, size = 5, color="black"))) +
  scale_color_manual(values = c("#000000","#D9D9D9","#003C30","#beb9db","#1D91C0","#A6CEE3",
                                "#50E991","#A6D854","#74C476","#C7E9B4","#00bfa0","#E5F5F9",
                                "#EDBF33","#E6D800","#FFF7BC","#ffee65","#C7EAE5",
                                "#67001F","#CB181D","#FD8D3C","#FC9272","#EF3B2C","#F16913",
                                "#9b19f5","#6A51A3","#762A83","#D4B9DA","#0bb4ff",
                                "#E60049","#DC0AB4","#AE017E","#DF65B0","#FDCCE5")) +
  scale_fill_manual(values = c("#000000","#D9D9D9","#003C30","#beb9db","#1D91C0","#A6CEE3",
                               "#50E991","#A6D854","#74C476","#C7E9B4","#00bfa0","#E5F5F9",
                               "#EDBF33","#E6D800","#FFF7BC","#ffee65","#C7EAE5",
                               "#67001F","#CB181D","#FD8D3C","#FC9272","#EF3B2C","#F16913",
                               "#9b19f5","#6A51A3","#762A83","#D4B9DA","#0bb4ff",
                               "#E60049","#DC0AB4","#AE017E","#DF65B0","#FDCCE5"))

ggsave(paste0(output_path, "UMAP_", s, ".tiff"))

# Save
openxlsx::addWorksheet(wb, sheetName = paste0("UMAP_", s))
openxlsx::writeData(wb, sheet = paste0("UMAP_", s), x = umap_val, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(output_path, "UMAP_Coords.xlsx"), overwrite = TRUE)

#******************************************************************************#
#       3D UMAP OF IMMUNE ENRICHED & DEPLETED SAMPLES ACROSS 33 CANCERS        #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/results/"

# Read pan cancer metadata
meta_data <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"))
proj <- make.names(unique(meta_data$Project_ID))
proj <- c(proj)

# Create workbook to save sample classification into Immune Enriched vs Depleted
wb <- openxlsx::createWorkbook()

for (s in c("gsva", "ssgsea", "z")){
  
  # Create an empty dataframe to store sample classifications
  group_data <- data.frame(Sample=c(""))
  
  for (f in proj){
    df <- read.xlsx(paste0(output_path, f, ".Scores.&.Heatmap.details.xlsx"),
                    sheet="NPEPPS") %>%
      dplyr::select(Sample, paste0(s, ".classification"), paste0(s, ".mean"), NPEPPS.x, NPEPPS.y)
    colnames(df) <- c("Sample", "Class", "Mean", "NPEPPS.x", "NPEPPS.y")
    
    group_data <- dplyr::bind_rows(group_data, df)
  }
  group_data <- group_data[-1,]
  
  # Save worksheet
  openxlsx::addWorksheet(wb, sheetName = paste0(s, ".classification"))
  openxlsx::writeData(wb, sheet = paste0(s, ".classification"), x = group_data, rowNames = FALSE)
}

openxlsx::saveWorkbook(wb, file = paste0(output_path, "Two_Group_Classification_from_Individual_cancers.xlsx"), overwrite = TRUE)

s <- "gsva"
mat <- read.xlsx(paste0(output_path, "UMAP_Coords.xlsx"), 
                 sheet = paste0("UMAP_", s))

group <- read.xlsx(paste0(output_path, "Two_Group_Classification_from_Individual_cancers.xlsx"), 
                   sheet = paste0(s, ".classification"))

mat <- mat %>% 
  dplyr::left_join(group, by=c("Sample"="Sample")) %>%
  dplyr::mutate(color = dplyr::case_when(Class == "Immune Depleted" ~ "#CB181D",
                                         TRUE ~ "#A6D854"))

# Save the plot
jpeg(file=paste0(output_path, "UMAP_", s, "_3D.tiff"))

# If you get "Error in check.plt(parplt) : figure margins too large,type 
# par("mar"). You will see values as [1] 5.1 4.1 4.1 2.1. Change like below
par(mar=c(1,1,1,1))
scatter3D(x=mat$UMAP1, y=mat$UMAP2, z=mat$UMAP3, 
          theta = 40, phi = 15,  # viewing angles
          bty = "b2",      # back panels and grid lines are visible
          colvar = NULL,   # avoid coloring by z variable
          col = mat$color, # variable used to color points
          pch = 19,        # shape of points
          cex = 0.25,       # size of points
          colkey = list(mat$Class),
          clab=c("Class"))

dev.off()

#******************************************************************************#
#HORIZONTAL STACKED BAR PLOT OF IMMUNE ENRICHED & DEPLETED SAMPLES ACROSS 33 CANCERS#
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/results/"

# Read pan cancer metadata
meta_data <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx")) %>% 
  dplyr::select(Sample_ID, Project_ID) %>%
  dplyr::mutate(Sample_ID = make.names(Sample_ID))

s <- "gsva"

# Read 2 group classification info
group <- read.xlsx(paste0(output_path, "Two_Group_Classification_from_Individual_cancers.xlsx"), 
                   sheet = paste0(s, ".classification")) %>%
  dplyr::left_join(meta_data, by=c("Sample"="Sample_ID")) %>%
  dplyr::group_by(Project_ID) %>%
  dplyr::count(Class) %>%
  dplyr::mutate(percent = round(100*n/sum(n, na.rm=TRUE), digits = 2), label_percent = paste0(percent,"%"))

# Plot stacked bar chart
ggplot2::ggplot(data = group, aes(x = percent, y = Project_ID, fill = Class)) +
  geom_bar(stat = "identity", position = "stack", width = 0.95) +
  theme_classic() + 
  scale_fill_manual(values = c("#CB181D","#A6D854")) +
  # scale_fill_brewer(palette = "Set1",
  #                   aesthetics = "fill") +
  geom_text(aes(x=percent, label=label_percent), 
            position=position_stack(vjust=0.5), 
            fontface="bold", colour="white", size=4, check_overlap=TRUE) +
  ggplot2::labs(title = "",
                fill = "Immune Class",
                x = "Percent composition",
                y = "")

ggsave(paste0(output_path, "zPercentages_", s, ".tiff"))

# Create workbook to save percentages of Immune Enriched vs Depleted
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = paste0(s, ".classification"))
openxlsx::writeData(wb, sheet = paste0(s, ".classification"), x = group_data, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(output_path, "zPercentages.xlsx"), overwrite = TRUE)


#******************************************************************************#
#     NPEPPS expression high in KRT13 population
#************************************************************************

proj <- "scRNASeq_Simon"
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")
source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/results/"

integrated_seurat <- readRDS(paste0(seurat_results, "Simon_integrated_seurat_snn.rds"))

Idents(integrated_seurat) <- "subtype"
Seurat::FeaturePlot(object=integrated_seurat, 
                    features = "NPEPPS", 
                    reduction = "umap",
                    cols = c("grey", viridis(n = 10, option = "C", direction = -1)),
                    min.cutoff='q10',
                    pt.size = 0.1,
                    order = TRUE, 
                    label = FALSE,
                    raster = FALSE,
                    combine = TRUE) +
  #NoLegend() +
  my_theme

ggsave("NPEPPS.Umap.tiff")

DotPlot(object=subset(integrated_seurat, celltype=="Epithelial"), 
        features = "NPEPPS", 
        assay = "RNA",
        dot.min = 0,
        dot.scale = 6,
        scale = TRUE,
        scale.by = "size",
        scale.min = 0,
        scale.max = 100) +
  #coord_flip() +
  ggplot2::geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.25) + #stroke is width of circle
  ggplot2::scale_colour_distiller(type = "div", palette = "RdYlGn", direction = -1) +
  ggplot2::guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white", stroke=0.75))) +
  ggplot2::theme(axis.text.x =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5, vjust = 0.5, angle = 45))

ggsave("NPEPPS.Dotplot.tiff")
#******************************************************************************#
#      PERFORM DEGs KRT13 (COLD) vs CDH12 (HOT) IN SIMON SINGLE CELL DATA      #
#          PERFORM DEGs NON-IMMUNE vs IMMUNE IN SIMON SINGLE CELL DATA         #
#******************************************************************************#

proj <- "scRNASeq_Simon"
source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")
source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/results/"

species <- "Homo sapiens"
annotations <- get_annotations(species)

integrated_seurat <- readRDS(paste0(seurat_results, "Simon_integrated_seurat_snn.rds"))

for (r in c("cold_hot", "non.immune_immune")){
  
  if (r == "cold_hot"){
    integrated_seurat@meta.data <- integrated_seurat@meta.data %>%
      dplyr::mutate(Condition = dplyr::case_when(subtype == "KRT13_Epithelial" ~ "KRT13_Epithelial",
                                                 subtype == "CDH12_Epithelial" ~ "CDH12_Epithelial",
                                                 TRUE ~ "Other_cell_types"))
    subtypes <- c("CDH12_Epithelial", "KRT13_Epithelial")
  } else if (r == "non.immune_immune"){
    integrated_seurat@meta.data <- integrated_seurat@meta.data %>%
      dplyr::mutate(Condition = dplyr::case_when(celltype %in% c("Lymphocyte", "Myeloid") ~ "Immune",
                                                 TRUE ~ "Non-Immune"))
    subtypes <- c("Immune", "Non-Immune")
  }
  
  # Prepare read_data and meta_data for each population and merge them
  meta_data_full <- data.frame()
  read_data_full <- data.frame(SYMBOL = rownames(integrated_seurat@assays$RNA$counts))
  
  for (s in subtypes){
    
    subset_seurat <- subset(x=integrated_seurat,
                            subset = Condition == s & !is.na(integrated_seurat@meta.data$Sample))
    
    # Subset metadata
    meta_data <- subset_seurat@meta.data %>%
      dplyr::distinct(Sample, .keep_all=TRUE) %>%
      dplyr::filter(!is.na(Sample)) %>%
      #dplyr::mutate(Batch=1, Condition = s, Sample_ID=paste0(Sample, "_", s)) %>%
      dplyr::mutate(Batch=1, Sample_ID=paste0(Sample, "_", s)) %>%
      dplyr::select(Sample_ID, Batch, Condition)
    
    #******************************IMPORT READ DATA******************************#                                           
    
    # The read data will have "the reads of all cells belonging to a single 
    # sample" merged together in each column. First, create a list of samples
    samples <- subset_seurat@meta.data %>% 
      dplyr::select(Sample) %>% 
      dplyr::filter(!is.na(Sample)) %>%
      unlist(., use.names=FALSE) %>% 
      unique()
    
    # Second, create an empty dataframe with rows=genes and columns=samples
    read_data <- data.frame(matrix(NA, nrow=nrow(subset_seurat@assays$RNA$counts), ncol=nrow(meta_data)))
    rownames(read_data) <- rownames(subset_seurat@assays$RNA$counts)
    colnames(read_data) <- samples
    
    # Thirdly, we will add row-wise, the counts of each gene for each sample
    for(i in samples){
      
      # Create a list of cells for each sample
      cells_subset <- rownames(subset_seurat@meta.data %>% dplyr::filter(Sample == i))
      
      # Use data.frame to convert "." in sparse matrix to "0"
      subset <- data.frame(subset_seurat@assays$RNA$counts[,cells_subset])
      read_data[,i]  <- rowSums(subset)
    }
    
    colnames(read_data) <- paste0(colnames(read_data), "_", s)
    read_data <- read_data %>% 
      tibble::rownames_to_column("SYMBOL")
    
    meta_data_full <- dplyr::bind_rows(meta_data_full, meta_data)
    read_data_full <- dplyr::left_join(read_data_full, read_data, by=c("SYMBOL"="SYMBOL"))
  }
  
  # Rename appropiate column in meta_data to "Sample_ID"
  # meta_data <- meta_data %>% dplyr::rename(Sample_ID = GEO_Accession__exp_)
  
  # Rename appropiate column in read_data to "SYMBOL"
  # read_data <- read_data %>% dplyr::rename(SYMBOL = GeneID)
  # colnames(read_data)[1] <- "SYMBOL"
  
  if (r == "cold_hot"){
    # Define comparisons
    Comparisons <- list(Variable    = c("Condition"),
                        Target      = c("KRT13_Epithelial"),
                        Reference   = c("CDH12_Epithelial"),
                        lfc.cutoff  = 0,
                        padj.cutoff = 0.1)
  } else if (r == "non.immune_immune"){
    # Define comparisons
    Comparisons <- list(Variable    = c("Condition"),
                        Target      = c("Non-Immune"),
                        Reference   = c("Immune"),
                        lfc.cutoff  = 0,
                        padj.cutoff = 0.1)
  }
  
  # Perform QC
  meta_data <- prep_metadata(meta_data_full, read_data_full)
  read_data <- prep_readdata(read_data_full, meta_data)
  l <- check_data(read_data, meta_data)
  meta_data <- l[[2]]
  read_data <- l[[1]]
  
  # Perform analysis for every comparison
  for (n in 1:length(Comparisons$Variable)){
    
    # This generates a new column "id" that has info on samples being comparared
    meta_data_comp <- meta_data %>%
      dplyr::mutate(id=get(Comparisons$Variable[n]))
    
    # Perform DESeq2() using in-built batch modelling
    approach <- "DESeq2"
    if (length(unique(meta_data_comp$Batch)) > 1){
      dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                            colData=meta_data_comp, 
                                            design=~ Batch+id)
    } else {
      dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                            colData=meta_data_comp, 
                                            design=~ id)
    }
    
    prefix <- "zSimon"
    dds <- run_deseq2(dds, meta_data_comp, annotations, Comparisons, n, approach, prefix, output_path)
  }
  
  # Alternative method of DEGs using FindMarkers
  if (r == "cold_hot"){
    subset_seurat <- subset(x=integrated_seurat,
                            subset = Condition %in% c("CDH12_Epithelial", "KRT13_Epithelial"))
  } else if(r == "non.immune_immune"){
    subset_seurat <- subset(x=integrated_seurat,
                            subset = Condition %in% c("Non-Immune", "Immune"))
  }
  
  DefaultAssay(subset_seurat) <- "RNA"
  Idents(object=subset_seurat) <- "Condition"
  
  all_markers <- Seurat::FindAllMarkers(object=subset_seurat,
                                        assay="RNA",
                                        features=NULL,
                                        logfc.threshold=0.25,
                                        test.use="wilcox",
                                        slot="data",
                                        min.pct=0.1,
                                        min.diff.pct=0.1,
                                        node=NULL,
                                        verbose=TRUE,
                                        only.pos=TRUE,
                                        max.cells.per.ident=Inf,
                                        random.seed=1,
                                        latent.vars=NULL,
                                        min.cells.feature=3,
                                        min.cells.group=1,
                                        pseudocount.use=1,
                                        mean.fxn=NULL,
                                        fc.name=NULL,
                                        base=2,
                                        return.thresh=0.01,
                                        densify=FALSE)
  
  
  t1 <- length(intersect( all_markers$gene, annotations$ENSEMBL_SYMBOL))
  t2 <- length(intersect( all_markers$gene, annotations$ENTREZ_SYMBOL))
  if (t1 >= t2){
    all_markers <- all_markers %>% 
      dplyr::mutate(pct.1=dplyr::if_else(pct.1 == 0, 0.001, pct.1),
                    pct.2=dplyr::if_else(pct.2 == 0, 0.001, pct.2),
                    ratio=pct.1/pct.2) %>%
      dplyr::left_join(y=unique(annotations[, c("ENSEMBL_SYMBOL", "CHR", "DESCRIPTION")]), by=c("gene"="ENSEMBL_SYMBOL")) %>%
      dplyr::relocate(cluster, gene, CHR, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, ratio, DESCRIPTION)
  } else {
    all_markers <- all_markers %>% 
      dplyr::mutate(pct.1=dplyr::if_else(pct.1 == 0, 0.001, pct.1),
                    pct.2=dplyr::if_else(pct.2 == 0, 0.001, pct.2),
                    ratio=pct.1/pct.2) %>%
      dplyr::left_join(y=unique(annotations[, c("ENTREZ_SYMBOL", "CHR", "DESCRIPTION")]), by=c("gene"="ENTREZ_SYMBOL")) %>%
      dplyr::relocate(cluster, gene, CHR, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, ratio, DESCRIPTION)
  }
  
  
  if (r == "cold_hot"){
    filename <- "zSimon.DEGs.KRT13_Epithelial_vs_CDH12_Epithelial_Seurat.xlsx"
  } else if (r == "non.immune_immune"){
    filename <- "zSimon.DEGs.Non-Immune_vs_Immune_Seurat.xlsx"
  }
  
  # Save all the markers
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb=wb, sheetName="All_Markers")
  openxlsx::writeData(wb=wb, sheet="All_Markers", x=all_markers)
  openxlsx::saveWorkbook(wb=wb, file=paste0(output_path, filename), overwrite=TRUE)
}

#******************************************************************************#
#                       HIT FILTERING AND IDENTIFICATION                       #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/results/"
supplemental_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/supplemental/"

#**********REACTOME**********#

# Get all reactome pathways from Reactome website
reactome_all_pathways <- read.xlsx(paste0(supplemental_path, "ReactomePathways.gmt.xlsx"),
                                   colNames= FALSE)
colnames(reactome_all_pathways) <- c("Pathway", "ID", paste0("Gene", seq(1, (ncol(reactome_all_pathways)-2))))

# Get all immune related pathway names from Reactome website
# Search immune in reactome website and choose pathways and copy paste all 191 
# pathways from search result into excel file and keep only the R-HSA id
reactome_search_immune <- read.xlsx(paste0(supplemental_path, "Reactome_Immune.xlsx"))

# Filter out non-immune related pathways as defined by reactome
custom_set1 <- reactome_all_pathways %>%
  dplyr::filter(ID %in% reactome_search_immune$ID)

# Create your own filters to identify immune related pathways
my_patterns <- c("Immune","Immuno", "Interferon","Interleukin", "Cytokine", 
                 "Toll", "Complement","Inflamma", "Antigen", "FC",
                 "Co-stimulation", "Co-inhibition")

# Filter out non-immune related pathways using keywords
custom_set2 <- reactome_all_pathways %>%
  dplyr::filter(grepl(pattern= paste0(my_patterns, collapse="|"), x = Pathway, ignore.case = TRUE)) %>%
  dplyr::distinct_at("Pathway", .keep_all = TRUE)

# Merge all immune pathways
reactome_immune_pathways <- dplyr::bind_rows(custom_set1, custom_set2) %>%
  dplyr::distinct_at("ID", .keep_all = TRUE)

# Melt into list of genes
immune_genes <- reactome_immune_pathways %>% 
  dplyr::select(everything(), -c(Pathway,ID)) %>%
  unlist(use.names=FALSE) %>%
  unique()

#**********TCGA**********#
# Get all DEGs that are upregulated in cold tumors
tcga_degs_up <- read.xlsx(paste0(output_path, "zSummary_of_DEGs.xlsx")) %>%
  dplyr::rename(SYMBOL = identity(1)) %>%
  dplyr::filter(n_UP>=10, n_DOWN <=1)

tcga_degs_down <- read.xlsx(paste0(output_path, "zSummary_of_DEGs.xlsx")) %>%
  dplyr::rename(SYMBOL = identity(1)) %>%
  dplyr::filter(n_UP<4, n_DOWN >1)

#**********DGIdb**********#
# Get all druggable genes from DGIdb
druggable_genes <- read.xlsx(paste0(supplemental_path, "DGIdb_categories.xlsx")) %>%
  dplyr::rename(Class = identity(2)) %>%
  dplyr::filter(Class == "DRUGGABLE GENOME")

#****Keep ONLY protein coding genes***********#
annotations <- get_annotations("Homo sapiens")
annotations <- annotations %>% 
  dplyr::filter(ENSEMBL_BIOTYPE == "protein_coding" | ENTREZ_BIOTYPE == "protein-coding")


tcga_degs_up <- tcga_degs_up %>% 
  dplyr::filter(SYMBOL %in% c(annotations$ENSEMBL_SYMBOL, annotations$ENTREZ_SYMBOL, annotations$ENSEMBL_ID))

tcga_degs_down <- tcga_degs_down %>% 
  dplyr::filter(SYMBOL %in% c(annotations$ENSEMBL_SYMBOL, annotations$ENTREZ_SYMBOL, annotations$ENSEMBL_ID))


# #**********Single Cell**********#
# # Get DEGs from single cell analysis
# krt_cdh12_deseq2 <- read.xlsx(paste0(output_path, "zSimon.DEGs.KRT13_Epithelial_vs_CDH12_Epithelial_DESeq2.xlsx")) %>%
#   dplyr::filter(padj < 0.05) %>% 
#   dplyr::select(SYMBOL, log2FoldChange) %>% 
#   dplyr::distinct(.keep_all = TRUE)
# 
# krt_cdh12_seurat <- read.xlsx(paste0(output_path, "zSimon.DEGs.KRT13_Epitheilal_vs_CDH12_Epithelial_Seurat.xlsx")) %>%
#   dplyr::mutate(avg_log2FC = dplyr::case_when(cluster == "CDH12_Epithelial" ~ -(avg_log2FC),
#                                               TRUE ~ avg_log2FC)) %>%
#   dplyr::filter(p_val_adj < 0.05) %>% 
#   dplyr::select(gene, avg_log2FC) %>% 
#   dplyr::distinct(.keep_all = TRUE)
# 
# non_immune_deseq2 <- read.xlsx(paste0(output_path, "zSimon.DEGs.Non-Immune_vs_Immune_DESeq2.xlsx")) %>%
#   dplyr::filter(padj < 0.05) %>% 
#   dplyr::select(SYMBOL, log2FoldChange) %>% 
#   dplyr::distinct(.keep_all = TRUE)
# 
# non_immune_seurat <- read.xlsx(paste0(output_path, "zSimon.DEGs.Non-Immune_vs_Immune_Seurat.xlsx")) %>%
#   dplyr::mutate(avg_log2FC = dplyr::case_when(cluster == "Immune" ~ -(avg_log2FC),
#                                               TRUE ~ avg_log2FC)) %>%
#   dplyr::filter(p_val_adj < 0.05) %>% 
#   dplyr::select(gene, avg_log2FC) %>% 
#   dplyr::distinct(.keep_all = TRUE)
# 
# deseq2_hits <- dplyr::full_join(krt_cdh12_deseq2,non_immune_deseq2, by=c("SYMBOL"="SYMBOL"))
# colnames(deseq2_hits) <- c("SYMBOL","log2FC_KRT13","log2FC_Non_Immune")
# seurat_hits <- dplyr::full_join(krt_cdh12_seurat,non_immune_seurat, by=c("gene"="gene"))
# colnames(seurat_hits) <- c("SYMBOL","log2FC_KRT13","log2FC_Non_Immune")
# 
# # Remove unwanted variables
# rm(custom_set1, custom_set2, krt_cdh12_deseq2, krt_cdh12_seurat, non_immune_deseq2, 
#    non_immune_seurat, reactome_all_pathways, reactome_search_immune, reactome_immune_pathways)
# 
# # Create a dataframe with all filtering so we can color genes appropriately
# deseq2_hits <- deseq2_hits %>%
#   #dplyr::mutate(Pass = dplyr::case_when(SYMBOL %in% tcga_degs_down$SYMBOL ~ "TCGA",
#   dplyr::mutate(Pass = dplyr::case_when(SYMBOL %in% tcga_degs_up$SYMBOL ~ "TCGA",
#                                         TRUE ~ "Others")) %>%
#   dplyr::mutate(Pass = dplyr::case_when(log2FC_KRT13 >= 1 & log2FC_Non_Immune >= 0.25 & Pass == "TCGA" ~ "TCGA+SC",
#                                         TRUE ~ Pass)) %>%
#   dplyr::mutate(Pass = dplyr::case_when(SYMBOL %in% immune_genes & Pass == "TCGA+SC" ~ "TCGA+SC+Reactome",
#                                         TRUE ~ Pass)) %>%
#   dplyr::mutate(Pass = dplyr::case_when(SYMBOL %in% druggable_genes$SYMBOL & Pass == "TCGA+SC+Reactome" ~ "TCGA+SC+Reactome+Druggable",
#                                         TRUE ~ Pass))
# 
# seurat_hits <- seurat_hits %>%
#   #dplyr::mutate(Pass = dplyr::case_when(SYMBOL %in% tcga_degs_down$SYMBOL ~ "TCGA",
#   dplyr::mutate(Pass = dplyr::case_when(SYMBOL %in% tcga_degs_up$SYMBOL ~ "TCGA",
#                                         TRUE ~ "Others")) %>%
#   dplyr::mutate(Pass = dplyr::case_when(log2FC_KRT13 >= 1 & log2FC_Non_Immune >= 0.25 & Pass == "TCGA" ~ "TCGA+SC",
#                                         TRUE ~ Pass)) %>%
#   dplyr::mutate(Pass = dplyr::case_when(SYMBOL %in% immune_genes & Pass == "TCGA+SC" ~ "TCGA+SC+Reactome",
#                                         TRUE ~ Pass)) %>%
#   dplyr::mutate(Pass = dplyr::case_when(SYMBOL %in% druggable_genes$SYMBOL & Pass == "TCGA+SC+Reactome" ~ "TCGA+SC+Reactome+Druggable",
#                                         TRUE ~ Pass))

hits_up <- tcga_degs_up[,1:4] %>%
  dplyr::mutate(Pass = "TCGA",
                Pass = dplyr::case_when(SYMBOL %in% immune_genes & Pass == "TCGA" ~ "TCGA+Reactome",
                                        TRUE ~ Pass),
                Pass = dplyr::case_when(SYMBOL %in% druggable_genes$SYMBOL & Pass == "TCGA+Reactome"~ "TCGA+Reactome+Druggable",
                                        TRUE ~ Pass))


hits_down <- tcga_degs_down[,1:4] %>%
  dplyr::mutate(Pass = "TCGA",
                Pass = dplyr::case_when(SYMBOL %in% immune_genes & Pass == "TCGA" ~ "TCGA+Reactome",
                                        TRUE ~ Pass),
                Pass = dplyr::case_when(SYMBOL %in% druggable_genes$SYMBOL & Pass == "TCGA+Reactome"~ "TCGA+Reactome+Druggable",
                                        TRUE ~ Pass))
# Save as excel
wb <- openxlsx::createWorkbook()

for (d in c("deseq2_hits", "seurat_hits")){
  
  data <- get(d)
  # data <- data %>% dplyr::mutate(size_col = dplyr::case_when(SYMBOL %in% hits$SYMBOL ~ 1.5,
  #                                                             TRUE ~ 0.25))
  
  ggplot(data = data, aes(x=log2FC_KRT13, y=log2FC_Non_Immune, color=Pass)) +
    geom_point(size=0.5) +
    geom_jitter() +
    ggplot2::theme_bw() +
    geom_text_repel(data = data %>% dplyr::filter(Pass == "TCGA+SC+Reactome+Druggable"),
                    mapping = aes(label = SYMBOL),
                    show.legend = FALSE,
                    size = 3,
                    force = 0.5,
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
                    position = position_quasirandom()) +
    ggplot2::labs(x = "High in Immune Enriched     High in Immune Depleted",
                  y = "High in Immune                             High in Non-Immune",
                  title = "",
                  fill = "log2FC") +
    ggplot2::theme(plot.title =   element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                   plot.caption = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0),
                   axis.title.x = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                   axis.title.y = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                   axis.text.x =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                   axis.text.y =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 1),
                   legend.title = element_text(family="sans", face="plain", colour="black", size=10, hjust = 0, vjust = 1),
                   legend.text =  element_text(family="sans", face="plain", colour="black", size=10, hjust = 0.5),
                   #legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
                   legend.position = "right",
                   legend.justification = "left",
                   legend.direction = "vertical",
                   legend.key.height = unit(0.5, 'cm'),
                   legend.key.width = unit(0.5, 'cm')) +
    #viridis::scale_color_viridis(option = "mako") +
    ggplot2::scale_color_manual(labels=c("Others", "TCGA", "TCGA+SC", "TCGA+SC+Reactome", "TCGA+SC+Reactome+Druggable"),
                                #values = c("#440154FF", "#2A788EFF", "#FDE725FF", "#7AD151FF", "#22A884FF")) +
                                values = c("#22A884FF", "#FDE725FF", "#FCA50AFF", "#7A0403FF","#440154FF")) +
    ggplot2::scale_x_continuous(limits = c(-ceiling(max(abs(data$log2FC_KRT13), na.rm=TRUE)),
                                           ceiling(max(abs(data$log2FC_KRT13),na.rm=TRUE))),
                                breaks = seq(-ceiling(max(abs(data$log2FC_KRT13),na.rm=TRUE)),
                                             ceiling(max(abs(data$log2FC_KRT13),na.rm=TRUE)),
                                             by = 1)) +                  # X-axis plot label location
    ggplot2::scale_y_continuous(limits = c(-ceiling(max(abs(data$log2FC_Non_Immune),na.rm=TRUE)),2),
                                #ceiling(max(abs(data$log2FC_Non_Immune)))),
                                breaks = seq(-ceiling(max(abs(data$log2FC_Non_Immune),na.rm=TRUE)),
                                             ceiling(max(abs(data$log2FC_Non_Immune),na.rm=TRUE)),
                                             by = 1)) +
    #ggplot2::scale_size(breaks = sapply(as.vector(quantile(c(min(data$log2FC.1), max(data$log2FC_Non_Immune)))), floor)) +
    geom_hline(yintercept = c(-0.25, 0.25), linetype = "dotted") +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted")
  #scale_size_continuous(range = c(0, 1.5))
  
  ggsave(paste0(output_path, "z", d, "_quadrant_up.tiff"), height = 7, width = 7, units = "in")
  #ggsave(paste0(output_path, "z", d, "_quadrant_down.tiff"), height = 7, width = 7, units = "in")
  data %>% dplyr::count(Pass) %>% dplyr::arrange(n) %>% dplyr::mutate(n= cumsum(n)) %>% dplyr::arrange(desc(n))
  
  openxlsx::addWorksheet(wb, sheetName = d)
  openxlsx::writeData(wb, sheet = d, x = data)
}
openxlsx::saveWorkbook(wb, file = paste0(output_path, "zFinal_Hits_up.xlsx"), overwrite = TRUE)
#openxlsx::saveWorkbook(wb, file = paste0(output_path, "zFinal_Hits_down.xlsx"), overwrite = TRUE)

#******************************************************************************#
#               (ICB DATASETS): SURVIVAL PLOTS & HR VALUES OF HITS             #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")
#source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/results/"

survival_params <- list(plot_by             = c(NA), 
                        split_by            = c(NA),
                        split_plot          = FALSE,
                        multiple_cutoff     = FALSE,
                        stratify_criteria   = c("q"),
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

gse <- c("GSE78220", "GSE91061", "GSE96619", "GSE115821", "GSE131521", "IMvigor210", "IMvigor010")

for (proj in gse){  
  
  meta_data <- read.xlsx(paste0(data_path, proj, ".Metadata.xlsx"))
  read_data <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, proj,".Normalized.counts.xlsx"))
  
  # Reformat metadata 
  meta_data <- meta_data %>% 
    dplyr::mutate(Sample_ID = make.names(names = Sample_ID)) %>%
    dplyr::mutate(Time = as.numeric(Time)) %>%
    dplyr::filter(Time > 0 & !is.na(Time)) %>%
    dplyr::distinct_at("Sample_ID", .keep_all = TRUE)
  
  if (proj == "IMvigor210"){
    meta_data <- meta_data %>% 
      dplyr::filter(Tissue == "bladder", Received.platinum == "Y", `Sample.collected.pre-platinum` == "N")
  }
  if (proj == "IMvigor010"){
    meta_data <- meta_data %>% 
      dplyr::filter(ARM == "Atezolizumab (MPDL3280A) 1200 mg", prior_neoadjuvant_chemotherapy == "YES")
  }
  
  # Reformat read data
  norm_counts <- read_data %>%
    dplyr::mutate(SYMBOL = make.names(names = SYMBOL, unique = TRUE)) %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    tibble::column_to_rownames(var = "SYMBOL") %>%
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
  colnames(norm_counts) <- base::make.names(names = colnames(norm_counts))
  norm_counts <- norm_counts[,intersect(make.names(meta_data$Sample_ID), colnames(norm_counts))]
  norm_counts <- norm_counts[!rowSums(norm_counts, na.rm=TRUE) == 0,]
  
  log_norm_counts <- log(1+norm_counts, base=2)
  t <- base::apply(X=log_norm_counts, MARGIN=1, FUN=median, na.rm=TRUE)
  log_norm_counts <- base::sweep(x=log_norm_counts, MARGIN=1, FUN="-", STATS=t)
  
  # List genes to plot
  hits <- read.xlsx(paste0(output_path, "zFinal_Hits.xlsx"),
                    sheet = "seurat_hits") %>%
    dplyr::filter(stringr::str_detect(string=Pass, pattern="Reactome"))
  plot_genes <- hits$SYMBOL
  plot_genes <- intersect(plot_genes, rownames(log_norm_counts))
  
  # Generate expr_df
  expr_df <- prep_expr_df(log_norm_counts, meta_data, plot_genes, survival_params)
  
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
}

#******************************************************************************#
#        (ImVigor): SURVIVAL PLOTS & HR VALUES OF KRT SCORE, CTL SCORE         #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/results/"

survival_params <- list(plot_by             = c(NA), 
                        split_by            = c(NA),
                        split_plot          = FALSE,
                        multiple_cutoff     = FALSE,
                        stratify_criteria   = c("q"), #m, #o, #t, #q, #th, #none
                        reference           = c("LOW"),
                        conf_interval       = FALSE,
                        plot_curve          = TRUE,
                        plot_risk_table     = TRUE,
                        legend_title        = "Expression",
                        legend_label        = c("High", "Low"),
                        color_palette       = c("#d73027","#0c2c84"),
                        plot_all_bins       = FALSE,
                        gene_sig_score      = TRUE)

gse <- c("IMvigor210", "IMvigor010")

for (proj in gse){  
  
  meta_data <- read.xlsx(paste0(data_path, proj, ".Metadata.xlsx"))
  read_data <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, proj,".Normalized.counts.xlsx"))
  
  # Reformat metadata 
  meta_data <- meta_data %>% 
    dplyr::mutate(Sample_ID = make.names(names = Sample_ID)) %>%
    dplyr::mutate(Time = as.numeric(Time)) %>%
    dplyr::filter(Time > 0 & !is.na(Time)) %>%
    dplyr::distinct_at("Sample_ID", .keep_all = TRUE)
  
  if (proj == "IMvigor210"){
    meta_data <- meta_data %>% 
      dplyr::filter(Tissue == "bladder", Received.platinum == "Y", `Sample.collected.pre-platinum` == "N")
  }
  if (proj == "IMvigor010"){
    meta_data <- meta_data %>% 
      dplyr::filter(ARM == "Atezolizumab (MPDL3280A) 1200 mg", prior_neoadjuvant_chemotherapy == "YES")
  }
  
  # Reformat read data
  norm_counts <- read_data %>%
    dplyr::mutate(SYMBOL = make.names(names = SYMBOL, unique = TRUE)) %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    tibble::column_to_rownames(var = "SYMBOL") %>%
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
  colnames(norm_counts) <- base::make.names(names = colnames(norm_counts))
  norm_counts <- norm_counts[,intersect(make.names(meta_data$Sample_ID), colnames(norm_counts))]
  norm_counts <- norm_counts[!rowSums(norm_counts, na.rm=TRUE) == 0,]
  
  log_norm_counts <- log(1+norm_counts, base=2)
  t <- base::apply(X=log_norm_counts, MARGIN=1, FUN=median, na.rm=TRUE)
  log_norm_counts <- base::sweep(x=log_norm_counts, MARGIN=1, FUN="-", STATS=t)
  
  # List genes to plot
  # krt13_sig <- read.xlsx(paste0(file_path, "41467_2021_25103_MOESM3_ESM.xlsx")) %>%
  #   dplyr::filter(group == "KRT_Epithelial", pvals_adj < 0.05) %>%
  #   dplyr::slice_max(order_by = logfoldchanges, n=200)
  # plot_genes <- krt13_sig$names
  
  plot_genes <- c("IFNG", "GZMA", "GZMB", "PRF1", "GZMK", "ZAP70", "GNLY", "FASLG", "TBX21", "EOMES", "CD8A", "CD8B")
  plot_genes <- intersect(plot_genes, rownames(log_norm_counts))
  
  # Generate expr_df
  expr_df <- prep_expr_df(log_norm_counts, meta_data, plot_genes, survival_params)
  
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
  
  if (survival_params$gene_sig_score == TRUE){
    
    gene <- "combined.exp"
    prefix <- paste0(proj, "_", gene)
    summary <- plot_survival(expr_df, gene, survival_params, prefix, output_path)
    
    # Merge classification info to parent dataframe
    classification_df <- summary[[1]] %>%
      dplyr::select(Sample_ID, combined.exp,	model, everything())
    
    # Merge all stats into a dataframe
    stats_df <- as.data.frame(summary[[2]])
    
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
    openxlsx::saveWorkbook(wb, file = paste0(output_path, proj, "CTL.score_q_Stats.xlsx"), overwrite = TRUE)
  }
}

#******************************************************************************#
#   (ImVigor): SURVIVAL PLOTS & HR VALUES OF IMMUNE ENRICHED/DEPLETED GROUPS   #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/results/"

survival_params <- list(plot_by             = c(NA), 
                        split_by            = c(NA),
                        split_plot          = FALSE,
                        multiple_cutoff     = FALSE,
                        stratify_criteria   = c("none"), #m
                        reference           = c("Immune Enriched"),
                        conf_interval       = FALSE,
                        plot_curve          = TRUE,
                        plot_risk_table     = TRUE,
                        legend_title        = "Immune Status",
                        legend_label        = c("Immune Depleted", "Immune Enriched"),
                        color_palette       = c("#CB181D","#A6D854"),
                        plot_all_bins       = FALSE,
                        plot_all_quartiles  = FALSE,      # If plotting expression and 
                        gene_sig_score      = FALSE)      # If plotting gene signature

meta_data <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"))
proj <- make.names(unique(meta_data$Project_ID))

gse <- c(proj, "GSE78220", "GSE91061", "GSE96619", "GSE115821", "GSE131521", "IMvigor210", "IMvigor010")

for (proj in gse){  
  
  meta_data <- read.xlsx(paste0(data_path, proj, ".Metadata.xlsx"))
  group_data <- read.xlsx(paste0(output_path, proj, ".Scores.&.Heatmap.details.xlsx"),
                          sheet="NPEPPS") %>%
    dplyr::select(Sample, gsva.classification, gsva.mean)
  
  # Merge meta_data and group_data
  meta_data <- meta_data %>%
    dplyr::mutate(Sample_ID = make.names(Sample_ID)) %>%
    dplyr::left_join(group_data, by=c("Sample_ID"="Sample")) %>% 
    dplyr::mutate(model = gsva.classification)
  
  # Reformat metadata 
  meta_data <- meta_data %>% 
    dplyr::mutate(Sample_ID = make.names(names = Sample_ID)) %>%
    dplyr::mutate(Time = as.numeric(Time), Status = as.numeric(Status)) %>%
    dplyr::filter(Time > 0 & !is.na(Time)) %>%
    dplyr::distinct_at("Sample_ID", .keep_all = TRUE)
  
  if (proj == "IMvigor210"){
    meta_data <- meta_data %>% 
      dplyr::filter(Tissue == "bladder", Received.platinum == "Y", `Sample.collected.pre-platinum` == "N")
  }
  if (proj == "IMvigor010"){
    meta_data <- meta_data %>% 
      dplyr::filter(ARM == "Atezolizumab (MPDL3280A) 1200 mg", prior_neoadjuvant_chemotherapy == "YES")
  }
  
  # Generate expr_df
  #expr_df <- prep_expr_df(log_norm_counts, meta_data, plot_genes, survival_params)
  expr_df <- meta_data
  
  gene <- ""
  prefix <- paste0(proj, "_", gene)
  summary <- plot_survival(expr_df, gene, survival_params, prefix, output_path)
  
  # Merge classification info to parent dataframe
  classification_df <- summary[[1]] %>%
    dplyr::select(Sample_ID, model, everything())
  
  # Merge all stats into a dataframe
  stats_df <- as.data.frame(summary[[2]])
  
  # Save the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Summary")
  openxlsx::writeData(wb, sheet = "Summary", x = stats_df)
  openxlsx::addWorksheet(wb, sheetName = "Classification")
  openxlsx::writeData(wb, sheet = "Classification", x = classification_df)
  openxlsx::saveWorkbook(wb, file = paste0(output_path, proj, "Immune_Stats.xlsx"), overwrite = TRUE)
}

#******************************************************************************#
#              STACKED BAR CHART IMMUNE ENRICHED/DEPLETED GROUPS               #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/results/"
#output_path <- "/hpc/home/kailasamms/"

meta_data <- read.xlsx(paste0(data_path, "IMvigor210.Metadata.xlsx"))
group_data <- read.xlsx(paste0(output_path, "IMvigor210.Scores.&.Heatmap.details.xlsx"),
                        sheet="NPEPPS") %>%
  dplyr::select(Sample, gsva.classification, gsva.mean)

meta_data <- meta_data %>%
  dplyr::mutate(Sample_ID = make.names(Sample_ID)) %>%
  dplyr::filter(Tissue == "bladder", Received.platinum == "Y", `Sample.collected.pre-platinum` == "N") %>%
  dplyr::left_join(group_data, by=c("Sample_ID"="Sample")) %>% 
  dplyr::mutate(model = gsva.classification) %>%
  dplyr::group_by(model) %>%
  dplyr::count(Best.Confirmed.Overall.Response) %>%
  dplyr::filter(Best.Confirmed.Overall.Response != "NE") %>%
  dplyr::mutate(percent = round(100*n/sum(n, na.rm=TRUE), digits = 2), label_percent = paste0(percent,"%")) %>%
  data.frame()

# Plot stacked bar chart
ggplot2::ggplot(data = meta_data, aes(x = factor(model), y = percent, fill = Best.Confirmed.Overall.Response)) +
  geom_bar(stat = "identity", width = 0.65) +
  theme_classic() + 
  #scale_fill_manual(values = c("#EDBF33","#A6D854","#DC0AB4","#C7EAE5", "#762A83")) +  
  scale_fill_brewer(palette = "Set1",
                    aesthetics = "fill") +
  geom_text(aes(y=percent, label=label_percent), position=position_stack(vjust=0.5), fontface="bold", colour="white", size=4, check_overlap=TRUE) +
  ggplot2::labs(title = "",
                fill = "Immune Response",
                y = "Percent composition",
                x = "")

ggsave(paste0(output_path, "Stacked_bar_Imvigor210.tiff"))

#******************************************************************************#
#ALL PROTEIN CODING DEGs, CD8A, CD8B, CTL SCORE IN NORMAL & TUMOR TCGA SAMPLES WITH CLASSIFICATION    #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/results/"

# hits <- openxlsx::read.xlsx(xlsxFile=paste0(output_path, "zFinal_Hits.xlsx"),
#                             sheet= "seurat_hits") %>%
#   dplyr::filter(grepl(pattern="Reactome", x=Pass)) %>%
#   dplyr::select(SYMBOL) %>%
#   unlist(use.names=FALSE)
# 
# hits <- c(hits, "CD8A", "CD8B")

# Read the median centered log norm counts
log_norm_counts <- openxlsx::read.xlsx(xlsxFile=paste0(output_path,"TCGA.PanAtlas.median_centered_log_norm_counts.xlsx"))
colnames(log_norm_counts)[1] <- "SYMBOL"

# Get clinical data of normal and tumor samples to identify their TCGA cancer type
meta_data_normal <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"),
                                        sheet = "Normal")
meta_data_tumor <- openxlsx::read.xlsx(xlsxFile=paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"),
                                       sheet = "Tumor")

# Calculate CTL score for each patient
plot_genes <- c("IFNG", "GZMA", "GZMB", "PRF1", "GZMK", "ZAP70", "GNLY", "FASLG", "TBX21", "EOMES", "CD8A", "CD8B")
plot_genes <- intersect(plot_genes, log_norm_counts$SYMBOL)
survival_params <- list(plot_by             = c(NA), 
                        split_by            = c(NA),
                        split_plot          = FALSE,
                        multiple_cutoff     = FALSE,
                        stratify_criteria   = c("none"), #m
                        reference           = c("Immune Enriched"),
                        conf_interval       = FALSE,
                        plot_curve          = TRUE,
                        plot_risk_table     = TRUE,
                        legend_title        = "Immune Status",
                        legend_label        = c("Immune Depleted", "Immune Enriched"),
                        color_palette       = c("#CB181D","#A6D854"),
                        plot_all_bins       = FALSE,
                        plot_all_quartiles  = FALSE,      # If plotting expression and 
                        gene_sig_score      = TRUE)      # If plotting gene signature
expr_df_tumor <- prep_expr_df(log_norm_counts %>% tibble::column_to_rownames("SYMBOL"), meta_data_tumor, plot_genes, survival_params)
expr_df_normal <- prep_expr_df(log_norm_counts %>% tibble::column_to_rownames("SYMBOL"), meta_data_normal, plot_genes, survival_params)
expr_df_tumor <- expr_df_tumor %>%
  dplyr::rename(CTL.score=combined.exp) %>%
  dplyr::select(Sample_ID, CTL.score)
expr_df_normal <- expr_df_normal %>%
  dplyr::rename(CTL.score=combined.exp) %>%
  dplyr::select(Sample_ID, CTL.score)

# Get classification of samples into immune enriched vs depleted
classification_df <- openxlsx::read.xlsx(xlsxFile=paste0(output_path, "Two_Group_Classification_from_Individual_cancers.xlsx"))
lower_cutoff <- quantile(classification_df$Mean, c(0.33,0.66))[[1]]
upper_cutoff <- quantile(classification_df$Mean, c(0.33,0.66))[[2]]

# Get all DEGs nUP>10, ndown<=1
DEGs_df <- openxlsx::read.xlsx(xlsxFile=paste0(output_path, "zSummary_of_DEGs.xlsx"))
colnames(DEGs_df)[1] <- "SYMBOL"
DEGs_df <- DEGs_df %>%
  dplyr::filter(n_UP > 10, n_DOWN <=1) %>%
  dplyr::select(SYMBOL) %>%
  unlist(use.names=FALSE)

# Keep ONLY protein coding DEGs
annotations <- get_annotations()
protein_coding <- annotations[[1]] %>% dplyr::filter(str_detect(string=ENSEMBL_BIOTYPE, pattern="protein") |
                                                       str_detect(string=ENTREZ_BIOTYPE,pattern="protein"))
DEGs_df <- intersect(DEGs_df, union(protein_coding$ENSEMBL_SYMBOL, protein_coding$ENTREZ_SYMBOL))
DEGs_df <- c(DEGs_df, "CD8A", "CD8B")

# Merge all info together
df_normal <- log_norm_counts %>% 
  dplyr::filter(SYMBOL %in% DEGs_df) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  t() %>%
  data.frame() %>% 
  tibble::rownames_to_column("Sample_ID") %>%
  dplyr::inner_join(meta_data_normal %>% dplyr::select(Sample_ID, Project_ID), by=c("Sample_ID"="Sample_ID")) %>%
  dplyr::inner_join(expr_df_normal) %>%
  dplyr::select(Project_ID, Sample_ID, CTL.score, everything())

df_tumor <- log_norm_counts %>% 
  dplyr::filter(SYMBOL %in% DEGs_df) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  t() %>%
  data.frame() %>% 
  tibble::rownames_to_column("Sample_ID") %>%
  dplyr::inner_join(meta_data_tumor %>% dplyr::select(Sample_ID, Project_ID), by=c("Sample_ID"="Sample_ID")) %>%
  dplyr::left_join(classification_df %>% dplyr::select(Sample, Class, Mean), by=c("Sample_ID"="Sample")) %>%
  dplyr::inner_join(expr_df_tumor) %>%
  dplyr::mutate(Class.33 = dplyr::case_when(Mean <= lower_cutoff ~ Class,
                                            Mean >= upper_cutoff ~ Class,
                                            TRUE ~ "")) %>%
  dplyr::select(Project_ID, Sample_ID, Class, Class.33, Mean, CTL.score, everything())

# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Normal")
openxlsx::writeData(wb, sheet = "Normal", x = df_normal)
openxlsx::addWorksheet(wb, sheetName = "Tumor")
openxlsx::writeData(wb, sheet = "Tumor", x = df_tumor)
openxlsx::saveWorkbook(wb, file = paste0(output_path, "zHit.Expression.Tumor.Normal.xlsx"), overwrite = TRUE)

# Calculate correlation between each gene and immune score; each gene vs CTL score
proj <- c(unique(meta_data_tumor$Project_ID))
DEGs <- intersect(make.names(DEGs_df), make.names(colnames(df_tumor)))
df <- data.frame()

for (f in proj){
  
  df_tumor_subset <- df_tumor %>%
    dplyr::filter(Project_ID == f)
  
  gene <- c()
  spearman.p.immune <- c()
  spearman.r.immune <- c()
  pearson.p.immune <- c()
  pearson.r.immune <- c()
  spearman.p.ctl <- c()
  spearman.r.ctl <- c()
  pearson.p.ctl <- c()
  pearson.r.ctl <- c()
  spearman.p.immune.33 <- c()
  spearman.r.immune.33 <- c()
  pearson.p.immune.33 <- c()
  pearson.r.immune.33 <- c()
  spearman.p.ctl.33 <- c()
  spearman.r.ctl.33 <- c()
  pearson.p.ctl.33 <- c()
  pearson.r.ctl.33 <- c()
  
  for (i in 1:length(DEGs)){
    
    gene <- c(gene, DEGs[i])
    cat(i,":", DEGs[i], "\n")
    
    # Correlation between gene z score and immune score
    immune.df <- df_tumor_subset %>%
      dplyr::select(Mean, all_of(DEGs[i]))
    
    res <- cor.test(x=immune.df[[1]], y=immune.df[[2]], method = "spearman")
    spearman.p.immune <- c(spearman.p.immune, res$p.value)
    spearman.r.immune <- c(spearman.r.immune, res$estimate)
    
    res <- cor.test(x=immune.df[[1]], y=immune.df[[2]], method = "pearson")
    pearson.p.immune <- c(pearson.p.immune, res$p.value)
    pearson.r.immune <- c(pearson.r.immune, res$estimate)
    
    # Correlation between gene z score and CTL score
    ctl.df <- df_tumor_subset %>%
      dplyr::select(CTL.score, all_of(DEGs[i]))
    
    res <- cor.test(x=ctl.df[[1]], y=ctl.df[[2]], method = "spearman")
    spearman.p.ctl <- c(spearman.p.ctl, res$p.value)
    spearman.r.ctl <- c(spearman.r.ctl, res$estimate)
    
    res <- cor.test(x=ctl.df[[1]], y=ctl.df[[2]], method = "pearson")
    pearson.p.ctl <- c(pearson.p.ctl, res$p.value)
    pearson.r.ctl <- c(pearson.r.ctl, res$estimate)
    
    # Correlation between gene z score and immune score using top 33% and bottom 33%
    immune.df <- df_tumor_subset %>%
      dplyr::filter(nchar(Class.33)!=0) %>% 
      dplyr::select(Mean, all_of(DEGs[i]))
    
    res <- cor.test(x=immune.df[[1]], y=immune.df[[2]], method = "spearman")
    spearman.p.immune.33 <- c(spearman.p.immune.33, res$p.value)
    spearman.r.immune.33 <- c( spearman.r.immune.33, res$estimate)
    
    res <- cor.test(x=immune.df[[1]], y=immune.df[[2]], method = "pearson")
    pearson.p.immune.33 <- c(pearson.p.immune.33, res$p.value)
    pearson.r.immune.33 <- c(pearson.r.immune.33, res$estimate)
    
    # Correlation between gene z score and ctl score using top 33% and bottom 33%
    ctl.df <- df_tumor_subset %>%
      dplyr::filter(nchar(Class.33)!=0) %>% 
      dplyr::select(CTL.score, all_of(DEGs[i]))
    
    res <- cor.test(x=ctl.df[[1]], y=ctl.df[[2]], method = "spearman")
    spearman.p.ctl.33 <- c(spearman.p.ctl.33, res$p.value)
    spearman.r.ctl.33 <- c(spearman.r.ctl.33, res$estimate)
    
    res <- cor.test(x=ctl.df[[1]], y=ctl.df[[2]], method = "pearson")
    pearson.p.ctl.33 <- c(pearson.p.ctl.33, res$p.value)
    pearson.r.ctl.33 <- c(pearson.r.ctl.33, res$estimate)
    
  }
  summary <- data.frame("Cancer" = f,
                        "Gene" = gene,
                        "spearman.p.immune" = spearman.p.immune,
                        "spearman.r.immune" = spearman.r.immune,
                        "pearson.p.immune" = pearson.p.immune,
                        "pearson.r.immune" = pearson.r.immune,
                        "spearman.p.ctl" = spearman.p.ctl,
                        "spearman.r.ctl" = spearman.r.ctl,
                        "pearson.p.ctl" = pearson.p.ctl,
                        "pearson.r.ctl" = pearson.r.ctl,
                        "spearman.p.immune.33" = spearman.p.immune.33, 
                        "spearman.r.immune.33" = spearman.r.immune.33, 
                        "pearson.p.immune.33" = pearson.p.immune.33,
                        "pearson.r.immune.33" = pearson.r.immune.33,
                        "spearman.p.ctl.33" = spearman.p.ctl.33,
                        "spearman.r.ctl.33" = spearman.r.ctl.33,
                        "pearson.p.ctl.33" = pearson.p.ctl.33,
                        "pearson.r.ctl.33" = pearson.r.ctl.33)
  
  df <- dplyr::bind_rows(df, summary)
}

# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Summary")
openxlsx::writeData(wb, sheet = "Summary", x = df)
openxlsx::saveWorkbook(wb, file = paste0(output_path, "Correlation.xlsx"), overwrite = TRUE)



#******************************************************************************#
#  ORA ANALYSIS oF UP AND DOWN REGULATED GENES IN ID vs IE IN EACH TCGA CANCER #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

file_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/original/"
data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/results/"

species <- "Homo sapiens"

# Read the gene set files
species <- "Human"
gmt_dir <- paste0("/hpc/home/kailasamms/projects/RNASeq/GSEA_genesets/", species)
gmt_files <- list.files(gmt_dir, full.names = TRUE)

# Create dummy dataframes to store pathway results from all cancers
up_df <- data.frame()
down_df <- data.frame()

# Read metadata for tumor samples ONLY
metadata <- read.xlsx(paste0(data_path, "TCGA.PanAtlas.Metadata.xlsx"),
                      sheet = "Tumor")

for (proj in make.names(unique(metadata$Project_ID))){
  
  # Read DEGs
  df <- read.xlsx(paste0(output_path, proj, ".DEGs.gsva.classification_Immune Depleted_vs_Immune Enriched_DESeq2_.xlsx"))
  
  # Define input genes (ONLY significant genes)
  input_genes_up <- df %>% 
    dplyr::filter(padj <= 0.05, log2FoldChange > 0.58) %>%
    dplyr::select(ENSEMBL_SYMBOL) %>%
    unlist(use.names=FALSE)
  
  input_genes_down <- df %>% 
    dplyr::filter(padj <= 0.05, log2FoldChange < -0.58) %>%
    dplyr::select(ENSEMBL_SYMBOL) %>%
    unlist(use.names=FALSE)
  
  # Define universe genes
  universe_genes <- df %>% 
    dplyr::select(ENSEMBL_SYMBOL) %>%
    unlist(use.names=FALSE)
  
  # Perform ORA analysis
  ora_result_up <- ora(input_genes_up, universe_genes, gmt_files)
  ora_result_down <- ora(input_genes_down, universe_genes, gmt_files)
  
  # Merge ora results from each cancer
  ora_result_up <- ora_result_up %>% 
    dplyr::mutate(ProjectID = proj)
  ora_result_down <- ora_result_down %>% 
    dplyr::mutate(ProjectID = proj)
  up_df <- dplyr::bind_rows(up_df, ora_result_up)
  down_df <- dplyr::bind_rows(down_df, ora_result_down)
  
  print(proj)
  print(dim(up_df))
  print(dim(down_df))
}

# Calculate number of cancers where pathway is up or down
up_df_sig <- up_df %>% 
  dplyr::filter(p.adjust <= 0.05) %>% 
  dplyr::add_count(Description)

down_df_sig <- down_df %>% 
  dplyr::filter(p.adjust <= 0.05) %>% 
  dplyr::add_count(Description)

common <- intersect(up_df_unique$Description, down_df_unique$Description)

up_df_unique <- up_df_sig %>% 
  dplyr::filter(!(Description %in% common)) %>%
  dplyr::select(n, ProjectID, Description, everything()) %>%
  dplyr::arrange(desc(n), Description)

down_df_unique <- down_df_sig %>% 
  dplyr::filter(!(Description %in% common)) %>%
  dplyr::select(n, ProjectID, Description, everything()) %>%
  dplyr::arrange(desc(n), Description)

# Create workbook to store results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName=paste0("ORA_UP_UNIQUE"))
openxlsx::writeData(wb, sheet=paste0("ORA_UP_UNIQUE"), x=up_df_unique, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName=paste0("ORA_DOWN_UNIQUE"))
openxlsx::writeData(wb, sheet=paste0("ORA_DOWN_UNIQUE"), x=down_df_unique, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName=paste0("ORA_UP_ALL"))
openxlsx::writeData(wb, sheet=paste0("ORA_UP_ALL"), x=up_df, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName=paste0("ORA_DOWN_ALL"))
openxlsx::writeData(wb, sheet=paste0("ORA_DOWN_ALL"), x=down_df, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(output_path, "zPanCancer_Pathway_Analysis_Results.xlsx"),
                       overwrite = TRUE)




# NOTE: You can convert the tsv to xlsx and also import into R
# identical(df1, df2) will be FALSE indicating there is difference in dataframes
# imported using read.table() vs read.xlsx().
# setdiff(df1,df2) will reveal no differences.
# However, sapply(df, class) will reveal that "size" column is stored as numeric
# by read.xlsx() but as integer by read.table().