# library(factoextra)
# library(igraph)
# library(mlbench)    # for Ionosphere data
# library(psych)      # for cor2dist
# data(Ionosphere)
# 
# scores1 <- scores[rownames(scores) %in% c(group_II, group_III),]
# t1 <- kmeans(t(scores1), centers=2, nstart=50)
# fviz_cluster(t1, t(scores1), labelsize=0)

# col_annotation <- scores1 %>%
#   data.frame() %>%
#   tibble::rownames_to_column("Pathway") %>%
#   dplyr::mutate(Group = dplyr::case_when(Pathway %in% group_I ~ "Group_I",
#                                          Pathway %in% group_II ~ "Group_II",
#                                          Pathway %in% group_III ~ "Group_III",
#                                          TRUE ~ "Group_IV")) %>%
#   tibble::column_to_rownames("Pathway") %>%
#   dplyr::group_by(Group) %>%
#   dplyr::summarise(across(.cols=everything(), .fns=mean)) %>%
#   dplyr::ungroup() %>%
#   tibble::column_to_rownames("Group") %>%
#   t() %>%
#   data.frame() %>%
#   dplyr::mutate(Class = dplyr::case_when(Group_I > 0 & Group_II > 0 & Group_III > 0 ~ "Immune Enriched Fibrotic",
#                                          Group_I < 0 & Group_II > 0 & Group_III > 0 ~ "Immune Enriched Non-Fibrotic",
#                                          Group_I > 0 & Group_II < 0 & Group_III < 0 ~ "Fibrotic",
#                                          Group_I < 0 & Group_II < 0 & Group_III < 0 ~ "Depleted",
#                                          TRUE ~ "Unclassified"),
#                 Dummy_col = NA) %>%  # ordering by rownames doesnt work for 1 column df
#   dplyr::select(Class, Dummy_col)

#*******Heatmap of individual cancers, pancancer, Imvigor010, Imvigor210*******#

#source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

parent_path <- "/hpc/home/kailasamms/scratch/TCGA_GDC/"
results_path <- "/hpc/home/kailasamms/scratch/proj_Saswat/"

meta_data <- read.xlsx(paste0(parent_path, "Meta_data_TCGA.xlsx"))
proj <- unique(meta_data$Project_ID)
proj <- c(proj, "IMvigor010", "IMvigor210", "PanCancer")

# Read signature gene sets for 29 pathways
sig <- read.xlsx(paste0(results_path, "Saswat_Sig_list.xlsx"))
sig <- sig[-1, -1] %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames(colnames(.)[1]) %>% 
  t() %>%
  data.frame()

# List the 4 groups 
group_I <- c("Angiogenesis","Endothelium","Cancer.associated.fibroblasts",
             "Matrix", "Matrix.remodeling", "Protumor.cytokines",
             "Neutrophil.signature", "Granulocyte.traffic")
group_II <- c("Tumor.associated.Macrophages", "Macrophage.and.DC.traffic",
              "Myeloid.cells.traffic","Immune.Suppression.by.Myeloid.Cells",
              "Th2.signature","Treg.and.Th2.traffic", "Treg", "M1.signature")
group_III <- c("MHCII", "Antitumor.cytokines", "Co.activation.molecules",
               "B.cells", "NK.cells", "Checkpoint.molecules", "Effector.cells",
               "T.cells", "Th1.signature", "Effector.cell.traffic", "MHCI")
group_IV <- c("EMT.signature", "Tumor.proliferation.rate")

# Keep pathways in same order as in the paper
ordering <- c("Angiogenesis","Endothelium","Cancer.associated.fibroblasts",
              "Matrix", "Matrix.remodeling", "Protumor.cytokines",
              "Neutrophil.signature", "Granulocyte.traffic",
              "Tumor.associated.Macrophages", "Macrophage.and.DC.traffic",
              "Myeloid.cells.traffic","Immune.Suppression.by.Myeloid.Cells",
              "Th2.signature","Treg.and.Th2.traffic", "Treg", "M1.signature",
              "MHCII", "Antitumor.cytokines", "Co.activation.molecules",
              "B.cells", "NK.cells", "Checkpoint.molecules", "Effector.cells",
              "T.cells", "Th1.signature", "Effector.cell.traffic", "MHCI",
              "EMT.signature", "Tumor.proliferation.rate")

# Keep only group II and group III as Immune Enriched and Immune Depleted are
# defined based on the pathways in these 2 groups
sig <- sig %>% dplyr::select(all_of(group_II), all_of(group_III))

for (f in proj){
  
  if (f == "PanCancer"){
    read_data <- read.table(paste0(parent_path, "Normalized_Counts_TCGA.tsv"), header=TRUE)
  } else if(f %in% c("IMvigor010", "IMvigor210")){
    read_data <- read.xlsx(paste0(results_path, "Normalized_Counts_", f, ".xlsx"))
    read_data <- read_data[,-c(2:4)]
  } else{
    read_data <- read.xlsx(paste0(parent_path, "Normalized_Counts_DESeq2_modelled", f, ".xlsx"))
    read_data <- read_data[,-c(2:5)]
  }
  
  read_data <- read_data %>% 
    base::replace(is.na(.), 0)                           # replace NA with 0
  read_data <- read_data[rowSums(read_data[,-1]) !=0,]   # remove genes with no expression in all samples
  colnames(read_data)[1] <- "SYMBOL"
  
  # log normalize and median center the counts
  normalized_counts <- read_data %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(across(.cols=everything(), .fns=mean)) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("SYMBOL")
  normalized_counts <- log(1+normalized_counts, base=2)
  t <- base::apply(X=normalized_counts, MARGIN=1, FUN=median, na.rm=TRUE)
  normalized_counts <- base::sweep(x=normalized_counts, MARGIN=1, FUN="-", STATS=t)
  
  if (f == "PanCancer"){
    # For Pancancer, save median centered log normalized counts as this step is time consuming
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "median centered log normalized")
    openxlsx::writeData(wb, sheet = "median centered log normalized", x = normalized_counts)
    openxlsx::saveWorkbook(wb, file = paste0(results_path, "TCGA_median_centered_log_normalized_counts.xlsx"), overwrite = TRUE)
  }
  
  # Create empty dataframe of samples as rownames and pathway as column names
  scores <- normalized_counts %>% 
    t() %>%
    data.frame() %>%
    dplyr::select(colnames(.)[1]) %>%
    tibble::rownames_to_column("Sample")
  
  # Calculate scores for each patient for each of the 29 signatures
  for (i in 1:ncol(sig)){
    
    plot_genes <- sig[,i][!is.na(sig[,i])] 
    
    expr_df <- as.data.frame(advanced_Z(plot_genes, normalized_counts))
    colnames(expr_df) <- colnames(sig)[i]
    expr_df <- expr_df %>%
      tibble::rownames_to_column("Sample")
    
    scores <- dplyr::left_join(scores, expr_df, by=c("Sample"="Sample"))
  }
  
  scores <- scores[,-2]
  scores <- scores %>% 
    tibble::column_to_rownames("Sample") %>%
    t()
  
  # scores <- scores[ordering,]
  rownames(scores) <- make.names(rownames(scores))
  colnames(scores) <- make.names(colnames(scores))
  
  if (f == "PanCancer"){
    # For Pancancer, save the transposed scores for each of the 29 pathways for each patient
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "scores")
    openxlsx::writeData(wb, sheet = "scores", x = t(scores), rowNames = TRUE)
    openxlsx::saveWorkbook(wb, file = paste0(results_path, "Pathway_29_TCGA_scores.xlsx"), overwrite = TRUE)
  }
  
  # # 4 way classification
  # col_annotation <- scores %>%
  #   data.frame() %>%
  #   tibble::rownames_to_column("Pathway") %>%
  #   dplyr::mutate(Group = dplyr::case_when(Pathway %in% group_I ~ "Group_I",
  #                                          Pathway %in% group_II ~ "Group_II",
  #                                          Pathway %in% group_III ~ "Group_III",
  #                                          TRUE ~ "Group_IV")) %>%
  #   tibble::column_to_rownames("Pathway") %>%
  #   dplyr::group_by(Group) %>%
  #   dplyr::summarise(across(.cols=everything(), .fns=mean)) %>%
  #   dplyr::ungroup() %>%
  #   tibble::column_to_rownames("Group") %>%
  #   t() %>%
  #   data.frame() %>%
  #   dplyr::mutate(Class = dplyr::case_when(Group_I > 0 & Group_II > 0 & Group_III > 0 ~ "Immune Enriched Fibrotic",
  #                                          Group_I < 0 & Group_II > 0 & Group_III > 0 ~ "Immune Enriched Non-Fibrotic",
  #                                          Group_I > 0 & Group_II < 0 & Group_III < 0 ~ "Fibrotic",
  #                                          Group_I < 0 & Group_II < 0 & Group_III < 0 ~ "Depleted",
  #                                          TRUE ~ "Unclassified"),
  #                 Dummy_col = NA) %>%  # ordering by rownames doesnt work for 1 column df
  #   dplyr::select(Class, Dummy_col)
  
  # 2 way classification
  scores <- scores[setdiff(rownames(scores), c(group_I, group_IV)),]
  col_annotation <- scores %>%
    t() %>%
    data.frame() %>%
    dplyr::mutate(Mean = rowMeans(as.matrix(.)),
                  Class = dplyr::case_when(Mean > 0 ~ "Immune Enriched",
                                           TRUE ~ "Immune Depleted")) %>%
    dplyr::select(Class, Mean)
  
  # Prepare row annotation for heatmap
  row_annotation <- data.frame(Signature=rownames(scores),
                               Group=c(#rep(x="Angiogenesis/Fibroblasts",times=8),
                                 rep(x="Pro-tumor Immune infiltrate",times=8),
                                 rep(x="Anti-tumor Immune infiltrate",times=11)),
                               #rep(x="EMT signature/Proliferation rate",times=2)),
                               Dummy_col = NA) %>%
    tibble::column_to_rownames("Signature")
  
  # # Within group column clustering
  # col_order <- c()
  # element_names <- sort(unique(col_annotation$Class))
  # for (g in element_names){
  #   temp_mat <- scores[,rownames(col_annotation)[col_annotation$Class == g]]
  #   colclust <- hclust(dist(t(temp_mat)))
  #   col_order <- c(col_order, colnames(temp_mat[,colclust$order]))
  # }
  
  col_order <- rownames(col_annotation %>% dplyr::arrange(Mean))
  
  # Within group row clustering
  row_order <- c()
  element_names <- sort(unique(row_annotation$Group))
  for (g in element_names){
    temp_mat <- scores[rownames(row_annotation)[row_annotation$Group == g],]
    rowclust <- hclust(dist(temp_mat))
    row_order <- c(row_order,rownames(temp_mat[rowclust$order,]))
  }
  reordered <- scores[row_order,col_order]
  col_annotation <- col_annotation[colnames(reordered),]
  row_annotation <- row_annotation[rownames(reordered),]
  
  # Get NPEPPS expression for each patient
  ctl_score <- reordered %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::select("Sample", "Effector.cells")
  
  npepps_zscore <- normalized_counts %>% 
    tibble::rownames_to_column("SYMBOL") %>%
    dplyr::filter(SYMBOL == "NPEPPS") %>% 
    tibble::column_to_rownames("SYMBOL") %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::mutate(Sample = make.names(Sample))
  
  npepps <- read_data %>% 
    dplyr::filter(SYMBOL == "NPEPPS") %>% 
    tibble::column_to_rownames("SYMBOL") %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    dplyr::mutate(Sample = make.names(Sample)) %>%
    dplyr::left_join(col_annotation %>% tibble::rownames_to_column("Sample"),
                     by=c("Sample"="Sample")) %>%
    dplyr::left_join(npepps_zscore, by=c("Sample"="Sample")) %>%
    dplyr::left_join(ctl_score, by=c("Sample"="Sample"))
  
  low_cutoff <- quantile(npepps$Mean, c(0.33))[[1]]
  high_cutoff <- quantile(npepps$Mean, c(0.66))[[1]]
  
  npepps <- npepps %>% 
    dplyr::mutate(Category = dplyr::case_when(Mean <= low_cutoff | Mean >= high_cutoff ~ Class,
                                              TRUE ~ NA))
  
  # Save heatmap details and NPEPPS expression
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Heatmap_matrix")
  openxlsx::writeData(wb, sheet = "Heatmap_matrix", x = reordered, rowNames = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "col_annotation")
  openxlsx::writeData(wb, sheet = "col_annotation", x = col_annotation, rowNames = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "row_annotation")
  openxlsx::writeData(wb, sheet = "row_annotation", x = row_annotation, rowNames = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "NPEPPS")
  openxlsx::writeData(wb, sheet = "NPEPPS", x = npepps, rowNames = TRUE)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, "Heatmap_details_", f, ".xlsx"), overwrite = TRUE)
  
  gaps_col <- sum(col_annotation$Class == "Immune Depleted")  #c(168,255)
  gaps_row <- sum(row_annotation$Group == "Anti-tumor Immune infiltrate") #c(8,16,27)
  mat <- scores
  
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
  
  # Plot heatmap
  pheatmap(reordered,
           scale = "row",
           border_color = "white",
           breaks = breaks,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           annotation_col = col_annotation %>% dplyr::select(Class),
           annotation_row = row_annotation %>% dplyr::select(Group),
           gaps_col = gaps_col,
           gaps_row = gaps_row,
           show_colnames=FALSE,
           show_rownames=TRUE,
           annotation_names_row = FALSE,
           annotation_names_col = FALSE,
           width = 20,
           filename = paste0(results_path, "Heatmap_", f, ".tiff"))
}

#**Perform DESeq2 between Immune enriched vs depleted patients in each cancer**#

#source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

parent_path <- "/hpc/home/kailasamms/scratch/proj_Saswat/"
results_path <- "/hpc/home/kailasamms/scratch/proj_Saswat/"

species <- "Homo sapiens"
padj.cutoff <- 0.1
lfc.cutoff <- 0
heatmap_plot <- FALSE
volcano_plot <- FALSE
cor_plot  <- FALSE

annotations <- get_annotations(species)
meta_data <- read.xlsx(paste0(parent_path, "Meta_data_TCGA.xlsx"))

proj <- unique(meta_data$Project_ID)
#proj <- c(proj, "PanCancer")

for (f in proj){
  
  meta_data <- read.xlsx(paste0(parent_path, "Meta_data_", f, ".xlsx"))
  read_data <- read.xlsx(paste0(parent_path, "Read_data_", f, ".xlsx"))
  read_data <- read_data[,-c(2)]
  
  read_data <- read_data %>% 
    base::replace(is.na(.), 0)                           # replace NA with 0
  read_data <- read_data[rowSums(read_data[,-1]) !=0,]   # remove genes with no expression in all samples
  colnames(read_data)[1] <- "SYMBOL"
  
  group_data <- read.xlsx(paste0(parent_path, "2 way classification mean 0/heatmaps 2 way mean cutoff 0/Heatmap_details_", f, ".xlsx"),
                          sheet="col_annotation")
  colnames(group_data) <- c("Sample_ID", "Class", "Mean")
  
  meta_data <- meta_data %>%
    dplyr::mutate(Sample_ID = make.names(Sample_ID)) %>%
    dplyr::left_join(group_data, by=c("Sample_ID"="Sample_ID"))
  
  Comparisons <- list(Variable=c("Class"),
                      Target=c("Immune Depleted"),
                      Reference=c("Immune Enriched"))
  suffix <- f
  meta_data <- prep_metadata(meta_data, Variable)
  read_data <- prep_readdata(read_data, meta_data)
  l <- check_data(read_data, meta_data)
  meta_data <- l[[2]]
  read_data <- l[[1]]
  
  for (n in 1:length(Comparisons$Target)){
    
    # This generates a new column "id" that has info on samples being comparared
    meta_data_comp <- meta_data %>%
      dplyr::mutate(id=get(Comparisons$Variable[n]))
    
    # Perform DESeq2() using in-built batch modelling
    approach <- "DESeq2_modelled"
    if (length(unique(meta_data_comp$Batch)) > 1){
      dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                            colData=meta_data_comp, 
                                            design=~ Batch+id)
    } else {
      dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                            colData=meta_data_comp, 
                                            design=~ id)
    }
    dds <- run_deseq2(dds, meta_data_comp, annotations, Comparisons, n, lfc.cutoff, padj.cutoff, approach, suffix)
    plot_qc(dds, meta_data_comp, approach, suffix)
  }
}

############ Read all DEG results and find frequency for each gene

meta_data <- read.xlsx(paste0(parent_path, "Meta_data_TCGA.xlsx"))
proj <- unique(meta_data$Project_ID)

deg_genes <- c()
for (f in proj){
  degs <- read.xlsx(paste0(parent_path, "2 way classification mean 0/degs/Results_id_Immune Depleted_vs_Immune Enriched_DESeq2_modelled_", f, "_DEGs.xlsx")) %>%
    dplyr::filter(padj <=0.05)
  deg_genes <- c(deg_genes, degs$SYMBOL)
}

deg_genes <- unique(deg_genes)

df <- data.frame(SYMBOL = deg_genes, 
                 DUMMY = 0)
for (f in proj){
  
  degs <- read.xlsx(paste0(parent_path, "2 way classification mean 0/degs/Results_id_Immune Depleted_vs_Immune Enriched_DESeq2_modelled_", f, "_DEGs.xlsx")) %>%
    dplyr::filter(padj <=0.05) %>%
    dplyr::select(SYMBOL, log2FoldChange) %>%
    dplyr::rename(!!rlang::sym(f) := "log2FoldChange") %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE)
  
  df <- df %>% 
    dplyr::left_join(degs, by=c("SYMBOL"="SYMBOL"))
}  

df <- df %>%
  dplyr::select(c(everything(), -c("DUMMY"))) %>%
  tibble::column_to_rownames("SYMBOL")

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
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Summary")
openxlsx::writeData(wb, sheet = "Summary", x = df, rowNames = TRUE)
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Summary.xlsx"), overwrite = TRUE)

#*****************Plot KM curves *******************

parent_path <- "/hpc/home/kailasamms/scratch/proj_Saswat/"
results_path <- "/hpc/home/kailasamms/scratch/proj_Saswat/"

subset_group <- NA    
subset_value <- NA

plot_by <- "Class"
split_by <- subset_group
combine_plot <- FALSE 
multiple_cutoff <- TRUE
stratify_criteria <- "o"
reference <- "Immune Enriched"
reference <- dplyr::if_else(plot_by == "Expression", "LOW", reference)
confidence_interval <- FALSE
legend_title <- "Class"
plot_risk_table <- TRUE
plot_curve <- TRUE
all_quartiles <- FALSE
gene_signature <- FALSE
legend_label <- c("Immune Depleted", "Immune Enriched")
color_palette <- c("#d73027","#0c2c84")
variable_x <- "Time"
variable_y <- "Status" 
gene <- ""
group <- ""

meta_data <- read.xlsx(paste0(parent_path, "Meta_data_TCGA.xlsx"))
proj <- c(unique(meta_data$Project_ID), "IMvigor010", "IMvigor010", "IMvigor210")

# Create a list to store survminer cutoffs, coxph stats, etc..
stats <- list("gene" = c(),
              "HR" = c(),
              "CI_lower" = c(),
              "CI_upper" = c(),
              "pvalue" = c())

n <- 0
for (f in proj){
  
  meta_data <- read.xlsx(paste0(parent_path, "Meta_data_", f, ".xlsx"))
  
  group_data <- read.xlsx(paste0(parent_path, "2 way classification mean 0/Heatmap_details_", f, ".xlsx"),
                          sheet="col_annotation")
  # group_data <- read.xlsx(paste0(parent_path, "4 way classification/Heatmap_matrix_", f, ".xlsx"),
  #                         sheet="col_annotation")
  
  colnames(group_data) <- c("Sample_ID", "Class", "Mean")
  
  group_data <- group_data %>%
    dplyr::mutate(Class = dplyr::case_when(grepl("Immune Enriched",Class) ~ "Immune Enriched",
                                           grepl("Depleted|Fibrotic",Class) ~ "Immune Depleted",
                                           TRUE ~ "Unclassified"))
  suffix <- ""
  
  # FOr Imvigor210 and Imvigor210_new
  if (f == "IMvigor210"){
    meta_data <- meta_data %>% 
      dplyr::filter(Tissue == "bladder", Received.platinum == "Y", `Sample.collected.pre-platinum` == "N")
    suffix <- "post_chemo"
  }
  
  
  if (f == "IMvigor010"){
    if (n == 0){
      # FOr Imvigor010
      meta_data <- meta_data %>%
        dplyr::filter(grepl("Atezolizumab", ARM), prior_neoadjuvant_chemotherapy == "YES")
      suffix <- "Neo_Atezolizumab"
      n <- n+1
    } else{
      # FOr Imvigor010
      meta_data <- meta_data %>%
        dplyr::filter(grepl("Observation", ARM), prior_neoadjuvant_chemotherapy == "YES")
      suffix <- "Neo_Observation"
    }
  }
  
  meta_data <- meta_data %>%
    dplyr::mutate(Sample_ID = make.names(Sample_ID)) %>%
    dplyr::left_join(group_data, by=c("Sample_ID"="Sample_ID")) %>% 
    dplyr::mutate(model = Class) %>%
    dplyr::filter(model %in% c("Immune Enriched", "Immune Depleted"))
  
  low_cutoff <- quantile(meta_data$Mean, c(0.33))[[1]]
  high_cutoff <- quantile(meta_data$Mean, c(0.66))[[1]]
  
  meta_data <- meta_data %>% 
    dplyr::filter(Mean <= low_cutoff | Mean >= high_cutoff)
  
  prefix <- paste0(f, "_", suffix)
  survival_data <- meta_data
  summary <- plot_survival(survival_data, group, prefix)
  
  print(f)
  print(summary[[2]])
  print(summary[[3]])
  print(summary[[4]])
  print(summary[[5]])
  stats$gene          <- c(stats$gene, f)
  stats$HR            <- c(stats$HR, summary[[2]])
  stats$CI_lower      <- c(stats$CI_lower, summary[[3]])
  stats$CI_upper      <- c(stats$CI_upper, summary[[4]])
  stats$pvalue        <- c(stats$pvalue, summary[[5]])
  
  stats_df <- data.frame(stats)
}  
# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Summary")
openxlsx::writeData(wb, sheet = "Summary", x = stats_df)
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Individual_gene_stats.xlsx"), overwrite = TRUE)

#*****************Correlation NPEPPS vs each of 9 pathways *******************

# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r

# Preleminary test to check the test assumptions
# (i) Is the covariation linear? 
# Yes, form the plot above, the relationship is linear. 
# In the situation where the scatter plots show curved patterns, we are dealing 
# with nonlinear association between the two variables. IF non-linear, use
# Kendall tau  or Spearman correlation

# (ii) Are the data from each of the 2 variables (x, y) follow a normal distribution?
# Use Shapiro-Wilk normality test shapiro.test()
# and then look at the normality plot ggpubr::ggqqplot()

# Shapiro-Wilk test can be performed as follow:
# Null hypothesis: the data are normally distributed
# Alternative hypothesis: the data are not normally distributed

# Shapiro-Wilk normality test for mpg and wt
# shapiro.test(my_data$mpg) # => p = 0.1229
# shapiro.test(my_data$wt) # => p = 0.09
# The two p-values are greater than the significance level 0.05 implying that the
# distribution of the data are not significantly different from normal 
# distribution. In other words, we can assume the normality.

# ggqqplot(my_data$mpg, ylab = "MPG")
# ggqqplot(my_data$wt, ylab = "WT")
# From the normality plots, we conclude that both populations may come from normal
# distributions.
parent_path <- "/hpc/home/kailasamms/scratch/proj_Saswat/"
results_path <- "/hpc/home/kailasamms/scratch/proj_Saswat/"

meta_data <- read.xlsx(paste0(parent_path, "Meta_data_TCGA.xlsx"))
proj <- c(unique(meta_data$Project_ID), "IMvigor010", "IMvigor010", "IMvigor210")

df <- data.frame("Cancer" = "",
                 "Pathway" = "",
                 "Spearman_p_norm_all" = 0, 
                 "Spearman_r_norm_all" = 0, 
                 "Pearson_p_norm_all" = 0, 
                 "Pearson_r_norm_all" = 0,
                 "Spearman_p_z_all" = 0, 
                 "Spearman_r_z_all" = 0, 
                 "Pearson_p_z_all" = 0, 
                 "Pearson_r_z_all" = 0,
                 "Spearman_p_norm_33" = 0,
                 "Spearman_r_norm_33" = 0,
                 "Pearson_p_norm_33" = 0,
                 "Pearson_r_norm_33" = 0,
                 "Spearman_p_z_33" = 0,
                 "Spearman_r_z_33" = 0,
                 "Pearson_p_z_33" = 0,
                 "Pearson_r_z_33" = 0)
for (f in proj){
  
  data_1 <- read.xlsx(paste0(parent_path, "2 way classification mean 0/heatmaps 2 way mean cutoff 0/Heatmap_details_", f, ".xlsx"),
                      sheet="NPEPPS")
  data_1 <- data_1[, -c(1,7)]
  data_1 <- data_1 %>% dplyr::rename("Immune.sig.score" = "Mean")
  
  data_2 <- read.xlsx(paste0(parent_path, "2 way classification mean 0/heatmaps 2 way mean cutoff 0/Heatmap_details_", f, ".xlsx"),
                      sheet="Heatmap_matrix")
  colnames(data_2)[1] <- "Pathway"
  data_2 <- data_2 %>%
    tibble::column_to_rownames("Pathway") %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column("Sample")
  
  data <- data_1 %>% 
    dplyr::left_join(data_2, by=c("Sample"="Sample"))
  
  var_y <- setdiff(c(colnames(data_2), "Immune.sig.score"), "Sample")
  pathway <- c()
  spearman_p_norm_all <- c()
  spearman_r_norm_all <- c()
  pearson_p_norm_all <- c()
  pearson_r_norm_all <- c()
  spearman_p_z_all <- c()
  spearman_r_z_all <- c()
  pearson_p_z_all <- c()
  pearson_r_z_all <- c()
  spearman_p_norm_33 <- c()
  spearman_r_norm_33 <- c()
  pearson_p_norm_33 <- c()
  pearson_r_norm_33 <- c()
  spearman_p_z_33 <- c()
  spearman_r_z_33 <- c()
  pearson_p_z_33 <- c()
  pearson_r_z_33 <- c()
  
  for (i in 1:length(var_y)){
    
    pathway <- c(pathway, var_y[i])
    
    # Correlation using NPEPPS norm counts
    comp_data <- data %>% 
      dplyr::select(NPEPPS.x, all_of(var_y[i]))
    
    res <- cor.test(x=comp_data[[1]], y=comp_data[[2]], method = "spearman")
    spearman_p_norm_all <- c(spearman_p, res$p.value)
    spearman_r_norm_all <- c(spearman_r, res$estimate)
    
    res <- cor.test(x=comp_data[[1]], y=comp_data[[2]], method = "pearson")
    pearson_p_norm_all <- c(pearson_p, res$p.value)
    pearson_r_norm_all <- c(pearson_r, res$estimate)
    
    # Correlation using NPEPPS z score
    comp_data <- data %>% 
      dplyr::select(NPEPPS.y, all_of(var_y[i]))
    
    res <- cor.test(x=comp_data[[1]], y=comp_data[[2]], method = "spearman")
    spearman_p_z_all <- c(spearman_p, res$p.value)
    spearman_r_z_all <- c(spearman_r, res$estimate)
    
    res <- cor.test(x=comp_data[[1]], y=comp_data[[2]], method = "pearson")
    pearson_p_z_all <- c(pearson_p, res$p.value)
    pearson_r_z_all <- c(pearson_r, res$estimate)
    
    # Correlation using NPEPPS norm counts on top 33% and bottom 33%
    comp_data <- data %>% 
      dplyr::filter(!is.na(Category)) %>%
      dplyr::select(NPEPPS.x, all_of(var_y[i]))
    
    res <- cor.test(x=comp_data[[1]], y=comp_data[[2]], method = "spearman")
    spearman_p_norm_33 <- c(spearman_p_norm_33, res$p.value)
    spearman_r_norm_33 <- c( spearman_r_norm_33, res$estimate)
    
    res <- cor.test(x=comp_data[[1]], y=comp_data[[2]], method = "pearson")
    pearson_p_norm_33 <- c(pearson_p_norm_33, res$p.value)
    pearson_r_norm_33 <- c(pearson_r_norm_33, res$estimate)
    
    # Correlation using NPEPPS z score on top 33% and bottom 33%
    comp_data <- data %>% 
      dplyr::filter(!is.na(Category)) %>%
      dplyr::select(NPEPPS.y, all_of(var_y[i]))
    
    res <- cor.test(x=comp_data[[1]], y=comp_data[[2]], method = "spearman")
    spearman_p_z_33 <- c(spearman_p_z_33, res$p.value)
    spearman_r_z_33 <- c(spearman_r_z_33, res$estimate)
    
    res <- cor.test(x=comp_data[[1]], y=comp_data[[2]], method = "pearson")
    pearson_p_z_33 <- c(pearson_p_z_33 , res$p.value)
    pearson_r_z_33 <- c(pearson_r_z_33, res$estimate)
    
  }
  summary <- data.frame("Cancer" = f,
                        "Pathway" = pathway,
                        "Spearman_p_norm_all" = spearman_p_norm_all, 
                        "Spearman_r_norm_all" = spearman_r_norm_all, 
                        "Pearson_p_norm_all" = pearson_p_norm_all, 
                        "Pearson_r_norm_all" = pearson_r_norm_all,
                        "Spearman_p_z_all" = spearman_p_z_all, 
                        "Spearman_r_z_all" = spearman_r_z_all, 
                        "Pearson_p_z_all" = pearson_p_z_all, 
                        "Pearson_r_z_all" = pearson_r_z_all,
                        "Spearman_p_norm_33" = spearman_p_norm_33,
                        "Spearman_r_norm_33" = spearman_r_norm_33,
                        "Pearson_p_norm_33" = pearson_p_norm_33,
                        "Pearson_r_norm_33" = pearson_r_norm_33,
                        "Spearman_p_z_33" = spearman_p_z_33,
                        "Spearman_r_z_33" = spearman_r_z_33,
                        "Pearson_p_z_33" = pearson_p_z_33,
                        "Pearson_r_z_33" = pearson_r_z_33)
  
  df <- df %>% dplyr::bind_rows(summary)
}
df <- df [-1,]

# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Summary")
openxlsx::writeData(wb, sheet = "Summary", x = df)
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Correlation.xlsx"), overwrite = TRUE)


