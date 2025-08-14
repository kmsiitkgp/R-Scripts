
#********************3D UMAP of NPEPPS in Immune enriched and depleted******
#**********DIDNT WORK..NPEPPS is all over the place not just in immune depleted patients

# for (s in c("gsva", "ssgsea", "z")){
#   
#   # Get the UMPA coordinates
#   mat <- read.xlsx(paste0(data_path, "UMAP_Coords.xlsx"), 
#                    sheet = paste0("UMAP_", s))
#   
#   # Get the immune classification along with NPEPPS expression
#   group <- read.xlsx(paste0(data_path, "Two_Group_Classification.xlsx"), 
#                      sheet = paste0("Classification_", s))
#   
#   # Combine UMAP coords with NPEPPS expression and immune classification 
#   mat <- mat %>% 
#     dplyr::left_join(group, by=c("Sample"="Sample"))
#   
#   # Create a dataframe with values and color codes
#   # col_df <- data.frame(col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(1000),
#   #                      breaks = c(seq(from = min(mat$NPEPPS.x), to = 0, length.out = 500),
#   #                                 seq(from = max(mat$NPEPPS.x)/100, to = max(mat$NPEPPS.x), length.out = 500))) %>%
#   #   dplyr::distinct_at("col", .keep_all = TRUE)
#   
#   col_df <- data.frame(col = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(1000),
#                        breaks = c(seq(from = -1, to = 0, length.out = 500),
#                                   seq(from = 1/100, to = 1, length.out = 500))) %>%
#     dplyr::distinct_at("col", .keep_all = TRUE)
#   
#   # Find the closest color code for NPEPPS expression
#   for (i in 1:nrow(mat)){
#     mat$color[i] <- col_df$col[which.min(abs(col_df$breaks - mat$NPEPPS.x[i]))]
#   } 
#   
#   # Save the plot
#   jpeg(file=paste0(data_path, "UMAP_3D_NPEPPS", s, ".tiff"))
#   
#   # If you get "Error in check.plt(parplt) : figure margins too large,type 
#   # par("mar"). You will see values as [1] 5.1 4.1 4.1 2.1. Change like below
#   par(mar=c(1,1,1,1))
#   scatter3D(x=mat$UMAP1, y=mat$UMAP2, z=mat$UMAP3, 
#             theta = 40, phi = 15,  # viewing angles
#             bty = "b2",      # back panels and grid lines are visible
#             colvar = NULL,   # avoid coloring by z variable
#             col = mat$color, # variable used to color points
#             pch = 19,        # shape of points
#             cex = 0.5,       # size of points
#             colkey = list(mat$Class),
#             clab=c("Class"))
#   
#   dev.off()
# }


#*****************Correlation NPEPPS vs each of 22 pathways *******************

# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r

# Preliminary test to check the test assumptions
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
data_path <- "/hpc/home/kailasamms/scratch/proj_Saswat/"
results_path <- "/hpc/home/kailasamms/scratch/proj_Saswat/"

meta_data <- read.xlsx(paste0(data_path, "Meta_data_TCGA.xlsx"))
proj <- c(unique(meta_data$Project_ID), "IMvigor010", "IMvigor210")

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

  data_1 <- read.xlsx(paste0(data_path, "2 way classification mean 0/heatmaps 2 way mean cutoff 0/Heatmap_details_", f, ".xlsx"),
                      sheet="NPEPPS")
  data_1 <- data_1[, -c(1,7)]
  data_1 <- data_1 %>% dplyr::rename("Immune.sig.score" = "Mean")

  data_2 <- read.xlsx(paste0(data_path, "2 way classification mean 0/heatmaps 2 way mean cutoff 0/Heatmap_details_", f, ".xlsx"),
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
openxlsx::saveWorkbook(wb, file = paste0(data_path, "Correlation.xlsx"), overwrite = TRUE)

#*****************Plot KM curves *******************

#source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")

#data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/"
data_path <- "/hpc/home/kailasamms/scratch/proj_Saswat/"
results_path <- "/hpc/home/kailasamms/"

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
color_palette <- c("#CB181D","#A6D854")
variable_x <- "Time"
variable_y <- "Status" 
gene <- ""
group <- ""

# meta_data <- read.xlsx(paste0(data_path, "Meta_data_TCGA.xlsx"))
# proj <- c(unique(meta_data$Project_ID), "IMvigor010", "IMvigor010", "IMvigor210")

proj <- c("IMvigor010", "IMvigor210")

# Create a list to store survminer cutoffs, coxph stats, etc..
stats <- list("gene" = c(),
              "HR" = c(),
              "CI_lower" = c(),
              "CI_upper" = c(),
              "pvalue" = c())

n <- 0
for (f in proj){
  
  meta_data <- read.xlsx(paste0(data_path, "Meta_data_", f, ".xlsx"))
  group_data <- read.xlsx(paste0(data_path, "Heatmap_details_", f, ".xlsx"),
                          sheet="NPEPPS") %>%
    dplyr::select(Sample, gsva_class, gsva.scores)
  
  suffix <- ""
  
  # For Imvigor210 and Imvigor210_new
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
    dplyr::left_join(group_data, by=c("Sample_ID"="Sample")) %>% 
    dplyr::mutate(model = gsva_class) %>%
    dplyr::filter(model %in% c("Immune Enriched", "Immune Depleted"))
  
  # low_cutoff <- quantile(meta_data$Mean, c(0.33))[[1]]
  # high_cutoff <- quantile(meta_data$Mean, c(0.66))[[1]]
  
  # meta_data <- meta_data %>% 
  #   dplyr::filter(Mean <= low_cutoff | Mean >= high_cutoff)
  
  prefix <- paste0(f, "_", suffix)
  survival_data <- meta_data
  summary <- plot_survival(survival_data, group, prefix, results_path)
  
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
openxlsx::saveWorkbook(wb, file = paste0(data_path, "KMplotter_Individual_gene_stats.xlsx"), overwrite = TRUE)

##### Hits in Imvigor210

hits <- c("NPEPPS", "MAOA", "IGF1R", "GAN", "KLF5", "FGFR3", "GCLC", "PTK2", 
          "VEGFA", "KDM1A", "PLXNB1", "SRPK1", "CXADR", "IL20RA", "SLK", 
          "ERBB2", "ULK1", "FASN")

meta_data <- read.xlsx(paste0(data_path, "Meta_data_", f, ".xlsx"))
meta_data <- meta_data %>% 
  dplyr::filter(Tissue == "bladder", Received.platinum == "Y", `Sample.collected.pre-platinum` == "N")
suffix <- "post_chemo"

prefix <- paste0(f, "_", suffix)
survival_data <- meta_data
wrangle_data(expr_df, stratify_criteria, prefix, results_path)

#****************Correlate tumor purity score to all genes
# unfortunately, NPEPPS doesnt correlate with tumor purity

# source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")
# #source("/hpc/home/kailasamms/projects/scRNASeq/scRNASeq_Seurat_Functions_Variables.R")
# data_path <- "/hpc/home/kailasamms/scratch/proj_Saswat/"
# library(GSVA)
# 
# # meta_data <- read.xlsx(paste0(data_path, "Meta_data_TCGA.xlsx"))
# # proj <- unique(meta_data$Project_ID)
# # proj <- c(proj, "IMvigor010", "IMvigor210", "PanCancer")
# 
# proj <- "PanCancer"
# purity_markers <- c("CSF2RB", "RHOH", "C1S", "CCDC69", "CCL22", "CYTIP", 
#                     "POU2AF1", "FGR", "CCL21", "IL7R")
# 
# # Import expression data that has already been normalized across samples
# for (f in proj){
#   
#   if (f == "PanCancer"){
#     read_data <- read.table(paste0(data_path, "Normalized_Counts_TCGA.tsv"), header=TRUE)
#     normalized_counts <- read.xlsx(paste0(data_path, "TCGA_median_centered_log_normalized_counts.xlsx"))
#   } else if(f %in% c("IMvigor010", "IMvigor210")){
#     read_data <- read.xlsx(paste0(data_path, "Normalized_Counts_", f, ".xlsx"))
#     read_data <- read_data[,-c(2:4)]
#   } else{
#     read_data <- read.xlsx(paste0(data_path, "Normalized_Counts_DESeq2_modelled", f, ".xlsx"))
#     read_data <- read_data[,-c(2:5)]
#   }
#   
#   read_data <- read_data %>% 
#     base::replace(is.na(.), 0)                           # replace NA with 0
#   read_data <- read_data[rowSums(read_data[,-1]) !=0,]   # remove genes with no expression in all samples
#   colnames(read_data)[1] <- "SYMBOL"
#   
#   if (f != "PanCancer"){
#     # log normalize and median center the counts
#     normalized_counts <- read_data %>%
#       dplyr::group_by(SYMBOL) %>%
#       dplyr::summarise(across(.cols=everything(), .fns=mean)) %>%
#       tibble::remove_rownames() %>%
#       tibble::column_to_rownames("SYMBOL")
#     normalized_counts <- log(1+normalized_counts, base=2)
#     t <- base::apply(X=normalized_counts, MARGIN=1, FUN=median, na.rm=TRUE)
#     normalized_counts <- base::sweep(x=normalized_counts, MARGIN=1, FUN="-", STATS=t)
#   }
#   
#   #***** Calculate purity scores for each patient *****#
#   
#   # # Calculate gsva and ssgsea scores (TAKING TOO MUCH TIME ON PANCANCER..SO SKIP)
#   # gs <- list()
#   # plot_genes <- list(purity_markers)
#   # names(plot_genes) <- "purity_markers"
#   # gs <- c(gs, plot_genes)
#   # 
#   # gsvaPar <- GSVA::gsvaParam(exprData = as.matrix(normalized_counts), 
#   #                            geneSets = as.list(gs))
#   # gsva.scores <- gsva(gsvaPar, 
#   #                     verbose=TRUE)
#   # 
#   # ssgseaPar <- GSVA::ssgseaParam(exprData = as.matrix(normalized_counts), 
#   #                                geneSets = as.list(gs))
#   # ssgsea.scores <- GSVA::gsva(ssgseaPar, 
#   #                             verbose=TRUE)
#   
#   # Calculate scores based on Levine et al 
#   # Create empty dataframe of samples as rownames and pathway as column names
#   z.scores <- normalized_counts %>% 
#     t() %>%
#     data.frame() %>%
#     dplyr::select(colnames(.)[1]) %>%  # Keep a dummy column to keep dataframe structure intact
#     tibble::rownames_to_column("Sample")
#   
#   plot_genes <- purity_markers
#   expr_df <- as.data.frame(advanced_Z(plot_genes, normalized_counts))
#   colnames(expr_df) <- "purity_markers"
#   expr_df <- expr_df %>%
#     tibble::rownames_to_column("Sample")
#   
#   z.scores <- dplyr::left_join(z.scores, expr_df, by=c("Sample"="Sample"))
#   
#   z.scores <- z.scores[,-2]  # remove dummy column
#   # z.scores <- z.scores %>% 
#   #   tibble::column_to_rownames("Sample") %>%
#   #   t()
#   
#   # rownames(z.scores) <- make.names(rownames(z.scores))
#   # colnames(z.scores) <- make.names(colnames(z.scores))
#   
#   # summary_df <- data.frame(Sample = rownames(t(z.scores)))
#   # for (s in c("gsva.scores", "ssgsea.scores", "z.scores")){
#   #   
#   #   scores <- get(s) %>%
#   #     t() %>%
#   #     data.frame() %>%
#   #     dplyr::rename(!!rlang::sym(gsub(".scores",".purity.scores",rlang::sym(s))) := purity_markers) %>%
#   #     tibble::rownames_to_column("Sample")
#   #   
#   #   summary_df <- summary_df %>% 
#   #     dplyr::left_join(scores, by=c("Sample"="Sample"))
#   # }
#   summary_df <- z.scores %>% 
#     dplyr::rename(z.purity.scores = purity_markers)
# }
# 
# # Identify genes upregulated in immune depleted patients across all cancers
# # NPEPPS was up in 13 cancers but down in 1 cancer. So, we set n_DOWN<2 instead 
# # of n_DOWN=0
# degs <- read.xlsx(paste0(data_path, "Summary.xlsx"),
#                   sheet="Summary_gsva")
# colnames(degs)[1] <- "Gene"
# degs <- degs %>%
#   dplyr::filter(n_DOWN < 2)
# degs <- degs$Gene
# degs <- make.names(degs)
# 
# # Get median centered values of DEGs
# degs_zscore <- normalized_counts %>% 
#   tibble::rownames_to_column("SYMBOL") %>%
#   dplyr::mutate(SYMBOL = make.names(SYMBOL)) %>%
#   dplyr::filter(SYMBOL %in% degs) %>% 
#   tibble::column_to_rownames("SYMBOL") %>%
#   t() %>%
#   data.frame() %>%
#   tibble::rownames_to_column("Sample") %>%
#   dplyr::mutate(Sample = make.names(Sample))
# 
# # Get normalized counts of DEGs
# degs_norm <- read_data %>% 
#   dplyr::mutate(SYMBOL = make.names(SYMBOL)) %>%
#   dplyr::filter(SYMBOL %in% degs) %>% 
#   tibble::column_to_rownames("SYMBOL") %>%
#   t() %>%
#   data.frame() %>%
#   tibble::rownames_to_column("Sample") %>%
#   dplyr::mutate(Sample = make.names(Sample))
# 
# # Combine purity scores with DEGs
# degs_zscore <- degs_zscore %>%
#   dplyr::right_join(summary_df, by=c("Sample"="Sample"))
# 
# degs_norm <- degs_norm %>%
#   dplyr::right_join(summary_df, by=c("Sample"="Sample"))
# 
# # Calculate correlation fo each of the hits with NPEPPS
# df <- data.frame("Gene" = "",
#                  "Spearman_p_z" = 0, 
#                  "Spearman_r_z" = 0, 
#                  "Pearson_p_z" = 0, 
#                  "Pearson_r_z" = 0,
#                  "Spearman_p_norm" = 0, 
#                  "Spearman_r_norm" = 0, 
#                  "Pearson_p_norm" = 0, 
#                  "Pearson_r_norm" = 0)
# gene <- c()
# spearman_p_z <-c() 
# spearman_r_z <-c() 
# pearson_p_z <-c() 
# pearson_r_z <-c()
# spearman_p_norm <-c() 
# spearman_r_norm <-c() 
# pearson_p_norm <-c() 
# pearson_r_norm <-c()
# 
# # Correlation using NPEPPS z score counts
# for (i in (2:(ncol(degs_zscore)-1))){   
#   
#   gene <- c(gene, colnames(degs_zscore)[i])
#   
#   comp_data1 <- degs_zscore %>%
#     dplyr::select(z.purity.scores, colnames(degs_zscore)[i])
#   
#   res <- cor.test(x=comp_data1[[1]], y=comp_data1[[2]], method = "spearman", exact=FALSE)
#   spearman_p_z <- c(spearman_p_z, res$p.value)
#   spearman_r_z <- c(spearman_r_z, res$estimate)
#   
#   res <- cor.test(x=comp_data1[[1]], y=comp_data1[[2]], method = "pearson", exact=FALSE)
#   pearson_p_z <- c(pearson_p_z, res$p.value)
#   pearson_r_z <- c(pearson_r_z, res$estimate)
#   
#   comp_data2 <- degs_norm %>%
#     dplyr::select(z.purity.scores, colnames(degs_zscore)[i])
#   
#   res <- cor.test(x=comp_data2[[1]], y=comp_data2[[2]], method = "spearman", exact=FALSE)
#   spearman_p_norm <- c(spearman_p_norm , res$p.value)
#   spearman_r_norm <- c(spearman_r_norm, res$estimate)
#   
#   res <- cor.test(x=comp_data2[[1]], y=comp_data2[[2]], method = "pearson", exact=FALSE)
#   pearson_p_norm <- c(pearson_p_norm, res$p.value)
#   pearson_r_norm <- c(pearson_r_norm, res$estimate)
#   
#   print(i)
# }
# 
# cor_results <- data.frame("Gene" = gene,
#                           "Spearman_p_z" = spearman_p_z, 
#                           "Spearman_r_z" = spearman_r_z, 
#                           "Pearson_p_z" = pearson_p_z, 
#                           "Pearson_r_z" = pearson_r_z,
#                           "Spearman_p_norm" = spearman_p_norm, 
#                           "Spearman_r_norm" = spearman_r_norm, 
#                           "Pearson_p_norm" = pearson_p_norm, 
#                           "Pearson_r_norm" = pearson_r_norm)
# 
# df <- df %>% dplyr::bind_rows(cor_results)
# df <- df [-1,]
# 
# df <- df %>% 
#   dplyr::mutate(Max = pmax(Spearman_r_z, Pearson_r_z, Spearman_r_norm, Pearson_r_norm),
#                 Hits = dplyr::case_when(Max > 0 ~ "YES", TRUE ~ "NO"))
# 
# # Save the results
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb, sheetName = "Summary")
# openxlsx::writeData(wb, sheet = "Summary", x = df)
# openxlsx::saveWorkbook(wb, file = paste0(data_path, "Correlation_with_purity.xlsx"), overwrite = TRUE)

###################

# library(factoextra)
# library(igraph)
# library(mlbench)    # for Ionosphere data
# library(psych)      # for cor2dist
# data(Ionosphere)
# 
# scores1 <- scores[rownames(scores) %in% c(pro_tumor, anti_tumor),]
# t1 <- kmeans(t(scores1), centers=2, nstart=50)
# fviz_cluster(t1, t(scores1), labelsize=0)

# col_annotation <- scores1 %>%
#   data.frame() %>%
#   tibble::rownames_to_column("Pathway") %>%
#   dplyr::mutate(Group = dplyr::case_when(Pathway %in% angio ~ "angio",
#                                          Pathway %in% pro_tumor ~ "pro_tumor",
#                                          Pathway %in% anti_tumor ~ "anti_tumor",
#                                          TRUE ~ "malignant_cell")) %>%
#   tibble::column_to_rownames("Pathway") %>%
#   dplyr::group_by(Group) %>%
#   dplyr::summarise(across(.cols=everything(), .fns=mean)) %>%
#   dplyr::ungroup() %>%
#   tibble::column_to_rownames("Group") %>%
#   t() %>%
#   data.frame() %>%
#   dplyr::mutate(Class = dplyr::case_when(angio > 0 & pro_tumor > 0 & anti_tumor > 0 ~ "Immune Enriched Fibrotic",
#                                          angio < 0 & pro_tumor > 0 & anti_tumor > 0 ~ "Immune Enriched Non-Fibrotic",
#                                          angio > 0 & pro_tumor < 0 & anti_tumor < 0 ~ "Fibrotic",
#                                          angio < 0 & pro_tumor < 0 & anti_tumor < 0 ~ "Depleted",
#                                          TRUE ~ "Unclassified"),
#                 Dummy_col = NA) %>%  # ordering by rownames doesnt work for 1 column df
#   dplyr::select(Class, Dummy_col)
