#!/usr/bin/env Rscript

#******************************************************************************#
#                   INFORMATION REGARDING SURVIVAL PARAMETERS                  #
#******************************************************************************#

## NOTE: Some predefined columns name MUST be present in the meta/clinical data
## "Sample_ID" : this column MUST contain sample names
## "Time:      : this column MUST contain survival duration in months
## "Status     : this column MUST contain Dead/Alive status (Alive = 0, Dead = 1)
## "Sex        : this column MUST contain "Male" or "Female"

# "plot_by" defines parameter based on which survival curves are plotted.
# By default, survival curves are plotted on gene expression of single/multiple
# genes i.e. plot_by = c(NA).  
# If "plot_by = c("Race")", a column named "Race" MUST exist in metadata.

# "split_by" defines parameter based on which patients are split into groups.
# To plot curves of gene X for male and female patients, a column named "Sex" 
# MUST exist in metadata.

# "split_plot" defines whether separate plots need to be made for patients
# grouped using "split_by" variable. 
# To plot curves of gene X for male and female patients SEPARATELY, set 
# "split_plot = TRUE"
# If "split_by == c(NA)", then split_plot variable is ignored

# "multiple_cutoff" defines whether to calculate a single cutoff for all 
# patients or multiple cutoffs for each group of patients grouped using 
# "split_by" variable

# "stratify_criteria" defines method used to determine cutoff for stratifying samples
# "m" : stratify using median cutoff (top 50% vs bottom 50%)
# "t" : stratify using tertile expression (top 33% vs bottom 33%)
# "q" : stratify using quartile expression (top 25% vs bottom 25%)
# "o" : stratify using optimal cutoff calculated by survminer
# "th": stratify using thirds cutoff

# "reference" defines the reference level for calculating Hazard Ratio (HR)
# By default, when plotting by expression, reference is set to "LOW".
# If plot_by = c("Race") etc, set appropriate reference

# "conf_interval" defines whether confidence interval are shown
# "plot_curve" defines whether to plot the survival curve
# "plot_risk_table" defines whether to plot the risk table

# "legend_title" indicate the title of legend

# "plot_all_bins" defines whether to plot all bins or ONLY 2 bins (high, low)
# stratify_criteria = "m" splits samples into 2 bins i.e. below 50%, above 50%
# stratify_criteria = "t" splits samples into 3 bins i.e. below 33%, 33%-67%, above 67%
# stratify_criteria = "q" splits samples into 4 bins i.e. below 25%, 25%-50%, 50%-75%, above 75%
# stratify_criteria = "o" splits samples into 2 bins i.e. above & below optimum cutoff
# stratify_criteria = "th" splits samples into 3 bins i.e. bottom 33%, middle33%, top 33% based on expression range

# "plot_all_quartiles" is relevant ONLY when stratify_criteria = "q" and defines
# whether to plot all 4 bins or ONLY 2 bins (high, low)

# "gene_sig_score" defines whether to calculate gene signature score

# Calculate z-score using the function described in Levine et al https://doi.org/10.1186/gb-2006-7-10-r93
# NOTE: z-score is calculated within each sample across all genes. So, 
# number of samples doesnt affect a sample's z-score.
# prep_expr_df() does this if survival_params$gene_sig_score = TRUE

# VERY VERY IMPORTANT: ggsurvplot() labels the groups in alphabetical order. So,
# when we want to use custom labels, initialize them in alphabetical order.
# Eg: c("High", "Low") instead of  c("Low, "High")
color_palette <- c("#F6D2E0", "#C8E7F5")      # 1st color~female, 2nd color~male
color_palette <- c("orchid2", "dodgerblue2")  # 1st color~female, 2nd color~male
color_palette <- c("#EE7AE9", "#1C86EE")      # 1st color~female, 2nd color~male
color_palette <- c("#F6D2E0", "#C8E7F5","#EE7AE9", "#1C86EE" )
legend_label <- c("High", "Low")              # Reference is Low
color_palette <- c("#DB6D00", "#490092")      # 1st color~high, 2nd color~low
color_palette <- c("#d73027","#0c2c84")
color_palette <- c(color_palette, 
                   colorspace::adjust_transparency(col = color_palette, alpha = 0.3), 
                   colorspace::adjust_transparency(col = color_palette, alpha = 0.6))

#******************************************************************************#
#                          DEFINE SURVIVAL PARAMETERS                          #
#******************************************************************************#

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")
source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

data_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"
output_path <- "/hpc/home/kailasamms/scratch/ICB_TCGA_datasets/processed/"

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

#******************************************************************************#
#                                 IMPORT DATA                                  #
#******************************************************************************#

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
    dplyr::select(all_of(intersect(make.names(meta_data$Sample_ID), colnames(norm_counts))))
  
  colnames(norm_counts) <- base::make.names(names = colnames(norm_counts))
  
  # Perform log1p transformation on the subsetted data if not done before
  # Perform median centering for each gene across samples as median expression 
  # is more robust to outliers as compared to mean expression.
  
  # NOTE: If plotting a single gene, using normalized counts or log1p transformed 
  # normalized counts or median centered log1p transformed normalized counts 
  # will yield same curves, pvalues, HR values, etc. 
  # Results will be very different in the 3 cases ONLY for gene signatures
  log_norm_counts <- log(1+norm_counts, base=2)
  t <- base::apply(X=log_norm_counts, MARGIN=1, FUN=median, na.rm=TRUE)
  log_norm_counts <- base::sweep(x=log_norm_counts, MARGIN=1, FUN="-", STATS=t)
  
  # List genes to plot
  plot_genes <- c("NPEPPS")
  plot_genes <- intersect(plot_genes, rownames(log_norm_counts))
  
  #******************************************************************************#
  #                           PERFORM SUVIVAL ANALYSIS                           #
  #******************************************************************************#
  
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
    openxlsx::saveWorkbook(wb, file = paste0(data_path, proj, "_Individual_stats.xlsx"), overwrite = TRUE)
  }
  
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
    openxlsx::saveWorkbook(wb, file = paste0(data_path, proj, "_Stats.xlsx"), overwrite = TRUE)
  }
}

#******************************************************************************#
#                      OPTIONAL: CHECK IF GENE SIGNATURE IS GOOD               #
#******************************************************************************#

# Create an AnnotatedDataFrame object for meta data
p_data <- meta_data
rownames(p_data) <- make.names(p_data$Sample_ID)
p_data <- Biobase::AnnotatedDataFrame(data = p_data,
                                      varMetadata = data.frame("labelDescription" = colnames(p_data))) 

# Create an AnnotatedDataFrame object for features
f_data <- data.frame("SYMBOL" = rownames(normalized_counts))
rownames(f_data) <- make.names(f_data$SYMBOL)
f_data <- Biobase::AnnotatedDataFrame(data = f_data,
                                      varMetadata = data.frame("labelDescription" = colnames(f_data))) 

# Create an ExpressionSet object for read data
e_data <- as.matrix(normalized_counts)
eset <- Biobase::ExpressionSet(assayData = e_data,
                               phenoData = p_data,
                               featureData = f_data,
                               annotation = "custom")

# NOTE: classes parameter in sigCheck() = Status column of varLabels(eset)
# NOTE: survival parameter in sigCheck() = Time column of varLabels(eset)
varLabels(eset)
eset$Status
eset$Time
# NOTE: annotation parameter in sigCheck() = SYMBOL column of fvarLabels(eset)
# NOTE: plot_genes and genes in SYMBOL column must have some overlap
fvarLabels(eset)   

# Perform anlysis for male and female
# sigCheck is not good at classifying samples into optimal groups. So, we 
# manually classify the samples using survminer and import the classification
# into sigCheck() using scoreMethod and threshold parameters

p <- c()

for (i in 1:100){
  
  plot_genes <- base::sample(x=rownames(normalized_counts), size=64, replace = FALSE)
  
  # Calculate z-score using the function described in above paper
  expr_df <- as.data.frame(advanced_Z(plot_genes, normalized_counts))
  
  # Merge expression data with survival data
  expr_df <- expr_df %>%
    data.frame() %>%
    dplyr::rename(combined.exp = identity(1)) %>%
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::inner_join(meta_data, by=c("Sample_ID"="Sample_ID"))
  
  gene <- "combined.exp"
  if (nrow(expr_df > 0)){
    summary <- wrangle_data(expr_df, stratify_criteria, prefix)
  }
  
  surv_df <- summary[[1]] %>%
    dplyr::filter(Status == 0 | Status == 1) %>%
    dplyr::mutate(Expression = stringi::stri_replace_all_regex(str = Expression,
                                                               pattern = c("HIGH", "LOW"),
                                                               replacement = c("1", "0"),
                                                               vectorise_all = FALSE)) %>%
    dplyr::filter(Sex == "Male")
  
  sigCheck_score <- function(eset){
    e <- surv_df  %>% 
      dplyr::select(Expression) %>%
      unlist(., use.names = FALSE) %>%
      as.numeric()
    
    return(e)
  }
  
  # Format the object to remove sample not present in surv_df
  # Also, remove samples which have status other than 0 or 1.
  eset_subset <- eset[, eset$Sample_ID %in% surv_df$Sample_ID]
  
  # Create a SigCheck object
  check <- sigCheck(expressionSet = eset_subset, 
                    classes = "Status", 
                    survival = "Time",
                    signature = plot_genes,
                    annotation = "SYMBOL",
                    scoreMethod = sigCheck_score(eset),
                    threshold = sum(sigCheck_score(eset))/length(sigCheck_score(eset)))
  
  p <- c(p, check@survivalPval)
  # sigCheckRandom(check = check,
  #                iterations=100)
}

#******************************************************************************#




























#******************************************************************************#
#                     STEP 2: PLOT SURVIVAL CURVES BY SEX                      #
#******************************************************************************#

# Read survival data
expr_df <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "BBN original figs/1A_data.xlsx"))

gse <- "Mice Expt"
sex <- "Male vs Female"
i <-  ""

# Plot
summary <- plot_survival(survival_data, group, prefix)

#******************************************************************************#
#                    STEP 3: PLOT SURVIVAL CURVES FOR GENES                    #
#******************************************************************************#

# If workign with huge datasets, try this
options(future.globals.maxSize = 30000 * 1024^2)

# Define the gse of the project
gse <- "TCGA_BLCA"         # DESeq2 normalized counts
gse <- "Blaveri"           # already log2 transformed & median centered
gse <- "mskcc"             # already log2 transformed
gse <- "gse13507"          # already log2 transformed
gse <- "gse31684"          # already log2 transformed
gse <- "gse32894"          # already log2 transformed
gse <- "Imvigor210"        # DESeq2 normalized counts
gse <- "Imvigor010"        # DESeq2 normalized counts
gse <- "Imvigor210_old" 


if (gse == "TCGA_BLCA" | gse == "Imvigor210" | gse == "Imvigor010"){
  normalized_counts <- log(1+normalized_counts, base=2)
}
if (gse != "Blaveri"){
  t <- base::apply(X=normalized_counts, MARGIN=1, FUN=median, na.rm=TRUE)
  normalized_counts <- base::sweep(x=normalized_counts, MARGIN=1, FUN="-", STATS=t)
}

# Human Y gene signature
plot_genes <- c("DDX3Y", "EIF1AY", "HSFY2", "KDM5D", "UTY", "NLGN4Y", 
                "PCDH11Y", "RPS4Y1", "TBL1Y", "TMSB4Y", "USP9Y", "ZFY", 
                "DAZ1", "DAZ2", "DAZ3", "DAZ4", "PRY2", "RBMY1A1")


