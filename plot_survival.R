#!/usr/bin/env Rscript

# BiocManager::install("SigCheck")

# Read this paper
# https://doi.org/10.1093/jncimonographs/lgu024

#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

# Data wrangling packages     
library("openxlsx")
library("dplyr")
library("tibble")

# Graph plotting packages
library("ggplot2")
library("grid")  

# Specialized Graph plotting packages
library("survival")
library("survminer")  #If this is not loaded %++% wont work
library("SigCheck")

#******************************************************************************#
#                           DECLARE GLOBAL VARIABLES                           #
#******************************************************************************#

## NOTE: Some predefined columns name MUST be present in the meta/clinical data
## "Sample_ID" : this column MUST contain sample names
## "Time:      : this column MUST contain survival duration in months
## "Status     : this column MUST contain Dead/Alive status (Alive = 0, Dead = 1)
## "Sex        : this column MUST contain sex/Sex "Male" or "Female"

## Declare column in metadata to subset samples and its values
## If plotting the whole data, set subset_group <- NA and subset_value <- NA
## If metadata has column named "Race" with "Asian", "African", "American" 
## values and you want to plot only Asian and African origin patients, define
## subset_group <- "Race" and subset_value <- c("Asian", "African").
subset_group <- "Sex"   #"Received.platinum" # "Sample.collected.pre.platinum" 
subset_value <- c("Male")  #c("Y", "N")  #c("Y", "N", "NE")
subset_group <- "Stage"   #"Received.platinum" # "Sample.collected.pre.platinum" 
subset_value <- c("T1", "T2", "T3", "T4")  #c("Y", "N")  #c("Y", "N", "NE")
subset_group <- "binaryResponse"
subset_value <- c("SD/PD","CR/PR")
subset_group <- "Condition"
subset_value <- c("PDL1", "Control")
subset_group <- "prior_neoadjuvant_chemotherapy"
subset_value <- c("YES", "NO")
subset_group <- "Sample.collected.pre.platinum"
subset_value <- c("Y", "N")
subset_group <- NA
subset_value <- NA
subset_group <- "Received.platinum"
subset_value <- c("Y", "N")

## "plot_by": define parameter to plot survival.
## To plot by gene expression, set as "Expression" (default). 
## To plot by "Race" etc, a column named "Race" MUST exist in metadata.
plot_by <- "Stage"
plot_by <- "Expression"

## "split_by": define parameter to split data. 
## If you want to plot curves of gene X for male and female patients, then a 
# column named "Sex" MUST exist in metadata.
split_by <- subset_group  #NA

## If you want to plot curves of gene X for male and female patients SEPARATELY,
## then set combine_plot <- FALSE
## If "split_by" <- NA, then combine_plot variable is ignored and its 
## value doesnt matter
combine_plot <- FALSE 

## If plotting by expression, expression cutoffs will be calculated to stratify
## patients into HIGH vs LOW groups. Define if you want to calculate a single 
## cutoff for all patients in subset_groups or individual cutoffs for each
## subset_group
multiple_cutoff <- TRUE

## Choose a stratify_criteria to set cutoff for stratifying samples
## "m" : stratify using median cutoff (top 50% vs bottom 50%)
## "t" : stratify using tertile expression (top 33% vs bottom 33%)
## "q" : stratify using quartile expression (top 25% vs bottom 25%)
## "o" : stratify using optimal cutoff calculated by survminer
## "th": stratify using thirds cutoff
stratify_criteria <- "o"

# Define reference level for calculating Hazard Ratio (HR)
# If you are plotting ONLY by expression, reference is set to "LOW" (default)
# If you are plotting by expression independent favtor like Race or Sex, set
# appropriate reference
reference <- "T1" #"NA"
reference <- dplyr::if_else(plot_by == "Expression", "LOW", reference)

## Indicate if you want to see confidence intervals in the plot
confidence_interval <- FALSE

## Indicate the title of legend
legend_title <- "Expression"  #Sex #Race

## indicate if you want to plot the risk table
plot_risk_table <- TRUE

## indicate if you want to plot the survival curve
plot_curve <- TRUE

## Indicate if you want to plot all quartiles or only HIGH vs LOW
all_quartiles <- FALSE

## Indicate if you are plotting a gene signature or single gene
gene_signature <- TRUE

# # VERY VERY IMPORTANT: ggsurvplot() labels the groups in alphabetical order. So, 
# # when we want to use custom labels, initialize them in alphabetical order. 
# # Eg: c("High", "Low") instead of  c("Low, "High")
legend_label <- c("Female", "Male")           # Reference is female
color_palette <- c("#F6D2E0", "#C8E7F5")      # 1st color~female, 2nd color~male
color_palette <- c("orchid2", "dodgerblue2")  # 1st color~female, 2nd color~male
color_palette <- c("#EE7AE9", "#1C86EE")      # 1st color~female, 2nd color~male
legend_label <- c("T1", "T2", "T3", "T4")
color_palette <- c("#F6D2E0", "#C8E7F5","#EE7AE9", "#1C86EE" )
legend_label <- c("High", "Low")              # Reference is Low
color_palette <- c("#DB6D00", "#490092")      # 1st color~high, 2nd color~low
color_palette <- c("#d73027","#0c2c84")
color_palette <- c(color_palette, 
                colorspace::adjust_transparency(col = color_palette, alpha = 0.3), 
                colorspace::adjust_transparency(col = color_palette, alpha = 0.6))

## Declare global variables for survival curves. 
variable_x <- "Time"
variable_y <- "Status" 

# Unique prefix to each plot
prefix <- ""

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


wb <- openxlsx::createWorkbook()
for (gse in c("TCGA_BLCA", "Blaveri")){
  
  ## Store path of parent directory i.e. root directory for the project
  parent_path <- parent_path <- "/hpc/home/kailasamms/"
  parent_path <- "C:/Users/KailasammS/Box/Saravana@cedars/10. Ongoing Projects/BBN project/BLCA_Cohorts_correct/"
  
  # Import read_data
  read_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, gse, "_Normalized.xlsx"))
  colnames(read_data)[1] <- "SYMBOL"
  
  # Import meta_data and subset meta_data if needed
  meta_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, gse, "_Metadata.xlsx"))  
  #%>% dplyr::filter(grepl(pattern = "T2|T3|T4", x = Stage))

  # Reformat metadata 
  meta_data <- meta_data %>% 
    dplyr::mutate(Sample_ID = make.names(names = Sample_ID)) %>%
    dplyr::mutate(Time = as.numeric(Time)) %>%
    dplyr::filter(Time > 0 & !is.na(Time)) %>%
    dplyr::distinct_at("Sample_ID", .keep_all = TRUE)
  
  if (!is.na(subset_group)){
    meta_data <- meta_data %>% 
      dplyr::filter(get(subset_group) %in% subset_value)
  }
  
  # meta_data <- meta_data %>%
  #   dplyr::filter(Sex == "Male", Received.platinum == "Y",
  #                 !is.na(`Sample.collected.pre.platinum`),
  #                 Tissue %in% c("bladder", "kidney", "ureter"))
  # 
  # meta_data <- meta_data %>%
  #   dplyr::filter(Sex == "Male", Received.platinum == "N",
  #                 Tissue %in% c("bladder", "kidney", "ureter"))
  # 
  # meta_data <- meta_data %>%
  #      dplyr::filter(Sex == "Male", Tissue == "bladder")
  # 
  # meta_data <- meta_data %>%
  #      dplyr::filter(Sex == "Male", Received.platinum == "Y",
  #                    !is.na(`Sample.collected.pre.platinum`))
  # 
  # meta_data <- meta_data %>%
  #      dplyr::filter(Sex == "Male", Received.platinum == "N")
  
  meta_data <- meta_data %>%
       dplyr::filter(Sex == "Male")
  
  # Reformat read data
  normalized_counts <- read_data %>%
    dplyr::mutate(SYMBOL = make.names(names = SYMBOL, unique = TRUE)) %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    tibble::column_to_rownames(var = "SYMBOL") %>%
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
  colnames(normalized_counts) <- base::make.names(names = colnames(normalized_counts))
  normalized_counts <- normalized_counts[,intersect(make.names(meta_data$Sample_ID), colnames(normalized_counts))]
  normalized_counts <- normalized_counts[!rowSums(normalized_counts, na.rm=TRUE) == 0,]
  
  # Perform log1p transformation on the subsetted data if not done before
  # Perform median centering for each gene across samples as median expression 
  # is more robust to outliers as compared to mean expression.
  
  # NOTE: If plotting a single gene, using normalized counts or log1p transformed 
  # normalized counts or median centered log1p transformed normalized counts 
  # will yield same curves, pvalues, HR values, etc. Results will be very 
  # different in the 3 cases ONLY for gene signatures
  
  if (gse == "TCGA_BLCA" | gse == "Imvigor210" | gse == "Imvigor010"){
    normalized_counts <- log(1+normalized_counts, base=2)
  }
  if (gse != "Blaveri"){
    t <- base::apply(X=normalized_counts, MARGIN=1, FUN=median, na.rm=TRUE)
    normalized_counts <- base::sweep(x=normalized_counts, MARGIN=1, FUN="-", STATS=t)
  }
  
  #*****************USE THIS SECTION FOR PLOTTING INDIVIDUAL GENES***************#
  
  if (gene_signature == FALSE){
    
    # Create a list of genes for which survival curves need to be plotted
    plot_genes <- c("ERK2")

    # Merge expression data with survival data
    expr_df <- normalized_counts %>%
      t() %>%
      data.frame() %>%
      tibble::rownames_to_column("Sample_ID") %>%
      dplyr::inner_join(meta_data, by=c("Sample_ID"="Sample_ID")) %>%
      dplyr::filter(Sample_ID %in% (meta_data %>% 
                                      #dplyr::filter(Received.platinum == "N") %>% 
                                      dplyr::filter(Received.platinum == "Y", Sample.collected.pre.platinum == "N") %>% 
                                      dplyr::select(Sample_ID) %>% 
                                      unlist(use.names=FALSE)))
    
    
    # Create a list to store survminer cutoffs, coxph stats, etc..
    stats <- list("gene" = c(),
                  "group" = c(),
                  "lower_cutoff" = c(),
                  "middle_cutoff" = c(),
                  "upper_cutoff" = c(),
                  "HR" = c(),
                  "CI_lower" = c(),
                  "CI_upper" = c(),
                  "pvalue" = c())
    
    classification_df <- expr_df %>% 
      dplyr::select(Sample_ID) %>%
      dplyr::mutate(Dummy_col = 0)
    
    # Plot survival curves
    for (gene in intersect(unique(plot_genes), rownames(normalized_counts))) {
      plot_curve <- TRUE
      summary <- wrangle_data(expr_df, stratify_criteria, prefix)
      
      class_df <- summary[[1]] %>%
        dplyr::select(Sample_ID, model) %>%
        dplyr::rename(!!gene := model)
      
      classification_df <- classification_df %>% 
        dplyr::left_join(class_df, by=("Sample_ID"="Sample_ID"))
      
      stats$gene          <- c(stats$gene, summary[[2]]$gene)
      stats$group         <- c(stats$group, summary[[2]]$group)
      stats$lower_cutoff  <- c(stats$lower_cutoff, summary[[2]]$lower)
      stats$middle_cutoff <- c(stats$middle_cutoff, summary[[2]]$middle)
      stats$upper_cutoff  <- c(stats$upper_cutoff, summary[[2]]$upper)
      stats$HR            <- c(stats$HR, summary[[2]]$HR )
      stats$CI_lower      <- c(stats$CI_lower, summary[[2]]$CI_lower)
      stats$CI_upper      <- c(stats$CI_upper, summary[[2]]$CI_upper)
      stats$pvalue        <- c(stats$pvalue, summary[[2]]$pvalue)
    }
    
    stats_df <- data.frame(stats)
    val_df <- normalized_counts[intersect(plot_genes, rownames(normalized_counts)),] %>%
      t() %>%
      data.frame()

    # Save the results
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "Summary")
    openxlsx::writeData(wb, sheet = "Summary", x = stats_df)
    openxlsx::addWorksheet(wb, sheetName = "Classification")
    openxlsx::writeData(wb, sheet = "Classification", x = classification_df)
    openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Individual_gene_stats.xlsx"), overwrite = TRUE)
  }
  
  #*****************USE THIS SECTION FOR PLOTTING GENE SIGNATURES**************#
  
  if (gene_signature == TRUE){
    
    # Human Y gene signature
    plot_genes <- c("DDX3Y", "EIF1AY", "HSFY2", "KDM5D", "UTY", "NLGN4Y", 
                    "PCDH11Y", "RPS4Y1", "TBL1Y", "TMSB4Y", "USP9Y", "ZFY", 
                    "DAZ1", "DAZ2", "DAZ3", "DAZ4", "PRY2", "RBMY1A1")
    
    # Calculate z-score using the function described in above paper
    # NOTE: z-score is calculated within each sample across all genes. So, 
    # number of samples doesnt affect a sample's z-score.
    expr_df <- as.data.frame(advanced_Z(plot_genes, normalized_counts))
    
    # Merge expression data with survival data
    expr_df <- expr_df %>%
      data.frame() %>%
      dplyr::rename(combined.exp = identity(1)) %>%
      tibble::rownames_to_column("Sample_ID") %>%
      dplyr::inner_join(meta_data, by=c("Sample_ID"="Sample_ID"))
    
    gene <- "combined.exp"
    
    #expr_df <- expr_df %>% dplyr::filter(Sex == "Male")
    if (nrow(expr_df > 0)){
      summary <- wrangle_data(expr_df, stratify_criteria, prefix)
    }
    
    surv_df <- summary[[1]] %>%
      dplyr::select(Sample_ID, combined.exp,	model, everything())
    stats_df <- as.data.frame(summary[[2]])
    
    # Save the results
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = paste0("Summary_",gse))
    openxlsx::writeData(wb, sheet = paste0("Summary_",gse), x = stats_df, rowNames = FALSE)
    openxlsx::addWorksheet(wb, sheetName = gse)
    openxlsx::writeData(wb, sheet = gse, x = surv_df, rowNames = FALSE)
    openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Summary.xlsx"), overwrite = TRUE)
  }
}
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Chao_summary.xlsx"), overwrite = TRUE)

#******************************************************************************#
#                      STEP 4: CHECK IF SIGNATURE IS GOOD                      #
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


#*****************USE THIS SECTION FOR PLOTTING GENE SIGNATURES****************#
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