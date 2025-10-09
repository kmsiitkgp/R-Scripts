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

full.df <- data.frame()
files.full <- list.files("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/custom", full.names = TRUE)
for (f in files.full){
  
  df <- read.xlsx(f) %>%
    dplyr::mutate(Comparison = gsub("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/custom/", "", f))
  
  full.df <- bind_rows(full.df, df)
}

full.df$Comparison <- gsub("\\.xlsx$|^Survival_Stats\\.", "",full.df$Comparison)
full.df <- full.df %>%
  tidyr::separate(
    col = Comparison,         # Column to split
    into = c("cutoff", "transform", "gene"),  # Names of new columns
    sep = "\\.",              # Delimiter (dot is special, so escape it)
    fill = "right",           # If fewer than 3 parts, fill with NA
    remove = TRUE             # Remove original column
  ) %>%
  dplyr::select(Gene, cutoff, transform, HR, pval, CI_lower, CI_upper, everything()) %>%
  arrange(Gene, cutoff, transform)

wb <- createWorkbook()
addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", full.df)
saveWorkbook(wb, "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/output.xlsx", overwrite = TRUE)
