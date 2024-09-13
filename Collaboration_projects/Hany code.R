#*********LOY score for each patient and correlation with ANPEPPS

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")
parent_path <- "/hpc/home/kailasamms/scratch/TCGA_GDC/"

read_data <- utils::read.table(file=paste0(parent_path, "Normalized_Counts_TCGA.tsv"),
                               header=TRUE, sep="\t", quote="", skip=0, fill=TRUE)
colnames(read_data)[1] <- "SYMBOL"

meta_data <- read.xlsx(paste0(parent_path, "Meta_data_TCGA.xlsx"))
proj <- unique(meta_data$Project_ID)

# Create dummy dataframe to store YScore and ANPEP expression from each cancer
df <- data.frame(Sample_ID=c(), YScore=c(), ANPEP=c(), TCGA=c())

for (f in proj){
  
  meta_data <- read.xlsx(paste0(parent_path, "Meta_data_", f, ".xlsx"))
  meta_data <- meta_data %>%
    dplyr::mutate(Sample_ID = make.names(Sample_ID, unique=TRUE)) %>%
    dplyr::filter(Sample_ID %in% colnames(read_data)) %>%
    dplyr::filter(Sex == "Male")
  
  raw_data <- read_data %>%
    dplyr::select(SYMBOL, meta_data$Sample_ID) %>% 
    dplyr::filter(SYMBOL %in% c("ANPEP")) %>% 
    tibble::column_to_rownames("SYMBOL") %>% 
    t() %>% 
    data.frame() %>% 
    tibble::rownames_to_column("Sample_ID")

  # log normalize and median center the counts
  normalized_counts <- read_data %>%
    dplyr::select(SYMBOL, meta_data$Sample_ID) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("SYMBOL")
    
  
  normalized_counts <- log(1+normalized_counts, base=2)
  t <- base::apply(X=normalized_counts, MARGIN=1, FUN=median, na.rm=TRUE)
  normalized_counts <- base::sweep(x=normalized_counts, MARGIN=1, FUN="-", STATS=t)
  
  plot_genes <- c("DDX3Y", "UTY", "KDM5D", "USP9Y", "ZFY", "RPS4Y1", "TMSB4Y", 
           "EIF1AY", "NLGN4Y", "TBL1Y", "HSFY2", "PCDH11Y", "PRY2", "DAZ1",
           "DAZ2", "DAZ3", "DAZ4", "RBMY1A1")
  
  expr_df <- as.data.frame(advanced_Z(plot_genes, normalized_counts))
  colnames(expr_df) <- "YScore"
  expr_df <- expr_df %>%
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::left_join(raw_data, by=c("Sample_ID"="Sample_ID")) %>%
    dplyr::mutate(TCGA=f)
  
  # Plot correlation
  library("ggpubr")
  df1 <- expr_df %>% 
    dplyr::mutate(log2ANPEP = log2(1+ANPEP))
  ggscatter(data = df1, x = "log2ANPEP", y = "YScore", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "ANPEP", ylab = "Y Score")
  
  ggsave(paste0(f, ".jpg"))
  
  #df <- dplyr::bind_rows(df, expr_df)
  print(f)
}

# Save Yscore and ANPEP expression for all patients
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "YScore")
openxlsx::writeData(wb, sheet = "YScore", x = df)
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Y_Scores.xlsx"), overwrite = TRUE)



