source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

#************HPA Fig 4A**********

# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3059-z
# https://www.proteinatlas.org/about/download

BiocManager::install("HPAanalyze")
library("HPAanalyze")

HPAanalyze::hpaListParam()
data <- HPAanalyze::hpaDownload(downloadList = c("Pathology"), 
                                version = "latest")

# Rename the cancer to match TCGA
unique(data$pathology$cancer)

data$pathology$cancer[data$pathology$cancer == "glioma"] <- "GBM"
data$pathology$cancer[data$pathology$cancer == "melanoma"] <- "UVM"
data$pathology$cancer[data$pathology$cancer == "lymphoma"] <- "DLBC"
data$pathology$cancer[data$pathology$cancer == "skin cancer"] <- "SKCM"
data$pathology$cancer[data$pathology$cancer == "ovarian cancer"] <- "OV"
data$pathology$cancer[data$pathology$cancer == "liver cancer"] <- "LIHC"
data$pathology$cancer[data$pathology$cancer == "breast cancer"] <- "BRCA"
data$pathology$cancer[data$pathology$cancer == "testis cancer"] <- "TGCT"
data$pathology$cancer[data$pathology$cancer == "stomach cancer"] <- "STAD"
data$pathology$cancer[data$pathology$cancer == "thyroid cancer"] <- "THCA"
data$pathology$cancer[data$pathology$cancer == "cervical cancer"] <- "CESC"
data$pathology$cancer[data$pathology$cancer == "prostate cancer"] <- "PRAD"
data$pathology$cancer[data$pathology$cancer == "lung cancer"] <- "LUAD/LUSC"
data$pathology$cancer[data$pathology$cancer == "pancreatic cancer"] <- "PAAD"
data$pathology$cancer[data$pathology$cancer == "urothelial cancer"] <- "BLCA"
data$pathology$cancer[data$pathology$cancer == "endometrial cancer"] <- "UCEC"
data$pathology$cancer[data$pathology$cancer == "head and neck cancer"] <- "HNSC"
data$pathology$cancer[data$pathology$cancer == "colorectal cancer"] <- "COAD/READ"
data$pathology$cancer[data$pathology$cancer == "renal cancer"] <- "KICH/KIRC/KIRP"
data$pathology$cancer[data$pathology$cancer == "carcinoid"] <- "Carcinoid"

# Check if all cancers have been renamed to match TCGA
unique(data$pathology$cancer)

HPAanalyze::hpaVisPatho(data = data,
                        targetGene = "MAPK1",
                        customTheme = TRUE) + 
  theme_classic() +
  ggplot2::theme(plot.title   = element_text(family="sans", face="bold",  colour="black", size=15,  hjust=0.5),
                 plot.subtitle= element_text(family="sans", face="bold",  colour="black", size=12,  hjust=0.5),
                 plot.tag     = element_text(family="sans", face="bold",  colour="black", size=10,  hjust=0.5),
                 plot.caption = element_text(family="sans", face="bold",  colour="black", size=10,  hjust=0.5),
                 strip.background = element_blank(),
                 strip.text   = element_text(family="sans", face="bold",  colour="black", size=10,  hjust=0.5),
                 legend.title = element_text(family="sans", face="bold",  colour="black", size=12, hjust=0,   vjust=0.5,   angle=0),
                 legend.text  = element_text(family="sans", face="plain", colour="black", size=10, hjust=0,   vjust=0.5,   angle=0), 
                 axis.line.x  = element_blank(),
                 axis.line.y  = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.text.x  = element_text(family="sans", face="plain", colour="black", size=10, hjust=0.5, vjust=0.5, angle=45),
                 axis.text.y  = element_text(family="sans", face="plain", colour="black", size=10, hjust=0.5, vjust=0.5, angle=0),
                 axis.title.x = element_text(family="sans", face="bold",  colour="black", size=12, hjust=0.5, vjust=0,   angle=0),
                 axis.title.y = element_text(family="sans", face="bold",  colour="black", size=12, hjust=0.5, vjust=1,   angle=90)) +
  labs(y = "Proportion of Patients",
       x= "",
       fill = "Expression Level")

ggsave(filename = "Patho.tiff",
       width = 7,
       height = 7,
       units = c("in"),
       dpi = 600)

# HPAanalyze::hpaVisSubcell(targetGene = "MAPK1")
# ggsave("Location.tiff")
# 
# HPAanalyze::hpaVisTissue(targetGene = "MAPK1")
# ggsave("Tissue.tiff", height = 20, width = 10)

#***********Get protein expression data from TCPA********************# 

# Import read data
read_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "TCGA-PANCAN32-L4_original.xlsx")) %>%
  tibble::column_to_rownames("Sample_ID") %>%
  dplyr::select(everything(), -c("Cancer_Type", "Sample_Type")) %>%
  t() %>%
  data.frame() %>%
  tibble::rownames_to_column("SYMBOL")
colnames(read_data)[1] <- "SYMBOL"

# Import meta data
meta_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Meta_data_TCGA.xlsx"))

# Rename column names to match with Sample_ID from meta_data
cols <- colnames(read_data)

for (i in 1:ncol(read_data)){
  
  if(sum(stringr::str_detect(string=colnames(read_data)[i],
                             pattern=make.names(meta_data$Sample_ID))) == 1){
    cols[i] <- meta_data$Sample_ID[stringr::str_detect(string=colnames(read_data)[i],
                                                       pattern=make.names(meta_data$Sample_ID))]
  }
}

# Replace with proper Sample_ID as column names
colnames(read_data) <- cols

# Check for duplicated columns
sum(duplicated(colnames(read_data)))

# Make them unique
colnames(read_data) <- make.names(colnames(read_data), unique=TRUE)

# Replace all NA values with 0
read_data <- read_data %>%
  base::replace(is.na(.), 0)

# Remove genes with no expression in all samples
read_data <- read_data[rowSums(read_data[,c(-1)]) !=0,]

# Save the reformatted read_data
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName="Normalized")
openxlsx::writeData(wb, sheet="Normalized", x=read_data)
openxlsx::saveWorkbook(wb, file=paste0(parent_path, "TCGA-PANCAN32-L4.xlsx"), overwrite=TRUE)

#***********Survival for Fig 4A using protein data**********

subset_group <- NA
subset_value <- NA
plot_by <- "Expression"
split_by <- subset_group
combine_plot <- FALSE 
multiple_cutoff <- TRUE
stratify_criteria <- "o"  # "m"
reference <- "T1" #"NA"
reference <- dplyr::if_else(plot_by == "Expression", "LOW", reference)
confidence_interval <- FALSE
legend_title <- "Expression"
plot_risk_table <- TRUE
plot_curve <- TRUE
all_quartiles <- FALSE
gene_signature <- FALSE
color_palette <- c("#d73027","#0c2c84")
color_palette <- c(color_palette, 
                   colorspace::adjust_transparency(col = color_palette, alpha = 0.3), 
                   colorspace::adjust_transparency(col = color_palette, alpha = 0.6))
variable_x <- "Time"
variable_y <- "Status"
parent_path <- "C:/Users/kailasamms/Desktop/Collaboration projects/Mukta/" 


# Import read_data
read_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "TCGA-PANCAN32-L4.xlsx"))
colnames(read_data)[1] <- "SYMBOL"

# NOTE: The normalized protein expression from TCPA already seems to be median 
# centered etc as some values are negative. Also, the values are log transformed
# given the max value is  below 100
max(read_data[,-1])

# Import meta data
meta_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Meta_data_TCGA.xlsx")) 

# Reformat metadata 
meta_data <- meta_data %>% 
  dplyr::mutate(Sample_ID = make.names(names = Sample_ID)) %>%
  dplyr::mutate(Time = as.numeric(Time)) %>%
  dplyr::filter(Time > 0 & !is.na(Time)) %>%
  dplyr::distinct_at("Sample_ID", .keep_all = TRUE)

# Reformat read data
normalized_counts <- read_data %>%
  dplyr::mutate(SYMBOL = make.names(names = SYMBOL, unique = TRUE)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  tibble::column_to_rownames(var = "SYMBOL") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
colnames(normalized_counts) <- base::make.names(names = colnames(normalized_counts))
normalized_counts <- normalized_counts[,intersect(make.names(meta_data$Sample_ID), colnames(normalized_counts))]
normalized_counts <- normalized_counts[!rowSums(normalized_counts, na.rm=TRUE) == 0,]


if (gene_signature == FALSE){
  
  # Create a list of genes for which survival curves need to be plotted
  plot_genes <- c("ERK2")
  
  # Merge expression data with survival data
  expr_df <- normalized_counts %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::inner_join(meta_data, by=c("Sample_ID"="Sample_ID"))
  
  # Plot survival curves for each cancer
  for (proj in unique(expr_df$Project_ID)){
    
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
    
    for (gene in intersect(unique(plot_genes), rownames(normalized_counts))) {
      
      prefix <- proj
      
      df <- expr_df %>%
        dplyr::filter(Project_ID == proj)
      
      summary <- wrangle_data(df, stratify_criteria, prefix)
      
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
    
    # # Save the results
    # wb <- openxlsx::createWorkbook()
    # openxlsx::addWorksheet(wb, sheetName = "Summary")
    # openxlsx::writeData(wb, sheet = "Summary", x = stats_df)
    # openxlsx::addWorksheet(wb, sheetName = "Classification")
    # openxlsx::writeData(wb, sheet = "Classification", x = classification_df)
    # openxlsx::saveWorkbook(wb, file = paste0(parent_path, prefix, "_Individual_gene_stats.xlsx"), overwrite = TRUE)
  }
}

#************************************

# PCA plots of CRISPR data from bowtie, bowtie2, mageck and kms counts

# Path to directory containing input files or folders
parent_path <- paste0("C:/Users/kailasamms/Box/Saravana@cedars/05. Bioinformatics/CRISPR/Jinfen/")
results_path <- parent_path

species <- "Mus musculus"
padj.cutoff <- 0.1
lfc.cutoff <- 0 
heatmap_plot <- FALSE
volcano_plot <- FALSE
cor_plot  <- FALSE
Variable <- "Culture"
Comparisons <- list(Target=c("3D"),
                    Reference=c("2D"))

annotations <- get_annotations(species)

for (suffix in c("bowtie", "bowtie2", "mageck", "kms")){
  
  meta_data <- openxlsx::read.xlsx(xlsxFile=paste0(parent_path, "Meta_data.xlsx"))
  
  #read_data <- openxlsx::read.xlsx(xlsxFile=paste0(parent_path, "Read_data.xlsx"))
  read_data <- read.table(file=paste0(parent_path, "count_results/", suffix, "/", suffix, ".count.txt"), header=TRUE)
  read_data <- read_data[,-2]
  colnames(read_data)[1] <- "SYMBOL"
  
  # Remove non-targeting guides
  read_data <- read_data %>% dplyr::filter(!grepl(pattern="non_targeting",x=SYMBOL))
  
  meta_data <- prep_metadata(meta_data, Variable)
  read_data <- prep_readdata(read_data, meta_data)
  l <- check_data(read_data, meta_data)
  meta_data <- l[[2]]
  read_data <- l[[1]]
  
  for (n in 1:length(Comparisons$Target)){
    
    # Perform DESeq2() using in-built batch modelling
    approach <- "DESeq2_modelled"
    
    if (length(unique(meta_data$Batch)) > 1){
      dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                            colData=meta_data, 
                                            design=~ Batch+id)
    } else {
      dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                            colData=meta_data, 
                                            design=~ id)
    }
    
    dds <- run_deseq2(dds, meta_data, annotations, Comparisons, n, lfc.cutoff, padj.cutoff, approach, suffix)
    deseq2_norm_counts(dds, annotations, approach, suffix) # batch corrected if you more than 1 batch
    plot_qc(dds, meta_data, approach, suffix)
  }
}

#**********Identifying good targets

# Read the common essential genes
common_essential <- read.xlsx("C:/Users/kailasamms/Box/Saravana@cedars/05. Bioinformatics/Depmap Gene Dependency Profile Summary.xlsx",
                              sheet = "Common Essential Human Genes")
common_essential <- common_essential$All.combined

for (suffix in c("bowtie", "bowtie2", "mageck", "kms")){
  for (comparison in c("Day7vs0_2D", "Day14vs0_2D", "Day7vs0_3D", "Day14vs0_3D")){
    
    day7_2D <- read.table(paste0(parent_path, "DEG_results/", suffix, "/", "Day7vs0_2D.gene_summary.txt"), 
                          header=TRUE)
    day7_3D <- read.table(paste0(parent_path, "DEG_results/", suffix, "/", "Day7vs0_3D.gene_summary.txt"), 
                          header=TRUE)
    day14_2D <- read.table(paste0(parent_path, "DEG_results/", suffix, "/", "Day14vs0_2D.gene_summary.txt"), 
                           header=TRUE)
    day14_3D <- read.table(paste0(parent_path, "DEG_results/", suffix, "/", "Day14vs0_3D.gene_summary.txt"), 
                           header=TRUE)
    day0 <- read.table(paste0(parent_path, "DEG_results/", suffix, "/", "Day0_3Dvs2D.gene_summary.txt"), 
                       header=TRUE)
    day7 <- read.table(paste0(parent_path, "DEG_results/", suffix, "/", "Day7_3Dvs2D.gene_summary.txt"), 
                       header=TRUE)
    day14 <- read.table(paste0(parent_path, "DEG_results/", suffix, "/", "Day14_3Dvs2D.gene_summary.txt"), 
                        header=TRUE)
    
    day7_2D <- day7_2D %>% 
      dplyr::filter(!(id %in% common_essential)) %>%
      dplyr::filter(neg.fdr < 0.05 | pos.fdr < 0.05)
    
    day7_3D <- day7_3D %>% 
      dplyr::filter(!(id %in% common_essential)) %>%
      dplyr::filter(neg.fdr < 0.05 | pos.fdr < 0.05)
    
    day14_2D <- day14_2D %>% 
      dplyr::filter(!(id %in% common_essential)) %>%
      dplyr::filter(neg.fdr < 0.05 | pos.fdr < 0.05)
    
    day14_3D <- day14_3D %>% 
      dplyr::filter(!(id %in% common_essential)) %>%
      dplyr::filter(neg.fdr < 0.05 | pos.fdr < 0.05)
    
    day0 <- day0 %>% 
      dplyr::filter(!(id %in% common_essential)) %>%
      dplyr::filter(neg.fdr < 0.05 | pos.fdr < 0.05)
    
    day7 <- day7 %>% 
      dplyr::filter(!(id %in% common_essential)) %>%
      dplyr::filter(neg.fdr < 0.05 | pos.fdr < 0.05)
    
    day14 <- day14 %>% 
      dplyr::filter(!(id %in% common_essential)) %>%
      dplyr::filter(neg.fdr < 0.05 | pos.fdr < 0.05)
    
    if ("MAPK1" %in% c(day7_2D$id, day7_3D$id, day14_2D$id, day14_3D$id, day0, day7, day14)){
      print("present")
    }
  }
}

# Unfortunately MAPK1 is not present in any of the comparisons.
    
#**************** Fig 5A

## Store path of parent directory i.e. root directory for the project
parent_path <- "C:/Users/KailasammS/Box/Saravana@cedars/10. Ongoing Projects/BBN project/BLCA_Cohorts_correct/"

gse <- "A"
group <- "NA"
survival_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "Mukta 17 gene expression.xlsx"))
plot_survival(survival_data, group)

#**************** Fig 5B signature

subset_group <- NA
subset_value <- NA
plot_by <- "Expression"
split_by <- subset_group
combine_plot <- FALSE 
multiple_cutoff <- FALSE
stratify_criteria <- "o"
reference <- "NA"
reference <- dplyr::if_else(plot_by == "Expression", "LOW", reference)
confidence_interval <- FALSE

legend_title <- "Expression"  #Sex #Race
plot_risk_table <- FALSE
plot_curve <- TRUE
all_quartiles <- FALSE
gene_signature <- TRUE
legend_label <- c("High", "Low")              # Reference is Low
color_palette <- c("#DB6D00", "#490092")      # 1st color~high, 2nd color~low
variable_x <- "Time"
variable_y <- "Status"

gse <- "Imvigor210" 

plot_genes <- c("USP18","IFIT5","IFI44L","PARP9","OAS2","CMPK2","OAS3",
                "CXCL16","SP110","IFIH1","HERC6","IFIT1","RSAD2","SP100",
                "XAF1","STAT1","IFI6")

expr_df <- as.data.frame(advanced_Z(plot_genes, normalized_counts))

# Merge expression data with survival data
expr_df <- expr_df %>%
  data.frame() %>%
  dplyr::rename(combined.exp = identity(1)) %>%
  tibble::rownames_to_column("GEO_ID") %>%
  dplyr::inner_join(meta_data, by=c("GEO_ID"="GEO_ID"))

gene <- "combined.exp"
if (nrow(expr_df > 0)){
  summary <- wrangle_data(stratify_criteria)
}

surv_df <- summary[[1]]
stats_df <- as.data.frame(summary[[2]])

# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = paste0("Summary_",gse))
openxlsx::writeData(wb, sheet = paste0("Summary_",gse), x = stats_df, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName = gse)
openxlsx::writeData(wb, sheet = gse, x = surv_df, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0("C:/Users/KailasammS/Desktop/", "a111.xlsx"), overwrite = TRUE)

#**************** Fig 5B MAPK

gse <- "Imvigor210" 
gene_signature <- FALSE
plot_genes <- "MAPK1"

#**************** Fig 2G  "MAPK1", "AKR1A1", "ZMYND8", "SLC16A1"

gse <- "TCGA_BLCA" 
gene_signature <- FALSE
plot_genes <- c("MAPK1", "AKR1A1", "ZMYND8", "SLC16A1")

# Merge expression data with survival data
expr_df <- normalized_counts %>%
  t() %>%
  data.frame() %>%
  tibble::rownames_to_column("GEO_ID") %>%
  dplyr::inner_join(meta_data, by=c("GEO_ID"="GEO_ID")) %>% 
  dplyr::filter(Time < 2500)

####### Fig 1G t score and plot scatter plot

calc_t_score <- function(data){
  
  data_control <- data %>%
    dplyr::filter(Gene %in% c("none", "safe"))
  
  median_ctrl <- median(data_control$LFC)
  sd_ctrl <- sd(data_control$LFC)
  
  data <- data %>%
    dplyr::mutate(pZ = (LFC-median_ctrl)/sd_ctrl)
  
  data_control <- data_control %>%
    dplyr::mutate(pZ = (LFC-median_ctrl)/sd_ctrl)
  
  U_ctrl <- median(data_control$pZ)
  Var_ctrl <- var(data_control$pZ)
  N_ctrl <- nrow(data_control)
  
  data <- data %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(U_gene = median(pZ),
                  Var_gene = var(pZ),
                  N_gene = n(),
                  S_gene = (Var_gene*(N_gene-1)) + (Var_ctrl*(N_ctrl-1)),
                  t_score = abs((U_gene - U_ctrl)/sqrt(S_gene/N_gene + S_gene/N_ctrl)),
                  U_ctrl = U_ctrl,
                  Var_ctrl = Var_ctrl,
                  N_ctrl = N_ctrl) %>%
    dplyr::select(Gene, U_gene, Var_gene, N_gene, U_ctrl, Var_ctrl, N_ctrl, S_gene, t_score) %>%
    dplyr::distinct_at("Gene", .keep_all = TRUE)
  
  return(data)
}

plot_scatter <- function(data, file){
  
  color_breaks <- c(-20,0,20)
  p <- ggplot2::ggplot(data = data,
                       aes(x = U_gene, 
                           y = t_score,
                           #size = t_score,
                           #color = pz,
                           fill = U_gene)) +
    # Plot dot plot
    ggplot2::geom_point(col="black", 
                        shape=21,
                        stroke=0.5,
                        position=position_jitter(h=0.01,w=0.01)) +
    # Define the theme of plot
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "U_gene") +
    coord_cartesian(xlim = c(-5,3), clip = "off") +
    scale_x_continuous(breaks = seq(-5, 3, by = 1)) +
    #scale_y_continuous(breaks = seq(0, 0.5, by = 0.1)) +
    ggplot2::guides(size = "none",
                    fill = guide_colourbar(theme = theme(legend.key.width  = unit(0.75, "lines"),
                                                         legend.key.height = unit(10, "lines"),
                                                         legend.ticks = element_blank(),
                                                         legend.frame = element_rect(colour = "Black",
                                                                                     linewidth = 0.5)))) +
    # Define the color of the dots
    ggplot2::scale_fill_viridis_c(option="C", limits =c(-5,5))
    # scale_fill_gradientn(colors=c("#007ba7", "Black","#FFFF00"), 
    #                      limits=c(-20, 20), 
    #                      values=c(0, scales::rescale(color_breaks, from = range(color_breaks)), 1))
    #scale_fill_gradient2(low="#007ba7", mid="Black", high="Yellow", midpoint = 0, limits=c(-5, 2))
    
  #scale_fill_continuous_diverging(palette = "Tofino")
  
  ggsave(paste0("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/", file, ".jpg"))
  return(p)
}

wb <- openxlsx::createWorkbook()
for (file in c("2D_Day7vs0.sgrna_summary.txt", "2D_Day14vs0.sgrna_summary.txt",
               "3D_Day7vs0.sgrna_summary.txt", "3D_Day14vs0.sgrna_summary.txt")){
               #"3D_Day14vs2D_Day14.sgrna_summary.txt")){
  
  data <- read.table(paste0("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/", file), header = TRUE) %>%
    dplyr::select(sgrna, Gene, control_mean, treat_mean, LFC)
  
  data <- calc_t_score(data)
  plot_scatter(data, file)
  
  t <- gsub(pattern=".sgrna_summary.txt", replacement="", x=file)
  openxlsx::addWorksheet(wb, sheetName = t)
  openxlsx::writeData(wb, sheet = t, x = data, rowNames = FALSE)  
}
openxlsx::saveWorkbook(wb, file = paste0("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/", "t_score.xlsx"), overwrite = TRUE)


#### Fig 2D Venndiagram vs upset plot

# Read file prepared fro venn diagram
d7_2d <- read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Mukta/2D_Day7vs0.gene_summary.xlsx")
d14_2d <- read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Mukta/2D_Day14vs0.gene_summary.xlsx")
d7_3d <- read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Mukta/3D_Day7vs0.gene_summary.xlsx")
d14_3d <- read.xlsx("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Mukta/3D_Day14vs0.gene_summary.xlsx")

colnames(d7_2d) <- make.names(colnames(d7_2d))
colnames(d7_3d) <- make.names(colnames(d7_3d))
colnames(d14_2d) <- make.names(colnames(d14_2d))
colnames(d14_3d) <- make.names(colnames(d14_3d))

d7_2d <- d7_2d %>%
  dplyr::filter(neg.goodsgrna >= 3, neg.lfc <= -0.1) %>%
  dplyr::select(id) %>%
  unlist(use.names=FALSE)
d7_2d <- setdiff(d7_2d, "none")
d7_3d <- d7_3d %>%
  dplyr::filter(neg.goodsgrna >=3, neg.lfc <= -0.1) %>%
  dplyr::select(id) %>%
  unlist(use.names=FALSE)
d7_3d <- setdiff(d7_3d, "none")
d14_2d <- d14_2d %>%
  dplyr::filter(neg.goodsgrna >=3, neg.lfc <= -0.1) %>%
  dplyr::select(id) %>%
  unlist(use.names=FALSE)
d14_2d <- setdiff(d14_2d, "none")
d14_3d <- d14_3d %>%
  dplyr::filter(neg.goodsgrna >=3, neg.lfc <= -0.1) %>%
  dplyr::select(id) %>%
  unlist(use.names=FALSE)
d14_3d <- setdiff(d14_3d, "none")

listInput <- list(Day.7.2D=d7_2d,
                  Day.14.2D=d14_2d,
                  Day.7.3D=d7_3d,
                  Day.14.3D=d14_3d)
                  
# Plot
upset(data = UpSetR::fromList(listInput),
      empty.intersections = "on",
      cutoff = 5,
      sets = c("Day.7.2D","Day.14.2D","Day.7.3D","Day.14.3D"),
      keep.order = TRUE,
      #nintersects = 5,                  # number of groups on X axis
      #nsets = 21,                         # number of groups on Y axis
      order.by = c("freq"),
      mb.ratio=c(0.7,0.3),
      point.size = 6,
      line.size=1,
      sets.bar.color=c("#440154FF", "#31688EFF","#35B779FF","#FDE725FF"))


max_l <- max(lengths(listInput))
for (i in 1:length(listInput)){
  listInput[[i]] <- c(listInput[[i]], rep(x="",times=max_l - length(listInput[[i]])))
}
path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/"
suffix <- ""
plot_venn(data.frame(listInput), path, suffix)




# scales::show_col(viridis(4))
"#440154FF", "#31688EFF","#35B779FF","#FDE725FF"
scales::show_col(plasma(4))
"#0D0887FF", "#9C179EFF", "#ED7953FF", "#F0F921FF" 


#************************Normalized counts

parent_path <- "C:/Users/KailasammS/Box/Saravana@cedars/10. Ongoing Projects/BBN project/BLCA_Cohorts_correct/"

read_data <- read.xlsx(paste0(parent_path, "TCGA_BLCA_Normalized.xlsx"))
colnames(read_data)[1] <- "SYMBOL"

meta_data <- read.xlsx(paste0(parent_path, "TCGA_BLCA_Metadata.xlsx"))
genes <- c("CDKL1", "CHRNA9", "DCLK1", "DYRK2", "GABRQ", "GRIK5",
           "HIPK4", "MAP2K2", "MAPK1", "MAPK7", "ROCK1")

# Reformat read data
normalized_counts <- read_data %>%
  dplyr::mutate(SYMBOL = make.names(names = SYMBOL, unique = TRUE)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  tibble::column_to_rownames(var = "SYMBOL") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
colnames(normalized_counts) <- base::make.names(names = colnames(normalized_counts))
normalized_counts <- normalized_counts[,intersect(make.names(meta_data$Sample_ID), colnames(normalized_counts))]
normalized_counts <- normalized_counts[!rowSums(normalized_counts, na.rm=TRUE) == 0,]
normalized_counts <- log(1+normalized_counts, base=2)
t <- base::apply(X=normalized_counts, MARGIN=1, FUN=median, na.rm=TRUE)
normalized_counts <- base::sweep(x=normalized_counts, MARGIN=1, FUN="-", STATS=t)

df1 <- normalized_counts[rownames(normalized_counts) %in% genes,] %>%
  t() %>%
  data.frame() %>%
  tibble::rownames_to_column("Sample_ID") %>%
  dplyr::left_join(meta_data %>% dplyr::select(Sample_ID, Stage) %>% dplyr::mutate(Sample_ID = make.names(Sample_ID)),
                   by = c("Sample_ID"="Sample_ID"))

df2 <- read_data %>%
  dplyr::filter(SYMBOL %in% genes) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  t() %>%
  data.frame() %>%
  tibble::rownames_to_column("Sample_ID") %>%
  dplyr::left_join(meta_data %>% dplyr::select(Sample_ID, Stage),
                   by = c("Sample_ID"="Sample_ID"))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "median centered log normalized")
openxlsx::writeData(wb, sheet = "median centered log normalized", x = df1)
openxlsx::addWorksheet(wb, sheetName = "normalized")
openxlsx::writeData(wb, sheet = "normalized", x = df2)
openxlsx::saveWorkbook(wb, file = paste0("TCGA_counts.xlsx"), overwrite = TRUE)






