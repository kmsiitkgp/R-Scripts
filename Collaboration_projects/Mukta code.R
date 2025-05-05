source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
parent_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Mukta/"

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
       path = parent_path,
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
    
    # Save the results
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "Summary")
    openxlsx::writeData(wb, sheet = "Summary", x = stats_df)
    openxlsx::addWorksheet(wb, sheetName = "Classification")
    openxlsx::writeData(wb, sheet = "Classification", x = classification_df)
    openxlsx::saveWorkbook(wb, file = paste0(parent_path, prefix, "_Individual_gene_stats.xlsx"), overwrite = TRUE)
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

############ Fig 1G t score

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Mukta/"

wb <- openxlsx::createWorkbook()
for (file in c("2D_Day7vs0.sgrna_summary.txt", "2D_Day14vs0.sgrna_summary.txt",
               "3D_Day7vs0.sgrna_summary.txt", "3D_Day14vs0.sgrna_summary.txt")){
  
  data <- read.table(paste0(data_path, file), header = TRUE) %>%
    dplyr::select(sgrna, Gene, control_mean, treat_mean, LFC)
  
  data <- calc_t_score(data)
  suffix <- gsub(pattern=".sgrna_summary.txt", replacement="", x=file)
  plot_t_score(data, data_path, suffix)
  
  openxlsx::addWorksheet(wb, sheetName = suffix)
  openxlsx::writeData(wb, sheet = suffix, x = data, rowNames = FALSE) 
}
openxlsx::saveWorkbook(wb, file = paste0(data_path, "t_score.xlsx"), overwrite = TRUE)


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

parent_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Mukta/"
suffix <- ""
plot_venn(data.frame(listInput), parent_path, suffix)

# scales::show_col(viridis(4))
"#440154FF", "#31688EFF","#35B779FF","#FDE725FF"
scales::show_col(plasma(4))
"#0D0887FF", "#9C179EFF", "#ED7953FF", "#F0F921FF" 


#************************Normalized counts

parent_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Mukta/"
input_path <- "C:/Users/KailasammS/Box/Saravana@cedars/10. Ongoing Projects/BBN project/BLCA_Cohorts_correct/"

read_data <- read.xlsx(paste0(input_path, "TCGA_BLCA_Normalized.xlsx"))
colnames(read_data)[1] <- "SYMBOL"

meta_data <- read.xlsx(paste0(input_path, "TCGA_BLCA_Metadata.xlsx"))
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
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "TCGA_counts.xlsx"), overwrite = TRUE)


###########Reviewer comments

#BiocManager::install("GEOquery")
library(GEOquery)
library(limma)
library(umap)

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

parent_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Mukta/"
dataset <- c("GSE3167") #, "GSE138118")

for (d in dataset){
  
  # Extract raw data from CEL files
  cel_files <- oligoClasses::list.celfiles(paste0(parent_path, "Microarray/", d, "_RAW"), full.names=TRUE, listGzipped=TRUE)
  raw_data <- oligo::read.celfiles(cel_files)
  
  # Perform RMA (background subtraction, quantile normalization, summarization)
  normalized_data <- oligo::rma(object = raw_data,
                                background = TRUE,
                                normalize = TRUE) 
  
  # View distribution before and after normalization
  # oligo::hist(raw_data)
  # oligo::hist(normalized_data)
  # oligo::boxplot(x = raw_data, which=c("pm"), transfo=log2)
  # par(mar=c(0,0,0,0)) # use it to avoid " figure margins too large" error
  # oligo::boxplot(x = normalized_data, which=c("pm"), transfo=log2)
  
  # Extract normalized counts
  normalized_counts <- exprs(normalized_data)
  colnames(normalized_counts) <- gsub(pattern="\\..*", replacement="", x= colnames(normalized_counts))
  colnames(normalized_counts) <- gsub(pattern="_.*", replacement="", x= colnames(normalized_counts))
  
  # metadata
  gene_data <- read.table(paste0(parent_path, "Microarray/", d, "_GENE.txt"), sep="\t", header=TRUE)
  meta_data <- read.xlsx(paste0(parent_path, "Microarray/", d, "_CLINICAL.xlsx"))
  
  #MAPK1
  probes <- gene_data %>% 
    dplyr::filter(Gene %in% c("MAPK1", "GABRQ", "CHRNA9", "GRIK5", "DYRK2", 
                              "MAP2K2", "HIPK4", "MAPK7", "ROCK1", "CDKL1",
                              "DCLK1")) %>% 
    dplyr::mutate(Gene = make.names(Gene, unique=TRUE))
    # dplyr::select(ID) %>% 
    # unlist(use.names=FALSE)
  
  # mapk1_probes <-  gene_data %>% 
  #   dplyr::filter(Gene %in% c("MAPK1")) %>% 
  #   dplyr::select(ID) %>% 
  #   unlist(use.names=FALSE)
  
  mapk1 <- normalized_counts[rownames(normalized_counts) %in% probes$ID,] %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column("GSM") %>%
    dplyr::left_join(meta_data %>% dplyr::select(GSM, Group), by=c("GSM"="GSM"))
  
  cols <- c("GSM")
  for (i in 2:ncol(mapk1)){
    cols <- c(cols, probes$Gene[which(make.names(probes$ID) == colnames(mapk1)[i])])
  }
  cols <- c(cols, "Group")
  
  colnames(mapk1) <- cols

  # cols <- c("GSM", "Group")
  # for (i in 1:nrow(probes)){
  #   cols <- c(cols, colnames(mapk1)[str_detect(colnames(mapk1), make.names(mapk1_probes[i]))])
  # }

  # mapk1 <- mapk1 %>%
  #   dplyr::select(all_of(cols))
    
  # Save
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "normalized counts")
  openxlsx::writeData(wb, sheet = "normalized counts", x = normalized_counts, rowNames = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "probe-gene mapping")
  openxlsx::writeData(wb, sheet = "probe-gene mapping", x = gene_data)
  openxlsx::addWorksheet(wb, sheetName = "sample-group mapping")
  openxlsx::writeData(wb, sheet = "sample-group mapping", x = meta_data)
  openxlsx::addWorksheet(wb, sheetName = "MAPK1")
  openxlsx::writeData(wb, sheet = "MAPK1", x = mapk1)
  openxlsx::saveWorkbook(wb, file = paste0(parent_path, d, ".xlsx"), overwrite = TRUE)
}

################## CCLE classification neuroendocrine vs non-neuroendocrine

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Mukta/"

# # (DO THIS ONLY ONCE) Generate DESeq2 normalized counts 
# read_data <- read.xlsx(paste0(data_path, "CCLE_Read_data.xlsx"))
# 
# # If there are duplicated genes, keep only data for highest expressing copy
# read_data <- read_data %>%
#   dplyr::select(everything(), -Name) %>%
#   dplyr::rename(SYMBOL = Description) %>%
#   dplyr::mutate(n = rowSums(.[,-1])) %>%
#   dplyr::group_by(SYMBOL) %>%
#   dplyr::slice_max(n) %>%
#   dplyr::ungroup() %>%
#   # Duplicated genes with 0 expression in all samples still remain, remove them
#   dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
#   dplyr::select(everything(), -n)
# colnames(read_data) <- gsub(pattern = "_.*", replacement = "", x=colnames(read_data))
# colnames(read_data) <- make.names(colnames(read_data), unique = TRUE)
# 
# meta_data <- read.xlsx(paste0(data_path, "CCLE_Meta_data.xlsx"))
# meta_data <- meta_data %>%
#   dplyr::mutate(Sample_ID = gsub(pattern = "_.*", replacement = "", x=meta_data$CCLE_ID)) %>%
#   dplyr::select(Sample_ID, everything())
# 
# # Get annotations
# species <- "Homo sapiens"
# annotations <- get_annotations(species)
# 
# 
# Comparisons <- list(Variable =c(NA),
#                     Target   =c(NA),
#                     Reference=c(NA)) 
# meta_data <- prep_metadata(meta_data, read_data)
# read_data <- prep_readdata(read_data, meta_data)
# l <- check_data(read_data, meta_data)
# meta_data <- l[[2]]
# read_data <- l[[1]]
# 
# # Normalize the raw read counts
# dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
#                                       colData=meta_data, 
#                                       design=~ 1)
# 
# approach <- ""
# suffix <- ""
# deseq2_norm_counts(dds, annotations, approach, suffix, data_path)

# Read normalized counts
read_data <- read.xlsx(paste0(data_path, "CCLE_Normalized_Counts.xlsx"))
plot_genes <- c("INSM1", "SYP", "MEIS2", "CHGA", "MKI67", "NCAM1")

# Reformat read data
normalized_counts <- read_data[, -c(2:5)] %>%
  dplyr::mutate(SYMBOL = make.names(names = SYMBOL, unique = TRUE)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  tibble::column_to_rownames(var = "SYMBOL") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
colnames(normalized_counts) <- base::make.names(names = colnames(normalized_counts))
normalized_counts <- normalized_counts[!rowSums(normalized_counts, na.rm=TRUE) == 0,]

normalized_counts <- log(1+normalized_counts, base=2) # log transform the normalized counts so it is not skewed
t <- base::apply(X=normalized_counts, MARGIN=1, FUN=median, na.rm=TRUE)  # calculate median for each gene across all samples
normalized_counts <- base::sweep(x=normalized_counts, MARGIN=1, FUN="-", STATS=t) # subtract median of each gene

expr_df <- as.data.frame(advanced_Z(plot_genes, normalized_counts))  # calculate mean of marker genes, subtract mean of all genes
expr_df <- expr_df %>%
  dplyr::rename(Score = identity(1)) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::left_join(normalized_counts %>% 
                     t() %>% 
                     data.frame() %>% 
                     tibble::rownames_to_column("Sample") %>%
                     dplyr::select(Sample, CDH12),
                   by=c("Sample"="Sample")) %>%
  dplyr::left_join(read_data %>%
                     dplyr::select(everything(), -contains(c("ENTREZ", "ENSEMBL"))) %>%
                     tibble::column_to_rownames("SYMBOL") %>%
                     t() %>%
                     data.frame() %>%
                     tibble::rownames_to_column("Sample") %>%
                     dplyr::select(Sample, CDH12),
                   by=c("Sample"="Sample")) %>%
  dplyr::rename(log2_CDH12 = identity(3), CDH12 = identity(4))

# Identify cutoffs and group the cell lines in neuro low and neuro high
pos_controls <- c("SKNAS", "SKNFI", "NCIH727", "TT.1")
neg_controls <- c("A427","A172","SW480")

low_cutoff <- max(expr_df %>% dplyr::filter(Sample %in% neg_controls) %>% dplyr::select(Score))
high_cutoff <- min(expr_df %>% dplyr::filter(Sample %in% pos_controls) %>% dplyr::select(Score))

neuro_low <- expr_df %>% dplyr::filter(Score <= low_cutoff) %>% dplyr::select(Sample) %>% unlist(use.names = FALSE)
neuro_high <- expr_df %>% dplyr::filter(Score >= high_cutoff) %>% dplyr::select(Sample) %>% unlist(use.names = FALSE)
neuro_mid <- setdiff(setdiff(expr_df$Sample, neuro_low), neuro_high)

neuro_low_df <- expr_df %>% dplyr::filter(Sample %in% neuro_low) 
neuro_high_df <- expr_df %>% dplyr::filter(Sample %in% neuro_high) 
neuro_mid_df <- expr_df %>% dplyr::filter(Sample %in% neuro_mid) 

# Save
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Neuro Low")
openxlsx::writeData(wb, sheet = "Neuro Low", x = neuro_low_df, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName = "Neuro High")
openxlsx::writeData(wb, sheet = "Neuro High", x = neuro_high_df, rowNames = FALSE)
openxlsx::addWorksheet(wb, sheetName = "Neuro Mid")
openxlsx::writeData(wb, sheet = "Neuro Mid", x = neuro_mid_df, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(data_path, "CCLE Neuro_Score-CDH12.xlsx"), overwrite = TRUE)


# Heatmap of 6 genes + CDH12
# mat <- read_data %>%
#   dplyr::select(everything(), -contains(c("ENTREZ", "ENSEMBL"))) %>%
#   tibble::column_to_rownames("SYMBOL") %>%
#   t() %>%
#   data.frame() %>%
#   dplyr::select(all_of(plot_genes), CDH12)
# 
# pheatmap(mat, scale = "column")

row_ann <- neuro_low_df %>% dplyr::mutate(Group = "LOW") %>% dplyr::select(Sample, Group) %>%
  dplyr::bind_rows(neuro_high_df %>% dplyr::mutate(Group = "HIGH") %>% dplyr::select(Sample, Group)) %>%
  dplyr::bind_rows(neuro_mid_df %>% dplyr::mutate(Group = "MID") %>% dplyr::select(Sample, Group)) %>%
  dplyr::mutate(Dummy = 0) %>%
  tibble::column_to_rownames("Sample")
row_elements <- unique(row_ann$Group)
  
mat <- normalized_counts %>%
  t() %>%
  data.frame() %>%
  dplyr::select(all_of(plot_genes), CDH12)

# Cluster within each group
row_order <- c()
for (g in row_elements){
  temp_mat <- mat[rownames(row_ann)[which(row_ann$Group == g)],]
  rowclust <- hclust(dist(temp_mat))
  row_order <- c(row_order, rownames(temp_mat[rowclust$order,]))
}

mat <- mat[rownames(row_ann), ]
mat <- mat[row_order, ]

pheatmap(mat, 
         cluster_cols = TRUE, cluster_rows = FALSE, 
         annotation_row = row_ann %>% dplyr::select(Group),
         filename = paste0(data_path, "CCLE Heatmap.tiff"))


################## TCGA classification neuroendocrine vs non-neuroendocrine

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Mukta/"

#remotes::install_github("cit-bioinfo/consensusMIBC", build_vignettes = TRUE)
library(consensusMIBC)

# # (DO THIS ONLY ONCE) Generate DESeq2 normalized counts
# read_data <- read.xlsx(paste0(data_path, "Read_data_TCGA-BLCA.xlsx"))
# meta_data <- read.xlsx(paste0(data_path, "Meta_data_TCGA-BLCA.xlsx"))
# 
# # If there are duplicated genes, keep only data for highest expressing copy
# read_data <- read_data %>%
#   dplyr::select(ENSEMBL_ID, make.names(meta_data$Sample_ID)) %>%
#   dplyr::mutate(n = rowSums(.[,-1])) %>%
#   dplyr::group_by(ENSEMBL_ID) %>%
#   dplyr::slice_max(n) %>%
#   dplyr::ungroup() %>%
#   # Duplicated genes with 0 expression in all samples still remain, remove them
#   dplyr::distinct_at("ENSEMBL_ID", .keep_all = TRUE) %>%
#   dplyr::select(everything(), -n)
# colnames(read_data) <- gsub(pattern = "_.*", replacement = "", x=colnames(read_data))
# colnames(read_data) <- make.names(colnames(read_data), unique = TRUE)
# colnames(read_data)[1] <- "SYMBOL"
# 
# # Get annotations
# species <- "Homo sapiens"
# annotations <- get_annotations(species)
# 
# Comparisons <- list(Variable =c(NA),
#                     Target   =c(NA),
#                     Reference=c(NA))
# meta_data <- prep_metadata(meta_data, read_data)
# read_data <- prep_readdata(read_data, meta_data)
# l <- check_data(read_data, meta_data)
# meta_data <- l[[2]]
# read_data <- l[[1]]
# 
# # Normalize the raw read counts
# dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
#                                       colData=meta_data,
#                                       design=~ 1)
# 
# approach <- ""
# suffix <- ""
# deseq2_norm_counts(dds, annotations, approach, suffix, data_path)

# CANNOT USE IT on CCLE as it has several non-bladder cell lines 
normalized_counts <- read.xlsx(paste0(data_path, "Normalized_Counts_TCGA-BLCA.xlsx"))
#normalized_counts <- read.xlsx(paste0(data_path, "CCLE_Normalized_Counts.xlsx"))

normalized_counts <- normalized_counts[, -c(1,3,4,5)] %>%
  dplyr::mutate(ENSEMBL_ID = make.names(names = ENSEMBL_ID, unique = TRUE)) %>%
  dplyr::distinct_at("ENSEMBL_ID", .keep_all = TRUE) %>%
  tibble::column_to_rownames(var = "ENSEMBL_ID") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
colnames(normalized_counts) <- base::make.names(names = colnames(normalized_counts))
normalized_counts <- normalized_counts[!rowSums(normalized_counts, na.rm=TRUE) == 0,]
normalized_counts_log <- log(1+normalized_counts, base=2)

res <- getConsensusClass(x=normalized_counts_log, gene_id = "ensembl_gene_id")

consensus_df <- res %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::left_join(normalized_counts_log %>%
                                 t() %>%
                                 data.frame() %>%
                                 tibble::rownames_to_column("Sample") %>%
                                 dplyr::select(Sample, ENSG00000154162),
                               by=c("Sample"="Sample")) %>%
  dplyr::left_join(normalized_counts %>%
                     t() %>%
                     data.frame() %>%
                     tibble::rownames_to_column("Sample") %>%
                     dplyr::select(Sample, ENSG00000154162),
                   by=c("Sample"="Sample")) %>%
  dplyr::rename(log2_CDH12 = identity(11), CDH12 = identity(12))

# Save
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "TCGA Consensus")
openxlsx::writeData(wb, sheet = "TCGA Consensus", x = consensus_df, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(data_path, "TCGA Consensus.xlsx"), overwrite = TRUE)
# openxlsx::addWorksheet(wb, sheetName = "CCLE Consensus")
# openxlsx::writeData(wb, sheet = "CCLE Consensus", x = consensus_df, rowNames = FALSE)
# openxlsx::saveWorkbook(wb, file = paste0(data_path, "CCLE Consensus.xlsx"), overwrite = TRUE)

######TCPA CDH12

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")
data_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Mukta/"

markers <- c("P53", "RB", "AKT", "AKT_pS473", "AKT_pT308")

protein_data <- read.xlsx(paste0(data_path, "TCGA-BLCA-L4.xlsx"))
protein_data <- protein_data[, -c(2:4)] %>%
  dplyr::mutate(Sample_ID = make.names(Sample_ID)) %>%
  dplyr::select(Sample_ID, all_of(markers))

read_data <- read.xlsx(paste0(data_path, "Normalized_Counts_TCGA-BLCA.xlsx"))
normalized_counts <- read_data[, -c(2:5)] %>%
  dplyr::mutate(SYMBOL = make.names(names = SYMBOL, unique = TRUE)) %>%
  dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
  tibble::column_to_rownames(var = "SYMBOL") %>%
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
colnames(normalized_counts) <- base::make.names(names = colnames(normalized_counts))
normalized_counts <- normalized_counts[!rowSums(normalized_counts, na.rm=TRUE) == 0,]

normalized_counts <- log(1+normalized_counts, base=2) # log transform the normalized counts so it is not skewed
t <- base::apply(X=normalized_counts, MARGIN=1, FUN=median, na.rm=TRUE)  # calculate median for each gene across all samples
normalized_counts <- base::sweep(x=normalized_counts, MARGIN=1, FUN="-", STATS=t) # subtract median of each gene
rna_data <- normalized_counts %>%
  t() %>%
  data.frame() %>%
  tibble::rownames_to_column("Sample_ID") %>%
  dplyr::select(Sample_ID, CDH12)

merge_data <- protein_data %>%
  dplyr::left_join(rna_data, by=c("Sample_ID"="Sample_ID")) %>%
  tibble::column_to_rownames("Sample_ID")

# Save
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "TCPA")
openxlsx::writeData(wb, sheet = "TCPA", x = merge_data, rowNames = TRUE)
openxlsx::saveWorkbook(wb, file = paste0(data_path, "TCPA.xlsx"), overwrite = TRUE)

breaks <- c(seq(from = -2, to = 0, length.out = 50), seq(from = 2/100, to = 2, length.out = 50))
pheatmap(merge_data, scale="row", breaks=breaks)

pheatmap(scale(merge_data), scale = "none", breaks=breaks)
