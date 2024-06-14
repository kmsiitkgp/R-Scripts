BiocManager::install("oligo")
BiocManager::install("oligoData")
BiocManager::install("illuminaHumanv4.db")
BiocManager::install("GEOquery")

#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

# Data analysis packages
library("ensembldb")
library("AnnotationHub")
library("GEOquery")
#library("affy")
#library("lumi")
library("oligo")
library("oligoData")
library("illuminaHumanv4.db")
#library("hgu133plus2.db")

# Data wrangling packages     
library("openxlsx")
library("dplyr")
library("tibble")
library("pheatmap")

#******************************************************************************#
#                           DECLARE GLOBAL VARIABLES                           #
#******************************************************************************#

# Save the gse of the project
gse <- "Blaveri"           # Data is already normalized, median centered. So, no analysis. 
# Blaveri data was obtained from Table S4 (10.1158/1078-0432.CCR-04-2409) and formatted. 
# UnigeneID was converted to gene_name using (http://resource.ibab.ac.in/GIDCON/geneid/home.html).
gse <- "GSE32894"
gse <- "GSE31684"
gse <- "GSE13507"
gse <- "mskcc"

parent_path <- "C:/Users/KailasammS/Desktop/"

# Choose species
species <- "Homo sapiens"
#species <- "Mus musculus"

# Choose platform and format of raw data
if (gse == "GSE32894"){
  array <- "Illumina"
  excel_files <- TRUE
  pvals <- TRUE
}
if (gse == "GSE31684"){
  array <- "Affy"
  excel_files <- FALSE
  pvals <- FALSE
}
if (gse == "GSE13507"){
  array <- "Illumina"
  excel_files <- TRUE
  pvals <- FALSE
}
if (gse == "mskcc"){
  array <- "Affy"
  excel_files <- TRUE
  pvals <- FALSE
}

#******************************************************************************#
#                         STEP 1: IMPORT CLINICAL DATA                         #
#******************************************************************************#

# Common abbreviations:
# AWD = alive with disease; DFS = diseasefree survival; 
# DOC = died of other cause; DOD = died of disease; 
# NED = no evidence of disease; NED II = no evidence of disease after recurrence;
# OS = overall survival; SCR = surgical complete remission. 

# Go to GEO (https://www.ncbi.nlm.nih.gov/geo/). 
# Search using GEO Accesssion (Eg: GSE32874)
# Once GSE32874 loads, scroll down and copy-paste sample info into 
# "sample_info.txt" file as per format below:

# GSM814052	UC_0001_1.a1.lbe1
# GSM814053	UC_0002_1.RNA.e1.a1.lbe1
# GSM814054	UC_0003_1.RNA.e1.a1.lbe1

# Name the file "sample_info.txt" and import it. Make sure the "Description"
# matches with "Sample_ID" in meta_data
sample_info <- utils::read.table(file = paste0(parent_path, gse, "_sample_info.txt"), 
                                 header = FALSE,
                                 sep = "\t") %>%
  dplyr::rename("GEO_ID" = V1, "Sample_id" = V2)

# Download patient data from the corresponding papers/GEO, save it as 
# "patient_data.xlsx" and import it. Make sure the there is a column which 
# matches with "Sample_id" in sample_info. 
# Also, some predefined columns name MUST be present in patient data
# "Sample_id" : this column MUST contain sample names
# "Time:      : this column MUST contain survival duration in days
# "Status     : this column MUST contain Dead/Alive status (Alive = 0, Dead = 1)
# "Sex        : this column MUST contain sex/Sex "Male" or "Female"
clinical_info <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, gse, "_clinical_info.xlsx"))

# Merge sample_info and meta_data
meta_data <-  dplyr::inner_join(x = sample_info, y = clinical_info, by = c("Sample_id"= "Sample_id")) %>%
  dplyr::select("GEO_ID", "Sample_id", "Sex", "Time", "Status", "Stage") %>%
  dplyr::distinct_at("Sample_id", .keep_all = TRUE)

#******************************************************************************#
#                          STEP 2: IMPORT COUNT DATA                           #
#******************************************************************************#

if (excel_files == "TRUE"){
  # Download background subtracted data and import it. Save it as xlsx file
  raw_data <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, gse, "_raw_data.xlsx"))
  colnames(raw_data)[1] <- "probe_id" 
  
  # (Optional) If p values are provided, ignore intensities where p val > 0.05
  if (pvals == "TRUE"){
    data <- raw_data[1]
    for (col in meta_data$Sample_id){
      temp <- raw_data %>%
        dplyr::select(probe_id, col, paste0(col,".detection.p.value")) %>%
        dplyr::filter(get(paste0(col,".detection.p.value")) < 0.05) %>%
        dplyr::select(probe_id, col)
      
      data <- dplyr::left_join(data, temp, by=c("probe_id"="probe_id"))
    }
    raw_data <- data
  }
  raw_data <- raw_data %>% 
    tibble::column_to_rownames(var="probe_id")
  
  # Rename the columns with appropriate GEO_ID if not already properly named
  for (i in 1:ncol(raw_data)){
    for (j in 1:nrow(meta_data)){
      if (colnames(raw_data)[i] == meta_data$Sample_id[j]){
        colnames(raw_data)[i] <- meta_data$GEO_ID[j]
      }
    }
  }
  
  # Remove probes which have no intensities in any sample
  raw_data <- raw_data[!(rowSums(raw_data, na.rm=TRUE) == 0),]
  
  # Log2 Transform the intensities if not already log2 trasnformed
  qx <- as.numeric(quantile(raw_data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE))
  LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { 
    raw_data <- raw_data %>% replace(. <=0, NA)
    raw_data_log <- log2(raw_data) 
  } else {
    raw_data_log <- raw_data
  }
  
  # Normalize data (assumes single-channel data)
  normalized_counts <- limma::normalizeBetweenArrays(object = as.matrix(raw_data_log), 
                                                     method = "quantile")
}else if (excel_files == "FALSE"){
  # NOTE: GEO cautions that values in series matrix files may not be comparable 
  # across samples. Moreover, the values in series matrix files are sometimes
  # normalized, sometimes normalized and median centered etc which complicates
  # the implementatio nof an universal code.
  # gset <- getGEO(gse)
  # gset <- gset[[1]]
  # raw_data <- exprs(gset) 
  
  # Extract raw data from CEL files
  cel_files <- oligoClasses::list.celfiles(paste0(parent_path, "Microarray"), full.names=TRUE)
  raw_data <- oligo::read.celfiles(cel_files)
  #raw_data <- ReadAffy(celfile.path = paste0(parent_path, gse, "_CEL"))
  
  # # (OPTIONAL) View one of the CEL files as image
  # oligo::image(x = raw_data[,1], transfo = rank)
  # oligo::image(x = raw_data[,1], transfo = log2)
  
  # (OPTIONAL) View MA plot
  # oligo::MAplot(object = raw_data[,1:3])
  # Getting phenotype data (raw_data@phenoData@data)
  #p_data <- pData(raw_data)
  
  # Perform RMA (background subtraction, quantile normalization, summarization)
  normalized_data <- oligo::rma(object = raw_data,
                                background = TRUE,
                                normalize = TRUE) 
  
  # View distribution before and after normalization
  oligo::hist(raw_data)
  oligo::hist(normalized_data)
  oligo::boxplot(x = raw_data, which=c("pm"), transfo=log2)
  oligo::boxplot(x = normalized_data, which=c("pm"), transfo=log2)
  
  # Extract normalized counts
  normalized_counts <- exprs(normalized_data)
  colnames(normalized_counts) <- gsub(pattern="\\..*", replacement="", x= colnames(normalized_counts))
}

# Filter out unwanted samples after normalization. If we perform this filtering
# before normalization, normalized values will be different
normalized_counts  <- as.data.frame(normalized_counts) %>% 
  dplyr::select(intersect(colnames(normalized_counts),meta_data$GEO_ID))

# Move rownames to 1st column
normalized_counts <- tibble::rownames_to_column(normalized_counts, "probe_id")

#******************************************************************************#
#                   STEP 4: GET PROBE ID-ENSEMBL ID MAPPINGS                   #
#******************************************************************************#

if (array == "Illumina"){
  # Get Ensembl IDs for probe IDs using "illuminaHumanv4.db" package
  xx <- as.data.frame(illuminaHumanv4ENSEMBL2PROBE)
} else if (array == "Affy"){
  # Get Ensembl IDs for probe IDs using "hgu133plus2.db" package
  x <- hgu133plus2.db::hgu133plus2ENSEMBL
  mapped_genes <- mappedkeys(x)
  xx <- as.data.frame(x[mapped_genes])
}

#******************************************************************************#
#                   STEP 5: MERGE GENE NAMES WITH PROBE IDS                    #
#******************************************************************************#

# NOTE: 7 probes may map to 1 gene in xx but only two of these probes might be 
# present in the expression data. So, we will remove probes lacking expressrion 
# and average the expression of the remaining probes so that each gene is 
# present ONLY ONCE in the data set 
normalized_counts <- xx %>% 
  dplyr::inner_join(normalized_counts, by=c("probe_id"="probe_id"))

# Get gene names for appropriate Ensembl IDs
normalized_counts <- annotations %>% 
  dplyr::inner_join(normalized_counts, by = c("ENSEMBL_ID"="ensembl_id"))

# Remove duplicated rows and probes with no gene names
normalized_counts <- normalized_counts %>%
  dplyr::filter(nchar(SYMBOL) > 0) %>%
  dplyr::select(everything(), -c(probe_id, ENSEMBL_ID, ENTREZ_ID, SYMBOL_ENTREZ, 
                                 START, END, CHR, STRAND, DESCRIPTION, BIOTYPE,
                                 BIOTYPE_ENTREZ)) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(across(.cols = everything(), .fns = ~mean(.x, na.rm=TRUE)))

# Save the normalized data
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = gse)
openxlsx::writeData(wb, sheet = gse, x = normalized_counts)
openxlsx::saveWorkbook(wb, file = paste0(parent_path, gse, "_Normalized.xlsx"), overwrite = TRUE)

# Save the meta data
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = gse)
openxlsx::writeData(wb, sheet = gse, x = meta_data)
openxlsx::saveWorkbook(wb, file = paste0(parent_path, gse, "_Metadata.xlsx"), overwrite = TRUE)

# #******************************************************************************#
# #                          STEP 6: VISUALIZE THE DATA                          #
# #******************************************************************************#
# 
# # Import genes of interest
# plot_genes <- openxlsx::read.xlsx(xlsxFile = paste0(parent_path, "hany paper/3B_venn_overlap_results.xlsx"))
# plot_genes <- unlist(unique(plot_genes[8]))
# plot_genes <- plot_genes[!is.na(plot_genes)]
# 
# #*******************************STEP 5A: HEATMAP*******************************# 
# 
# # Plot heatmap of normalized counts based on gender for genes of interest
# # Since some genes have multiple probes, we will average their expression. 
# # IT IS NOT RECOMMENDED TO AVERAGE BEFORE MAKING SURE THE PROBES TARGET SAME TRANSCRIPT
# # However, for sake of easiness, we are averaging here blindly.
# 
# # Keep only MIBC patients indicated by invasiveness = 2
# normalized_counts <- normalized_counts %>% 
#   dplyr::filter(gene_name %in% plot_genes) %>%
#   dplyr::group_by(gene_name) %>%
#   dplyr::summarize_at(vars(meta_data$case_submitter_id[which(meta_data$invasiveness>1)]), list(mean)) %>%
#   as.data.frame()
# 
# rownames(normalized_counts) <- normalized_counts[,1]
# normalized_counts <- normalized_counts[,-1]
# 
# # Create a dataframe to annotate samples of heatmap
# # The rownames of this data frame must match the colnames of the input matrix 
# sample_annotation <- data.frame(sex = meta_data$gender[which(meta_data$invasiveness>1)])
# if (identical(colnames(normalized_counts), meta_data$case_submitter_id[which(meta_data$invasiveness>1)])){
#   rownames(sample_annotation) <- colnames(normalized_counts)
# }
# 
# # Perform scaling of each gene across samples
# scaled_normalized_counts <- t(scale(t(normalized_counts)))
# 
# # If you try to plot heatmap with just 1 gene, you will get error
# pheatmap(mat = scaled_normalized_counts,
#          color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100),
#          kmeans_k = NA,
#          breaks = NA, 
#          border_color = "grey60",
#          cellwidth = NA, 
#          cellheight = NA, 
#          scale = "none", # since we are comparing same gene across samples,
#          # we will do row wise scaling
#          cluster_rows = TRUE,  #cluster the rows
#          cluster_cols = TRUE,  #cluster the columns
#          clustering_distance_rows = "euclidean", #"correlation",
#          clustering_distance_cols = "euclidean", #"correlation",
#          clustering_method = "complete", #average",
#          #clustering_callback = identity2,
#          cutree_rows = NA, 
#          cutree_cols = NA,
#          #treeheight_row = ifelse((class(cluster_rows) == "hclust") || cluster_rows,50, 0),
#          #treeheight_col = ifelse((class(cluster_cols) == "hclust") || cluster_cols, 50, 0),
#          legend = TRUE, 
#          legend_breaks = NA,
#          legend_labels = NA, 
#          #annotation_row = gene_anno,
#          annotation_col = sample_annotation,
#          annotation = NA,   # deprecated parameter
#          annotation_colors = NA, 
#          annotation_legend = TRUE,
#          annotation_names_row = FALSE,
#          annotation_names_col = FALSE,
#          drop_levels = TRUE, 
#          show_rownames = if_else(nrow(scaled_normalized_counts)<60, TRUE, FALSE, missing = NULL),
#          show_colnames = TRUE,
#          main = NA,
#          fontsize = 8, 
#          fontsize_row = 8, 
#          fontsize_col = 8,
#          angle_col = c("270", "0", "45", "90", "315"),
#          display_numbers = FALSE,
#          number_format = "%.2f", 
#          number_color = "grey30",
#          fontsize_number = 0.8*fontsize, 
#          gaps_row = NULL, 
#          gaps_col = NULL, 
#          labels_row = NULL,
#          labels_col = NULL,
#          filename = paste0(parent_path, "DEGs_Heatmap.pdf"),
#          width = NA, 
#          height = NA,
#          silent = FALSE,
#          na_col = "#DDDDDD")
# 
# #*******************************STEP 5B: BOXPLOT*******************************#
# 
# for (gene in plot_genes){
#   
#   if(gene %in% normalized_counts$gene_name){
#     plot_data <- normalized_counts %>%
#       dplyr::filter(gene_name == gene) %>%
#       dplyr::group_by(gene_name) %>%
#       dplyr::summarize_at(vars(meta_data$case_submitter_id), list(mean)) %>%
#       t() %>%
#       as.data.frame() %>% 
#       tibble::rownames_to_column("case_submitter_id")
#     
#     colnames(plot_data)[2] <- gene
#     
#     # Convert to numeric. Use :=  VERY IMPORTANT
#     plot_data <- plot_data %>% dplyr::mutate(!!rlang::sym(gene) := as.numeric(!!rlang::sym(gene)))
#     
#     # Keep only MIBC patients indicated by invasiveness = 2
#     plot_data <- dplyr::inner_join(meta_data, plot_data, by=c("case_submitter_id"="case_submitter_id")) %>%
#       dplyr::filter(invasiveness > 1)
#     
#     # Assumption 1: Are the two samples independents?
#     # Yes, since the samples from men and women are not related.
#     
#     # Assumtion 2: Are the data from each of the 2 groups follow a normal distribution?
#     # Null hypothesis ALWAYS assumes the opposite of what we want. 
#     # We hope that the data follows normal distribution. So,
#     # Null hypothesis: The distribution of the data is different from normal distribution
#     # If p < 0.05, then we accept null hypothesis =>  data DOESNT follow normal distribution
#     # If p > 0.05, we reject null hypothesis => data DOES follow normal distribution
#     
#     # Shapiro-Wilk normality test for Men
#     res.shapiro_1 <- stats::shapiro.test(x = as.list(plot_data %>% dplyr::filter(gender == "male") %>% dplyr::select(gene))[[1]])
#     # Shapiro-Wilk normality test for Women
#     res.shapiro_2 <- stats::shapiro.test(x = as.list(plot_data %>% dplyr::filter(gender == "female") %>% dplyr::select(gene))[[1]])
#     
#     # Assumption 3: Do the two populations have the same variances?
#     # IT IS ALWAYS RECOMMENDED TO USE WELCH'S T TEST
#     # NULL hypothesis: Two groups have different variance
#     # If p < 0.05, then we accept null hypothesis =>  2 groups have DIFFERENT variance. We MUST use WELCH's T TEST
#     # If p > 0.05, we reject null hypothesis => 2 groups have SAME variance & we can use STUDENT'S T TEST
#     
#     # Student's t-test
#     # res <- t.test(get(gene) ~ gender, data = plot_data, var.equal = TRUE)
#     
#     # Welch's t-test
#     res <- t.test(get(gene) ~ gender, data = plot_data, var.equal = FALSE)
#     
#     if (res$p.value < 0.05){
#       
#       cat("\n", gene, "\t\t","p-value: ", res$p.value, "\t\t", 
#           names(res$estimate[1]), " ", res$estimate[[1]], "\t\t", 
#           names(res$estimate[2]), " ", res$estimate[[2]])
#       
#       if (res$estimate[[1]] > res$estimate[[2]]){
#         if (stringr::str_detect(names(res$estimate[1]), "group male")){
#           signature_genes$male_high <- c(signature_genes$male_high, gene)
#         } else{
#           signature_genes$female_high <- c(signature_genes$female_high, gene)
#         }
#       } else if (res$estimate[[2]] > res$estimate[[1]]){
#         if (stringr::str_detect(names(res$estimate[2]), "group male")){
#           signature_genes$male_high <- c(signature_genes$male_high, gene)
#         } else{
#           signature_genes$female_high <- c(signature_genes$female_high, gene)
#         }
#       } else {}
#       
#       # Plot violin plot
#       lower_limit <- floor(min(plot_data[,gene]))
#       upper_limit <- ceiling(max(plot_data[,gene]))
#       
#       # fill=<> has to be categorical. If it is numerical like 1 and 2, it wort work
#       ggplot(data = plot_data, aes(x=gender, y=get(gene), fill=gender)) +
#         geom_violin(trim=FALSE) +         #display violin plot
#         theme_classic() +       
#         labs(x = "Gender", y = "Normalized Counts", title = gene) +
#         coord_cartesian(ylim = c(lower_limit, upper_limit)) +
#         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10), 
#               axis.text.y = element_text(size = 10),                                 
#               axis.title = element_text(size = 14),                                  
#               plot.title = element_text(hjust=0.5, size = 16, face="bold")) +  
#         geom_boxplot(width=0.1)
#       
#       # Save the plot
#       ggplot2::ggsave(filename = paste0("Box Plot for gene ", gene, ".tiff"),
#                       plot = last_plot(),
#                       device = "jpeg",
#                       path = paste0(parent_path, "microarray/"),
#                       scale = 1,
#                       width = 8,
#                       height = 8,
#                       units = c("in"),
#                       dpi = 600,
#                       limitsize = TRUE,
#                       bg = NULL)
#     }
#   }
# }
# 
# #**********************STEP 5C: SAVE THE SIGNATURE GENES***********************#
# 
# # Calculate maximum length to use for the dataframe
# j <- 1
# for (i in 2:length(signature_genes)){
#   if (length(signature_genes[[i]]) > length(signature_genes[[j]])){
#     nrow <- length(signature_genes[[i]])
#     j <- i
#     print(nrow)
#   } else{
#     nrow <- length(signature_genes[[j]])
#   }
# }
# 
# signature_dataframe <- data.frame(matrix(NA, nrow=nrow, ncol=length(signature_genes)))
# for (i in 1:length(signature_genes)){
#   colnames(signature_dataframe)[i] <- names(signature_genes)[i]
#   for (j in 1:length(signature_genes[[i]])){
#     signature_dataframe[j, i] <- signature_genes[[i]][j]
#     
#   } 
# }
# 
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb, sheetName = "signature_genes")
# openxlsx::writeData(wb, sheet = "signature_genes", x = signature_dataframe, rowNames = FALSE)
# 
# openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Signature_", gse, ".xlsx"),
#                        overwrite = TRUE,  returnValue = FALSE)
# 
