#!/usr/bin/env Rscript

#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

# Data analysis packages
library("ensembldb")            # Needed for annotating genes
library("AnnotationHub")        # Needed for annotating genes
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("DESeq2")               # Needed for Differential Expression analysis
library("sva")
library("preprocessCore")
#library("MAGeCKFlute")

# Data wrangling packages
library("openxlsx")             # Needed for reading, writing xlsx files
library("dplyr")                # Needed for data wrangling
library("tibble")               # Needed for data wrangling
library("stringr")              # Needed for data wrangling
library("purrr")                # Needed for data wrangling

# Graph plotting packages
library("ggplot2")              # Needed for making graphs
library("cowplot")              # Needed for merging multiple graphs
library("viridis")              # Needed for nice graph coloring
library("RColorBrewer")         # Needed for nice graph coloring
library("ggrepel")              # Needed for making graphs prettier
library("ggbeeswarm")           # Needed for proper positioning of labels in scatter plots
library("colorspace")
library("UpSetR")

# Specialized Graph plotting packages
library("pheatmap")             # Needed for making heatmap
library("ggridges")             # Needed for making ridgeplots
library("VennDiagram")          # Needed for making Venn diagram
library("survival")             # Needed for making survival curves
library("survminer")            # Needed for making survival curves, to handle %++% in survival function

# Define axis font etc to use in all plots
my_theme <- ggplot2::theme(plot.title=  element_text(family="sans", face="bold",  colour="black", size=15, hjust=0.5),
                           legend.title=element_text(family="sans", face="bold",  colour="black", size=12, hjust=0,   vjust=1,   angle=0),
                           axis.title.x=element_text(family="sans", face="bold",  colour="black", size=12, hjust=0.5, vjust=0,   angle=0),
                           axis.title.y=element_text(family="sans", face="bold",  colour="black", size=12, hjust=0.5, vjust=1,   angle=90),
                           legend.text= element_text(family="sans", face="plain", colour="black", size=10, hjust=0.5),
                           axis.text.x= element_text(family="sans", face="plain", colour="black", size=10, hjust=0.5, vjust=0.5, angle=45),
                           axis.text.y= element_text(family="sans", face="plain", colour="black", size=10, hjust=0.5, vjust=0.5, angle=0))
# strip.text.x=element_text(family="sans", face="bold",  colour="black", size=10, hjust=0.5),
# legend.background=element_rect(fill="lightblue", size=0.5, linetype="solid", colour ="darkblue"),
# legend.position="right",
# legend.justification="left",
# legend.direction="vertical",
# legend.key.height=unit(0.5, 'cm'),
# legend.key.width =unit(0.5, 'cm'), 
# legend.text.align=0)


#******************************************************************************#
#                               GENE ANNOTATIONS                               #
#******************************************************************************#

# This function returns a dataframe with following columns:
# "ENSEMBL_ID", "ENTREZ_ID", "SYMBOL", "SYMBOL_ENTREZ", "BIOTYPE", 
# "BIOTYPE_ENTREZ", "START", "END", "CHR", "STRAND", "DESCRIPTION"

get_annotations <- function(species){
  
  # AnnotationHub has SYMBOL-ENSEMBL_ID info ONLY.
  # AnnotationDbi has SYMBOL-ENSEMBL_ID as well as SYMBOL-ENTREZ_ID info.
  # hubCache(AnnotationHub()) to find location where cache is stored and delete
  # it and strart fresh if you get errors like "Error: failed to load resource"
  
  #setAnnotationHubOption("CACHE", "C:/Users/kailasamms/AppData/Local/R/win-library/4.4/AnnotationHub")
  
  #**************************GET ENSEMBL ANNOTATIONS***************************#
  # Connect to AnnotationHub
  ah <- AnnotationHub::AnnotationHub()
  
  # Access the Ensembl database for organism
  ahDb <- AnnotationHub::query(x=ah,
                               pattern=c(species, "EnsDb"),
                               ignore.case=TRUE)
  
  # Acquire the latest annotation files
  id <- ahDb %>%
    mcols() %>%
    rownames() %>%
    tail(n=1)
  
  # Download the appropriate Ensembldb database
  edb <- ah[[id]]
  
  # Extract gene-level information from database
  ensembl <- ensembldb::genes(x=edb,
                              return.type="data.frame")
  
  # Select annotations of interest
  ensembl <- ensembl %>%
    dplyr::rename(ENSEMBL_ID=gene_id, SYMBOL=gene_name, 
                  BIOTYPE=gene_biotype, START=gene_seq_start, END=gene_seq_end, 
                  CHR=seq_name, STRAND=seq_strand, DESCRIPTION=description,
                  ENSEMBL_TRANSCRIPT = canonical_transcript) %>%
    dplyr::mutate(SYMBOL=dplyr::case_when(nchar(SYMBOL) == 0 ~ NA,
                                          TRUE ~ SYMBOL)) %>%
    dplyr::select(ENSEMBL_ID, ENSEMBL_TRANSCRIPT, SYMBOL, BIOTYPE, START, END, CHR, STRAND, DESCRIPTION)
  
  #***************************GET ENTREZ ANNOTATIONS***************************# 
  
  # NOTE: mapIds can ONLY retrieve one of "EMSEMBL/SYMBOL/GENETYPE" at a time
  # mapping <- AnnotationDbi::mapIds(x=org.Hs.eg.db, 
  #                                  keys=keys(org.Hs.eg.db),
  #                                  keytype="ENTREZID", 
  #                                  column="SYMBOL") %>%
  #   as.data.frame(do.call(cbind, list(.))) %>%
  #   tibble::rownames_to_column("ENTREZID") %>%
  #   dplyr::rename(ENTREZID=identity(1), SYMBOL=identity(2))
  
  if (species == "Homo sapiens"){
    entrez <- AnnotationDbi::select(x=org.Hs.eg.db, keys=keys(org.Hs.eg.db),
                                    columns=c("ENSEMBL", "SYMBOL","GENETYPE"))
  } else if (species == "Mus musculus"){
    entrez <- AnnotationDbi::select(x=org.Mm.eg.db, keys=keys(org.Mm.eg.db),
                                    columns=c("ENSEMBL", "SYMBOL","GENETYPE"))
  }
  
  colnames(entrez) <- c("ENTREZ_ID", "ENSEMBL_ID", "SYMBOL_ENTREZ", "BIOTYPE_ENTREZ")
  
  # Merge ensembl and entrez dataframes
  annotations <- dplyr::full_join(ensembl, entrez, by=c("ENSEMBL_ID"="ENSEMBL_ID")) %>%
    dplyr::select(ENSEMBL_ID, ENSEMBL_TRANSCRIPT, ENTREZ_ID, SYMBOL, 
                  SYMBOL_ENTREZ, BIOTYPE, BIOTYPE_ENTREZ, START, END, CHR, 
                  STRAND, DESCRIPTION)
  
  #DBI::dbDisconnect(conn=ah)
  return(annotations)
}

#******************************************************************************#
#                COMPILE RAW COUNTS FROM STAR/HTSEQ COUNT FILES                #                                           
#******************************************************************************#

# Read individual count files.txt and merge all counts into "Read_data.xlsx"
compile_raw_counts <- function(parent_path, method){
  
  # Create a list of all txt files within folder that will be analyzed
  count_folder <- paste0(parent_path, "count_results/")
  files <- list.files(path=count_folder)
  
  # Create an empty dataframe with 0's
  read_data <- data.frame(0)
  
  # Create the reads table 
  for (i in 1:length(files)){
    
    # Read the txt file
    temp_file <- read.table(file=paste0(count_folder, files[i]), header=FALSE, sep="\t")  
    
    # Remove top 4 rows in STAR count output:  
    # N_unmapped, N_multimapping, N_noFeature, N_ambiguous 
    if (method == "STAR"){
      temp_file <- temp_file[5:nrow(temp_file),]
    }
    
    # Remove last 5 rows in HTSeq count output: __no_feature,  __ambiguous, 
    # __too_low_aQual, __not_aligned, __alignment_not_unique
    if (method == "HTSEQ"){
      temp_file <- temp_file[1:(nrow(temp_file)-5),]
    }
    
    # The 1st Column will have Ensembl ids. 
    # For HTSeq, gene counts may be in 2nd or 3rd column. 
    # For STAR, gene counts are in 2nd (unstranded), 3rd (+), 4th (-) column. 
    # Append appropriate column.
    if (method == "HTSEQ" & sum(temp_file[2], na.rm=TRUE) == 0 & sum(temp_file[3], na.rm=TRUE) > 0){
      read_data <- bind_cols(read_data, temp_file[,3])
    } else if (method == "HTSEQ" &  sum(temp_file[2], na.rm=TRUE) > 0 & sum(temp_file[3], na.rm=TRUE) ==0){
      read_data <- bind_cols(read_data, temp_file[,2])                  
    } else if (method == "STAR" & abs((sum(temp_file[2])/sum(temp_file[3])) - (sum(temp_file[2])/sum(temp_file[4]))) < 2){
      print("Unstranded")
      read_data <- bind_cols(read_data, temp_file[,2])                  
    } else if (method == "STAR" & sum(temp_file[3]) > 3*sum(temp_file[4])){
      print("Pos stranded")
      read_data <- bind_cols(read_data, temp_file[,3])                  
    } else if (method == "STAR" & sum(temp_file[4]) > 3*sum(temp_file[3])){
      print("Neg stranded")
      read_data <- bind_cols(read_data, temp_file[,4])                  
    } else{
      print("Error: Gene counts NOT PRESENT in either column 2 or 3 of count file")
    }
    
    # Rename the column names to sample names
    colnames(read_data)[i+1] <- gsub(pattern="\\..*$", replacement="", x=files[i])
  }
  
  # Check if all count files have same order of genes in the rows so that files can be merged together
  temp_file <- read.table(file=paste0(count_folder, files[1]), header=FALSE, sep="\t")
  gene_list <- temp_file[, 1]
  for (i in 1:length(files)){
    
    temp_file <- read.table(file=paste0(count_folder, files[i]), header=FALSE, sep="\t")
    genes <- temp_file[, 1]
    
    if (!identical(gene_list, genes)){ 
      print("Gene order is different between the count files")
    }
  }
  
  if (method == "STAR"){
    gene_list <- gene_list[5:length(gene_list)]
  } else if(method == "HTSEQ"){
    gene_list <-gene_list[1:(length(gene_list)-5)]
  }
  
  # Add gene names to 1st column
  read_data[, 1] <- gene_list
  colnames(read_data)[1] <- "SYMBOL"
  colnames(read_data) <- gsub(pattern="ReadsPerGene", replacement="", x=colnames(read_data))
  
  # Remove genes with 0 counts in all samples
  read_data <- read_data[rowSums(read_data[,-1]) != 0,]
  
  # Remove samples with 0 counts in all genes
  read_data <- read_data[,colSums(read_data[,-1]) != 0]
  
  # Save the results as xlsx file
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Raw_counts")
  openxlsx::writeData(wb, sheet="Raw_counts", x=read_data)
  openxlsx::saveWorkbook(wb, file=paste0(parent_path, "Read_data.xlsx"),
                         overwrite=TRUE)
  
  return(read_data)
}

#****************************************************************************#
#                           STEP 1: REFORMAT META DATA                         #
#****************************************************************************#

# Remove rows with missing sample names
# Add sample names as rownames
# Create a new column "id" which has groups that need to be compared
prep_metadata <- function(meta_data, Variable){  
  
  meta_data <- meta_data %>%
    dplyr::filter(!is.na(Sample_ID)) %>%
    dplyr::mutate(Sample_ID = make.names(Sample_ID, unique=TRUE)) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames(var="Sample_ID") #%>%
  #  dplyr::mutate(id=get(Variable))
  # dplyr::mutate(id=dplyr::if_else(get(Variable) %in% c(Comparisons$Target, Comparisons$Reference),
  #                                   get(Variable), "ignore")) %>%
  # dplyr::mutate(id=dplyr::if_else(is.na(id), "ignore", id))
  
  return(meta_data)
}  

#****************************************************************************#
#                          STEP 2: REFORMAT READ DATA                          #                                           
#****************************************************************************#

# Remove rows with missing gene names
# Add gene names as rownames
# Keep ONLY samples present in metadata
# Replace all blank values i.e. counts with 0
prep_readdata <- function(read_data, meta_data){
  
  read_data <- read_data %>%
    dplyr::filter(!is.na(SYMBOL)) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames(var="SYMBOL") %>%
    dplyr::select(intersect(colnames(.), rownames(meta_data))) %>%
    replace(is.na(.), 0)
  
  return (read_data)
}  

#****************************************************************************#
#              STEP 3: PREPARE META DATA and READ DATA FOR DESEQ2            #
#****************************************************************************#

# Arrange columns of read_data to match meta_data
# Remove samples that have NA for the variable being compared for DEG
# Remove any column names "sizeFactor" in meta_data if present
# Check if all columns in metadata are factors. If not, convert to factors
# Check if read_data and meta_Data are dataframes. If not, convert to dataframe
# Check if column names of read_data are present in same order as rownames of meta_data
# Check if ALL column names of read_data are present in rownames of metadata
check_data <- function(read_data, meta_data){  
  
  # DESeq2 automatically removes genes that have 0 counts but doesnt remove 
  # samples that have 0 counts for all genes (unlikely scenario). So, remove
  # such samples, else, the geometric mean will be 0 for all genes and DESeq2 
  # will halt execution.
  
  if (TRUE %in% (colSums(read_data)==0)){
    read_data <- read_data[, -c(which(colSums(read_data)==0))]
    meta_data <- meta_data[-c(which(colSums(read_data)==0)),]
  }
  
  # Rearrange read_data to match meta_data
  read_data <- read_data[,rownames(meta_data)]
  
  # Remove samples that have "NA" for variable being analyzed in metadata
  #NA_sample_list <- which(is.na(meta_data[Variable]))
  NA_sample_list <- c()
  for (i in 1:length(Comparisons$Variable)){
    NA_sample_list <- c(NA_sample_list, which(is.na(meta_data[Comparisons$Variable[i]])))
  }
  NA_sample_list <- unique(NA_sample_list)

  if (length(NA_sample_list)!=0){
    meta_data <- meta_data[-c(NA_sample_list),]
    read_data <- read_data[,-c(NA_sample_list)]
  }
  
  # Remove any column named "sizeFactor" in meta_data as it interferes with DESeq2()
  if ("sizeFactor" %in% colnames(meta_data)){
    meta_data <- meta_data %>% dplyr::select(!sizeFactor)
  }
  
  # Check if all columns of meta table are factors. Else, convert to factors.
  print("meta_data before factor conversion")
  str(meta_data)
  for (i in 1:ncol(meta_data)){
    meta_data[,i] <- as.factor(meta_data[,i])
  }
  print("meta_data after factor conversion")
  str(meta_data)
  
  # Check (i) read_data & meta_data are dataframes
  # (ii) column names of read_data (i.e.sample name) are PRESENT in row names of meta_data (i.e. sample name)
  # (iii) column names of read_data (i.e.sample name) are PRESENT in same order as row names of meta_data (i.e. sample name)
  if (class(read_data) == "data.frame" & class(meta_data) == "data.frame" & 
      base::all(colnames(read_data) %in% rownames(meta_data)) &
      base::all(colnames(read_data) == rownames(meta_data))){
    print("ALL OK")
  } else if (class(read_data) != "data.frame"){
    print("read_data is NOT a dataframe")
  } else if (class(meta_data) != "data.frame"){
    print("meta_data is NOT a dataframe")
  } else if (!(all(colnames(read_data) %in% rownames(meta_data)))){
    print("Missing samples")
  } else if (all(colnames(read_data) != rownames(meta_data))){
    print("Varying sample order")
  } else{}
  
  return(list(read_data, meta_data)) 
}

#****************************************************************************#
#                      STEP 4: PERFORM BATCH CORRECTION                      #
#****************************************************************************#

# NOTE: It is RECOMMENDED to perform batch correction ONLY if you know the batch
# information for all the samples in meta_data.

# https://support.bioconductor.org/p/76099/
# https://support.bioconductor.org/p/9149116/
# https://support.bioconductor.org/p/133222/#133225/
# https://support.bioconductor.org/p/125386/#125387
# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#using-sva-with-deseq2

# There are multiple approaches to batch correction:

# (i) Modelling batch effect using DESeq2 design [RECOMMENDED]
# Add Batch to design of DESEq2 like design=~Batch + id

# (ii) Using combatseq to remove batch effects
# Use combatseq to get batch corrected read counts and then use the corrected 
# counts as input for DESeq2() and a design without batch variable (design=~ id)

# NOTE: I have tested (i) and (ii) and found that 
# (a) Almost all DEGs identified as significant in (i) are present in (ii) 
# (b) All significant DEGs with -ve log2FC from (i) also have -ve log2FC in (ii)
# (c) All significant DEGs with +ve log2FC from (i) also have +ve log2FC in (ii)
# (d) (ii) identifies many other significant DEGs missing in (i)
# (e) log2FC from (ii) differs from (i) for all genes but difference is minimal
# for most of the significant DEGs

# (iii) Using sva to identify hidden batch variables and model them in DeSeq2
# Use sva to find upto 2 surrogate variables that could be causing batch effects
# and include them in the design of DESeq2 like design=~SV1 + SV2 + id

# NOTE: I have tested (i) and (iii) and found that
# (a) Majority of significant DEGs in (iii) match with (i) but there are many
# significant DEGs in (i) absent in (iii) and vice versa
# (b) Significant DEGs with -ve log2FC from (i) also have -ve log2FC in (iii)
# (c) Significant DEGs with +ve log2FC from (i) also have +ve log2FC in (iii) 
# (d) (iii) identifies many other significant DEGs missing in (i) but (iii)
# fails to identify many DEGs identified by (i)
# (e) log2FC from (iii) differs from (i) for MOST genes to a great extent,
# HOWEVER, DDR2KO samples had log2FC(DDR2) of -0.7 in (iii) but only -0.37 in 
# (i) and (ii). So, sva might infact be removing batch effects and detecting 
# true biological effects !!! 

# combatseq batch correction:
# NOTE: This batch correction of known factors is done on raw counts
combatseq_batch <- function(read_data, meta_data){  
  
  if (length(unique(as.numeric(meta_data$Batch))) > 1){
    
    corrected_read_data <- sva::ComBat_seq(counts=as.matrix(read_data), 
                                           batch=as.numeric(meta_data$Batch), 
                                           group=as.numeric(as.factor(meta_data$id)),
                                           full_mod = TRUE)
  } else{
    corrected_read_data <- read_data
  }
  
  return(corrected_read_data) 
}

# svaseq batch correction:
# NOTE: This batch correction of unknown factors is done on normalized counts
# NOTE: svaseq() can find n number of surrogate variables. If we model for all 
# of them there could be over correction. Hence, we limit batch correction to
# only the top 2 surrogate variables.
# NOTE: Lowly expressed genes are removed before finding surrogate variables
# in "rowMeans(dat) > 1".
# Here, we just create a new object sva_dds with sva design variables
svaseq_batch <- function(read_data, meta_data){
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~1)
  dds <- DESeq2::estimateSizeFactors(dds) 
  dat  <- DESeq2::counts(dds, normalized = TRUE)
  idx  <- rowMeans(dat) > 1
  dat  <- dat[idx, ]
  mod  <- stats::model.matrix(~ id, colData(dds))
  mod0 <- stats::model.matrix(~ 1, colData(dds))
  svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0)
  
  if (svseq$n.sv > 2){
    svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0, n.sv=2)
  }
  ddssva <- dds
  
  # ddssva$SV1 <- svseq$sv[,1]
  # ddssva$SV2 <- svseq$sv[,2]
  # design(ddssva) <- ~ SV1 + SV2 + id
  
  for (i in 1:ncol(svseq$sv)){
    var <- paste0("SV",i)
    ddssva[[var]] <- svseq$sv[,i]
  }
  design(ddssva) <- as.formula(paste0("~", paste0("SV", seq(1:ncol(svseq$sv)), collapse = "+"), "+id"))
  
  return(ddssva)
}

#****************************************************************************#
#                     STEP 5: CALCULATE NORMALIZED COUNTS                    #
#****************************************************************************#

add_annotation <- function(normalized_counts, annotations){
  
  ensembl <- length(intersect(normalized_counts$ID, annotations$ENSEMBL_ID))
  entrez <- length(intersect(normalized_counts$ID, annotations$ENTREZ_ID))
  symbol <- length(intersect(normalized_counts$ID, annotations$SYMBOL))
  
  if (ensembl > entrez & ensembl > symbol){
    normalized_counts <- annotations %>%
      dplyr::right_join(normalized_counts, by=c("ENSEMBL_ID"="ID"), multiple="all") %>%
      dplyr::select(SYMBOL, everything(), -c(CHR, DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
      dplyr::mutate(SYMBOL=dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                            is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                            TRUE ~ SYMBOL)) %>%
      dplyr::distinct_at("ENSEMBL_ID", .keep_all=TRUE)
  }
  
  if (entrez > ensembl & entrez > symbol){
    normalized_counts <- annotations %>%
      dplyr::right_join(normalized_counts, by=c("ENTREZ_ID"="ID"), multiple="all") %>%
      dplyr::select(SYMBOL, everything(), -c(CHR, DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
      dplyr::mutate(SYMBOL=dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                            is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                            TRUE ~ SYMBOL)) %>%
      dplyr::distinct_at("ENTREZ_ID", .keep_all=TRUE)
  }
  
  if (symbol > ensembl & symbol > entrez){
    normalized_counts <- annotations %>%
      dplyr::right_join(normalized_counts, by=c("SYMBOL"="ID"), multiple="all") %>%
      dplyr::select(SYMBOL, everything(), -c(CHR, DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
      dplyr::mutate(SYMBOL=dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                            is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                            TRUE ~ SYMBOL)) %>%
      dplyr::distinct_at("SYMBOL", .keep_all=TRUE)
  }
  
  return(normalized_counts)
}

# Normalized counts are influenced by sizeFactors.
# sizeFactors are affected by number of samples (all samples vs subset of samples)
# sizeFactors are NOT affected by design formula. So, design ~1 or 
# design ~Batch+1 or design ~id will give the same normalized values.
# sizeFactors MUST be estimated first before normalization.

# Normalized counts from dds object are NOT batch corrected. We do this below.
# https://www.biostars.org/p/490181/
deseq2_norm_counts <- function(dds, annotations, approach, suffix){  
  
  dds <- DESeq2::estimateSizeFactors(dds) 
  
  normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
  
  if (length(unique(meta_data$Batch)) > 1){
    normalized_counts <- limma::removeBatchEffect(x=log2(normalized_counts+1), 
                                                  dds$Batch)
  }
  
  normalized_counts <- normalized_counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID")
  
  normalized_counts <- add_annotation(normalized_counts, annotations)
  
  # Save batch corrected normalized counts for entire dataset
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="normalized")
  openxlsx::writeData(wb, sheet="normalized", x=normalized_counts)
  openxlsx::saveWorkbook(wb, file=paste0(results_path, "Normalized_Counts_", approach, "_", suffix, ".xlsx"),
                         overwrite=TRUE)
}

combatseq_norm_counts <- function(dds, annotations, approach, suffix){
  
  dds <- DESeq2::estimateSizeFactors(dds) 
  
  normalized_counts <- DESeq2::counts(dds, normalized=TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID")
  
  normalized_counts <- add_annotation(normalized_counts, annotations) 
  
  # Save batch corrected normalized counts for entire dataset
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="normalized")
  openxlsx::writeData(wb, sheet="normalized", x=normalized_counts)
  openxlsx::saveWorkbook(wb, file=paste0(results_path, "Normalized_Counts_", approach, "_", suffix, ".xlsx"),
                         overwrite=TRUE)
  
}

# svaseq corrected normalized counts:
# NOTE: ddssva object from svaseq_batch has the top 2 surrogate variables that 
# will be used in DESeq2() but the normalized counts from ddssva object are NOT 
# batch corrected. We do this below.
# https://www.biostars.org/p/121489/

# NOTE: Lowly expressed genes are removed before finding surrogate variables
# in "rowMeans(dat) > 1". As a result, the number of genes in normalized counts
# excel file is lower than DeSeq2.
svaseq_norm_counts <- function(dds, annotations, approach, suffix){
  
  dat  <- DESeq2::counts(dds, normalized = TRUE)
  idx  <- rowMeans(dat) > 1
  dat  <- dat[idx, ]
  mod  <- stats::model.matrix(~ id, colData(dds))
  mod0 <- stats::model.matrix(~ 1, colData(dds))
  svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0)
  
  if (svseq$n.sv > 2){
    svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0, n.sv=2)
  }
  
  # %*% indicates Matrix multiplication
  X <- base::cbind(mod, svseq$sv)
  Hat <- base::solve(t(X) %*% X) %*% t(X)
  beta <- (Hat %*% t(dat))
  P <- ncol(mod)
  corrected_data <- dat - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])
  
  normalized_counts <- corrected_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID")
  
  normalized_counts <- add_annotation(normalized_counts, annotations) 
  
  # Save batch corrected normalized counts for entire dataset
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="normalized")
  openxlsx::writeData(wb, sheet="normalized", x=normalized_counts)
  openxlsx::saveWorkbook(wb, file=paste0(results_path, "Normalized_Counts_", approach, "_", suffix, ".xlsx"),
                         overwrite=TRUE)
}

#****************************************************************************#
#                             STEP 6: RUN DESEQ2                             #
#****************************************************************************#

# This function runs DESeq2 and plots volcano plot, heatmap, PCA plot and
# correlation plot.
run_deseq2 <- function(dds, meta_data, annotations, Comparisons, n, lfc.cutoff, padj.cutoff, approach, suffix){  
  
  # Run analysis using both local and parameteric fit
  # betaPrior: default=FALSE, shrunken LFCs are obtained later using lfcShrink
  # fitType: "local" may be better sometimes
  dds_para <- DESeq2::DESeq(object=dds, 
                            test="Wald",
                            fitType="parametric",
                            betaPrior=FALSE,
                            minReplicatesForReplace=7)
  dds_local <- DESeq2::DESeq(object=dds, 
                             test="Wald",
                             fitType="local",
                             betaPrior=FALSE,
                             minReplicatesForReplace=7)
  
  # Determine whether to use "local" or "parametric" fit
  residual_para <- mcols(dds_para)$dispGeneEst - mcols(dds_para)$dispFit
  residual_local <- mcols(dds_local)$dispGeneEst - mcols(dds_local)$dispFit
  
  if (median(residual_para^2, na.rm=TRUE) <= median(residual_local^2, na.rm=TRUE)){
    fit <- "parametric"
    dds <- dds_para
  } else{
    fit <- "local"
    dds <- dds_local
  }
  
  # Define contrast and coeff
  DE_levels <- as.vector(unique(meta_data$id))
  contrast <- c("id", 
                DE_levels[DE_levels==Comparisons$Target[n]],
                DE_levels[DE_levels==Comparisons$Reference[n]])
  
  
  coeff <- paste(contrast[1], contrast[2], "vs", contrast[3], sep='_')
  cat("Coeff is ", coeff, "\n")
  
  # Calculate results
  res <- DESeq2::results(object=dds,
                         contrast=contrast,
                         lfcThreshold=lfc.cutoff,
                         altHypothesis="greaterAbs",
                         cooksCutoff=TRUE,
                         independentFiltering=TRUE,
                         alpha=padj.cutoff,
                         pAdjustMethod="BH")
  
  # Perform lfcshrinkage to account for variability between replicates
  # For ashr, if res is provided, then coef and contrast are ignored.
  # lfcshrinkage will not change the number of DEGs and affects only logFC
  res <- DESeq2::lfcShrink(dds=dds, res=res, type="ashr")
  
  # Summarize results
  summary(res)
  
  # Export results as dataframe and replace ensembl IDs with Gene names
  results <- res %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") 
  
  results <- add_annotation(results, annotations)
  
  # Add unique identifier to each result file
  file_prefix <- gsub(pattern="/", replacement="", x=coeff)
  
  # Save the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="DEGs")
  openxlsx::writeData(wb, sheet="DEGs", x=results, rowNames=FALSE)
  openxlsx::saveWorkbook(wb, file=paste0(results_path, "Results_", file_prefix , "_", approach, "_", suffix, "_DEGs.xlsx"), 
                         overwrite=TRUE)
  
  return(dds)
  
  #**************************************************************************#
  #             CREATE HEATMAPS, VOLCANO PLOTS, PCA PLOTS            #
  #**************************************************************************#
  
  if (heatmap_plot == TRUE){
    
    # Refer https://www.biostars.org/p/339065/ where the authors of DESeq2 recommend
    # using vst or rld rather than log1p transformation while plotting DESeq2
    # normalized counts in heatmap.
    # NOTE: We scale the vst values and replace all NAs with 0 before feeding it to
    # pheatmap(). If we perform the scaling in pheatmap, the NAs generated during
    # scaling will give error.
    # NOTE: While use of log1p transformation is NOT recommended, I notice only
    # minor changes in color
    
    # Define any filename you want added to final file
    file_prefix <- paste0(Variable2_value, "_", coeff)
    file_prefix <- gsub(pattern="/", replacement="", x=file_prefix)
    file_prefix <- paste0(file_prefix, "_", f_suffix)
    
    # (I) Read metadata
    metadata <- meta_data %>%
      tibble::rownames_to_column("Sample")
    
    # (II) Define genes to plot
    plot_genes <- results %>%
      dplyr::filter(padj < 0.05) %>%
      dplyr::select(SYMBOL) %>%
      unlist(use.names=FALSE)
    
    # (III) Genes to display in heatmap
    disp_genes <- c()
    
    # (IV) Read expr data
    #normalized_counts <- normalized_counts
    normalized_counts <- assay(DESeq2::vst(dds, blind=FALSE)) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("ENSEMBL_ID") %>%
      dplyr::left_join(annotations, by=c("ENSEMBL_ID"="ENSEMBL_ID"), multiple="all") %>%
      dplyr::filter(nchar(SYMBOL) > 0) %>%
      dplyr::select(SYMBOL, everything(), -c(ID, CHR, DESCRIPTION, START, END, STRAND))
    colnames(normalized_counts)[1] <- "SYMBOL"
    
    # Run the function
    plot_heatmap(normalized_counts, metadata, plot_genes, disp_genes, file_prefix)
  }
  
  if (volcano_plot == TRUE){
    
    # Define any filename you want added to final file
    file_prefix <- paste0(Variable2_value, "_", coeff)
    file_prefix <- gsub(pattern="/", replacement="", x=file_prefix)
    file_prefix <- paste0(file_prefix, "_", f_suffix)
    
    # (I) Import expression data with log2FC and pval
    volcano_df <- results %>%
      dplyr::rename(log2FC=log2FoldChange, Gene=SYMBOL)
    
    # (II) Define or import metadata
    # NOTE: Metadata is a dataframe with "Sample" and "Condition" columns
    metadata <- meta_data %>%
      tibble::rownames_to_column("Sample")
    
    # (III) Define any genes you want to mark in volcano plot
    disp_genes <- c()
    
    # Make volcano plots
    plot_volcano(volcano_df, disp_genes, file_prefix)
  }
  
  if (cor_plot == TRUE){
    
    #Determine number of cuts in heatmap
    cuts <- 20
    
    # Use vst as well as rld
    vst <- DESeq2::vst(dds)
    sig_genes <- results %>%
      dplyr::filter(abs(log2FoldChange) >= 0.58 & padj < 0.05 & padj > 0 & !is.na(padj)) %>%
      dplyr::select(SYMBOL) %>%
      unlist(use.names=FALSE)
    vst <- vst[rownames(vst) %in% sig_genes]
    head(assay(vst))
    sampleDists <- dist(assay(vst))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- dplyr::left_join(data.frame("ENSEMBL_ID"=rownames(sampleDistMatrix)), annotations, by=c("ENSEMBL_ID"="ENSEMBL_ID")) %>%
      dplyr::select(SYMBOL) %>%
      unlist(use.names=FALSE)
    colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
    
    # cluster and re-order rows
    rowclust <- hclust(dist(sampleDistMatrix))
    
    # cluster and re-order columns
    colclust <- hclust(dist(t(sampleDistMatrix)))
    reordered <- sampleDistMatrix[rowclust$order,colclust$order]
    
    # Save batch corrected normalized counts for entire dataset
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName="Correlation")
    openxlsx::writeData(wb, sheet="Correlation", x=reordered , rowNames=TRUE)
    openxlsx::saveWorkbook(wb,
                           file=paste0(results_path, "Results_", Variable2_value, "_Correlation_Matrix.xlsx"),
                           overwrite=TRUE)
    
    pheatmap::pheatmap(mat=reordered,
                       scale="none",
                       cellwidth=3,
                       cellheight=2,
                       cutree_rows=cuts,
                       cutree_cols=cuts,
                       cluster_rows=TRUE,   #cluster the rows
                       cluster_cols=TRUE,   #cluster the columns
                       fontsize=8,
                       fontsize_row=8,
                       fontsize_col=8,
                       show_rownames=FALSE,
                       show_colnames=FALSE,
                       angle_col=c("270", "0", "45", "90", "315"),
                       fontsize_number=0.8*fontsize,
                       width=40,
                       height=40,
                       filename=paste0(results_path, "Correlation_vst.jpg"))
  }
}

#****************************************************************************#
#                          STEP 7: PREPARE QC PLOTS                          #
#****************************************************************************#

plot_qc <- function(dds, meta_data, approach, suffix){  
  
  #****************************************************************************#
  #                                   PLOT PCA                                 #  
  #****************************************************************************#
  
  # NOTE: DESeq2::vst() only needs DESeqDataSet or matrix of counts as input
  # NOTE: DESeq2::plotDispEsts needs DeSeqDataSet with dispersion estimated.
  # So, better to plot these qc plots after running DESeq2()
  
  vst_mat <- DESeq2::vst(dds, blind=FALSE)
  vst_mat <- SummarizedExperiment::assay(vst_mat)
  
  # Calculate PCs
  pca <- prcomp(t(vst_mat))
  
  # Create data frame with metadata, PC1 & PC2 values for ggplot
  df <- cbind(meta_data, pca$x)      
  
  ggplot(df, aes(x=PC1, y=PC2, color=id)) + 
    geom_point() +
    geom_text_repel(label=rownames(df))
  ggplot2::ggsave(filename=paste0("PCA_Plot_overall_", approach, "_", suffix, ".jpg"),
                  plot=last_plot(),
                  device="jpeg",
                  path=results_path,
                  width=7,
                  height=7,
                  units=c("in"),	 
                  dpi=300,
                  limitsize=TRUE,
                  bg="white")
  
  #****************************************************************************#
  #                         PLOT DISPERSION ESTIMATES                          #
  #****************************************************************************#
  
  # Expected results: Higher the mean, lower the dispersion
  # NOTE: The output of plotDispEsts() is not a ggplot2 object.
  
  grDevices::jpeg(filename=paste0(results_path, "Dispersion_Plot_overall_", approach, "_", suffix, ".jpg"),
                  width=8.5,
                  height=11,
                  units=c("in"),
                  quality=75,
                  bg="white",
                  res=300)
  
  DESeq2::plotDispEsts(object=dds,
                       genecol="black",
                       fitcol="red",
                       finalcol="dodgerblue",
                       legend=TRUE,
                       xlab="Mean of Normalized counts",
                       ylab="Dispersion",
                       log="xy",
                       cex=0.45)
  
  grDevices::dev.off()
} 

#******************************************************************************#
#
#******************************************************************************#

# Values should be raw i.e. untransformed. DO NOT use log transformed values etc.
# DO NOT replace NA with 0 etc. Leave NA as they are.
# NA values will be imputed/replaced with average of non-NA values.
# If all values are NA, they will be set to 0.
# NO duplicated genes MUST be present.
impute_with_mean <- function(raw_counts){
  
  # Replace NA with average
  for (j in 1:nrow(raw_counts)){
    
    data <- raw_counts[j,] %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Sample") %>%
      dplyr::left_join(metadata, by = c("Sample" = "Sample"))
    colnames(data) <- c("Sample", "values", "Condition")
    
    # If all values for Reference == NA or Target == NA, set them to 0. 
    # Replace NA with average wherever possible.
    data <- data %>% 
      dplyr::mutate(values = as.numeric(values)) %>%
      dplyr::group_by(Condition) %>% 
      dplyr::mutate(average = mean(values, na.rm=TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(average = dplyr::if_else(is.na(average), 0, average),
                    values = dplyr::if_else(is.na(values), average, values))
    
    data <- data %>% 
      dplyr::select(Sample, values) %>% 
      tibble::column_to_rownames("Sample") %>% 
      t()
    
    if(all(colnames(data) == colnames(raw_counts))){
      raw_counts[j,] <- data[,colnames(raw_counts[j,])] %>% unlist(use.names=FALSE)
    }
  }
  
  imputed_counts <- raw_counts %>% 
    dplyr::mutate(across(.cols=everything(), .fns = as.numeric))
  
  return(imputed_counts)
}

#******************************************************************************#
#
#******************************************************************************#

# Perform normalization before imputation so counts can be compared across samples
# NOTE: meta_data MUST have column "Sample" that matches with column names of
# raw_counts 
# NOTE: meta_data MUST have column "Condition" that defines the groups
quantile_norm <- function(raw_counts, meta_data, quant_norm){
  
  # Perform quantile normalization
  if (quant_norm == TRUE){
    
    # https://doi.org/10.1038/s41598-020-72664-6
    # The above paper recommends to quantile normalize within each group rather
    # than whole dataset
    quant_norm_counts <- data.frame(matrix(data=NA, nrow=nrow(raw_counts), ncol=0))
    
    # Subset raw_counts for each group
    for (c in unique(meta_data$Condition)){
      
      samples <- meta_data %>% 
        dplyr::filter(Condition %in% c) %>% 
        dplyr::select(Sample) %>% 
        unlist(use.names=FALSE)
      
      counts <- raw_counts[,samples]
      counts <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(counts)))
      rownames(counts) <- rownames(raw_counts)
      colnames(counts) <- samples
      
      quant_norm_counts <- dplyr::bind_cols(quant_norm_counts, counts)
    }
  } else {
    quant_norm_counts <- raw_counts
  }
  
  return(quant_norm_counts)
}

#******************************************************************************#
#                              CALCULATE padj AND log2FC
#******************************************************************************#

# Function to calculate pval and log2FoldChange
# norm_counts is matrix with log2 transformed values or non-log transformed values
# DO NOT use log10 transformed values
calc_stats <- function(norm_counts, metadata, Target, Reference,log2_transformed_already){
  
  # Perform t.test
  SYMBOL <- c()
  expt <- c()
  control <- c()
  pval <- c()
  for (j in 1:nrow(norm_counts)){
    
    data <- norm_counts[j,] %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Sample") %>%
      dplyr::left_join(metadata, by = c("Sample" = "Sample"))
    colnames(data) <- c("Sample", "values", "Condition")
    
    # Use F.test() to determine if variances are equal or not.
    # NOTE: p < 0.05 => unequal variance
    
    # NOTE: If Iso <- c(NA, NA, 18) and IP <- c(17, NA, 18),
    # var.test and t.test with throw error since Iso has only 1 value.
    
    # NOTE: If Iso <- c(1,1,1) and IP <- c(1,1,1),
    # t.test will throw error "data are essentially constant".
    
    # NOTE: If Iso <- c(0,0,0) and IP <- c(1,1,1),
    # t.test will throw error "data are essentially constant".
    
    if (sum(!is.na(data[data$Condition == Reference, ]$values)) > 1 &
        sum(!is.na(data[data$Condition == Target, ]$values)) > 1 & 
        (length(unique(data[data$Condition == Reference, ]$values)) + 
         length(unique(data[data$Condition == Target, ]$values)) > 2)){
      #   f_test <- var.test(formula = values ~ Condition, 
      #                      data = data,
      #                      alternative = "two.sided")
      #   
      #   if(!is.na(f_test$p.value)){
      # if (f_test$p.value < 0.05){
      # t_test <- t.test(formula = values ~ Condition,
      #                  data = data,
      #                  alternative = "two.sided",
      #                  var.equal = FALSE)
      # }
      
      # # Remove outliers
      # ref_data <- data[data$Condition == Reference, ]$values
      # low_ref <- quantile(ref_data, na.rm=TRUE)[2] - 1.5*IQR(ref_data, na.rm=TRUE)
      # high_ref <- quantile(ref_data, na.rm=TRUE)[3] + 1.5*IQR(ref_data, na.rm=TRUE)
      # 
      # target_data <- data[data$Condition == Target, ]$values
      # low_tar <- quantile(target_data, na.rm=TRUE)[2] - 1.5*IQR(target_data, na.rm=TRUE)
      # high_tar <- quantile(target_data, na.rm=TRUE)[3] + 1.5*IQR(target_data, na.rm=TRUE)
      # 
      # data <- data %>%
      #   dplyr::filter(!(Condition == Reference & (values > high_ref | values < low_ref))) %>%
      #   dplyr::filter(!(Condition == Target & (values > high_tar | values < low_tar)))
      
      
      # Calculate p values, mean expression
      t_test <- t.test(formula = values ~ Condition, 
                       data = data,
                       alternative = "two.sided",
                       var.equal = FALSE)
      
      if (grepl(Reference, names(t_test$estimate[1]))){
        SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[1]])
        expt <- c(expt, t_test$estimate[[2]])
      } else if (grepl(Reference, names(t_test$estimate[2]))) {
        SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[2]])
        expt <- c(expt, t_test$estimate[[1]])
      }
    } else{
      SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
      pval <- c(pval, 1) # Note: DO NOT SET to NA. It will increase padj.
      control <- c(control, mean(data[data$Condition == Reference, ]$values, na.rm=TRUE))
      expt <- c(expt, mean(data[data$Condition == Target, ]$values, na.rm=TRUE))
    }
  }
  
  stats_df <- data.frame(SYMBOL, expt, control, pval)
  stats_df$padj <- stats::p.adjust(p = stats_df$pval, method = "fdr", n = length(stats_df$pval))
  if(log2_transformed_already){
    stats_df$log2FoldChange <- stats_df$expt - stats_df$control  # if data is already log transformed
  }else{
    stats_df$log2FoldChange <- log(stats_df$expt/stats_df$control, base=2)
  }
  
  result <- norm_counts %>% 
    tibble::rownames_to_column(var = "SYMBOL") %>% 
    dplyr::left_join(stats_df, by = c("SYMBOL" = "SYMBOL"))
  
  return(result)
}

#******************************************************************************#
#                          GSEA ANALYSIS USING FGSEA                           #
#******************************************************************************#

# Enrichment analysis takes differential data from every measured gene and looks
# for pathways displaying significantly coordinated shifts in those values.

# https://www.biostars.org/p/12182/
# https://www.biostars.org/p/17628/
# https://www.reddit.com/r/bioinformatics/comments/11o7mrv/gene_set_enrichment_analysis_do_you_separate_out/?rdt=59356
# https://groups.google.com/g/gsea-help/c/oXsBOAUYnH4
# https://www.biostars.org/p/132575/

# DEGs_df MUST contain columns "SYMBOL", "padj" & "log2FoldChange"
# gmt_file MUST be the original unmodified filename of gmt_file with full 
# path as downloaded from GSEA

fgsea <- function(DEGs_df, gmt_file, annotations, results_path){
  
  #****************************************************************************#
  #     PREPARE A RANKED GENE LIST OF ALL GENES EXPRESSED IN YOUR SAMPLES      #
  #****************************************************************************#
  
  # NOTE: ALL genes MUST be used for this analysis, NOT just DEGs. 
  # NOTE: Genes MUST be ranked i.e. sorted in descending fold change. You can 
  # also rank based on log2FC & p value like: sign(df$log2fc)*(-log10(df$pval)))
  # NOTE: Genes MUST be stored in list format, not as a dataframe.
  # NOTE: No NA MUST be present in SYMBOL column. Else, fgsea::collapsePathways()
  # will give error 
  # "Error in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  : 
  # Not all stats values are finite numbers"
  # You can figure using table(is.na(names(DEGs_list))), 
  # is.infinite(DEGs_df$log2FoldChange) or sapply(DEGs_df, class) to make sure
  # log2FoldChange and padj are numeric
  # NOTE: If you excel file has "inf" in padj or log2FoldChange columns, the
  # column will be read into R as character column instead of numeric. So, remove
  # text values from log2FoldChange and padj columns
  
  DEGs_df <- DEGs_df %>%
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    dplyr::filter(!is.na(padj), !is.na(SYMBOL)) %>%
    dplyr::mutate(log2FoldChange = as.numeric(log2FoldChange),
                  padj = as.numeric(padj)) %>%
    dplyr::arrange(desc(log2FoldChange))
  
  DEGs_list <- DEGs_df$log2FoldChange
  names(DEGs_list) <- DEGs_df$SYMBOL
  
  #****************************************************************************#
  #                         DEFINE scoreType PARAMETER                         #
  #****************************************************************************#
  
  # Define score type in fgseaMultilevel() based on fold change values.
  # NOTE: Use "pos", if you are ONLY interested in activated pathways.
  # NOTE: Use "neg", if you are ONLY interested in inhibited pathways. 
  # NOTE: Else, use "std" for both activated & inhibited pathways. 
  score_type <- dplyr::if_else(max(DEGs_list) > 0 & min(DEGs_list) < 0, "std", 
                               dplyr::if_else(max(DEGs_list) < 0 & min(DEGs_list) < 0, "neg", "pos"))
  
  #****************************************************************************#
  #                      IMPORT YOUR GENE SETS OF INTEREST                     #
  #****************************************************************************#
  
  # NOTE: Unlike clusterProfiler::GSEA(), fgsea::fgseaMultilevel() needs gene 
  # sets in a specific list format
  
  # NOTE: If you run multiple gene sets like C5 and C2 together, padj will not 
  # be significant as there will be too many multiple comparisons.
  
  gmt <- fgsea::gmtPathways(gmt_file)
  gmt_name <- gsub(pattern="^.*/|v2023.*$", replacement="", x=gmt_file)
  
  # From each gene set, remove genes that are absent in your DEGs_list
  for (i in 1:length(gmt)){
    gmt[[i]] <- gmt[[i]][gmt[[i]] %in% names(DEGs_list)]
  }
  
  #****************************************************************************#
  #                                  RUN fGSEA                                 #
  #****************************************************************************#
  
  fgsea <- fgsea::fgseaMultilevel(pathways = gmt,
                                  stats = DEGs_list,
                                  scoreType = score_type,
                                  sampleSize = 101,
                                  minSize = 1,
                                  maxSize = length(DEGs_list) - 1,
                                  eps = 1e-50,
                                  nproc = 0,
                                  gseaParam = 1,
                                  BPPARAM = NULL,
                                  nPermSimple = 10000)
  
  #****************************************************************************#
  #                               FORMAT RESULTS                               #
  #****************************************************************************#
  
  # NOTE: Output of fgsea is a data.table & data.frame. 
  # "leadingEdge" column is a list of genes. 
  # So, DO NOT FORCE the output of fgsea to a dataframe as this will lead to 
  # data loss from "leadingEdge" column & affect plotting using fgsea::plotEnrichment()
  
  # Reformat the output
  gsea_results <- fgsea %>%
    #dplyr::filter()
    dplyr::mutate(abs_NES = abs(NES)) %>%
    dplyr::arrange(desc(abs_NES)) %>%
    dplyr::mutate(Direction = dplyr::if_else(NES > 0, "Upregulated", "Downregulated"))
  
  # If you ordered your gene list based on fold change, then +ve NES indicates
  # that the genes in this gene set are mostly at the top of your gene list
  # (hence, most of them are upregulated) and -ve NES indicates that the genes
  # in this gene set are mostly at the bottom of your gene list (hence, most 
  # of them are downregulated)
  
  # Identify overlapping pathways and collapse them into major pathways
  concise_gsea_results <- fgsea::collapsePathways(fgseaRes = gsea_results,
                                                  pathways = gmt,
                                                  stats = DEGs_list)
  # Filter out overlapping pathways
  concise_gsea_results <- gsea_results %>% 
    dplyr::filter(pathway %in% concise_gsea_results$mainPathways)
  
  # Function to convert genelist to df
  convert_genelist_to_df <- function(gsea_results){  
    
    # Create a dataframe containing genes for each pathway
    max_len <- max(unlist(lapply(X=stringr::str_split(string = gsea_results$leadingEdge, pattern = ","), FUN=length)))
    df <- data.frame(matrix(NA, nrow=max_len))
    
    # Add genes from each pathway to a separate column
    for (i in 1:nrow(gsea_results)){
      l1 <- unlist(stringr::str_split(string = unlist(gsea_results$leadingEdge[[i]]), pattern = ","))
      
      # If leading edge column has ENSEMBL_IDs, convert them to SYMBOL
      if (length(intersect(l1, annotations$ENSEMBL_ID)) > length(intersect(l1, annotations$SYMBOL))){
        l1 <- annotations %>% 
          dplyr::filter(ENSEMBL_ID %in% l1) %>% 
          dplyr::select(SYMBOL) %>% 
          unlist(., use.names=FALSE)
      }
      
      l1 <- c(l1, rep(x=NA, times=max_len-length(l1)))
      df <- dplyr::bind_cols(df, l1)
      colnames(df)[i+1] <- gsea_results$pathway[i]
    }
    
    # Remove the dummy column containing NAs
    df <- df[,-1]
    
    return(df)
  }
  
  if(nrow(gsea_results) > 0){
    gsea_gene_df <- convert_genelist_to_df(gsea_results)
  } else{
    gsea_gene_df <- data.frame(matrix(NA, nrow=1))
  }
  if(nrow(concise_gsea_results) > 0){
    concise_gsea_gene_df <- convert_genelist_to_df(concise_gsea_results)
  } else{
    concise_gsea_gene_df <- data.frame(matrix(NA, nrow=1))
  }
  
  # Save the results in excel file
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="gsea_full")
  openxlsx::writeData(wb, sheet="gsea_full", x=gsea_results, rowNames=FALSE)
  openxlsx::addWorksheet(wb, sheetName="gsea_concise")
  openxlsx::writeData(wb, sheet="gsea_concise", x=concise_gsea_results, rowNames=FALSE)
  openxlsx::addWorksheet(wb, sheetName="gsea_full_genes")
  openxlsx::writeData(wb, sheet="gsea_full_genes", x=gsea_gene_df, rowNames=FALSE)
  openxlsx::addWorksheet(wb, sheetName="gsea_concise_genes")
  openxlsx::writeData(wb, sheet="gsea_concise_genes", x=concise_gsea_gene_df, rowNames=FALSE)
  openxlsx::saveWorkbook(wb, file = paste0(results_path, gmt_name, "fGSEA_Results.xlsx"),
                         overwrite = TRUE)
}

#******************************************************************************#
#                                 VOLCANO PLOT                                 #
#******************************************************************************#

# Function to plot volcano plots
# Input dataframe MUST have columns "log2FC", "padj", "SYMBOL"

plot_volcano <- function(volcano_df, disp_genes, Target, Reference, log2_cutoff, padj_cutoff, file_suffix, output_path){
  
  # Categorize the data points
  volcano_df <- volcano_df %>%
    dplyr::mutate(Direction = dplyr::case_when(padj < padj_cutoff & log2FC > log2_cutoff ~ paste0("Up in ", Target),
                                               padj < padj_cutoff & log2FC < -log2_cutoff ~ paste0("Up in ", Reference),
                                               TRUE ~ "Not Significant"),
                  padj = dplyr::case_when(is.na(padj) ~ 0,
                                          TRUE ~ padj),
                  Significance = dplyr::case_when(abs(log2FC) >= log2_cutoff & padj <= 0.05 & padj > 0.01 ~ "FDR < 0.05",
                                                  abs(log2FC) >= log2_cutoff & padj <= 0.01 & padj > 0.001 ~ "FDR < 0.01",
                                                  abs(log2FC) >= log2_cutoff & padj <= 0.001  ~ "FDR < 0.001",
                                                  TRUE ~ "Not Significant"))
  
  # Define the limits of the x-axis
  x_limits <- as.vector(quantile(volcano_df$log2FC, probs = c(0,1)))
  if (is.infinite(max(x_limits))){
    x_limits <- as.vector(quantile(volcano_df$log2FC, probs = c(0,0.999)))
  }
  
  # Define the limits of the y-axis
  y_limits <- as.vector(quantile(-log10(volcano_df$padj), na.rm = TRUE, probs = c(0,1)))
  if (is.infinite(max(y_limits))){
    y_limits <- as.vector(quantile(-log10(volcano_df$padj), na.rm = TRUE, probs = c(0,0.999)))
  }
  
  
  
  # Define the standard colors
  volcano_palette <- c("#808080", RColorBrewer::brewer.pal(11, "RdBu")[c(11)], RColorBrewer::brewer.pal(11, "RdBu")[c(1)])
  names(volcano_palette) <- c("Not Significant", paste0("Up in ", Reference), paste0("Up in ", Target))
  
  # if (color_by == "Significance"){
  #   volcano_palette <- c(viridis(4)[4], viridis(4)[3], viridis(4)[2], viridis(4)[1])
  #   names(volcano_palette) <- c("Not Significant", "FDR < 0.05", "FDR < 0.01", "FDR < 0.001")
  # } else if (color_by == "Direction" & same_color == "TRUE"){
  #   x <- unique(volcano_df$Direction)
  #   volcano_palette <- dplyr::case_when(grepl(Target, x) ~ "orange", 
  #                                       grepl(Reference, x) ~ "purple", 
  #                                       TRUE ~ "grey")
  #   names(volcano_palette) <- unique(volcano_df$Direction)
  # } else if (color_by == "Direction" & same_color == "FALSE"){
  #   x <- unique(volcano_df$Direction)
  #   volcano_palette <- c("grey", my_palette[1:(length(x)-1)])
  #   names(volcano_palette) <- unique(volcano_df$Direction)
  # }
  
  # NOTE: DO NOT USE labels for defining colors due to reasons below. 
  # RECOMMEND using a named vector.
  # NOTE: If using labels, sort labels in alphabetical order and then assign 
  # color because R by default will arrange the labels in alphabetical order 
  # first and then match them to colors indicated in values vector and then 
  # color the plot. The coloring in the legend is however dependent on the 
  # order of labels vector and values vector. To understand, create a plot first 
  # using the sort and then without the sort(). 
  
  # We set the max size of a point as 0.2
  max_val <- max(abs(volcano_df$log2FC))*max(-log10(volcano_df$padj))
  
  # Plot the volcano plot
  p <-  ggplot2::ggplot(data = volcano_df, 
                  aes(x = log2FC, 
                      y = -log10(padj),
                      #fill = Direction,
                      fill = log2FC,
                      size = abs(log2FC)*-log10(padj)*0.2/max_val)) +
    
    # Plot dot plot
    # If you want all points to be fixed size, use size=0.2 within geom_point()
    # If you want to adjust point size based on column in volcano_df, declare in global aes() within ggplot()
    ggplot2::geom_point(col="black", 
                        shape=21,
                        stroke=0.5,
                        position=position_jitter(h=0.01,w=0.01)) +
    
    # Define the theme of plot
    ggplot2::theme_classic() +
    
    # Adjust size of symbols in legend 
    ggplot2::guides(shape = guide_legend(override.aes = list(size = 3)),
                    size = "none", 
                    fill = guide_colourbar(theme = theme(legend.key.width  = unit(0.75, "lines"),
                                                         legend.key.height = unit(10, "lines"),
                                                         legend.ticks = element_blank(),
                                                         legend.frame = element_rect(colour = "Black", linewidth = 1)))) +
    
    # Define the color of the dots
    scale_fill_gradient2(low="#007ba7", mid="Black", high="Yellow") +
    # scale_color_viridis_d() +
    #scale_fill_manual(values = volcano_palette) +
    
    # Define the axis, plot headings
    ggplot2::labs(x = expression("log"[2]*"FC"),
                  y = expression("-log"[10]*"padj"),
                  #fill = "Direction",
                  fill = "log2FC",
                  size = "-log10(padj)*log2FC",
                  title = paste0(Target, " vs ", Reference)) +
    
    # Draw line to mark the cutoffs
    geom_vline(xintercept = c(-log2_cutoff,log2_cutoff), color = "black", linetype = "dotted", linewidth = 0.5) +
    geom_hline(yintercept = -log10(padj_cutoff),         color = "black", linetype = "dotted", linewidth = 0.5) +
    
    # Define x-axis start and end
    coord_cartesian(xlim = c(-ceiling(max(abs(x_limits))), ceiling(max(abs(x_limits))))) +
    
    # Define the axis tick marks
    #scale_x_continuous(breaks = seq(-ceiling(max(abs(x_limits))), ceiling(max(abs(x_limits))), 1)) +
    #scale_y_continuous(breaks = seq(0, ceiling(y_limits[2]/10)*10, 20)) +
    
    # Adjust font size, style
    my_theme +
    
    # Add gene labels
    geom_text_repel(data = volcano_df %>% dplyr::filter(SYMBOL %in% disp_genes),
                    mapping = aes(label = SYMBOL),
                    #size = 2,
                    force = 0.5,
                    point.size = 1,
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
                    position = position_quasirandom())
  
  # # Save the plot
  # ggplot2::ggsave(filename =  paste0("Volcano_Plot",  file_suffix, ".pdf"),
  #                 plot = last_plot(),
  #                 device = "pdf",
  #                 path = output_path,
  #                 width = 7,
  #                 height = 7,
  #                 units = c("in"),	 
  #                 dpi = 300,
  #                 limitsize = TRUE,
  #                 bg = "white")
  
  return(p)
}


#******************************************************************************#
#                                 VENN DIAGRAM                                 #
#******************************************************************************#

# Input dataframe MUST have maximum of 5 columns

plot_venn <- function(data, path, suffix){
  
  plot_title <- suffix
  
  # Number of columns being plotted
  ncol <- ncol(data)
  
  # Fix column names
  colnames(data) <- base::gsub("_", " ", colnames(data))
  colnames(data) <- base::gsub("\\.", " ", colnames(data))
  
  # Declare cat.pos and cat.dis
  if (ncol == 4){
    pos <- c(330, 15, 330, 15)
    dist <- c(0.27, 0.25, 0.15, 0.13)
    cex = 2
    palette1 <- c("#C8E7F5", "#00008C", "#F6D2E0", "#E75480")        # for 5B, 6B
  } else if (ncol == 3){
    pos <- c(0, 0, 180)
    dist <- c(0.1, 0.1, 0.1)
    cex = 2
    palette1 <- c("#C8E7F5", "#F6D2E0", "#db6d00")    # c("#00008C", "#E75480", "#db6d00") # for 7A
  } else if (ncol == 2){
    pos <- c(0, 0)
    dist <- c(0.05, 0.05)
    cex = 2.75
    palette1 <- c("#C8E7F5", "#db6d00")                              # for 7B
  } else if (ncol == 1){
    pos <- c(0)
    dist <- c(0.1)
    cex = 2.75
    palette1 <- c("#F6D2E0")                                         # for 7B
  }
  
  # Create a dataframe where we store the wrapped column names
  annotation <- data.frame(colnames(data))
  colnames(annotation) <- c("Labels")
  annotation <- annotation %>% 
    dplyr::mutate(Labels = stringr::str_wrap(Labels, 6))  # Set 10 for #5C. Normally, 15
  
  # Convert the data frame to a named list
  genes <- base::vector(mode = "list", length = ncol(data))
  names(genes) <- annotation$Labels
  
  for (i in 1:ncol(data)){
    
    # remove NA values and create a list of genes for each label
    genes[[i]] <- data[!is.na(data[i]),i]
  }
  
  # Plot the venn diagram
  VennDiagram::venn.diagram(x = genes,
                            main = plot_title, 
                            category.names = annotation$Labels,
                            filename = paste0(path, "Venn_Diagram_", suffix, ".tiff"),
                            output = TRUE,
                            scaled =FALSE,
                            imagetype = "tiff",
                            height = 11, 
                            width = 11,
                            units = "in",
                            resolution = 600,
                            compression = "lzw",
                            margin = 0.3,    #amount of white space around Venn Diagram in grid units
                            
                            # Formatting the shapes of venn diagram
                            lwd = 1.5,                 #thickness of line
                            lty = 1,                   #type of line
                            col = "black",             #color of line
                            
                            # Formatting numbers inside venn diagram
                            cex = cex,                 #font size (2 or 2.75)
                            fontface = "bold",         #font style
                            fontfamily = "sans",       #font type
                            
                            # Formatting title of venn diagram
                            main.cex = 2,              #font size
                            main.fontface = "bold",    #font style
                            main.fontfamily = "sans",  #font type
                            main.col = "black",        #font color
                            
                            # Formatting category of venn diagram
                            cat.cex = 2,               #font size
                            cat.fontface = "bold",     #font style
                            cat.fontfamily = "sans",   #font type
                            cat.col = palette1,  #"black",  #c("#00008C", "#00008C", "#E75480", "#E75480"),
                            
                            # Formatting colors of venn diagram
                            fill = palette1,
                            alpha = rep(0.5, ncol), #0.5=50% transparency, 1=0% transparency
                            #cat.default.pos = "outer",    
                            
                            cat.pos = pos,    
                            cat.dist = dist, 
                            disable.logging = TRUE,
                            ext.text = TRUE)
  
  #******************************************************************************#
  #                          SAVE THE OVERLAPPING GENES                          #
  #******************************************************************************#
  
  # Save the list of overlapping genes. NOTE: You need to manually figure out
  # which genes belong to which overlap based on number of genes overlapping
  overlap <- VennDiagram::calculate.overlap(x = genes)
  
  # Identify maximum number of genes present in any overlap
  max = max(lengths(overlap))
  
  # Create an dataframe of size length(overlap), max with NAs
  results = data.frame(matrix("", nrow = max, ncol = length(overlap)))
  
  rownames(results) <- paste0("Gene#", seq(max))
  colnames(results) <- paste0("Intersection#", seq(length(overlap)))
  
  # Populate the dataframe with gene names
  for (i in 1:length(overlap)){
    if (length(overlap[[i]]) > 0){
      for (j in 1:length(overlap[[i]])){
        results[[j,i]] <- overlap[[i]][j]
      }
    }
  }
  
  # Save the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Output")
  openxlsx::writeData(wb, sheet = "Output", x = results, rowNames = TRUE)
  openxlsx::addWorksheet(wb, sheetName = "Input")
  openxlsx::writeData(wb, sheet = "Input", x = data, rowNames = FALSE)
  openxlsx::saveWorkbook(wb, file = paste0(path, "Overlap_", suffix, ".xlsx"), overwrite = TRUE,  returnValue = FALSE)
}

#******************************************************************************#
#                       FUNCTIONS TO PLOT SURVIVAL CURVE                       #
#******************************************************************************#

wrangle_data <- function(expr_df, stratify_criteria, prefix){
  
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
                "Tarone_Ware.early" = c(), 
                "Peto_Peto.early" = c(),  
                "modified_Peto_Peto" = c(), 
                "Fleming_Harrington" = c())
  
  # STEP 1: Calculate cutoffs
  # If cutoffs need to be calculated for each group, subset the expr_df and pass
  # it to calculate_cutoffs(). Else, pass entire expr_df to calculate_cutoffs()
  if (multiple_cutoff == TRUE & !is.na(split_by)) {
    
    # Create empty dataframe to store results of calculate_cutoffs() for each group
    survival_data <- data.frame(model = " ")
    
    # Calculate cutoffs for each group in split_by variable
    for (x in (expr_df %>% dplyr::distinct(get(split_by)))[[1]]) {
      mat <- expr_df %>% 
        dplyr::filter(get(split_by) == x) %>% 
        dplyr::select("Sample_ID", "Sex", "Time", "Status", all_of(split_by), all_of(gene))
      mat <- calculate_cutoffs(mat, x)
      
      # Save the data from output of calculate_cutoffs()
      survival_data        <- dplyr::bind_rows(survival_data, mat[[1]])
      stats$gene           <- c(stats$gene,          mat[[2]]$gene)
      stats$group          <- c(stats$group,         mat[[2]]$group)
      stats$lower_cutoff   <- c(stats$lower_cutoff,  mat[[2]]$lower)
      stats$middle_cutoff  <- c(stats$middle_cutoff, mat[[2]]$middle)
      stats$upper_cutoff   <- c(stats$upper_cutoff,  mat[[2]]$upper)
    }
    
    # Populate the model variable by concatenating "Expression" and "split_by"
    survival_data <- survival_data %>%
      dplyr::mutate(model = paste0(Expression, "_", get(split_by))) %>%
      dplyr::filter(!is.na(get(split_by)))
  } else{
    mat <- expr_df
    x <- "NA"
    mat <- calculate_cutoffs(mat, x)
    
    # Save the data from output of calculate_cutoffs()
    survival_data <- mat[[1]]
    stats$gene <- c(stats$gene, mat[[2]]$gene)
    stats$group <- c(stats$group, mat[[2]]$group)
    stats$lower_cutoff <- c(stats$lower_cutoff, mat[[2]]$lower)
    stats$middle_cutoff <- c(stats$middle_cutoff, mat[[2]]$middle)
    stats$upper_cutoff <- c(stats$upper_cutoff, mat[[2]]$upper)
    
    # Rename the column "Expression" to "model"
    survival_data <- survival_data %>%
      dplyr::rename(model = Expression)
  }
  
  # STEP 2: Calculate survival stats
  # If each group has to be plotted in separate plots, subset the survival_data
  # and pass it to plot_survival(). Else, pass entire survival_data to 
  # plot_survival().
  if (combine_plot == "FALSE") {
    if (!is.na(split_by)){
      for (x in (expr_df %>% dplyr::distinct(get(split_by)))[[1]]) {
        s_data <- survival_data %>% dplyr::filter(get(split_by) == x)
        cox_stats <- plot_survival(s_data, x, prefix)
        
        # Save the data from output of plot_survival()
        stats$HR                  <- c(stats$HR,                  cox_stats$HR )
        stats$CI_lower            <- c(stats$CI_lower,            cox_stats$CI_lower)
        stats$CI_upper            <- c(stats$CI_upper,            cox_stats$CI_upper)
        stats$logrank             <- c(stats$logrank,             cox_stats$logrank)
        stats$reg_logrank.late    <- c(stats$reg_logrank.late,    cox_stats$reg_logrank.late)
        stats$Gehan_Breslow.early <- c(stats$Gehan_Breslow.early, cox_stats$Gehan_Breslow.early)
        stats$Tarone_Ware.early   <- c(stats$Tarone_Ware.early,   cox_stats$Tarone_Ware.early)
        stats$Peto_Peto.early     <- c(stats$Peto_Peto.early,     cox_stats$Peto_Peto.early)
        stats$modified_Peto_Peto  <- c(stats$modified_Peto_Peto,  cox_stats$modified_Peto_Peto)
        stats$Fleming_Harrington  <- c(stats$Fleming_Harrington,  cox_stats$Fleming_Harrington)
      }
    } else {
      s_data <- survival_data
      x <- "NA"
      cox_stats <- plot_survival(s_data, x, prefix)
      
      # Save the data from output of plot_survival()
      stats$HR                  <- c(stats$HR,                  cox_stats$HR )
      stats$CI_lower            <- c(stats$CI_lower,            cox_stats$CI_lower)
      stats$CI_upper            <- c(stats$CI_upper,            cox_stats$CI_upper)
      stats$logrank             <- c(stats$logrank,             cox_stats$logrank)
      stats$reg_logrank.late    <- c(stats$reg_logrank.late,    cox_stats$reg_logrank.late)
      stats$Gehan_Breslow.early <- c(stats$Gehan_Breslow.early, cox_stats$Gehan_Breslow.early)
      stats$Tarone_Ware.early   <- c(stats$Tarone_Ware.early,   cox_stats$Tarone_Ware.early)
      stats$Peto_Peto.early     <- c(stats$Peto_Peto.early,     cox_stats$Peto_Peto.early)
      stats$modified_Peto_Peto  <- c(stats$modified_Peto_Peto,  cox_stats$modified_Peto_Peto)
      stats$Fleming_Harrington  <- c(stats$Fleming_Harrington,  cox_stats$Fleming_Harrington)
    }
  }  else{
    s_data <- survival_data
    x <- "NA"
    cox_stats <- plot_survival(s_data, x, prefix)
    
    # Save the data from output of plot_survival()
    stats$HR                  <- c(stats$HR,                  cox_stats$HR )
    stats$CI_lower            <- c(stats$CI_lower,            cox_stats$CI_lower)
    stats$CI_upper            <- c(stats$CI_upper,            cox_stats$CI_upper)
    stats$logrank             <- c(stats$logrank,             cox_stats$logrank)
    stats$reg_logrank.late    <- c(stats$reg_logrank.late,    cox_stats$reg_logrank.late)
    stats$Gehan_Breslow.early <- c(stats$Gehan_Breslow.early, cox_stats$Gehan_Breslow.early)
    stats$Tarone_Ware.early   <- c(stats$Tarone_Ware.early,   cox_stats$Tarone_Ware.early)
    stats$Peto_Peto.early     <- c(stats$Peto_Peto.early,     cox_stats$Peto_Peto.early)
    stats$modified_Peto_Peto  <- c(stats$modified_Peto_Peto,  cox_stats$modified_Peto_Peto)
    stats$Fleming_Harrington  <- c(stats$Fleming_Harrington,  cox_stats$Fleming_Harrington)
  }
  
  return(list(survival_data, stats))
}

calculate_cutoffs <- function(df, group){
  
  # Identify upper & lower cutoffs based on stratify_criteria
  #*************************Split samples by median**************************#
  if(stratify_criteria == "m"){
    quartiles <- stats::quantile(x = df[[gene]], 
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    
    cutoff_lower_end <- quartiles[[3]]
    cutoff_upper_end <- quartiles[[3]]
    cutoff_middle <- "NA"
  }
  
  #****************Split samples into top and bottom tertiles****************#
  if(stratify_criteria == "t"){
    tertiles <- stats::quantile(x = df[[gene]],
                                probs = c(0, 0.33, 0.66, 1),
                                na.rm=TRUE)
    
    cutoff_lower_end <- tertiles[[2]]
    cutoff_upper_end <- tertiles[[3]]
    cutoff_middle <- "NA"
  }
  
  #***************Split samples into top and bottom quartiles****************#
  if(stratify_criteria == "q"){
    quartiles <- stats::quantile(x = df[[gene]], 
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    
    cutoff_lower_end <- quartiles[[2]]
    cutoff_upper_end <- quartiles[[4]]
    cutoff_middle <- quartiles[[3]]
  }
  
  #*********************Split expression range by thirds*********************#
  if(stratify_criteria == "th"){
    quartiles <- stats::quantile(x = df[[gene]], 
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    iqr <- stats::IQR(x = df[[gene]],
                      na.rm=TRUE)
    
    # Normal range of expression values lie between cutoff_lower & cutoff_upper
    cutoff_upper <- quartiles[[4]]+1.5*iqr
    cutoff_lower <- dplyr::if_else(quartiles[[1]]-1.5*iqr > 0, quartiles[[1]]-1.5*iqr, 0)
    
    # Based on normal range of expression, identify onethird & twothird cutoff
    cutoff_lower_end <- cutoff_lower + (cutoff_upper-cutoff_lower)/3
    cutoff_upper_end <- cutoff_lower + (cutoff_upper-cutoff_lower)*2/3
    cutoff_middle <- "NA"
  }
  
  #***************Split expression range using optimal cutoff****************#
  if(stratify_criteria == "o"){
    
    quartiles <- stats::quantile(x = df[[gene]], 
                                 probs = c(0, 0.25, 0.50, 0.75, 1),
                                 na.rm=TRUE)
    
    # Sometimes quartiles will look like: 
    # 0%       25%      50%      75%     100% 
    # 0.000000 0.000000 0.000000 0.000000 3.495493 
    # In such cases, surv_cutpoint() will fail. So, we add extra if() here.
    if (quartiles[[4]] > quartiles[[2]]){
      res.cut <- survminer::surv_cutpoint(data = df,
                                          time = "Time",
                                          event = "Status",
                                          variables = gene)
      
      cutoff_lower_end <- res.cut$cutpoint$cutpoint
      cutoff_upper_end <- res.cut$cutpoint$cutpoint
      cutoff_middle <- "NA"
    } else{
      #cat("Surv cutpoint unable to detect optimum cutoff")
      cutoff_lower_end <- "NA"
      cutoff_upper_end <- "NA"
      cutoff_middle <- "NA"
    }
  }
  
  # Categorize the data based on above cutoffs
  if (all_quartiles == "TRUE"){
    df <- df %>% 
      dplyr::mutate(Expression = dplyr::if_else(get(gene) > cutoff_upper_end, "HIGH", 
                                                dplyr::if_else(get(gene) < cutoff_lower_end, "LOW",
                                                               dplyr::if_else(get(gene) < cutoff_middle, "MED_LOW",
                                                                              "MED_HIGH"))))
  } else{
    df <- df %>% 
      dplyr::mutate(Expression = dplyr::if_else(get(gene) > cutoff_upper_end, "HIGH", 
                                                dplyr::if_else(get(gene) <= cutoff_lower_end, "LOW", "MID"))) %>%
      dplyr::filter(Expression != "MID")
  }
  
  # # Print the cutoffs
  # cat("\nGene         :", gene)
  # cat("\nGroup        :", group)
  # cat("\nLower cutoff :", cutoff_lower_end)
  # cat("\nUpper cutoff :", cutoff_upper_end)
  # cat("\nMiddle cutoff:", cutoff_middle)
  
  # Create a list to store cutoff values
  ls <- list("group" = c(), 
             "gene" = c(), 
             "lower" = c(), 
             "upper" = c(), 
             "middle" = c())
  
  ls$group <- c(group)
  ls$gene <- c(gene)
  ls$lower <- c(cutoff_lower_end)
  ls$upper <- c(cutoff_upper_end)
  ls$middle <- c(cutoff_middle)
  
  # Return the df
  return(list(df, ls))
  
}

plot_survival <- function(survival_data, group, prefix){
  
  # If all samples belong to one group (like LOW or HIGH or males or female),
  # then quit the function as comparison cannot be done
  if (nrow(survival_data %>% dplyr::count(model)) > 1){
    
    # Create a survival object where Alive = 0, Dead = 1
    survival_object <- survival::Surv(time = survival_data$Time,
                                      event = survival_data$Status,
                                      type = "right",
                                      origin = 0)
    
    # Create a formula for plotting survival curve
    survival_formula <- survival_object ~ model
    
    # Create a fit for survival curve. survfit() gives error in ggsurvplot(). Use surv_fit()
    survival_curve <- survminer::surv_fit(formula = survival_formula,
                                          data = survival_data,
                                          type = "kaplan-meier",
                                          group.by = NULL,
                                          match.fd = FALSE)
    
    # Check summary of the survival curve with time duration of our interest
    #cat("\nRange of survival (months):", range(survival_data$Time, na.rm=TRUE), "\n")
    # base::summary(survival_curve, times = base::seq(from = floor(range(survival_data$Time, na.rm=TRUE)[[1]]),
    #                                                 to = ceiling(range(survival_data$Time, na.rm=TRUE)[[2]]),
    #                                                 by = 3))
    
    # Create a Cox model for the survival curve and calculate stats
    cox_model <- survival::coxph(formula = survival_formula,
                                 data = survival_data)
    #print(summary(cox_model))
    #cat("\n")
    
    # Calculate HR, 95% CI for HR, p-val
    # NOTE: Variable mentioned here is the numerator in h1(t)/h0(t).
    # The reference variable h0(t) will not be mentioned in co-efficients.
    # Make sure this is not the reference level i.e. low expression. If this is
    # the reference, then we need to reverse the HR ratio, legend labels
    #print(names(cox_model$coefficients))  
    
    # If all samples belong to more than 2 groups (like LOW, MID, HIGH), then
    # we cannot have survival stats. SO, we set them to 0.
    if (nrow(survival_data %>% dplyr::count(model)) == 2){
      # Store HR and CI
      if (stringr::str_detect(names(cox_model$coefficients), reference)){
        HR <- round(exp(-cox_model$coefficients[[1]]), 2)
        CI <- round(exp(-confint(cox_model)), 2)
        CI_1 <- CI[1]
        CI[1] <- CI[2]
        CI[2] <- CI_1
      } else {
        HR <- round(exp(cox_model$coefficients[[1]]),2)
        CI <- round(exp(confint(cox_model)),2)
      }
      
      # Store pvalues by different methods
      pvals <- c()
      for (test.method in c("survdiff", "1", "n", "sqrtN", "S1","S2", "FH_p=1_q=1")){
        
          p_val <- survminer::surv_pvalue(fit = survival_curve,
                                          method = test.method,
                                          test.for.trend = FALSE,
                                          combine = FALSE)
          pvals <- c(pvals, p_val[[2]])
      } 
    } else {
      HR <- 0
      CI <- c(0, 0)
      pvals <- c(0,0,0,0,0,0,0) # 7 pvalues for each method
    }
    
    # Plot the survival curve using survminer::ggsurvplot() instead of base::plot()
    # ggsurvplot() produces a list of ggplot objects: a survival curve and a risk table
    # Saving it using cowplot() first and then using ggsave() works nicely as
    # compared to saving directly using ggsave()
    if(plot_curve){
      
      # Plot the survival curve
      legend_label <- survival_data %>% 
        dplyr::count(model) %>% 
        dplyr::select(model) %>% 
        unlist(.,use.names=FALSE)
      
      # We identify proper breaks based on max duration of the dataset
      # We want a maximum of 10 timepoint intervals that are multiples of 12
      max_time <- max(survival_data$Time,na.rm=TRUE)
      n <- floor(max_time/10/12)*12
      if(max_time %/% n <= 10){
        breaks <- n
      } else{
        breaks <- n+12
      }
      
      survival_plot <- survminer::ggsurvplot(fit = survival_curve,
                                             pval = FALSE,
                                             conf.int = confidence_interval,
                                             risk.table = plot_risk_table,
                                             legend.title = legend_title,
                                             legend.labs = legend_label,
                                             palette = color_palette,
                                             break.time.by = breaks,           # break X axis in time intervals by 3 months
                                             xlab = "Time (Months)",           # customize X axis label
                                             ylab = "Survival Probability",    # customize Y axis label
                                             title = dplyr::if_else(gene == "combined.exp", "", gene),
                                             risk.table.y.text.col = T,        # colour risk table text annotations
                                             size = 2,
                                             censor.size = 9)
        survival_plot$table <- survival_plot$table + 
        coord_cartesian(x=c(0,ceiling(max_time/breaks)*breaks), clip = "off")
        survival_plot$plot <- survival_plot$plot + 
          coord_cartesian(x=c(0,ceiling(max_time/breaks)*breaks), clip = "off")
        
      
      # Plot p and HR value
      method_plot <- "log-rank"
      p_plot <- pvals[1]  
      
      grob1 <- grobTree(textGrob(label = paste0("p = ", formatC(p_plot, format = "e", digits = 1),
                                                "\nHR = ", round(HR,1), " [", round(CI[1],1), ", ", round(CI[2],1), "]",
                                                "\nMethod = ", method_plot),
                                 x = 0.50,
                                 y = 0.90,
                                 hjust = 0,
                                 gp = grid::gpar(fontfamily="Times", fontface="bold", col="black", fontsize=15)))
      
      # Add p values and HR values to plot
      survival_plot$plot <- survival_plot$plot %++%
        ggplot2::annotation_custom(grob1)
      
      cowplot::plot_grid(plotlist = survival_plot,
                         align = "hv",
                         axis = "tblr",
                         nrow = 2,  
                         ncol = 1, 
                         rel_widths = 1,
                         rel_heights = c(1,0.45),
                         labels = NULL,
                         label_size = 14,
                         label_fontfamily = NULL,
                         label_fontface = "bold",
                         label_colour = NULL,
                         label_x = 0,
                         label_y = 1,
                         hjust = -0.5,
                         vjust = 1.5)
      
      # Save the plot
      ggplot2::ggsave(filename = gsub("/", "-", x=paste0(prefix, "_", group, "_", stratify_criteria, "_", gene, ".tiff")),
                      plot = last_plot(),
                      device = "tiff",
                      path = parent_path,
                      width = 5,
                      height = 7,
                      units = c("in"),
                      dpi = 300,
                      limitsize = TRUE,
                      bg = NULL)
    }
  } else {
    HR <- 0
    CI <- c(0, 0)
    pvals <- c(0,0,0,0,0,0,0) # 7 pvalues for each method
  }
  
  # Create a list to store survival stats
  ls <- list("group" = c(), "HR" = c(), "CI_lower" = c(), "CI_upper" = c(), "pvalue" =c(), 
             "logrank" = c(), "reg_logrank.late" = c(), "Gehan_Breslow.early" = c(),
             "Tarone_Ware.early" = c(), "Peto_Peto.early" = c(),  
             "modified_Peto_Peto" = c(), "Fleming_Harrington" = c())
  
  ls$group <- c(group)
  ls$HR <- c(HR)
  ls$CI_lower <- c(CI[1])
  ls$CI_upper <- c(CI[2])
  ls$logrank <- c(pvals[1])
  ls$reg_logrank.late <- c(pvals[2])
  ls$Gehan_Breslow.early <- c(pvals[3])
  ls$Tarone_Ware.early <- c(pvals[4])
  ls$Peto_Peto.early <- c(pvals[5])
  ls$modified_Peto_Peto <- c(pvals[6])
  ls$Fleming_Harrington <-c(pvals[7])
  
  return(ls)
}

# We need to calculate a combined expression value of all genes in the gene
# signature for each sample. Next, we use surv_cutpoint() to classify samples 
# into high and low groups and plot survival curves.

# There are several approaches to calculate the combined expression value. While
# normal z-scaling is logical, https://doi.org/10.1186/gb-2006-7-10-r93 
# recommends a more advanced z-scaling which is implemented below. Refer
# "Measuring gene set expression" in "Materials & stratify_criteria" section.

advanced_Z <- function(gset, eset) {
  
  # Compute z-score for gene set of interest
  ix <- is.element(toupper(rownames(eset)), toupper(gset))
  cat(sum(ix), "\n")
  
  if (sum(ix)>0){
    avg_gset <- base::apply(X=eset[ix,], MARGIN=2, FUN=mean, na.rm=TRUE)
    avg_all <- base::apply(X=eset, MARGIN=2, FUN=mean, na.rm=TRUE)
    sd_all <- base::apply(X=eset, MARGIN=2, FUN=sd, na.rm=TRUE)
    z <- (avg_gset - avg_all)*sqrt(sum(ix))/sd_all
  } else{
    z <- NA
  }
  
  return(z)
}

normal_Z <- function(gset, eset) {
  
  # Compute z-score for gene set of interest
  eset <- eset[gset,]
  a <- t(scale(t(eset)))
  
  z <- colSums(a, na.rm=TRUE) 
  
  return(z)
}

#******************************************************************************#
#                                 HEATMAP PLOT                                 #
#******************************************************************************#

# Function to plot heatmap
# NOTE: The input normalized_counts MUST be a dataframe.  
# (i) First column of normalized_counts MUST be "SYMBOL". Duplicates ok.
# (ii) plot_genes is a vector of genes to be plotted
# (iii) disp_genes is a vector of plotted gene names that must be labelled
# (iv) metadata_column MUST contain the column names specified in anno_columns
# metadata_column MUST have rownames.
# (v) metadata_row MUST contain a column named "SYMBOL" containing gene names 
# without duplicated gene names. metadata_row MUST have rownames.
# (vi) If plotting gene signature score, name the column as "Score". The script
# has been adjusted to use same color as survival curves for HIGH LOW samples

# NOTE: Other variables that MUST be defined are:
# perform_log_transform; perform_scaling; row_clustering; col_clustering;
# row_clustering_alphabetical; col_clustering_alphabetical;
# gaps_in_row; gap_rows; gaps_in_col; gap_columns anno_columns; anno_rows; 
# my_palette; color_by_cols, color_by_rows 
# columns (this defines which column of metadata will be plotted as columns, 
# usually the column containing the sample names)
plot_heatmap <- function(normalized_counts, metadata_column, metadata_row, 
                         plot_genes, disp_genes, file_suffix, file_format, 
                         results_path, bar_width, bar_height, expr_legend){
  
  #****************************************************************************#
  # Format the matrix for heatmap
  #****************************************************************************#
  
  mat <- normalized_counts %>%
    # Keep only genes that need to be plotted
    dplyr::filter(str_to_upper(SYMBOL) %in% str_to_upper(plot_genes)) %>%
    # If there are duplicated genes, keep only data for highest expressing copy
    dplyr::mutate(n = rowSums(.[,-1])) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_max(n) %>%
    dplyr::ungroup() %>%
    # Duplicated genes with 0 expression in all samples still remain, remove them
    dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
    dplyr::select(everything(), -n) %>%
    # Move gene names to rownames
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("SYMBOL") %>%
    # Make sure all values are numeric
    dplyr::mutate(across(.cols = everything(), .fns = as.numeric))
  # mat[, unlist(lapply(mat, is.numeric))]    #alternative to mutate(across())
  
  # Perform log transform if needed. Count data is usually skewed right or left.
  # So, it will mostly be red or blue. log transform to make it less skewed.
  if (perform_log_transform == TRUE){
    mat <- log(1+mat, base = 2)
  }
  
  # Perform scaling for each gene across samples/cells if needed. Since scale() 
  # performs only column scaling, we transpose the dataframe first, so that 
  # genes are on columns & cells are on rows & then perform scaling.
  if (perform_scaling == TRUE){
    mat <- mat %>% t() %>% scale() %>% t()
  }
  
  # Replace NA values with 0
  mat[is.na(mat)] <- 0
  
  # Remove rows which have 0 in all samples
  mat <- mat[rowSums(mat) != 0,]
  
  # Keep only genes that need to be plotted
  mat <- mat[intersect(plot_genes, rownames(mat)), ]
  
  # Keep ONLY samples common in metadata_column and mat
  metadata_column <- metadata_column %>% 
    dplyr::filter(make.names(get(columns)) %in% make.names(colnames(mat)))
  
  # Arrange samples in mat in the same order as in metadata_column. 
  # NOTE: This is important because in the next step we assign rownames to 
  # col_annotation assuming identical sample order between metadata_column & mat
  mat <- mat[,metadata_column[,columns]]
  
  #****************************************************************************#
  # Define column and row annotations
  #****************************************************************************#
  
  if (gtools::invalid(anno_columns)){
    col_annotation <- NA
  } else {
    col_annotation <- metadata_column %>% dplyr::select(all_of(anno_columns))
    rownames(col_annotation) <- colnames(mat)
  }
  
  # Define row annotation for genes
  if (gtools::invalid(anno_rows)){
    row_annotation <- NA
  } else {
    row_annotation <- dplyr::left_join(x = mat %>% as.data.frame() %>% tibble::rownames_to_column("SYMBOL"),
                                       y = metadata_row,
                                       by=c("SYMBOL"="SYMBOL")) %>%
      dplyr::select(everything(), -colnames(mat)) %>%
      dplyr::distinct_at("SYMBOL", .keep_all = TRUE) %>%
      tibble::column_to_rownames("SYMBOL")
  }
  
  #****************************************************************************#
  # Define colors column annotation
  #****************************************************************************#
  
  ann_colors_col <- list()
  ann_colors_row <- list()
  colors <- c("#BF812D", "#35978F", "#C51B7D", "#7FBC41", "#762A83",
              "#E08214", "#542788", "#D6604D", "#4393C3", "#878787",
              "#E41A1C", "#F781BF", "#4DAF4A", "#FFFFBF", "#377EB8", 
              "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#999999",
              "#66C2A5", "#FC8D62", "#000000", "#9E0142", "#1A1A1A",
              "#74006F", "#FFC606", "#F6D2E0", "#C8E7F5")
  #colors <- c("#74006F", "#FFC606")    # Female: Purple, Male:Gold
  #colors <- c("#F6D2E0", "#C8E7F5")    # Female: Pink,   Male:Blue (BBN paper)
  
  if (color_by_cols == TRUE){
    for (y in 1:length(anno_columns)){
      
      # Find number of different elements in each column annotation
      # Eg: Sex has 2: Male, Female; Race has 3: White, Asian, African
      elements <- sort(unique(metadata_column[,anno_columns[y]]))
      palette_n <- length(elements)
      
      if (length(anno_columns) > 1){
        # Create a vector of transparency values. 
        alphas <- base::seq(from=0, to=1, by=1/palette_n)
        
        #Remove alpha=0 as it is always white
        alphas <- setdiff(alphas, 0)
        
        # Create a color palette with different transparencies 
        palette_colors <- colorspace::adjust_transparency(col = colors[y], alpha = alphas[1])
        for (n in 2:palette_n){
          palette_colors <- c(palette_colors, 
                              colorspace::adjust_transparency(col = colors[y], alpha = alphas[n]))
        }
      } else{
        palette_colors <- colors[1:palette_n]
      }
      palette_colors <- rev(palette_colors)
      names(palette_colors) <- rev(elements)  # reversed so order is "LOW","HIGH"
      ann_colors_col <- c(ann_colors_col, list(palette_colors))
    }
    names(ann_colors_col) <- gsub(pattern="\\.", replacement=" ", x=anno_columns)
  }
  
  ann_colors_col$Score[[1]] <- "#0c2c84"
  ann_colors_col$Score[[2]] <- "#d73027"
  
  if (color_by_rows == TRUE){
    for (y in 1:length(unique(metadata_row[,anno_rows]))){
      
      # Find number of different elements in each row annotation
      # Eg: Sex has 2: Male, Female; Race has 3: White, Asian, African
      elements <- sort(unique(unlist(mat[y,], use.names=FALSE)))
      palette_n <- length(elements)
      
      # Create a vector of transparency values. 
      alphas <- base::seq(from=0, to=1, by=1/palette_n)
      
      #Remove alpha=0 as it is always white
      alphas <- setdiff(alphas, 0)
      
      # Create a color palette with different transparencies 
      palette_colors <- colorspace::adjust_transparency(col = colors[y], alpha = alphas[1])
      for (n in 2:palette_n){
        palette_colors <- c(palette_colors, 
                            colorspace::adjust_transparency(col = colors[y], alpha = alphas[n]))
      }
      names(palette_colors) <- elements
      ann_colors_row <- c(ann_colors_row, list(palette_colors))
    }
    names(ann_colors_row) <- gsub(pattern="\\.", replacement=" ", x=rownames(mat))
  }  
  
  # $Sample
  # FB1       FB2       FB3       FB4       FB5        FC       MB1     
  # "#BF812D" "#35978F" "#C51B7D" "#7FBC41" "#762A83" "#E08214" "#542788"
  # $Sex
  # Female      Male
  # "#9E0142" "#E41A1C"
  # $Condition
  # Tumor    Normal
  # "#377EB8" "#4DAF4A"
  
  #****************************************************************************#
  # Determine breaks for heatmap color scale.
  # breaks correspond to numerical ranges for color palette's bins 
  # i.e. 0 to length(my_palette)
  #****************************************************************************#
  
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
  
  #****************************************************************************#
  # Define vectors indicating positions where you want to have gaps in heatmap
  #****************************************************************************#
  
  if (gaps_in_row == TRUE){
    
    # Count() automatically arranges by alphabetical order unfortunately
    # So, we do this manually.
    element_names <- c()
    element_counts <- c()
    c <- 0
    for (i in 1:nrow(row_annotation)){
      if (!(row_annotation[,anno_rows][i] %in% element_names)){
        element_names <- c(element_names, row_annotation[,anno_rows][i])
        element_counts <- c(element_counts, c)
        c <- 1
      } else{
        c <- c+1
      }
    }
    element_counts <- c(element_counts,c)
    
    gaps_row <- data.frame("Description" = element_names, n = element_counts[-1]) %>%
      dplyr::mutate(n = cumsum(n)) %>%
      dplyr::select(n) %>%
      unlist(use.names = FALSE)
  } else{
    gaps_row <- NULL
  }
  
  if (gaps_in_col == TRUE){
    
    # Count() automatically arranges by alphabetical order unfortunately
    # So, we do this manually.
    element_names <- c()
    element_counts <- c()
    c <- 0
    for (i in 1:nrow(col_annotation)){
      if (!(col_annotation[,gap_columns][i] %in% element_names)){
        element_names <- c(element_names, col_annotation[,gap_columns][i])
        element_counts <- c(element_counts, c)
        c <- 1
      } else{
        c <- c+1
      }
    }
    element_counts <- c(element_counts,c)
    
    gaps_col <- data.frame("Description" = element_names, n = element_counts[-1]) %>%
      dplyr::mutate(n = cumsum(n)) %>%
      dplyr::select(n) %>%
      unlist(use.names = FALSE)
  } else {
    gaps_col <- NULL
  }
  
  #****************************************************************************#
  # Define how samples will be arranged in the heatmap
  #****************************************************************************#
  
  if (row_clustering_alphabetical == TRUE){
    mat <- mat[sort(rownames(mat)),]
  } 
  
  if (col_clustering_alphabetical == TRUE){
    mat <- mat[,sort(colnames(mat))]
  } 
  
  if(row_clustering == TRUE){
    # cluster and re-order rows
    rowclust <- hclust(dist(mat))
    reordered <- mat[rowclust$order,]
  } else{
    reordered <- mat
  }
  
  if(col_clustering == TRUE){
    # cluster and re-order columns
    colclust <- hclust(dist(t(mat)))
    reordered <- reordered[, colclust$order]
  } else{
    reordered <- reordered
  }
  
  # If you have 4 groups and want to cluster samples within each of these groups
  if(col_clustering_within_group == TRUE){
    element_names <- unique(col_annotation$Group)
    col_order <- c()
    for (g in element_names){
      temp_mat <- reordered[,rownames(col_annotation)[which(col_annotation == g)]] 
      colclust <- hclust(dist(t(temp_mat)))
      col_order <- c(col_order, colnames(temp_mat[,colclust$order]))
    }
  } else{
    col_order <- colnames(reordered)
  }
  
  if(row_clustering_within_group == TRUE){
    element_names <- unique(row_annotation$Group)
    row_order <- c()
    for (g in element_names){
      temp_mat <- reordered[rownames(row_annotation)[which(row_annotation == g)],]
      rowclust <- hclust(dist(temp_mat))
      row_order <- c(row_order,rownames(temp_mat[rowclust$order,]))
    }
  } else{
    row_order <- rownames(reordered)
  }
  reordered <- reordered[row_order,col_order]
  
  #****************************************************************************#
  # List genes and samples you want to display in the plot
  #****************************************************************************#
  
  display_col <- colnames(reordered)
  display_row <- data.frame("gene" = rownames(reordered)) %>%
    dplyr::mutate(gene = dplyr::case_when(gene %in% disp_genes ~ gene, 
                                          TRUE ~ "")) %>% 
    unlist(use.names=FALSE)
  
  if (ncol(reordered) > 500){
    
    # Save the clustered scores in xlsx
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "Heatmap_matrix")
    openxlsx::writeData(wb, sheet = "Heatmap_matrix", x = t(reordered), rowNames = TRUE)
    openxlsx::saveWorkbook(wb, file = paste0(results_path, "Heatmap_matrix", file_suffix, ".xlsx"), 
                           overwrite = TRUE)
  } else{
    # Save the clustered scores in xlsx
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, sheetName = "Heatmap_matrix")
    openxlsx::writeData(wb, sheet = "Heatmap_matrix", x = reordered, rowNames = TRUE)
    openxlsx::saveWorkbook(wb, file = paste0(results_path, "Heatmap_matrix", file_suffix, ".xlsx"), 
                           overwrite = TRUE)
  }
  
  #****************************************************************************#
  # Plot heatmap
  #****************************************************************************#
  pheatmap::pheatmap(mat = as.matrix(reordered),
                     color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
                     breaks = breaks, 
                     border_color = NA, #"grey90", "white"
                     cellwidth = bar_width,
                     cellheight = bar_height,
                     # cellwidth = dplyr::case_when(ncol(reordered) > 25 ~ NA,
                     #                              ncol(reordered) <= 25 ~ 6), 
                     # cellheight = dplyr::case_when(nrow(reordered) > 25 ~ NA,
                     #                               nrow(reordered) <= 25 ~ 6), 
                     scale = "none",   
                     cluster_rows = row_clustering,   #cluster the rows
                     cluster_cols = col_clustering,   #cluster the columns
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method = "complete",
                     legend = expr_legend,  # set FALSE if it overlaps with annotation legends
                     legend_breaks = NA,
                     legend_labels = NA, 
                     annotation_row = row_annotation,
                     annotation_col = col_annotation,
                     annotation_colors = ann_colors_col,
                     annotation_legend = TRUE,
                     annotation_names_row = FALSE,
                     annotation_names_col = FALSE, #set to TRUE if more than 1 
                     show_rownames = dplyr::if_else(length(disp_genes) < 80, TRUE, FALSE, missing = NULL), 
                     show_colnames = dplyr::if_else(length(unique(display_col)) < 50, TRUE, FALSE, missing = NULL),
                     fontsize = 5, 
                     fontsize_row = 5, 
                     fontsize_col = 5,
                     gaps_row = gaps_row,
                     gaps_col = gaps_col,
                     angle_col = "45",
                     fontsize_number = 0.8*fontsize, 
                     labels_row = display_row, 
                     labels_col = display_col,
                     #width = 11,
                     #height = 11,
                     filename = paste0(results_path, "Heatmap_", file_suffix, ".", file_format))
}

#******************************************************************************#