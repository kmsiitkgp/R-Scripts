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

# Specialized Graph plotting packages
library("pheatmap")             # Needed for making heatmap
library("ggridges")             # Needed for making ridgeplots
library("VennDiagram")          # Needed for making Venn diagram
library("survival")             # Needed for making survival curves
library("survminer")            # Needed for making survival curves, to handle %++% in survival function

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
  
  setAnnotationHubOption("CACHE", "C:/Users/kailasamms/AppData/Local/R/win-library/4.4/AnnotationHub")
  
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
    dplyr::filter(!is.na(Sample)) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames(var="Sample") %>%
    dplyr::mutate(id=get(Variable))
  # dplyr::mutate(id=dplyr::if_else(get(Variable) %in% c(Comparisons$Target, Comparisons$Reference),
  #                                   get(Variable), "ignore")) %>%
  # dplyr::mutate(id=dplyr::if_else(is.na(id), "ignore", id))
  
  return(meta_data)
}  

#****************************************************************************#
#                          STEP 2: REFORMAT READ DATA                          #                                           
#****************************************************************************#

# Remove rows with missng gene names
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
  NA_sample_list <- which(is.na(meta_data[Variable]))
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
deseq2_norm_counts <- function(dds, annotations, approach){  
  
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
  openxlsx::saveWorkbook(wb, file=paste0(results_path, "Normalized_Counts_", approach, ".xlsx"),
                         overwrite=TRUE)
}

combatseq_norm_counts <- function(dds, annotations, approach){
  
  dds <- DESeq2::estimateSizeFactors(dds) 
  
  normalized_counts <- DESeq2::counts(dds, normalized=TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID")
  
  normalized_counts <- add_annotation(normalized_counts, annotations) 
  
  # Save batch corrected normalized counts for entire dataset
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="normalized")
  openxlsx::writeData(wb, sheet="normalized", x=normalized_counts)
  openxlsx::saveWorkbook(wb, file=paste0(results_path, "Normalized_Counts_", approach, ".xlsx"),
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
svaseq_norm_counts <- function(dds, annotations, approach){
  
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
  openxlsx::saveWorkbook(wb, file=paste0(results_path, "Normalized_Counts_", approach, ".xlsx"),
                         overwrite=TRUE)
}

#****************************************************************************#
#                             STEP 6: RUN DESEQ2                             #
#****************************************************************************#

# This function runs DESeq2 and plots volcano plot, heatmap, PCA plot and
# correlation plot.
run_deseq2 <- function(dds, meta_data, annotations, Comparisons, n, lfc.cutoff, padj.cutoff, approach){  
  
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
                DE_levels[grepl(pattern=Comparisons$Target[n], x=DE_levels)],
                DE_levels[grepl(pattern=Comparisons$Reference[n], x=DE_levels)])
  
  
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
  
  ensembl <- length(intersect(results$ID, annotations$ENSEMBL_ID))
  entrez <- length(intersect(results$ID, annotations$ENTREZ_ID))
  symbol <- length(intersect(results$ID, annotations$SYMBOL))
  
  if (ensembl > entrez & ensembl > symbol){
    results <- annotations %>%
      dplyr::right_join(results, by=c("ENSEMBL_ID"="ID"), multiple="all") %>%
      dplyr::select(SYMBOL, CHR, everything(), -c(DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
      dplyr::mutate(SYMBOL=dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                            is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                            TRUE ~ SYMBOL)) %>%
      dplyr::distinct_at("ENSEMBL_ID", .keep_all=TRUE)
  }
  
  if (entrez > ensembl & entrez > symbol){
    results <- annotations %>%
      dplyr::right_join(results, by=c("ENTREZ_ID"="ID"), multiple="all") %>%
      dplyr::select(SYMBOL, CHR, everything(), -c(DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
      dplyr::mutate(SYMBOL=dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                            is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                            TRUE ~ SYMBOL)) %>%
      dplyr::distinct_at("ENTREZ_ID", .keep_all=TRUE)
  }
  
  if (symbol > ensembl & symbol > entrez){
    results <- annotations %>%
      dplyr::right_join(results, by=c("SYMBOL"="ID"), multiple="all") %>%
      dplyr::select(SYMBOL, CHR, everything(), -c(DESCRIPTION, START, END, STRAND, BIOTYPE, BIOTYPE_ENTREZ)) %>%
      dplyr::mutate(SYMBOL=dplyr::case_when(is.na(SYMBOL) & !is.na(ENSEMBL_ID) ~ ENSEMBL_ID,
                                            is.na(SYMBOL) & !is.na(ENTREZ_ID) ~ ENTREZ_ID,
                                            TRUE ~ SYMBOL)) %>%
      dplyr::distinct_at("SYMBOL", .keep_all=TRUE)
  }
  
  # Add unique identifier to each result file
  file_prefix <- gsub(pattern="/", replacement="", x=coeff)
  
  # Save the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="DEGs")
  openxlsx::writeData(wb, sheet="DEGs", x=results, rowNames=FALSE)
  openxlsx::saveWorkbook(wb, file=paste0(results_path, "Results_", file_prefix, "_", approach, "_DEGs.xlsx"), 
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
  
  if (volcano_plot ==TRUE){
    
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
  
  if (cor_plot==TRUE){
    
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

plot_qc <- function(dds, meta_data, approach){  
  
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
  ggplot2::ggsave(filename=paste0("PCA_Plot_overall_", approach, ".jpg"),
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
  
  grDevices::jpeg(filename=paste0(results_path, "Dispersion_Plot_overall_", approach, ".jpg"),
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
