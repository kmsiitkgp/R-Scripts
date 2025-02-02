#!/usr/bin/env Rscript

source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")
source("C:/Users/kailasamms/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

#******************************************************************************#
#                          INSTALL NECESSARY PACKAGES                          #
#******************************************************************************#

# if (!require("BiocManager", quietly=TRUE)) 
#   install.packages("BiocManager")       # Needed to download BioManager packages

# # Data analysis packages
# BiocManager::install("TCGAbiolinks")    # Needed for TCGA data analysis

# # Data wrangling packages
# install.packages("openxlsx")            # Needed for reading, writing xlsx files
# install.packages("dplyr")               # Needed for data wrangling
# install.packages("tibble")              # Needed for advanced data wrangling
# install.packages("stringr")             # Needed for advanced data wrangling

#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

# Data analysis packages
#library("TCGAbiolinks")

# Data wrangling packages     
library("openxlsx")
library("dplyr")
library("tibble")
library("stringr")

#******************************************************************************#
#                           DECLARE GLOBAL VARIABLES                           #
#******************************************************************************#

# Store path of parent directory i.e. root directory for the project
parent_path <- "C:/Users/KailasammS/Desktop/"
parent_path <- "/hpc/home/kailasamms/scratch/TCGA_GDC/"

#******************************************************************************#
#              STEP 1: DOWNLOAD TCGA EXPRESSION AND CLINICAL DATA              #
#******************************************************************************#

# There are multiple ways to download TCGA data. 
# (i) Directly download from GDC Portal (SAFEST, BEST & EASIEST) 
# (ii) TCGAbiolinks (OK as it is updated frequently)
# (iii) RTCGA (NOT RECOMMENDED as it is not updated frequently)

# This script covers ONLY (i)
# Go to https://portal.gdc.cancer.gov/ ; click "projects" and filter as below:
# "Program"               -> "TCGA" ; 
# "Data Category"         -> "sequencing reads"
# "Experimental Strategy" -> "RNA-Seq"
# "Primary Site"          -> "Bladder" (or) choose all 33 TCGA projects
# Click "Save New cohort", give a name and select your cohort.

# Now, click "Repository" and filter as below:
# "Experimental Strategy" -> "RNA-Seq"
# "Workflow Type"         -> "STAR-Counts"
# "Access"                -> "open"

# Verify that filenames are star_gene_counts.tsv
# (i) Click "Manifest" to download manifest info
# (ii) Click "Download Associated Data" and download "Sample sheet" to get info 
# linking sample name to manifest info.
# NOTE: "Sample Sheet" is NECESSARY to match the downloaded count files with 
# correct patient as the downloaded files have random names.

# (iii) Click on "11428 cases", "Clinical" and "tsv" to download clinical data 
# for all 33 TCGA projects. The clinical data will be in tar.gz format. 
# Extract the clinical.tsv file and save it as xlsx file. Discard other files.

# Next, upload the manifest txt file to "TCGA_GDC" folder in HPC cluster. 
# "TCGA_GDC" folder also has "TCGA_GDC.sh" file. Adjust the manifest file name 
# in this sh file and run it to download the count data. Once download from GDC 
# portal is complete, move all tsv files containing counts and then import
# them into R.

#******************************************************************************#
#                       STEP 2: LOAD TCGA EXPRESSION DATA                      #
#******************************************************************************#

# Read sample sheet
sample_sheet <- read.xlsx(xlsxFile=paste0(parent_path, "gdc_sample_sheet.2024-06-13.xlsx"))
# read.table(file=, header=TRUE, sep="\t", quote="", skip=0, fill=TRUE) 

# Read the count files for each Project.ID and save xlsx file
# NOTE: We use the value in "unstranded" column i.e. column 4 as counts since
# the kits used for library preparation may or may not be strand specific.
# https://www.biostars.org/p/9536035/
for (proj in unique(sample_sheet$Project.ID)){
  
  files <- sample_sheet %>% 
    dplyr::filter(Project.ID == proj) %>% 
    dplyr::select(File.Name) %>%
    unlist(use.names=FALSE)
  
  # Load one of the files to extract the genes
  temp_file <- read.table(file=paste0(parent_path, "counts/", files[1]),
                          header=TRUE,
                          sep="\t",
                          quote="",
                          fill=TRUE)
  
  df <- data.frame(ENSEMBL_ID=temp_file[,1],
                   SYMBOL=temp_file[,2])
  
  for (i in 1:length(files)){
    print(proj)
    temp_data <- read.table(file=paste0(parent_path, "counts/", files[i]),
                            header=TRUE,
                            sep="\t",
                            quote="",
                            fill=TRUE) %>%
      dplyr::select(gene_id, unstranded)
    
    colnames(temp_data) <- c("ENSEMBL_ID", files[i])
    
    df <- df %>% 
      dplyr::left_join(temp_data, by=c("ENSEMBL_ID"="ENSEMBL_ID"))
  }  
  
  # Remove rows containing counts for N_unmapped, N_multimapping, N_noFeature, etc
  df <- df[-c(1,2,3,4),]
  
  # Remove version number from Ensembl_IDs
  df$ENSEMBL_ID <- base::gsub(pattern=".[0-9]+$", replacement="", x=df$ENSEMBL_ID) 
  
  # Make sure there are no duplicate Case.ID or Sample.ID
  sample_sheet$'Case.ID' <- make.names(sample_sheet$'Case.ID', unique = TRUE)
  sample_sheet$'Sample.ID' <- make.names(sample_sheet$'Sample.ID', unique = TRUE)
  
  # Rename columns names with appropriate sample names from sample_sheet
  for (i in 1:ncol(df)){
    colnames(df)[i] <- dplyr::if_else(colnames(df)[i] %in% sample_sheet$File.Name,
                                      sample_sheet$Case.ID[which(sample_sheet$File.Name == colnames(df)[i], arr.ind=TRUE)[1]],
                                      colnames(df)[i])
  }
  
  # Save the counts data and sample sheet in xlsx file
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="raw_counts")
  openxlsx::writeData(wb, sheet="raw_counts", x=df)
  openxlsx::addWorksheet(wb, sheetName="sample_sheet")
  openxlsx::writeData(wb, sheet="sample_sheet", x=sample_sheet)
  openxlsx::saveWorkbook(wb, file=paste0(parent_path, "Read_data_", proj, ".xlsx"), 
                         overwrite=TRUE)
}

#******************************************************************************#
#                        STEP 3: LOAD TCGA CLINICAL DATA                       #
#******************************************************************************#

# https://gdc.cancer.gov/about-data/gdc-data-processing/clinical-data-standardization
# Check above link to understand clinical data columns
# prior treatment: 
# indicator related to the administration of therapeutic agents received before the body specimen was collected.
# treatment_or_therapy: 
# indicator related to the administration of therapeutic agents received.
# NOS = Not otherwise specified

# Read metadata
clinical_data <- read.xlsx(xlsxFile=paste0(parent_path, "Meta_data_TCGA_original.xlsx"))
#utils::read.table(file=,header=TRUE, sep="\t", quote="", skip=0, fill=TRUE) 

# Wrangle the data to keep important columns and reformat the metadata to be
# usable for survival analysis in future if needed.
# NOTE: We plot "days_to_follow_up" for alive patients and "days_to_death" for
# dead patients. Alive=0, Dead=1

# Identify blank entries
clinical_data[clinical_data == "'--"] <- ""

# Identify all columns that are blank for all patients and remove them
t <- apply(X=clinical_data, MARGIN=2,FUN=nchar)
clinical_data <- clinical_data[,colnames(t)[colSums(t) != 0]]

clinical_data <- original %>%
  dplyr::rename(Sample_ID=case_submitter_id,
                Project_ID=project_id,
                Sex=gender,
                Race=race,
                Ethnicity=ethnicity,
                Disease=primary_diagnosis,
                Stage=ajcc_pathologic_stage,
                T=ajcc_pathologic_t,
                N=ajcc_pathologic_n,
                M=ajcc_pathologic_m,
                Tissue.Organ=tissue_or_organ_of_origin,
                Resection.Biopsy_Site=site_of_resection_or_biopsy,
                Prior_Malignancy=prior_malignancy,
                Prior_Treatment=prior_treatment) %>%
  dplyr::group_by(Sample_ID) %>%
  dplyr::mutate(Time=dplyr::case_when(vital_status=="Alive" ~ days_to_last_follow_up,
                                      vital_status=="Dead" ~ days_to_death,
                                      TRUE ~ NA),
                Time=as.numeric(Time)/30, 
                Status=dplyr::case_when(vital_status=="Alive" ~ 0,
                                        vital_status=="Dead" ~ 1,
                                        TRUE ~ NA),
                Sex=stringr::str_to_title(Sex),
                Race=stringr::str_to_title(Race),
                Ethnicity=stringr::str_to_title(Ethnicity),
                Prior_Malignancy=stringr::str_to_title(Prior_Malignancy),
                Prior_Treatment=stringr::str_to_title(Prior_Treatment),
                Treatment=dplyr::case_when(treatment_or_therapy=="no" ~ "No",
                                           treatment_or_therapy=="yes" ~ treatment_type,
                                           treatment_or_therapy=="not reported" ~ "Not Reported",
                                           TRUE ~ "Not available"),
                Treatment = paste0(Treatment, collapse = "_"),
                Treatment=dplyr::case_when(grepl("Radiation", x=Treatment) & grepl("Pharmaceutical", x=Treatment) ~ "Radiation & Pharmaceutical",
                                           grepl("Radiation", x=Treatment) & grepl("No|Not Reported", x=Treatment) ~ "Radiation",
                                           grepl("Pharmaceutical", x=Treatment) & grepl("No|Not Reported", x=Treatment) ~ "Pharmaceutical",
                                           grepl("No", x=Treatment) & !grepl("Not Reported|Not available", x= Treatment) ~ "No Treatment",
                                           grepl("Not Reported", x= Treatment) ~ "Not Reported",
                                           TRUE ~ "Not available")) %>%
  dplyr::ungroup() %>%
  data.frame() %>%
  dplyr::distinct_at(c("Sample_ID"), .keep_all = TRUE) %>%
  dplyr::select(Project_ID, Sample_ID, Sex, Time, Status, Race, Ethnicity, 
                Stage, T, N, M, Prior_Malignancy, Prior_Treatment, Treatment, 
                Disease, Tissue.Organ, Resection.Biopsy_Site, everything(), 
                -c(treatment_or_therapy, treatment_type)) %>%
  dplyr::arrange(Project_ID, Sample_ID)

# Save clinical data in xlsx file
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName="Clinical")
openxlsx::writeData(wb, sheet="Clinical", x=clinical_data)
openxlsx::saveWorkbook(wb, file=paste0(parent_path, "Meta_data_TCGA.xlsx"), overwrite=TRUE)


# # (ii) TCGAbiolinks (OK as it is updated frequently)
# # Check GDC server status using the api https://api.gdc.cancer.gov/status
# GDC_info <- TCGAbiolinks::getGDCInfo()
# GDC_info
# 
# # Check all the available projects at TCGA
# GDC_projects <- TCGAbiolinks::getGDCprojects()
# 
# # Check project summary
# GDC_project_summary <- TCGAbiolinks:::getProjectSummary(project="TCGA-BLCA",
#                                                         legacy=FALSE)
# 
# # Download clinical data using TCGAbiolinks
# TCGAbiolinks_clinical <- TCGAbiolinks::GDCquery_clinic(project="TCGA-BLCA",
#                                                        type="clinical",
#                                                        save.csv=TRUE)
# 
# # (iii) RTCGA (NOT RECOMMENDED as it is not updated frequently) 
# # Download clinical data using RTCGA
# RTCGA_clinical <- RTCGA::survivalTCGA(BLCA.clinical,
#                                       extract.cols="admin.disease_code",
#                                       extract.names=FALSE,
#                                       barcode.name="patient.bcr_patient_barcode",
#                                       event.name="patient.vital_status",
#                                       days.to.followup.name="patient.days_to_last_followup",
#                                       days.to.death.name="patient.days_to_death")
#
# # You will notice the clinical data is very different. Two main differences:
# # (i) "RTCGA_clinical$patient.vital_status" is same as
# # "TCGAbiolinks_clinical$year_of_death". Clearly, RTCGA is wrong because you can
# # notice some patients are dead based on "TCGAbiolonks_clinical$vital_status"
# # (ii) "RTCGA_clinical$times" is combination of
# # "TCGAbiolonks_clinical$days_to_last_follow_up" and "TCGAbiolonks_clinical$days_to_death".

#******************************************************************************#
#        STEP 4: NORMALIZE PAN CANCER DATA & CORRECT FOR BATCH EFFECTS         #
#******************************************************************************#

# Download batch corrected normalized pan cancer RNA expression data from below
# https://gdc.cancer.gov/about-data/publications/pancanatlas
# RNA (Final) - EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv
# The clinical data from previous step can be used 

parent_path <- "/hpc/home/kailasamms/scratch/TCGA_GDC/"
meta_data <- openxlsx::read.xlsx(xlsxFile=paste0(parent_path, "Meta_data_TCGA.xlsx"))
#read_data <- utils::read.table(file=paste0(parent_path, "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"),
read_data <- utils::read.table(file=paste0(parent_path, "Normalized_Counts_TCGA_original.tsv"),
                               header=TRUE, sep="\t", quote="", skip=0, fill=TRUE) 

# Rename column names to match with Sample_ID from meta_data
cols <- colnames(read_data)

# NOTE: str_detect is better than grepl
# grepl(pattern, replacement) where both pattern and replacement MUST be string
# str_detect(string, pattern) either string or pattern can be a vector but not BOTH.

for (i in 1:ncol(read_data)){
  
  if(sum(stringr::str_detect(string=colnames(read_data)[i],
                             pattern=make.names(meta_data$Sample_ID))) == 1){
    cols[i] <- meta_data$Sample_ID[stringr::str_detect(string=colnames(read_data)[i],
                                                       pattern=make.names(meta_data$Sample_ID))]
    print(i)
  }
}

# Replace with proper Sample_ID as column names
colnames(read_data) <- cols
colnames(read_data)[1] <- "ENTREZ_ID"

# There are 795 duplicated columns
sum(duplicated(colnames(read_data)))
which(colnames(read_data) == "TCGA-06-0125")
read_data[1:5, 6695]
read_data[1:5, 6696]

# Make them unique
colnames(read_data) <- make.names(colnames(read_data), unique=TRUE)

# Separate the first column which contains SYMBOL and ENTREZ_ID joined together
# using piping character | into 2 columns
read_data <- read_data %>% 
  tidyr::separate(col=ENTREZ_ID, into=c("SYMBOL", "ENTREZ_ID"), sep="\\s*\\|\\s*")

# Remove the \" from SYMBOL and ENTREZ_ID
# cat(read_data$ENTREZ_ID[100]) will show 5826"
# print(read_data$ENTREZ_ID[100]) will show "5826\""
read_data$SYMBOL <- gsub(pattern="[\"\"]", replacement="", x=read_data$SYMBOL)
read_data$ENTREZ_ID <- gsub(pattern="[\"\"]", replacement="", x=read_data$ENTREZ_ID)

# Replace all NA values with 0
read_data <- read_data %>%
  base::replace(is.na(.), 0)

# Remove genes with no expression in all samples
read_data <- read_data[rowSums(read_data[,c(-1,-2)]) !=0,]

# Replace symbols with ? with ENTREZ_ID
read_data <- read_data %>% 
  dplyr::mutate(SYMBOL = dplyr::case_when(SYMBOL == "?" ~ ENTREZ_ID, TRUE ~ SYMBOL))

# Remove duplicated genes
read_data <- read_data %>% dplyr::distinct_at("SYMBOL", .keep_all = TRUE)

# Save the reformatted read_data(DO NOT SAVE xlsx as file is too large)
write.table(x=read_data, file=paste0(parent_path, "Normalized_Counts_TCGA.tsv"), 
            quote=FALSE, sep='\t', row.names = FALSE)

#******************************************************************************#
#              STEP 4A: NORMALIZE EACH CANCER DATA INDIVIDUALLY                #
#******************************************************************************#

parent_path <- "/hpc/home/kailasamms/scratch/TCGA_GDC/"
results_path <- parent_path
species <- "Homo sapiens"
Variable <- "Treatment"
Comparisons <- list(Target=c("Radiation"),
                    Reference=c("No Treatment"))
padj.cutoff <- 0.1
lfc.cutoff <- 0  
annotations <- get_annotations(species)

meta_data_full <- openxlsx::read.xlsx(xlsxFile=paste0(parent_path, "Meta_data_TCGA.xlsx"))

for(proj in unique(meta_data_full$Project_ID)){
  
  read_data <- openxlsx::read.xlsx(xlsxFile=paste0(parent_path, "Read_data_", proj, ".xlsx")) %>% 
    dplyr::mutate(SYMBOL = make.names(SYMBOL, unique=TRUE))
  
  meta_data <- meta_data_full %>% 
    dplyr::filter(make.names(Sample_ID) %in% colnames(read_data)) %>%
    dplyr::mutate(Batch=Project_ID)
  
  # Save the cancer specific meta_data
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Metadata")
  openxlsx::writeData(wb, sheet="Metadata", x=meta_data)
  openxlsx::saveWorkbook(wb, file=paste0(parent_path, "Meta_data_", proj, ".xlsx"), overwrite=TRUE)
  
  meta_data <- prep_metadata(meta_data, Variable)
  read_data <- prep_readdata(read_data, meta_data)
  l <- check_data(read_data, meta_data)
  meta_data <- l[[2]]
  read_data <- l[[1]]
  sva_dds <- svaseq_batch(read_data, meta_data)
  n <- 1
  
  #Perform DESeq2() using in-built batch modelling
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
  dds <- run_deseq2(dds, meta_data, annotations, Comparisons, n, lfc.cutoff, padj.cutoff, approach)
  deseq2_norm_counts(dds, annotations, approach) # batch corrected if you more than 1 batch
  plot_qc(dds, meta_data, approach)
  
  # Perform DESeq2() using sva modelled surrogate variables SV1 and SV2
  approach <- "sva_modelled"
  sva_dds <- run_deseq2(sva_dds, meta_data, annotations, Comparisons, n, lfc.cutoff, padj.cutoff, approach)
  # calc_norm_counts(sva_dds, annotations, approach)   # uncorrected
  svaseq_norm_counts(sva_dds, annotations, approach)   # sva seq batch corrected
  plot_qc(sva_dds, meta_data, approach)
}


#*****************************************************************************#



# The below 2 posts recommend using all samples from a single experiment for 
# normalizing but avoiding using all samples from multiple experiments.
# https://support.bioconductor.org/p/92879/, https://www.biostars.org/p/9560478/
# https://support.bioconductor.org/p/59711/
# https://support.bioconductor.org/p/98765/
# https://www.biostars.org/p/336298/
# https://support.bioconductor.org/p/75260/

# LINK #1: In general, it is recommended to use all the available samples from 
# an experiment in the analysis, even if you are not interested in differential 
# expression for some of those samples. The main reason for this is that having
# more samples allows more accurate estimation of the gene dispersion values.

# LINK #2: If you are not going to test for differential expression between
# experiments, then there is no purpose in normalizing them together. The more
# worrying problem with analyzing all your experiments as a single data set is 
# that a single dispersion value will be estimated for each gene across all 
# experiments. This is only ok if you believe that every gene has equal 
# biological variability in all your experiments, which is unlikely to be case.

# Since, it is not appropriate to merge all samples from 33 different TCGA 
# projects and normalize them together due to issues explained above, we will
# perform normalization on each TCGA project individually ASSUMING there are
# no batch effects within each experiment.