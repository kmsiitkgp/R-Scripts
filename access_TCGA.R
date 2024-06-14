#!/usr/bin/env Rscript

#******************************************************************************#
#                          INSTALL NECESSARY PACKAGES                          #
#******************************************************************************#

# if (!require("BiocManager", quietly = TRUE)) 
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
library("TCGAbiolinks")

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
# Click "Manifest" to download manifest info
# Click "Download Associated Data" and download "Sample sheet" to get info linking sample name to manifest info 
# Click "Download Associated Data" and download "Metadata"

# Next, upload the manifest txt file to "TCGA_GDC" folder in HPC cluster. 
# "TCGA_GDC" folder also has "TCGA_GDC.sh" file. Adjust the manifest file name 
# in this sh file and run it to download the count data. Once download is 
# complete, download the count folder to desktop

# NOTE: "Sample Sheet" is NECESSARY to match the downloaded count files with 
# correct patient as the downloaded files have random names.

# NOTE: The downloaded clinical file will be in tar.gz format. 
# Extract using Peazip or other softwares and you will see clinical.tsv file.
# Discard rest of extracted files.

#******************************************************************************#
#                        STEP 2: LOAD TCGA CLINICAL DATA                       #
#******************************************************************************#

# Import clinical data from clinical.tsv file.
GDC_clinical <- utils::read.table(file = paste0(parent_path, "clinical.tsv"),
                                  header = TRUE,
                                  sep = "\t",
                                  quote="",
                                  skip = 0,
                                  fill = TRUE) %>%
  dplyr::distinct(case_id, .keep_all=TRUE)


# Wrangle the data to keep interesting columns and format the metadata to be
# usable for survival analysis in future if needed.
# NOTE: Manually format clinical data as columns vary for each dataset.
# NOTE: We plot "days_to_follow_up" for alive patients and "days_to_death" for
# dead patients. Alive = 0, Dead = 1

GDC_clinical <- GDC_clinical %>%
  dplyr::rename(Sample_id = case_submitter_id,
                AJCC_stage = ajcc_pathologic_stage,
                Race = race,
                Disease = primary_diagnosis) %>%
  dplyr::mutate(Time = dplyr::if_else(vital_status=="Alive", days_to_last_follow_up,
                                      dplyr::if_else(vital_status=="Dead", days_to_death, "NA")),
                Status = dplyr::if_else(vital_status=="Alive", 0,
                                        dplyr::if_else(vital_status=="Dead", 1, 1000)),
                Sex = dplyr::if_else(gender=="male", "Male",
                                     dplyr::if_else(gender=="female", "Female", "NA"))) %>%
  dplyr::select(Sample_id, Sex, Time, Status, vital_status, Race,
                AJCC_stage, Disease)

# Save clinical data in xlsx file
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Clinical")
openxlsx::writeData(wb, sheet = "Clinical", x = GDC_clinical)
openxlsx::saveWorkbook(wb,
                       file = paste0(parent_path, "Metadata.xlsx"),
                       overwrite = TRUE)


# # (ii) TCGAbiolinks (OK as it is updated frequently)
# # Check GDC server status using the api https://api.gdc.cancer.gov/status
# GDC_info <- TCGAbiolinks::getGDCInfo()
# GDC_info
# 
# # Check all the available projects at TCGA
# GDC_projects <- TCGAbiolinks::getGDCprojects()
# 
# # Check project summary
# GDC_project_summary <- TCGAbiolinks:::getProjectSummary(project = "TCGA-BLCA",
#                                                         legacy = FALSE)
# 
# # Download clinical data using TCGAbiolinks
# TCGAbiolinks_clinical <- TCGAbiolinks::GDCquery_clinic(project = "TCGA-BLCA",
#                                                        type = "clinical",
#                                                        save.csv = TRUE)
# 
# # (iii) RTCGA (NOT RECOMMENDED as it is not updated frequently) 
# # Download clinical data using RTCGA
# RTCGA_clinical <- RTCGA::survivalTCGA(BLCA.clinical,
#                                       extract.cols="admin.disease_code",
#                                       extract.names = FALSE,
#                                       barcode.name = "patient.bcr_patient_barcode",
#                                       event.name = "patient.vital_status",
#                                       days.to.followup.name = "patient.days_to_last_followup",
#                                       days.to.death.name = "patient.days_to_death")
#
# # You will notice the clinical data is very different. Two main differences:
# # (i) "RTCGA_clinical$patient.vital_status" is same as
# # "TCGAbiolinks_clinical$year_of_death". Clearly, RTCGA is wrong because you can
# # notice some patients are dead based on "TCGAbiolonks_clinical$vital_status"
# # (ii) "RTCGA_clinical$times" is combination of
# # "TCGAbiolonks_clinical$days_to_last_follow_up" and "TCGAbiolonks_clinical$days_to_death".

#******************************************************************************#
#                       STEP 2: LOAD TCGA EXPRESSION DATA                      #
#******************************************************************************#

# Read sample sheet and keep ONLY tumor samples (VERY IMPORTANT)
# You will notice some count files belong to normal tissues
sample_sheet <- read.table(file = paste0(parent_path, "gdc_sample_sheet.2023-04-17.tsv"),
                           header = TRUE,
                           sep = "\t",
                           quote="",
                           skip = 0,
                           fill = TRUE) %>%
  dplyr::filter(Sample.Type == "Primary Tumor")

# Read all count files
# List all files in count folder
files <- list.files(path = paste0(parent_path, "counts/"))

# Load 1 file to determine row position where the Ensembl genes start
temp_file <- read.table(file = paste0(parent_path, "counts/", files[1]),
                        header = TRUE,
                        sep = "\t",
                        quote="",
                        fill = TRUE)

#Find the row position where the Ensembl genes start
junk_rows <- which(grepl(pattern = "^EN", x = temp_file[,1]), arr.ind = TRUE)[1]
count_column <- which(colnames(temp_file) == "unstranded", arr.ind = TRUE)[1]

# Read the 1st count file. 
# NOTE: The top 6 lines of count files have junk data. So, skip them.
# NOTE: You need to adjust the value of skip parameter based on your data
# NOTE: STAR has counts for unstranded, stranded_first and stranded_second
# Decide which column to use.
# NOTE: We use the value in "unstranded" column i.e. column 4 as counts as the kits 
# used for library preparation may or may not be strand specific.

data <- read.table(file = paste0(parent_path, "counts/", files[1]),
                   header = TRUE,
                   sep = "\t",
                   quote="",
                   fill = TRUE)
data <- data[c(junk_rows:nrow(data)),c(1,count_column)]
colnames(data) <- c("Ensembl_ID", files[1])

# Append rest of count files to the previous count file
for (i in 2:length(files)){
  temp_data <- read.table(file = paste0(parent_path, "counts/", files[i]),
                          header = TRUE,
                          sep = "\t",
                          quote="",
                          fill = TRUE)
  temp_data <- temp_data[c(junk_rows:nrow(data)),c(1,count_column)]
  colnames(temp_data) <- c("Ensembl_ID", files[i])
  
  data <- data %>% dplyr::left_join(temp_data, by=c("Ensembl_ID" = "Ensembl_ID"))
}

# Remove version number from Ensembl_IDs
data$Ensembl_ID <- base::gsub(pattern = ".[0-9]+$", replacement = "", x = data$Ensembl_ID) 

# Rename columns names with appropriate sample names from sample_sheet
for (i in 2:ncol(data)){
  colnames(data)[i] <- dplyr::if_else(colnames(data)[i] %in% sample_sheet$File.Name,
                                      sample_sheet$Case.ID[which(sample_sheet$File.Name == colnames(data)[i], arr.ind = TRUE)[1]],
                                      colnames(data)[i])
}

# Remove samples without clincal data
data <- data %>% 
  dplyr::select("Ensembl_ID", GDC_clinical$Sample_id)

# Save the counts data in xlsx file
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "raw_counts")
openxlsx::writeData(wb, sheet = "raw_counts", x = data)
openxlsx::saveWorkbook(wb, 
                       file = paste0(parent_path, "Readdata.xlsx"), 
                       overwrite = TRUE)

#******************************************************************************#