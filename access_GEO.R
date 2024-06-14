#!/usr/bin/env Rscript

#******************************************************************************#
#                          INSTALL NECESSARY PACKAGES                          #
#******************************************************************************#

# # Data wrangling packages
# install.packages("openxlsx")            # Needed for reading, writing xlsx files
# install.packages("dplyr")               # Needed for data wrangling
# install.packages("tibble")              # Needed for advanced data wrangling
# install.packages("stringr")             # Needed for advanced data wrangling

#******************************************************************************#
#                           LOAD NECESSARY PACKAGES                            #
#******************************************************************************#

# Data wrangling packages     
library("openxlsx")
library("dplyr")
library("tibble")
library("stringr")

#******************************************************************************#
#                                  IMPORT DATA                                 #
#******************************************************************************#

# NOTE: SRA Run Selector displays a table with several columns
# SRR ids              : Corresponds to each fastq file
# SAMN, SRX, GSM ids   : Corresponds to each sample
# We are interested ONLY in getting the SRR ids (from SRA Run Selector) and 
# corresponding sample description (from sample_info.txt).

# Go to GEO (https://www.ncbi.nlm.nih.gov/geo/). 
# Search using GEO Accesssion (Eg: GSE169379)
# Once GSE169379 loads, scroll down and copy-paste sample info into 
# "sample_info.txt" file as per format below:

# GSM5199001	B1246-GEX: MIBC_rxn1246
# GSM5199003	B1246-HTO: MIBC_rxn1246
# .......     ........

results_path <- "C:/Users/KailasammS/Desktop/"

# Name the file "sample_info.txt" and import it.
# Delete samples you don't want to download from this file.
sample_info <- utils::read.table(file = paste0(results_path, "sample_info.txt"), 
                                 header = FALSE,
                                 sep = "\t") %>%
  dplyr::rename("GEO_Accession__exp_" = V1, "Description" = V2)

# Click SRA Run Selector at bottom of page, then download metadata and import it.
meta_data <- utils::read.table(file = paste0(results_path, "SraRunTable.txt"), 
                               header = TRUE,
                               sep = ",") %>%
  dplyr::rename_with(.fn = ~gsub(pattern = "\\.", replacement = "_", x = .x), .cols = everything()) 
  #dplyr::select(GEO_Accession__exp_, Run)

# Merge sample_info and meta_data
meta_data <- sample_info %>% 
  dplyr::left_join(meta_data, by=c("GEO_Accession__exp_"= "GEO_Accession__exp_")) %>%
  dplyr::select(Run, Description, everything())

# Save meta_data
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Metadata")
openxlsx::writeData(wb, sheet = "Metadata", x = meta_data, rowNames = FALSE)
openxlsx::saveWorkbook(wb, 
                       file = paste0(results_path, "Metadata.xlsx"), 
                       overwrite = TRUE)

#******************************************************************************#
#                            CHOOSE SRR OF INTEREST                            #
#******************************************************************************#

# Choose runs of interest and save it to "SRR_Acc_List.txt"
srr <- meta_data %>% 
  dplyr::select(Run)

# Save to a txt file.
utils::write.table(x = srr, 
                   file = "C:/Users/KailasammS/Desktop/SRR_Acc_List.txt",
                   quote = FALSE,
                   col.names = FALSE,
                   row.names = FALSE)

#******************************************************************************#
# MOST IMPORTANT: Open "SRR_Acc_List.txt" in Notepad++.                        #
# Right click "Windows (CR LF)" at right bottom and change to "Unix (LF)".     #
# Next, upload "SRR_Acc_List.txt" to HPC cluster and use 02a_Fastq_dump.sh to  #
# download the SRR files to HPC cluster. If you do not change this, arrays     #  
# created in Linux by reading the txt file wont work properly and only the     #
# last value will be stored in the array.                                      #
#******************************************************************************#
