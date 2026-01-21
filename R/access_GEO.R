#!/usr/bin/env Rscript



output_path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/"

#******************************************************************************#
#                           PREPARE sample_info.txt                            #
#******************************************************************************#

# Step 1: Go to GEO (https://www.ncbi.nlm.nih.gov/geo/)
# Step 2: Search using the GEO accession (e.g. GSE169379)
# Step 3: Scroll down to the “Samples” section
# Step 4: Copy the sample information into "sample_info.txt" in this format:

# GSM5199001	B1246-GEX: MIBC_rxn1246
# GSM5199003	B1246-HTO: MIBC_rxn1246
# GSM5199005	B1247-GEX: MIBC_rxn1247
# GSM5199007	B1247-HTO: MIBC_rxn1247
# ...and so on

# Delete any rows for samples you do not wish to download.

sample_info <- utils::read.table(file = paste0(output_path, "sample_info.txt"),
                                 header = FALSE,
                                 sep = "\t") %>%
  dplyr::rename("GEO_Accession" = V1,
                "Description" = V2)

#******************************************************************************#
#                           PREPARE SraRunTable.csv                            #
#******************************************************************************#

# Step 5: Download SRA metadata via the "SRA Run Selector"
# Step 6: Save the CSV as "SraRunTable.csv" in your working directory

meta_data <- utils::read.table(file = paste0(output_path, "SraRunTable.csv"),
                               header = TRUE,
                               sep = ",") %>%
  dplyr::rename_with(~ gsub("\\.", "_", .x), everything())

#******************************************************************************#
#                     MERGE sample_info.txt & Metadata.xlsx                    #
#******************************************************************************#

# Identify which column in meta_data contains the GEO accession numbers
match_col <- colnames(meta_data)[
  as.numeric(which(meta_data == sample_info$GEO_Accession[1], arr.ind = TRUE)[1, 2])]

# Merge sample_info with SRA metadata
meta_merged <- sample_info %>%
  dplyr::left_join(meta_data, by = c("GEO_Accession" = match_col)) %>%
  dplyr::select(Run, Description, everything())

#******************************************************************************#
#                     SAVE Metadata.xlsx & SRR_Acc_List.txt                    #
#******************************************************************************#

# Save the merged metadata as Excel
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Metadata")
openxlsx::writeData(wb, sheet = "Metadata", x = meta_merged, rowNames = FALSE)
openxlsx::saveWorkbook(wb, file = paste0(output_path, "Metadata_temp.xlsx"), overwrite = TRUE)

# Save SRR Run IDs to text file
srr <- meta_merged %>% dplyr::select(Run)

utils::write.table(x = srr,
                   file = paste0(output_path, "SRR_Acc_List.txt"),
                   quote = FALSE,
                   col.names = FALSE,
                   row.names = FALSE)

#******************************************************************************#
#                                  IMPORTANT                                   #
#******************************************************************************#

# After saving "SRR_Acc_List.txt":
# 1. Open it in Notepad++.
# 2. Change line endings from "Windows (CR LF)" → "Unix (LF)".
# 3. Upload to HPC cluster to use with 01_Fastqer_dump.sh.
# Failure to do this may result in only the last value being read in Linux arrays.
