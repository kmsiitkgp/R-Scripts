library(readr)
library(dplyr)
library(openxlsx)

# Number of CpG sites per region ( all 3 files MUST be identical)
count_hmc <- read_tsv("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/biomodal reports/count_hmc.tsv")
count_mc <- read_tsv("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/biomodal reports/count_mc.tsv")
count_total_c <- read_tsv("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/biomodal reports/count_total_c.tsv")

all.equal(as.data.frame(count_hmc),  as.data.frame(count_mc), check.attributes = FALSE)
all.equal(as.data.frame(count_hmc),  as.data.frame(count_total_c),  check.attributes = FALSE)


frac_hmc <- read_tsv("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/biomodal reports/frac_hmc.tsv")
frac_mc <- read_tsv("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/biomodal reports/frac_mc.tsv")

# Number of CpG sites * reads per region
sum_hmc <- read_tsv("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/biomodal reports/sum_hmc.tsv")
sum_mc <- read_tsv("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/biomodal reports/sum_mc.tsv")
sum_total_c <- read_tsv("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/biomodal reports/sum_total_c.tsv")
frac_hmc <- read_tsv("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/biomodal reports/frac_hmc.tsv")
frac_mc <- read_tsv("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/biomodal reports/frac_mc.tsv")

sum(round(sum_hmc[4:27]/sum_total_c[4:27] - frac_hmc[4:27],2), na.rm=TRUE)
sum(round(sum_mc[4:27]/sum_total_c[4:27] - frac_mc[4:27],2), na.rm=TRUE)

# Filter out weak signals
dmr_hmc <- read_tsv("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/biomodal reports/DMR_20251124_130348_DMR_hmc_Pre__Post_20251124_130348.tsv")
dmr_mc <- read_tsv("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/biomodal reports/DMR_20251124_130348_DMR_mc_Pre__Post_20251124_130348.tsv")

sig_dmr_hmc <- dmr_hmc %>% 
  dplyr::mutate(mean_mod_group_1 = ifelse(mean_mod_group_1 == 0, 1e-6, mean_mod_group_1),
                mean_mod_group_2 = ifelse(mean_mod_group_2 == 0, 1e-6, mean_mod_group_2),
                mod_fold_change = mean_mod_group_2 / mean_mod_group_1,
                mod_fold_change = log2(mod_fold_change)) %>%
  dplyr::filter(dmr_qvalue <= 0.05)
  #dplyr::filter(dmr_qvalue <= 0.05, !is.na(Name), abs(mod_fold_change) >= 0.58)

sig_dmr_mc <- dmr_mc %>% 
  dplyr::mutate(mean_mod_group_1 = ifelse(mean_mod_group_1 == 0, 1e-6, mean_mod_group_1),
                mean_mod_group_2 = ifelse(mean_mod_group_2 == 0, 1e-6, mean_mod_group_2),
                mod_fold_change = mean_mod_group_2 / mean_mod_group_1,
                mod_fold_change = log2(mod_fold_change)) %>% 
  dplyr::filter(dmr_qvalue <= 0.05)
  #dplyr::filter(dmr_qvalue <= 0.05, !is.na(Name), abs(mod_fold_change) >= 0.58)

# Create a new workbook
wb <- createWorkbook()
addWorksheet(wb, "Original_5mc")
writeData(wb, sheet = "Original_5mc", dmr_mc)
addWorksheet(wb, "Processed_5mc")
writeData(wb, sheet = "Processed_5mc", sig_dmr_mc)
addWorksheet(wb, "Original_5hmc")
writeData(wb, sheet = "Original_5hmc", dmr_hmc)
addWorksheet(wb, "Processed_5hmc")
writeData(wb, sheet = "Processed_5hmc", sig_dmr_hmc)
saveWorkbook(wb, "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/biomodal reports/Biomodal_results.xlsx", overwrite = TRUE)



library(GenomicRanges)
library(rtracklayer)

# Load normal prostate 5hmC peaks
normal_5hmC_file <- "/mnt/data/GSM4290243_prostate1_peaks.narrowPeak"
normal_peaks <- import(normal_5hmC_file, format = "narrowPeak")

# Load your cfDNA peaks (replace with your file path and format)
cfDNA_peaks <- import("your_cfDNA_peaks.bed", format = "bed")

# Find overlaps between cfDNA and normal peaks
overlaps <- findOverlaps(cfDNA_peaks, normal_peaks)

# Get cfDNA peaks NOT overlapping normal 5hmC
cfDNA_unique <- cfDNA_peaks[-queryHits(overlaps)]

# cfDNA_unique now contains peaks unique (non-overlapping) to your cfDNA samples

