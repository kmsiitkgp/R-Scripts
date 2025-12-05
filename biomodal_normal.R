# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/

# ---- 0️⃣ Setup ----
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(openxlsx)

# ---- Input files ----
path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/biomodal reports"
gtf_file <- file.path(path, "hg19.refGene.gtf")         # RefSeq GTF
peak_files <- c(file.path(path, "prostate1_peaks.narrowPeak"), 
                file.path(path, "prostate2_peaks.narrowPeak"),
                file.path(path, "prostate3_peaks.narrowPeak"),
                file.path(path, "prostate4_peaks.narrowPeak"),
                file.path(path, "prostate5_peaks.narrowPeak"))

# ---- 1️⃣ Import GTF ----
gtf <- rtracklayer::import(con = gtf_file, format = "gtf")


# ---- 2️⃣ Define regions to match with Biomodal  ----

transcripts <- gtf[GenomicRanges::mcols(gtf)$type == "transcript"]

# GRange object for gene body
gene_body_gr <- transcripts

# GRange object for promoters
promoter_gr <- GenomicRanges::promoters(x = transcripts, 
                                        upstream = 1000, 
                                        downstream = 0)

# GRange object around TSS
tss_region_gr <- GenomicRanges::promoters(x = transcripts, 
                                          upstream = 200, 
                                          downstream = 200)

# # GRange object at TSS
# tss_pos <- dplyr::if_else(condition = as.logical(GenomicRanges::strand(transcripts) == "+"),
#                           true = GenomicRanges::start(transcripts),
#                           false = GenomicRanges::end(transcripts))
# tss_gr <- GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(transcripts),
#                                  ranges = IRanges::IRanges(start = tss_pos, end = tss_pos),
#                                  strand = GenomicRanges::strand(transcripts),
#                                  gene_id = transcripts$gene_id)

# ---- 3️⃣ Annotate peaks for each file ----

# Helper function to get collapsed gene names for a single Hits object
get_collapsed_genes <- function(peaks, hits, subject_gr) {
  # We use seq_along(peaks) here, but define peaks inside lapply
  sapply(seq_along(peaks), function(i) {
    # Extract gene names for the current peak 'i'
    unique_names <- unique(subject_gr$gene_name[S4Vectors::subjectHits(hits)[S4Vectors::queryHits(hits) == i]])
    
    # Return "NA" if no names are found, else return semicolon-separated string
    if (length(unique_names) == 0) {
      return("NA")
    } else {
      return(paste(unique_names, collapse = ";"))
    }
  })
}


annotated_list <- lapply(peak_files, function(peak_file) {
  
  # Read the peak file
  peaks <- rtracklayer::import(con = peak_file, format = "narrowPeak")
  
  # Overlaps
  promoter_hits   <- GenomicRanges::findOverlaps(query = peaks, subject = promoter_gr)
  tss_hits        <- GenomicRanges::findOverlaps(query = peaks, subject = tss_region_gr)
  gene_body_hits <- GenomicRanges::findOverlaps(query = peaks, subject = gene_body_gr)
  
  # Build annotation table
  peak_df <- data.frame(peak_name     = GenomicRanges::mcols(peaks)$name,
                        peak_chr      = GenomicRanges::seqnames(peaks),
                        peak_start    = GenomicRanges::start(peaks),
                        peak_end      = GenomicRanges::end(peaks),
                        in_promoter   = GenomicRanges::countOverlaps(query = peaks, subject = promoter_gr) > 0,
                        in_tss        = GenomicRanges::countOverlaps(query = peaks, subject = tss_region_gr) > 0,
                        in_gene_body  = GenomicRanges::countOverlaps(query = peaks, subject = gene_body_gr) > 0, 
                        promoter_gene_names  = get_collapsed_genes(peaks, promoter_hits, promoter_gr),
                        tss_gene_names       = get_collapsed_genes(peaks, tss_hits, tss_region_gr),
                        gene_body_gene_names = get_collapsed_genes(peaks, gene_body_hits, gene_body_gr)
                        )
  return(peak_df)
})

# ---- 4️⃣ Merge genes across samples ----

# Get all gene IDs in promoter/TSS/gene body for each sample
promoter_genes <- lapply(annotated_list, function(df) {
  df %>%
    dplyr::filter(in_promoter & promoter_gene_names != "NA") %>%
    tidyr::separate_longer_delim(cols = promoter_gene_names, delim = ";") %>%
    dplyr::pull(promoter_gene_names) %>%
    unique()
})

tss_genes <- lapply(annotated_list, function(df) {
  df %>%
    dplyr::filter(in_tss & tss_gene_names != "NA") %>%
    tidyr::separate_longer_delim(cols = tss_gene_names, delim = ";") %>%
    dplyr::pull(tss_gene_names) %>%
    unique()
})

gene_body_genes <- lapply(annotated_list, function(df) {
  df %>%
    dplyr::filter(in_gene_body & gene_body_gene_names != "NA") %>%
    tidyr::separate_longer_delim(cols = gene_body_gene_names, delim = ";") %>%
    dplyr::pull(gene_body_gene_names) %>%
    unique()
})

# Genes present in ALL 5 samples
common_promoter_genes <- purrr::reduce(.x = promoter_genes, .f = intersect)
common_tss_genes      <- purrr::reduce(.x = tss_genes, .f = intersect)
common_gene_body_genes <- purrr::reduce(.x = gene_body_genes, .f = intersect)

# ---- 5️⃣ Save results ----

# Create data frames for each region with a column indicating the region
# Using Name, Annotation, Promoter, TSS, Gene to match biomodal DMR xlsx
promoter_df   <- data.frame(Name = common_promoter_genes, Annotation = "Promoter")
tss_df        <- data.frame(Name = common_tss_genes, Annotation = "TSS")
gene_body_df  <- data.frame(Name = common_gene_body_genes, Annotation = "Gene")

# Merge all into one data frame
all_genes_df <- bind_rows(promoter_df, tss_df, gene_body_df) %>%
  distinct()  # Remove duplicates if a gene is in multiple regions

# Write to Excel
wb <- createWorkbook()
addWorksheet(wb, "Genes_by_region")
writeData(wb, sheet = "Genes_by_region", all_genes_df)
saveWorkbook(wb, file.path(path, "all_genes_by_region.xlsx"), overwrite = TRUE)

# ************************** #


# Read Biomodal results
biomodal_5hmc <- read.xlsx(file.path(path, "Biomodal_results.xlsx"), sheet = "Processed_5hmc")
biomodal_5mc <- read.xlsx(file.path(path, "Biomodal_results.xlsx"), sheet = "Processed_5mc")
biomodal_5hmc_raw <- read.xlsx(file.path(path, "Biomodal_results.xlsx"), sheet = "Original_5hmc")
biomodal_5mc_raw <- read.xlsx(file.path(path, "Biomodal_results.xlsx"), sheet = "Original_5mc")
normal_peaks <- read.xlsx(file.path(path,  "all_genes_by_region.xlsx"))

# Remove normal genes
biomodal_tumor <- biomodal_5hmc %>%
  anti_join(normal_peaks, by = c("Name"="Name", "Annotation"="Annotation")) %>%
  arrange(Annotation, Name)



# Get summary
# Create the summary table
dmr_summary <- biomodal_tumor %>%
  # Filter for significant DMRs
  dplyr::filter(dmr_qvalue < 0.05) %>%
  # Create a column to categorize the change direction
  dplyr::mutate(Direction = dplyr::case_when(mod_difference > 0 ~ "Hyper-modified (Up)",
                                             mod_difference < 0 ~ "Hypo-modified (Down)",
                                             TRUE ~ "No Change")) %>% # Should be rare after filtering
  # Count unique genes by Annotation and Direction
  dplyr::group_by(Annotation, Direction) %>%
  dplyr::summarise(Unique_Genes = n_distinct(Name),
                   .groups = "drop") %>%
  # Pivot wider for a clean summary table format
  tidyr::pivot_wider(names_from = Direction,
                     values_from = Unique_Genes,
                     values_fill = 0)

# Display the summary table
print(dmr_summary)


#  Comparison directionof 5mc and 5hmc
comparison_df <- dplyr::inner_join(x = biomodal_5hmc %>% dplyr::rename(hmc_diff = mod_difference, hmc_logFC = mod_fold_change) %>% dplyr::select(Name, Annotation, hmc_diff, hmc_logFC ),
                                   y = biomodal_5mc %>% dplyr::rename(mc_diff = mod_difference, mc_logFC = mod_fold_change) %>% dplyr::select(Name, Annotation, mc_diff, mc_logFC ),
                                   by = c("Name"= "Name", "Annotation"="Annotation")) %>%
  dplyr::mutate(Trend = dplyr::case_when((mc_diff > 0 & hmc_diff < 0) | (mc_diff < 0 & hmc_diff > 0) ~ "Inverse",
                                         (mc_diff > 0 & hmc_diff > 0) | (mc_diff < 0 & hmc_diff < 0) ~ "Concordant",
                                         TRUE ~ "Unknown")) %>%
  dplyr::arrange(Trend, Annotation, Name) %>%
  dplyr::select(Name, Annotation, hmc_diff, mc_diff, hmc_logFC, mc_logFC, Trend)

comparison_df %>% dplyr::count(Trend)

# Write to Excel
wb <- createWorkbook()
addWorksheet(wb, "5hmc_Tumor")
writeData(wb, sheet = "5hmc_Tumor", biomodal_tumor)
addWorksheet(wb, "5hmc.5mc.Trend")
writeData(wb, sheet = "5hmc.5mc.Trend", comparison_df)
addWorksheet(wb, "5hmc_Processed")
writeData(wb, sheet = "5hmc_Processed", biomodal_5hmc)
addWorksheet(wb, "5mc_Processed")
writeData(wb, sheet = "5mc_Processed", biomodal_5mc)
addWorksheet(wb, "5hmc_Raw")
writeData(wb, sheet = "5hmc_Raw", biomodal_5hmc_raw)
addWorksheet(wb, "5mc_Raw")
writeData(wb, sheet = "5mc_Raw", biomodal_5mc_raw)
saveWorkbook(wb, file.path(path, "Biomodal_tumor.xlsx"), overwrite = TRUE)



