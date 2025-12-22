library(readr)
library(dplyr)
library(openxlsx)

path <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Biomodal"

# ---- Cross check biomodal outputs for accuracy ----

# (i) Number of CpG sites per region ( all 3 files MUST be identical)
count_hmc <- readr::read_tsv(file.path(path, "biomodal_output", "count_hmc.tsv"), show_col_types = FALSE)
count_mc <- readr::read_tsv(file.path(path, "biomodal_output", "count_mc.tsv"), show_col_types = FALSE)
count_total_c <- readr::read_tsv(file.path(path, "biomodal_output", "count_total_c.tsv"), show_col_types = FALSE)

all.equal(as.data.frame(count_hmc),  as.data.frame(count_mc), check.attributes = FALSE)
all.equal(as.data.frame(count_hmc),  as.data.frame(count_total_c),  check.attributes = FALSE)

# (ii) Fraction of 5mc and 5hmc
frac_hmc <- readr::read_tsv(file.path(path, "biomodal_output", "frac_hmc.tsv"), show_col_types = FALSE)
frac_mc <- readr::read_tsv(file.path(path, "biomodal_output", "frac_mc.tsv"), show_col_types = FALSE)

# Number of CpG sites * reads per region
sum_hmc <- readr::read_tsv(file.path(path, "biomodal_output", "sum_hmc.tsv"), show_col_types = FALSE)
sum_mc <- readr::read_tsv(file.path(path, "biomodal_output", "sum_mc.tsv"), show_col_types = FALSE)
sum_total_c <- readr::read_tsv(file.path(path, "biomodal_output", "sum_total_c.tsv"), show_col_types = FALSE)

sum(round(sum_hmc[4:27]/sum_total_c[4:27] - frac_hmc[4:27],2), na.rm=TRUE)
sum(round(sum_mc[4:27]/sum_total_c[4:27] - frac_mc[4:27],2), na.rm=TRUE)

# ---- Get genes that already have 5mc and 5hmc in normal samples ----

# http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/

#  0️⃣ Setup
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(openxlsx)

# Input files
normal_5hmc_files <- list.files(file.path(path, "Benign_prostate_5hmc"), full.names = TRUE)
normal_5mc_files <- list.files(file.path(path, "Healthy_plasma_5mc"), full.names = TRUE)

# 1️⃣ Import GTF
gtf_file <- file.path(path, "biomodal_output", "hg19.refGene.gtf.gz")         # RefSeq GTF
gtf <- rtracklayer::import(con = gtf_file, format = "gtf")

# 2️⃣ Define regions to match with Biomodal 
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

# 3️⃣ Annotate peaks for each file 

# Helper function to get collapsed gene names for a single Hits object
get_collapsed_genes <- function(peaks, hits, subject_gr) {
  # We use seq_along(peaks) here, but define peaks inside lapply
  sapply(seq_along(peaks), function(i) {
    # Extract gene names for the current peak 'i'
    #unique_names <- unique(subject_gr$gene_name[S4Vectors::subjectHits(hits)[S4Vectors::queryHits(hits) == i]])
    unique_names <- unique(S4Vectors::mcols(subject_gr)$gene_name[S4Vectors::subjectHits(hits)[S4Vectors::queryHits(hits) == i]])
    
    
    # Return "NA" if no names are found, else return semicolon-separated string
    if (length(unique_names) == 0) {
      return("NA")
    } else {
      return(paste(unique_names, collapse = ";"))
    }
  })
}

annotate_peaks <- function(peak_file, format) {
  
  # Read the peak file
  peaks <- rtracklayer::import(con = peak_file, format = format)
  
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
}

normal_5hmc <- lapply(normal_5hmc_files, function(f) {annotate_peaks(f, "narrowPeak")})
normal_5mc <- lapply(normal_5mc_files, function(f) {annotate_peaks(f, "BED")})

annotated_list <- list()
annotated_list[["normal_5hmc"]] <- normal_5hmc 
annotated_list[["normal_5mc"]] <- normal_5mc

for (i in names(annotated_list)){
  
  # 4️⃣ Merge genes across samples
  
  # Get all gene IDs in promoter/TSS/gene body for each sample
  promoter_genes <- lapply(annotated_list[[i]], function(df) {
    df %>%
      dplyr::filter(in_promoter & promoter_gene_names != "NA") %>%
      tidyr::separate_longer_delim(cols = promoter_gene_names, delim = ";") %>%
      dplyr::pull(promoter_gene_names) %>%
      unique()
  })
  
  tss_genes <- lapply(annotated_list[[i]], function(df) {
    df %>%
      dplyr::filter(in_tss & tss_gene_names != "NA") %>%
      tidyr::separate_longer_delim(cols = tss_gene_names, delim = ";") %>%
      dplyr::pull(tss_gene_names) %>%
      unique()
  })
  
  gene_body_genes <- lapply(annotated_list[[i]], function(df) {
    df %>%
      dplyr::filter(in_gene_body & gene_body_gene_names != "NA") %>%
      tidyr::separate_longer_delim(cols = gene_body_gene_names, delim = ";") %>%
      dplyr::pull(gene_body_gene_names) %>%
      unique()
  })
  
  # Genes present in ALL normal samples
  common_promoter_genes <- purrr::reduce(.x = promoter_genes, .f = intersect)
  common_tss_genes      <- purrr::reduce(.x = tss_genes, .f = intersect)
  common_gene_body_genes <- purrr::reduce(.x = gene_body_genes, .f = intersect)
  
  # 5️⃣ Save results
  
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
  addWorksheet(wb, i)
  writeData(wb, sheet = i, all_genes_df)
  saveWorkbook(wb, file.path(path, paste0(i, ".xlsx")), overwrite = TRUE)
}


# ---- Identify DMR in tumor ----

# DMR in normal DNA
normal_dmr_hmc <- read.xlsx(file.path(path,  "normal_5hmc.xlsx"))
normal_dmr_mc <- read.xlsx(file.path(path,  "normal_5mc.xlsx"))

# DMR in cfDNA ( = ctDNA + normal DNA)
#dmr_hmc_24 <- readr::read_tsv(file.path(path, "biomodal_output", "DMR_20251124_130348_DMR_hmc_Pre__Post_20251124_130348.tsv"), show_col_types = FALSE)
#dmr_mc_24 <- readr::read_tsv(file.path(path, "biomodal_output", "DMR_20251124_130348_DMR_mc_Pre__Post_20251124_130348.tsv"), show_col_types = FALSE)

dmr_hmc_21 <- readr::read_tsv(file.path(path, "biomodal_output", "DMR_20251208_144543_DMR_hmc_Pre__Post_20251208_144543.tsv"), show_col_types = FALSE)
dmr_mc_21 <- readr::read_tsv(file.path(path, "biomodal_output", "DMR_20251208_144543_DMR_mc_Pre__Post_20251208_144543.tsv"), show_col_types = FALSE)

sig_dmr_hmc <- dmr_hmc_21 %>% 
  dplyr::mutate(mean_mod_group_1 = ifelse(mean_mod_group_1 == 0, 1e-6, mean_mod_group_1),
                mean_mod_group_2 = ifelse(mean_mod_group_2 == 0, 1e-6, mean_mod_group_2),
                mod_fold_change = mean_mod_group_2 / mean_mod_group_1,
                mod_fold_change = log2(mod_fold_change)) %>%
  dplyr::filter(dmr_qvalue <= 0.05, num_contexts >= 10, abs(mod_fold_change) >= 0.58, !is.na(Name))

sig_dmr_mc <- dmr_mc_21 %>% 
  dplyr::mutate(mean_mod_group_1 = ifelse(mean_mod_group_1 == 0, 1e-6, mean_mod_group_1),
                mean_mod_group_2 = ifelse(mean_mod_group_2 == 0, 1e-6, mean_mod_group_2),
                mod_fold_change = mean_mod_group_2 / mean_mod_group_1,
                mod_fold_change = log2(mod_fold_change)) %>% 
  dplyr::filter(dmr_qvalue <= 0.05, num_contexts >= 10, abs(mod_fold_change) >= 0.58, !is.na(Name))

# Remove normal DMR 
tumor_DMR_hmc <- sig_dmr_hmc %>%
  anti_join(normal_dmr_hmc, by = c("Name"="Name", "Annotation"="Annotation")) %>%
  arrange(Annotation, Name) %>%
  dplyr::mutate(key = paste(Chromosome, Start, End, Name, Annotation, sep = "_"))

tumor_DMR_mc <- sig_dmr_mc %>%
  inner_join(normal_dmr_mc, by = c("Name"="Name", "Annotation"="Annotation")) %>%
  arrange(Annotation, Name) %>%
  dplyr::mutate(key = paste(Chromosome, Start, End, Name, Annotation, sep = "_"))

# ---- UMAP ----

metadata <- read.xlsx(file.path(path, "Metadata.xlsx"))
  
# Some genes like CD99, MKSS have more than 1 TSS. So, multiple entries exist.
# Collapse into single entry using max signal.
# logit transform the fraction data
hmc_mat <- frac_hmc[,4:ncol(frac_hmc)] %>% 
  dplyr::group_by(Name, Annotation) %>%
  dplyr::summarize(across(.cols = everything(), .fns = max), .groups = "drop") %>%
  dplyr::inner_join(sig_dmr_hmc %>% dplyr::select(Name, Annotation),
                    by = c("Name", "Annotation")) %>%
  dplyr::mutate(ID = paste0(Name, ".", Annotation)) %>%
  dplyr::select(-Name, -Annotation) %>%
  tibble::column_to_rownames("ID") %>%
  dplyr::rename_with(.fn = function(x){ gsub("_num_hmc_region_frac", "", x)})

# Replace NAs and 0/1 extremes and logit transform
epsilon <- 1e-6
hmc_mat[is.na(hmc_mat)] <- epsilon
hmc_mat[hmc_mat == 0] <- epsilon
hmc_mat[hmc_mat == 1] <- 1 - epsilon
logit_hmc <- log(hmc_mat / (1 - hmc_mat)) %>% as.matrix()

mc_mat <- frac_mc[,4:ncol(frac_mc)] %>% 
  dplyr::group_by(Name, Annotation) %>%
  dplyr::summarize(across(.cols = everything(), .fns = max), .groups = "drop") %>%
  dplyr::inner_join(sig_dmr_mc %>% dplyr::select(Name, Annotation),
                    by = c("Name", "Annotation")) %>%
  dplyr::mutate(ID = paste0(Name, ".", Annotation)) %>%
  dplyr::select(-Name, -Annotation) %>%
  tibble::column_to_rownames("ID") %>%
  dplyr::rename_with(.fn = function(x){ gsub("_num_mc_region_frac", "", x)})

# Replace NAs and 0/1 extremes and logit transform
epsilon <- 1e-6
mc_mat[is.na(mc_mat)] <- epsilon
mc_mat[mc_mat == 0] <- epsilon
mc_mat[mc_mat == 1] <- 1 - epsilon
logit_mc <- log(mc_mat / (1 - mc_mat)) %>% as.matrix()

filename <- "PCA_hmc"
plot_pca(expr_mat = logit_hmc, metadata, top_n_genes = 5000, skip_plot = FALSE, filename, output_dir = path)
filename <- "PCA_mc"
plot_pca(expr_mat = logit_mc, metadata, top_n_genes = 5000, skip_plot = FALSE, filename, output_dir = path)

filename <- "UMAP_hmc"
plot_umap(expr_mat = logit_hmc, metadata, n_pcs = 50, n_neighbors = NULL, filename, output_dir = path)
filename <- "UMAP_mc"
plot_umap(expr_mat = logit_mc, metadata, n_pcs = 50, n_neighbors = NULL, filename, output_dir = path)

# ---- Pie Chart ----

sig_dmr_hmc_all <- dmr_hmc_21 %>% 
  dplyr::mutate(mean_mod_group_1 = ifelse(mean_mod_group_1 == 0, 1e-6, mean_mod_group_1),
                mean_mod_group_2 = ifelse(mean_mod_group_2 == 0, 1e-6, mean_mod_group_2),
                mod_fold_change = mean_mod_group_2 / mean_mod_group_1,
                mod_fold_change = log2(mod_fold_change)) %>%
  dplyr::filter(dmr_qvalue <= 0.05, num_contexts >= 10, abs(mod_fold_change) >= 0.58)

sig_dmr_mc_all <- dmr_mc_21 %>% 
  dplyr::mutate(mean_mod_group_1 = ifelse(mean_mod_group_1 == 0, 1e-6, mean_mod_group_1),
                mean_mod_group_2 = ifelse(mean_mod_group_2 == 0, 1e-6, mean_mod_group_2),
                mod_fold_change = mean_mod_group_2 / mean_mod_group_1,
                mod_fold_change = log2(mod_fold_change)) %>% 
  dplyr::filter(dmr_qvalue <= 0.05, num_contexts >= 10, abs(mod_fold_change) >= 0.58)


plot_piechart(metadata = sig_dmr_hmc_all, segment_col = "Annotation", filename = "5hmc_all", output_dir = path, split_col = NULL)
plot_piechart(metadata = sig_dmr_mc_all, segment_col = "Annotation", filename = "5mc_all", output_dir = path, split_col = NULL)


# ---- Patient wise analysis ----

metadata <- read.xlsx(file.path(path, "Metadata.xlsx"))

# Get all possible comparisons between controls and experiments
samples <- metadata %>%
  dplyr::pull(Sample_ID) %>%
  unique()

combns <- utils::combn(x = samples, m = 2)
controls <- c()
expts <- c()
comparisons <- list()
for (i in 1:ncol(combns)){
  
  a <- gsub(pattern = "C1|C2|C3|EOT", "", x = combns[1, i])
  b <- gsub(pattern = "C1|C2|C3|EOT", "", x = combns[2, i])
  
  if (a == b){
    if(grepl("C1", combns[1,i]) & !grepl("C1", combns[2,i])){
      control <- combns[1, i]
      expt <- combns[2, i]
      controls <- c(controls, control)
      expts <- c(expts, expt)
    } 
  }
}

comparisons[["control"]] <- controls
comparisons[["expt"]] <- expts

# Get fraction of 5mc and 5hmc
frac_hmc <- readr::read_tsv(file.path(path, "biomodal_output", "frac_hmc.tsv"), show_col_types = FALSE)
frac_mc <- readr::read_tsv(file.path(path, "biomodal_output", "frac_mc.tsv"), show_col_types = FALSE)

colnames(frac_hmc) <- gsub("_num_hmc_region_frac", "", colnames(frac_hmc))
colnames(frac_mc) <- gsub("_num_mc_region_frac", "", colnames(frac_mc))

frac_hmc <- frac_hmc %>% 
  dplyr::mutate(across(.cols = everything(), .fns = function(x) { replace(x, is.na(x), 0) })) %>%
  dplyr::mutate(n_UP = 0, n_DOWN = 0)
frac_mc <- frac_mc %>% 
  dplyr::mutate(across(.cols = everything(), .fns = function(x) { replace(x, is.na(x), 0) })) %>%
  dplyr::mutate(n_UP = 0, n_DOWN = 0)

for (i in seq_along(comparisons$control)){
  
  ctrl <- comparisons$control[i]
  expt <- comparisons$expt[i]
  
  # Get regions between Control and Experiment
  frac_hmc <- frac_hmc %>%
    dplyr::mutate(n_UP   = n_UP   + as.integer(.data[[expt]] > .data[[ctrl]]),
                  n_DOWN = n_DOWN + as.integer(.data[[expt]] < .data[[ctrl]]))
  
  # Get regions between Control and Experiment
  frac_mc <- frac_mc %>%
    dplyr::mutate(n_UP   = n_UP   + as.integer(.data[[expt]] > .data[[ctrl]]),
                  n_DOWN = n_DOWN + as.integer(.data[[expt]] < .data[[ctrl]]))
}

frac_hmc <- frac_hmc %>% 
  dplyr::mutate(across(.cols = everything(), .fns = function(x) { replace(x, x == 0, NA) })) %>%
  dplyr::filter(n_UP != n_DOWN, !is.na(Name))
frac_mc <- frac_mc %>% 
  dplyr::mutate(across(.cols = everything(), .fns = function(x) { replace(x, x == 0, NA) })) %>%
  dplyr::filter(n_UP != n_DOWN, !is.na(Name))

# Get pvalues and other stats from biomodal DMR analysis
dmr_hmc_21 <- readr::read_tsv(file.path(path, "biomodal_output", "DMR_20251208_144543_DMR_hmc_Pre__Post_20251208_144543.tsv"), show_col_types = FALSE)
dmr_mc_21 <- readr::read_tsv(file.path(path, "biomodal_output", "DMR_20251208_144543_DMR_mc_Pre__Post_20251208_144543.tsv"), show_col_types = FALSE)

frac_hmc <- frac_hmc %>%
  dplyr::left_join(dmr_hmc_21, by=c("Chromosome", "Start", "End", "Name", "Annotation"))
frac_mc <- frac_mc %>%
  dplyr::left_join(dmr_mc_21, by=c("Chromosome", "Start", "End", "Name", "Annotation"))

# Add "Type" column indicating if region was tumor specific
frac_hmc <- frac_hmc %>%
  mutate(key = paste(Chromosome, Start, End, Name, Annotation, sep = "_")) %>%
  mutate(Type = ifelse(key %in% tumor_DMR_hmc$key, "Tumor only", "Normal")) %>%
  select(-key)  # remove temporary key

frac_mc <- frac_mc %>%
  mutate(key = paste(Chromosome, Start, End, Name, Annotation, sep = "_")) %>%
  mutate(Type = ifelse(key %in% tumor_DMR_mc$key, "Tumor only", "Normal")) %>%
  select(-key)  # remove temporary key

# Reformat and save to excel
frac_hmc <- frac_hmc %>%
  dplyr::mutate(n_Diff = n_UP - n_DOWN) %>%
  dplyr::rename(n_CpGs = num_contexts, FDR = dmr_qvalue, log2FC = mod_fold_change) %>%
  dplyr::select(Chromosome, Start, End, Name, Annotation, Type, n_UP, n_DOWN, n_Diff, n_CpGs, log2FC, FDR, everything())

frac_mc <- frac_mc %>%
  dplyr::mutate(n_Diff = n_UP - n_DOWN) %>%
  dplyr::rename(n_CpGs = num_contexts, FDR = dmr_qvalue, log2FC = mod_fold_change) %>%
  dplyr::select(Chromosome, Start, End, Name, Annotation, Type, n_UP, n_DOWN, n_Diff, n_CpGs, log2FC, FDR, everything())
 
# Write to Excel
wb <- createWorkbook()
addWorksheet(wb, "hmc")
writeData(wb, sheet = "hmc", frac_hmc)
addWorksheet(wb, "mc")
writeData(wb, sheet = "mc", frac_mc)
saveWorkbook(wb, file.path(path, "Biomodal_final_results.xlsx"), overwrite = TRUE)





# dmr_summary <- df %>%
#   dplyr::mutate(Direction = dplyr::case_when(mod_difference > 0 ~ "Up in Carotuximab",
#                                              mod_difference < 0 ~ "Down in Carotuximab",
#                                              TRUE ~ "No Change")) %>% # Should be rare after filtering
#   # Count unique genes by Annotation and Direction
#   dplyr::group_by(Annotation, Direction) %>%
#   dplyr::summarise(Unique_Genes = n_distinct(Name),
#                    .groups = "drop") %>%
#   # Pivot wider for a clean summary table format
#   tidyr::pivot_wider(names_from = Direction,
#                      values_from = Unique_Genes,
#                      values_fill = 0)
# 
# # Display the summary table
# print(dmr_summary)
# 
# #  Comparison directionof 5mc and 5hmc
# comparison_df <- dplyr::inner_join(x = tumor_DMR_hmc %>% dplyr::rename(hmc_diff = mod_difference, hmc_logFC = mod_fold_change) %>% dplyr::select(Name, Annotation, hmc_diff, hmc_logFC ),
#                                    y = tumor_DMR_mc %>% dplyr::rename(mc_diff = mod_difference, mc_logFC = mod_fold_change) %>% dplyr::select(Name, Annotation, mc_diff, mc_logFC ),
#                                    by = c("Name"= "Name", "Annotation"="Annotation")) %>%
#   dplyr::mutate(Trend = dplyr::case_when((mc_diff > 0 & hmc_diff < 0) | (mc_diff < 0 & hmc_diff > 0) ~ "Inverse",
#                                          (mc_diff > 0 & hmc_diff > 0) | (mc_diff < 0 & hmc_diff < 0) ~ "Concordant",
#                                          TRUE ~ "Unknown")) %>%
#   dplyr::arrange(Trend, Annotation, Name) %>%
#   dplyr::select(Name, Annotation, hmc_diff, mc_diff, hmc_logFC, mc_logFC, Trend)
# 
# comparison_df %>% dplyr::count(Trend)

