# # Different ways to access the data
# identical(probeData@featureData@data[1:5,1:5], fData(probeData)[1:5, 1:5]) 
# identical(probeData@phenoData@data[1:5,1:5], pData(probeData)[1:5, 1:5])
# identical(probeData@protocolData@data[1:5,1:5], protocolData(probeData)@data[1:5, 1:5])
# identical(probeData@assayData$exprs[1:5,1:5], exprs(probeData)[1:5, 1:5])
# identical(sData(probeData), dplyr::bind_cols(pData(probeData), protocolData(probeData)@data))

# geomx_dir <- system.file("extdata", "WTA_NGS_Example", package="GeoMxWorkflows")
# DCC_files <- dir(file.path(geomx_dir, "dccs"), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
# PKC_files <- list.files(path = file.path(geomx_dir), pattern = ".pkc", full.names = TRUE, recursive = TRUE)
# meta_file <-  dir(file.path(geomx_dir, "annotation"), pattern = ".xlsx$", full.names = TRUE, recursive = TRUE)

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/Custom_Functions.R")
#source("/hpc/home/kailasamms/projects/scRNASeq/Custom_Functions.R")

# ---- Default Thresholds and Project Parameters ----
# Purpose: Define all default QC thresholds and project-specific parameters
# in a single, centralized location for consistency and easy modification.

# 1. Segment-level QC thresholds
DEFAULT_MIN_SEGMENT_READS     <- 1000  # Minimum reads per segment
DEFAULT_MIN_TRIMMED           <- 80    # Minimum % trimmed reads per segment
DEFAULT_MIN_STITCHED          <- 80    # Minimum % stitched reads per segment
DEFAULT_MIN_ALIGNED           <- 80    # Minimum % aligned reads per segment
DEFAULT_MIN_SATURATION        <- 50    # Minimum % sequencing saturation per segment
DEFAULT_MIN_NEGATIVE_COUNT    <- 10    # Minimum negative control reads per segment
DEFAULT_MAX_NTC_COUNT         <- 1000  # Maximum NTC reads per segment
DEFAULT_MIN_NUCLEI            <- 100   # Minimum nuclei per segment
DEFAULT_MIN_AREA              <- 5000  # Minimum area per segment

# 2. Probe-level QC thresholds
DEFAULT_MIN_PROBE_RATIO       <- 0.1   # Minimum probe geometric mean / target geometric mean
DEFAULT_PERCENT_FAIL_GRUBBS   <- 20    # % of segments where probe fails Grubb's test

# 3. Limit of Quantification (LOQ) thresholds
DEFAULT_LOQ_CUTOFF            <- 2     # Number of geometric SDs above mean
DEFAULT_MIN_LOQ               <- 2     # Minimum LOQ per segment

# 4. Detection thresholds
DEFAULT_MIN_GENE_DETECTION    <- 0.1  # Minimum fraction of genes above LOQ per segment
DEFAULT_MIN_SEGMENT_DETECTION <- 0.1  # Fraction of segments where gene must be detected

# 5. Project parameters
proj <- "Neil_GeoMx"
proj.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data"

proj.params <- list(
  proj              = proj,
  species           = "Homo sapiens",
  data.dir          = file.path(proj.dir, proj, "data"),
  geomx.dir         = file.path(proj.dir, proj, "geomx"),
  
  # Segment-level QC
  minSegmentReads   = DEFAULT_MIN_SEGMENT_READS,
  minTrimmed        = DEFAULT_MIN_TRIMMED,
  minStitched       = DEFAULT_MIN_STITCHED,
  minAligned        = DEFAULT_MIN_ALIGNED,
  minSaturation     = DEFAULT_MIN_SATURATION,
  minNegativeCount  = DEFAULT_MIN_NEGATIVE_COUNT,
  maxNTCCount       = DEFAULT_MAX_NTC_COUNT,
  minNuclei         = DEFAULT_MIN_NUCLEI,
  minArea           = DEFAULT_MIN_AREA,
  
  # Probe-level QC
  minProbeRatio     = DEFAULT_MIN_PROBE_RATIO,
  percentFailGrubbs = DEFAULT_PERCENT_FAIL_GRUBBS,
  
  # LOQ thresholds
  LOQ_cutoff        = DEFAULT_LOQ_CUTOFF,
  minLOQ            = DEFAULT_MIN_LOQ,
  
  # Gene Detection
  minGeneDetection     = DEFAULT_MIN_GENE_DETECTION,
  minSegmentDetection  = DEFAULT_MIN_SEGMENT_DETECTION
)

cat("Project parameters initialized for:", proj, "\n")

# ---- Locate and Load Input Files ----
# Purpose: Collect DCC, PKC, and metadata files and load into a GeoMx object

# 1. Locate input files
DCC_files <- list.files(
  path = proj.params$data.dir, pattern = "\\.dcc$", full.names = TRUE, recursive = TRUE)

PKC_files <- list.files(
  path = proj.params$data.dir, pattern = "\\.pkc$", full.names = TRUE, recursive = TRUE)

meta_file <- list.files(
  path = proj.params$data.dir, pattern = "\\.xlsx$", full.names = TRUE, recursive = TRUE)

# 2. Load Metadata
metadata <- openxlsx::read.xlsx(meta_file)

# 3. Load GeoMx Data
geomxData <- GeomxTools::readNanoStringGeoMxSet(
  dccFiles             = DCC_files,
  pkcFiles             = PKC_files,
  phenoDataFile        = meta_file,              # Excel file with phenotypic data
  phenoDataSheet       = "Template",             # Sheet name with phenotypic data
  phenoDataDccColName  = "Sample_ID",            # Sample ID column
  protocolDataColNames = c("aoi", "roi"),        # Protocol/sequencing columns
  experimentDataColNames = c("panel")            # Panel column
)

# ---- Remove unwanted samples ----
probeData <- geomxData
samples_to_keep <- rownames(pData(probeData))[pData(probeData)$TissueOrigin != "Pancreas"]
probeData <- probeData[, sampleNames(probeData) %in% samples_to_keep]

# ---- Identify Spike-ins, Modules, Negative Controls, and NTCs ----
# Purpose: Extract key assay components for downstream QC and normalization.

# 1. Identify Spike-in Probes
# Spike-ins are synthetic controls added to samples; exclude endogenous and negative controls
spike_in <- setdiff(
  unique(fData(probeData)$CodeClass), 
  c("Endogenous", "Negative"))
cat("Spike-in probes detected:", paste(spike_in, collapse = ", "), "\n")

# 2. Identify Modules
modules <- fData(probeData)[["Module"]] %>% unique()
cat("Modules detected:", paste(modules, collapse = ", "), "\n")

# 3. Identify Negative Control Probes
neg_probes <- fData(probeData) %>% 
  dplyr::filter(CodeClass == "Negative") %>% 
  dplyr::pull(RTS_ID) %>%
  unique()
cat("Negative control probe RTS_IDs:", paste(head(neg_probes, 10), collapse = ", "), "...\n")

# 4. Identify Negative Control Genes
neg_genes <- fData(probeData) %>% 
  dplyr::filter(CodeClass == "Negative") %>% 
  dplyr::pull(TargetName) %>%
  unique()
cat("Negative control target genes:", paste(neg_genes, collapse = ", "), "\n")

# 5. Identify No Template Controls (NTCs)
# Segments labeled as "NTC" in the 'segment' column
ntc_samples <- sData(probeData) %>% 
  dplyr::filter(segment == "NTC") %>% 
  rownames()
cat("NTC segments:", paste(ntc_samples, collapse = ", "), "\n")




# ---- Sankey Diagram ----
sankeyCols <- c("source", "target", "value")

# link1 <- count(pData(probeData), `slide name`, class)
# link2 <- count(pData(probeData),  class, region)
# link3 <- count(pData(probeData),  region, segment)

link1 <- pData(probeData) %>% dplyr::count(TissueOrigin, TissueStatus)
link2 <- pData(probeData) %>% dplyr::count(TissueStatus, TissueType)
link3 <- pData(probeData) %>% dplyr::count(TissueType, segment)

colnames(link1) <- sankeyCols
colnames(link2) <- sankeyCols
colnames(link3) <- sankeyCols

links <- rbind(link1,link2,link3)
nodes <- unique(data.frame(name=c(links$source, links$target)))

# sankeyNetwork is 0 based, not 1 based
links$source <- as.integer(match(links$source,nodes$name)-1)
links$target <- as.integer(match(links$target,nodes$name)-1)


sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name",
              units = "TWh", fontSize = 12, nodeWidth = 30)



# ---- Segment-Level QC ----
# Prior to downstream analysis, we assess sequencing quality and tissue adequacy for each ROI/AOI segment.
# This is an outlier-removal and quality-assessment step applied at the segment level.

# Notes:
# - Segment QC evaluates raw reads, alignment metrics, sequencing saturation, negative probe counts,
#   no-template controls, nuclei counts, and segment area.
# - Segments failing one or more criteria may be excluded or flagged for additional review.
# - Cutoffs may vary depending on tissue type and study design; consistency across segments is key.

# Segment QC Criteria:

# 1. Raw sequencing reads:
#    - Segments with < 1,000 raw reads are removed.

# 2. % Aligned, % Trimmed, or % Stitched reads:
#    - Segments with values < ~80% for any of these metrics are removed.

# 3. % Sequencing saturation ([1 − deduplicated reads / aligned reads] × 100):
#    - Segments < ~50% may require additional sequencing.
#    - Typically, segments below this threshold are not analyzed until improved.

# 4. Negative count:
#    - Calculated as the geometric mean of unique negative probes in the GeoMx panel.
#    - Represents background count per segment.
#    - Low negative counts (1–10) do not automatically exclude segments but indicate
#      potential low endogenous signal or insufficient tissue.

# 5. No Template Control (NTC) count:
#    - Values > 1,000 may indicate contamination.
#    - Segments with NTC counts between 1,000–10,000 can be retained if NTC counts are
#      consistently low across all probes (e.g., 0–2 counts).

# 6. Nuclei count:
#    - Generally, > 100 nuclei per segment is recommended.
#    - Optimal cutoff is tissue- and study-dependent; consistency across segments is more important than a strict threshold.

# 7. Segment area:
#    - Typically correlates with nuclei count.
#    - No strict cutoff is generally applied.

# Shift Zero Counts to One
zero_counts_before <- sum(assayData(probeData)$exprs == 0)
cat("Number of zero counts (before shift):", zero_counts_before, "\n")

probeData <- GeomxTools::shiftCountsOne(
  object = probeData,
  elt = "expr",
  useDALogic = TRUE
)

zero_counts_after <- sum(assayData(probeData)$exprs == 0)
cat("Number of zero counts (after shift):", zero_counts_after, "\n")

# Define Segment QC Parameters
segment_QC_params <- list(
  minSegmentReads   = proj.params$minSegmentReads,
  percentTrimmed    = proj.params$minTrimmed,
  percentStitched   = proj.params$minStitched,    
  percentAligned    = proj.params$minAligned,
  percentSaturation = proj.params$minSaturation,
  minNegativeCount  = proj.params$minNegativeCount,    
  maxNTCCount       = proj.params$maxNTCCount,  
  minNuclei         = proj.params$minNuclei,
  minArea           = proj.params$minArea
)

# Apply Segment QC Flags
probeData <- GeomxTools::setSegmentQCFlags(
  object    = probeData, 
  qcCutoffs = segment_QC_params
)

cat("Segment QC flags added to probeData@protocolData.\n")

# Calculate geometric means of negative control probes per segment for each module
negativeGeoMeans_list <- list()
for (module in modules){
  
  # Get negative probes for each module
  neg_probes <- fData(probeData) %>%
    dplyr::filter(Module == module, CodeClass == "Negative") %>%
    dplyr::pull(RTS_ID)
  
  # Calculate geometric mean of negative control probes for each segment
  neg_geo_means <- apply(X = exprs(probeData)[neg_probes, , drop = FALSE], 
                            MARGIN = 2, 
                            FUN =  function(x){exp(mean(log(x[x>0])))})
  
  # Save to a named list
  negativeGeoMeans_list[[paste0("NegGeoMean_", module)]] <- neg_geo_means
}

# Convert list to dataframe and add segment names
negativeGeoMeans <- dplyr::bind_cols(negativeGeoMeans_list) %>% as.matrix()
rownames(negativeGeoMeans) <- colnames(exprs(probeData))

# negativeGeoMeans1 <- NanoStringNCTools::esBy(
#   X     = negativeControlSubset(x = probeData), 
#   GROUP = "Module", 
#   FUN   = function(x) {
#     assayDataApply(X = x, MARGIN = 2, FUN = ngeoMean, elt = "exprs")
#   }
# )
# colnames(negativeGeoMeans1) <- paste0("NegGeoMean_", colnames(negativeGeoMeans1))

# Add negativeGeoMeans to protocolData
for (col in colnames(negativeGeoMeans)) {
  protocolData(probeData)[[col]] <- as.data.frame(negativeGeoMeans)[[col]]
}

# ---- Plot and tabulate Segment QC ----

# QC Histogram Function
QC_histogram <- function(assay_data, x_axis, fill_by, cutoff, scale_trans = NULL) {
  
  # Skip plotting if x_axis column not found
  if (!(x_axis %in% colnames(assay_data))) {
    warning(paste("Skipping plot:", x_axis, "not found in assay_data."))
    return(NULL)
  }
  
  p <- ggplot2::ggplot(
    data    = assay_data,
    mapping = aes_string(x = paste0("unlist(`", x_axis, "`)"), fill = fill_by)
  ) +
    ggplot2::geom_histogram(bins = 50, color = "black") +
    ggplot2::geom_vline(xintercept = cutoff, lty = "dashed", color = "red") +
    ggplot2::theme_bw() +
    ggplot2::guides(fill = "none") +
    ggplot2::facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    ggplot2::labs(x = x_axis, y = "Segments (#)", title = x_axis)
  
  if (!is.null(scale_trans)) {
    p <- p + ggplot2::scale_x_continuous(trans = scale_trans)
  }
  return(p)
}

# Generate QC Plots
p1 <- QC_histogram(sData(probeData), "Trimmed (%)",   "segment", segment_QC_params$percentTrimmed)
p2 <- QC_histogram(sData(probeData), "Stitched (%)",  "segment", segment_QC_params$percentStitched)
p3 <- QC_histogram(sData(probeData), "Aligned (%)",   "segment", segment_QC_params$percentAligned)
p4 <- QC_histogram(sData(probeData), "Saturated (%)", "segment", segment_QC_params$percentSaturation) +
  ggplot2::labs(title = "Sequencing Saturation (%)", x = "Sequencing Saturation (%)")
p5 <- QC_histogram(sData(probeData), "area", "segment", segment_QC_params$minArea, scale_trans = "log10")
p6 <- QC_histogram(sData(probeData), "nuclei", "segment", segment_QC_params$minNuclei)

plot_list <- list(p1, p2, p3, p4, p5, p6)

# Remove NULL plots (avoids blank space in final figure)
plot_list <- Filter(Negate(is.null), plot_list)

# Add plots for Negative Control Geometric Means
for (col in colnames(negativeGeoMeans)) {
  p <- QC_histogram(sData(probeData), col, "segment", 2, scale_trans = "log10")
  plot_list <- c(plot_list, list(p))
}

# Combine and Save Plots
combined_plot <- cowplot::plot_grid(
  plotlist = plot_list,
  ncol     = 2
)

ggplot2::ggsave(
  filename = file.path(proj.params$geomx.dir, "QC_histograms.jpg"),
  plot     = combined_plot,
  width    = 2 * 8, 
  height   = ceiling(length(plot_list) / 2) * 6,
  units    = "in",
  dpi      = 300
)

# Summarize Segment QC Flags
segment_qc_summary <- protocolData(probeData)[["QCFlags"]] %>%
  dplyr::mutate(QCStatus = dplyr::case_when(rowSums(across(everything())) == 0 ~ FALSE,
                                            TRUE                              ~ TRUE)) %>%
  dplyr::summarise(across(everything(), list(PASS    = ~sum(!.x, na.rm = TRUE),
                                             WARNING = ~sum(.x, na.rm = TRUE)))) %>%
  tidyr::pivot_longer(everything(),
                      names_to  = c("Column", ".value"),
                      names_sep = "_") %>%
  tibble::column_to_rownames("Column")

cat("\n--- Segment QC Summary ---\n")
print(segment_qc_summary)

# Summarize NTC Counts (if available)
if ("NTC" %in% colnames(sData(probeData))) {
  ntc_df <- as.data.frame(table(NTC_Count = sData(probeData)$NTC), stringsAsFactors = FALSE)
  colnames(ntc_df) <- c("NTC Count", "# of Segments")
  
  cat("\n--- NTC Count Summary ---\n")
  print(ntc_df)
} else {
  warning("No 'NTC' column found in sData(probeData). Skipping NTC summary.")
}

# ---- Remove Flagged Segments ----
# Compute QC Status
protocolData(probeData)[["QCStatus"]] <- protocolData(probeData)[["QCFlags"]] %>%
  dplyr::mutate(QCStatus = dplyr::case_when(rowSums(across(everything())) == 0 ~ "PASS",
                                            TRUE                              ~ "WARNING")) %>%
  dplyr::pull(QCStatus)

cat("\n--- Dimensions Before Segment QC Filtering ---\n")
print(dim(probeData))

# Remove Flagged Segments
failed_segments <- rownames(sData(probeData))[sData(probeData)[["QCStatus"]] != "PASS"]
probeData <- probeData[, protocolData(probeData)[["QCStatus"]] == "PASS"]

cat("\n--- Dimensions After Segment QC Filtering ---\n")
print(dim(probeData))



# ---- Probe-Level QC ----
# Prior to summarizing data into gene-level counts, we remove low-performing probes.
# This is an outlier-removal step applied to probes representing negative controls.

# Notes:
# - Endogenous genes in WTA data typically have a single probe per target, so Probe QC
#   does not apply to them.
# - Negative control probes, which do not target any genomic sequence, often have multiple
#   probes per target. These enable calculation of segment-specific background and are
#   important for downstream gene detection.
# - Probe QC only removes probes; every gene target retains at least one probe. Genes are
#   never removed.

# Global Probe Removal Criteria:
# A probe is removed from the dataset entirely if either:
# 1. Its geometric mean across all segments, divided by the geometric mean of all probes
#    representing the same target, is < 0.1.
# 2. It is identified as an outlier via Grubb’s test in ≥ 20% of segments.

# Local Probe Removal Criteria:
# A probe is removed from a specific segment if it is an outlier according to Grubb’s test
# in that segment. Set `removeLocalOutliers = FALSE` to skip local outlier removal.

# Define Probe QC Parameters
probe_QC_params <- list(
  minProbeRatio     = proj.params$minProbeRatio,
  percentFailGrubbs = proj.params$percentFailGrubbs
)

# Apply Probe QC Flags
probeData <- GeomxTools::setBioProbeQCFlags(
  object             = probeData, 
  qcCutoffs          = probe_QC_params, 
  removeLocalOutliers = TRUE
)

cat("Probe QC flags added to probeData@featureData.\n")

# ---- Tabulate Probe QC ----
# Summarize Probe QC Flags
probe_qc_summary <- fData(probeData)[["QCFlags"]] %>%
  dplyr::mutate(LowRatio = LowProbeRatio,
                Global   = GlobalGrubbsOutlier,
                Local    = rowSums(across(-c(LowProbeRatio, GlobalGrubbsOutlier))) > 0) %>%
  dplyr::summarise(across(c(LowRatio, Global, Local), list(PASS    = ~sum(!.x, na.rm = TRUE),
                                             WARNING = ~sum(.x, na.rm = TRUE)))) %>%
  tidyr::pivot_longer(everything(),
                      names_to  = c("Column", ".value"),
                      names_sep = "_") %>%
  tibble::column_to_rownames("Column")

cat("\n--- Probe QC Summary ---\n")
print(probe_qc_summary)

# ---- Remove Flagged Probes ----
cat("\n--- Dimensions Before Probe QC Filtering ---\n")
print(dim(probeData))

failed_probes <- rownames(fData(probeData))[fData(probeData)[["QCFlags"]]$LowProbeRatio | fData(probeData)[["QCFlags"]]$GlobalGrubbsOutlier]
probeData <- probeData[!fData(probeData)[["QCFlags"]]$LowProbeRatio & 
                         !fData(probeData)[["QCFlags"]]$GlobalGrubbsOutlier, ]

cat("\n--- Dimensions After Probe QC Filtering ---\n")
print(dim(probeData))


# ---- Create Gene-Level Count Data ----
# For genes represented by multiple probes per segment, counts are summarized
# using the geometric mean of the probes.

cat("\n--- Dimensions Before Aggregation ---\n")
print(dim(probeData))

geneData <- GeomxTools::aggregateCounts(
  object = probeData,
  FUN    = ngeoMean
)

cat("\n--- Dimensions After Aggregation ---\n")
print(dim(geneData))



# ---- Limit of Quantification (LOQ) ----
# The LOQ is a segment-specific and module-specific threshold.
# LOQ defines the minimum expression level at which a probe’s signal can be 
# reliably distinguished from background noise. It is calculated based on 
# negative control probes (probes that should not hybridize to any target sequence).

# Notes:
# - Targets with expression below the LOQ are considered unreliable or 
#   indistinguishable from background.
# - LOQ calculations are more stable in larger segments, which contain more 
#   negative control probes. Segments with very few negative control probes 
#   (e.g., <2) may yield less reliable LOQs.
# - Formula for the LOQ in the ith segment:
#     LOQ_i = geomean(NegProbe_i) * (geoSD(NegProbe_i) ^ n)
# - We typically use n = 2 (2 geometric standard deviations above the geometric mean).
# - A minimum LOQ of 2 is recommended for segments where the calculated LOQ falls below this threshold.

# Calculate Limit of Quantification (LOQ) per Segment for each module using
# NegGeoMean and NegGeoSD we calculated earlier
LOQ <- data.frame(row.names = colnames(geneData))

for(module in modules) {
  mean_col <- paste0("NegGeoMean_", module)
  sd_col   <- paste0("NegGeoSD_", module)
  LOQ_col <- paste0("LOQ_", module)
  
  if(all(c(mean_col, sd_col) %in% colnames(pData(geneData)))) {
    LOQ[[LOQ_col]] <-  pmax(proj.params$minLOQ, 
                            pData(geneData)[[mean_col]] * pData(geneData)[[sd_col]] ^ proj.params$LOQ_cutoff)
  }
}

# Add LOQ values to geneData protocolData
pData(geneData)$LOQ <- LOQ
cat("LOQ values added to geneData@protocolData.\n")

# Generate LOQ Logical Matrix
# This creates a logical matrix indicating whether each gene's expression
# in each segment is above the segment-specific LOQ.
# NOTE: A given gene cannot be present in multiple modules.
LOQ_Mat <- NULL 

for (module in modules) {
  # Construct the LOQ column name for this module
  LOQ_col <- paste0("LOQ_", module)
  
  # Get feature (gene) names belonging to this module
  genes <- fData(geneData) %>% 
    dplyr::filter(Module == module) %>% 
    rownames() 
  
  # Extract expression values for these genes across all segments
  expr_mat <- exprs(geneData)[genes, , drop = FALSE]
  
  # Extract LOQ values for each segment corresponding to this module
  # Ensure segments (columns) are aligned with expr_mat
  loq_values <- pData(geneData)$LOQ[colnames(expr_mat), LOQ_col]
  
  # Compare each gene's expression to the segment-specific LOQ
  # Result is a logical matrix: TRUE if expression > LOQ, FALSE otherwise
  log_mat <- base::apply(
    X      = expr_mat,
    MARGIN = 1,  # Apply over rows (genes)
    FUN    = function(x) { x > loq_values }
  ) %>% t()
  
  # Append results to the cumulative LOQ matrix
  # rbind() is fine since there are no common genes between modules
  LOQ_Mat <- base::rbind(LOQ_Mat, log_mat)
}


# ---- Gene Detection Rate per Segment QC ----
# Segments with exceptionally low signal are filtered out, as they typically have 
# a small fraction of panel genes detected above the LOQ compared to other segments.
#
# Key points:
# - A typical detection threshold for segment filtering is 5–10% of genes above LOQ.
# - Thresholds may need adjustment based on experimental design factors 
#   (e.g., segment type, area, number of nuclei) and tissue characteristics 
#   (e.g., type, age).
# - Filtering low-detection segments ensures that downstream analyses are 
#   based on reliable, informative data.

# Count number of genes above LOQ per segment
pData(geneData)[["GenesDetected"]] <- colSums(LOQ_Mat, na.rm = TRUE)

# Calculate detection rate (fraction of panel genes above LOQ)
pData(geneData)[["GeneDetectionRate"]] <- pData(geneData)[["GenesDetected"]] / nrow(geneData)

# Categorize detection rate into bins for plotting
pData(geneData)[["DetectionThreshold"]] <- pData(geneData) %>%
  dplyr::mutate(
    DetectionThreshold = dplyr::case_when(
      GeneDetectionRate > 0.15 ~ ">15%",
      GeneDetectionRate > 0.10 ~ "10-15%",
      GeneDetectionRate > 0.05 ~ "5-10%",
      GeneDetectionRate > 0.01 ~ "1-5%",
      TRUE                      ~ "<1%"), 
    DetectionThreshold = factor(DetectionThreshold,
                                levels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))) %>%
  dplyr::pull(DetectionThreshold)

# ---- Plot and tabulate Gene Detection Rate per Segment QC ----

# Stacked bar plot of different cut points
ggplot2::ggplot(data = pData(geneData), 
                mapping = aes(x = DetectionThreshold)) +
  ggplot2::geom_bar(aes(fill = segment)) +          # stacked bars by segment type
  ggplot2::geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +  # counts on top
  ggplot2::theme_bw() +
  ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ggplot2::labs(x    = "Gene Detection Rate",
    y    = "Segments, #",
    fill = "Segment Type")

ggsave(filename = file.path(proj.params$geomx.dir, "QC_Gene_Detection_Rate.jpg"),
  plot     = last_plot(),
  width    = 6,
  height   = 6,
  units    = "in",
  dpi      = 300
)

# ---- Remove low-detection segments ----
cat("\n--- Dimensions Before Gene Detection Filtering ---\n")
print(dim(geneData))

failed_segments <- c(failed_segments, rownames(sData(geneData))[pData(geneData)[["GeneDetectionRate"]] < proj.params$minGeneDetection])
geneData <- geneData[, pData(geneData)[["GeneDetectionRate"]] >= proj.params$minGeneDetection]

cat("\n--- Dimensions After Gene Detection Filtering ---\n")
print(dim(geneData))


# ---- Overall Gene Detection Rate QC ----
# This analysis visualizes the total number of genes detected across segments,
# helping to assess global gene detection in the study and guide low-detection
# gene filtering.
#
# Key points:
# - Filtering genes detected in very few segments improves downstream statistical
#   performance and reduces noise, enhancing true biological signals.
# - Typical % segment cutoffs for filtering range from 5–20%, depending on the
#   biological diversity of the dataset.
# - If a key gene is naturally expressed in very few segments (<5%), it can be
#   retained manually to ensure important biological information is preserved.
#
# This visualization allows informed selection of gene-level filtering thresholds
# and ensures that filtering decisions are data-driven rather than arbitrary.

# Calculate detection rate for each gene across segments:

# Ensure LOQ_Mat columns match geneData segments
LOQ_Mat <- LOQ_Mat[, colnames(geneData)]

# Count the number of segments each gene is detected in
fData(geneData)[["DetectedSegments"]] <- rowSums(LOQ_Mat, na.rm = TRUE)

# Calculate segment-level detection rate (fraction of segments where gene is detected)
fData(geneData)[["SegmentDetectionRate"]] <- fData(geneData)[["DetectedSegments"]] / nrow(pData(geneData))

# Categorize detection rate into bins for plotting
fData(geneData)[["SegmentDetectionThreshold"]] <- fData(geneData) %>%
  dplyr::mutate(
    SegmentDetectionThreshold = dplyr::case_when(
      SegmentDetectionRate > 0.5  ~ ">50%",
      SegmentDetectionRate > 0.3  ~ "30-50%",
      SegmentDetectionRate > 0.2  ~ "20-30%",
      SegmentDetectionRate > 0.1  ~ "10-20%",
      SegmentDetectionRate > 0.05 ~ "5-10%",
      TRUE                         ~ "<1%"),
    SegmentDetectionThreshold = factor(
      SegmentDetectionThreshold,
      levels = c("<1%", "5-10%", "10-20%", "20-30%", "30-50%", ">50%"))) %>%
  dplyr::pull(SegmentDetectionThreshold)

# ---- Plot and tabulate Gene Detection Rate across segments ----

# Stacked bar plot of different cut points
ggplot2::ggplot(data = fData(geneData), 
                mapping = aes(x = SegmentDetectionThreshold)) +
  ggplot2::geom_bar() +
  ggplot2::geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  ggplot2::theme_bw() +
  ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ggplot2::labs(
    x = "% of Segments",
    y = "Genes, #",
    fill = "Segment Type")

ggsave(filename = file.path(proj.params$geomx.dir, "QC_Segment_Detection_Rate.jpg"),
  plot     = last_plot(),
  width    = 6,
  height   = 6,
  units    = "in",
  dpi      = 300)

# ---- Remove low-detected genes ----

cat("\n--- Dimensions Before Gene Detection Filtering ---\n")
print(dim(geneData))

# Keep genes detected in at least proj.params$minGeneDetection fraction of segments
# Always include negative control genes
failed_probes <- c(failed_probes, rownames(fData(geneData))[fData(geneData)[["SegmentDetectionRate"]] < proj.params$minSegmentDetection])
geneData <- geneData[fData(geneData)[["SegmentDetectionRate"]] >= proj.params$minSegmentDetection  |
                                     fData(geneData)$TargetName %in% neg_genes, 
]

cat("\n--- Dimensions After Gene Detection Filtering ---\n")
print(dim(geneData))


# ---- Normalization ----

# Q3 normalization is the preferred strategy for most DSP-NGS RNA studies. 
# It normalizes counts based on the upper quartile (75th percentile) of counts 
# per segment. Background normalization is avoided here due to low negative 
# probe counts, which can make background estimates unstable.

# Assess Q3 vs. Negative Control Signals
# Before performing normalization, explore the relationship between:
# 1) Q3 of counts per segment
# 2) Geometric mean of negative control probes per segment

# Goal: ensure separation between these two values for stable normalization.
# If the values are too close, it may indicate low-signal segments, which 
# could benefit from additional filtering of low-quality segments or genes.

# Visual assessment can guide whether stricter filtering is required prior 
# to normalization.

# Q3 vs Negative Control Probe Assessment
# Purpose: Examine relationship between segment Q3 counts and geometric mean
# of negative control probes (NegProbe) to ensure stable normalization.

ann_of_interest <- "segment" # column in metadata

# Prepare data frame
Stat_data <- data.frame(
  row.names   = colnames(exprs(geneData)),
  Segment     = colnames(exprs(geneData)),
  Annotation  = pData(geneData)[, ann_of_interest],
  Q3          = unlist(apply(X = exprs(geneData), MARGIN = 2, FUN = function(x) quantile(x, probs = 0.75, na.rm = TRUE))),
  NegProbe    = exprs(geneData)[neg_genes, ])

# Pivot for histogram plotting
Stat_data_m <- Stat_data %>%
  tidyr::pivot_longer(
    cols       = c(Q3, NegProbe),   # columns to pivot
    names_to   = "Statistic",       # new column for original column names
    values_to  = "Value"            # new column for the values
  ) %>%
  # Compute range per Statistic (optional, for ordering)
  dplyr::group_by(Statistic) %>%
  dplyr::mutate(Range = max(Value) - min(Value)) %>%
  dplyr::ungroup() %>%
  # Order factor by descending range
  dplyr::mutate(
    Statistic = factor(Statistic,
                       levels = names(sort(tapply(Value, Statistic, function(x) max(x)-min(x)), decreasing = TRUE)))) %>%
  dplyr::select(-Range) %>%  # remove temporary Range column
  as.data.frame()

# It is important to set levels as c("Q3", "NegProbe") because Q3 values
# are usually larger than NegProbe values. This ensures proper color assignment
# and plotting order.

# Histogram of Q3 vs NegProbe per segment
plt1 <- ggplot(Stat_data_m, aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) +
  theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) +
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #", title = "Histogram of Q3 vs NegProbe Counts")

# Scatter plot: NegProbe vs Q3 per segment
plt2 <- ggplot(Stat_data, aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() +
  guides(color = "none") +
  theme_bw() +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts", title = "Q3 vs NegProbe Counts")

# Scatter plot: Q3 / NegProbe ratio
plt3 <- ggplot(Stat_data, aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() +
  theme_bw() +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 / NegProbe Value", title = "Q3 / NegProbe Ratio")

# Combine plots using cowplot
btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""), rel_widths = c(0.43, 0.57))
combined_plot <- plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))
ggsave(
  filename = file.path(proj.params$geomx.dir, "QC_Q3_vs_NegProbe.jpg"),
  plot     = combined_plot,
  width    = 10,
  height   = 6,
  units    = "in",
  dpi      = 300
)

# 1. Background Normalization
# This normalizes counts relative to negative controls (neg_norm).
# This is robust for DSP-NGS RNA studies without spike-ins.
geneData <- GeomxTools::normalize(
  object     = geneData,
  norm_method = "neg",
  fromElt     = "exprs",
  toElt       = "neg_norm"
)

# 2. Q3 Normalization
# Normalizes counts to the 75th percentile per segment (q_norm)
# This is robust for DSP-NGS RNA studies, with or without spike-ins.
geneData <- GeomxTools::normalize(
  object        = geneData,
  norm_method   = "quant",
  desiredQuantile = 0.75,
  fromElt       = "exprs",      # Can also use "neg_norm" if background-normalized counts preferred
  toElt         = "q_norm"
)

# Visualize the first 10 segments with each normalization method
plot_box <- function(data, cols, main_title, ylab_title, names = 1:10) {
  boxplot(
    data[, cols],
    col   = cols,
    main  = main_title,
    log   = "y",            # log scale for better visualization of wide-range counts
    names = names,
    xlab  = "Segment",
    ylab  = ylab_title
  )
}

# Prepare inputs
plot_list <- list(
  list(data = assayDataElement(geneData, "exprs"),    col = "#9EDAE5", main = "Raw Counts", ylab = "Counts, Raw"),
  list(data = assayDataElement(geneData, "neg_norm"), col = "#FF7F0E", main = "Neg Norm Counts", ylab = "Counts, Neg. Normalized"),
  list(data = assayDataElement(geneData, "q_norm"),   col = "#2CA02C", main = "Q3 Norm Counts", ylab = "Counts, Q3 Normalized")
)

# Generate combined boxplots
jpeg(file.path(proj.params$geomx.dir, "QC_Normalizations.jpg"), width = 2000, height = 2500, res = 300)  # high-res output
par(mfrow = c(3, 1), mar = c(4, 4, 3, 2))  # 3 rows, 1 column, adjust margins

for (p in plot_list) {
  boxplot(
    p$data, 
    col   = p$col, 
    main  = p$main, 
    log   = "y",           # log scale for wide-range count distributions
    names = NULL,          # optional: can supply segment names
    xlab  = "Segment", 
    ylab  = p$ylab)
}

dev.off()  # close the plotting device

# ---- Plot UMAP ----
# Set UMAP parameters
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42  # for reproducibility

# Run UMAP on log2-transformed Q3-normalized counts
umap_out <- umap(
  t(log2(assayDataElement(geneData, elt = "q_norm"))), 
  config = custom_umap)

# Store UMAP coordinates in pData
pData(geneData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, 1:2]

# Generate UMAP plot
ggplot(data = pData(geneData),
       aes(x = UMAP1, y = UMAP2, color = TissueType)) + # shape = TissueOrigin)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = gsub("DSP-1001660006572-H-|\\.dcc", "", rownames(pData(geneData)))),
    size = 3,  # Adjust text size
    max.overlaps = 50, # Allow more overlaps to be resolved
    box.padding = 0.3, # Space around text
    point.padding = 0.2, # Space between text and points
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  # geom_text(
  #   aes(label = gsub("DSP-1001660006572-H-|\\.dcc", "", rownames(pData(geneData)))),
  #   hjust = 0, nudge_x = 0.1, show.legend = FALSE) +
  scale_color_manual(values = custom_palette) +  # pre-defined color palette
  theme_bw() +
  labs(title = "UMAP of Segments (Q3-normalized counts)",
       x = "UMAP1", y = "UMAP2")

# Save UMAP plot
ggsave(filename = file.path(proj.params$geomx.dir, "QC_UMAP.jpg"),
  plot     = last_plot(),
  width    = 11,
  height   = 11,
  units    = "in",
  dpi      = 300)

# ---- Plot heatmap of Highly Variable Genes (HVGs) ----
# Purpose: Use coefficient of variation (CV) to find genes with high variability 
# across segments. These genes often show strong biological differences between 
# conditions or tissues.

# Log2 transformation of Q3-normalized counts
assayDataElement(geneData, elt = "log_q") <- 
  assayDataApply(geneData, 2, FUN = log, base = 2, elt = "q_norm")

# Define CV function
calc_CV <- function(x) {
  sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}

# Calculate CV per gene
CV_dat <- assayDataApply(geneData, elt = "log_q", MARGIN = 1, FUN = calc_CV)

# Display top 5 most variable genes
top_CV_genes <- sort(CV_dat, decreasing = TRUE)[1:5]
print(top_CV_genes)

# Select genes in the top 20% CV for heatmap
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]

# Plot heatmap of highly variable genes
pheatmap::pheatmap(
  assayDataElement(geneData[GOI, ], elt = "log_q"),
  scale = "row",                          # row-wise z-score normalization
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  clustering_method = "average",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  breaks = seq(-3, 3, 0.05),             # fine-grained color breaks
  color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
  annotation_col = pData(geneData)[, c("TissueType", "segment")],
  filename = file.path(proj.params$geomx.dir, "QC_Heatmap.jpg"))






# ---- DEG analysis using LMM ----

# [OPTIONAL] Remove segments based on PCA plot
remove_patterns <- c("F08", "B12", "B01")
keep_idx <- !grepl(paste(remove_patterns, collapse = "|"),
                   sampleNames(geneData))
geneData <- geneData[, keep_idx]

# Apply log2 transformation to the expression data
assayDataElement(object = geneData, elt = "log_q") <- 
  assayDataApply(geneData, 2, FUN = log, base = 2, elt = "q_norm")

# Define test_col to contain the groups to be compared
pData(geneData)[["test_col"]] <- factor(pData(geneData)[["TissueType"]])

# Design, expression matrix, metadata
expr <- assayData(geneData)[["log_q"]]
meta <- pData(geneData)
design <- ~ test_col + (1 | CaseNumber)
# design <- ~ test_col + (1 + test_col | slide)   # if tubule and glomerulus are present on same slide
# design <- ~ test_col + (1 | slide)             # if tubule and glomerulus are present on different slides 

# ---- Direct lmer method [slow for many genes] ----

# Function: run lmer + emmeans for one gene
run_lmer_gene <- function(dat, design, gene_id) {
  
  # Skip if not enough levels
  if (length(unique(dat[["test_col"]])) < 2) return(NULL)
  
  model_formula <- stats::update(design, y ~ .) 
  
  # Fit model
  fit <- tryCatch(
    lme4::lmer(formula = model_formula, data = dat),
    error = function(e) return(NULL)
  )
  if (is.null(fit)) return(NULL)
  
  # Pairwise contrasts
  emm <- emmeans::emmeans(object = fit, specs = ~ test_col)
  res <- emmeans::contrast(object = emm, method = "pairwise") %>%
    base::summary(infer = TRUE, adjust = "fdr")
  
  # Add gene name
  res_tbl <- as.data.frame(res) %>%
    mutate(Gene = gene_id) %>%
    relocate(Gene)
  
  return(res_tbl)
}

if (all(rownames(meta) == colnames(expr))){
  
  dat <- meta
  results <- vector("list", nrow(expr))
  for(i in seq_len(nrow(expr))) {
    
    gene_id <- rownames(expr)[i]
    dat$y <- expr[gene_id, ]    # just overwrite the 'y' column
  
    results[[i]] <- run_lmer_gene(dat, design, gene_id)
  }
}

results <- dplyr::bind_rows(results) %>%
  dplyr::mutate(FDR = p.adjust(.data[["p.value"]], method = "fdr")) %>%
  dplyr::rename(SYMBOL = Gene,
                log2FoldChange = estimate,
                lfcSE = SE,
                stat = t.ratio,
                pvalue = p.value,
                padj = FDR) %>%
  dplyr::select(SYMBOL, contrast, log2FoldChange, lfcSE, stat, pvalue, padj, everything())


# ---- Default GeoMx method [fast for many genes] ----

# If you want to subset data, define column to use for subsetting
results <- data.frame()
subset_col <- NULL #"class"
subset_vals <- if(is.null(subset_col)) NA else unique(pData(geneData)[[subset_col]])

# Subset geomx object
subset_obj <- function(geomx_obj, subset_col = NULL, subset_val = NULL){
  
  if(is.null(subset_col) || is.na(subset_val)){
    # No subsetting, return full object
    return(geomx_obj)
  }
  
  # Subset data
  ind <- pData(geomx_obj)[[subset_col]] == subset_val
  sub_obj <- geomx_obj[, ind]
  
  return(sub_obj)
}

run_lme4 <- function(geomx_obj, design, subset_val = NA) {
  
  # Run mixed model
  mixedOut <- GeomxTools::mixedModelDE(object = geomx_obj,
                                       elt = "log_q",
                                       modelFormula = design,
                                       groupVar = "test_col",
                                       nCores = 1,
                                       multiCore = FALSE)
  
  # Extract results
  r_test <- do.call(rbind, mixedOut["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- r_test %>%
    as.data.frame() %>%
    dplyr::mutate(Contrast = tests,
                  Gene = unlist(lapply(colnames(mixedOut), rep, nrow(mixedOut["lsmeans", ][[1]]))),
                  Subset = ifelse(is.null(subset_val), NA, subset_val),
                  FDR = p.adjust(.data$`Pr(>|t|)`, method = "fdr")) %>%
    dplyr::select(Gene, Subset, Contrast, Estimate, `Pr(>|t|)`, FDR) %>%
    tibble::remove_rownames()
  
  return(r_test)
}

# Run LMM for each subset
for (subset_val in subset_vals){
  
  sub_obj <- subset_obj(geneData, subset_col, subset_val)
  res <- run_lme4(sub_obj, design, subset_val)
  results <- dplyr::bind_rows(results, res)
}

results <- results %>%
  dplyr::rename(SYMBOL = Gene,
                log2FoldChange = Estimate,
                pvalue = `Pr(>|t|)`,
                padj = FDR) %>%
  dplyr::select(SYMBOL, Subset, Contrast, log2FoldChange, pvalue, padj, everything())


# --- Save as xlsx ----

wb <- createWorkbook()
addWorksheet(wb, "All_Results")
writeData(wb, sheet = "All_Results", results)
saveWorkbook(wb, file = file.path(proj.params$geomx.dir, "lmer_results.xlsx"), overwrite = TRUE)


q_norm_df   <- Biobase::assayDataElement(geneData, "q_norm") %>% 
  as.data.frame()

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Q3_norm")
openxlsx::writeData(wb, sheet = "Q3_norm", x = q_norm_df , rowNames = TRUE)
writeData(wb, sheet = "Q3_norm", x = "SYMBOL", startCol = 1, startRow = 1)
openxlsx::saveWorkbook(wb, file = file.path(proj.params$geomx.dir, "Q3_Norm_Counts.xlsx"), overwrite = TRUE)
openxlsx::saveWorkbook(wb, file = file.path("Q3_Norm_Counts.xlsx"), overwrite = TRUE)

# ---- ----
