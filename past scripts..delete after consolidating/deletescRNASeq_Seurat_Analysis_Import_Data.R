#!/usr/bin/env Rscript

# Read and store variables from command line interface (CLI)
cli <- base::commandArgs(trailingOnly = TRUE) 
args <- base::strsplit(x = cli, split = "=", fixed = TRUE)

for (e in args){
  argname <- e[1]
  argval <- e[2]
  assign(argname, argval)
}

#***********************STEP 2A: IDENTIFY EMPTY DROPLETS***********************#

# Empty droplets often contain RNA from the ambient solution, resulting in 
# non-zero counts. These empty droplets MUST be removed before identifying 
# doublets. The emptyDrops function is designed to distinguish between empty 
# droplets and cells by testing each barcode’s expression profile for 
# significant deviation from the ambient profile. The null hypothesis is that
# "transcript molecules are included into droplets by multinomial sampling from
# the ambient profile" i.e. "droplets contain ambient RNA". So, if FDR < 0.05, 
# we reject null hypothesis.

# NOTE: The ouput of emptyDrops() is a dataframe with columns 
# Total   : total UMI count for each barcode (droplet)
# LogProb : log-probability of observing the barcode's count vector under the null model.
# PValue  : the Monte Carlo p-value against the null model
# Limited : whether a lower p-value could be obtained by increasing ‘niters’
# FDR     : Non-empty droplets have FDR < 0.05 (i.e. 5%)
# glimpse(e.out) shows metadata associated with the dataframe also
# ?emptyDrops shows explanation of why there are NA values in FDR

#**************************STEP 2B: IDENTIFY DOUBLETS**************************#

# NOTE: Doublet finding algorithms need good quality data as starting point for
# correct estimation of parameters used in identifying doublets. So, we perform 
# this step after removing empty droplets.

# cell_type has classification based on UMAP clusters
# ucell_class has classification based on UCell scores
# seurat_class has classification based on seurat scores

integ <- add_module_scores(integ, "All Markers")
integ <- annotate_data_score(integ, celltype)
save_data(integ, celltype)


# for (reduc in c("CCA", "RPCA", "Harmony", "JointPCA")){
#   res <- 1.4
#   plot_conserved_modules(res, reduc, celltype, "All Markers")
# }
#plot_pre_integration(sct)

#******************************************************************************#

# # Finding number of common barcodes between samples in Simon snRNASeq
# for (i in samples){
#   barcodes <- rownames(get(i)@meta.data)
#   t <- paste0(i, ".cell")
#   assign(t, barcodes)
# }
# 
# samples.cell <- paste0(samples, ".cell")
# # Identify all possible combinations of 2 samples
# combinations <- utils::combn(x=samples.cell, m=2)
# combinations_result <- c()
# total_sample1 <- c()
# total_sample2 <- c()
# for (i in 1:ncol(combinations)){
#   
#   # Find number of common barcodes for any two samples
#   combinations_result <- c(combinations_result, 
#                            length(intersect(get(combinations[1,i]), get(combinations[2,i]))))
#   total_sample1 <- c(total_sample1, length(get(combinations[1,i])))
#   total_sample2 <- c(total_sample2, length(get(combinations[2,i])))
# }
# 
# results <- data.frame(t(combinations), combinations_result, total_sample1, total_sample2)
# colnames(results) <- c("Sample1", "Sample2", "Number of common barcodes", 
#                        "Total barcodes Sample1", "Total barcodes Sample2")


#******************************************************************************#

# # Check QC metrics for singlets and doublets
# VlnPlot(object = sample.seurat,
#         group.by = "Sample",
#         split.by = "DF.Class",
#         features = c("nUMIs", "nGenes", "MitoRatio", "RiboRatio", "HemeRatio", "Novelty"))

# ggplot2::ggplot(data = filtered_seurat@meta.data, aes(x=nUMIs, y=Sample, fill= after_stat(x))) +
#   geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
#   viridis::scale_fill_viridis(name = "nUMIs", alpha = 1, begin = 0, end = 1, 
#                               direction = 1, discrete = FALSE, option = "D") +
#   labs(title = 'UMI Distribution') +
#   coord_cartesian(xlim = c(100, 10000)) +
#   scale_x_continuous(breaks = c(100, 1000, 2500, 5000, 10000))+
#   #hrbrthemes::theme_ipsum() +
#   theme(legend.position="right",
#         panel.spacing = unit(0.1, "lines"),
#         strip.text.x = element_text(size = 8))
# 
# ggsave("1.jpg", bg="white")  