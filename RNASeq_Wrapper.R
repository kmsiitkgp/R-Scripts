#!/usr/bin/env Rscript

# Read and store variables from command line interface (CLI)
cli <- base::commandArgs(trailingOnly = TRUE) 
args <- base::strsplit(x = cli, split = "=", fixed = TRUE)

for (e in args){
  argname <- e[1]
  argval <- e[2]
  assign(argname, argval)
}

# Run the Custom_Functions.R script
path1 <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/Custom_Functions.R"
path2 <- "/hpc/home/kailasamms/projects/scRNASeq/Custom_Functions.R"
if (file.exists(path1)) {
  source(path1)
} else if (file.exists(path2)) {
  source(path2)
}

# ---- PROJECT SET UP ----

parent.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data"
gmt.dir    <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/GSEA_genesets"
scripts.dir <- NULL

# parent.dir  <- "/hpc/home/kailasamms/scratch"
# gmt.dir     <- "/hpc/home/kailasamms/projects/GSEA_genesets"
# scripts.dir <- "/hpc/home/kailasamms/projects/scRNASeq"

proj.params <- setup_project(proj             = proj,
                             species          = species,  #"Mus musculus", "Homo sapiens"
                             contrasts        = contrasts,
                             parent.dir       = parent.dir,
                             gmt.dir          = gmt.dir,
                             scripts.dir      = scripts.dir,
                             deseq2.override  = deseq2.override,
                             heatmap.override = heatmap.override,
                             volcano.override = volcano.override)

# ---- PRE-ANALYSIS ----

meta_data <- read.xlsx(file.path(proj.params$proj.dir, paste0(proj, "_Metadata.xlsx")))
read_data <- NULL
trial <- TRUE
main_analysis(meta_data, read_data, proj.params, trial)

# ---- FULL-ANALYSIS ----

meta_data <- read.xlsx(file.path(proj.params$proj.dir, paste0(proj, "_Metadata.xlsx")))
read_data <- read.xlsx(file.path(proj.params$proj.dir, paste0(proj, "_Raw_counts.xlsx")))
trial <- FALSE

## Remove bad samples [CASE by CASE BASIS]
# read_data <- read_data %>%
#   dplyr::select(everything(), -c("SBQuadFc2", "SBQuadFc4"))
main_analysis(meta_data, read_data, proj.params, trial)