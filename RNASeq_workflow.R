#!/usr/bin/env Rscript

# Read and store variables from command line interface (CLI)
cli <- base::commandArgs(trailingOnly = TRUE) 
args <- base::strsplit(x = cli, split = "=", fixed = TRUE)

for (e in args){
  argname <- e[1]
  argval <- e[2]
  assign(argname, argval)
}

# ---- PROJECT SET UP ----

if (.Platform$OS.type == "windows") {
  parent_dir  <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data"
  gmt_dir     <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GSEA_genesets"
  scripts_dir <- NULL
  script_file <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/Custom_Functions.R"
} else {  # Linux/macOS (e.g., HPC)
  parent_dir  <- "/hpc/home/kailasamms/scratch"
  gmt_dir     <- "/hpc/home/kailasamms/projects/GSEA_genesets"
  scripts_dir <- "/hpc/home/kailasamms/projects/scRNASeq"
  script_file <- "/hpc/home/kailasamms/projects/scRNASeq/Custom_Functions.R"
}

if (file.exists(script_file)) {
  source(script_file)
} else{
  stop(paste("Custom_Functions.R not found at:", script_file))
}

# Define these 4 variables in project specific script
# contrasts <- c()
# deseq2.override <- list()
# heatmap.override <- list()
# volcano.override <- list()

proj.params <- setup_project(proj             = proj,
                             species          = species,  #"Mus musculus", "Homo sapiens"
                             contrasts        = contrasts,
                             parent_dir       = parent_dir,
                             gmt_dir          = gmt_dir,
                             scripts_dir      = scripts_dir,
                             deseq2.override  = deseq2.override,
                             heatmap.override = heatmap.override,
                             volcano.override = volcano.override)

# ---- PRE-ANALYSIS ----

meta_data <- read.xlsx(file.path(proj.params$proj_dir, paste0(proj, "_Metadata.xlsx")))
read_data <- NULL
trial <- TRUE
main_analysis(meta_data, read_data, proj.params, trial)

# ---- BULK RNA SEQ WORKFLOW ----

meta_data <- read.xlsx(file.path(proj.params$proj_dir, paste0(proj, "_Metadata.xlsx")))
read_data <- read.xlsx(file.path(proj.params$proj_dir, paste0(proj, "_Raw_counts.xlsx")))
trial <- FALSE

## Remove bad samples [CASE by CASE BASIS]
# read_data <- read_data %>%
#   dplyr::select(everything(), -c("SBQuadFc2", "SBQuadFc4"))
main_analysis(meta_data, read_data, proj.params, trial)