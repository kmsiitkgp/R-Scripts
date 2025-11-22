
proj.dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Collaboration projects data/Kamya_IMC"
# df1 <- openxlsx::read.xlsx(file.path(proj.dir, "Simian Files/name_mapping.xlsx"))
# df2 <- openxlsx::read.xlsx(file.path(proj.dir, "KS_LiverMets_Cohort2_clinicalinfo 08.15.25 - de-identified.xlsx"))
# df <- dplyr::full_join(df1, df2, by=c("PatientID"="PatientID"))


df <- openxlsx::read.xlsx(file.path(proj.dir, "Metadata.xlsx"))
  
df$DateBirth <- openxlsx::convertToDate(df$DateBirth)
df$DateDiagnosis <- openxlsx::convertToDate(df$DateDiagnosis)
df$DateSurgery <- openxlsx::convertToDate(df$DateSurgery)
df$DateDeath <- openxlsx::convertToDate(df$DateDeath)
df$DateLastFollowUp <- openxlsx::convertToDate(df$DateLastFollowUp)
df$DateImmunotherapyStart <- openxlsx::convertToDate(df$DateImmunotherapyStart)

df <- df %>% dplyr::mutate(OS = dplyr::case_when(!is.na(DateDeath) ~ DateDeath - DateDiagnosis,
                                                 !is.na(DateLastFollowUp) ~ DateLastFollowUp -DateDiagnosis,
                                                 TRUE ~ NA),
                           DPS = dplyr::case_when(!is.na(DateDeath) ~ DateDeath - DateSurgery,
                                                  !is.na(DateLastFollowUp) ~ DateLastFollowUp - DateSurgery,
                                                  TRUE ~ NA)) %>%
  dplyr::select(SlideName, DateName, SurgicalAccession, SampleID, PatientID, SlideID, 
                AnonymID, DateBirth, DateDeath, DateLastFollowUp, DateDiagnosis,
                DateSurgery, Status, OS, DPS, AgeDiagnosis, StageDiagnosis, 
                everything())


# wb <- createWorkbook()
# addWorksheet(wb, "Metadata")
# writeData(wb, sheet = "Metadata", df)
# saveWorkbook(wb, file.path(proj.dir, "Metadata.xlsx"), overwrite = TRUE)


# ---- Correlation ---- 
freq_subtype <- openxlsx::read.xlsx(file.path(proj.dir, "sample_counts.xlsx"))
freq_celltype <- openxlsx::read.xlsx(file.path(proj.dir, "sample_counts_meta.xlsx"))

freq_celltype <- freq_celltype %>%
  tidyr::separate(col = txt, 
                  into = c("DateName", "AnonymID", "SlideID", "ROI"), 
                  sep = "_", 
                  extra = "merge") %>%
  dplyr::left_join(clinical, by=c("DateName"="DateName",
                                  "AnonymID"="AnonymID")) %>%
  rowwise() %>%
  mutate(
    total_cells = sum(c_across(Endothelial:Tumor)),
    Endothelial = Endothelial / total_cells,
    Innate = Innate / total_cells,
    Lymphoid = Lymphoid / total_cells,
    Stromal = Stromal / total_cells,
    Tumor = Tumor / total_cells
  ) %>%
  ungroup() %>%
  select(-total_cells) %>%
  group_by(AnonymID) %>%
  mutate(
    Endothelial = mean(Endothelial, na.rm = TRUE),
    Innate      = mean(Innate, na.rm = TRUE),
    Lymphoid    = mean(Lymphoid, na.rm = TRUE),
    Stromal     = mean(Stromal, na.rm = TRUE),
    Tumor       = mean(Tumor, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  distinct_at("AnonymID", .keep_all = TRUE)
  

# Select only relevant columns
cell_cols <- c("Endothelial","Innate","Lymphoid","Stromal","Tumor")
clinical_cols <- c("Metastasis","Sex","Histology","ECOGPSDiagnosis","Smoking","Immunotherapy")

# Split clinical variables into numeric vs categorical
numeric_clinical <- clinical_cols[sapply(freq_celltype[clinical_cols], is.numeric)]
categorical_clinical <- setdiff(clinical_cols, numeric_clinical)

# 1) Numeric vars: Spearman correlation
results_num <- expand.grid(Clinical = numeric_clinical, CellType = cell_cols) %>%
  rowwise() %>%
  mutate(
    test = list(cor.test(df_corr[[Clinical]], df_corr[[CellType]],
                         method = "spearman", exact = FALSE,
                         use = "complete.obs")),
    rho  = test$estimate,
    pval = test$p.value,
    n    = sum(complete.cases(df_corr[[Clinical]], df_corr[[CellType]])),
    method = "Spearman"
  ) %>%
  select(-test) %>%
  ungroup()

# 2) Categorical vars: group tests
results_cat <- expand.grid(Clinical = categorical_clinical, CellType = cell_cols) %>%
  rowwise() %>%
  mutate(
    groups = list(df_corr[[Clinical]]),
    values = list(df_corr[[CellType]]),
    
    test = list({
      dat <- tibble(g = groups[[1]], v = values[[1]]) %>% 
        filter(!is.na(g), !is.na(v))
      ng <- n_distinct(dat$g)
      
      if (ng < 2) {
        NULL   # not testable
      } else if (ng == 2) {
        wilcox.test(v ~ g, data = dat)
      } else {
        kruskal.test(v ~ g, data = dat)
      }
    }),
    
    rho  = NA_real_,
    pval = if (!is.null(test)) test$p.value else NA_real_,
    n    = sum(complete.cases(groups[[1]], values[[1]])),
    method = case_when(
      is.null(test)            ~ "Not testable",
      n_distinct(na.omit(groups[[1]])) == 2 ~ "Wilcoxon",
      n_distinct(na.omit(groups[[1]])) > 2  ~ "Kruskal–Wallis"
    )
  ) %>%
  select(-test, -groups, -values) %>%
  ungroup()

# 3) Combine & adjust p-values
results <- bind_rows(results_num, results_cat) %>%
  mutate(padj = p.adjust(pval, method = "fdr")) %>%
  arrange(pval)

results %>% filter(padj <= 0.05)


# ---- Survival ----

# remove lung samples, OS/DPS <= 30 days
meta_data <- df %>% 
  dplyr::mutate(Sample_ID = as.character(PatientID), 
                OS_Time = as.numeric(OS)/30.44, DPS_Time = as.numeric(DPS)/30.44) %>%
  dplyr::filter(!is.na(Status), !is.na(SlideName), SourceSite == "Liver", (OS > 30 | DPS > 30)) %>%
  distinct(Sample_ID, .keep_all = TRUE)

# Keep only samples retained in meta_data  
expr <- openxlsx::read.xlsx(file.path(proj.dir, "cell_norm_expression.xlsx")) %>%
  tidyr::separate(col = txt, 
                  into = c("DateName", "AnonymID", "SlideID", "ROI"), 
                  sep = "_", 
                  extra = "merge")

expr <- expr %>%
  dplyr::filter(AnonymID %in% meta_data$AnonymID) %>%
  dplyr::left_join(meta_data %>% dplyr::select(AnonymID, PatientID) %>% dplyr::mutate(AnonymID = as.character(AnonymID)), 
                   by=c("AnonymID"="AnonymID"))

marker_cols <- colnames(expr)[16:57]

# Merge counts for retained patients
expr_metapheno <- expr %>% 
  dplyr::group_by(PatientID, metapheno) %>%
  dplyr::summarise(across(all_of(marker_cols), mean, na.rm = TRUE),.groups = "drop")
metapheno_groups <- unique(expr_metapheno$metapheno)

# Merge counts for retained patients
expr_pheno <- expr %>% 
  dplyr::group_by(PatientID, pheno) %>%
  dplyr::summarise(across(all_of(marker_cols), mean, na.rm = TRUE),.groups = "drop")
pheno_groups <- unique(expr_pheno$pheno)

# Merge counts for retained patients
expr_custom <- expr %>% 
  dplyr::filter(pheno %in% c("T helper", "T cytotoxic")) %>%
  dplyr::group_by(PatientID) %>%
  dplyr::summarise(across(all_of(marker_cols), mean, na.rm = TRUE),.groups = "drop") %>%
  dplyr::mutate(pheno = "Tcell")
custom_groups <- unique(expr_custom$pheno)
  
# Prepare a list of datasets with their grouping column and output folder prefix
grouped_data <- list(
  list(data = expr_metapheno, group_col = "metapheno", output_dir = "metapheno"),
  list(data = expr_pheno,     group_col = "pheno",     output_dir = "pheno"),
  list(data = expr_custom,    group_col = "pheno",     output_dir = "custom")
)

for (item in grouped_data) {
  # Extract components
  df <- item$data
  group_col <- item$group_col
  output_dir <- item$output_dir
  
  # Get unique groups for this dataset
  groups <- unique(df[[group_col]])
  
  for (group in groups) {
    expr_data <- df %>%
      dplyr::filter(!!rlang::sym(group_col) == group) %>%
      dplyr::select(-!!rlang::sym(group_col)) %>%
      tibble::column_to_rownames("PatientID") %>%
      t()
    
    survival_params <- list(
    # ---- Stratification (Expression + Metadata-based survival) ----
    stratify_var     = marker_cols,          # one or more genes or metadata columns
    sig_score        = FALSE,          # TRUE = combine genes into one signature score
    substratify_var  = NULL,          # optional metadata column for sub-stratification
    facet_var        = NULL,          # optional faceting variable
    
    # ---- Cutoff settings (ONLY for Expression-based survival) ----
    cutoff_method    = "optimal",      # median, quartile, tertile, optimal, thirds
    show_all_bins    = FALSE,          # TRUE = plot all bins (LOW, HIGH, MID/MED_HIGH/MED_LOW)
    multiple_cutoff  = FALSE,          # TRUE = compute cutoffs separately for substratify_var
    
    # ---- Plot settings ----
    conf_interval    = FALSE,          # TRUE = show confidence interval in survival curve
    plot_curve       = TRUE,           # TRUE = plot the survival curve
    plot_risk_table  = TRUE,           # TRUE = plot the risk table below the curve
    color_palette    = custom_palette, # vector of colors for groups c("#d73027","#0c2c84")
    
    # ---- Survival data columns ----
    time_col         = "OS_Time", #"DPS_Time",        # metadata column containing Time values
    status_col       = "Status",       # metadata column containing Status values
    
    # ---- Output ----
    prefix           = "",
    output_path      = file.path("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop", output_dir, group)
    )
    
    survival_analysis(meta_data, expr_data, survival_params)
  }
}

# Read all excel file and combine the results
base_dir <- "C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/OS"
categories <- c("metapheno", "pheno", "custom")

merged_data <- map_dfr(categories, function(cat) {
  cat_path <- file.path(base_dir, cat)
  
  # Get all subdirectories inside this category
  subdirs <- list.dirs(cat_path, recursive = FALSE, full.names = TRUE)
  
  # For each subdir, read the one Excel file, add source column
  map_dfr(subdirs, function(subdir) {
    excel_file <- list.files(subdir, pattern = "\\.xlsx?$", full.names = TRUE)[1]
    
    if (is.na(excel_file)) return(NULL)
    
    df <- openxlsx::read.xlsx(excel_file)
    
    # Add source column like "metapheno_tcell"
    source_name <- paste0(cat, "_", basename(subdir))
    df$source <- source_name
    
    return(df)
  })
})

# Save merged data
output_file <- file.path(base_dir, "merged_all.xlsx")
openxlsx::write.xlsx(merged_data, output_file)

message("Merged data saved at: ", output_file)



