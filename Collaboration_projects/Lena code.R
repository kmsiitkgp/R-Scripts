source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

parent_path <- "C:/Users/kailasamms/Box/Saravana@cedars/05. Bioinformatics/CRISPR/CRISPR_Lena/"
count_path <- paste0(parent_path, "count_results/")
DEG_path <- paste0(parent_path, "DEG_results/")
  
meta_data <- read.xlsx(paste0(parent_path, "CRISPR_Lena_Meta_data.xlsx")) %>%
  dplyr::mutate(Sample_ID = make.names(Sample_ID)) %>%
  dplyr::filter(Original_Library == Actual_Library)

#### Calculate scores from count files
for (batch in unique(meta_data$Batch)){
  
  meta_data_sub <- meta_data %>% 
    dplyr::filter(Batch == batch) %>% 
    dplyr::select(Sample_ID, Cell_line, Batch, Actual_Library)
  
  # library
  library <- unique(meta_data_sub$Actual_Library)
  
  #cell line
  cell_line <- unique(meta_data_sub$Cell_line)
  
  # Read the counts file
  counts <- read.table(paste0(count_path, library, "/mageck.total.count.txt"), header=TRUE) %>%
    dplyr::select(sgRNA, Gene, meta_data_sub$Sample_ID)
  
  # If first cell has "Byte-Order-Mark" (BOM), correct it
  #counts[1,1] <- stringr::str_replace(count_table[1,1], "ï»¿","")
  
  # control samples
  control_samples <- colnames(counts)[grepl(pattern="DO", x=colnames(counts))]
  
  # treatment samples
  treatment_samples <- setdiff(make.names(meta_data_sub$Sample_ID), control_samples)
  
  cat(control_samples)
  cat(treatment_samples)
  
  # Find Seq depth of each sample
  meta_data_sub$sample_depth <- as.numeric(colSums(counts[,-c(1,2)]))
  
  # sgRNAs that have 10 or less raw reads are most likely absent
  meta_data_sub <- meta_data_sub %>% 
    dplyr::mutate(sample_specific_cutoff = 10/sample_depth*100)
  
  # Remove sgRNAs from Day 0 samples that have less than 10 reads as these were
  # absent at t=0
  day_zero_cutoff <- 10/(2*10^6)*100
  
  # We will only focus on sgRNAs that had at least 100 or 200 reads at t=0 for depletion hits
  depletion_cutoff <- 100/(2*10^6)*100
  #depletion_cutoff <- 200/(2*10^6)*100
  
  # Normalized counts
  norm_counts <- counts %>%
    dplyr::mutate(across(.cols=c(everything(), -c(sgRNA,Gene)), ~ (100*.)/sum(.))) %>%
    dplyr::mutate(across(.cols=c(everything(), -c(sgRNA,Gene)), ~ dplyr::case_when(get(control_samples) >= day_zero_cutoff ~ .,
                                                                                   TRUE ~ NA)))
  # Score matrix
  score_matrix <- norm_counts
  
  for (row in 1:nrow(norm_counts)){
    for (sample in treatment_samples){
      
      cutoff <- meta_data_sub %>% 
        dplyr::filter(Sample_ID == sample) %>% 
        dplyr::select(sample_specific_cutoff) %>%
        unlist(use.names=FALSE)

      # Score based on expression of sgRNA
      # If norm count <= sample specific cutoff, then 0 (surely depleted)
      # If sample specific cutoff < norm count < Day 0, then 1 (possibly depleted)
      # If Day 0 <= norm_count, then 2 (possibly enriched) 
      # OLD SCORING:
      # if (is.na(norm_counts[row,sample])){
      #   score_matrix[row,sample] <- NA
      # } else if(norm_counts[row,sample] <= cutoff){
      #   score_matrix[row,sample] <- 0
      # } else if (norm_counts[row,sample] < norm_counts[row,control_samples]){
      #   score_matrix[row,sample] <- 1
      # } else{
      #   score_matrix[row,sample] <- 2
      # }
      
      # NEW SCORING:
      if (is.na(norm_counts[row,sample])){
        score_matrix[row,sample] <- NA
      } else if(norm_counts[row,sample] <= cutoff){
        score_matrix[row,sample] <- 0
      } else if (norm_counts[row,sample] <= norm_counts[row,control_samples]/2){
        score_matrix[row,sample] <- 1
      } else if (norm_counts[row,sample] <= norm_counts[row,control_samples]){
        score_matrix[row,sample] <- 2
      } else if (norm_counts[row,sample] <= norm_counts[row,control_samples]*3/2){
        score_matrix[row,sample] <- 3
      } else{
        score_matrix[row,sample] <- 4
      }
      
    }
  }
  
  score_matrix$sgRNA_score <- rowMeans(score_matrix[,treatment_samples])
  #score_matrix$sgRNA_score <- rowMedians(as.matrix(score_matrix[,treatment_samples]))
  #score_matrix$sgRNA_score <- rowSums(score_matrix[,treatment_samples])
  
  # Enrichment score and n
  score_matrix1 <- score_matrix %>%
    dplyr::filter(!is.na(get(control_samples))) %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(n=n(),
                  Enrichment_Score = mean(sgRNA_score, na.rm=TRUE),
                  #Enrichment_Score = median(sgRNA_score, na.rm=TRUE),
                  #Enrichment_Score = sum(sgRNA_score, na.rm=TRUE)/n,
                  Enrichment_n=length(treatment_samples)*n) %>%
    ungroup() %>%
    dplyr::select(sgRNA, Enrichment_Score, Enrichment_n)
  
  # Depletion score and n
  score_matrix2 <- score_matrix %>%
    dplyr::filter(!is.na(get(control_samples))) %>%
    dplyr::filter(get(control_samples) > depletion_cutoff) %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(n=n(),
                  Depletion_Score = mean(sgRNA_score, na.rm=TRUE),
                  #Depletion_Score = median(sgRNA_score, na.rm=TRUE),
                  #Depletion_Score = sum(sgRNA_score, na.rm=TRUE)/n,
                  Depletion_n=length(treatment_samples)*n) %>%
    ungroup() %>%
    dplyr::select(sgRNA, Depletion_Score, Depletion_n)
  
  score_matrix <- score_matrix %>%
    dplyr::left_join(score_matrix1, by=c("sgRNA"="sgRNA")) %>%
    dplyr::left_join(score_matrix2, by=c("sgRNA"="sgRNA"))
  
  # Identify cutoffs
  control_matrix <- score_matrix %>%
    dplyr::filter(Gene == "negative_control" | grepl(pattern="Olfr",x=Gene))
  
  treatment_matrix <- score_matrix[!(score_matrix$sgRNA %in% control_matrix$sgRNA),]
  
  cutoff_depletion <- mean(control_matrix$Depletion_Score, na.rm=TRUE)
  control_n <- max(control_matrix$Enrichment_n, na.rm=TRUE)
  treatment_n <-  max(treatment_matrix$Enrichment_n, na.rm=TRUE)
 
  cat(library, ":", cell_line, ";", treatment_n, ":", control_n, "\n")
  
  #Final score and n
  score_matrix <- score_matrix %>%
    dplyr::mutate(Cutoff = cutoff_depletion,
                  SD = sd(Depletion_Score, na.rm=TRUE),
                  Score = dplyr::case_when(Depletion_Score < Cutoff ~ Depletion_Score,
                                           Depletion_Score >= Cutoff ~ Enrichment_Score),
                  Score = Score - cutoff_depletion,
                  n = dplyr::case_when(Depletion_Score < cutoff_depletion ~ Depletion_n,
                                       Depletion_Score >= cutoff_depletion ~ Enrichment_n),
                  Percent_n = dplyr::case_when(Gene == "negative_control" | grepl(pattern="Olfr",x=Gene) ~ n*100/control_n,
                                       TRUE ~ n*100/treatment_n),
                  Hits = dplyr::case_when(Score < -2*SD & Percent_n >= 60 ~ "Depletion 2SD",
                                          Score < -SD & Percent_n >= 60 ~ "Depletion 1SD",
                                          Score > 2*SD & Percent_n >= 60 ~ "Enrichment 2SD",
                                          Score >  SD & Percent_n >= 60 ~ "Enrichment 1SD",
                                          TRUE ~ "Not Significant")) %>%
    dplyr::arrange(sgRNA) %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(Score = median(Score,na.rm=TRUE), n=median(n, na.rm=TRUE), Percent_n=median(Percent_n,na.rm=TRUE)) %>%
    dplyr::ungroup()
  
  # Save output in excel
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Scores")
  openxlsx::writeData(wb, sheet = "Scores", x = score_matrix)
  openxlsx::addWorksheet(wb, sheetName = "Normalized")
  openxlsx::writeData(wb, sheet = "Normalized", x = norm_counts)
  openxlsx::addWorksheet(wb, sheetName = "Cutoffs")
  openxlsx::writeData(wb, sheet = "Cutoffs", x = meta_data_sub)
  #openxlsx::saveWorkbook(wb, file = paste0(parent_path,"Scores_Mean100/", cell_line, "_", library, "_Batch_", batch , ".xlsx"), overwrite=TRUE)
  #openxlsx::saveWorkbook(wb, file = paste0(parent_path,"Scores_Median100/", cell_line, "_", library, "_Batch_", batch , ".xlsx"), overwrite=TRUE)
  #openxlsx::saveWorkbook(wb, file = paste0(parent_path,"Scores_Sum100/", cell_line, "_", library, "_Batch_", batch , ".xlsx"), overwrite=TRUE)
  openxlsx::saveWorkbook(wb, file = paste0(parent_path,"Scores_New100/", cell_line, "_", library, "_Batch_", batch , ".xlsx"), overwrite=TRUE)
}

################# Merge batch 6 and 7 mean/new

data1 <- read.xlsx(paste0(parent_path, "Scores_New100/pre-merge/MB83F_CDH12ACT_Batch_6.xlsx")) %>%
  dplyr::select(everything(), -c(Score, n, Percent_n,SD,Hits)) %>%
  dplyr::rename(sgRNA_score_6 = sgRNA_score, 
                Enrichment_Score_6 = Enrichment_Score,
                Depletion_Score_6 = Depletion_Score,
                Cutoff_6 = Cutoff,
                Enrichment_n_6 = Enrichment_n,
                Depletion_n_6 = Depletion_n)
data2 <- read.xlsx(paste0(parent_path, "Scores_New100/pre-merge/MB83F_CDH12ACT_Batch_7.xlsx")) %>%
  dplyr::select(everything(), -c(Score, n, Percent_n,SD,Hits,Gene)) %>%
  dplyr::rename(sgRNA_score_7 = sgRNA_score, 
                Enrichment_Score_7 = Enrichment_Score,
                Depletion_Score_7 = Depletion_Score,
                Cutoff_7 = Cutoff,
                Enrichment_n_7 = Enrichment_n,
                Depletion_n_7 = Depletion_n)

score_matrix <- data1 %>% 
  dplyr::left_join(data2, by=c("sgRNA"="sgRNA")) %>%
  rowwise() %>%
  dplyr::mutate(Enrichment_Score = mean(c_across(starts_with("Enrichment_Score")), na.rm=TRUE),
                Depletion_Score = mean(c_across(starts_with("Depletion_Score")), na.rm=TRUE),
                Enrichment_n = sum(c_across(starts_with("Enrichment_n")), na.rm=TRUE),
                Depletion_n = sum(c_across(starts_with("Depletion_n")), na.rm=TRUE)) %>%
  ungroup()

control_matrix <- score_matrix %>%
  dplyr::filter(Gene == "negative_control" | grepl(pattern="Olfr",x=Gene))

treatment_matrix <- score_matrix[!(score_matrix$sgRNA %in% control_matrix$sgRNA),]

cutoff_depletion <- mean(control_matrix$Depletion_Score, na.rm=TRUE)
control_n <- max(control_matrix$Enrichment_n, na.rm=TRUE)
treatment_n <-  max(treatment_matrix$Enrichment_n, na.rm=TRUE)

score_matrix <- score_matrix %>%
  dplyr::mutate(Cutoff = cutoff_depletion,
                SD = sd(Depletion_Score, na.rm=TRUE),
                Score = dplyr::case_when(Depletion_Score < Cutoff ~ Depletion_Score,
                                         Depletion_Score >= Cutoff ~ Enrichment_Score),
                Score = Score - Cutoff,
                n = dplyr::case_when(Depletion_Score < Cutoff ~ Depletion_n,
                                     Depletion_Score >= Cutoff ~ Enrichment_n),
                Percent_n = dplyr::case_when(Gene == "negative_control" | grepl(pattern="Olfr",x=Gene) ~ n*100/control_n,
                                             TRUE ~ n*100/treatment_n),
                Hits = dplyr::case_when(Score < -2*SD & Percent_n >= 60 ~ "Depletion 2SD",
                                        Score < -SD & Percent_n >= 60 ~ "Depletion 1SD",
                                        Score > 2*SD & Percent_n >= 60 ~ "Enrichment 2SD",
                                        Score >  SD & Percent_n >= 60 ~ "Enrichment 1SD",
                                        TRUE ~ "Not Significant")) %>%
  dplyr::arrange(sgRNA) %>%
  dplyr::group_by(Gene) %>%
  dplyr::mutate(Score = median(Score,na.rm=TRUE), n=median(n, na.rm=TRUE), Percent_n=median(Percent_n,na.rm=TRUE)) %>%
  dplyr::ungroup()

# Save output in excel
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Scores")
openxlsx::writeData(wb, sheet = "Scores", x = score_matrix)
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Scores_New100/MB83F_CDH12ACT_Batch_6+7.xlsx"), overwrite=TRUE)

################# Merge batch 6 and 7 median

data1 <- read.xlsx(paste0(parent_path, "Scores_Median200/pre-merge/MB83F_CDH12ACT_Batch_6.xlsx")) %>%
  dplyr::select(sgRNA, Gene, contains("PCR"), Enrichment_n, Depletion_n) %>%
  dplyr::rename(Enrichment_n_6 = Enrichment_n,
                Depletion_n_6 = Depletion_n)

data2 <- read.xlsx(paste0(parent_path, "Scores_Median200/pre-merge/MB83F_CDH12ACT_Batch_7.xlsx")) %>%
  dplyr::select(sgRNA, Gene, contains("PCR"),  Enrichment_n, Depletion_n) %>%
  dplyr::rename(Enrichment_n_7 = Enrichment_n,
                Depletion_n_7 = Depletion_n)

score_matrix <- data1 %>% 
  dplyr::left_join(data2, by=c("sgRNA"="sgRNA", "Gene"="Gene"))

score_matrix$sgRNA_score <- matrixStats::rowMedians(as.matrix(score_matrix[,c("PCR.46", "PCR.47", "PCR.48", "PCR.50", "PCR.51")]))

score_matrix <- score_matrix %>%
  rowwise %>%
  dplyr::mutate(Enrichment_n = sum(c_across(starts_with("Enrichment_n")), na.rm=TRUE),
                Depletion_n = sum(c_across(starts_with("Depletion_n")), na.rm=TRUE)) 

# control samples
control_samples <- c("PCR.06DO", "PCR.06.2DO")

# treatment samples
treatment_samples <- c("PCR.46", "PCR.47", "PCR.48", "PCR.50", "PCR.51")

depletion_cutoff <- 100/(2*10^6)*100
depletion_cutoff <- 200/(2*10^6)*100

# Enrichment score and n
score_matrix1 <- score_matrix %>%
  dplyr::group_by(Gene) %>%
  dplyr::mutate(Enrichment_Score = median(sgRNA_score, na.rm=TRUE)) %>%
  ungroup() %>%
  dplyr::select(sgRNA, Enrichment_Score)

# Depletion score and n
score_matrix2 <- score_matrix %>%
  dplyr::group_by(Gene) %>%
  dplyr::mutate(Depletion_Score = median(sgRNA_score, na.rm=TRUE)) %>%
  ungroup() %>%
  dplyr::select(sgRNA, Depletion_Score)

score_matrix <- score_matrix %>%
  dplyr::left_join(score_matrix1, by=c("sgRNA"="sgRNA")) %>%
  dplyr::left_join(score_matrix2, by=c("sgRNA"="sgRNA"))

# Identify cutoffs
control_matrix <- score_matrix %>%
  dplyr::filter(Gene == "negative_control" | grepl(pattern="Olfr",x=Gene))

treatment_matrix <- score_matrix[!(score_matrix$sgRNA %in% control_matrix$sgRNA),]

cutoff_depletion <- mean(control_matrix$Depletion_Score, na.rm=TRUE)
control_n <- max(control_matrix$Enrichment_n, na.rm=TRUE)
treatment_n <-  max(treatment_matrix$Enrichment_n, na.rm=TRUE)

#Final score and n
score_matrix <- score_matrix %>%
  dplyr::mutate(Cutoff = cutoff_depletion,
                SD = sd(Depletion_Score, na.rm=TRUE),
                Score = dplyr::case_when(Depletion_Score < Cutoff ~ Depletion_Score,
                                         Depletion_Score >= Cutoff ~ Enrichment_Score),
                Score = Score - cutoff_depletion,
                n = dplyr::case_when(Depletion_Score < cutoff_depletion ~ Depletion_n,
                                     Depletion_Score >= cutoff_depletion ~ Enrichment_n),
                Percent_n = dplyr::case_when(Gene == "negative_control" | grepl(pattern="Olfr",x=Gene) ~ n*100/control_n,
                                             TRUE ~ n*100/treatment_n),
                Hits = dplyr::case_when(Score < -2*SD & Percent_n >= 60 ~ "Depletion 2SD",
                                        Score < -SD & Percent_n >= 60 ~ "Depletion 1SD",
                                        Score > 2*SD & Percent_n >= 60 ~ "Enrichment 2SD",
                                        Score >  SD & Percent_n >= 60 ~ "Enrichment 1SD",
                                        TRUE ~ "Not Significant")) %>%
  dplyr::arrange(sgRNA) %>%
  dplyr::group_by(Gene) %>%
  dplyr::mutate(Score = median(Score,na.rm=TRUE), n=median(n, na.rm=TRUE), Percent_n=median(Percent_n,na.rm=TRUE)) %>%
  dplyr::ungroup()

# Save output in excel
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Scores")
openxlsx::writeData(wb, sheet = "Scores", x = score_matrix)
openxlsx::saveWorkbook(wb, file = paste0(parent_path,"Scores_Median200/MB83F_CDH12ACT_Batch_6+7.xlsx"), overwrite=TRUE)

################# Merge batch 6 and 7 Sum

data1 <- read.xlsx(paste0(parent_path, "Scores_Sum100/pre-merge/MB83F_CDH12ACT_Batch_6.xlsx")) %>%
  dplyr::select(everything(), -c(Score, n, Percent_n,SD,Hits)) %>%
  dplyr::rename(sgRNA_score_6 = sgRNA_score, 
                Enrichment_Score_6 = Enrichment_Score,
                Depletion_Score_6 = Depletion_Score,
                Cutoff_6 = Cutoff,
                Enrichment_n_6 = Enrichment_n,
                Depletion_n_6 = Depletion_n)
data2 <- read.xlsx(paste0(parent_path, "Scores_Sum100/pre-merge/MB83F_CDH12ACT_Batch_7.xlsx")) %>%
  dplyr::select(everything(), -c(Score, n, Percent_n,SD,Hits,Gene)) %>%
  dplyr::rename(sgRNA_score_7 = sgRNA_score, 
                Enrichment_Score_7 = Enrichment_Score,
                Depletion_Score_7 = Depletion_Score,
                Cutoff_7 = Cutoff,
                Enrichment_n_7 = Enrichment_n,
                Depletion_n_7 = Depletion_n)

score_matrix <- data1 %>% 
  dplyr::left_join(data2, by=c("sgRNA"="sgRNA")) %>%
  rowwise() %>%
  dplyr::mutate(Enrichment_Score = mean(c_across(starts_with("Enrichment_Score")), na.rm=TRUE),
                Depletion_Score = mean(c_across(starts_with("Depletion_Score")), na.rm=TRUE),
                Enrichment_n = sum(c_across(starts_with("Enrichment_n")), na.rm=TRUE),
                Depletion_n = sum(c_across(starts_with("Depletion_n")), na.rm=TRUE)) %>%
  ungroup()

control_matrix <- score_matrix %>%
  dplyr::filter(Gene == "negative_control" | grepl(pattern="Olfr",x=Gene))

treatment_matrix <- score_matrix[!(score_matrix$sgRNA %in% control_matrix$sgRNA),]

cutoff_depletion <- mean(control_matrix$Depletion_Score, na.rm=TRUE)
control_n <- max(control_matrix$Enrichment_n, na.rm=TRUE)
treatment_n <-  max(treatment_matrix$Enrichment_n, na.rm=TRUE)

score_matrix <- score_matrix %>%
  dplyr::mutate(Cutoff = cutoff_depletion,
                SD = sd(Depletion_Score, na.rm=TRUE),
                Score = dplyr::case_when(Depletion_Score < Cutoff ~ Depletion_Score,
                                         Depletion_Score >= Cutoff ~ Enrichment_Score),
                Score = Score - Cutoff,
                n = dplyr::case_when(Depletion_Score < Cutoff ~ Depletion_n,
                                     Depletion_Score >= Cutoff ~ Enrichment_n),
                Percent_n = dplyr::case_when(Gene == "negative_control" | grepl(pattern="Olfr",x=Gene) ~ n*100/control_n,
                                             TRUE ~ n*100/treatment_n),
                Hits = dplyr::case_when(Score < -2*SD & Percent_n >= 60 ~ "Depletion 2SD",
                                        Score < -SD & Percent_n >= 60 ~ "Depletion 1SD",
                                        Score > 2*SD & Percent_n >= 60 ~ "Enrichment 2SD",
                                        Score >  SD & Percent_n >= 60 ~ "Enrichment 1SD",
                                        TRUE ~ "Not Significant")) %>%
  dplyr::arrange(sgRNA) %>%
  dplyr::group_by(Gene) %>%
  dplyr::mutate(Score = median(Score,na.rm=TRUE), n=median(n, na.rm=TRUE), Percent_n=median(Percent_n,na.rm=TRUE)) %>%
  dplyr::ungroup()

# Save output in excel
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Scores")
openxlsx::writeData(wb, sheet = "Scores", x = score_matrix)
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Scores_Sum100/MB83F_CDH12ACT_Batch_6+7.xlsx"), overwrite=TRUE)


#######MB49 MAGECK mean, total, control, median, alphamedian, secondbest comparisons

# Normalization: 
# Median (with control sgRNAs), 
# Median (without control sgRNAs),
# total (with control sgRNAs), 
# total (without control sgRNAs),
# control (with control sgRNAs)
# For each of the above normalizations,  LFC was calculated using 3 methods: 
# median, aplhamedian, secondbest
# We now want to see which approach gives maximum overlap with our manual method
libraries <- c("CDH12KO", "CDH12ACT", "ImmuneKO", "ImmuneACT")

for (library in libraries){

  files <- list.files(paste0(DEG_path, library, "/"), full.names = TRUE)
  files <- files[grepl(pattern="MB49",x=files)]
  
  p_removed <- gsub(pattern=paste0(DEG_path,library, "/"), replacement="",x=files)
  cell_lines <- gsub(pattern="\\..*", replacement="",x=p_removed)
  norm_methods <- gsub(pattern=".*mageck.", replacement="",x=p_removed)
  norm_methods <- gsub(pattern="..txt", replacement="",x=norm_methods)
  
  final_l <- list()
  for (i in 1:length(files)){
    data <- read.table(files[i], header=TRUE)
    
    # Identify cutoffs
    cutoff <- mean(data %>%
                     dplyr::filter(id == "negative_control" | grepl(pattern="Olfr",x=id)) %>%
                     dplyr::select(neg.score) %>%
                     unlist(use.names=FALSE), na.rm=TRUE)
    
    stdev <- sd(data$neg.score, na.rm=TRUE)
    cat(cell_lines[i], ":\t", library, ":\t", norm_methods[i], ":\t", cutoff, ":\t", stdev, "\n")
    
    # Depletion hits
    dl <- data %>% 
      dplyr::filter(neg.score < cutoff-stdev) %>%
      dplyr::select(id) %>%
      unlist(use.names=FALSE) %>%
      list()
    names(dl) <- paste0(library, ".", norm_methods[i], ".depletion")
    
    # Enrichment hits
    el <- data %>% 
      dplyr::filter(neg.score > cutoff+stdev) %>%
      dplyr::select(id) %>%
      unlist(use.names=FALSE) %>%
      list()
    names(el) <- paste0(library, ".", norm_methods[i], ".enrich")
    
    final_l <- c(final_l, dl, el)
  }
  
  max_l <- max(lengths(final_l))
  final_l <- lapply(X=final_l, FUN=function(x) { c(x, rep(x="", times=max_l-length(x)))})
  final_l <- data.frame(final_l)
  
  # Save output in excel
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Scores")
  openxlsx::writeData(wb, sheet = "Scores", x = final_l)
  openxlsx::saveWorkbook(wb, file = paste0(parent_path, "MB49c5", library, "hits.xlsx"), overwrite=TRUE)

}

# ######################## Compile the hits
# 
# parent_path <- "C:/Users/kailasamms/Box/Saravana@cedars/05. Bioinformatics/CRISPR/CRISPR_Lena/"
# count_path <- paste0(parent_path, "count_results/")
# 
# meta_data <- read.xlsx(paste0(parent_path, "CRISPR_Lena_Meta_data.xlsx")) %>%
#   dplyr::filter(Original_Library == Actual_Library) %>%
#   dplyr::mutate(Sample_ID = make.names(Sample_ID))
# 
# depletion_df <- c()
# enrichment_df <- c()
# for (batch in unique(meta_data$Batch)){
#   
#   meta_data_sub <- meta_data %>% 
#     dplyr::filter(Batch == batch) %>% 
#     dplyr::select(Sample_ID, Cell_line, Batch, Actual_Library)
#   
#   # library
#   library <- unique(meta_data_sub$Actual_Library)
#   
#   # Read the counts file
#   counts <- read.table(paste0(count_path, library, "/mageck.total.count.txt"), header=TRUE) %>%
#     dplyr::select(sgRNA, Gene, meta_data_sub$Sample_ID)
#   
#   depletion_df <- c(depletion_df, counts$Gene)
#   enrichment_df <- c(enrichment_df, counts$Gene)
# }
# 
# depletion_df <- unique(depletion_df) %>% data.frame()
# enrichment_df <- unique(enrichment_df) %>% data.frame()
# colnames(depletion_df) <- "Gene"
# colnames(enrichment_df) <- "Gene"
# 
# files <- list.files(paste0(parent_path, "Scores/"), full.names = TRUE)
# files <- files[grepl(pattern=".xlsx",x=files)]
# for (f in files){
#   
#   score_matrix <- read.xlsx(f)
#   filename <- gsub(pattern=paste0(parent_path,"Scores/"), replacement="",x=f)
#   cell_line <- gsub(pattern="_.*", replacement="",x=filename)
#   library <- gsub(pattern=paste0(cell_line, "_"), replacement="",x=filename)
#   library <- gsub(pattern="_.*", replacement="",x=library)
#   
#   depletion_hits <- score_matrix %>%
#     dplyr::filter(Hits %in% c("Depletion 2SD", "Depletion 1SD")) %>%
#     dplyr::select(Gene) %>%
#     dplyr::distinct_at("Gene", .keep_all = TRUE) %>%
#     dplyr::mutate(!!paste0(cell_line, ".", library) := Gene)
#   
#   enrichment_hits <- score_matrix %>%
#     dplyr::filter(Hits %in% c("Enrichment 2SD", "Enrichment 1SD")) %>%
#     dplyr::select(Gene) %>%
#     dplyr::distinct_at("Gene", .keep_all = TRUE) %>%
#     dplyr::mutate(!!paste0(cell_line, ".", library) := Gene)
#   
#   depletion_df <- depletion_df %>%
#     dplyr::left_join(depletion_hits,by=c("Gene"="Gene"))
#   
#   enrichment_df <- enrichment_df %>%
#     dplyr::left_join(enrichment_hits,by=c("Gene"="Gene"))
#   
# }
# 
# #Save output in excel
# wb <- openxlsx::createWorkbook()
# openxlsx::addWorksheet(wb, sheetName = "Depletion_Hits")
# openxlsx::writeData(wb, sheet = "Depletion_Hits", x = depletion_df)
# openxlsx::addWorksheet(wb, sheetName = "Enrichment_Hits")
# openxlsx::writeData(wb, sheet = "Enrichment_Hits", x = enrichment_df)
# openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Hits.xlsx"), overwrite=TRUE)


#######  Plot PCA

# Read all count files
df <- read.table(paste0(count_path, "CDH12KO/mageck.total.count.txt"), header=TRUE) %>%
  dplyr::mutate(Gene = make.names(Gene, unique=TRUE)) %>%
  dplyr::select("Gene") 

libraries <- c("CDH12KO", "CDH12ACT", "ImmuneKO", "ImmuneACT")
for (library in libraries){
  
  # Read the counts file
  counts <- read.table(paste0(count_path, library, "/mageck.total.count.txt"), header=TRUE) %>%
    dplyr::mutate(Gene = make.names(Gene, unique=TRUE)) %>%
    dplyr::select(everything(), -sgRNA)
  
  # df <- df %>% 
  #   full_join(counts, by=c("Gene"="Gene"))
  df <- counts
  
  # Normalized counts
  norm_counts <- df %>%
    dplyr::mutate(across(.cols=c(everything(), -c(Gene)), ~ (100*.)/sum(., na.rm=TRUE)))
  
  norm_counts[is.na(norm_counts)] <- 0
  norm_counts <- norm_counts[rowSums(norm_counts[,-1]) != 0,]
  
  norm_counts <- norm_counts %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("Gene") %>%
    t() %>%
    data.frame() %>%
    tibble::rownames_to_column("Sample_ID")
  
  pca <- prcomp(norm_counts[,-1], 
                center = TRUE, 
                scale. = TRUE)
  library("factoextra")
  fviz_eig(pca, addlabels = TRUE)
  ggsave(paste0(library, "_Scree.jpeg"),
         height = 11,
         width = 8)
  
  df <- pca$x %>%
    data.frame()
  
  rownames(df) <- norm_counts$Sample_ID
  
  df <- df %>%
    tibble::rownames_to_column("Sample_ID") %>%
    dplyr::left_join(meta_data %>% 
                       dplyr::mutate(Cell_line = dplyr::if_else(Original_Library != Actual_Library, "Excluded", Cell_line)) %>%
                       dplyr::select(Sample_ID, Cell_line, Original_Library),
                     by=c("Sample_ID"="Sample_ID")) %>%
    dplyr::select(Sample_ID, Cell_line, Original_Library, everything())
  
  ggplot(df, aes(x=PC1, y=PC2, color=Cell_line, shape=Original_Library, label = Sample_ID)) + 
    geom_point() + 
    geom_text_repel()
  
  ggsave(paste0(library, "_PCA.jpeg"),
         height = 11,
         width = 8)
  # geom_text_repel(label=df$Group)
}

########## Plot hits for each library

libraries <- c("CDH12KO", "CDH12ACT", "ImmuneKO", "ImmuneACT")
for (library in libraries){

  files <- list.files(paste0(parent_path, "Scores/"), full.names = TRUE)
  files <- files[grepl(pattern=".xlsx", x= files)]
  files <- files[grepl(pattern=library, x= files)]

  cell_lines <- gsub(pattern=paste0(parent_path,"Scores/"), replacement="",x=files)
  cell_lines <- gsub(pattern="_.*", replacement="",x=cell_lines)

  data1 <- read.xlsx(files[1]) %>%
    dplyr::filter(!is.na(Score)) %>%
    dplyr::distinct_at(c("Gene","Score"), .keep_all = TRUE) %>%
    dplyr::select(Gene, Score, Percent_n, Cutoff, SD, Hits) %>%
    dplyr::mutate(Cell = cell_lines[1],
                  Control = dplyr::case_when(Gene == "negative_control" | grepl(pattern="Olfr",x=Gene) ~ "Control",
                                             TRUE ~ ""))

  data2 <- read.xlsx(files[2]) %>%
    dplyr::filter(!is.na(Score)) %>%
    dplyr::distinct_at(c("Gene","Score"), .keep_all = TRUE) %>%
    dplyr::select(Gene, Score, Percent_n, Cutoff, SD,Hits) %>%
    dplyr::mutate(Cell = cell_lines[2],
                  Control = dplyr::case_when(Gene == "negative_control" | grepl(pattern="Olfr",x=Gene) ~ "Control",
                                             TRUE ~ ""))

  data3 <- read.xlsx(files[3]) %>%
    dplyr::filter(!is.na(Score)) %>%
    dplyr::distinct_at(c("Gene","Score"), .keep_all = TRUE) %>%
    dplyr::select(Gene, Score, Percent_n, Cutoff, SD,Hits) %>%
    dplyr::mutate(Cell = cell_lines[3],
                  Control = dplyr::case_when(Gene == "negative_control" | grepl(pattern="Olfr",x=Gene) ~ "Control",
                                             TRUE ~ ""))

  data <- dplyr::bind_rows(data1, data2) %>%
    dplyr::bind_rows(data3) %>%
    dplyr::mutate(Size = Percent_n*abs(Score))

  palette <- c("#808080", "#66A61E", "#E7298A", "#E6AB02", "#D95F02" )
  names(palette) <- c("Not Significant", "Depletion 2SD", "Depletion 1SD", "Enrichment 1SD", "Enrichment 2SD")

  ggplot(data=data, mapping=aes(x=Score, y=Percent_n, color=Hits, shape=Cell, size=Percent_n)) +
    geom_point(alpha = 0.5, stroke = 0) +
    theme_classic() +
    coord_cartesian(ylim=c(0,100)) + 
    geom_hline(yintercept = 60, linetype="dashed") +
    scale_color_manual(values=palette) + 
    scale_size_continuous(range  = c(0.1, 3), 
                          limits = c(0, 100), 
                          breaks = c(0, 20, 40, 60, 80, 100))

  ggsave(paste0(library,"_Hits_1.tiff"),
         width = 8,
         height = 5,
         path = parent_path)

  ggplot(data=data, mapping=aes(x=Gene, y=Score, color=Hits, shape=Cell, size=Percent_n)) +
    geom_point(alpha = 0.5, stroke = 0) +
    theme_classic() +
    scale_color_manual(values=palette) + 
    scale_size_continuous(range  = c(0.1, 3), 
                          limits = c(0, 100), 
                          breaks = c(0, 20, 40, 60, 80, 100))

  ggsave(paste0(library,"_Hits_2.tiff"),
         width = 8,
         height = 5,
         path = parent_path)
  
  ggplot(data=data, mapping=aes(x=reorder(Score, Gene), y=Score, color=Hits, shape=Cell, size=Percent_n)) +
    geom_point(alpha = 0.5, stroke = 0) +
    theme_classic() +
    scale_color_manual(values=palette) +
    scale_size_continuous(range  = c(0.1, 3), 
                          limits = c(0, 100), 
                          breaks = c(0, 20, 40, 60, 80, 100))
  
  ggsave(paste0(library,"_Hits_3.tiff"),
         width = 8,
         height = 5,
         path = parent_path)
}

files <- list.files(paste0(parent_path, "Scores/"), full.names = TRUE)
files <- files[grepl(pattern=".xlsx", x= files)]
for (f in files){  
  
  filename <- gsub(pattern=paste0(parent_path,"Scores/"), replacement="",x=f)
  cell_line <- gsub(pattern="_.*", replacement="",x=filename)
  library <- gsub(pattern=paste0(cell_line, "_"), replacement="",x=filename)
  library <- gsub(pattern="_.*", replacement="",x=library)
  
  data <- read.xlsx(files[1]) %>%
    dplyr::filter(!is.na(Score)) %>%
    dplyr::distinct_at(c("Gene","Score"), .keep_all = TRUE) %>%
    dplyr::select(Gene, Score, Percent_n, Cutoff, SD, Hits) %>%
    dplyr::mutate(Cell = cell_lines[1],
                  Control = dplyr::case_when(Gene == "negative_control" | grepl(pattern="Olfr",x=Gene) ~ "Control",
                                             TRUE ~ "")) %>%
    dplyr::mutate(Size = Percent_n*abs(Score))
  
  palette <- c("#808080", "#66A61E", "#E7298A", "#E6AB02", "#D95F02" )
  names(palette) <- c("Not Significant", "Depletion 2SD", "Depletion 1SD", "Enrichment 1SD", "Enrichment 2SD")
  
  ggplot(data=data, mapping=aes(x=Score, y=Percent_n, color=Hits, shape=Cell, size=Percent_n)) +
    geom_point(alpha = 0.5, stroke = 0) +
    theme_classic() +
    coord_cartesian(ylim=c(0,100)) +  
    geom_hline(yintercept = 60, linetype="dashed") +
    scale_color_manual(values=palette) +
    scale_size_continuous(range  = c(0.1, 3), 
                          limits = c(0, 100), 
                          breaks = c(0, 20, 40, 60, 80, 100))
  
  ggsave(paste0(cell_line, "_", library,"_Hits_1.tiff"),
         width = 8,
         height = 5,
         path = parent_path)
  
  ggplot(data=data, mapping=aes(x=Gene, y=Score, color=Hits, shape=Cell, size=Percent_n)) +
    geom_point(alpha = 0.5, stroke = 0) +
    theme_classic() +
    scale_color_manual(values=palette) +
    scale_size_continuous(range  = c(0.1, 3), 
                          limits = c(0, 100), 
                          breaks = c(0, 20, 40, 60, 80, 100))
 
  ggsave(paste0(cell_line, "_", library,"_Hits_2.tiff"),
         width = 11,
         height = 5,
         path = parent_path)
         
  
  ggplot(data=data, mapping=aes(x=reorder(Score, Gene), y=Score, color=Hits, shape=Cell, size=Percent_n)) +
    geom_point(alpha = 0.5, stroke = 0) +
    theme_classic() +
    scale_color_manual(values=palette) + 
    scale_size_continuous(range  = c(0.1, 3), 
                          limits = c(0, 100), 
                          breaks = c(0, 20, 40, 60, 80, 100))
  
  ggsave(paste0(cell_line, "_", library,"_Hits_3.tiff"),
         #width = 11,
         #height = 5,
         path = parent_path)
  
}

######Plot venn diagram

libraries <- c("CDH12KO", "CDH12ACT", "ImmuneKO", "ImmuneACT")
pattern <- c("Depletion 2SD", "Depletion","Enrichment 2SD","Enrichment")
for (library in libraries){
  for (p in pattern){
  
  files <- list.files(paste0(parent_path, "Scores_New100/"), full.names = TRUE)
  files <- files[grepl(pattern=".xlsx", x= files)]
  files <- files[grepl(pattern=library, x= files)]
  
  cell_lines <- gsub(pattern=paste0(parent_path,"Scores_New100/"), replacement="",x=files)
  cell_lines <- gsub(pattern="_.*", replacement="",x=cell_lines)
  
  data1 <- read.xlsx(files[1]) %>%
    dplyr::filter(!is.na(Score)) %>%
    dplyr::distinct_at(c("Gene","Score"), .keep_all = TRUE) %>%
    dplyr::filter(grepl(pattern=p, x=Hits)) %>%
    dplyr::select(Gene) %>%
    unlist(use.names=FALSE)
  
  data1 <- data1[!grepl(pattern="negative_control|Olfr", x=data1)]
  
  data2 <- read.xlsx(files[2]) %>%
    dplyr::filter(!is.na(Score)) %>%
    dplyr::distinct_at(c("Gene","Score"), .keep_all = TRUE) %>%
    dplyr::filter(grepl(pattern=p, x=Hits)) %>%
    dplyr::select(Gene) %>%
    unlist(use.names=FALSE)
  
  data2 <- data2[!grepl(pattern="negative_control|Olfr", x=data2)]
  
  data3 <- read.xlsx(files[3]) %>%
    dplyr::filter(!is.na(Score)) %>%
    dplyr::distinct_at(c("Gene","Score"), .keep_all = TRUE) %>%
    dplyr::filter(grepl(pattern=p, x=Hits)) %>%
    dplyr::select(Gene) %>%
    unlist(use.names=FALSE)
  
  data3 <- data3[!grepl(pattern="negative_control|Olfr", x=data3)]
  
  l <-  max(length(data1), length(data2), length(data3))
  data1 <- c(data1, rep(x="",times=l-length(data1)))
  data2 <- c(data2, rep(x="",times=l-length(data2)))
  data3 <- c(data3, rep(x="",times=l-length(data3)))
  
  data <- dplyr::bind_cols(data1, data2,data3) %>%
    data.frame()
  
  colnames(data) <- cell_lines
  suffix <- paste0(library, "_", p)
  
  plot_venn(data, parent_path, suffix)
  }
}


############ Calculate missing guides

source("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

parent_path <- "C:/Users/kailasamms/Box/Saravana@cedars/05. Bioinformatics/CRISPR/CRISPR_Lena/"
count_path <- paste0(parent_path, "count_results/kms/")
DEG_path <- paste0(parent_path, "DEG_results/")

libraries <- c("CDH12KO", "CDH12ACT", "ImmuneKO", "ImmuneACT")

files <- list.files(paste0(count_path, libraries[1]), full.names = TRUE)
samples <- gsub("^.*/|\\.csv", "", files)
summary <- data.frame(Sample = samples)

for (library in libraries){
  
  files <- list.files(paste0(count_path, library), full.names = TRUE)
  samples <- gsub("^.*/|\\.csv", "", files)
  zero <- c()
  total <- c()
  
  print("Here")
  for (f in files){
    data <- read.table(f, header=TRUE, sep=",") 
    zero <- c(zero,sum(data$COUNT == 0))
    total <- c(total, nrow(data))
  }
  
  df <- data.frame(Sample = samples,
                   Missing = zero,
                   Total = total) %>%
    dplyr::mutate(Percent = round(100*Missing/Total,2)) %>%
    #dplyr::select(Sample, Percent) %>%
    dplyr::rename(!!rlang::sym(library) := "Percent")
  
  summary <- summary %>%
    dplyr::left_join(df, by=c("Sample"="Sample"))
}

#Save output in excel
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Missing")
openxlsx::writeData(wb, sheet = "Missing", x = summary)
openxlsx::saveWorkbook(wb, file = paste0(parent_path, "Missing_Guides.xlsx"), overwrite=TRUE)
