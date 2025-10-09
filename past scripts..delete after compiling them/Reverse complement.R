data <- read.xlsx("C:/Users/KailasammS/Desktop/CRISPR_Jinfen(from Addgene).xlsx", colNames = FALSE)
data <- data[,1:3]
colnames(data) <- c("ID", "Seq","Gene")
data <- data %>% dplyr::mutate(Rev_Seq = "")

for (index in 1:nrow(data)){
  
  string = strsplit(data[index,2],"")[[1]]
  rc_string = string
  for (i in seq(from=length(string), to=1, by=-1)){
    rc_string[(length(string)-i+1)] = dplyr::case_when(string[i] == "A" ~ "T",
                                                       string[i] == "T" ~ "A",
                                                       string[i] == "G" ~ "C",
                                                       string[i] == "C" ~ "G",
                                                       TRUE ~ string[i])
  }
  rc_string = paste(rc_string, collapse = "")
  data[index,4] <- rc_string
}

data <- data %>% dplyr::mutate(Mageck_format = paste(ID, Seq, Gene, sep = ","),
                               Mageck_format_rev = paste(ID, Rev_Seq, Gene, sep = ","))

# Save the results
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "Guides")
openxlsx::writeData(wb, sheet = "Guides", x = data, rowNames = FALSE)
openxlsx::saveWorkbook(wb, 
                       file = paste0("C:/Users/KailasammS/Desktop/CRISPR_Jinfen(from Addgene).xlsx"), 
                       overwrite = TRUE)

