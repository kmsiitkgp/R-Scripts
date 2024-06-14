# Mouse annotations from ensembl
mouse_full <- read.xlsx("C:/Users/KailasammS/Desktop/Mouse.xlsx", sheet="MOUSE")
# Human annotations from ensembl
human_full <- read.xlsx("C:/Users/KailasammS/Desktop/Mouse.xlsx", sheet="HUMAN")

#Get human_mouse homologs from https://www.informatics.jax.org/downloads/reports/index.html
# Human and Mouse Homology Classes with Sequence information (tab-delimited)
# Separate them into human and mouse on different sheets
mouse_homology <- read.xlsx("C:/Users/KailasammS/Desktop/homology.xlsx", sheet="MOUSE")
human_homology <- read.xlsx("C:/Users/KailasammS/Desktop/homology.xlsx", sheet="HUMAN")

# Filter genes without gene symbols, remove duplicates
mouse_full <- mouse_full %>% 
  dplyr::filter(nchar(MOUSE_SYMBOL) > 0) %>%
  dplyr::distinct(pick(MOUSE_SYMBOL), .keep_all = TRUE)

# Filter genes without gene symbols, remove duplicates
human_full <- human_full %>% 
  dplyr::filter(nchar(HUMAN_SYMBOL) > 0) %>%
  dplyr::distinct(pick(HUMAN_SYMBOL), .keep_all = TRUE)

# Filter duplicated homology pairs
mouse_homology <- mouse_homology %>%
  dplyr::distinct(pick(M_Symbol, DB.Class.Key), .keep_all = TRUE)

# Filter duplicated homology pairs
human_homology <- human_homology %>%
  dplyr::distinct(pick(H_Symbol, DB.Class.Key), .keep_all = TRUE)

# # Find duplicates
# nrow(mouse_full %>% dplyr::count(MOUSE_SYMBOL) %>% dplyr::filter(n>1))
# nrow(mouse_homology %>% dplyr::count(M_Symbol) %>% dplyr::filter(n>1))
# length(unique(intersect(mouse_full$MOUSE_SYMBOL, mouse_homology$M_Symbol)))

mouse_final <- dplyr::inner_join(mouse_full, mouse_homology, by=c("MOUSE_SYMBOL"="M_Symbol")) %>%
  dplyr::select(DB.Class.Key, MOUSE_ID, MOUSE_SYMBOL, M_Genetic.Location, everything()) %>%
  dplyr::rename(MOUSE_GENETIC_LCATION =M_Genetic.Location)

# # Find duplicates
# nrow(human_full %>% dplyr::count(HUMAN_SYMBOL) %>% dplyr::filter(n>1))
# nrow(human_homology %>% dplyr::count(H_Symbol) %>% dplyr::filter(n>1))
# length(unique(intersect(human_full$HUMAN_SYMBOL, human_homology$H_Symbol)))

human_final <- dplyr::inner_join(human_full, human_homology, by=c("HUMAN_SYMBOL"="H_Symbol")) %>%
  dplyr::select(DB.Class.Key, HUMAN_ID, HUMAN_SYMBOL, H_Genetic.Location, everything()) %>%
  dplyr::rename(HUMAN_GENETIC_LCATION =H_Genetic.Location)


homologs1 <- dplyr::inner_join(mouse_homology, human_homology, by=c("DB.Class.Key"="DB.Class.Key")) %>%
  dplyr::select(DB.Class.Key, MOUSE_ID, MOUSE_SYMBOL, HUMAN_ID, HUMAN_SYMBOL, everything(), -c(MOUSE_DESCRIPTION, HUMAN_DESCRIPTION))


