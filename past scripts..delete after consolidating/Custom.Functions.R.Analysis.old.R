# NOTE: 'proj' must be defined by the parent script
# Metadata Excel file should be named <proj>_Metadata.xlsx
# Must contain a column "Unique_ID" matching Seurat object metadata

#******************************************************************************#
#                             MARKER IDENTIFICATION                            #
#******************************************************************************#

find_markers_across_datasets <- function(datasets){
  
  # Rather than trying to identify clusters that are similar across datasets 
  # and then trying to find marker genes, we try to find which genes are 
  # co-expressed frequently across datasets and then assign cell types based on 
  # such co-expressed gene sets. Since our approach is focused on identifying 
  # markers that are frequently co-expressed across all clusters from all 
  # datasets, it is robust to:
  # (i) sequencing depth (low UMIs dataset vs high UMIs dataset)
  # (ii) cell composition (immune enriched vs whole tumor dataset)
  # (iii) experiment type (single cell vs single nuclei dataset)
  # (iv) purity (high ambient RNA vs low ambient RNA dataset)
  # (v) clustering resolution
  
  datasets <- c("scRNASeq_BBN_C57BL6", "scRNASeq_BBN_Rag", "scRNASeq_GSE217093", 
                "scRNASeq_Jinfen", "scRNASeq_Jyoti", "scRNASeq_GSE164557",
                "scRNASeq_Chen", "scRNASeq_GSE222315", "scRNASeq_HRA003620")
  
  # Get human to mouse ortholog mapping
  ortho <- get_orthologs()
  
  # Create empty dataframe to store markers from all datasets
  markers <- data.frame(cluster = c(""))
  
  # Read markers identified at resolution Harmony 0.8 from each dataset
  for (proj in datasets){
    
    seurat_results <- paste0("/hpc/home/kailasamms/scratch/", proj, "/results_seurat/")
    
    # Genes present repeatedly due to DESCRIPTION column. Remove them.
    df <- read.xlsx(paste0(seurat_results, proj, ".Markers.All.cluster.0.8.harmony.Full.xlsx")) %>%
      dplyr::mutate(Proj = gsub(pattern="scRNASeq_", replacement="", x=proj),
                    Proj.cluster = paste0(Proj, ".", cluster),
                    Avg_Expr = avg_log2FC*pct.1) %>%
      dplyr::distinct_at(c("gene", "cluster", "Proj", "p_val_adj", "avg_log2FC", "pct.1", "pct.2"), .keep_all = TRUE) %>%
      dplyr::select(Proj, cluster, Proj.cluster, gene, p_val_adj, avg_log2FC, pct.1, pct.2, ratio, Avg_Expr)
    
    # Merge markers from all datasets
    markers <- dplyr::bind_rows(markers, df) %>%
      dplyr::filter(!is.na(Proj))
    cat(nrow(markers), "\n")
  }
  
  # Add Human orthologs after removing poor markers
  markers <- markers %>%
    dplyr::filter(p_val_adj <= 0.05) %>% #, pct.1 >= 0.4, ratio >= 2) %>%
    dplyr::left_join(ortho, by=c("gene"="Mouse")) %>%
    dplyr::mutate(avg_log2FC = round(avg_log2FC, 2),
                  ratio = round(ratio, 2),
                  Avg_Expr = round( Avg_Expr, 2),
                  p_val_adj = round(p_val_adj, 2),
                  SYMBOL = dplyr::case_when(is.na(Human) ~ gene,
                                            TRUE ~ Human))
  
  # Get top 100 markers based on avg_log2FC, ratio, Avg_Expr, pct.1 for each cluster
  markers_log2FC <- markers %>%
    dplyr::group_by(Proj.cluster) %>%
    dplyr::slice_max(n=100, order_by = avg_log2FC)
  
  markers_ratio <- markers %>%
    dplyr::group_by(Proj.cluster) %>%
    dplyr::slice_max(n=100, order_by = ratio) %>%
    dplyr::filter(ratio > 1)
  
  markers_pct1 <- markers %>%
    dplyr::group_by(Proj.cluster) %>%
    dplyr::slice_max(n=100, order_by = pct.1) %>%
    # if you exclude this filter, you will get specific but sparsely expressed genes
    dplyr::filter(pct.1 >= 0.4)  
  
  markers_expr <- markers %>%
    dplyr::group_by(Proj.cluster) %>%
    dplyr::slice_max(n=100, order_by = Avg_Expr)
  
  # Merge top 100 markers for each cluster and remove duplicates
  markers_top <- dplyr::bind_rows(markers_log2FC, markers_ratio, 
                                  markers_pct1, markers_expr) %>%
    dplyr::distinct_all(.keep_all = TRUE)
  
  # Get all possible combinations of markers
  marker.pair.df <- data.frame()
  count <- 0
  for (i in unique(markers_top$Proj.cluster)){
    
    # Get markers from each proj and cluster
    markers.subset <- markers_top %>%
      dplyr::filter(Proj.cluster == i)
    
    count <- count+1
    cat(count, ":", nrow(markers.subset), "\t")
    
    if (nrow(markers.subset) >= 2){
      df <- utils::combn(x=markers.subset$SYMBOL, m=2) %>%
        #df <- mixtools::perm(n=length(markers.subset$SYMBOL), r=2, v=markers.subset$SYMBOL)
        t() %>%
        data.frame()
      
      marker.pair.df <- dplyr::bind_rows(marker.pair.df, df)
    }
  }  
  
  # Count all possible combinations of markers from all datasets
  marker.pair.df <- marker.pair.df %>%
    dplyr::rename(PairA = identity(1), PairB = identity(2)) %>%
    dplyr::add_count(PairA, PairB) %>%
    dplyr::distinct_at(c("PairA", "PairB"), .keep_all = TRUE) %>%
    dplyr::rename(n_clusters = n) %>%
    dplyr::filter(n_clusters >=3)   # remove combinations not observed in even 3 clusters
  
  # Calculate overlapping number of genes between PairA and PairB
  marker.pair.df$n_common <- 0  # number of genes commonly coexpressed between A & B
  marker.pair.df$n_PairA <- 0   # number of genes coexpressed with A
  marker.pair.df$n_PairB <- 0   # number of genes coexpressed with B
  marker.pair.df$n_ratio <- 0   # n_common/(n_PairA+n_PairB-n_Common)
  for (i in 1:nrow(marker.pair.df)){
    
    # Find all genes coexpressed with PairA
    coexp_A1 <- marker.pair.df %>% 
      dplyr::filter(PairA == marker.pair.df$PairA[i]) %>% 
      dplyr::select(PairB) %>% unlist(use.names=FALSE)
    coexp_A2 <- marker.pair.df %>% 
      dplyr::filter(PairB ==  marker.pair.df$PairA[i]) %>% 
      dplyr::select(PairA) %>% unlist(use.names=FALSE)
    coexp_A <- unique(coexp_A1, coexp_A2)
    
    # Find all genes coexpressed with PairB
    coexp_B1 <- marker.pair.df %>% 
      dplyr::filter(PairA ==  marker.pair.df$PairB[i]) %>% 
      dplyr::select(PairB) %>% unlist(use.names=FALSE)
    coexp_B2 <- marker.pair.df %>% 
      dplyr::filter(PairB == marker.pair.df$PairB[i]) %>% 
      dplyr::select(PairA) %>% unlist(use.names=FALSE)
    coexp_B <- unique(coexp_B1, coexp_B2)
    
    # Calculate stats
    marker.pair.df$n_common[i] <- length(intersect(coexp_A, coexp_B))
    marker.pair.df$n_PairA[i] <- length(coexp_A)
    marker.pair.df$n_PairB[i] <- length(coexp_B)
    marker.pair.df$n_ratio[i] <- marker.pair.df$n_common[i]/(marker.pair.df$n_PairA[i]+marker.pair.df$n_PairB[i]-marker.pair.df$n_common[i])
    
    cat(i, "\t")
  }
  
  # Remove all pairs that have poor overlap (n_ratio < 0.5) 
  top.marker.pair.df <- marker.pair.df %>% 
    dplyr::filter(n_ratio >= 0.5)
  
  final.markers <- list(T.NK.cell        = c("CD3D"),
                        B.Plasma.cell    = c("CD79A"),
                        Erythrocyte      = c("HBB"),
                        Mast.cell        = c("KIT"),   # tissue resident granule producing cell
                        #Granulocyte      = c(),        # blood resident granule producing cell (Basophil, Eosinophil, Neutrophil)                     
                        Monocyte         = c("GOS2"),  # blood resident phagocyte
                        Macrophage       = c("C1QA"),  # tissue resident phagocyte
                        #Dendritic.cell   = c(),
                        Endothelial.cell = c("VWF"),
                        Myocyte          = c("MYH11"),
                        Neurons          = c("KCNA1"))
  
  # Fill the marker list
  for (i in 1:length(final.markers)){
    
    final.markers[[i]] <- c(final.markers[[i]], marker.pair.df %>% 
                              dplyr::filter(PairA == final.markers[[i]]) %>% 
                              dplyr::select(PairB) %>%
                              unlist(use.names=FALSE))
  }
  
  # Convert to dataframe
  max_l <- max(lengths(final.markers)) 
  final.markers.df <- lapply(X=final.markers, FUN=function(x){c(x, base::rep(x="", times=max_l-length(x)))})
  
  # Save the clustered similarity matrix
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "Markers")
  openxlsx::writeData(wb, sheet = "Markers", x = final.markers.df, rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName = "Datasets")
  openxlsx::writeData(wb, sheet = "Datasets", x = data.frame(datasets), rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName = "All")
  openxlsx::writeData(wb, sheet = "All", x = marker.pair.df, rowNames = FALSE)
  openxlsx::addWorksheet(wb, sheetName = "Top")
  openxlsx::writeData(wb, sheet = "Top", x = top.marker.pair.df, rowNames = FALSE)
  openxlsx::saveWorkbook(wb = wb, file = "Markers.Compiled.xlsx", overwrite = TRUE)
}

plot_markers_across_datasets <- function(datasets, markers){
  
  # Read integrated seurat object for each dataset
  for (d in datasets){
    seurat_results <- paste0("/hpc/home/kailasamms/scratch/", d, "/results_seurat/")
    integrated.seurat <- readRDS(paste0(seurat_results, "integrated.seurat.rds"))
    Idents(integrated.seurat) <- "cluster.0.8.harmony"
    assign(d, integrated.seurat)
  }
  
  #*****************Plot UMAP at Harmony 0.8 for each dataset******************#
  
  purrr::map(.x = datasets,
             .f = function(x){
               obj <- get(x)
               Seurat::DimPlot(object=obj,
                               reduction="umap.harmony",
                               cols=custom_palette,
                               pt.size=0.2,
                               label.size=1,
                               order = TRUE,  # plot doublets on above rest of cells
                               label = TRUE,
                               raster = FALSE,
                               combine = TRUE) +
                 NoLegend() +
                 ggplot2::labs(title = x,  x="UMAP_1", y="UMAP_2") +
                 custom_theme}) %>%
    cowplot::plot_grid(plotlist=.,
                       align="hv",
                       axis="tblr",
                       #nrow=,
                       ncol=3,
                       rel_widths=1,
                       rel_heights=1,
                       greedy=TRUE,
                       byrow=TRUE)
  
  # Save the plot
  ggplot2::ggsave(filename = "Marker_UMAPs.tiff",
                  plot = last_plot(),
                  device = "jpeg",
                  #path = ,
                  scale = 1,
                  width = 4*3,
                  height = 4*ceiling(length(datasets)/3),
                  units = c("in"),
                  dpi = 300,
                  limitsize = TRUE,
                  bg = "white")
  
  #*****************Plot UMAP of identified markers for each dataset******************#
  
  # Plot these putative markers on UMAP
  for (f in markers){
    
    purrr::map(.x = datasets,
               .f = function(x){
                 obj <- get(x)
                 Seurat::FeaturePlot(object = obj,
                                     features = intersect(rownames(obj@assays$RNA$counts), c(f, stringr::str_to_title(f))),
                                     reduction = "umap.harmony",
                                     cols = c("grey", viridis(n = 10, option = "C", direction = -1)),
                                     pt.size = 0.2,
                                     label.size = 1,
                                     min.cutoff='q10',
                                     order = TRUE,  # plot doublets on above rest of cells
                                     label = TRUE,
                                     raster = FALSE,
                                     combine = TRUE) +
                   NoLegend() +
                   # scale_colour_gradientn(colours=rev(brewer.pal(n=11, name="RdBu"))[5:11])
                   ggplot2::labs(title = x, x="UMAP_1", y="UMAP_2") +
                   custom_theme}) %>% cowplot::plot_grid(plotlist=.,
                                                     align="hv",
                                                     axis="tblr",
                                                     #nrow=,
                                                     ncol=3,
                                                     rel_widths=1,
                                                     rel_heights=1,
                                                     greedy=TRUE,
                                                     byrow=TRUE)
    
    ggplot2::ggsave(filename = paste0(f, ".tiff"),
                    plot = last_plot(),
                    device = "jpeg",
                    #path = ,
                    scale = 1,
                    width = 4*3,
                    height = 4*ceiling(length(datasets)/3),
                    units = c("in"),
                    dpi = 300,
                    limitsize = TRUE,
                    bg = "white")
    
  }
}

#******************************************************************************#
#                               GENE ANNOTATIONS                               #
#******************************************************************************#

### Get human-mouse orthologs
# Output is a dataframe with columns DB.Class.Key, Human, Mouse
# NOTE: Mouse genes (H2-Q10,H2-Q8,H2-Q7,H2-Q6,H2-Q4,H2-Q2,H2-Q1,H2-T23,H2-K1,
# H2-D1) are othologs of the same human gene HLA-A. Similarly, human genes 
# ZNG1A, ZNG1B, ZNG1C, ZNG1E, ZNG1F) are orthologs of the same mouse gene (Zng1).
# So, we make them unique as well syntactically valid (hyphens repalced with .)
# to avoid errors in data analysis.
get_orthologs <- function(){
  
  # This website has a list of human mouse orthologs
  df <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  # Get human genes
  df_h <- df %>% 
    dplyr::select(DB.Class.Key, Common.Organism.Name, Symbol) %>%
    dplyr::filter(Common.Organism.Name == "human") %>%
    dplyr::select(DB.Class.Key, Symbol) %>%
    dplyr::distinct_at(c("DB.Class.Key","Symbol"),.keep_all = TRUE) %>%
    dplyr::rename(Human = Symbol)
  
  # Get mouse genes 
  df_m <- df %>% 
    dplyr::select(DB.Class.Key, Common.Organism.Name, Symbol) %>%
    dplyr::filter(Common.Organism.Name == "mouse, laboratory") %>%
    dplyr::select(DB.Class.Key, Symbol) %>%
    dplyr::distinct_at(c("DB.Class.Key","Symbol"),.keep_all = TRUE) %>%
    dplyr::rename(Mouse = Symbol)
  
  # Get human-mouse orthologs & remove genes that dont have orthologs
  df_h_m <- dplyr::full_join(df_h, df_m, by=c("DB.Class.Key"="DB.Class.Key")) %>%
    base::replace(is.na(.), "None") %>%
    dplyr::filter(Human != "None", Mouse != "None")
  #df_h_m[is.na(df_h_m)] <- "None"
  
  # Similar orthologs (mouse and human gene names are identical)
  conf_h_m <- df_h_m %>% dplyr::filter(Human == base::toupper(Mouse))
  
  # Dissimilar orthologs (mouse and human gene names are NOT identical)
  fix_h_m <- df_h_m %>% 
    dplyr::filter(!(Human %in% conf_h_m$Human)) %>%
    dplyr::filter(!(Mouse %in% conf_h_m$Mouse))
  
  # Merge
  df <- dplyr::bind_rows(conf_h_m, fix_h_m) %>%
    dplyr::mutate(Human = make.names(Human, unique=TRUE),
                  Mouse = make.names(Mouse, unique=TRUE))
  
  # Save the excel file  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb=wb, sheetName="Homologs")
  openxlsx::writeData(wb=wb, sheet="Homologs", x=df)
  openxlsx::saveWorkbook(wb=wb, file="Human.Mouse.Homologs.xlsx", overwrite=TRUE)
  
  return (df)
  
  # # Get mouse and human annotations
  # df <- get_annotations()
  # h <- df[[1]]
  # m <- df[[2]]
  # 
  # # Identify identical genes
  # common <- base::intersect(h$SYMBOL, base::toupper(m$SYMBOL))
  # common_h <- h %>% dplyr::filter(h$SYMBOL %in% common)
  # common_m <- m %>% dplyr::filter(base::toupper(m$SYMBOL) %in% common)
  # 
  # # Identify non-common genes
  # uniq_h <- h %>% dplyr::filter(!(h$SYMBOL %in% common))
  # uniq_m <- m %>% dplyr::filter(!(base::toupper(m$SYMBOL) %in% common))
  # 
  # # Identify similar genes (TIME CONSUMING. SO, save excel and use it in future)
  # # Eg: CD3 and CD3D etc
  # similar_h <- c()
  # similar_m <- c()
  # for (p in unique(uniq_h$SYMBOL)){
  #   
  #   if (sum(grepl(pattern=paste0(p, "."), x=toupper(unique(uniq_m$SYMBOL)))) > 0){
  #     similar_h <- c(similar_h, p)
  #     similar_m <- c(similar_m, unique(uniq_m$SYMBOL)[grepl(pattern=paste0(p, "."), x=toupper(unique(uniq_m$SYMBOL)))][1])
  #                  
  #     #cat(p, ":", unique(uniq_m$SYMBOL)[grepl(pattern=paste0(p, "."), x=toupper(unique(uniq_m$SYMBOL)))][1], "\n")
  #   }
  # }
  # 
  # # Generate final dataframe of human-mouse homolog
  # common_df <- data.frame(HUMAN = unique(sort(common_h$SYMBOL)), 
  #                  MOUSE = unique(sort(common_m$SYMBOL)))
  # similar_df <- data.frame(HUMAN = similar_h,
  #                          MOUSE = similar_m)
  # df1 <- dplyr::bind_rows(common_df, similar_df)
}

#******************************************************************************#
#                      RNA SEQ ANALYSIS RELATED FUNCTIONS                      #       
#******************************************************************************#

### Read individual count files.txt and merge all counts into "Read_data.xlsx"
# Input is path to counts file, project and method (STAR, HTSeq)
compile_raw_counts <- function(count_dir, proj, method){
  
  data_path <- gsub(pattern = paste0("counts/", method, "_raw_counts/"), 
                    replacement = "", x = count_dir)
  
  # Create a list of all txt files within folder that will be analyzed
  files <- list.files(path=count_dir)
  
  # Create an empty dataframe with 0's
  read_data <- data.frame(0)
  
  # Create the reads table 
  for (i in 1:length(files)){
    
    # Read the txt file
    temp_file <- read.table(file=paste0(count_dir, files[i]), header=FALSE, sep="\t")  
    
    if (method == "HTSEQ"){
      # Remove last 5 rows in HTSeq count output: __no_feature,  __ambiguous, 
      # __too_low_aQual, __not_aligned, __alignment_not_unique
      temp_file <- temp_file[1:(nrow(temp_file)-5),]
    } else if (method == "STAR"){
      # Remove top 4 rows in STAR count output:  
      # N_unmapped, N_multimapping, N_noFeature, N_ambiguous
      temp_file <- temp_file[5:nrow(temp_file),]
    }
    
    # The 1st Column will have Ensembl ids. 
    # For HTSeq, gene counts may be in 2nd or 3rd column. 
    # For STAR, gene counts are in 2nd (unstranded), 3rd (+), 4th (-) column. 
    # Append appropriate column.
    if (method == "HTSEQ" & sum(temp_file[2], na.rm=TRUE) == 0 & sum(temp_file[3], na.rm=TRUE) > 0){
      read_data <- bind_cols(read_data, temp_file[,3])
    } else if (method == "HTSEQ" &  sum(temp_file[2], na.rm=TRUE) > 0 & sum(temp_file[3], na.rm=TRUE) ==0){
      read_data <- bind_cols(read_data, temp_file[,2])                  
    } else if (method == "STAR" & abs((sum(temp_file[2])/sum(temp_file[3])) - (sum(temp_file[2])/sum(temp_file[4]))) < 2){
      print("Unstranded")
      read_data <- bind_cols(read_data, temp_file[,2])                  
    } else if (method == "STAR" & sum(temp_file[3]) > 3*sum(temp_file[4])){
      print("Pos stranded")
      read_data <- bind_cols(read_data, temp_file[,3])                  
    } else if (method == "STAR" & sum(temp_file[4]) > 3*sum(temp_file[3])){
      print("Neg stranded")
      read_data <- bind_cols(read_data, temp_file[,4])                  
    } else{
      print("Error: Gene counts NOT PRESENT in either column 2 or 3 of count file")
    }
    
    # Rename the column names to sample names
    colnames(read_data)[i+1] <- gsub(pattern="\\..*$|ReadsPerGene.out.tab", replacement="", x=files[i])
  }
  
  # Check if all count files have same order of genes in the rows so that files can be merged together
  temp_file <- read.table(file=paste0(count_dir, files[1]), header=FALSE, sep="\t")
  gene_list <- temp_file[, 1]
  for (i in 1:length(files)){
    
    temp_file <- read.table(file=paste0(count_dir, files[i]), header=FALSE, sep="\t")
    genes <- temp_file[, 1]
    
    if (!identical(gene_list, genes)){ 
      print("Gene order is different between the count files")
    }
  }
  
  if (method == "STAR"){
    gene_list <- gene_list[5:length(gene_list)]
  } else if(method == "HTSEQ"){
    gene_list <-gene_list[1:(length(gene_list)-5)]
  }
  
  # Add gene names to 1st column
  read_data[, 1] <- gene_list
  colnames(read_data)[1] <- "SYMBOL"
  colnames(read_data) <- gsub(pattern="ReadsPerGene", replacement="", x=colnames(read_data))
  
  # Remove genes with 0 counts in all samples
  read_data <- read_data[rowSums(read_data[,-1]) != 0,]
  
  # Remove samples with 0 counts in all genes
  read_data <- read_data[,colSums(read_data[,-1]) != 0]
  
  # Save the results as xlsx file
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Raw_counts")
  openxlsx::writeData(wb, sheet="Raw_counts", x=read_data)
  openxlsx::saveWorkbook(wb, file=paste0(data_path, proj, ".raw.counts.xlsx"),
                         overwrite=TRUE)
  
  return(read_data)
}

# Example: If samples were isolated on different days & prepared using different
# kits, "Batch" column must have values like "1_Ribo", "1_Poly-A", "2_Ribo", etc
# NOTE: Make sure there are no white spaces in the Target and Reference columns
# in excel file. R will change any white space (" ") to dot ("."). So, replace 
# white space (" ") with underscore ("_") before importing into R.

# DESeq2 automatically removes genes that have 0 counts but doesnt remove 
# samples that have 0 counts for all genes (unlikely scenario). So, remove
# such samples, else, the geometric mean will be 0 for all genes and DESeq2 
# will halt execution.

# Built-in PCA plot in DESeq2
# Use getMethod("plotPCA", "DESeqTransform") to understand how DESeq2 makes 
# PCA plot, it uses top 500 genes with highest variance and uses scale=FALSE
# pcaData <- plotPCA(object = vsd, intgroup = "Sample_ID", returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))

# Calculate Principal components
# (i) PCA is performed across columns. So, have variables i.e. genes on columns.
# (ii) "Error cannot rescale a constant/zero column to unit variance" is due  
# to genes that have 0 expression across all samples. Remove them.
# (iii) All NA must be replaced with 0
# (iv) prcomp() is better than princomp()

# Normalized counts are influenced by sizeFactors.
# sizeFactors are affected by number of samples (all samples vs subset of samples)
# sizeFactors are NOT affected by design formula.
# sizeFactors MUST be estimated first before normalization.
# Normalized counts from dds object are NOT batch corrected. We do this below.
# https://www.biostars.org/p/490181/

# design doesnt affect size factors. Hence, normalized counts are not affected by design
# but vst counts are affected by design blind=TRUE vs blind=FALSE

# AnnotationHub has SYMBOL-ENSEMBL_ID info ONLY.
# AnnotationDbi has SYMBOL-ENSEMBL_ID as well as SYMBOL-ENTREZ_ID info.
# hubCache(AnnotationHub()) to find location where cache is stored and delete
# it and start fresh if you get errors like "Error: failed to load resource"
# NOTE: mapIds can ONLY retrieve one of "EMSEMBL/SYMBOL/GENETYPE" at a time
# mapping <- AnnotationDbi::mapIds(x=org.Hs.eg.db, 
#                                  keys=keys(org.Hs.eg.db),
#                                  keytype="ENTREZID", 
#                                  column="SYMBOL") %>%
#   as.data.frame(do.call(cbind, list(.))) %>%
#   tibble::rownames_to_column("ENTREZID") %>%
#   dplyr::rename(ENTREZID=identity(1), SYMBOL=identity(2))

# Two types of analysis: 
# (i) Gene Set Enrichment Analysis (GSEA)
# (ii) Over Representation Analysis (ORA)
# GSEA uses all DEGs ordered by fold change or other parameter
# ORA uses ONLY significant DEGs and ignores fold change etc

# NOTE: Genes MUST be ranked i.e. sorted in descending fold change. You can 
# also rank based on log2FC & p value like: sign(df$log2fc)*(-log10(df$pval)))

# NOTE: Genes MUST be stored in list format, not as a dataframe.

# NOTE: No NA MUST be present in SYMBOL column. Else, fgsea::collapsePathways()
# will give "Error in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  : 
# Not all stats values are finite numbers"
# You can figure out using table(is.na(names(DEGs_list))), 
# is.infinite(DEGs_df$log2FoldChange) or sapply(DEGs_df, class) to make sure
# log2FoldChange and padj are numeric

# NOTE: If your excel file has "inf" in padj or log2FoldChange columns, the
# column will be read into R as character column instead of numeric. So, remove
# text values from log2FoldChange and padj columns

# Define score type in fgseaMultilevel() based on fold change values.
# NOTE: Use "pos", if you are ONLY interested in activated pathways.
# NOTE: Use "neg", if you are ONLY interested in inhibited pathways. 
# NOTE: Else, use "std" for both activated & inhibited pathways. 

# NOTE: If you run multiple gene sets like C5 and C2 together, padj will not 
# be significant as there will be too many multiple comparisons. So, run
# each gene set separately and merge results

# If you ordered your gene list based on fold change, then +ve NES indicates
# that the genes in this gene set are mostly at the top of your gene list
# (hence, most of them are upregulated) and -ve NES indicates that the genes
# in this gene set are mostly at the bottom of your gene list (hence, most 
# of them are downregulated)

# NOTE: Output of fgsea is a data.table & data.frame. 
# "leadingEdge" column is a list of genes. 
# So, DO NOT FORCE the output of fgsea to a dataframe as this will lead to 
# data loss from "leadingEdge" column & affect plotting using fgsea::plotEnrichment()

# NOTE: DO NOT USE labels for defining colors due to reasons below. 
# RECOMMEND using a named vector.
# NOTE: If using labels, sort labels in alphabetical order and then assign 
# color because R by default will arrange the labels in alphabetical order 
# first and then match them to colors indicated in values vector and then 
# color the plot. The coloring in the legend is however dependent on the 
# order of labels vector and values vector. To understand, create a plot first 
# using the sort and then without the sort(). 

### combatseq Batch correction:
# NOTE: This batch correction of known factors is done on raw counts
# The batch corrected raw reads are used in DESeq2
batch_correct_combat <- function(meta_data, read_data, formula_string){ 
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~1)
  dds <- DESeq2::estimateSizeFactors(dds) 
  dat  <- DESeq2::counts(dds, normalized = TRUE)

  # Full model matrix with the variable of interest
  mod  <- stats::model.matrix(as.formula(formula_string), colData(dds))
  
  if (length(unique(as.numeric(meta_data$Batch))) > 1){
    # Instead of using group & full_mod=TRUE, use covar_mod
    read_data_combat <- sva::ComBat_seq(counts=as.matrix(read_data), 
                                        batch=as.numeric(meta_data$Batch), 
                                        #group=as.numeric(as.factor(meta_data$Condition)),
                                        #full_mod = TRUE,
                                        group = NULL,
                                        covar_mod = mod)
  } else{
    read_data_combat <- read_data
  }
  return(read_data_combat) 
}

### svaseq Batch correction:
# NOTE: This batch correction of unknown factors is done on normalized counts
# NOTE: svaseq() can find n number of surrogate variables. If we model for all 
# of them there could be over correction. Hence, we limit batch correction to
# only the top 3 surrogate variables.
# Here, we just create a new object sva_dds with sva design variables
batch_correct_sva <- function(meta_data, read_data, formula_string){
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~1)
  dds <- DESeq2::estimateSizeFactors(dds) 
  dat  <- DESeq2::counts(dds, normalized = TRUE)
  
  # Remove lowly expressed genes before finding surrogate variables
  idx  <- rowMeans(dat) > 1
  dat  <- dat[idx, ]
  
  # Full model matrix with the variable of interest #Reduced/null model matrix with only an intercept term
  mod  <- stats::model.matrix(as.formula(formula_string), colData(dds))
  # Reduced/null model matrix with only an intercept term
  mod0 <- stats::model.matrix(~ 1, colData(dds))
  # Estimate all surrogate variables
  svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0)
  
  # If there are more than 3 SV, estimate upto 3 SVs
  # svseq$sv values will change based on n.sv
  if (svseq$n.sv > 3){
    svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0, n.sv=3)
  }
  
  # Add these SVs as columns to the DESeqDataSet and then add them to the design
  # ddssva$SV1 <- svseq$sv[,1]
  # ddssva$SV2 <- svseq$sv[,2]
  # design(ddssva) <- ~ SV1 + SV2 + id
  
  ddssva <- dds
  for (i in 1:ncol(svseq$sv)){
    var <- paste0("SV",i)
    ddssva[[var]] <- svseq$sv[,i]
  }
  
  design(ddssva) <- as.formula(paste0("~", 
                                      paste0("SV", seq(1:ncol(svseq$sv)), collapse = "+"), 
                                      "+", 
                                      gsub(pattern="~", replacement="",x=formula_string)))
  
  return(ddssva)
}

norm_counts_combat <- function(meta_data, read_data_combat, output_path){
  
  # design doesnt affect size factors. Hence, normalized counts are not affected by design
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data_combat,
                                        colData=meta_data, 
                                        design=~ 1)
  # Estimate sizefactors
  dds <- DESeq2::estimateSizeFactors(object = dds)
  
  # EXtract normalized counts from dds object
  normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
  
  normalized_counts <- normalized_counts %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID")
  
  # Add gene names
  normalized_counts <- add_annotation(normalized_counts, annotations)
  
  # Save the results
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Norm_counts")
  openxlsx::writeData(wb, sheet="Norm_counts", x=normalized_counts, rowNames=FALSE)
  openxlsx::saveWorkbook(wb, file=paste0(output_path, "Normalized.counts.combat.xlsx"), 
                         overwrite=TRUE)
  
  return(normalized_counts)
}

# svaseq corrected normalized counts:
# NOTE: ddssva object from svaseq_batch has the top 2 surrogate variables that 
# will be used in DESeq2() but the normalized counts from ddssva object are NOT 
# batch corrected. We do this below.  https://www.biostars.org/p/121489/
# NOTE: Lowly expressed genes are removed before finding surrogate variables.
# So, number of genes is lower than number of DESeq2 normalized counts excel.
norm_counts_sva <- function(meta_data, read_data, formula_string, output_path){
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=read_data,
                                        colData=meta_data, 
                                        design=~1)
  dds <- DESeq2::estimateSizeFactors(dds) 
  dat  <- DESeq2::counts(dds, normalized = TRUE)
  
  # Remove lowly expressed genes before finding surrogate variables
  idx  <- rowMeans(dat) > 1
  dat  <- dat[idx, ]
  
  # Full model matrix with the variable of interest #Reduced/null model matrix with only an intercept term
  mod  <- stats::model.matrix(as.formula(formula_string), colData(dds))
  # Reduced/null model matrix with only an intercept term
  mod0 <- stats::model.matrix(~ 1, colData(dds))
  # Estimate all surrogate variables
  svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0)
  
  # If there are more than 3 SV, estimate upto 3 SVs
  # svseq$sv values will change based on n.sv
  if (svseq$n.sv > 3){
    svseq <- sva::svaseq(dat=dat, mod=mod, mod0=mod0, n.sv=3)
  }
  
  # %*% indicates Matrix multiplication
  X <- base::cbind(mod, svseq$sv)
  Hat <- base::solve(t(X) %*% X) %*% t(X)
  beta <- (Hat %*% t(dat))
  P <- ncol(mod)
  corrected_data <- dat - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])
  
  normalized_counts <- corrected_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID")
  
  normalized_counts <- add_annotation(normalized_counts, annotations) 
  
  # Save batch corrected normalized counts for entire dataset
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName="Norm_counts_batch_corrected")
  openxlsx::writeData(wb, sheet="Norm_counts_batch_corrected", x=normalized_counts_batch, rowNames=FALSE)
  openxlsx::saveWorkbook(wb, file=paste0(output_path, "Normalized.counts.SVA.xlsx"), 
                         overwrite=TRUE)
}

#******************************************************************************#
#                              CALCULATE padj AND log2FC
#******************************************************************************#

# Function to calculate pval and log2FoldChange
# norm_counts is matrix with log2 transformed values or non-log transformed values
# DO NOT use log10 transformed values
calc_stats <- function(norm_counts, metadata, Target, Reference, log2_transformed_already){
  
  # Perform t.test
  SYMBOL <- c()
  expt <- c()
  control <- c()
  pval <- c()
  for (j in 1:nrow(norm_counts)){
    
    data <- norm_counts[j,] %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Sample") %>%
      dplyr::left_join(metadata, by = c("Sample" = "Sample"))
    colnames(data) <- c("Sample", "values", "Condition")
    
    # Use F.test() to determine if variances are equal or not.
    # NOTE: p < 0.05 => unequal variance
    
    # NOTE: If Iso <- c(NA, NA, 18) and IP <- c(17, NA, 18),
    # var.test and t.test with throw error since Iso has only 1 value.
    
    # NOTE: If Iso <- c(1,1,1) and IP <- c(1,1,1),
    # t.test will throw error "data are essentially constant".
    
    # NOTE: If Iso <- c(0,0,0) and IP <- c(1,1,1),
    # t.test will throw error "data are essentially constant".
    
    if (sum(!is.na(data[data$Condition == Reference, ]$values)) > 1 &
        sum(!is.na(data[data$Condition == Target, ]$values)) > 1 & 
        (length(unique(data[data$Condition == Reference, ]$values)) + 
         length(unique(data[data$Condition == Target, ]$values)) > 2)){
      #   f_test <- var.test(formula = values ~ Condition, 
      #                      data = data,
      #                      alternative = "two.sided")
      #   
      #   if(!is.na(f_test$p.value)){
      # if (f_test$p.value < 0.05){
      # t_test <- t.test(formula = values ~ Condition,
      #                  data = data,
      #                  alternative = "two.sided",
      #                  var.equal = FALSE)
      # }
      
      # # Remove outliers
      # ref_data <- data[data$Condition == Reference, ]$values
      # low_ref <- quantile(ref_data, na.rm=TRUE)[2] - 1.5*IQR(ref_data, na.rm=TRUE)
      # high_ref <- quantile(ref_data, na.rm=TRUE)[3] + 1.5*IQR(ref_data, na.rm=TRUE)
      # 
      # target_data <- data[data$Condition == Target, ]$values
      # low_tar <- quantile(target_data, na.rm=TRUE)[2] - 1.5*IQR(target_data, na.rm=TRUE)
      # high_tar <- quantile(target_data, na.rm=TRUE)[3] + 1.5*IQR(target_data, na.rm=TRUE)
      # 
      # data <- data %>%
      #   dplyr::filter(!(Condition == Reference & (values > high_ref | values < low_ref))) %>%
      #   dplyr::filter(!(Condition == Target & (values > high_tar | values < low_tar)))
      
      
      # Calculate p values, mean expression
      t_test <- stats::t.test(formula = values ~ Condition, 
                              data = data,
                              alternative = "two.sided",
                              var.equal = FALSE)
      
      if (grepl(Reference, names(t_test$estimate[1]))){
        SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[1]])
        expt <- c(expt, t_test$estimate[[2]])
      } else if (grepl(Reference, names(t_test$estimate[2]))) {
        SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
        pval <- c(pval, t_test$p.value)
        control <- c(control, t_test$estimate[[2]])
        expt <- c(expt, t_test$estimate[[1]])
      }
    } else{
      SYMBOL <- c(SYMBOL, rownames(norm_counts)[j])
      pval <- c(pval, 1) # Note: DO NOT SET to NA. It will increase padj.
      control <- c(control, mean(data[data$Condition == Reference, ]$values, na.rm=TRUE))
      expt <- c(expt, mean(data[data$Condition == Target, ]$values, na.rm=TRUE))
    }
  }
  
  stats_df <- data.frame(SYMBOL, expt, control, pval)
  stats_df$padj <- stats::p.adjust(p = stats_df$pval, method = "fdr", n = length(stats_df$pval))
  if(log2_transformed_already){
    stats_df$log2FoldChange <- stats_df$expt - stats_df$control  # if data is already log transformed
  }else{
    stats_df$log2FoldChange <- log(stats_df$expt/stats_df$control, base=2)
  }
  
  result <- norm_counts %>% 
    tibble::rownames_to_column(var = "SYMBOL") %>% 
    dplyr::left_join(stats_df, by = c("SYMBOL" = "SYMBOL"))
  
  return(result)
}

flattenCorrMatrix_pmatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  df <- data.frame(row = rownames(cormat)[row(cormat)[ut]],
                   column = rownames(cormat)[col(cormat)[ut]],
                   cor  =(cormat)[ut],
                   p = pmat[ut])
  
  return(df)
}

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  df <- data.frame(row = rownames(cormat)[row(cormat)[ut]],
                   column = rownames(cormat)[col(cormat)[ut]],
                   cor  =(cormat)[ut])
  return(df)
}

# Compare two dataframe and output similarity
compare_deg_results <- function(df1, df2, file_suffix, output_path){
  
  df1 <- df1 %>% 
    dplyr::select(SYMBOL, ENSEMBL_ID, ENSEMBL_SYMBOL, log2FoldChange, padj) %>%
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::mutate(padj = round(padj, 2), log2FoldChange = round(log2FoldChange, 2)) 
  
  df2 <- df2 %>% 
    dplyr::select(SYMBOL, ENSEMBL_ID, ENSEMBL_SYMBOL, log2FoldChange, padj) %>%
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::mutate(padj = round(padj, 2), log2FoldChange = round(log2FoldChange, 2)) 
  
  merged_df <- dplyr::full_join(df1, df2,by=c("ENSEMBL_ID"="ENSEMBL_ID")) %>%
    dplyr::mutate(SYMBOL = dplyr::case_when(!is.na(SYMBOL.x) ~ SYMBOL.x,
                                            !is.na(SYMBOL.y) ~ SYMBOL.y,
                                            TRUE ~ ENSEMBL_ID)) %>%
    dplyr::select(SYMBOL, ENSEMBL_ID, log2FoldChange.x,  log2FoldChange.y, padj.x, padj.y) %>%
    dplyr::mutate(Group = dplyr::case_when(padj.x <= 0.05 & padj.y <= 0.05 & log2FoldChange.x >= 0.58 & log2FoldChange.y >= 0.58 ~ "Up in both",
                                           padj.x <= 0.05 & padj.y <= 0.05 & log2FoldChange.x <= -0.58 & log2FoldChange.y <= -0.58 ~ "Down in both",
                                           padj.x <= 0.05 & log2FoldChange.x >= 0.58 & (padj.y > 0.05 | is.na(padj.y)) ~ "Up x",
                                           padj.x <= 0.05 & log2FoldChange.x <= -0.58 & (padj.y > 0.05 | is.na(padj.y)) ~ "Down x",
                                           (padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y >= 0.58 ~ "Up y",
                                           (padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y <= -0.58 ~ "Down y",
                                           TRUE ~ "Not Significant in both"))
  
  up.x.y <- nrow(merged_df %>% dplyr::filter(log2FoldChange.x >= 0.58, log2FoldChange.y >= 0.58, padj.x <=0.05, padj.y <= 0.05))
  down.x.y <- nrow(merged_df %>% dplyr::filter(log2FoldChange.x <= -0.58, log2FoldChange.y <= -0.58, padj.x <=0.05, padj.y <= 0.05))
  up.x <- nrow(merged_df %>% dplyr::filter(padj.x <= 0.05 & log2FoldChange.x >= 0.58 & (padj.y > 0.05 | is.na(padj.y))))
  down.x <- nrow(merged_df %>% dplyr::filter(padj.x <= 0.05 & log2FoldChange.x <= -0.58 & (padj.y > 0.05 | is.na(padj.y))))
  up.y <- nrow(merged_df %>% dplyr::filter((padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y >= 0.58))
  down.y <- nrow(merged_df %>% dplyr::filter((padj.x > 0.05 | is.na(padj.x)) & padj.y <= 0.05 & log2FoldChange.y <= -0.58))               
  
  
  ggplot2::ggplot(data = merged_df, 
                  mapping = aes(x=log2FoldChange.x, y = log2FoldChange.y, color = Group)) +
    geom_point(size=0.75) +
    theme_classic() +
    coord_cartesian(clip = "off") +
    geom_hline(yintercept = c(-0.58, 0.58), color="black", linetype="dashed") +
    geom_vline(xintercept = c(-0.58, 0.58), color="black", linetype="dashed") +
    annotate(geom="text", label=up.x.y, x=10, y=10, col="black", size=5) +
    annotate(geom="text", label=down.x.y , x=-10, y=-10, col="black", size=5) +
    annotate(geom="text", label=up.x, x=10, y=0, col="black", size=5) +
    annotate(geom="text", label=down.x, x=-10, y=0, col="black", size=5) +
    annotate(geom="text", label=up.y, x=0, y=10, col="black", size=5) +
    annotate(geom="text", label=down.y, x=0, y=-10, col="black", size=5)
  
  ggplot2::ggsave(filename = paste0(file_suffix, ".jpg"),
                  plot = last_plot(),
                  device = "tiff",
                  path = output_path,
                  width = 7,
                  height = 7,
                  units = c("in"),
                  dpi = 300,
                  limitsize = TRUE,
                  bg = NULL)
}

#******************************************************************************#
#                         MICROARRAY RELATED FUNCTIONS                         #                       
#******************************************************************************#

# Values should be raw i.e. untransformed. DO NOT use log transformed values etc.
# DO NOT replace NA with 0 etc. Leave NA as they are.
# NA values will be imputed/replaced with average of non-NA values.
# If all values are NA, they will be set to 0.
# NO duplicated genes MUST be present.
impute_with_mean <- function(raw_counts){
  
  # Replace NA with average
  for (j in 1:nrow(raw_counts)){
    
    data <- raw_counts[j,] %>% 
      t() %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "Sample") %>%
      dplyr::left_join(metadata, by = c("Sample" = "Sample"))
    colnames(data) <- c("Sample", "values", "Condition")
    
    # If all values for Reference == NA or Target == NA, set them to 0. 
    # Replace NA with average wherever possible.
    data <- data %>% 
      dplyr::mutate(values = as.numeric(values)) %>%
      dplyr::group_by(Condition) %>% 
      dplyr::mutate(average = mean(values, na.rm=TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(average = dplyr::if_else(is.na(average), 0, average),
                    values = dplyr::if_else(is.na(values), average, values))
    
    data <- data %>% 
      dplyr::select(Sample, values) %>% 
      tibble::column_to_rownames("Sample") %>% 
      t()
    
    if(all(colnames(data) == colnames(raw_counts))){
      raw_counts[j,] <- data[,colnames(raw_counts[j,])] %>% unlist(use.names=FALSE)
    }
  }
  
  imputed_counts <- raw_counts %>% 
    dplyr::mutate(across(.cols=everything(), .fns = as.numeric))
  
  return(imputed_counts)
}

# Perform normalization before imputation so counts can be compared across samples
# NOTE: meta_data MUST have column "Sample" that matches with column names of
# raw_counts 
# NOTE: meta_data MUST have column "Condition" that defines the groups
quantile_norm <- function(raw_counts, meta_data, quant_norm){
  
  # Perform quantile normalization
  if (quant_norm == TRUE){
    
    # https://doi.org/10.1038/s41598-020-72664-6
    # The above paper recommends to quantile normalize within each group rather
    # than whole dataset
    quant_norm_counts <- data.frame(matrix(data=NA, nrow=nrow(raw_counts), ncol=0))
    
    # Subset raw_counts for each group
    for (c in unique(meta_data$Condition)){
      
      samples <- meta_data %>% 
        dplyr::filter(Condition %in% c) %>% 
        dplyr::select(Sample) %>% 
        unlist(use.names=FALSE)
      
      counts <- raw_counts[,samples]
      counts <- as.data.frame(preprocessCore::normalize.quantiles(as.matrix(counts)))
      rownames(counts) <- rownames(raw_counts)
      colnames(counts) <- samples
      
      quant_norm_counts <- dplyr::bind_cols(quant_norm_counts, counts)
    }
  } else {
    quant_norm_counts <- raw_counts
  }
  
  return(quant_norm_counts)
}

#******************************************************************************#
#                          NOTES ON PATHWAY ANALYSIS                           #
#******************************************************************************#

# Enrichment analysis takes differential data from every measured gene and looks
# for pathways displaying significantly coordinated shifts in those values.
# https://www.biostars.org/p/12182/
# https://www.biostars.org/p/17628/
# https://www.reddit.com/r/bioinformatics/comments/11o7mrv/gene_set_enrichment_analysis_do_you_separate_out/?rdt=59356
# https://groups.google.com/g/gsea-help/c/oXsBOAUYnH4
# https://www.biostars.org/p/132575/
# https://support.bioconductor.org/p/85681/

#******************************************************************************#
#               NOTES ON DESEQ2 & RNA SEQ BATCH CORRECTION                     #
#******************************************************************************#

# NOTE: coeff MUST match one of columns in resultsNames(dds)
# betaPrior: default=FALSE, shrunken LFCs are obtained later using lfcShrink
# Perform lfcshrinkage to account for variability between replicates
# For ashr, if res is provided, then coef and contrast are ignored.
# lfcshrinkage will not change the number of DEGs and affects only logFC

# NOTE: It is RECOMMENDED to perform batch correction ONLY if you know the batch
# information for all the samples in meta_data.
# https://support.bioconductor.org/p/76099/
# https://support.bioconductor.org/p/9149116/
# https://support.bioconductor.org/p/133222/#133225/
# https://support.bioconductor.org/p/125386/#125387
# https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#using-sva-with-deseq2

# There are multiple approaches to batch correction:
# (i) Modelling batch effect using DESeq2 design [RECOMMENDED]
# Add Batch to design of DESEq2 like design=~Batch + id
# (ii) Using combatseq to remove batch effects
# Use combatseq to get batch corrected read counts and then use the corrected 
# counts as input for DESeq2() and a design without batch variable (design=~ id)
# (iii) Using sva to identify hidden batch variables and model them in DeSeq2
# Use sva to find upto 2 surrogate variables that could be causing batch effects
# and include them in the design of DESeq2 like design=~SV1 + SV2 + id

# NOTE: I have tested (i) and (ii) and found that 
# (a) Almost all DEGs identified as significant in (i) are present in (ii) 
# (b) All significant DEGs with -ve log2FC from (i) also have -ve log2FC in (ii)
# (c) All significant DEGs with +ve log2FC from (i) also have +ve log2FC in (ii)
# (d) (ii) identifies many other significant DEGs missing in (i)
# (e) log2FC from (ii) differs from (i) for all genes but difference is minimal
# for most of the significant DEGs

# NOTE: I have tested (i) and (iii) and found that
# (a) Majority of significant DEGs in (iii) match with (i) but there are many
# significant DEGs in (i) absent in (iii) and vice versa
# (b) Significant DEGs with -ve log2FC from (i) also have -ve log2FC in (iii)
# (c) Significant DEGs with +ve log2FC from (i) also have +ve log2FC in (iii) 
# (d) (iii) identifies many other significant DEGs missing in (i) but (iii)
# fails to identify many DEGs identified by (i)
# (e) log2FC from (iii) differs from (i) for MOST genes to a great extent,
# HOWEVER, DDR2KO samples had log2FC(DDR2) of -0.7 in (iii) but only -0.37 in 
# (i) and (ii). So, sva might infact be removing batch effects and detecting 
# true biological effects !!! 

#******************************************************************************#
#                         NOTES ON SINGLE CELL ANALYSIS                        #
#******************************************************************************#

# https://ggplot2.tidyverse.org/reference/guide_colourbar.html
# https://stackoverflow.com/questions/56777529/how-to-pass-bash-variable-into-r-script
# https://github.com/satijalab/seurat/issues/4082
# Major changes in Seurat v5 
# https://satijalab.org/seurat/articles/seurat5_integration
# https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html
# https://satijalab.org/seurat/articles/integration_mapping
# https://satijalab.org/seurat/articles/integration_introduction
# https://satijalab.org/seurat/reference/integratelayers
# https://rdrr.io/r/base/split.html
# https://satijalab.org/seurat/reference/aggregateexpression
# https://github.com/satijalab/seurat/issues/2101

# # Indicate if data is from human or mice. We will adjust gene names accordingly.
# species <- dplyr::case_when(proj %in% c("scRNASeq_Chen",
#                                         "scRNASeq_Simon",
#                                         "visium_GSE171351",
#                                         "scRNASeq_HRA003620",
#                                         "scRNASeq_GSE222315") ~ "Homo sapiens",
#                             TRUE ~ "Mus musculus")

#******************************************************************************#
#        NORMALIZE DATA, IDENTIFY HIGHLY VARIABLE FEATURES, SCALE DATA,        # 
#                PERFORM DIMENSIONAL REDUCTION USING PCA & UMAP                #
#******************************************************************************#

# Use the sctransform method as a more accurate method of normalizing, 
# estimating the variance of the filtered data, and identifying the most 
# variable genes. By default, sctransform accounts for cellular sequencing 
# depth (i.e. nUMIs). Also, we can regress out variation from cell cycle genes
# and mitochondrial genes if needed. 

# Refer https://satijalab.org/seurat/articles/sctransform_vignette.html
# The residuals (normalized values) are stored in pbmc[["SCT"]]@scale.data and 
# used directly as input to PCA. Please note that this matrix is non-sparse, and
# can therefore take up a lot of memory if stored for all genes. To save memory,
# we store these values only for variable genes, by setting the 
# return.only.var.genes=TRUE by default in the SCTransform().

# To assist with visualization and interpretation, we also convert Pearson 
# residuals back to ‘corrected’ UMI counts. You can interpret these as the UMI 
# counts we would expect to observe if all cells were sequenced to the same depth.
# The ‘corrected’ UMI counts are stored in pbmc[["SCT"]]@counts. 

# The log-normalized versions of these corrected counts are stored in 
# pbmc[["SCT"]]@data, which are very helpful for visualization.

# You can use the corrected log-normalized counts for differential expression
# and integration. However, in principle, it would be most optimal to perform
# these calculations directly on the residuals (stored in the scale.data slot) 
# themselves.

#******************************************************************************#
#               PREPARE THE DATA FOR INTEGRATION & INTEGRATE DATA              #
#******************************************************************************#

# As you see from the UMAP, the cells cluster differently in each sample. 
# To find the same cell population (say macrophages) between 2 samples,
# it is necessary for both samples to have similar clustering pattern in UMAP.
# So, we have to integrate the samples. 

# The goal of integration is to ensure that cell types of one 
# condition/dataset align with the same cell types of the other 
# conditions/datasets (e.g. macrophages in control samples align with 
# macrophages in stimulated condition).

# To integrate, we will use the shared highly variable genes from each 
# condition identified using SCTransform, then, we will "integrate" or 
# "harmonize" the conditions to overlay cells that are similar or have a 
# "common set of biological features" between groups. 

# STEP 7A: DECLARE REFERENCE SAMPLES FOR INTEGRATING THE DATA
# STEP 7B: SELECT 3000 MOST VARIABLE GENES TO USE FOR INTEGRATING THE DATA
# STEP 7C: FIND RESIDUALS FOR MISSING GENES
# Each sample has different 3000 most variable genes. Gene X which is most 
# variable among cells of "sample A" may not be one of the top 3000 most 
# variable genes in "sample B". PrepSCTIntegration() will calculate Pearson 
# residuals for missing genes so that all samples have the same 3000 genes

# STEP 7D: FIND COMMON ANCHORS BETWEEN SAMPLES TO INTEGRATE THE DATA
# NOTE: Data must be scaled & PCA must have been run before doing cca or rpca
# in this step. cca is computationally intensive if more than 2 samples are 
# integrated. In such cases, use "rpca". Also, using reference based integration
# is faster.

# (i) Perform canonical correlation analysis (CCA):
# CCA identifies shared sources of variation between the conditions/groups. It
# is a form of PCA, in that it identifies the greatest sources of variation in
# the data, but only if it is shared or conserved across the conditions/groups
# (using the 3000 most variant genes from each sample). This step roughly aligns
# the cells using the greatest shared sources of variation.

# NOTE: The shared highly variable genes are used because they are the most 
# likely to represent those genes distinguishing the different cell types 
# present.

# (ii) Identify anchors or mutual nearest neighbors (MNNs) across datasets 
# (sometimes incorrect anchors are identified): MNNs can be thought of as 
# 'best buddies'. For each cell in one condition:   
# (a) The cell's closest neighbor in the other condition is identified based on
# gene expression values - it's 'best buddy'.
# (b) The reciprocal analysis is performed, and if the two cells are 'best 
# buddies' in both directions, then those cells will be marked as anchors to 
# 'anchor' the two datasets together.

# NOTE: The difference in expression values between cells in an MNN pair 
# provides an estimate of the batch effect, which is made more precise by 
# averaging across many such pairs. A correction vector is obtained and applied
# to the expression values to perform batch correction."
# 
# (iii) Filter anchors to remove incorrect anchors:
# Assess the similarity between anchor pairs by the overlap in their local 
# neighborhoods (incorrect anchors will have low scores)

# STEP 7E: FIND OPTIMUM k.weight FOR USE IN Seurat::IntegrateData()
# k.weight MUST be less than number of anchors. Else, error will be thrown.

# STEP 7F: INTEGRATE THE DATA
# Use anchors and corresponding scores to transform the cell expression values,
# allowing for the integration of the conditions/datasets (different samples, 
# conditions, datasets, modalities)

# NOTE: Transformation of each cell uses a weighted average of the two cells of 
# each anchor across anchors of the datasets. Weights determined by cell 
# similarity score (distance between cell and k nearest anchors) and anchor 
# scores, so cells in the same neighborhood should have similar correction values.

# If cell types are present in one dataset, but not the other, then the cells 
# will still appear as a separate sample-specific cluster.

# STEP 7G: RUN PCA USING 3000 INTEGRATION FEATURES & UMAP USING FIRST 40 PCs
# You need to run PCA and UMAP after integration in order to visualize correctly
# because IntegrateData() uses a different set of 3000 variable genes. So, new
# PCs will need to be calculated.
# Note: If you used SCTransform() before integration, you don't need to run 
# ScaleData() after integration. However, if you ONLY used NormalizeData() 
# before integration, you need to use ScaleData() after integration.

#************************STEP 7B: INTEGRATE THE DATA*************************#

# NOTE: The work of SelectIntegrationFeatures(), PrepSCTIntegration(), 
# FindIntegrationAnchors() and IntegrateData() are done by IntegrateLayers().
# Additionally, a new reduction which is equivalent of RunPCA() is also 
# created after integration.

# NOTE: RPCA needs proper kweight. Else, it throws error. I have not yet found
# a way to calculate optimal kweight unlike seurat v3. If script gives error
# regarding kweight, use the kweight it recommends in the error and re-run.

#******************************************************************************#
#                 CLUSTER THE CELLS & REMOVE SCARCE CLUSTERS                   #
#******************************************************************************#

# FindNeighbors() uses the user indicated "reduction" to calculate the k-nearest
# neighbors and construct the SNN graph.
# FindClusters() then performs graph-based clustering on the SNN graph. 

# NOTE: It is recommended to adjust k.param of FindNeighbors() [default=20] to 
# the same value as n.neighbors of UMAP() [default=30] 
# https://github.com/satijalab/seurat/issues/2152

#**************************STEP 8C: MERGE ALL LAYERS*************************#

# Once integrative analysis is complete, you can rejoin the layers - which 
# collapses the individual datasets together and recreates the original 
# counts and data layers. You will need to do this before performing any 
# differential expression analysis. However, you can always resplit the 
# layers in case you would like to reperform integrative analysis.

#******************************************************************************#


# ### Visualize the number of UMIs per cell
# umi_qc <- function(meta){
#   
#   umi_cutoff <- 500
#   # metadata <- meta %>%
#   #   dplyr::mutate(QC = factor(QC, levels = qc_levels))
#   
#   # p <- ggplot(data = metadata, aes(x = Sample, y = nUMIs, fill = QC)) +
#     # geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) + 
#     # geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) +
#     # theme_classic() +  
#     # custom_theme +
#     #labs(x = "Sample", y = "Number of UMIs", title = "Distribution of UMIs") +
#     #coord_cartesian(ylim = c(1,1000000), clip = "off") +
#     #scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000)) +
#     # scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
#     #                              "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
#     #geom_hline(yintercept = umi_cutoff, linetype = 2)
#   
#   return(p)
# }

# ### Visualize the number of genes per cell
# gene_qc <- function(meta){
#   
#   gene_cutoff <- 250
#   # metadata <- meta %>%
#   #   dplyr::mutate(QC = factor(QC, levels = c("Empty Droplet","Doublet", "Low Quality", "Singlet")))
#   
#   # p <- ggplot(data = metadata, aes(x = Sample, y = nGenes, fill = QC)) +
#   #   geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) + 
#   #   geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) +
#   #   theme_classic() + 
#   #   custom_theme + 
#   #  labs(x = "Sample", y = "Number of Genes", title = "Distribution of Genes") +
#     coord_cartesian(ylim = c(1, 30000), clip = "off") +
#     scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000)) + 
#     # scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
#     #                              "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
#     #geom_hline(yintercept = gene_cutoff, linetype = 2)
#   
#   return(p)
# }

# # Visualize the MitoRatio of each cell
# mito_qc <- function(meta){
#   
#   mito_cutoff <- 0.2
#   # metadata <- meta %>%
#   #   dplyr::mutate(QC = factor(QC, levels = c("Empty Droplet","Doublet", "Low Quality", "Singlet")))
#   
#   # p <- ggplot(data = metadata, aes(x = Sample, y = MitoRatio, fill = QC)) +
#   #   geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) + 
#   #   geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) + 
#   #   theme_classic() + 
#   #   custom_theme +
#     #labs(x = "Sample", y = "MitoRatio", title = "Distribution of MitoRatio") +
#     coord_cartesian(ylim = c(0.00001, 1), clip = "off") +
#     scale_y_log10(breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1)) + 
#     # scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
#     #                              "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
#     #geom_hline(yintercept = mito_cutoff, linetype = 2)
#   
#   return(p)
# }

# ### Visualize the RiboRatio of each cell
# ribo_qc <- function(meta){
#   
#   ribo_cutoff <- 0.05
#   # metadata <- meta %>%
#   #   dplyr::mutate(QC = factor(QC, levels = c("Empty Droplet","Doublet", "Low Quality", "Singlet")))
#   
#   # p <- ggplot(data = metadata, aes(x = Sample, y = RiboRatio, fill = QC)) +
#   #   geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) + 
#   #   geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) + 
#   #   theme_classic() + 
#   #   custom_theme +
#    # labs(x = "Sample", y = "RiboRatio", title = "Distribution of RiboRatio") +
#     coord_cartesian(ylim = c(0.0001, 1), clip = "off") +
#     scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1)) + 
#     # scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
#     #                              "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
#     #geom_hline(yintercept = ribo_cutoff, linetype = 2)
#   
#   return(p)
# }

# ### Visualize the novelty or complexity of each cell
# novelty_qc <- function(meta){
#   
#   novelty_cutoff <- 0.8
#   # metadata <- meta %>%
#   #   dplyr::mutate(QC = factor(QC, levels = c("Empty Droplet","Doublet", "Low Quality", "Singlet")))
#   
#   # p <- ggplot(data = metadata, aes(x = Sample, y = Novelty, fill = QC)) +
#   #   geom_violin(position = position_dodge(0.9), scale = "width", drop = FALSE) + 
#   #   geom_boxplot(position = position_dodge(0.9), width = 0.15, outlier.size = 0.5, drop = FALSE) + 
#   #   theme_classic() + 
#   #   custom_theme +
#     # labs(x = "Sample", y = "Novelty", title = "Distribution of Novelty Score") +
#     coord_cartesian(ylim = c(0.3, 1), clip = "off") +
#     scale_y_log10(breaks = c(0.3, 1)) + 
#     # scale_fill_manual(values = c("Empty Droplet" = "#FFC61E", "Singlet" = "#2CA02C",
#     #                              "Doublet" = "#1F77B4","Low Quality" = "#D62728")) +
#     #geom_hline(yintercept = novelty_cutoff, linetype = 2)
#   
#   return(p)
# }

# # Plot all QC metrics before and after QC
# funcs <- c("cell_qc", "umi_qc", "gene_qc", "mito_qc", "ribo_qc", "novelty_qc",
#            "gene_umi_mito_qc")
# 
# filenames <- c("Cell_Counts", "UMI_Distribution", "Gene_Distribution",
#                "MitoRatio_Distribution", "RiboRatio_Distribution", 
#                "Novelty_Score_Distribution", "Genes_UMI_MitoRatio_Distribution")
# 
#for (i in 1:length(funcs)){

# # Plot QC metrics
# purrr::map(.x = c("raw_metadata"), .f = get(funcs[i])) %>% 
#   cowplot::plot_grid(plotlist = ., align = "hv", axis = "tblr", nrow = 1, ncol = 1)

#p <- get(funcs[i])(raw_metadata)

#   # Save the plot
#   ggplot2::ggsave(filename = paste0("QC_", filenames[i], ".pdf"),
#                   plot = p, #last_plot(),
#                   device = "pdf",
#                   path = output_path,
#                   #scale = 1,
#                   width = 11,
#                   height = 8,
#                   units = c("in"),	 
#                   dpi = 600,
#                   limitsize = TRUE,
#                   bg = NULL)
# }

#cat("QC plots generated\n")  
#}

#******************************************************************************#

# Create a list of sample names which will be added to each barcode.
# Since folder names correspond to sample name, we just use list.files()
samples <- list.files(path = raw_matrix_path) 

# Create empty dataframe to store raw metadata
raw_metadata <- data.frame(Cell = c("")) #Sample = as.factor(1), nUMIs = c(0)

for (s in samples){
  s.obj <- read_cellranger(s, raw_matrix_path)
  s.obj <- mark_empty_droplets_dropletutils(s.obj)
  s.obj <- mark_empty_droplets_cellranger(s.obj)
  s.obj <- doublet_finder(s.obj)
  s.obj <- scdbl_finder(s.obj)
  s.obj <- calc_qc_metrics_sc_sp(s.obj, "RNA")
  s.obj <- mark_low_quality_sc_sp(s.obj)
  raw_metadata <- generate_plotdata(s.obj, raw_metadata)
  s.obj <- filter_singlets_sc_sp(s.obj)
  assign(s, s.obj)
}

plot_qc(raw_metadata, seurat_results)
#xpectr::suppress_mw(generate_whitelist(filt, seurat_results))

# Merge all samples
# NOTE: seurat objects MUST have been loaded into R prior to this step
# Import any other meta data associated with data set
# NOTE: This xslx file should have column named "Unique_ID" whose values matches 
# with  column "Unique_ID" of seurat object's metadata.
filt.obj <- merge_filtered_sc_sp(samples, "RNA", extra_metadata, seurat_results)
sct.obj  <- sctransform_sc_sp(filt.obj, "RNA", seurat_results)
reference.samples <- NULL
kweight <- min(sct.obj@meta.data %>% dplyr::count(Sample) %>% dplyr::select(n) %>% min()/2, 100) 
integ <- integrate_sc_sp(sct.obj, "RNA", reference.samples, kweight, seurat_results)
integ.clust <- cluster_sc_sp(integ, "RNA", seurat_results)
integ.final <- remove_sparse_clusters_sc_sp(integ.clust, "RNA", seurat_results)


suffix <- "Full"
resolution <- 0.8
reduction <- "Harmony"
identify_markers_sc_sp(integ.final, "RNA", resolution, reduction, suffix, seurat_results)
plot_metrics_post_integration_sc_sp(integ.final, "RNA", diagnostics_path)

# Decide which reduction gives (i) clear separation of clusters 
# (ii) clusters cells together i.e. B cells are not split between clusters.
# In all cases so far, RPCA and Harmony seemed to work best.