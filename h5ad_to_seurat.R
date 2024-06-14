# We export metadata from adata
# We export cell names i.e. barcodes from adata
# We import expr data using SeuratDisk so that we dont have to save expr data in
# csv files

qrsh
conda activate scVelo
conda install -c anaconda openpyxl
python

import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad("/hpc/home/kailasamms/scratch/scRNASeq_GSE132042/scRNASeq_GSE132042.h5ad")

# Info on expr 
adata.X
adata.to_df()

# Info on cells
adata.obs

# Info on genes
adata.var 

# Subsetting adata
# bdata = adata[adata.obs.cell_ontology_class == "bladder urothelial cell"]

# Convert expr data to df and save it. to_excel() is too slow. So, use to_csv()
# Recommend using SeuratDisk rather to import expr data
# df = pd.DataFrame(adata.X.toarray()).transpose()
# df.index = adata.var_names
# df.columns = adata.obs_names
# df.to_csv("/hpc/home/kailasamms/Bladder_droplet_counts.csv", index=True)

# Save metadata to df
adata.obs.to_excel("/hpc/home/kailasamms/scratch/scRNASeq_GSE132042/scRNASeq_GSE132042_Meta.xlsx", index=True)


#############################################################

# count_mat <- read.csv("/hpc/home/kailasamms/Bladder_droplet_counts.csv", header = TRUE) %>%
#   tibble::column_to_rownames("index")
# colnames(count_mat) <- gsub(pattern= "\\.", replacement= "-", x=colnames(count_mat))
# colnames(count_mat) <- gsub(pattern= "X10X", replacement= "10X", x=colnames(count_mat))
#rownames(metadata_df) <- metadata_df$Cell

# Load h5ad
SeuratDisk::Convert(source = paste0(parent_path,"scRNASeq_GSE132042.h5ad"),
                    dest = "h5seurat",
                    assay="RNA",
                    overwrite = FALSE)

raw_seurat <- SeuratDisk::LoadH5Seurat(file = paste0(parent_path,"scRNASeq_GSE132042.h5seurat"))

# This group has messed up metadata. There is no metadata.
# Correct the excel file to contain columns orig.ident, nCount_RNA, nFeature_RNA, 
# Sex, Condition, Age.
# Also, correct the barcodes like "AAACCTGAGTACGTTC-1-24-0-0", 
# "10X_P7_7_TTGCGTCGTAGGACAC-1" to follow the format "FR1_AAACCTGAGTACGTTC-1"
metadata <- openxlsx::read.xlsx(paste0(parent_path, "scRNASeq_GSE132042_Meta.xlsx")) %>% 
  dplyr::rename(orig.ident = mouse.id,
                nCount_RNA = n_counts,
                nFeature_RNA = n_genes) %>%
  dplyr::mutate(index = paste0(orig.ident, "_", index)) %>%
  tibble::column_to_rownames("index") %>%
  dplyr::select(orig.ident, nCount_RNA, nFeature_RNA)

# Ideally counts must be  integers but these idiots dont share this info.
# So, we have normalized values in counts
readdata <- raw_seurat@assays$RNA@counts
colnames(readdata) <- rownames(metadata)
# readdata <- (exp(readdata)-1)/10000
# collapse::TRA(x=as.matrix(readdata), STATS=metadata$nCount_RNA, FUN="*")
  
# Import these into a seurat object
raw_seurat <- CreateSeuratObject(counts = readdata,
                                 project = "GSE132042",
                                 assay = "RNA",
                                 meta.data = metadata,
                                 row.names = rownames(readdata))

# Follow steps from Step 3 Quality control making only the a small change
# while defining Sample in transmute()
raw_metadata <- raw_metadata %>% 
  dplyr::transmute(Cell = rownames(raw_metadata),
                   Folder = orig.ident,
                   Sample = Folder,
                   ...)
