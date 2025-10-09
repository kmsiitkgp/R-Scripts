library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)


# User defined
organism <- "human"  #"mouse"

DEGs_df <- read.xlsx(.....)    # MUST have columsn SYMBOL, padj, log2FOldChange
vst_counts <- read.xlsx(...)   # log transformed normalized counts preferably vst or rlog
samples <- c()                 # samples belonging to groups being compared



tf_net <- decoupleR::get_collectri(organism = organism, 
                                   split_complexes = FALSE)
stats <- c("mlm", "ulm", "wsum")

# TF analysis using t-statistics
t_stats_mat <- DEGs_df %>% 
  as.data.frame() %>%
  dplyr::mutate(t = -log10(padj) * log2FoldChange) %>%
  dplyr::filter(!is.na(t)) %>%
  dplyr::select(SYMBOL, t) %>%
  tibble::column_to_rownames("SYMBOL") %>%
  as.matrix()

# TF analysis using counts
norm_counts_sub <- vst_counts[, samples]

# Run TF analysis 
# @group level (using t_stats_mat) or
# @sample level (using norm_counts_sub)
tf_df <- decoupleR::decouple(mat = ,         #t_stats_mat or norm_counts_sub, 
                             network = tf_net,
                             statistics = stats, 
                             minsize = 5)

# Plot bar plots for TF analysis @group level
ggplot(tf_df, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

# Plot heatmap for TF analysis @sample level
sample_acts_mat <- tf_df %>%
  pivot_wider(id_cols = 'condition', names_from = 'source', values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

# Scale per feature
sample_acts_mat <- scale(sample_acts_mat)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap(sample_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 

