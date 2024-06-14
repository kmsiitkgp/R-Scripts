# This is uploaded to Github.

library(RColorBrewer)
library(colorspace)

# Extract the default ggplot2 color scheme
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

# Save the default color scheme for 12 color palette
palette <- ggplotColours(n=12)
palette <- c("#FF1F5B", "#F28522", "#009ADE", "#AF58BA", "#3C005A", "#00B000", "#FFC61E", "#808080", "#A6761D", 
             "#F6D2E0", "#E75480", "#C8E7F5", "#2E5984", "#FFFF99", "#B15928", "#A8A9AD", "#000000") 
palette <- RColorBrewer::brewer.pal(n=8, name="Dark2")
palette <- RColorBrewer::brewer.pal(n=12, name="Paired")
palette <- RColorBrewer::brewer.pal(n=11, name="RdYlBu")
#palette <- RColorBrewer::brewer.pal(11, "RdYlBu")[c(1,11)]

# Print the color codes
palette

# Visualize the colors
pie(rep(1,length(palette)), col=palette)


library(RColorBrewer)
n <- 25
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'seq',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))


for (i in 1:ceiling(length(col_vector)/10)){
  setwd("C:/Users/KailasammS/Desktop/pals/")
  jpeg(paste0(i, ".jpg")) 
  
  start <- 10*i-9 
  end <- 10*i
  palette <- col_vector[start:end]
  pie(rep(1,length(palette)), col=palette)
  
  dev.off() 
}

# Go through the plots and create palette using picture# and color#
l <- data.frame(file_num = c(1, 1, 2, 2, 3, 4, 5, 5, 6, 7, 7, 8, 9, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 23),
                col_num =  c(3, 9, 3, 10, 4, 6, 3, 7, 3, 4, 6, 2, 9, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 7))

palette <- c()
for (i in 1:nrow(l)){
  
    palette <- c(palette, col_vector[(l$file_num[i]-1)*10+l$col_num[i]])
}

# Visualize the colors
pie(rep(1,length(palette)), col=palette)

# Print the color codes
palette
#"#F28522"
my_palette <- c("#FF1F5B", "#F28522", "#009ADE", "#AF58BA", "#00B000", 
                   "#FFC61E", "#808080", "#BF812D", "#35978F", "#C51B7D", 
                   "#7FBC41", "#762A83", "#D6604D", "#4393C3", "#FFFFBF", 
                   "#9E0142", "#E41A1C", "#4DAF4A", "#FF7F00", "#FFFF33", 
                   "#A65628", "#F781BF", "#66C2A5", "#FC8D62", "#000000")

my_palette <- c("#9E0142", "#C51B7D", "#009ADE", "#AF58BA", "#00B000", 
                   "#FFC61E", "#762A83", "#7FBC41", "#C51B7D", "#FF7F00",
                   "#FF1F5B", "#D6604D", "#808080", "#BF812D", "#000000", 
                   "#4393C3", "#FFFFBF", "#FFFF33", "#E41A1C", "#4DAF4A", 
                   "#35978F", "#A65628", "#F781BF", "#66C2A5", "#FC8D62")
#scanpy tab20
my_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", 
                   "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                   "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
                   "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5")

#modified scanpy tab20
my_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", 
                   "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                   "#FFC61E", "#762A83", "#333333", "#FF1F5B", "#B8E80C",
                   "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
                   "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5",
                   "#FFFFBF", "#C51B7D", 
                   "#DE77AE", "#F1B6DA", "#FDE0EF", "#F7F7F7", "#E6F5D0",   
                   "#B8E186")  #"#FFC3FF", "#C51B7D", "#FFFF72"


my_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", 
                "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
                "#FFC61E", "#762A83", "#333333", "#FF1F5B", "#B8E80C",
                "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
                "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5",
                "#FFFFBF", "#C51B7D", "#DE77AE", "#F1B6DA", "#FDE0EF", 
                "#F7F7F7", "#E6F5D0", "#B8E186")

my_palette1 <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#8C564B", 
                 "#E377C2", "#BCBD22", "#17BECF", "#FFC61E", "#762A83",
                 "#333333", "#FF1F5B", "#B8E80C")

my_palette1 <- c(my_palette1, 
                 colorspace::adjust_transparency(col = my_palette1, alpha = 0.6), 
                 colorspace::adjust_transparency(col = my_palette1, alpha = 0.3))
               

pie(rep(1,length(my_palette1)), col=my_palette1)
my_palette2 <- c()
for(i in 1:13){
  my_palette2 <- c(my_palette2, my_palette1[i], my_palette1[13+i], my_palette1[26+i])
}
pie(rep(1,length(my_palette2)), col=my_palette2)

my_palette2 <-c("grey", viridis(n = 10, option = "C", direction = -1))

my_palette <- c("#BF812D", "#35978F", "#C51B7D", "#7FBC41", "#762A83",
                "#E08214", "#542788", "#D6604D", "#4393C3", "#878787",
                "#1A1A1A", "#FFFFBF", "#377EB8", "#4DAF4A", "#984EA3",
                "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
                "#66C2A5", "#FC8D62", "#000000", "#9E0142", "#E41A1C")
pie(rep(1,length(my_palette)), col=my_palette)

library(colorspace)
hcl_palettes(plot = TRUE)
