install.packages("UpSetR")

library("UpSetR")

plot_example <- function(){ 
# The names of list will be on Y axis of bottom graph
# The intersection size based on Y values of top graph
# All possible intersections will be displayed in bottom graph
listInput <- list(DMPK=c(1,2,3,4,5,6,7,8,9,11,16,17),
              MAPK1=c(1,2,3,4,5,6,7,8,9,16,17,18),
              RAF1=c(1,2,3,4,5,6,7,8,9,16,17,18),
              ROCK1=c(1,2,3,4,5,6,7,8,9,17,18),
              PIM1=c(1,2,3,4,5,6,7,8,9,17),
              DYRK2=c(1,2,3,4,5,6,7,8,9,17),
              STK26=c(1,2,3,4,5,6,7,8,9,17),
              MAP2K2=c(1,2,3,4,5,6,7,8,9),
              LIMK1=c(1,2,3,4,5,6,7,8,9),
              MYO3A=c(1,2,3,4,5,6,7,8,9),
              TSSK1B=c(1,2,3,4,5,6,7,8,9),
              HIPK4=c(1,2,3,4,5,6,7,8,9),
              DCLK1=c(1,2,3,4,5,6,7,8,9),
              NEK3=c(1,2,3,4,5,6,7,8,9),
              CDKL1=c(1,2,3,4,5,6,7,8,9),
              TRPM7=c(1,2,3,4,5,6,7,8,10,11,12,13,14,15),
              MAPK7=c(1,2,3,4,5,6,7,8,17),
              IKBKE=c(1,2,4,5,6,7,8,9,17),
              PRKAG3=c(1,2,4,5,6,7,8,9),
              TAF1L=c(1,2,4,5),
              JAK2=c(2,4,5,6,7,8,9,16,17,18))


# Plot
ggplotify::as.ggplot(upset(data = UpSetR::fromList(listInput),
      empty.intersections = "on",
      cutoff = 5,
      mb.ratio = c(0.5,0.5),
      sets = c("DCLK1","MAPK1","CDKL1","ROCK1","MAPK7","HIPK4","MAP2K2","DYRK2"),
      #nintersects = 5,                  # number of groups on X axis
      #nsets = 21,                         # number of groups on Y axis
      order.by = c("freq")))

ggsave("C:/Users/kailasamms/OneDrive - Cedars-Sinai Health System/Desktop/Example.jpg",
       height = 11,
       width=11)
}
plot_example()


