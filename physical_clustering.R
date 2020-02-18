library(SOMbrero)
library(tidyverse)
library(oce)
library(lubridate)
library(vegan)


physical_data <- read.csv(file = "data/NCOG_sample_metadata.csv")

simple_phys <- physical_data[,c(33,34,37:39)]

simple_phys <- apply(simple_phys, 2, scale)

physical_data <- physical_data[complete.cases(simple_phys),]
simple_phys <- simple_phys[complete.cases(simple_phys),]

eco.som <- trainSOM(x.data = simple_phys, dimension = c(5, 5), nb.save = 10, maxit = 2000, 
                    scaling = "none")

prop_mat <- matrix(NA,14,12)

# running SOM
for(i in 2:15){
  
  eco.clust <- superClass(eco.som, k = i)
  
  clusters <- eco.clust$cluster
  
  ids <- eco.clust$som$clustering
  
  som_ids <- clusters[ids]
  
  prop_table <- (table(som_ids)/length(som_ids))[order(table(som_ids),decreasing = TRUE)]
  
  for (k in 2:7) {
    
    prop_mat[(i-1),(k-1)] <- sum(prop_table[1:k])
    
    prop_mat[(i-1),(6+(k-1))] <- k/i
    
  }
  
}

prop_mat[which(prop_mat > 1)] <- NA

par(mar = c(4,6,4,2), mfrow = c(2,3))
for (i in 1:6) {
  plot(2:15,prop_mat[,i], type = "l", xlab = "Clusters", ylab = "Proportion of Clusters", main = paste0("Physical Data", " ", (i+1), " Clusters"), ylim = c(0,1))
  points(2:15,prop_mat[,(i+6)], type = "l", lty = 2, col = "red")
  
}

som_kmeans_cascade <- cascadeKM(simple_phys, inf.gr = 2,
                                sup.gr = 15, 
                                iter = 10000,
                                criterion = "ssi")

summary(som_kmeans_cascade)
plot(som_kmeans_cascade, sortg = TRUE)


som_kmeans_cascade <- cascadeKM(eco.som$prototypes, inf.gr = 2,
                                sup.gr = 15, 
                                iter = 10000,
                                criterion = "ssi")

summary(som_kmeans_cascade)
plot(som_kmeans_cascade, sortg = TRUE)

