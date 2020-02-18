# Why two clusters 

library(SOMbrero)
library(tidyverse)
library(oce)
library(rworldmap)
library(ggplot2)
library(ggmap)
library(cowplot)
library(scales)
library(leaps)
library(MASS)
library(randomForest)
library(grid)
library(gridExtra)
library(lubridate)
library(vegan)
library(spatialEco)
library(geosphere)
library(viridis)




result_tables <- function(input = "data/16s_bacteria.Rdata", 
                          physical_data = "data/NCOG_sample_metadata.csv",
                          full_data_file = "output/bacteria_16s_full_data.Rdata",
                          sample_regime = "both", name = "bact"){



load(input)

asv_copy <- asv_table

physical_dat <- read.csv(physical_data, header = TRUE, stringsAsFactors = FALSE)

# check variables

data <- physical_dat[colnames(physical_dat)[c(1,33:41,46:47)]]

data$eco_name <- paste0("X20",substr(data$Sample.Name,start = 3,stop = 100))

full_dat <- physical_dat[!is.na(match(physical_dat$Sample.Name, data$Sample.Name)),]

full_dat$eco_name <- paste0("X20",substr(full_dat$Sample.Name,start = 3,stop = 100))

# physical data match to biological data

full_dat <- full_dat[which(!is.na(match(full_dat$eco_name, rownames(scaled_inputs)))),]

data <- data[which(!is.na(match(data$eco_name, rownames(scaled_inputs)))),]

# partition data

if(sample_regime == "both"){}
if(sample_regime == "surface"){
  
  eco_vals <- full_dat$eco_name[which(full_dat$Depthm < 15)]
  
  scaled_inputs <- scaled_inputs[eco_vals,]  
  asv_copy <- asv_copy[eco_vals,]   
  asv_table <- asv_table[eco_vals,] 
  
  full_dat <- full_dat[which(full_dat$Depthm < 15),]
  
}

if(sample_regime == "depth"){
  
  eco_vals <- full_dat$eco_name[which(full_dat$Depthm > 14)]
  
  scaled_inputs <- scaled_inputs[eco_vals,]  
  asv_copy <- asv_copy[eco_vals,]   
  asv_table <- asv_table[eco_vals,] 
  
  full_dat <- full_dat[which(full_dat$Depthm > 14),]
  
}

som_kmeans_cascade <- cascadeKM(scaled_inputs, inf.gr = 2,
                                sup.gr = 15, 
                                iter = 1000,
                                criterion = "ssi")

summary(som_kmeans_cascade)
plot(som_kmeans_cascade, sortg = TRUE)


eco.som <- trainSOM(x.data = scaled_inputs, dimension = c(5, 5), nb.save = 10, maxit = 2000, 
                    scaling = "none")

som_kmeans_cascade <- cascadeKM(eco.som$prototypes, inf.gr = 2,
                                sup.gr = 15, 
                                iter = 1000,
                                criterion = "ssi")

summary(som_kmeans_cascade)
plot(som_kmeans_cascade, sortg = TRUE)


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
  plot(2:15,prop_mat[,i], type = "l", xlab = "Clusters", ylab = "Proportion of Clusters", main = paste0(name, " ", (i+1), " Clusters"), ylim = c(0,1))
  points(2:15,prop_mat[,(i+6)], type = "l", lty = 2, col = "red")

  }

# par(mfrow = c(1,1))
# 
# plot(2:15, prop_mat[,1]-prop_mat[,(1+6)], type = "l", xlab = "Clusters", 
#      ylab = "Difference", ylim = c(0,0.5))
# for (i in 2:6) {
#   points(2:15, prop_mat[,i]-prop_mat[,(i+6)], type = "l", xlab = "Clusters", 
#          ylab = "Difference")
# }


}





# 16s cyanos
set.seed(5)
result_tables(input = "data/16s_cyanos.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              full_data_file = "output/cyano_16s_full_data.Rdata",
              sample_regime = "both", name = "Cyanobacteria")

# 16s minus eukaryotes and plastids
set.seed(192)
result_tables(input = "data/16s_bacteria_m_euks.Rdata", 
              physical_data = "data/NCOG_sample_metadata.csv",
              full_data_file = "output/bacteria_m_euks_16s_full_data.Rdata",
              sample_regime = "both", name = "Bacteria/Archaea")

# 18sv9 Autotrophic Eukaryotes
set.seed(15)
result_tables(input = "data/18s_autotrophic_euks.Rdata", 
              physical_data = "data/NCOG_sample_metadata.csv",
              full_data_file = "output/euks_auto_18sv9_full_data.Rdata",
              sample_regime = "both", name = "Eukaryotic Phytoplankton")

# 18sv9 Eukaryotes
set.seed(987)
result_tables(input = "data/18s_heterotrophic_euks.Rdata", 
              physical_data = "data/NCOG_sample_metadata.csv",
              full_data_file = "output/euks_hetero_18sv9_full_data.Rdata",
              sample_regime = "both", name = "Heterotrophic Eukaryotes")
 