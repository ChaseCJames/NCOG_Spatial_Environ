library(SOMbrero)
library(tidyverse)
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
library(GUniFrac)
library(reldist)

result_tables <- function(input = "data/18s_diatom_S.Rdata",
                          som_file = "output/diatom_18sv9_som_S.Rdata",
                          physical_data = "data/NCOG_sample_log_DNA_stvx_meta_2014-2020.csv",
                          map_file = "output/diatom_18sv9_map_S.Rdata",
                          regression_file = "output/diatom_18sv9_glm_S.Rdata",
                          dissimmilar_matrix = "output/diatom_18sv9_dissimilar_S.Rdata", 
                          full_data_file = "output/diatom_18sv9_full_data_S.Rdata",
                          sample_regime = "both"){
  
  load(input)
  
  asv_copy <- asv_table
  
  physical_dat <- read.csv(physical_data, header = TRUE, stringsAsFactors = FALSE) 
  
  # remove samples if there is no data for a groupt
  
  # if(length((which(rowSums(asv_copy) == 0))) > 0){
  #   asv_copy <- asv_copy[-(which(rowSums(asv_copy) == 0)),]
  #   asv_table <- asv_table[-(which(rowSums(asv_table) == 0)),]
  #   scaled_inputs <- scaled_inputs[-(which(rowSums(scaled_inputs) == 0)),]
  # }

  # filter samples

  full_dat <- physical_dat %>% filter(paste0("X",Sample.Name) %in% rownames(asv_copy))
  
  # focus on main grid
  full_dat <- full_dat %>% filter(as.numeric(substr(Sta_ID,2,3)) > 75)
  
  asv_copy <- asv_copy[match(paste0("X",full_dat$Sample.Name), rownames(asv_copy)),]
  asv_table <- asv_table[match(paste0("X",full_dat$Sample.Name), rownames(asv_table)),]
  scaled_inputs <- scaled_inputs[match(paste0("X",full_dat$Sample.Name), rownames(scaled_inputs)),] %>% as.matrix()
  
  full_dat$eco_name <- paste0("X",full_dat$Sample.Name)
  
  # distance to coast
  
  full_dat$dist_to_coast <- abs(full_dat$Distance)
  
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
  
  # running SOM
  
  eco.som <- trainSOM(x.data = scaled_inputs, dimension = c(6, 6), nb.save = 10, maxit = 2000, 
                      scaling = "none")
  
  save(eco.som, file = som_file)
  
  # summary(eco.som)
  
  eco.clust <- superClass(eco.som, k = 2)
  
  clusters <- eco.clust$cluster
  
  ids <- eco.clust$som$clustering
  
  som_ids <- clusters[ids]
  
  asv_copy <- asv_copy[which(!is.na(match(rownames(asv_copy), rownames(scaled_inputs)))),]
  asv_table <- asv_table[which(!is.na(match(rownames(asv_table), rownames(scaled_inputs)))),]
  
  asv_copy$som_id <- som_ids
  
  full_dat$som_id <- asv_copy$som_id[match(full_dat$eco_name, rownames(asv_copy))]
  
  # shannon diversity
  
  asv_copy$shannon_index <- vegan::diversity(asv_table, MARGIN = 1, index = "shannon")
  full_dat$shannon <- asv_copy$shannon_index[match(full_dat$eco_name, rownames(asv_copy))]
  
  # evenness
  
  S <- apply(asv_table>0,1,sum)
  asv_copy$evenness <- vegan::diversity(asv_table, index="shannon")/log(S)
  full_dat$evenness <- asv_copy$evenness[match(full_dat$eco_name, rownames(asv_copy))]
  
  # richness
  
  richness <- apply(asv_table, 1, function(x) length(which(x != 0)))
  
  asv_copy$richness <- richness
  full_dat$richness <- asv_copy$richness[match(full_dat$eco_name, rownames(asv_copy))]
  
  # chao1
  chao1 <- estimateR(asv_table)
  
  asv_copy$chao1 <- chao1[2,]
  full_dat$chao1 <- asv_copy$chao1[match(full_dat$eco_name, rownames(asv_copy))]
  
  # gini
  
  gini <- apply(asv_table, 1,gini)
  
  asv_copy$gini <- gini
  full_dat$gini <- asv_copy$gini[match(full_dat$eco_name, rownames(asv_copy))]
  
  # inv_simp

  asv_copy$inv_simp <- vegan::diversity(asv_table, MARGIN = 1, index = "invsimpson")
  full_dat$inv_simp <- asv_copy$inv_simp[match(full_dat$eco_name, rownames(asv_copy))]
  
  
  # Dissimilarity
  
  dissimilar <- vegdist(asv_table, method = "bray", binary = FALSE)
  dissimilar <- as.matrix(dissimilar)
  
  save(dissimilar, file = dissimmilar_matrix)
  
  # maps
  
  som_maps <- full_dat %>% 
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som_id == 1, na.rm = TRUE)/n(), som_2 = sum(som_id == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE),
              temp_mean = mean(T_degC, na.rm = TRUE), sal_mean = mean(Salnty, na.rm = TRUE),
              PO4_mean = mean(PO4ug, na.rm = TRUE), NO3_mean = mean(NO3ug, na.rm = TRUE),
              SiO3_mean = mean(SiO3ug, na.rm = TRUE), 
              MLD_mean = mean(MLD_Sigma, na.rm = TRUE), NC_mean = mean(NCDepth, na.rm = TRUE),
              Dist_mean = mean(dist_to_coast, na.rm = TRUE), Chl_mean = mean(ChlorA, na.rm = TRUE),
              temp_coeff = sd(T_degC, na.rm = TRUE)/mean(T_degC, na.rm = TRUE),
              sal_coeff = sd(Salnty, na.rm = TRUE)/mean(Salnty, na.rm = TRUE),
              PO4_coeff =  sd(PO4ug, na.rm = TRUE)/mean(PO4ug, na.rm = TRUE),
              NO3_coeff = sd(NO3ug, na.rm = TRUE)/mean(NO3ug, na.rm = TRUE),
              SiO3_coeff = sd(SiO3ug, na.rm = TRUE)/mean(SiO3ug, na.rm = TRUE),
              MLD_coeff = sd(MLD_Sigma, na.rm = TRUE)/mean(MLD_Sigma, na.rm = TRUE),
              NC_coeff = sd(NCDepth, na.rm = TRUE)/mean(NCDepth, na.rm = TRUE),
              Chl_coeff = sd(ChlorA, na.rm = TRUE)/mean(ChlorA, na.rm = TRUE),
              evenness = mean(evenness, na.rm = TRUE), shannon = mean(shannon, na.rm = TRUE),
              richness = mean(richness, na.rm = TRUE), chao1 = mean(chao1, na.rm = TRUE),
              gini = mean(gini, na.rm = TRUE), inv_simp = mean(inv_simp, na.rm = TRUE)
    )
  
  
  full_dat$Date <- as.Date(full_dat$Date, format = "%m/%d/%Y")
  
  save(full_dat, file = full_data_file)
  
  save(som_maps, file = map_file)
  
  # look at stations with more than X samples
  
  # x <- 2
  
  # som_maps <- som_maps[which(som_maps$n_samps > 2),]
  
  som_plots <- gather(som_maps, cluster, freq, 2:3)

  save(som_plots, file = regression_file)

  }






