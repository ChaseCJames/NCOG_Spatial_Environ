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



result_tables <- function(input = "data/16s_bacteria.Rdata", som_file = "output/bacteria_16s_som.Rdata",
                          SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
                          physical_data = "data/NCOG_sample_metadata.csv",
                          map_file = "output/bacteria_16s_map.Rdata", regression_file = "output/bacteria_16s_glm.Rdata",
                          dissimmilar_matrix = "output/bacteria_16s_dissimilar.Rdata", 
                          full_data_file = "output/bacteria_16s_full_data.Rdata",
                          sample_regime = "both"){
  
  load(SST_data)
  sst_coeff_table <- coeff_table
  sst_mean_table <- mean_table
  
  load(SLA_data)
  sla_coeff_table <- coeff_table
  sla_mean_table <- mean_table
  
  load(input)
  
  asv_copy <- asv_table
  
  physical_dat <- read.csv(physical_data, header = TRUE, stringsAsFactors = FALSE)
  
  # check variables
  
  data <- physical_dat[colnames(physical_dat)[c(1,33:41,46:47)]]
  
  # data <- data[which(complete.cases(data)==TRUE),]
  
  data$spice <- swSpice(data$Salnty, data$T_degC)
  
  data$eco_name <- paste0("X20",substr(data$Sample.Name,start = 3,stop = 100))
  
  full_dat <- physical_dat[!is.na(match(physical_dat$Sample.Name, data$Sample.Name)),]
  
  full_dat$eco_name <- paste0("X20",substr(full_dat$Sample.Name,start = 3,stop = 100))
  
  full_dat$spice <- swSpice(full_dat$Salnty, full_dat$T_degC)
  
  # physical data match to biological data
  
  full_dat <- full_dat[which(!is.na(match(full_dat$eco_name, rownames(scaled_inputs)))),]
  
  data <- data[which(!is.na(match(data$eco_name, rownames(scaled_inputs)))),]
  
  # distance to coast
  
  coast_calc <- vector()
  
  map <- map_data("world")  
  
  for (i in 1:nrow(full_dat)) {
    
    long_lat <-  c(full_dat$Lon_Dec[i], full_dat$Lat_Dec[i])
    
    if (length(which(!is.na(long_lat))) == 2) {
      
      map_dist <- map[,1:2]
      
      distances <- distm(long_lat, map_dist, fun = distGeo)
      
      coast_calc[i] <- min(distances, na.rm = TRUE)
      
    }
    
    else{coast_calc[i] <- NA}
    
  }
  
  full_dat$dist_to_coast <- coast_calc/1000
  
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
  
  eco.som <- trainSOM(x.data = scaled_inputs, dimension = c(5, 5), nb.save = 10, maxit = 2000, 
                      scaling = "none")
  
  save(eco.som, file = som_file)
  
  # summary(eco.som)
  
  eco.clust <- superClass(eco.som, k = 2)
  
  clusters <- eco.clust$cluster
  
  ids <- eco.clust$som$clustering
  
  som_ids <- clusters[ids]
  
  asv_copy$som_id <- som_ids
  
  full_dat$som_id <- asv_copy$som_id[match(full_dat$eco_name, rownames(asv_copy))]
  
  # shannon diversity
  
  asv_copy$shannon_index <- diversity(asv_table, MARGIN = 1, index = "shannon")
  full_dat$shannon <- asv_copy$shannon_index[match(full_dat$eco_name, rownames(asv_copy))]
  
  # evenness
  
  S <- apply(asv_table>0,1,sum)
  asv_copy$evenness <- diversity(asv_table, index="shannon")/log(S)
  full_dat$evenness <- asv_copy$evenness[match(full_dat$eco_name, rownames(asv_copy))]
  
  # richness
  
  richness <- apply(asv_table, 1, function(x) length(which(x != 0)))
  
  asv_copy$richness <- richness
  full_dat$richness <- asv_copy$richness[match(full_dat$eco_name, rownames(asv_copy))]
  
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
              SiO3_mean = mean(SiO3ug, na.rm = TRUE), C14_mean = mean(IntC14, na.rm = TRUE),
              MLD_mean = mean(MLD_Sigma, na.rm = TRUE), NC_mean = mean(NCDepth, na.rm = TRUE),
              temp_coeff = sd(T_degC, na.rm = TRUE)/mean(T_degC, na.rm = TRUE),
              sal_coeff = sd(Salnty, na.rm = TRUE)/mean(Salnty, na.rm = TRUE),
              PO4_coeff =  sd(PO4ug, na.rm = TRUE)/mean(PO4ug, na.rm = TRUE),
              NO3_coeff = sd(NO3ug, na.rm = TRUE)/mean(NO3ug, na.rm = TRUE),
              SiO3_coeff = sd(SiO3ug, na.rm = TRUE)/mean(SiO3ug, na.rm = TRUE),
              C14_coeff = sd(IntC14, na.rm = TRUE)/mean(IntC14, na.rm = TRUE),
              MLD_coeff = sd(MLD_Sigma, na.rm = TRUE)/mean(MLD_Sigma, na.rm = TRUE),
              NC_coeff = sd(NCDepth, na.rm = TRUE)/mean(NCDepth, na.rm = TRUE),
              evenness = mean(evenness, na.rm = TRUE), shannon = mean(shannon, na.rm = TRUE),
              richness = mean(richness, na.rm = TRUE)
    )
  
  # FOR NOW REMOVE NORTHERN TRANSECTS
  
  som_maps <- som_maps[-c(1:12),]
  
  full_dat$Date <- as.Date(full_dat$Date, format = "%m/%d/%Y")
  
  save(full_dat, file = full_data_file)
  
  # regression data sla
  
  coeff_val <- vector()
  mean_val <- vector()
  
  sst_coeff_val <- vector()
  sst_mean_val <- vector()
  
  coast_distance <- vector()
  
  for (i in 1:nrow(som_maps)) {
    
    long_lat <-  c(som_maps$long[i], som_maps$lat[i])
    sla_dat <- sla_coeff_table[,1:2]
    sst_dat <- sst_coeff_table[,1:2]
    map_dat <- map[,1:2]
    
    sla_dist <- distm(long_lat, sla_dat, fun = distGeo)
    sst_dist <- distm(long_lat, sst_dat, fun = distGeo)
    map_dist <- distm(long_lat, map_dat, fun = distGeo)
    
    coeff_val[i] <- sla_coeff_table$coeff_var[which.min(sla_dist)]
    mean_val[i] <- sla_mean_table$Mean[which.min(sla_dist)]
    
    sst_coeff_val[i] <- sst_coeff_table$coeff_var[which.min(sst_dist)]
    sst_mean_val[i] <- sst_mean_table$Mean[which.min(sst_dist)]
    
    # find distance between coast and location
    coast_distance[i] <- min(map_dist, na.rm = TRUE)
    
  }
  
  som_maps$sla_mean <- mean_val
  som_maps$sla_coeff <- coeff_val
  som_maps$sst_mean <- sst_mean_val
  som_maps$sst_coeff <- sst_coeff_val
  
  som_maps$dist_to_coast <- coast_distance/1000
  
  save(som_maps, file = map_file)
  
  # look at stations with more than 7 samples
  
  som_maps <- som_maps[which(som_maps$n_samps > 7),]
  
  som_plots <- gather(som_maps, cluster, freq, 2:3)

  save(som_plots, file = regression_file)

  }


# both surface and deep

# 16s plastids
set.seed(202)
result_tables(input = "data/16s_plastids.Rdata", som_file = "output/plastid_16s_som.Rdata",
              SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              map_file = "output/plastid_16s_map.Rdata", regression_file = "output/plastid_16s_glm.Rdata",
              dissimmilar_matrix = "output/plastid_16s_dissimilar.Rdata", 
              full_data_file = "output/plastid_16s_full_data.Rdata",
              sample_regime = "both")

# 16s cyanos
set.seed(5)
result_tables(input = "data/16s_cyanos.Rdata", som_file = "output/cyano_16s_som.Rdata",
              SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              map_file = "output/cyano_16s_map.Rdata", regression_file = "output/cyano_16s_glm.Rdata",
              dissimmilar_matrix = "output/cyano_16s_dissimilar.Rdata", 
              full_data_file = "output/cyano_16s_full_data.Rdata",
              sample_regime = "both")

# 16s minus eukaryotes and plastids
set.seed(192)
result_tables(input = "data/16s_bacteria_m_euks.Rdata", som_file = "output/bacteria_m_euks_16s_som.Rdata",
              SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              map_file = "output/bacteria_m_euks_16s_map.Rdata", regression_file = "output/bacteria_m_euks_16s_glm.Rdata",
              dissimmilar_matrix = "output/bacteria_m_euks_16s_dissimilar.Rdata", 
              full_data_file = "output/bacteria_m_euks_16s_full_data.Rdata",
              sample_regime = "both")

# 18sv9 Autotrophic Eukaryotes
set.seed(15)
result_tables(input = "data/18s_autotrophic_euks.Rdata", som_file = "output/euks_auto_18sv9_som.Rdata",
              SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              map_file = "output/euks_auto_18sv9_map.Rdata", regression_file = "output/euks_auto_18sv9_glm.Rdata",
              dissimmilar_matrix = "output/euks_auto_18sv9_dissimilar.Rdata", 
              full_data_file = "output/euks_auto_18sv9_full_data.Rdata",
              sample_regime = "both")

# 18sv9 Eukaryotes
set.seed(987)
result_tables(input = "data/18s_heterotrophic_euks.Rdata", som_file = "output/euks_hetero_18sv9_som.Rdata",
              SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              map_file = "output/euks_hetero_18sv9_map.Rdata", regression_file = "output/euks_hetero_18sv9_glm.Rdata",
              dissimmilar_matrix = "output/euks_hetero_18sv9_dissimilar.Rdata", 
              full_data_file = "output/euks_hetero_18sv9_full_data.Rdata",
              sample_regime = "both")

# all 16s minus plastids
# set.seed(928)
# result_tables(input = "data/16s_bacteria.Rdata", som_file = "output/bacteria_16s_som.Rdata",
#               SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
#               physical_data = "data/NCOG_sample_metadata.csv",
#               map_file = "output/bacteria_16s_map.Rdata", regression_file = "output/bacteria_16s_glm.Rdata",
#               dissimmilar_matrix = "output/bacteria_16s_dissimilar.Rdata", 
#               full_data_file = "output/bacteria_16s_full_data.Rdata",
#               sample_regime = "both")


# Surface Only

# 16s plastids
set.seed(202)
result_tables(input = "data/16s_plastids.Rdata", som_file = "output/plastid_16s_surf_som.Rdata",
              SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              map_file = "output/plastid_16s_surf_map.Rdata", regression_file = "output/plastid_16s_surf_glm.Rdata",
              dissimmilar_matrix = "output/plastid_16s_surf_dissimilar.Rdata", 
              full_data_file = "output/plastid_16s_surf_full_data.Rdata",
              sample_regime = "surface")

# 16s cyanos
set.seed(5)
result_tables(input = "data/16s_cyanos.Rdata", som_file = "output/cyano_16s_surf_som.Rdata",
              SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              map_file = "output/cyano_16s_surf_map.Rdata", regression_file = "output/cyano_16s_surf_glm.Rdata",
              dissimmilar_matrix = "output/cyano_16s_surf_dissimilar.Rdata", 
              full_data_file = "output/cyano_16s_surf_full_data.Rdata",
              sample_regime = "surface")

# 16s minus eukaryotes and plastids
set.seed(192)
result_tables(input = "data/16s_bacteria_m_euks.Rdata", som_file = "output/bacteria_m_euks_16s_surf_som.Rdata",
              SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              map_file = "output/bacteria_m_euks_16s_surf_map.Rdata", regression_file = "output/bacteria_m_euks_16s_surf_glm.Rdata",
              dissimmilar_matrix = "output/bacteria_m_euks_16s_surf_dissimilar.Rdata", 
              full_data_file = "output/bacteria_m_euks_16s_surf_full_data.Rdata",
              sample_regime = "surface")

# 18sv9 Autotrophic Eukaryotes
set.seed(15)
result_tables(input = "data/18s_autotrophic_euks.Rdata", som_file = "output/euks_auto_18sv9_surf_som.Rdata",
              SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              map_file = "output/euks_auto_18sv9_surf_map.Rdata", regression_file = "output/euks_auto_18sv9_surf_glm.Rdata",
              dissimmilar_matrix = "output/euks_auto_18sv9_surf_dissimilar.Rdata", 
              full_data_file = "output/euks_auto_18sv9_surf_full_data.Rdata",
              sample_regime = "surface")

# 18sv9 Eukaryotes
set.seed(987)
result_tables(input = "data/18s_heterotrophic_euks.Rdata", som_file = "output/euks_hetero_18sv9_surf_som.Rdata",
              SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              map_file = "output/euks_hetero_18sv9_surf_map.Rdata", regression_file = "output/euks_hetero_18sv9_surf_glm.Rdata",
              dissimmilar_matrix = "output/euks_hetero_18sv9_surf_dissimilar.Rdata", 
              full_data_file = "output/euks_hetero_18sv9_surf_full_data.Rdata",
              sample_regime = "surface")

# Depth Only

# 16s plastids
set.seed(202)
result_tables(input = "data/16s_plastids.Rdata", som_file = "output/plastid_16s_deep_som.Rdata",
              SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              map_file = "output/plastid_16s_deep_map.Rdata", regression_file = "output/plastid_16s_deep_glm.Rdata",
              dissimmilar_matrix = "output/plastid_16s_deep_dissimilar.Rdata", 
              full_data_file = "output/plastid_16s_deep_full_data.Rdata",
              sample_regime = "depth")

# 16s cyanos
set.seed(5)
result_tables(input = "data/16s_cyanos.Rdata", som_file = "output/cyano_16s_deep_som.Rdata",
              SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              map_file = "output/cyano_16s_deep_map.Rdata", regression_file = "output/cyano_16s_deep_glm.Rdata",
              dissimmilar_matrix = "output/cyano_16s_deep_dissimilar.Rdata", 
              full_data_file = "output/cyano_16s_deep_full_data.Rdata",
              sample_regime = "depth")

# 16s minus eukaryotes and plastids
set.seed(192)
result_tables(input = "data/16s_bacteria_m_euks.Rdata", som_file = "output/bacteria_m_euks_16s_deep_som.Rdata",
              SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              map_file = "output/bacteria_m_euks_16s_deep_map.Rdata", regression_file = "output/bacteria_m_euks_16s_deep_glm.Rdata",
              dissimmilar_matrix = "output/bacteria_m_euks_16s_deep_dissimilar.Rdata", 
              full_data_file = "output/bacteria_m_euks_16s_deep_full_data.Rdata",
              sample_regime = "depth")

# 18sv9 Autotrophic Eukaryotes
set.seed(15)
result_tables(input = "data/18s_autotrophic_euks.Rdata", som_file = "output/euks_auto_18sv9_deep_som.Rdata",
              SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              map_file = "output/euks_auto_18sv9_deep_map.Rdata", regression_file = "output/euks_auto_18sv9_deep_glm.Rdata",
              dissimmilar_matrix = "output/euks_auto_18sv9_deep_dissimilar.Rdata", 
              full_data_file = "output/euks_auto_18sv9_deep_full_data.Rdata",
              sample_regime = "depth")

# 18sv9 Eukaryotes
set.seed(987)
result_tables(input = "data/18s_heterotrophic_euks.Rdata", som_file = "output/euks_hetero_18sv9_deep_som.Rdata",
              SST_data = "output/CALCOFI_temp_tables.Rdata", SLA_data = "output/CALCOFI_sla_tables.Rdata",
              physical_data = "data/NCOG_sample_metadata.csv",
              map_file = "output/euks_hetero_18sv9_deep_map.Rdata", regression_file = "output/euks_hetero_18sv9_deep_glm.Rdata",
              dissimmilar_matrix = "output/euks_hetero_18sv9_deep_dissimilar.Rdata", 
              full_data_file = "output/euks_hetero_18sv9_deep_full_data.Rdata",
              sample_regime = "depth")

