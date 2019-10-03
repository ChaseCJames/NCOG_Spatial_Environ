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


########## Data Processing ############

load("output/CALCOFI_temp_tables.Rdata")

physical_dat <- read.csv("data/NCOG_sample_metadata.csv", header = TRUE, stringsAsFactors = FALSE)

excludes <- read.delim("data/exclude_samples.txt", header = FALSE)

# check variables

data <- physical_dat[colnames(physical_dat)[c(1,33:41)]]

# data <- data[which(complete.cases(data)==TRUE),]

data$spice <- swSpice(data$Salnty, data$T_degC)

data$eco_name <- paste0("X20",substr(data$Sample.Name,start = 3,stop = 100))

full_dat <- physical_dat[!is.na(match(physical_dat$Sample.Name, data$Sample.Name)),]

full_dat$eco_name <- paste0("X20",substr(full_dat$Sample.Name,start = 3,stop = 100))

# read in 18s
eighteen_s <- read.csv("data/18-asv_count_tax_final_update.csv", stringsAsFactors = FALSE)

eight_id_names <- eighteen_s$ï..Feature.ID

eight_tax_id <- eighteen_s[,c(1,ncol(eighteen_s)-1, ncol(eighteen_s))]

eighteen_s <- eighteen_s[,-c(1,ncol(eighteen_s)-1, ncol(eighteen_s))]

eighteen_s <- apply(eighteen_s, 2, as.numeric)

eighteen_s <- t(eighteen_s)

eighteen_s <- as.data.frame(eighteen_s)

eight_tp <- rownames(eighteen_s)

colnames(eighteen_s) <- eight_id_names

# remove mocks from data and rows with no reads
eighteen_s <- eighteen_s[-c(479:491,669:680,873:878),]

# remove bad samples
eighteen_s <- eighteen_s[-which(!is.na(match(rownames(eighteen_s), paste0("X",excludes$V1)))),]

r_sums <- rowSums(eighteen_s, na.rm = TRUE)

eighteen_s <- eighteen_s[-which(r_sums == 0),]

r_sums <- r_sums[-which(r_sums == 0)]

eighteen_s <- as.matrix(eighteen_s)

for (i in 1:nrow(eighteen_s)){
  
  eighteen_s[i,] <- eighteen_s[i,]/r_sums[i]
  
}

eighteen_s <- as.data.frame(eighteen_s)

eight_copy <- eighteen_s

scaled_inputs <- eighteen_s

scaled_inputs <- as.matrix(scaled_inputs)

# physical data match to biological data

full_dat <- full_dat[which(!is.na(match(full_dat$eco_name, rownames(eighteen_s)))),]

data <- data[which(!is.na(match(data$eco_name, rownames(eighteen_s)))),]


############## Run 18sv9 SOM ##############

set.seed(12)

eco.som <- trainSOM(x.data = scaled_inputs, dimension = c(5, 5), nb.save = 10, maxit = 2000, 
                    scaling = "none")

# summary(eco.som)

eco.clust <- superClass(eco.som, k = 2)

clusters <- eco.clust$cluster

ids <- eco.clust$som$clustering

som_ids <- clusters[ids]

eighteen_s$som_id <- som_ids

full_dat$som_id <- eighteen_s$som_id[match(full_dat$eco_name, rownames(eighteen_s))]

# diversity evenness

eighteen_s$shannon_index <- diversity(eight_copy, MARGIN = 1, index = "shannon")
full_dat$shannon <- eighteen_s$shannon_index[match(full_dat$eco_name, rownames(eighteen_s))]

S <- apply(eight_copy>0,1,sum)
eighteen_s$evenness <- diversity(eight_copy, index="shannon")/log(S)
full_dat$evenness <- eighteen_s$evenness[match(full_dat$eco_name, rownames(eighteen_s))]

richness <- apply(eight_copy, 1, function(x) length(which(x != 0)))

eighteen_s$richness <- richness
full_dat$richness <- eighteen_s$richness[match(full_dat$eco_name, rownames(eighteen_s))]

# Dissimilarity

dissimilar <- vegdist(eight_copy, method = "jaccard", binary = TRUE)

dissimilar <- as.matrix(dissimilar)


############ 18sv9 SOM Plots #############

som_maps <- full_dat %>% 
  group_by(Sta_ID) %>%
  summarise(som_1 = sum(som_id == 1, na.rm = TRUE)/n(), som_2 = sum(som_id == 2, na.rm = TRUE)/n(),
             n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE),
            PO4_mean = mean(PO4ug, na.rm = TRUE), NO3_mean = mean(NO3ug, na.rm = TRUE),
            SiO3_mean = mean(SiO3ug, na.rm = TRUE),
            PO4_coeff =  sd(PO4ug, na.rm = TRUE)/mean(PO4ug, na.rm = TRUE),
            NO3_coeff = sd(NO3ug, na.rm = TRUE)/mean(NO3ug, na.rm = TRUE),
            SiO3_coeff = sd(SiO3ug, na.rm = TRUE)/mean(SiO3ug, na.rm = TRUE),
            evenness = mean(evenness, na.rm = TRUE), shannon = mean(shannon, na.rm = TRUE),
            richness = mean(richness, na.rm = TRUE)
            )


som_list <- list()

names_vect <- c("Offshore Cluster", "Nearshore Cluster")
color_vect <- c("darkblue", "darkred")
centroid_vect <- c("red","blue")

# FOR NOW REMOVE NORTHERN TRANSECTS

som_maps <- som_maps[-c(1:12),]

map <- map_data("world")    

# find centroids
centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)

for (i in 1:2) {

p <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",i)), color = "black", size =5, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = color_vect[i], limits = c(0,1)) +
    ggtitle(paste0("Rel. % ",names_vect[i]," per Station")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          title = element_text(hjust = 0.5))

centroid <- wt.centroid(x =centroid_df , p = i+1)

p <- p + geom_point(aes(x = centroid@coords[1], y = centroid@coords[2]), color = centroid_vect[i], size = 5, pch = 10)

print(p)

som_list[[i]] <- p

}

# shannon
shannon <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "shannon"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient(low = "white", high = "blue", aes(min = min(shannon, na.rm = TRUE, max = max(shannon, na.rm = TRUE)), limits = c(min,max))) +
  ggtitle(paste0("Shannon Diversity")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        title = element_text(hjust = 0.5))

evenness <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "evenness"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient(low = "white", high = "red", aes(min = min(evenness, na.rm = TRUE, max = max(evenness, na.rm = TRUE)), limits = c(min,max))) +
  ggtitle(paste0("Evenness")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        title = element_text(hjust = 0.5))

richness <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "richness"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient(low = "white", high = "darkgreen", aes(min = min(richness, na.rm = TRUE, max = max(richness, na.rm = TRUE)), limits = c(min,max))) +
  ggtitle(paste0("Species Richness")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        title = element_text(hjust = 0.5))

# distance to coast

coast_calc <- vector()

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

############## HAB ID ############

taxas <- eight_tax_id[which(!is.na(match(eight_tax_id$ï..Feature.ID, colnames(eighteen_s)))),]

split_taxa <- separate(taxas, Taxon, sep = ";", into = c("A","B","C", "D", "E", "F", "G", "H", "I"))

##### Pseudo-nitschia #####

pseudo_nitz <- which(split_taxa$G == "Pseudo-nitzschia")

pseudo_ids <- split_taxa$id[pseudo_nitz]

neuron_vals <- eco.som$prototypes

pseudo_neuron_vals <- as.data.frame(neuron_vals[,which(!is.na(match(colnames(neuron_vals), pseudo_ids)))])

colnames(pseudo_neuron_vals) <- make.names(split_taxa$H[match(colnames(pseudo_neuron_vals),
                                                              split_taxa$id)], unique = TRUE)

pseudo_neuron_vals$clusters <- clusters

pseudo_nitz_cluster_vals <- pseudo_neuron_vals %>%
  group_by(clusters) %>%
  summarise_all(mean) %>% 
  gather(id, "mean_val",2:37)

pseudo_nitz_cluster_vals <- pseudo_nitz_cluster_vals %>%
  group_by(id) %>%
  mutate(prop=(mean_val/sum(mean_val))) %>%
  ungroup()

pseudo_nitz_cluster_vals$clusters[which(pseudo_nitz_cluster_vals$clusters == 1)] = "Offshore"
pseudo_nitz_cluster_vals$clusters[which(pseudo_nitz_cluster_vals$clusters == 2)] = "Nearshore"

pn1 <- ggplot(pseudo_nitz_cluster_vals, aes(x = id, y = mean_val, fill = clusters)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -90),
        title = element_text(hjust = 0.5)) + 
  ylab("Mean Rel. Abun per Cluster") +
  xlab("Species") + ggtitle("Pseudo-nitzschia")

pn2 <- ggplot(pseudo_nitz_cluster_vals, aes(x = id, y = prop, fill = clusters)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -45),
        title = element_text(hjust = 0.5)) + 
  ylab("Rel. Appearance in each Cluster") +
  xlab("Species") + ggtitle("Pseudo-nitzschia") + geom_hline(yintercept = 0.5, linetype = "dashed")

###### Alexandrium #####

alex <- which(split_taxa$H == "Alexandrium_01" | split_taxa$H == "Alexandrium_02")

alex_ids <- split_taxa$id[alex]

alex_neuron_vals <- as.data.frame(neuron_vals[,which(!is.na(match(colnames(neuron_vals), alex_ids)))])

colnames(alex_neuron_vals) <- make.names(split_taxa$I[match(colnames(alex_neuron_vals),
                                                            split_taxa$id)], unique = TRUE)

alex_neuron_vals$clusters <- clusters

alex_cluster_vals <- alex_neuron_vals %>%
  group_by(clusters) %>%
  summarise_all(mean) %>% 
  gather(id, "mean_val",2:34)

alex_cluster_vals <- alex_cluster_vals %>%
  group_by(id) %>%
  mutate(prop=(mean_val/sum(mean_val))) %>%
  ungroup()

alex_cluster_vals$clusters[which(alex_cluster_vals$clusters == 1)] = "Offshore"
alex_cluster_vals$clusters[which(alex_cluster_vals$clusters == 2)] = "Nearshore"

a1 <- ggplot(alex_cluster_vals, aes(x = id, y = mean_val, fill = clusters)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -90),
        title = element_text(hjust = 0.5)) + 
  ylab("Mean Rel. Abun per Cluster") +
  xlab("Species") + ggtitle("Alexandrium")

a2 <- ggplot(alex_cluster_vals, aes(x = id, y = prop, fill = clusters)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -45),
        title = element_text(hjust = 0.5)) + 
  ylab("Rel. Appearance in each Cluster") +
  xlab("Species") + ggtitle("Alexandrium") + geom_hline(yintercept = 0.5, linetype = "dashed")

###### Dinophysis #####

dino <- which(split_taxa$H == "Dinophysis_03" | split_taxa$H == "Dinophysis_04")

dino_ids <- split_taxa$id[dino]

dino_neuron_vals <- as.data.frame(neuron_vals[,which(!is.na(match(colnames(neuron_vals), dino_ids)))])

colnames(dino_neuron_vals) <- make.names(split_taxa$I[match(colnames(dino_neuron_vals),
                                                            split_taxa$id)], unique = TRUE)

dino_neuron_vals$clusters <- clusters

dino_cluster_vals <- dino_neuron_vals %>%
  group_by(clusters) %>%
  summarise_all(mean) %>% 
  gather(id, "mean_val",2:12)

dino_cluster_vals <- dino_cluster_vals %>%
  group_by(id) %>%
  mutate(prop=(mean_val/sum(mean_val))) %>%
  ungroup()

dino_cluster_vals$clusters[which(dino_cluster_vals$clusters == 1)] = "Offshore"
dino_cluster_vals$clusters[which(dino_cluster_vals$clusters == 2)] = "Nearshore"

d1 <- ggplot(dino_cluster_vals, aes(x = id, y = mean_val, fill = clusters)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -90),
        title = element_text(hjust = 0.5)) + 
  ylab("Mean Rel. Abun per Cluster") +
  xlab("Species") + ggtitle("Dinophysis")

d2 <- ggplot(dino_cluster_vals, aes(x = id, y = prop, fill = clusters)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -45),
        title = element_text(hjust = 0.5)) + 
  ylab("Rel. Appearance in each Cluster") +
  xlab("Species") + ggtitle("Dinophysis") + geom_hline(yintercept = 0.5, linetype = "dashed")

###### Lingulodinium #####

ling <- which(split_taxa$H == "Lingulodinium")

ling_ids <- split_taxa$id[ling]

ling_neuron_vals <- as.data.frame(neuron_vals[,which(!is.na(match(colnames(neuron_vals), ling_ids)))])

colnames(ling_neuron_vals) <- make.names(split_taxa$I[match(colnames(ling_neuron_vals),
                                                            split_taxa$id)], unique = TRUE)

ling_neuron_vals$clusters <- clusters

ling_cluster_vals <- ling_neuron_vals %>%
  group_by(clusters) %>%
  summarise_all(mean) %>% 
  gather(id, "mean_val",2:11)

ling_cluster_vals <- ling_cluster_vals %>%
  group_by(id) %>%
  mutate(prop=(mean_val/sum(mean_val))) %>%
  ungroup()

ling_cluster_vals$clusters[which(ling_cluster_vals$clusters == 1)] = "Offshore"
ling_cluster_vals$clusters[which(ling_cluster_vals$clusters == 2)] = "Nearshore"

l1 <- ggplot(ling_cluster_vals, aes(x = id, y = mean_val, fill = clusters)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -90),
        title = element_text(hjust = 0.5)) + 
  ylab("Mean Rel. Abun per Cluster") +
  xlab("Species") + ggtitle("Lingulodinium polyedrum")

l2 <- ggplot(ling_cluster_vals, aes(x = id, y = prop, fill = clusters)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = -45),
        title = element_text(hjust = 0.5)) + 
  ylab("Rel. Appearance in each Cluster") +
  xlab("Species") + ggtitle("Lingulodinium polyedrum") + geom_hline(yintercept = 0.5, linetype = "dashed")

##### HAB plots ########

pdf("figures/pseudo_nitzschia_clusters.pdf", width = 8, height = 4)
plot_grid(pn2, nrow = 1)
dev.off()

pdf("figures/alexandrium_clusters.pdf", width = 8, height = 4)
plot_grid(a2, nrow = 1)
dev.off()

pdf("figures/dinophysis_clusters.pdf", width = 6, height = 4)
plot_grid(d2, nrow = 1)
dev.off()

########## NetCDF OI SST Data ################

# source("coeff_var_temp.R")

# coeff_center <- SpatialPointsDataFrame(coords = coeff_table[,1:2], data = coeff_table)

sst <- ggplot() + 
  geom_tile(data = coeff_table, aes(x = lon, y = lat, fill = coeff_var), width =0.26, height = 0.26) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", limits = c(0.09,0.12), oob = squish, midpoint = 0.1066851) +geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("SST Coeff. Var.") +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
        title = element_text(hjust = 0.5), axis.line = element_blank())

sst_mean <- ggplot() + 
  geom_tile(data = mean_table, aes(x = lon, y = lat, fill = Mean), width =0.26, height = 0.26) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", limits = c(15,18), oob = squish, midpoint = 16.5) +geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("SST Mean (°C)") +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
        title = element_text(hjust = 0.5), axis.line = element_blank())


######## Combined Plots ##############

p1 <-  ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",1)), color = "black", size =6, stroke = 0.1, shape = 21) +
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0,1)) +
  ggtitle(paste0("C. Offshore Cluster")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
        title = element_text(hjust = 0.5), axis.line = element_blank())

centroid1 <- wt.centroid(x =centroid_df , p = 2)
p1 <- p1 + geom_point(aes(x = centroid1@coords[1], y = centroid1@coords[2]), color = centroid_vect[1], size = 5, pch = 10)

p2 <-  ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",2)), color = "black", size =6, stroke = 0.1, shape = 21) +
  scale_fill_gradient(low = "white", high = "darkred", limits = c(0,1)) +
  ggtitle(paste0("D. Nearshore Cluster")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
        title = element_text(hjust = 0.5), axis.line = element_blank())

centroid2 <- wt.centroid(x =centroid_df , p = 3)
p2 <- p2 + geom_point(aes(x = centroid2@coords[1], y = centroid2@coords[2]), color = centroid_vect[2], size = 5, pch = 10)

theme_set(theme_cowplot(font_size=8))

pdf("figures/3panel_NCOG_expanded_region_0604.pdf", width = 8, height = 3)
plot_grid(sst,p1,p2, ncol = 3, labels = c("A","B","C"), label_size = 10)
dev.off()

theme_set(theme_cowplot(font_size=10))

pdf("figures/4panel_NCOG_expand_ASV_0619.pdf", width = 8, height = 8)
plot_grid(sst,sst_mean, p1,p2, ncol = 2)
dev.off()


########## Time Series ##################

som_map_ts <- full_dat %>% 
  group_by(Sta_ID,Cruise) %>%
  summarise(som_1 = sum(som_id == 1)/n(), som_2 = sum(som_id == 2)/n(),
            lat = mean(Lat_Dec), long = mean(Lon_Dec))


som_map_ts$Cruise <- paste0(substr(som_map_ts$Cruise,1,4),"/",substr(som_map_ts$Cruise,5,6),"/01")

som_map_ts$Cruise <- as.Date(som_map_ts$Cruise, format = "%Y/%m/%d")

dates <- unique(som_map_ts$Cruise)[order(unique(som_map_ts$Cruise))]

for (i in 1:length(dates)) {
  
  subsamp <- som_map_ts[which(som_map_ts$Cruise==dates[i]),]
  
  p1 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-125.3206, -117.45),ylim= c(29.5,36.08734), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = subsamp, aes_string(x = "long", y = "lat",
               fill = paste0("som_",1)), color = "black",
               size =5, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "darkblue", limits = c(0,1)) +
    ggtitle(paste0("Offshore Cluster ", dates[i])) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          title = element_text(hjust = 0.5), axis.line = element_blank())
  
  p2 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-125.3206, -117.45),ylim= c(29.5,36.08734), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = subsamp, aes_string(x = "long", y = "lat",
                fill = paste0("som_",2)), color = "black",
                size =5, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "darkred", limits = c(0,1)) +
    ggtitle(paste0("Nearshore Cluster ",dates[i])) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          title = element_text(hjust = 0.5), axis.line = element_blank())
  
  a <- plot_grid(p1, p2)
  
  print(a)
  
}

############ Temp Plots (Using Station Data) ##################

# surf_dat <- full_dat[which(full_dat$CC_Depth < 40),]
surf_dat <- full_dat


temp_summary <- surf_dat %>%
  group_by(Sta_ID) %>%
  summarise(mean = mean(Temp, na.rm = TRUE), sd = sd(Temp, na.rm = TRUE), coeff_var = sd(Temp, na.rm = TRUE)/mean(Temp, na.rm = TRUE), lat = mean(Lat_Dec, na.rm = TRUE), long = mean(Lon_Dec, na.rm = TRUE))

# for places with only one sample
temp_summary[is.na(temp_summary)] = 0

map <- map_data("world")    

# sst <- ggplot() + 
#   geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
#   coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,36.08734), 1.3) +
#   xlab("Longitude") + ylab("Latitude") + 
#   geom_point(data = temp_summary, aes(x = long, y = lat, fill = coeff_var), color = "black", size =5, stroke = 0.1, shape = 21) +
#   scale_fill_gradient(low = "white", high = "darkblue") +
#   ggtitle("Coefficient of Variation for Temperature") +
#   theme(legend.title = element_blank(),
#         panel.background = element_blank(),
#         panel.border = element_rect(fill = NA),
#         title = element_text(hjust = 0.5))

# print(sst)

############# Physical Chemical SOM ############

data_bits <- data[,5:11]

phys.som <- trainSOM(x.data = data_bits, dimension = c(5, 5), nb.save = 10, maxit = 2000)

summary(phys.som)

phys.clust <- superClass(phys.som, k = 2)

plot(phys.clust, type = "dendrogram")
plot(phys.clust, type = "mds")

plot(phys.som, what = "obs", type = "hitmap")

plot(phys.som, what = "obs", type = "radar", print.title = TRUE)

plot(phys.som, what = "prototypes", type = "smooth.dist")

plot(phys.som, what = "prototypes", type = "lines")

for (i in 1:ncol(data_bits)) {
  plot(phys.som, what = "obs", type = "color", var = i, main = colnames(data_bits)[i])
}

clusters <- phys.clust$cluster

ids <- phys.clust$som$clustering

som_ids <- clusters[ids]

data$som_id_phys <- som_ids

full_dat$som_id_phys <- data$som_id_phys[match(full_dat$eco_name, data$eco_name)]

########## Physical Chemical SOM Plots ##############


som_maps_phys <- full_dat %>% 
  group_by(Sta_ID) %>%
  summarise(som_1 = sum(som_id_phys == 1)/n(), som_2 = sum(som_id_phys == 2)/n(),
            som_3 = sum(som_id_phys == 3)/n(), som_4 = sum(som_id_phys == 4)/n(),
            som_5 = sum(som_id_phys == 5)/n(), lat = mean(Lat_Dec), long = mean(Lon_Dec))


som_list_phys <- list()

for (i in 1:2) {
  
  p <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-125.3206, -117.45),ylim= c(29.5,36.08734), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",i)), color = "black", size =5, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "darkred", limits = c(0,1)) +
    ggtitle(paste0("Proportional Dominance of SOM Cluster ",i," per Station")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          title = element_text(hjust = 0.5))
  
  print(p)
  
  som_list_phys[[i]] <- p
  
}

####### Regression bw Temp and Clusters ########

# look at stations with more than 5 samples

som_maps <- som_maps[which(som_maps$n_samps > 7),]

# Identify closest grid to station 

# som_maps <- som_maps[-which(is.na(som_maps$lat)),]

# Identify closest grid to station 

coeff_val <- vector()
mean_val <- vector()

coast_distance <- vector()

for (i in 1:nrow(som_maps)) {
  
  long_lat <-  c(som_maps$long[i], som_maps$lat[i])
  sst_dat <- coeff_table[,1:2]  
  map_dat <- map[,1:2]
  
  sst_dist <- distm(long_lat, sst_dat, fun = distGeo)
  map_dist <- distm(long_lat, map_dat, fun = distGeo)
  
  coeff_val[i] <- coeff_table$coeff_var[which.min(sst_dist)]
  mean_val[i] <- mean_table$Mean[which.min(sst_dist)]
  
  # find distance between coast and location
  coast_distance[i] <- min(map_dist, na.rm = TRUE)
  
  
}

som_maps$temp_mean <- mean_val
som_maps$temp_coeff <- coeff_val
som_maps$dist_to_coast <- coast_distance/1000


# mean

summary_fit <- summary(lm(som_maps$som_1~som_maps$temp_mean))

p_val <- summary_fit$coefficients[2,4]
r_sq_mean <- summary_fit$r.squared

# text(x = 19, y = 0.2, labels = paste0("R^2 = ", round(r_sq, digits = 4)))
# text(x = 19, y = 0.1, labels = paste0("p-val = ", round(p_val, digits = 4)))

# coeff

summary_fit <- summary(lm(som_maps$som_1~som_maps$temp_coeff))

p_val <- summary_fit$coefficients[2,4]
r_sq_coeff <- summary_fit$r.squared

# text(x = .155, y = 0.2, labels = paste0("R^2 = ", round(r_sq, digits = 4)))
# text(x = .155, y = 0.1, labels = paste0("p-val = ", round(p_val, digits = 8)))

# ggplot

som_plots <- gather(som_maps, cluster, freq, 2:3)

som_plots$cluster[som_plots$cluster == "som_1"] <- "Offshore"
som_plots$cluster[som_plots$cluster == "som_2"] <- "Nearshore"

mean_plot <- ggplot(som_plots, aes(x = temp_mean, y = freq, color = cluster)) +
  geom_point(pch = 20, size = 4) +
  theme(legend.title = element_blank(),
        panel.border = element_rect(color = "black", linetype = "solid"),
        axis.line = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid")) +
  stat_smooth(method = "lm", level = 0.95) +
  scale_color_manual(values = c("red", "blue")) +
  xlab("Mean SST (°C)") + ylab("Frequency") +
  annotate("text", x = 16.2, y = 0.5, label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 4) + ggtitle("F. SOM Freq ~ Mean SST")

coeff_plot <- ggplot(som_plots, aes(x = temp_coeff, y = freq, color = cluster)) +
  geom_point(pch = 20, size = 4) +
  theme(legend.title = element_blank(),
        panel.border = element_rect(color = "black", linetype = "solid"),
        axis.line = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid"),
        legend.position = "none") +
  stat_smooth(method = "lm", level = 0.95) +
  scale_color_manual(values = c("red", "blue")) +
  xlab("SST Coeff. Var") + ylab("Frequency") +
  annotate("text", x = 0.1, y = 0.5, label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 4) + ggtitle("E. SOM Freq ~ Coeff. Var. SST")

temp_som_plots <- plot_grid(coeff_plot, mean_plot, ncol = 2,rel_widths = c(0.8,1))

theme_set(theme_cowplot(font_size=8))

pdf("figures/lm_mean_coeff_asv_multiple_samps_0606.pdf", width = 8, height = 3)
plot_grid(coeff_plot, mean_plot, ncol = 2,rel_widths = c(0.8,1))
dev.off()

# Diversity

summary_fit <- summary(lm(som_maps$shannon~som_maps$temp_mean))

p_val <- summary_fit$coefficients[2,4]
r_sq_mean <- summary_fit$r.squared

summary_fit <- summary(lm(som_maps$shannon~som_maps$temp_coeff))

p_val <- summary_fit$coefficients[2,4]
r_sq_coeff <- summary_fit$r.squared

mean_plot <- ggplot(som_plots, aes(x = temp_mean, y = shannon)) +
  geom_point(pch = 20, size = 4) +
  theme(legend.title = element_blank(),
        panel.border = element_rect(color = "black", linetype = "solid"),
        axis.line = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid")) +
  stat_smooth(method = "lm", level = 0.95) +
  xlab("Mean SST (°C)") + ylab("Shannon Index") +
  annotate("text", x = 16, y = 3.8, label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 4)

coeff_plot <- ggplot(som_plots, aes(x = temp_coeff, y = shannon)) +
  geom_point(pch = 20, size = 4) +
  theme(legend.title = element_blank(),
        panel.border = element_rect(color = "black", linetype = "solid"),
        axis.line = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid"),
        legend.position = "none") +
  stat_smooth(method = "lm", level = 0.95) +
  xlab("SST Coeff. Var") + ylab("Shannon Index") +
  annotate("text", x = 0.1, y = 4, label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 4)

theme_set(theme_cowplot(font_size=8))
temp_diversity_plots <- plot_grid(coeff_plot, mean_plot, ncol = 1,rel_widths = c(1,1))

# Evenness

summary_fit <- summary(lm(som_maps$evenness~som_maps$temp_mean))

p_val <- summary_fit$coefficients[2,4]
r_sq_mean <- summary_fit$r.squared

summary_fit <- summary(lm(som_maps$evenness~som_maps$temp_coeff))

p_val <- summary_fit$coefficients[2,4]
r_sq_coeff <- summary_fit$r.squared

mean_plot <- ggplot(som_plots, aes(x = temp_mean, y = evenness)) +
  geom_point(pch = 20, size = 4) +
  theme(legend.title = element_blank(),
        panel.border = element_rect(color = "black", linetype = "solid"),
        axis.line = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid")) +
  stat_smooth(method = "lm", level = 0.95) +
  xlab("Mean SST (°C)") + ylab("Evenness") +
  annotate("text", x = 15.3, y = 0.75, label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 4)

coeff_plot <- ggplot(som_plots, aes(x = temp_coeff, y = evenness)) +
  geom_point(pch = 20, size = 4) +
  theme(legend.title = element_blank(),
        panel.border = element_rect(color = "black", linetype = "solid"),
        axis.line = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid"),
        legend.position = "none") +
  stat_smooth(method = "lm", level = 0.95) +
  xlab("SST Coeff. Var") + ylab("Evenness") +
  annotate("text", x = 0.1, y = 0.65, label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 4)

theme_set(theme_cowplot(font_size=8))
temp_evenness_plots <- plot_grid(coeff_plot, mean_plot, ncol = 1,rel_widths = c(1,1))

# Richness

summary_fit <- summary(lm(som_maps$richness~som_maps$temp_mean))

p_val <- summary_fit$coefficients[2,4]
r_sq_mean <- summary_fit$r.squared

summary_fit <- summary(lm(som_maps$richness~som_maps$temp_coeff))

p_val <- summary_fit$coefficients[2,4]
r_sq_coeff <- summary_fit$r.squared

mean_plot <- ggplot(som_plots, aes(x = temp_mean, y = richness)) +
  geom_point(pch = 20, size = 4) +
  theme(legend.title = element_blank(),
        panel.border = element_rect(color = "black", linetype = "solid"),
        axis.line = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid")) +
  stat_smooth(method = "lm", level = 0.95) +
  xlab("Mean SST (°C)") + ylab("Evenness") +
  annotate("text", x = 15.3, y = 800, label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 4)

coeff_plot <- ggplot(som_plots, aes(x = temp_coeff, y = richness)) +
  geom_point(pch = 20, size = 4) +
  theme(legend.title = element_blank(),
        panel.border = element_rect(color = "black", linetype = "solid"),
        axis.line = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid"),
        legend.position = "none") +
  stat_smooth(method = "lm", level = 0.95) +
  xlab("SST Coeff. Var") + ylab("Evenness") +
  annotate("text", x = 0.1, y = 800, label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 4)

theme_set(theme_cowplot(font_size=8))
temp_richness_plots <- plot_grid(coeff_plot, mean_plot, ncol = 1,rel_widths = c(1,1))

var_names <- c("PO4", "NO3", "SiO3", "temp")

for (i in 1:4) {
  
  mean_plot <- ggplot(som_plots, aes_string(x = paste0(var_names[i],"_mean"), y = "freq", color = "cluster")) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid")) +
    stat_smooth(method = "lm", level = 0.95) +
    xlab(paste0("Mean ",var_names[i])) + ylab("Frequency") 
  
  coeff_plot <- ggplot(som_plots, aes_string(x = paste0(var_names[i],"_coeff"), y = "freq", color = "cluster")) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none") +
    stat_smooth(method = "lm", level = 0.95) +
    xlab(paste0("Coeff ",var_names[i])) + ylab("Frequency") 
  
  a <- plot_grid(coeff_plot, mean_plot, ncol = 2,rel_widths = c(0.8,1))
  
  print(a)
  
}


summary_fit <- summary(lm(som_maps$som_1~som_maps$PO4_mean))

p_val <- summary_fit$coefficients[2,4]
r_sq_coeff <- summary_fit$r.squared

# trying to fit a glm

som_glm <- som_maps[,c(1:9,16,10:12,17)]

lm_mean_var <- lm(som_1 ~ lat + long + PO4_mean + NO3_mean + SiO3_mean + temp_mean + 
                    PO4_coeff + NO3_coeff + SiO3_coeff + temp_coeff, data = som_glm)

response_df <- som_glm$som_1
predictor_df <- som_glm[,7:14]

regsubsetsObj <- regsubsets(x=predictor_df ,y=response_df, nbest = 2, really.big = T)
regsubsetsObj$xnames <- c("(Intercept)", "Mean PO4", "Mean NO3", "Mean SiO3",
                          "Mean Temp", "Coeff. Var. PO4","Coeff. Var. NO3",
                          "Coeff. Var. SiO3", "Coeff. Var. Temp")

pdf(file = "figures/18sv9_stepwise_regression_leaps_0723.pdf", width = 6, height = 5)
plot(regsubsetsObj, scale = "adjr2", main = "Backwards Subset Selection")
dev.off()

# 
# glm_mean_var <- glm(som_1 ~ lat + long + PO4_mean + NO3_mean + SiO3_mean + temp_mean + 
#       PO4_coeff + NO3_coeff + SiO3_coeff + temp_coeff, data = som_maps)
# 
# stepAIC(glm_mean_var)

model_AIC <- as.data.frame(matrix(NA,10,2))
colnames(model_AIC) <- c("Model","AIC")

# Temp

glm_mean_temp <- glm(som_1 ~  temp_mean, data = som_glm)
mt_sum <- summary(glm_mean_temp)
model_AIC[1,2] <- mt_sum$aic
model_AIC[1,1] <- "Mean Temp"

glm_coeff_temp <- glm(som_1 ~  temp_coeff, data = som_glm)
ct_sum <- summary(glm_coeff_temp)
model_AIC[2,2] <- ct_sum$aic
model_AIC[2,1] <- "Coeff Var. Temp"

# NO3

glm_mean_no3 <- glm(som_1 ~  NO3_mean, data = som_glm)
mn_sum <- summary(glm_mean_no3)
model_AIC[3,2] <- mn_sum$aic
model_AIC[3,1] <- "Mean NO3"

glm_coeff_no3 <- glm(som_1 ~  NO3_coeff, data = som_glm)
cn_sum <- summary(glm_coeff_no3)
model_AIC[4,2] <- cn_sum$aic
model_AIC[4,1] <- "Coeff Var. NO3"

# PO4

glm_mean_po4 <- glm(som_1 ~  PO4_mean, data = som_glm)
mp_sum <- summary(glm_mean_po4)
model_AIC[5,2] <- mp_sum$aic
model_AIC[5,1] <- "Mean PO4"

glm_coeff_po4 <- glm(som_1 ~  PO4_coeff, data = som_glm)
cp_sum <- summary(glm_coeff_po4)
model_AIC[6,2] <- cp_sum$aic
model_AIC[6,1] <- "Coeff Var. PO4"

# SiO3

glm_mean_sio3 <- glm(som_1 ~  SiO3_mean, data = som_glm)
ms_sum <- summary(glm_mean_sio3)
model_AIC[7,2] <- ms_sum$aic
model_AIC[7,1] <- "Mean SiO3"

glm_coeff_sio3 <- glm(som_1 ~  SiO3_coeff, data = som_glm)
cs_sum <- summary(glm_coeff_sio3)
model_AIC[8,2] <- cs_sum$aic
model_AIC[8,1] <- "Coeff Var. SiO3"

# Everything

glm_mean_var <- glm(som_1 ~ PO4_mean + NO3_mean + SiO3_mean + temp_mean + 
                      PO4_coeff + NO3_coeff + SiO3_coeff + temp_coeff, data = som_glm)

all_sum <- summary(glm_mean_var)
model_AIC[9,2] <- all_sum$aic
model_AIC[9,1] <- "Full Model"

# best fit

stepAIC(glm_mean_var)

glm_simple <- glm(som_1 ~ PO4_mean + NO3_mean + temp_mean +
                    PO4_coeff + NO3_coeff + SiO3_coeff + temp_coeff, data = som_glm)

simple_sum <- summary(glm_simple)
model_AIC[10,2] <- simple_sum$aic
model_AIC[10,1] <- "Mean PO4 + Mean NO3 +\n Mean Temp + Coeff Var. PO4 +\n Coeff Var. NO3 + Coeff Var. SiO3 +\n Coeff Var. Temp"

AIC_table <- model_AIC[order(model_AIC$AIC, decreasing = TRUE),]

AIC_table$AIC <- round(AIC_table$AIC, 3)

pdf(file = "figures/18sv9_AIC_table_0723.pdf", height = 5, width = 7)
grid.table(AIC_table, rows = NULL)
dev.off()

####### Random Forests ##########

train <- sample(nrow(som_maps), 0.7*nrow(som_maps), replace = FALSE)


model1 <- randomForest(som_1 ~ PO4_mean + NO3_mean +
                         SiO3_mean + temp_mean + 
                         PO4_coeff + NO3_coeff + SiO3_coeff + temp_coeff,
                       data = som_maps, importance = TRUE)

varImpPlot(model1)

###### Full Data Test ############

full_dat$Date <- as.Date(full_dat$Date, format = "%m/%d/%Y")
full_dat$Julian <- as.numeric(format(full_dat$Date, "%j"))

forest_dat <- full_dat[,-which(is.na(rowSums(full_dat[,c(32,33,37,38,39,59)])))]

forest_dat$som_id <- as.factor(forest_dat$som_id)

train <- sample(nrow(forest_dat), 0.7*nrow(forest_dat), replace = FALSE)

forest_out <- randomForest(som_id ~  CC_Depth + T_degC + PO4ug + SiO3ug +
                             NO3ug, data = forest_dat, importance = TRUE,
                           na.action = na.omit, subset = train, mtry = 4)


# Predicting on validation set
predTrain <- predict(forest_out, forest_dat[train,], type = "class")
# Checking classification accuracy
table(predTrain, forest_dat$som_id[train]) 

# Predicting on validation set
predValid <- predict(forest_out, forest_dat[-train,], type = "class")
# Checking classification accuracy
table(predValid, forest_dat$som_id[-train])  

varImpPlot(forest_out, type = 2, main = paste0("18sv9 Random Forest"))

pdf(file = "figures/18sv9_random_forest_importance_0723.pdf", width = 5, height = 4)
varImpPlot(forest_out, type = 2, main = paste0("18sv9 Random Forest"))
dev.off()


### Other Shannon Evenness plots

shannon_depth_dist <- ggplot(full_dat, aes(x = dist_to_coast, y = CC_Depth, fill = shannon)) + 
  geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
  scale_fill_viridis(name = "Shannon Diversity", option = "D", limits = c(3.0,5.2), oob = scales::squish) +
  ylab("Depth (m)") + xlab("Distance to Coast (km)") 

even_depth_dist <- ggplot(full_dat, aes(x = dist_to_coast, y = CC_Depth, fill = evenness)) + 
  geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
  scale_fill_viridis(name = "Evenness", option = "D", limits = c(0.3,0.8), oob = scales::squish) +
  ylab("Depth (m)") + xlab("Distance to Coast (km)") 

rich_depth_dist <- ggplot(full_dat, aes(x = dist_to_coast, y = CC_Depth, fill = richness)) + 
  geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
  scale_fill_viridis(name = "Richness", option = "D", limits = c(300,900), oob = scales::squish) +
  ylab("Depth (m)") + xlab("Distance to Coast (km)") 

full_dat$Date <- as.Date(full_dat$Date, format = "%m/%d/%Y")

shannon_timeplot <- ggplot(full_dat, aes(x = Date, y = dist_to_coast, fill = shannon)) + 
  geom_jitter(size = 4, pch = 21, color = "black",height = 0, width = 10) + 
  scale_fill_viridis(name = "Shannon", option = "D", limits = c(3.0,5.2), oob = scales::squish) +
  ylab("Distance to Coast") + xlab("Time") 

even_timeplot <- ggplot(full_dat, aes(x = Date, y = dist_to_coast, fill = evenness)) + 
  geom_jitter(size = 4, pch = 21, color = "black",height = 0, width = 10) + 
  scale_fill_viridis(name = "Evenness", option = "D", limits = c(0.3,0.8), oob = scales::squish) +
  ylab("Distance to Coast") + xlab("Time") 

rich_timeplot <- ggplot(full_dat, aes(x = Date, y = dist_to_coast, fill = richness)) + 
  geom_jitter(size = 4, pch = 21, color = "black",height = 0, width = 10) + 
  scale_fill_viridis(name = "Richness", option = "D", limits = c(300,900), oob = scales::squish) +
  ylab("Distance to Coast") + xlab("Time") 

title <- ggdraw() + draw_label("All 18sv9", fontface='bold')

pdf(file = "figures/all18sv9_som_summary_0717.pdf", width = 8, height = 10)
plot_grid(title, plot_grid(sst,sst_mean, ncol = 2), 
          plot_grid(p1,p2, ncol = 2),
          temp_som_plots, ncol = 1, nrow = 4, rel_heights = c(0.1,1,1,1))
dev.off()

title1 <- ggdraw() + draw_label("All 18sv9 Diversity/Evenness/Richness", fontface='bold')
pdf(file = "figures/all18sv9_diversity_summary_0717.pdf", width = 12, height = 12)
plot_grid(title1, plot_grid(shannon,evenness,richness, ncol = 3, labels = c("A","B", "C")),
          plot_grid(temp_diversity_plots, temp_evenness_plots, temp_richness_plots, ncol = 3, labels = c("D","E", "F")),
          plot_grid(shannon_depth_dist, even_depth_dist, rich_depth_dist, ncol = 3, labels = c("G","H", "I")),
          ncol = 1, nrow = 4, rel_heights = c(0.1,1,2,1))
dev.off()

title2 <- ggdraw() + draw_label("All 18sv9 Diversity/Evenness/Richness over Time", fontface='bold')
pdf(file = "figures/all18sv9_diversity_time_0717.pdf", width = 8, height = 12)
plot_grid(title2, shannon_timeplot, even_timeplot, rich_timeplot,
          ncol = 1, nrow = 4, rel_heights = c(0.1,1,1,1))
dev.off()

####### Dissimilarity ######

full_dat$Month_Day1 <- format(full_dat$Date, format="%m-%d-2000")
full_dat$Month_Day1 <- as.Date(full_dat$Month_Day1, format = "%m-%d-%Y")

full_dat$Month_Day2 <- format(full_dat$Date, format="%m-%d-1999")
full_dat$Month_Day2 <- as.Date(full_dat$Month_Day2, format = "%m-%d-%Y")


temp_mat <- matrix(NA, nrow = NROW(dissimilar), ncol = NCOL(dissimilar))
chl_mat <- matrix(NA, nrow = NROW(dissimilar), ncol = NCOL(dissimilar))
nitrate_mat <- matrix(NA, nrow = NROW(dissimilar), ncol = NCOL(dissimilar))
depth_mat <- matrix(NA, nrow = NROW(dissimilar), ncol = NCOL(dissimilar))
date_mat <- matrix(NA, nrow = NROW(dissimilar), ncol = NCOL(dissimilar))
season_mat <- matrix(NA, nrow = NROW(dissimilar), ncol = NCOL(dissimilar))
distance_mat <- matrix(NA, nrow = NROW(dissimilar), ncol = NCOL(dissimilar))

# mean

temp_mean <- matrix(NA, nrow = NROW(dissimilar), ncol = NCOL(dissimilar))
chl_mean <- matrix(NA, nrow = NROW(dissimilar), ncol = NCOL(dissimilar))
nitrate_mean <- matrix(NA, nrow = NROW(dissimilar), ncol = NCOL(dissimilar))
depth_mean <- matrix(NA, nrow = NROW(dissimilar), ncol = NCOL(dissimilar))

for (i in 1:nrow(dissimilar)) {
  for (j in 1:ncol(dissimilar)) {
    
    # delta t
    t_1 <- full_dat$T_degC[which(full_dat$eco_name == rownames(dissimilar)[i])]
    t_2 <- full_dat$T_degC[which(full_dat$eco_name == colnames(dissimilar)[j])]
    temp_mat[i,j] <- abs(t_1-t_2)
    temp_mean[i,j] <- mean(c(t_1,t_2), na.rm = TRUE)
    
    
    # delta chl
    chl_1 <- full_dat$ChlorA[which(full_dat$eco_name == rownames(dissimilar)[i])]
    chl_2 <- full_dat$ChlorA[which(full_dat$eco_name == colnames(dissimilar)[j])]
    chl_mat[i,j] <- abs(chl_1-chl_2)
    chl_mean[i,j] <- mean(c(chl_1,chl_2), na.rm = TRUE)
    
    # delta nitrate
    n_1 <- full_dat$NO3ug[which(full_dat$eco_name == rownames(dissimilar)[i])]
    n_2 <- full_dat$NO3ug[which(full_dat$eco_name == colnames(dissimilar)[j])]
    nitrate_mat[i,j] <- abs(n_1-n_2)
    nitrate_mean[i,j] <- mean(c(n_1,n_2), na.rm = TRUE)
    
    # delta depth
    d_1 <- full_dat$CC_Depth[which(full_dat$eco_name == rownames(dissimilar)[i])]
    d_2 <- full_dat$CC_Depth[which(full_dat$eco_name == colnames(dissimilar)[j])]
    depth_mat[i,j] <- abs(d_1-d_2)
    depth_mean[i,j] <- mean(c(d_1,d_2), na.rm = TRUE)
    
    # delta date
    date_1 <- full_dat$Date[which(full_dat$eco_name == rownames(dissimilar)[i])]
    date_2 <- full_dat$Date[which(full_dat$eco_name == colnames(dissimilar)[j])]
    date_mat[i,j] <- abs(as.numeric(date_1-date_2))
    
    # delta season
    diff_vect <- vector()
    
    season_1 <- full_dat$Month_Day1[which(full_dat$eco_name == rownames(dissimilar)[i])]
    season_2 <- full_dat$Month_Day1[which(full_dat$eco_name == colnames(dissimilar)[j])]
    season_3 <- full_dat$Month_Day2[which(full_dat$eco_name == colnames(dissimilar)[j])]
    season_4 <- full_dat$Month_Day2[which(full_dat$eco_name == rownames(dissimilar)[i])]
    diff_vect[1] <- abs(as.numeric(season_1-season_2))
    diff_vect[2] <- abs(as.numeric(season_1-season_3))
    diff_vect[3] <- abs(as.numeric(season_2-season_4))
    
    season_mat[i,j] <- min(diff_vect, na.rm = TRUE)
    
    # delta distance
    
    long_lat_1 <- c(full_dat$Lon_Dec[which(full_dat$eco_name == rownames(dissimilar)[i])],
                    full_dat$Lat_Dec[which(full_dat$eco_name == rownames(dissimilar)[i])])
    
    long_lat_2 <- c(full_dat$Lon_Dec[which(full_dat$eco_name == rownames(dissimilar)[j])],
                    full_dat$Lat_Dec[which(full_dat$eco_name == rownames(dissimilar)[j])])
    
    distance <- distm(long_lat_1, long_lat_2, fun = distGeo)/1000
    distance_mat[i,j] <- distance
    
  }
}

diag(dissimilar) <- NA
dissimilar[upper.tri(dissimilar)] <- NA

dissimilar_plots <- as.data.frame(matrix(NA,length(as.vector(dissimilar)),12))
colnames(dissimilar_plots) <- c("Dissimilarity", "Delta_T", "Delta_Chl", "Delta_N03",
                                "Delta_Depth", "Delta_Date", "Delta_Season", "Temp_Mean",
                                "Chl_Mean", "N03_Mean", "Depth_Mean", "Distance")

dissimilar_plots$Dissimilarity <- as.vector(dissimilar)
dissimilar_plots$Delta_T <- as.vector(temp_mat)
dissimilar_plots$Delta_Chl <- as.vector(chl_mat)
dissimilar_plots$Delta_N03 <- as.vector(nitrate_mat)
dissimilar_plots$Delta_Depth <- as.vector(depth_mat)
dissimilar_plots$Delta_Date <- as.vector(date_mat)
dissimilar_plots$Delta_Season <- as.vector(season_mat)
dissimilar_plots$Temp_Mean <- as.vector(temp_mean)
dissimilar_plots$Chl_Mean <- as.vector(chl_mean)
dissimilar_plots$N03_Mean <- as.vector(nitrate_mean)
dissimilar_plots$Depth_Mean <- as.vector(depth_mean)
dissimilar_plots$Distance <- as.vector(distance_mat)

dissimilar_plots <- dissimilar_plots[-which(is.na(dissimilar_plots$Dissimilarity)),]

dissimilar_plots <- dissimilar_plots[-which(is.na(dissimilar_plots$Delta_T)),]

dissimilar_plots$Similarity <- 1 - dissimilar_plots$Dissimilarity

# temp_dissimilar <- ggplot(dissimilar_plots,
#                           aes(x = Delta_T, y = Dissimilarity, color = Delta_Season)) +
#   geom_point(alpha = 0.3) + xlab("\u0394 Temperature (°C)") + ylab("Jaccard's Dissimilarity") +
#   scale_color_viridis(name = "\u0394 Days (Excluding Year)", option = "B")
# 
# chl_dissimilar <- ggplot(dissimilar_plots,
#                          aes(x = Delta_Chl, y = Dissimilarity, color = Delta_Season)) +
#   geom_point(alpha = 0.3) + xlab("\u0394 Chl-a") + ylab("Jaccard's Dissimilarity") +
#   scale_color_viridis(name = "\u0394 Days (Excluding Year)", option = "B")
# 
# n03_dissimilar <- ggplot(dissimilar_plots,
#                          aes(x = Delta_N03, y = Dissimilarity, color = Delta_Season)) +
#   geom_point(alpha = 0.3) + xlab(expression(Delta~"N03 ("*mu*"g)")) + ylab("Jaccard's Dissimilarity")  +
#   scale_color_viridis(name = "\u0394 Days (Excluding Year)", option = "B")
# 
# depth_dissimilar <- ggplot(dissimilar_plots,
#                            aes(x = Delta_Depth, y = Dissimilarity, color = Delta_Season)) +
#   geom_point(alpha = 0.3) + xlab("\u0394 Depth (m)") + ylab("Jaccard's Dissimilarity") +
#   scale_color_viridis(name = "\u0394 Days (Excluding Year)", option = "B")
# 
# date_dissimilar <- ggplot(dissimilar_plots,
#                           aes(x = Delta_Date, y = Dissimilarity, color = Delta_Season)) +
#   geom_point(alpha = 0.3) + xlab("\u0394 Date (days)") + ylab("Jaccard's Dissimilarity") +
#   scale_color_viridis(name = "\u0394 Days (Excluding Year)", option = "B")
# 
# season_dissimilar <- ggplot(dissimilar_plots,
#                             aes(x = Delta_Season, y = Dissimilarity, color = Delta_Season)) +
#   geom_point(alpha = 0.3) + xlab("\u0394 Days (Excluding Year)") + ylab("Jaccard's Dissimilarity") +
#   scale_color_viridis(name = "\u0394 Days (Excluding Year)", option = "B")

# hexbins

temp_hexbin <- hexbin(dissimilar_plots$Delta_T, dissimilar_plots$Similarity, IDs = TRUE, xbins = 15)
plot(temp_hexbin)
temp_hexbin@cAtt <- as.vector(hexTapply(temp_hexbin, dissimilar_plots$Delta_Season, mean))
hexbin_df <- data.frame(hcell2xy(temp_hexbin),  hexID = temp_hexbin@cell, counts = temp_hexbin@cAtt)
temp_hex <- ggplot(hexbin_df, aes(x = x, y = y, fill = counts, hexID = hexID)) + geom_hex(stat= "identity") +
  scale_fill_viridis(name = expression(Delta~"Days (Excluding Year)"), option = "B") + xlab(expression(Delta~"Temperature (°C)")) +
  ylab("Jaccard's Similarity")

chl_hexbin <- hexbin(dissimilar_plots$Delta_Chl, dissimilar_plots$Similarity, IDs = TRUE, xbins = 15)
plot(chl_hexbin)
chl_hexbin@cAtt <- as.vector(hexTapply(chl_hexbin, dissimilar_plots$Delta_Season, mean))
hexbin_df <- data.frame(hcell2xy(chl_hexbin),  hexID = chl_hexbin@cell, counts = chl_hexbin@cAtt)
chl_hex <- ggplot(hexbin_df, aes(x = x, y = y, fill = counts, hexID = hexID)) + geom_hex(stat= "identity") +
  scale_fill_viridis(name = expression(Delta~"Days (Excluding Year)"), option = "B") + xlab(expression(Delta~"Chl-a")) +
  ylab("Jaccard's Similarity")

n03_hexbin <- hexbin(dissimilar_plots$Delta_N03, dissimilar_plots$Similarity, IDs = TRUE, xbins = 15)
plot(n03_hexbin)
n03_hexbin@cAtt <- as.vector(hexTapply(n03_hexbin, dissimilar_plots$Delta_Season, mean))
hexbin_df <- data.frame(hcell2xy(n03_hexbin),  hexID = n03_hexbin@cell, counts = n03_hexbin@cAtt)
n03_hex <- ggplot(hexbin_df, aes(x = x, y = y, fill = counts, hexID = hexID)) + geom_hex(stat= "identity") +
  scale_fill_viridis(name = expression(Delta~"Days (Excluding Year)"), option = "B") + xlab(expression(Delta~"N03 ("*mu*"g)")) +
  ylab("Jaccard's Similarity")

depth_hexbin <- hexbin(dissimilar_plots$Delta_Depth, dissimilar_plots$Similarity, IDs = TRUE, xbins = 15)
plot(depth_hexbin)
depth_hexbin@cAtt <- as.vector(hexTapply(depth_hexbin, dissimilar_plots$Delta_Season, mean))
hexbin_df <- data.frame(hcell2xy(depth_hexbin),  hexID = depth_hexbin@cell, counts = depth_hexbin@cAtt)
depth_hex <- ggplot(hexbin_df, aes(x = x, y = y, fill = counts, hexID = hexID)) + geom_hex(stat= "identity") +
  scale_fill_viridis(name = expression(Delta~"Days (Excluding Year)"), option = "B") + xlab(expression(Delta~"Depth (m)")) +
  ylab("Jaccard's Similarity")

date_hexbin <- hexbin(dissimilar_plots$Delta_Date, dissimilar_plots$Similarity, IDs = TRUE, xbins = 25)
plot(date_hexbin)
date_hexbin@cAtt <- as.vector(hexTapply(date_hexbin, dissimilar_plots$Delta_Season, mean))
hexbin_df <- data.frame(hcell2xy(date_hexbin),  hexID = date_hexbin@cell, counts = date_hexbin@cAtt)
date_hex <- ggplot(hexbin_df, aes(x = x, y = y, fill = counts, hexID = hexID)) + geom_hex(stat= "identity") +
  scale_fill_viridis(name = expression(Delta~"Days (Excluding Year)"), option = "B") + xlab(expression(Delta~"Days")) +
  ylab("Jaccard's Similarity")

season_hexbin <- hexbin(dissimilar_plots$Delta_Season, dissimilar_plots$Similarity, IDs = TRUE, xbins = 15)
plot(season_hexbin)
season_hexbin@cAtt <- as.vector(hexTapply(season_hexbin, dissimilar_plots$Delta_Season, mean))
hexbin_df <- data.frame(hcell2xy(season_hexbin),  hexID = season_hexbin@cell, counts = season_hexbin@cAtt)
season_hex <- ggplot(hexbin_df, aes(x = x, y = y, fill = counts, hexID = hexID)) + geom_hex(stat= "identity") +
  scale_fill_viridis(name = expression(Delta~"Days (Excluding Year)"), option = "B") + xlab(expression(Delta~"Days (Excluding Year)")) +
  ylab("Jaccard's Similarity")

distance_hexbin <- hexbin(dissimilar_plots$Distance, dissimilar_plots$Similarity, IDs = TRUE, xbins = 25)
plot(distance_hexbin)
# depth_hexbin@cAtt <- as.vector(hexTapply(depth_hexbin, dissimilar_plots$Delta_Season, mean))
hexbin_df <- data.frame(hcell2xy(distance_hexbin),  hexID = distance_hexbin@cell, counts = distance_hexbin@count)
distance_hex <- ggplot(hexbin_df, aes(x = x, y = y, fill = counts, hexID = hexID)) + geom_hex(stat= "identity") + 
  xlab("Distance (km)") + ylab("Jaccard's Similarity")

season_hexbin <- hexbin(dissimilar_plots$Delta_Season, dissimilar_plots$Similarity, IDs = TRUE, xbins = 25)
plot(season_hexbin)
hexbin_df <- data.frame(hcell2xy(season_hexbin),  hexID = season_hexbin@cell, counts = season_hexbin@count)
season_count_hex <- ggplot(hexbin_df, aes(x = x, y = y, fill = counts, hexID = hexID)) + geom_hex(stat= "identity") +
  xlab(expression(Delta~"Days (Excluding Year)")) + ylab("Jaccard's Similarity")

# main plot

title_dissim <- ggdraw() + draw_label("18sv9 Dissimilarity Plots", fontface='bold')

theme_set(theme_cowplot(font_size=8))

pdf(file = "figures/18sv9_dissimilarity_0807.pdf", width = 12, height = 6)
plot_grid(title_dissim,
          plot_grid(distance_hex, season_count_hex,
                    nrow = 1, ncol = 2,
                    labels = c("A", "B")),
          nrow = 2, ncol = 1, rel_heights = c(0.1,1))
dev.off()

# Individual Plots

pdf(file = "figures/18sv9_temp_dissim_0722.pdf", width = 5, height = 5)
temp_dissimilar
dev.off()

pdf(file = "figures/18sv9_chl_dissim_0722.pdf", width = 5, height = 5)
chl_dissimilar
dev.off()

pdf(file = "figures/18sv9_n03_dissim_0722.pdf", width = 5, height = 5)
n03_dissimilar
dev.off()

pdf(file = "figures/18sv9_depth_dissim_0722.pdf", width = 5, height = 5)
depth_dissimilar
dev.off()

pdf(file = "figures/18sv9_date_dissim_0722.pdf", width = 5, height = 5)
date_dissimilar
dev.off()

pdf(file = "figures/18sv9_season_dissim_0722.pdf", width = 5, height = 5)
season_dissimilar
dev.off()







####### Total Diversity/Evenness/Richness #############

# Mean Alpha Diversity

summary_fit <- summary(lm(som_maps$shannon~som_maps$temp_mean))

alpha_mean_p_val <- summary_fit$coefficients[2,4]
alpha_r_sq_mean <- summary_fit$r.squared

summary_fit <- summary(lm(som_maps$shannon~som_maps$temp_coeff))

alpha_coeff_p_val <- summary_fit$coefficients[2,4]
alpha_r_sq_coeff <- summary_fit$r.squared

alpha_mean_plot <- ggplot(som_plots, aes(x = temp_mean, y = shannon)) +
  geom_point(pch = 20, size = 4) +
  theme(legend.title = element_blank(),
        panel.border = element_rect(color = "black", linetype = "solid"),
        axis.line = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid")) +
  stat_smooth(method = "lm", level = 0.95) +
  xlab("Mean SST (°C)") + ylab("Shannon Index") +
  annotate("text", x = 16, y = 3.8, label = paste0("R Squared = ", round(alpha_r_sq_mean, digits = 3)), size = 4)

alpha_coeff_plot <- ggplot(som_plots, aes(x = temp_coeff, y = shannon)) +
  geom_point(pch = 20, size = 4) +
  theme(legend.title = element_blank(),
        panel.border = element_rect(color = "black", linetype = "solid"),
        axis.line = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid"),
        legend.position = "none") +
  stat_smooth(method = "lm", level = 0.95) +
  xlab("SST Coeff. Var") + ylab("Shannon Index") +
  annotate("text", x = 0.1, y = 4, label = paste0("R Squared = ", round(alpha_r_sq_coeff, digits = 3)), size = 4)

# Gamma Diversity

eight_sums <- eight_copy

eight_sums$station <- full_dat$Sta_ID[match(rownames(eight_copy), full_dat$eco_name)]

station_sums <- eight_sums %>% 
  group_by(station) %>% 
  summarise_all(mean, na.rm = TRUE)

station_sums <- as.data.frame(station_sums)

rownames(station_sums) <- station_sums$station

station_sums$station <- NULL

station_plot_df <- as.data.frame(matrix(NA,NROW(station_sums),6))

colnames(station_plot_df) <- c("Station", "Diversity","Evenness","Richness", "Latitude", "Longitude")

station_plot_df$Station <- rownames(station_sums)

# diversity evenness

station_plot_df$Diversity <- diversity(station_sums, MARGIN = 1, index = "shannon")

S <- apply(station_sums>0,1,sum)
station_plot_df$Evenness <- diversity(station_sums, index="shannon")/log(S)

# richness

station_plot_df$Richness <- apply(station_sums, 1, function(x) length(which(x != 0)))

station_position <- full_dat %>% 
  group_by(Sta_ID) %>% 
  summarise(Lat = mean(Lat_Dec, na.rm = TRUE), Lon = mean(Lon_Dec, na.rm = TRUE))

station_plot_df$Latitude <- station_position$Lat[match(station_plot_df$Station, station_position$Sta_ID)]
station_plot_df$Longitude <- station_position$Lon[match(station_plot_df$Station, station_position$Sta_ID)]

# remove northern transects

line <- as.numeric(substr(station_plot_df$Station,2,3))
n_stations <- which(line < 76)
if(length(n_stations) > 0){station_plot_df <- station_plot_df[-n_stations,]}

# shannon
shannon_stat <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = station_plot_df, aes_string(x = "Longitude", y = "Latitude", fill = "Diversity"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient(low = "white", high = "blue", aes(min = min(Diversity, na.rm = TRUE, max = max(Diversity, na.rm = TRUE)), limits = c(min,max))) +
  ggtitle(paste0("Shannon Diversity")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        title = element_text(hjust = 0.5))

evenness_stat <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = station_plot_df, aes_string(x = "Longitude", y = "Latitude", fill = "Evenness"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient(low = "white", high = "red", aes(min = min(Evenness, na.rm = TRUE, max = max(Evenness, na.rm = TRUE)), limits = c(min,max))) +
  ggtitle(paste0("Evenness")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        title = element_text(hjust = 0.5))

richness_stat <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = station_plot_df, aes_string(x = "Longitude", y = "Latitude", fill = "Richness"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient(low = "white", high = "darkgreen", aes(min = min(Richness, na.rm = TRUE, max = max(Richness, na.rm = TRUE)), limits = c(min,max))) +
  ggtitle(paste0("Species Richness")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        title = element_text(hjust = 0.5))


title_stat <- ggdraw() + draw_label("18sv9 Diversity/Evenness/Richness by Station over Time", fontface='bold')

theme_set(theme_cowplot(font_size=8))

pdf(file = "figures/18sv9_station_div_0722.pdf", width = 10, height = 6)
plot_grid(title_stat,
          plot_grid(shannon_stat, evenness_stat, richness_stat, nrow = 1,
                    labels = c("A", "B", "C")),
          nrow = 2, ncol = 1, rel_heights = c(0.1,1))
dev.off()

# Alpha Gamma Plots

som_maps$Gamma_Diversity <- station_plot_df$Diversity[match(som_maps$Sta_ID, station_plot_df$Station)]

summary_fit <- summary(lm(som_maps$Gamma_Diversity~som_maps$temp_coeff))

gamma_coeff_p_val <- summary_fit$coefficients[2,4]
gamma_r_sq_coeff <- summary_fit$r.squared



div_plot <- ggplot(som_maps, aes(x = temp_coeff)) + 
  geom_point(aes(y = shannon, color = "Alpha Diversity")) +
  stat_smooth(aes(y = shannon), method = "lm", level = 0.95, color = "black") +
  geom_point(aes(y = Gamma_Diversity, color = "Gamma Diversity")) +
  stat_smooth(aes(y = Gamma_Diversity), method = "lm", level = 0.95, color = "black") +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Gamma Diversity")) +
  ylab("Alpha Diversity") + xlab("Coeff. Var. SST") +
  theme(legend.title = element_blank(),
        legend.box.background = element_rect(color = "black"),
        panel.background = element_blank()) +
  ggtitle("18sv9 Mean Alpha Diversity vs Gamma Diversity per Station") +
  scale_color_manual(values = c("royalblue2", "seagreen3")) +
  annotate("text", x = 0.097, y = 4, label = paste0("R Squared = ", round(alpha_r_sq_coeff, digits = 3)), size = 4) +
  annotate("text", x = 0.125, y = 6.2, label = paste0("R Squared = ", round(gamma_r_sq_coeff, digits = 3)), size = 4)

theme_set(theme_cowplot(font_size=10))

pdf(file = "figures/18sv9_alpha_v_gamma_0723.pdf", width = 8, height = 6)
div_plot
dev.off()


