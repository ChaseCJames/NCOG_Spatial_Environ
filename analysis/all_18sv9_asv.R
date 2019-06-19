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


########## Data Processing ############

load("output/CALCOFI_temp_tables.Rdata")

physical_dat <- read.csv("data/NCOG_sample_metadata.csv", header = TRUE, stringsAsFactors = FALSE)

# check variables

data <- physical_dat[colnames(physical_dat)[c(1,33:41)]]

# data <- data[which(complete.cases(data)==TRUE),]

data$spice <- swSpice(data$Salnty, data$T_degC)

data$eco_name <- paste0("X20",substr(data$Sample.Name,start = 3,stop = 100))

full_dat <- physical_dat[!is.na(match(physical_dat$Sample.Name, data$Sample.Name)),]

full_dat$eco_name <- paste0("X20",substr(full_dat$Sample.Name,start = 3,stop = 100))

# read in 18s
eighteen_s <- read.csv("data/asv_count_tax_final_update.csv", stringsAsFactors = FALSE)

eight_id_names <- eighteen_s$Feature.ID

eight_tax_id <- eighteen_s[,c(1,ncol(eighteen_s)-1, ncol(eighteen_s))]

eighteen_s <- eighteen_s[,-c(1,ncol(eighteen_s)-1, ncol(eighteen_s))]

eighteen_s <- apply(eighteen_s, 2, as.numeric)

eighteen_s <- t(eighteen_s)

eighteen_s <- as.data.frame(eighteen_s)

eight_tp <- rownames(eighteen_s)

colnames(eighteen_s) <- eight_id_names

# remove mocks from data and rows with no reads
eighteen_s <- eighteen_s[-c(479:491,669:680,873:878),]

r_sums <- rowSums(eighteen_s, na.rm = TRUE)

eighteen_s <- eighteen_s[-which(r_sums == 0),]

r_sums <- r_sums[-which(r_sums == 0)]

eighteen_s <- as.matrix(eighteen_s)

for (i in 1:nrow(eighteen_s)){
  
  eighteen_s[i,] <- eighteen_s[i,]/r_sums[i]
  
}

# cutoff

# remove_vect <- vector()
# 
# for (i in 1:ncol(eighteen_s)) {
# 
#  vals <- which(eighteen_s[,i] > 0.001)
# 
#  if(length(vals) < 1){remove_vect[i] <- 1}else{remove_vect[i] <- NA}
# 
# }
# 
# eighteen_s <- eighteen_s[,-which(!is.na(remove_vect))]

# make dataframes

eighteen_s <- as.data.frame(eighteen_s)

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

# plot(eco.clust, type = "dendrogram")
# plot(eco.clust, type = "mds")

# plot(eco.som, what = "obs", type = "hitmap")

# for (i in 1:length(colnames(scaled_inputs))) {
#   plot(eco.som, what = "obs", type = "color", var = i, print.title = TRUE, main = colnames(scaled_inputs)[i])
# }

# plot(eco.som, what = "obs", type = "radar", print.title = TRUE)

# plot(eco.som, what = "prototypes", type = "smooth.dist")

# plot(eco.clust, type = "lines")

clusters <- eco.clust$cluster

ids <- eco.clust$som$clustering

som_ids <- clusters[ids]

eighteen_s$som_id <- som_ids

full_dat$som_id <- eighteen_s$som_id[match(full_dat$eco_name, rownames(eighteen_s))]

# diversity evenness

eighteen_s$shannon_index <- diversity(eighteen_s, MARGIN = 1, index = "shannon")
full_dat$shannon <- eighteen_s$shannon_index[match(full_dat$eco_name, rownames(eighteen_s))]

S <- apply(eighteen_s>0,1,sum)
eighteen_s$evenness <- diversity(eighteen_s, index="shannon")/log(S)
full_dat$evenness <- eighteen_s$evenness[match(full_dat$eco_name, rownames(eighteen_s))]


############ 18sv9 SOM Plots #############

som_maps <- full_dat %>% 
  group_by(Sta_ID) %>%
  summarise(som_1 = sum(som_id == 1)/n(), som_2 = sum(som_id == 2)/n(),
             n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE),
            PO4_mean = mean(PO4ug, na.rm = TRUE), NO3_mean = mean(NO3ug, na.rm = TRUE),
            SiO3_mean = mean(SiO3ug, na.rm = TRUE),
            PO4_coeff =  sd(PO4ug, na.rm = TRUE)/mean(PO4ug, na.rm = TRUE),
            NO3_coeff = sd(NO3ug, na.rm = TRUE)/mean(NO3ug, na.rm = TRUE),
            SiO3_coeff = sd(SiO3ug, na.rm = TRUE)/mean(SiO3ug, na.rm = TRUE),
            evenness = mean(evenness, na.rm = TRUE), shannon = mean(shannon, na.rm = TRUE)
            )


som_list <- list()

names_vect <- c("Offshore Cluster", "Nearshore Cluster")
color_vect <- c("darkblue", "darkred")

# FOR NOW REMOVE NORTHERN TRANSECTS

som_maps <- som_maps[-c(1:12),]

map <- map_data("world")    

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
  scale_fill_gradient(low = "white", high = "red", aes(min = min(evenness, na.rm = TRUE, max = max(shannon, na.rm = TRUE)), limits = c(min,max))) +
  ggtitle(paste0("Evenness")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        title = element_text(hjust = 0.5))

############## HAB ID ############

taxas <- eight_tax_id[which(!is.na(match(eight_tax_id$Feature.ID, colnames(eighteen_s)))),]

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

sst <- ggplot() + 
  geom_tile(data = coeff_table, aes(x = lon, y = lat, fill = coeff_var), width =0.26, height = 0.26) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", limits = c(0.09,0.12), oob = squish, midpoint = 0.1066851) +geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("A. SST Coeff. Var.") +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
        title = element_text(hjust = 0.5), axis.line = element_blank())

sst_mean <- ggplot() + 
  geom_tile(data = mean_table, aes(x = lon, y = lat, fill = Mean), width =0.26, height = 0.26) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", limits = c(15,18), oob = squish, midpoint = 16.5) +geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("B. SST Mean (°C)") +
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

theme_set(theme_cowplot(font_size=8))

pdf("figures/3panel_NCOG_expanded_region_0604.pdf", width = 8, height = 3)
plot_grid(sst,p1,p2, ncol = 3, labels = c("A","B","C"), label_size = 10)
dev.off()

theme_set(theme_cowplot(font_size=10))

pdf("figures/4panel_NCOG_expand_ASV_0606.pdf", width = 8, height = 8)
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

som_maps <- som_maps[-which(is.na(som_maps$lat)),]


coeff_val <- vector()
mean_val <- vector()

for (i in 1:nrow(som_maps)) {
  
  lat_1 <- som_maps$lat[i]
  lon_1 <- som_maps$long[i]
  
  distance_vect <- sqrt((lat_1 - coeff_table$lat)^2 + (lon_1 - coeff_table$lon)^2)
   
  coeff_val[i] <- coeff_table$coeff_var[which.min(distance_vect)]
  mean_val[i] <- mean_table$Mean[which.min(distance_vect)]

  
}

som_maps$temp_mean <- mean_val
som_maps$temp_coeff <- coeff_val


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

mean_plot <- ggplot(som_plots, aes(x = mean, y = freq, color = cluster)) +
  geom_point(pch = 20, size = 4) +
  theme(legend.title = element_blank(),
        panel.border = element_rect(color = "black", linetype = "solid"),
        axis.line = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid")) +
  stat_smooth(method = "lm", level = 0.95) +
  xlab("Mean SST (°C)") + ylab("Frequency") +
  annotate("text", x = 15, y = 0.08, label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 2)

coeff_plot <- ggplot(som_plots, aes(x = coeff, y = freq, color = cluster)) +
  geom_point(pch = 20, size = 4) +
  theme(legend.title = element_blank(),
        panel.border = element_rect(color = "black", linetype = "solid"),
        axis.line = element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid"),
        legend.position = "none") +
  stat_smooth(method = "lm", level = 0.95) +
  xlab("SST Coeff. Var") + ylab("Frequency") +
  annotate("text", x = 0.095, y = 0.3, label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 2)

theme_set(theme_cowplot(font_size=8))

pdf("figures/lm_mean_coeff_asv_multiple_samps_0606.pdf", width = 8, height = 3)
plot_grid(coeff_plot, mean_plot, ncol = 2,rel_widths = c(0.8,1))
dev.off()

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

som_maps <- som_maps[,c(1:6,13,7:9,14,10:12)]

lm_mean_var <- lm(som_1 ~ lat + long + PO4_mean + NO3_mean + SiO3_mean + temp_mean + 
                    PO4_coeff + NO3_coeff + SiO3_coeff + temp_coeff, data = som_maps)

response_df <- som_maps$som_1
predictor_df <- som_maps[,5:14]

regsubsetsObj <- regsubsets(x=predictor_df ,y=response_df, nbest = 2, really.big = T)
regsubsetsObj$xnames <- c("(Intercept)", "Latitude", "Longitude", "Mean Temp", "Mean PO4",
                          "Mean NO3", "Mean SiO3", "Coeff. Var. Temp", "Coeff. Var. PO4",
                          "Coeff. Var. NO3", "Coeff. Var. SiO3")

pdf(file = "figures/stepwise_regression_leaps_0612.pdf", width = 6, height = 5)
plot(regsubsetsObj, scale = "adjr2", main = "Backwards Subset Selection")
dev.off()


glm_mean_var <- glm(som_1 ~ lat + long + PO4_mean + NO3_mean + SiO3_mean + temp_mean + 
      PO4_coeff + NO3_coeff + SiO3_coeff + temp_coeff, data = som_maps)

stepAIC(glm_mean_var)

model_AIC <- as.data.frame(matrix(NA,12,2))
colnames(model_AIC) <- c("Model","AIC")

# Temp

glm_mean_temp <- glm(som_1 ~  temp_mean, data = som_maps)
mt_sum <- summary(glm_mean_temp)
model_AIC[1,2] <- mt_sum$aic
model_AIC[1,1] <- "Mean Temp"

glm_coeff_temp <- glm(som_1 ~  temp_coeff, data = som_maps)
ct_sum <- summary(glm_coeff_temp)
model_AIC[2,2] <- ct_sum$aic
model_AIC[2,1] <- "Coeff Var. Temp"

# NO3

glm_mean_no3 <- glm(som_1 ~  NO3_mean, data = som_maps)
mn_sum <- summary(glm_mean_no3)
model_AIC[3,2] <- mn_sum$aic
model_AIC[3,1] <- "Mean NO3"

glm_coeff_no3 <- glm(som_1 ~  NO3_coeff, data = som_maps)
cn_sum <- summary(glm_coeff_no3)
model_AIC[4,2] <- cn_sum$aic
model_AIC[4,1] <- "Coeff Var. NO3"

# PO4

glm_mean_po4 <- glm(som_1 ~  PO4_mean, data = som_maps)
mp_sum <- summary(glm_mean_po4)
model_AIC[5,2] <- mp_sum$aic
model_AIC[5,1] <- "Mean PO4"

glm_coeff_po4 <- glm(som_1 ~  PO4_coeff, data = som_maps)
cp_sum <- summary(glm_coeff_po4)
model_AIC[6,2] <- cp_sum$aic
model_AIC[6,1] <- "Coeff Var. PO4"

# SiO3

glm_mean_sio3 <- glm(som_1 ~  SiO3_mean, data = som_maps)
ms_sum <- summary(glm_mean_sio3)
model_AIC[7,2] <- ms_sum$aic
model_AIC[7,1] <- "Mean SiO3"

glm_coeff_sio3 <- glm(som_1 ~  SiO3_coeff, data = som_maps)
cs_sum <- summary(glm_coeff_sio3)
model_AIC[8,2] <- cs_sum$aic
model_AIC[8,1] <- "Coeff Var. SiO3"

# Everything

glm_mean_var <- glm(som_1 ~ lat + long + PO4_mean + NO3_mean + SiO3_mean + temp_mean + 
                      PO4_coeff + NO3_coeff + SiO3_coeff + temp_coeff, data = som_maps)

all_sum <- summary(glm_mean_var)
model_AIC[9,2] <- all_sum$aic
model_AIC[9,1] <- "Full Model"

# best fit

glm_simple <- glm(som_1 ~ long + temp_mean + PO4_coeff + NO3_coeff + 
                      SiO3_coeff + temp_coeff, data = som_maps)

simple_sum <- summary(glm_simple)
model_AIC[10,2] <- simple_sum$aic
model_AIC[10,1] <- "Longitude + Mean Temp +\n Coeff Var. PO4 + Coeff Var. NO3 +\n Coeff Var. SiO3 + Coeff Var. temp"

# Lat

glm_lat <- glm(som_1 ~ lat, data = som_maps)

lat_sum <- summary(glm_lat)
model_AIC[11,2] <- lat_sum$aic
model_AIC[11,1] <- "Latitude"

# Lon

glm_lon <- glm(som_1 ~ long, data = som_maps)

lon_sum <- summary(glm_lon)
model_AIC[12,2] <- lon_sum$aic
model_AIC[12,1] <- "Longitude"

AIC_table <- model_AIC[order(model_AIC$AIC, decreasing = TRUE),]

pdf(file = "figures/AIC_table_0612.pdf", height = 5, width = 7)
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

forest_dat <- full_dat[,-which(is.na(rowSums(full_dat[,c(27,28,33,34,37,38,39,59,60)])))]

forest_dat$som_id <- as.factor(forest_dat$som_id)

forest_out <- randomForest(som_id ~ Lat_Dec + Lon_Dec + CC_Depth + T_degC + Salnty +
                             PO4ug + SiO3ug + NO3ug + Julian, data = forest_dat, importance = TRUE, na.action = na.omit)

pdf(file = "figures/random_forest_importance_0617.pdf", width = 5, height = 4)
varImpPlot(forest_out, type = 1, main = "Random Forest Variable Importance")
dev.off()



