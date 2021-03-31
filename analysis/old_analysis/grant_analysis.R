### Grant Analysis

library(tidyverse)
library(SOMbrero)
library(ggmap)
library(oce)
library(sp)
library(spatialEco)
library(patchwork)

NCOG_dat <- read.csv("data/NCOG_sample_log_DNA_meta_2014-2019.csv")

NCOG_dat$Line <- as.numeric(substr(NCOG_dat$Sta_ID,1,3))

# only look at southern cruises

NCOG_dat <- NCOG_dat %>%
  filter(Line > 75, Depthm < 15)

NCOG_dat <- NCOG_dat[,c(3,6,12,25,28,29,34:35,37:43,48,49)]

NCOG_dat$Spice <- swSpice(NCOG_dat$Salnty, NCOG_dat$T_degC)

NCOG_dat$Year <- as.numeric(substr(NCOG_dat$Cruise,1,4))

# run different soms

# first som is t,s

set.seed(15)

som_dat_1 <- NCOG_dat[,c(7:8)]
som_1_order <- which(complete.cases(som_dat_1) == TRUE)
som_dat_1 <- som_dat_1[complete.cases(som_dat_1),]

som.out.1 <- trainSOM(x.data = som_dat_1, dimension = c(6, 6), nb.save = 10, maxit = 2000)

som.1.clust <- superClass(som.out.1, k = 2)
clusters <- som.1.clust$cluster
ids <- som.1.clust$som$clustering
som_1_ids <- clusters[ids]

som_1_map <- NCOG_dat[som_1_order,]
som_1_map$som <- som_1_ids

# second Som is temp, salinity, 02, nutrients

som_dat_2 <- NCOG_dat[,c(7:13)]
som_2_order <- which(complete.cases(som_dat_2) == TRUE)
som_dat_2 <- som_dat_2[complete.cases(som_dat_2),]

som.out.2 <- trainSOM(x.data = som_dat_2, dimension = c(6, 6), nb.save = 10, maxit = 2000)

som.2.clust <- superClass(som.out.2, k = 2)
clusters <- som.2.clust$cluster
ids <- som.2.clust$som$clustering
som_2_ids <- clusters[ids]

som_2_map <- NCOG_dat[som_2_order,]
som_2_map$som <- som_2_ids

# third Som is spice, 02, nutrients

som_dat_3 <- NCOG_dat[,c(18,9:13)]
som_3_order <- which(complete.cases(som_dat_3) == TRUE)
som_dat_3 <- som_dat_3[complete.cases(som_dat_3),]

som.out.3 <- trainSOM(x.data = som_dat_3, dimension = c(6, 6), nb.save = 10, maxit = 2000)

som.3.clust <- superClass(som.out.3, k = 2)
clusters <- som.3.clust$cluster
ids <- som.3.clust$som$clustering
som_3_ids <- clusters[ids]

som_3_map <- NCOG_dat[som_3_order,]
som_3_map$som <- som_3_ids

som_list <- list(som1 = som_1_map, som2 = som_2_map, som3 = som_3_map)

map <- map_data("world") 
som_plots <- list()

for (i in 1:length(som_list)) {

  som_map <- som_list[[i]]

  warm_map <- som_map %>%
    filter(Year < 2017) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som == 1, na.rm = TRUE)/n(),
              som_2 = sum(som == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE),
              long = mean(Lon_Dec, na.rm = TRUE))

  cool_map <- som_map %>%
    filter(Year > 2016, Year < 2019) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som == 1, na.rm = TRUE)/n(),
              som_2 = sum(som == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE),
              long = mean(Lon_Dec, na.rm = TRUE))

  centroid_df <- SpatialPointsDataFrame(coords = cool_map[,c(6,5)], data = cool_map)

  wt_1 <- wt.centroid(x = centroid_df, p = 2)
  wt_2 <- wt.centroid(x = centroid_df, p = 3)

  clust1 <- which.max(c(wt_1@coords[1], wt_2@coords[1]))
  clust2 <- which.min(c(wt_1@coords[1], wt_2@coords[1]))



  if(clust1 == 1){
    cool_map$nearshore <- cool_map$som_1
    cool_map$offshore <- cool_map$som_2
    warm_map$nearshore <- warm_map$som_1
    warm_map$offshore <- warm_map$som_2
  }

  if(clust1 == 2){
    cool_map$nearshore <- cool_map$som_2
    cool_map$offshore <- cool_map$som_1
    warm_map$nearshore <- warm_map$som_2
    warm_map$offshore <- warm_map$som_1
  }

  warm <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = warm_map, aes(x = long, y = lat, fill = nearshore),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(low = "darkred", high = "darkblue",
                         mid = "white", limits = c(0,1), midpoint = 0.5) +
    ggtitle("a. Warm phase") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(size = 16), axis.line = element_blank(),
          axis.text = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    labs(fill = "Habitat probability")

  cool <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = cool_map, aes(x = long, y = lat, fill = nearshore),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(low = "darkred", high = "darkblue",
                         mid = "white", limits = c(0,1), midpoint = 0.5) +
    ggtitle("b. Cool phase") +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(size = 16), axis.line = element_blank(),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.position = "none")

  som_plots[[i]] <- warm + cool + plot_layout(guides = "collect", ncol = 1)

  }

pdf(file = "figures/t_s_nutrient_phase.pdf", width = 6, height = 10)
print(som_plots[[2]])
dev.off()



# by cruise

cruise_plots <- list()

for (i in 1:length(som_list)) {
  
  som_map <- som_list[[i]]
  
  c_2014 <- som_map %>%
    filter(Year == 2014) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som == 1, na.rm = TRUE)/n(),
              som_2 = sum(som == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE),
              long = mean(Lon_Dec, na.rm = TRUE))
  
  c_2015 <- som_map %>%
    filter(Year == 2015) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som == 1, na.rm = TRUE)/n(),
              som_2 = sum(som == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE),
              long = mean(Lon_Dec, na.rm = TRUE))
  
  c_2016 <- som_map %>%
    filter(Year == 2016) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som == 1, na.rm = TRUE)/n(),
              som_2 = sum(som == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE),
              long = mean(Lon_Dec, na.rm = TRUE))
  
  c_2017 <- som_map %>%
    filter(Year == 2017) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som == 1, na.rm = TRUE)/n(),
              som_2 = sum(som == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE),
              long = mean(Lon_Dec, na.rm = TRUE))
  
  c_2018 <- som_map %>%
    filter(Year == 2018) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som == 1, na.rm = TRUE)/n(),
              som_2 = sum(som == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE),
              long = mean(Lon_Dec, na.rm = TRUE))
  
  c_2019 <- som_map %>%
    filter(Year == 2019) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som == 1, na.rm = TRUE)/n(),
              som_2 = sum(som == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE),
              long = mean(Lon_Dec, na.rm = TRUE))
  
  all_map <- som_map %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som == 1, na.rm = TRUE)/n(),
              som_2 = sum(som == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE),
              long = mean(Lon_Dec, na.rm = TRUE))
  
  centroid_df <- SpatialPointsDataFrame(coords = all_map[,c(6,5)], data = all_map)
  
  wt_1 <- wt.centroid(x = centroid_df, p = 2)
  wt_2 <- wt.centroid(x = centroid_df, p = 3)
  
  clust1 <- which.max(c(wt_1@coords[1], wt_2@coords[1]))
  clust2 <- which.min(c(wt_1@coords[1], wt_2@coords[1]))
  
  
  
  if(clust1 == 1){
    c_2014$nearshore <- c_2014$som_1
    c_2014$offshore <- c_2014$som_2
    c_2015$nearshore <- c_2015$som_1
    c_2015$offshore <- c_2015$som_2
    c_2016$nearshore <- c_2016$som_1
    c_2016$offshore <- c_2016$som_2
    c_2017$nearshore <- c_2017$som_1
    c_2017$offshore <- c_2017$som_2
    c_2018$nearshore <- c_2018$som_1
    c_2018$offshore <- c_2018$som_2
    c_2019$nearshore <- c_2019$som_1
    c_2019$offshore <- c_2019$som_2
  }
  
  if(clust1 == 2){
    c_2014$nearshore <- c_2014$som_2
    c_2014$offshore <- c_2014$som_1
    c_2015$nearshore <- c_2015$som_2
    c_2015$offshore <- c_2015$som_1
    c_2016$nearshore <- c_2016$som_2
    c_2016$offshore <- c_2016$som_1
    c_2017$nearshore <- c_2017$som_2
    c_2017$offshore <- c_2017$som_1
    c_2018$nearshore <- c_2018$som_2
    c_2018$offshore <- c_2018$som_1
    c_2019$nearshore <- c_2019$som_2
    c_2019$offshore <- c_2019$som_1
  }
  
  p_2014 <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = c_2014, aes(x = long, y = lat, fill = nearshore),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(low = "darkred", high = "darkblue",
                         mid = "white", limits = c(0,1), midpoint = 0.5) +
    ggtitle("a. 2014") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(size = 16), axis.line = element_blank(),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    labs(fill = "Habitat probability")
  
  p_2015 <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = c_2015, aes(x = long, y = lat, fill = nearshore),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(low = "darkred", high = "darkblue",
                         mid = "white", limits = c(0,1), midpoint = 0.5) +
    ggtitle("a. 2015") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(size = 16), axis.line = element_blank(),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + labs(fill = "Habitat probability")
  
  p_2016 <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = c_2016, aes(x = long, y = lat, fill = nearshore),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(low = "darkred", high = "darkblue",
                         mid = "white", limits = c(0,1), midpoint = 0.5) +
    ggtitle("c. 2016") +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(size = 16), axis.line = element_blank(),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")
  
  p_2017 <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = c_2017, aes(x = long, y = lat, fill = nearshore),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(low = "darkred", high = "darkblue",
                         mid = "white", limits = c(0,1), midpoint = 0.5) +
    ggtitle("d. 2017") +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(size = 16), axis.line = element_blank(),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")
  
  p_2018 <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = c_2018, aes(x = long, y = lat, fill = nearshore),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(low = "darkred", high = "darkblue",
                         mid = "white", limits = c(0,1), midpoint = 0.5) +
    ggtitle("b. 2018") +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(size = 16), axis.line = element_blank(),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.position = "none")
  
  p_2019 <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = c_2019, aes(x = long, y = lat, fill = nearshore),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(low = "darkred", high = "darkblue",
                         mid = "white", limits = c(0,1), midpoint = 0.5) +
    ggtitle("f. 2019") +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(size = 16), axis.line = element_blank(),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")
  
  cruise_plots[[i]] <- p_2015 + p_2018 + plot_layout(ncol = 1, guides = "collect")
  
}

pdf(file = "figures/t_s_nutrient_year.pdf", width = 6, height = 10)
print(cruise_plots[[2]])
dev.off()






