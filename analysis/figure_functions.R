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
library(gganimate)
library(gifski)
library(chron)
library(purrr)
library(magick)
library(reshape2)
library(metR)
library(ncdf4)
library(patchwork)
library(MASS)
library(standardize)
library(mgcv)


#### SOM Figure ####

som_figure <- function(map_file = "output/bacteria_m_euks_16s_map.Rdata",
                       figure_name = paste0("figures/som_maps/bacteria_16s_som_",Sys.Date(),".pdf"),
                       main = "16s Bacteria", cluster1 = "Nearshore", cluster2 = "Offshore",
                       tsize = 12, psize = 4){
  
  
  
  map <- map_data("world")   
  
  load(map_file)
  
  som_maps <- som_maps %>% filter(substr(Sta_ID,2,3) > 75)
  
  # find centroids
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  wt_1 <- wt.centroid(x = centroid_df, p = 2)
  wt_2 <- wt.centroid(x = centroid_df, p = 3)
  
  clust1 <- which.max(c(wt_1@coords[1], wt_2@coords[1]))
  clust2 <- which.min(c(wt_1@coords[1], wt_2@coords[1]))
  
  p1 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",clust1)), color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "darkblue", limits = c(0,1)) +
    ggtitle(paste0(cluster1)) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank(),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize))
  
  if(clust1 == 1){p1 <- p1 + geom_point(aes(x = wt_1@coords[1], y = wt_1@coords[2]), color = "red", size = 5, pch = 10)}
  if(clust1 == 2){p1 <- p1 + geom_point(aes(x = wt_2@coords[1], y = wt_2@coords[2]), color = "red", size = 5, pch = 10)}
  
  p2 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",clust2)), color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "darkred", limits = c(0,1)) +
    ggtitle(paste0(cluster2)) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank(),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize))
  
  if(clust2 == 1){p2 <- p2 + geom_point(aes(x = wt_1@coords[1], y = wt_1@coords[2]), color = "blue", size = 5, pch = 10)}
  if(clust2 == 2){p2 <- p2 + geom_point(aes(x = wt_2@coords[1], y = wt_2@coords[2]), color = "blue", size = 5, pch = 10)}
  
  title <- ggdraw() + draw_label(main, fontface='bold')
  a <- plot_grid(title,plot_grid(p2,p1), nrow = 2,
                 rel_heights = c(0.1,1))
  
  pdf(file = figure_name, width = 6, height = 4)
  print(a)
  dev.off()
  
  blue_plot <- p2
  red_plot <- p1
  
  blue_plot <- blue_plot + ggtitle("")
  red_plot <- red_plot + ggtitle("")
  
  return(list(bp = blue_plot, rp = red_plot))
  
}

#### Regression Models ####

regression_figure <- function(glm_file = "output/bacteria_m_euks_16s_glm.Rdata",
                              map_file = "output/bacteria_m_euks_16s_map.Rdata",   
                              figure_name = "figures/glm_plots/bacteria_m_euks_16s_som_",
                              main = "16s Bacteria", cluster1 = "Nearshore", cluster2 = "Offshore",
                              var = "NC_mean", var_name = "Nitracline Depth (m)"){
  
  
  load(glm_file)
  load(map_file)
  
  # find centroids
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  wt_1 <- wt.centroid(x = centroid_df, p = 2)
  wt_2 <- wt.centroid(x = centroid_df, p = 3)
  
  clust1 <- which.max(c(wt_1@coords[1], wt_2@coords[1]))
  clust2 <- which.min(c(wt_1@coords[1], wt_2@coords[1]))
  
  if(clust1 == 1){som_plots$cluster[which(som_plots$cluster == "som_1")] <- "Offshore"}
  if(clust1 == 2){som_plots$cluster[which(som_plots$cluster == "som_2")] <- "Offshore"}
  
  if(clust2 == 1){som_plots$cluster[which(som_plots$cluster == "som_1")] <- "Nearshore"}
  if(clust2 == 2){som_plots$cluster[which(som_plots$cluster == "som_2")] <- "Nearshore"}
  
  
  reg_plot <- ggplot(som_plots, aes_string(x = var, y = "freq", color = "cluster")) + geom_point() +
    stat_smooth(method="glm", method.args = list(family = "binomial")) + 
    theme(legend.title = element_blank(), 
          panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    xlab(var_name) + ylab("Frequency") + scale_color_manual(values = c("blue", "red")) + ggtitle(main)
  
  
  
  pdf(file = paste0(figure_name,var,".pdf"), width = 5, height = 3)
  print(reg_plot)
  dev.off()
  
  return(reg_plot)
  
}

#### Diversity Maps ####

diveristy_figure <- function(map_file = "output/bacteria_m_euks_16s_map.Rdata",
                             full_dat = "output/bacteria_m_euks_16s_full_data.Rdata",
                             figure_start = "figures/diversity/bacteria_16s_",
                             main = "16s Bacteria"){
  
  
  map <- map_data("world")   
  
  load(map_file)
  load(full_dat)
  
  even <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes(x = long, y = lat, fill = evenness), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = mean(som_maps$evenness)) +
    ggtitle(paste0(main,"\nMean Evenness")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  rich <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes(x = long, y = lat, fill = richness), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = mean(som_maps$richness)) +
    ggtitle(paste0(main,"\nMean Richness")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  shannon <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes(x = long, y = lat, fill = shannon), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_viridis() +
    ggtitle(paste0(main,"\nMean Shannon Diversity")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  wt_1 <- wt.centroid(x = centroid_df, p = 2)
  wt_2 <- wt.centroid(x = centroid_df, p = 3)
  
  clust1 <- which.max(c(wt_1@coords[1], wt_2@coords[1]))
  clust2 <- which.min(c(wt_1@coords[1], wt_2@coords[1]))
  
  if(clust1 == 1){full_dat$som_id[which(full_dat$som_id == 1)] <- "Nearshore"}
  if(clust1 == 2){full_dat$som_id[which(full_dat$som_id == 2)] <- "Nearshore"}
  
  if(clust2 == 1){full_dat$som_id[which(full_dat$som_id == 1)] <- "Offshore"}
  if(clust2 == 2){full_dat$som_id[which(full_dat$som_id == 2)] <- "Offshore"}
  
  
  vio_even <- ggplot(full_dat, aes(x = som_id, y = evenness, fill = som_id)) +
    geom_violin() + geom_boxplot(width=0.1, fill = "white") +
    scale_fill_manual(values = c("red", "blue")) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),
          plot.title = element_text(hjust=0.5)) + xlab("") + ylab("Evenness") +
    ggtitle(paste0(main,"\nEvenness"))
  
  vio_rich <- ggplot(full_dat, aes(x = som_id, y = richness, fill = som_id)) +
    geom_violin() + geom_boxplot(width=0.1, fill = "white") +
    scale_fill_manual(values = c("red", "blue")) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),
          plot.title = element_text(hjust=0.5)) + xlab("") + ylab("Richness") +
    ggtitle(paste0(main,"\nRichness"))
  
  vio_shannon <- ggplot(full_dat, aes(x = som_id, y = shannon, fill = som_id)) +
    geom_violin() + geom_boxplot(width=0.1, fill = "white") +
    scale_fill_manual(values = c("red", "blue")) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),
          plot.title = element_text(hjust=0.5)) + xlab("") + ylab("Shannon Diversity") +
    ggtitle(paste0(main,"\nShannon Diversity"))
  
  pdf(file = paste0(figure_start,"even_map.pdf"), width = 4, height = 4)
  print(even)
  dev.off()
  
  pdf(file = paste0(figure_start,"rich_map.pdf"), width = 4, height = 4)
  print(rich)
  dev.off()
  
  pdf(file = paste0(figure_start,"shannon_map.pdf"), width = 4, height = 4)
  print(shannon)
  dev.off()
  
  pdf(file = paste0(figure_start,"even_vio.pdf"), width = 6, height = 4)
  print(vio_even)
  dev.off()
  
  pdf(file = paste0(figure_start,"rich_vio.pdf"), width = 6, height = 4)
  print(vio_rich)
  dev.off()
  
  pdf(file = paste0(figure_start,"shannon_vio.pdf"), width = 6, height = 4)
  print(vio_shannon)
  dev.off()
  
  return(shannon)
  
}

alpha_versus_gamma_figure <- function(full_data_file = "output/bacteria_m_euks_16s_full_data.Rdata",
                                      raw_data_file = "data/16s_bacteria_m_euks.Rdata",
                                      map_file = "output/bacteria_m_euks_16s_map.Rdata", minimum_tp = 8,
                                      figure_name = paste0("figures/bacteria_m_euks_16s_alpha_gamma_",Sys.Date(),".pdf"),
                                      main = "16s Bacteria"){
  
  load(full_data_file)
  load(raw_data_file)  
  load(map_file)
  
  som_maps <- som_maps[which(som_maps$n_samps > (minimum_tp-1)),]
  
  # Gamma Diversity
  
  asv_sums <- asv_table
  
  asv_sums <- as.data.frame(asv_sums)
  
  asv_sums$station <- full_dat$Sta_ID[match(rownames(asv_sums), full_dat$eco_name)]
  
  if(length(which(is.na(asv_sums$station))) > 0){asv_sums <- asv_sums[which(!is.na(asv_sums$station)),]}
  
  station_sums <- asv_sums %>% 
    group_by(station) %>% 
    summarise(across(.cols = everything(), .fns = ~sum(.x,na.rm = TRUE)))
  
  station_sums <- as.data.frame(station_sums)
  
  rownames(station_sums) <- station_sums$station
  
  station_sums$station <- NULL
  
  station_plot_df <- as.data.frame(matrix(NA,NROW(station_sums),4))
  
  colnames(station_plot_df) <- c("Station", "Diversity", "Latitude", "Longitude")
  
  station_plot_df$Station <- rownames(station_sums)
  
  # diversity
  
  station_plot_df$total_rich <- apply(station_sums, 1, function(x) length(which(x != 0)))
  
  station_plot_df$Diversity <- diversity(station_sums, MARGIN = 1, index = "shannon")
  
  # station location
  
  station_position <- full_dat %>% 
    group_by(Sta_ID) %>% 
    summarise(Lat = mean(Lat_Dec, na.rm = TRUE), Lon = mean(Lon_Dec, na.rm = TRUE))
  
  station_plot_df$Latitude <- station_position$Lat[match(station_plot_df$Station, station_position$Sta_ID)]
  station_plot_df$Longitude <- station_position$Lon[match(station_plot_df$Station, station_position$Sta_ID)]
  
  # remove northern transects
  
  line <- as.numeric(substr(station_plot_df$Station,2,3))
  n_stations <- which(line < 76)
  if(length(n_stations) > 0){station_plot_df <- station_plot_df[-n_stations,]}
  
  # Alpha Gamma Plots
  som_maps$Gamma_Diversity <- station_plot_df$Diversity[match(som_maps$Sta_ID, station_plot_df$Station)]
  som_maps$Total_Rich <- station_plot_df$total_rich[match(som_maps$Sta_ID, station_plot_df$Station)]
  
  # div_plot <- ggplot(som_maps, aes(x = Dist_mean)) + 
  #   geom_point(aes(y = richness, color = "Mean Richnness")) +
  #   stat_smooth(aes(y = richness), method = "loess", level = 0.95, color = "black") +
  #   geom_point(aes(y = Total_Rich, color = "Total Richness")) +
  #   stat_smooth(aes(y = Total_Rich), method = "loess", level = 0.95, color = "black") +
  #   scale_y_continuous(sec.axis = sec_axis(~., name = "Total Richness")) +
  #   ylab("Mean Richness") + xlab("Distance to Coast (km)") +
  #   theme(legend.title = element_blank(),
  #         panel.background = element_blank(),
  #         panel.border = element_rect(fill = NA, color = "black"),
  #         legend.position = "bottom", plot.title = element_text(hjust =0.5)) +
  #   ggtitle(paste0(main,"\nMean Richness vs\nTotal Richness per Station")) +
  #   scale_color_manual(values = c("royalblue2", "seagreen3")) 
  
  div_plot <- ggplot(som_maps, aes(x = Dist_mean)) + 
    geom_point(aes(y = shannon, color = "Mean Alpha Diversity")) +
    stat_smooth(aes(y = shannon), method = "loess", level = 0.95, color = "royalblue2") +
    geom_point(aes(y = Gamma_Diversity, color = "Gamma Diversity")) +
    stat_smooth(aes(y = Gamma_Diversity), method = "loess", level = 0.95, color = "seagreen3") +
    ylab("Shannon Index") + xlab("Distance to Coast (km)") +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.position = "bottom", plot.title = element_text(hjust =0.5)) +
    ggtitle(paste0(main,"\nAlpha Diversity vs\nGamma Diversity per Station")) +
    scale_color_manual(values = c("seagreen3","royalblue2")) 
  
  pdf(file = figure_name, width = 8, height = 6)
  print(div_plot)
  dev.off()
  
  return(div_plot)
  
}


beta_diversity_figure2 <- function(full_data_file = "output/bacteria_m_euks_16s_full_data.Rdata",
                                  bc_data_file = "output/bacteria_m_euks_16s_dissimilar.Rdata",
                                  map_file = "output/bacteria_m_euks_16s_map.Rdata", minimum_tp = 8,
                                  figure_name = paste0("figures/bacteria_m_euks_16s_beta_",Sys.Date(),".pdf"),
                                  main = "16s Bacteria"){
  
  
  load(full_data_file)
  load(bc_data_file)  
  load(map_file)
  
  # remove northern sites
  
  ns <- which(substr(colnames(dissimilar), 10,11) < 76)
  
  dissimilar <- dissimilar[-ns,-ns]
  
  # remove lower tri
  
  dissimilar[lower.tri(dissimilar, diag = TRUE)] <- NA
  
  # find cluster location
  
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  wt_1 <- wt.centroid(x = centroid_df, p = 2)
  wt_2 <- wt.centroid(x = centroid_df, p = 3)
  
  clust1 <- which.max(c(wt_1@coords[1], wt_2@coords[1]))
  clust2 <- which.min(c(wt_1@coords[1], wt_2@coords[1]))
  
  if(clust1 == 1){full_dat$som_id[which(full_dat$som_id == 1)] <- "Nearshore"}
  if(clust1 == 2){full_dat$som_id[which(full_dat$som_id == 2)] <- "Nearshore"}
  
  if(clust2 == 1){full_dat$som_id[which(full_dat$som_id == 1)] <- "Offshore"}
  if(clust2 == 2){full_dat$som_id[which(full_dat$som_id == 2)] <- "Offshore"}
  
  # som id
  
  dissimilar <- as.data.frame(dissimilar)
  
  dissimilar$som_id <- full_dat$som_id[match(colnames(dissimilar), full_dat$eco_name)]
  
  diss_df <- pivot_longer(dissimilar, -som_id, values_drop_na = TRUE)
  
  diss_df$name1 <- rownames(diss_df)
  
  diss_df$comp <- paste0(diss_df$som_id,"-",full_dat$som_id[match(diss_df$name, full_dat$eco_name)])
  
  diss_class <- diss_df[which(diss_df$comp == "Offshore-Offshore" | diss_df$comp == "Nearshore-Nearshore"),]
  
  diss_class$comp[diss_class$comp == "Nearshore-Nearshore"] = "Nearshore"
  diss_class$comp[diss_class$comp == "Offshore-Offshore"] = "Offshore"
  
  bc_plot <- ggplot(diss_class, aes(x = comp, y = value, fill = comp)) + geom_violin() +
    geom_boxplot(fill = "white", width = 0.05) + labs(fill = "") +
    scale_fill_manual(values = c("red","blue")) + xlab("") +
    ylab("Bray-Curtis\nDissimilarity") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
    ggtitle(main)
  
  diss_near <- diss_class[which(diss_class$comp == "Nearshore"),]
  diss_off <- diss_class[which(diss_class$comp == "Offshore"),]
  
  t_test <- t.test(diss_near$value, diss_off$value)
  
  print(t_test)
  
  pdf(file = figure_name, width = 4, height = 5)
  print(bc_plot)
  dev.off()
  
  return(bc_plot)
  
  
}

beta_diversity_figure <- function(full_data_file = "output/bacteria_m_euks_16s_full_data.Rdata",
                                  bc_data_file = "output/bacteria_m_euks_16s_dissimilar.Rdata",
                                  raw_data_file = "data/16s_bacteria_m_euks.Rdata",
                                  map_file = "output/bacteria_m_euks_16s_map.Rdata", minimum_tp = 8,
                                  figure_name = paste0("figures/bacteria_m_euks_16s_beta_",Sys.Date(),".pdf"),
                                  main = "16s Bacteria"){
  
  
  load(full_data_file)
  load(bc_data_file)  
  load(map_file)
  load(raw_data_file)
  
  map <- map_data("world")   
  
  colnames(dissimilar) <- rownames(asv_table)
  rownames(dissimilar) <- rownames(asv_table)
  
  # remove northern sites
  
  ns <- which(substr(colnames(dissimilar), 10,11) < 76)
  
  dissimilar <- dissimilar[-ns,-ns]
  
  # remove lower tri
  
  dissimilar[lower.tri(dissimilar, diag = TRUE)] <- NA
  
  # som id
  
  dissimilar <- as.data.frame(dissimilar)
  
  dissimilar$Sta_ID <- full_dat$Sta_ID[match(colnames(dissimilar), full_dat$eco_name)]
  
  diss_df <- pivot_longer(dissimilar, -Sta_ID, values_drop_na = TRUE)
  
  stationsplit <- strsplit(diss_df$name,"_")
  
  station_func <- function(df) paste0(df[2]," ",df[3])

  diss_df$name <- unlist(lapply(stationsplit, station_func))
  
  diss_sta <- diss_df[which(diss_df$Sta_ID == diss_df$name),]
  
  stations <- diss_sta %>%
    group_by(Sta_ID) %>%
    summarise(mean_beta = mean(value, na.rm = TRUE))
  
  som_maps$mean_beta <- (1 - stations$mean_beta[match(som_maps$Sta_ID, stations$Sta_ID)])
  
  bc_plot <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps %>% filter(n_samps > 30), aes(x = long, y = lat, fill = mean_beta), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = mean(som_maps$mean_beta)) +
    ggtitle(paste0(main,"\nMean Bray-Curtis Similarity")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = test, aes(x = long, y = lat, fill = mean_beta), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = mean(test$mean_beta)) +
    ggtitle(paste0(main,"\nMean Bray-Curtis Similarity")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  pdf(file = figure_name, width = 4, height = 5)
  print(bc_plot)
  dev.off()
  
  return(bc_plot)
  
  
}

### AIC Tables ####

aic_table_func <- function(som_maps = cyano_plots){  
  
  som_glm <- som_maps
  
  model_AIC <- matrix(NA,19,3)
  
  # temperature
  
  glm_mean_temp <- glm(som_1 ~  temp_mean, data = som_glm, family = binomial)
  mt_sum <- summary(glm_mean_temp)
  model_AIC[1,2] <- mt_sum$aic
  model_AIC[1,1] <- "Mean Temp"
  model_AIC[1,3] <- coef(mt_sum)[2,4]
  
  glm_coeff_temp <- glm(som_1 ~  temp_coeff, data = som_glm, family = binomial)
  ct_sum <- summary(glm_coeff_temp)
  model_AIC[10,2] <- ct_sum$aic
  model_AIC[10,1] <- "Coeff. Var. Temp"
  model_AIC[10,3] <- coef(ct_sum)[2,4]
  
  # sea surface temperature
  
  glm_mean_sst <- glm(som_1 ~  sst_mean, data = som_glm, family = binomial)
  mt_sum <- summary(glm_mean_sst)
  model_AIC[2,2] <- mt_sum$aic
  model_AIC[2,1] <- "Mean SST"
  model_AIC[2,3] <- coef(mt_sum)[2,4]
  
  glm_coeff_sst <- glm(som_1 ~  sst_coeff, data = som_glm, family = binomial)
  ct_sum <- summary(glm_coeff_sst)
  model_AIC[11,2] <- ct_sum$aic
  model_AIC[11,1] <- "Coeff. Var. SST"
  model_AIC[11,3] <- coef(ct_sum)[2,4]
  
  # salinity
  
  glm_mean_sal <- glm(som_1 ~  sal_mean, data = som_glm, family = binomial)
  mt_sum <- summary(glm_mean_sal)
  model_AIC[3,2] <- mt_sum$aic
  model_AIC[3,1] <- "Mean Salinity"
  model_AIC[3,3] <- coef(mt_sum)[2,4]
  
  glm_coeff_sal <- glm(som_1 ~  sal_coeff, data = som_glm, family = binomial)
  ct_sum <- summary(glm_coeff_sal)
  model_AIC[12,2] <- ct_sum$aic
  model_AIC[12,1] <- "Coeff. Var. Salinity"
  model_AIC[12,3] <- coef(ct_sum)[2,4]
  
  # NO3
  
  glm_mean_no3 <- glm(som_1 ~  NO3_mean, data = som_glm, family = binomial)
  mn_sum <- summary(glm_mean_no3)
  model_AIC[4,2] <- mn_sum$aic
  model_AIC[4,1] <- "Mean NO3"
  model_AIC[4,3] <- coef(mn_sum)[2,4]
  
  glm_coeff_no3 <- glm(som_1 ~  NO3_coeff, data = som_glm, family = binomial)
  cn_sum <- summary(glm_coeff_no3)
  model_AIC[13,2] <- cn_sum$aic
  model_AIC[13,1] <- "Coeff. Var. NO3"
  model_AIC[13,3] <- coef(cn_sum)[2,4]
  
  # PO4
  
  glm_mean_po4 <- glm(som_1 ~  PO4_mean, data = som_glm, family = binomial)
  mp_sum <- summary(glm_mean_po4)
  model_AIC[5,2] <- mp_sum$aic
  model_AIC[5,1] <- "Mean PO4"
  model_AIC[5,3] <- coef(mp_sum)[2,4]
  
  glm_coeff_po4 <- glm(som_1 ~  PO4_coeff, data = som_glm, family = binomial)
  cp_sum <- summary(glm_coeff_po4)
  model_AIC[14,2] <- cp_sum$aic
  model_AIC[14,1] <- "Coeff. Var. PO4"
  model_AIC[14,3] <- coef(cp_sum)[2,4]
  
  # SiO3
  
  glm_mean_sio3 <- glm(som_1 ~  SiO3_mean, data = som_glm, family = binomial)
  ms_sum <- summary(glm_mean_sio3)
  model_AIC[6,2] <- ms_sum$aic
  model_AIC[6,1] <- "Mean SiO4"
  model_AIC[6,3] <- coef(ms_sum)[2,4]
  
  glm_coeff_sio3 <- glm(som_1 ~  SiO3_coeff, data = som_glm, family = binomial)
  cs_sum <- summary(glm_coeff_sio3)
  model_AIC[15,2] <- cs_sum$aic
  model_AIC[15,1] <- "Coeff. Var. SiO4"
  model_AIC[15,3] <- coef(cs_sum)[2,4]
  
  # C14
  
  # glm_mean_c14 <- glm(som_1 ~  C14_mean, data = som_glm)
  # ms_sum <- summary(glm_mean_c14)
  # model_AIC[7,2] <- ms_sum$aic
  # model_AIC[7,1] <- "Mean C14"
  # 
  # glm_coeff_c14 <- glm(som_1 ~  C14_coeff, data = som_glm)
  # cs_sum <- summary(glm_coeff_c14)
  # model_AIC[17,2] <- cs_sum$aic
  # model_AIC[17,1] <- "Coeff. Var. C14"
  
  # SLA
  
  # glm_mean_sla <- glm(som_1 ~  sla_mean, data = som_glm, family = binomial)
  # ms_sum <- summary(glm_mean_sla)
  # model_AIC[7,2] <- ms_sum$aic
  # model_AIC[7,1] <- "Mean SLA"
  # 
  # glm_coeff_sla <- glm(som_1 ~  sla_coeff, data = som_glm, family = binomial)
  # cs_sum <- summary(glm_coeff_sla)
  # model_AIC[16,2] <- cs_sum$aic
  # model_AIC[16,1] <- "Coeff. Var. SLA"
  
  # MLD
  
  glm_mean_mld <- glm(som_1 ~  MLD_mean, data = som_glm, family = binomial)
  ms_sum <- summary(glm_mean_mld)
  model_AIC[7,2] <- ms_sum$aic
  model_AIC[7,1] <- "Mean MLD"
  model_AIC[7,3] <- coef(ms_sum)[2,4]
  
  glm_coeff_mld <- glm(som_1 ~  MLD_coeff, data = som_glm, family = binomial)
  cs_sum <- summary(glm_coeff_mld)
  model_AIC[16,2] <- cs_sum$aic
  model_AIC[16,1] <- "Coeff. Var. MLD"
  model_AIC[16,3] <- coef(cs_sum)[2,4]
  
  # NC Depth
  
  glm_mean_nc <- glm(som_1 ~  NC_mean, data = som_glm, family = binomial)
  ms_sum <- summary(glm_mean_nc)
  model_AIC[8,2] <- ms_sum$aic
  model_AIC[8,1] <- "Mean NCD"
  model_AIC[8,3] <- coef(ms_sum)[2,4]
  
  glm_coeff_nc <- glm(som_1 ~  NC_coeff, data = som_glm, family = binomial)
  cs_sum <- summary(glm_coeff_nc)
  model_AIC[17,2] <- cs_sum$aic
  model_AIC[17,1] <- "Coeff. Var. NCD"
  model_AIC[17,3] <- coef(cs_sum)[2,4]
  
  # Distance to Coast
  
  glm_mean_dc <- glm(som_1 ~  Dist_mean, data = som_glm, family = binomial)
  ms_sum <- summary(glm_mean_dc)
  model_AIC[18,2] <- ms_sum$aic
  model_AIC[18,1] <- "Distance to Coast"
  model_AIC[18,3] <- coef(ms_sum)[2,4]
  
  # Chlorophyll
  
  glm_mean_chl <- glm(som_1 ~  Chl_mean, data = som_glm, family = binomial)
  ms_sum <- summary(glm_mean_chl)
  model_AIC[9,2] <- ms_sum$aic
  model_AIC[9,1] <- "Mean Chl-a"
  model_AIC[9,3] <- coef(ms_sum)[2,4]
  
  glm_coeff_chl <- glm(som_1 ~  Chl_coeff, data = som_glm, family = binomial)
  cs_sum <- summary(glm_coeff_chl)
  model_AIC[19,2] <- cs_sum$aic
  model_AIC[19,1] <- "Coeff. Var. Chl-a"
  model_AIC[19,3] <- coef(cs_sum)[2,4]
  
  
  return(model_AIC)
  
}

aic_table_func_diveristy <- function(som_maps = cyano_plots, col_num = 2){  
  
  som_glm <- som_maps
  
  colnames(som_glm)[col_num] <- "response"
  
  model_AIC <- matrix(NA,19,3)
 
  
  
  # temperature
  
  glm_mean_temp <- glm(response ~  temp_mean, data = som_glm)
  mt_sum <- summary(glm_mean_temp)
  model_AIC[1,2] <- mt_sum$aic
  model_AIC[1,1] <- "Mean Temp"
  model_AIC[1,3] <- coef(mt_sum)[2,4]
  
  glm_coeff_temp <- glm(response ~  temp_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_temp)
  model_AIC[10,2] <- ct_sum$aic
  model_AIC[10,1] <- "Coeff. Var. Temp"
  model_AIC[10,3] <- coef(ct_sum)[2,4]
  
  # sea surface temperature
  
  glm_mean_sst <- glm(response ~  sst_mean, data = som_glm)
  mt_sum <- summary(glm_mean_sst)
  model_AIC[2,2] <- mt_sum$aic
  model_AIC[2,1] <- "Mean SST"
  model_AIC[2,3] <- coef(mt_sum)[2,4]
  
  glm_coeff_sst <- glm(response ~  sst_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_sst)
  model_AIC[11,2] <- ct_sum$aic
  model_AIC[11,1] <- "Coeff. Var. SST"
  model_AIC[11,3] <- coef(ct_sum)[2,4]
  
  # salinity
  
  glm_mean_sal <- glm(response ~  sal_mean, data = som_glm)
  mt_sum <- summary(glm_mean_sal)
  model_AIC[3,2] <- mt_sum$aic
  model_AIC[3,1] <- "Mean Salinity"
  model_AIC[3,3] <- coef(mt_sum)[2,4]
  
  glm_coeff_sal <- glm(response ~  sal_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_sal)
  model_AIC[12,2] <- ct_sum$aic
  model_AIC[12,1] <- "Coeff. Var. Salinity"
  model_AIC[12,3] <- coef(ct_sum)[2,4]
  
  # NO3
  
  glm_mean_no3 <- glm(response ~  NO3_mean, data = som_glm)
  mn_sum <- summary(glm_mean_no3)
  model_AIC[4,2] <- mn_sum$aic
  model_AIC[4,1] <- "Mean NO3"
  model_AIC[4,3] <- coef(mn_sum)[2,4]
  
  glm_coeff_no3 <- glm(response ~  NO3_coeff, data = som_glm)
  cn_sum <- summary(glm_coeff_no3)
  model_AIC[13,2] <- cn_sum$aic
  model_AIC[13,1] <- "Coeff. Var. NO3"
  model_AIC[13,3] <- coef(cn_sum)[2,4]
  
  # PO4
  
  glm_mean_po4 <- glm(response ~  PO4_mean, data = som_glm)
  mp_sum <- summary(glm_mean_po4)
  model_AIC[5,2] <- mp_sum$aic
  model_AIC[5,1] <- "Mean PO4"
  model_AIC[5,3] <- coef(mp_sum)[2,4]
  
  glm_coeff_po4 <- glm(response ~  PO4_coeff, data = som_glm)
  cp_sum <- summary(glm_coeff_po4)
  model_AIC[14,2] <- cp_sum$aic
  model_AIC[14,1] <- "Coeff. Var. PO4"
  model_AIC[14,3] <- coef(cp_sum)[2,4]
  
  # SiO3
  
  glm_mean_sio3 <- glm(response ~  SiO3_mean, data = som_glm)
  ms_sum <- summary(glm_mean_sio3)
  model_AIC[6,2] <- ms_sum$aic
  model_AIC[6,1] <- "Mean SiO4"
  model_AIC[6,3] <- coef(ms_sum)[2,4]
  
  glm_coeff_sio3 <- glm(response ~  SiO3_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_sio3)
  model_AIC[15,2] <- cs_sum$aic
  model_AIC[15,1] <- "Coeff. Var. SiO4"
  model_AIC[15,3] <- coef(cs_sum)[2,4]
  
  # C14
  
  # glm_mean_c14 <- glm(som_1 ~  C14_mean, data = som_glm)
  # ms_sum <- summary(glm_mean_c14)
  # model_AIC[7,2] <- ms_sum$aic
  # model_AIC[7,1] <- "Mean C14"
  # 
  # glm_coeff_c14 <- glm(som_1 ~  C14_coeff, data = som_glm)
  # cs_sum <- summary(glm_coeff_c14)
  # model_AIC[17,2] <- cs_sum$aic
  # model_AIC[17,1] <- "Coeff. Var. C14"
  
  # SLA
  
  # glm_mean_sla <- glm(response ~  sla_mean, data = som_glm)
  # ms_sum <- summary(glm_mean_sla)
  # model_AIC[7,2] <- ms_sum$aic
  # model_AIC[7,1] <- "Mean SLA"
  # 
  # glm_coeff_sla <- glm(response ~  sla_coeff, data = som_glm)
  # cs_sum <- summary(glm_coeff_sla)
  # model_AIC[16,2] <- cs_sum$aic
  # model_AIC[16,1] <- "Coeff. Var. SLA"
  
  # MLD
  
  glm_mean_mld <- glm(response ~  MLD_mean, data = som_glm)
  ms_sum <- summary(glm_mean_mld)
  model_AIC[7,2] <- ms_sum$aic
  model_AIC[7,1] <- "Mean MLD"
  model_AIC[7,3] <- coef(ms_sum)[2,4]
  
  glm_coeff_mld <- glm(response ~  MLD_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_mld)
  model_AIC[16,2] <- cs_sum$aic
  model_AIC[16,1] <- "Coeff. Var. MLD"
  model_AIC[16,3] <- coef(cs_sum)[2,4]
  
  # NC Depth
  
  glm_mean_nc <- glm(response ~  NC_mean, data = som_glm)
  ms_sum <- summary(glm_mean_nc)
  model_AIC[8,2] <- ms_sum$aic
  model_AIC[8,1] <- "Mean NCD"
  model_AIC[8,3] <- coef(ms_sum)[2,4]
  
  glm_coeff_nc <- glm(response ~  NC_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_nc)
  model_AIC[17,2] <- cs_sum$aic
  model_AIC[17,1] <- "Coeff. Var. NCD"
  model_AIC[17,3] <- coef(cs_sum)[2,4]
  
  # Distance to Coast
  
  glm_mean_dc <- glm(response ~  Dist_mean, data = som_glm)
  ms_sum <- summary(glm_mean_dc)
  model_AIC[18,2] <- ms_sum$aic
  model_AIC[18,1] <- "Distance to Coast"
  model_AIC[18,3] <- coef(ms_sum)[2,4]
  
  # Chlorophyll
  
  glm_mean_chl <- glm(response ~  Chl_mean, data = som_glm)
  ms_sum <- summary(glm_mean_chl)
  model_AIC[9,2] <- ms_sum$aic
  model_AIC[9,1] <- "Mean Chl-a"
  model_AIC[9,3] <- coef(ms_sum)[2,4]
  
  glm_coeff_chl <- glm(response ~  Chl_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_chl)
  model_AIC[19,2] <- cs_sum$aic
  model_AIC[19,1] <- "Coeff. Var. Chl-a"
  model_AIC[19,3] <- coef(cs_sum)[2,4]
  
  return(model_AIC)
  
}


full_aic_table_figure <- function(in_group_list = c("bacteria_m_euks_16s", "cyano_16s","plastid_16s",
                                                    "euks_auto_18sv9", "euks_hetero_18sv9", "flavo_16s",
                                                    "rhodo_16s", "sar_16s", "archaea_16s", "diatom_18sv9",
                                                    "dino_18sv9", "syndin_18sv9","hapto_18sv9", "metazoa_18sv9"),
                                  in_group_names = c("Heterotrophic\n Bacteria", "Cyanobacteria", "Eukaryotic\n Plastids",
                                                     "Eukaryotic\n Phytoplankton", "Heterotrophic\n Eukaryotes", 
                                                     "Flavobacteriales","Rhodobacterales", "Sar Clade", "Archaea", "Diatoms",
                                                     "Dinoflagellates", "Syndiniales", "Haptophytes", "Metazoans"),
                                  minimum_tp = 4, width_plot = 15,
                                  figure_name = paste0("figures/full_aic_table_logit_",Sys.Date(),".pdf"),
                                  figure_name_2 = paste0("figures/full_aic_plot_logit_",Sys.Date(),".pdf"),
                                  title_name = "Variable Importance", tsize = 12){
  # load data
  
  map_list <- list()
  p_val_list <- list()
  
  for (i in 1:length(in_group_list)) {
    
    load(paste0("output/",in_group_list[i], "_map.Rdata"))
    som_maps2 <- som_maps[which(som_maps$n_samps > (minimum_tp-1)),]
    AIC_table <- aic_table_func(som_maps = som_maps2)
    AIC_table <- as.data.frame(AIC_table)
    colnames(AIC_table) <- c("Variables","AIC","P_val")
    AIC_table <- AIC_table[c(1,3,4,5,6,9,8,10,12,13,14,15,19,17,18),]
    AIC_table[,2] <- as.numeric(as.character(AIC_table[,2]))
    AIC_table[,2] <- round(AIC_table[,2], digits = 2)
    
    AIC_table$sig <- "significant"
    AIC_table$sig[which(as.numeric(AIC_table$P_val) > 0.05)] <- "not significant"
    
    AIC_table$P_val <- NULL
    
    map_list[[i]] <- AIC_table[,1:2]
    p_val_list[[i]] <- AIC_table[,c(1,3)]
    
  }
  
  AIC_full <- full_join(map_list[[1]], map_list[[2]], by = "Variables")
  p_full <- full_join(p_val_list[[1]], p_val_list[[2]], by = "Variables")
  
  for (i in 3:length(map_list)) {
    
    AIC_full <- full_join(AIC_full, map_list[[i]], by = "Variables")
    p_full <- full_join(p_full, p_val_list[[i]], by = "Variables")
    
  }
  
  colnames(AIC_full) <- c("Variables", in_group_names)
  colnames(p_full) <- c("Variables", in_group_names)
  
  colfunc <- colorRampPalette(c("white", "red"))
  
  
  # phyto_col <- round(seq(max(phyto_AIC$AIC, na.rm = TRUE),min(phyto_AIC$AIC, na.rm = TRUE), by = -0.01), digits = 2)
  # scale <- colfunc(length(phyto_col))
  # phyto_fill <- scale[match(round(phyto_AIC$AIC, digits = 2), phyto_col)]
  # 
  # 
  # t0 <- tableGrob(AIC_full["Variables"], 
  #                 theme=ttheme_default(
  #                   core=list(bg_params = list(fill="grey90", col = "black"),
  #                             fg_params = list(fontface="bold")),
  #                   colhead = list(bg_params=list(fill="white", col="black"))), 
  #                 rows = NULL)
  # t1 <- tableGrob(AIC_full["Cyanobacteria\n AIC"],
  #                 theme=ttheme_default(
  #                   core=list(bg_params = list(fill=cyano_fill, col = "black")),
  #                   colhead = list(bg_params=list(fill="white", col="black"))),
  #                 rows = NULL)
  # 
  # # join tables
  # tab <- gtable_combine(t0,t1,t2,t3,t4)
  # 
  # pdf(file = figure_name, width = 10, height = 6)
  # grid.arrange(tab)
  # dev.off()
  
  AIC_scaled <- AIC_full
  
  for (i in 2:ncol(AIC_scaled)) {
    zero_one_scale <- 1-(AIC_scaled[,i]-min(AIC_scaled[,i], na.rm = TRUE))/
      abs(min(AIC_scaled[,i], na.rm = TRUE) - max(AIC_scaled[,i], na.rm = TRUE))
    AIC_scaled[,i] <- 50^zero_one_scale
    
  }
  
  p_df <- melt(p_full, id.vars = "Variables")
  colnames(p_df)[3] <- "p_val"
  
  plot_df <- melt(AIC_scaled)
  
  plot_df <- right_join(plot_df,p_df)

  colnames(plot_df) <- c("Variables", "Group", "AIC","P_Val")
  
  # # removing plastids from plot
  # plot_df <- plot_df[-which(plot_df$Group == "Eukaryotic\nPlastid\n AIC"),]
  
  plot_df$Variables <- as.factor(plot_df$Variables)
  plot_df$Variables <- factor(plot_df$Variables, levels = c("Distance to Coast",
                                                            "Coeff. Var. NCD","Coeff. Var. Chl-a",
                                                            "Coeff. Var. SiO4","Coeff. Var. PO4",
                                                            "Coeff. Var. NO3", "Coeff. Var. Salinity",
                                                            "Coeff. Var. Temp",
                                                            "Mean NCD", "Mean Chl-a",
                                                            "Mean SiO4", "Mean PO4",
                                                            "Mean NO3","Mean Salinity",
                                                            "Mean Temp" ))
  
  pdf(figure_name_2, width = width_plot, height = 8)
  
  print(ggplot(data = plot_df, aes(x = Group, y = Variables, size = AIC, fill = P_Val)) + 
          geom_point(color = "black", alpha = 0.6, shape = 21, show.legend = FALSE) +
          scale_fill_manual(values = c("grey90", "red")) +
          labs(size = "Variable\n Importance") + ylab("Variable") +
          theme(panel.background = element_blank(),
                panel.border = element_rect(color = "black", fill = NA),
                legend.position = "none",
                panel.grid.major.y = element_line(color = "grey", linetype = 2),
                axis.text.x = element_text(angle = 0),
                axis.text = element_text(size = tsize),
                axis.title = element_text(size = tsize)) +
          scale_size_continuous(range = c(1,18)) + xlab("") + ylab(""))
  
  dev.off()
  
  a <- ggplot(data = plot_df, aes(x = Group, y = Variables, size = AIC, fill = P_Val)) + 
    geom_point(color = "black", alpha = 0.6, shape = 21) +
    scale_fill_manual(values = c("grey90", "red")) +
    labs(size = "Variable\n Importance") + ylab("Variable") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          legend.position = "none",
          panel.grid.major.y = element_line(color = "grey", linetype = 2),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 0)) + ggtitle(title_name) +
    scale_size_continuous(range = c(1,18)) + xlab("") + ylab("")
  
  return(a)
  
}

full_aic_table_figure_diversity <- function(in_group_list = c("cyano_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                                                              "plastid_16s", "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                              "hapto_18sv9", "metazoa_18sv9"),
                                            in_group_names = c("Cyanobacteria", "Flavobacteriales","Rhodobacterales", "Sar Clade", "Archaea",
                                                               "Eukaryotic\n Plastids", "Diatoms",
                                                               "Dinoflagellates", "Syndiniales", "Haptophytes", "Metazoans"),
                                            figure_name = paste0("figures/full_aic_table_logit_even_",Sys.Date(),".pdf"),
                                            figure_name_2 = paste0("figures/group_aic_plot_logit_even_",Sys.Date(),".pdf"),
                                            title_name = "Variable Importance Evenness", # col 2 = even, col 3 = shan col 4 = rich
                                            col = 2, color_fill = "purple", width_plot = 15){
  # load data
  
  map_list <- list()
  p_val_list <- list()
  
  for (i in 1:length(in_group_list)) {
    
    load(paste0("output/",in_group_list[i], "_map.Rdata"))
    AIC_table <- aic_table_func_diveristy(som_maps = som_maps, col_num = col)
    AIC_table <- as.data.frame(AIC_table, stringsAsFactors = FALSE)
    colnames(AIC_table) <- c("Variables","AIC","P_val")
    AIC_table <- AIC_table[c(1,3,4,5,6,9,8,10,12,13,14,15,19,17,18),]
    AIC_table[,2] <- as.numeric(AIC_table[,2])
    AIC_table[,2] <- round(AIC_table[,2], digits = 2)
    
    AIC_table$sig <- "significant"
    AIC_table$sig[which(AIC_table$P_val > 0.05)] <- "not significant"
    
    AIC_table$P_val <- NULL
    
    map_list[[i]] <- AIC_table[,1:2]
    p_val_list[[i]] <- AIC_table[,c(1,3)]
    
  }
  
  AIC_full <- full_join(map_list[[1]], map_list[[2]], by = "Variables")
  p_full <- full_join(p_val_list[[1]], p_val_list[[2]], by = "Variables")
  
  for (i in 3:length(map_list)) {
    
    AIC_full <- full_join(AIC_full, map_list[[i]], by = "Variables")
    p_full <- full_join(p_full, p_val_list[[i]], by = "Variables")
    
  }
  
  colnames(AIC_full) <- c("Variables", in_group_names)
  colnames(p_full) <- c("Variables", in_group_names)
  
  
  
  # colfunc <- colorRampPalette(c("white", color_fill))
  
  # plastid_col <- round(seq(max(plastid_AIC$AIC, na.rm = TRUE),min(plastid_AIC$AIC, na.rm = TRUE), by = -0.01), digits = 2)
  # scale <- colfunc(length(plastid_col))
  # plastid_fill <- scale[match(round(plastid_AIC$AIC, digits = 2), plastid_col)]
  
  # cyano_col <- round(seq(max(cyano_AIC$AIC, na.rm = TRUE),min(cyano_AIC$AIC, na.rm = TRUE), by = -0.01), digits = 2)
  # scale <- colfunc(length(cyano_col))
  # cyano_fill <- scale[match(round(cyano_AIC$AIC, digits = 2), cyano_col)]
  # cyano_fill[5] <- scale[3903]
  # 
  # t0 <- tableGrob(AIC_full["Variables"], 
  #                 theme=ttheme_default(
  #                   core=list(bg_params = list(fill="grey90", col = "black"),
  #                             fg_params = list(fontface="bold")),
  #                   colhead = list(bg_params=list(fill="white", col="black"))), 
  #                 rows = NULL)
  # 
  # t1 <- tableGrob(AIC_full["Cyanobacteria\n AIC"],
  #                 theme=ttheme_default(
  #                   core=list(bg_params = list(fill=cyano_fill, col = "black")),
  #                   colhead = list(bg_params=list(fill="white", col="black"))),
  #                 rows = NULL)
  # 
  # # join tables
  # tab <- gtable_combine(t0,t1,t2,t3,t4)
  # 
  # pdf(file = figure_name, width = 10, height = 6)
  # grid.arrange(tab)
  # dev.off()
  
  AIC_scaled <- AIC_full
  
  for (i in 2:ncol(AIC_scaled)) {
    zero_one_scale <- 1-(AIC_scaled[,i]-min(AIC_scaled[,i], na.rm = TRUE))/
      abs(min(AIC_scaled[,i], na.rm = TRUE) - max(AIC_scaled[,i], na.rm = TRUE))
    AIC_scaled[,i] <- 50^zero_one_scale
    
  }
  
  plot_df <- melt(AIC_scaled)
  
  p_df <- melt(p_full, id.vars = "Variables")
  colnames(p_df)[3] <- "p_val"
  
  plot_df <- right_join(plot_df,p_df)
  
  colnames(plot_df) <- c("Variables", "Group", "AIC","P_Val")
 
  
  
  # # removing plastids from plot
  # plot_df <- plot_df[-which(plot_df$Group == "Eukaryotic\nPlastid\n AIC"),]
  
  plot_df$Variables <- as.factor(plot_df$Variables)
  plot_df$Variables <- factor(plot_df$Variables, levels = c("Distance to Coast",
                                                            "Coeff. Var. NCD","Coeff. Var. Chl-a",
                                                            "Coeff. Var. SiO4","Coeff. Var. PO4",
                                                            "Coeff. Var. NO3", "Coeff. Var. Salinity",
                                                            "Coeff. Var. Temp",
                                                            "Mean NCD", "Mean Chl-a",
                                                            "Mean SiO4", "Mean PO4",
                                                            "Mean NO3","Mean Salinity",
                                                            "Mean Temp" ))
  
  pdf(figure_name_2, width = width_plot, height = 8)
  print(ggplot(data = plot_df, aes(x = Group, y = Variables, size = AIC)) + 
          geom_point(fill = color_fill, color = "black", alpha = 0.6, shape = 21) +
          labs(size = "Variable\n Importance") + ylab("Variable") +
          theme(panel.background = element_blank(),
                panel.border = element_rect(color = "black", fill = NA),
                legend.position = "none",
                panel.grid.major.y = element_line(color = "grey", linetype = 2),
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 0)) + ggtitle(title_name) +
          scale_size_continuous(range = c(1,18)) + xlab("") + ylab(""))
  dev.off()
  
  
}


#### AIC Tables Sign ####

aic_table_func_diveristy_sign <- function(som_maps = cyano_plots, col_num = 2){  
  
  som_glm <- som_maps
  
  colnames(som_glm)[col_num] <- "response"
  
  model_AIC <- matrix(NA,19,4)
  # temperature
  
  glm_mean_temp <- glm(response ~  temp_mean, data = som_glm)
  mt_sum <- summary(glm_mean_temp)
  model_AIC[1,2] <- mt_sum$aic
  model_AIC[1,1] <- "Mean Temp"
  model_AIC[1,3] <- cor(som_glm$temp_mean, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[1,4] <- coef(mt_sum)[2,4]
  
  glm_coeff_temp <- glm(response ~  temp_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_temp)
  model_AIC[10,2] <- ct_sum$aic
  model_AIC[10,1] <- "Coeff. Var. Temp"
  model_AIC[10,3] <- cor(som_glm$temp_coeff, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[10,4] <- coef(ct_sum)[2,4]
  
  # sea surface temperature
  
  glm_mean_sst <- glm(response ~  sst_mean, data = som_glm)
  mt_sum <- summary(glm_mean_sst)
  model_AIC[2,2] <- mt_sum$aic
  model_AIC[2,1] <- "Mean SST"
  model_AIC[2,3] <- cor(som_glm$sst_mean, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[2,4] <- coef(mt_sum)[2,4]
  
  glm_coeff_sst <- glm(response ~  sst_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_sst)
  model_AIC[11,2] <- ct_sum$aic
  model_AIC[11,1] <- "Coeff. Var. SST"
  model_AIC[11,3] <- cor(som_glm$sst_coeff, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[11,4] <- coef(ct_sum)[2,4]
  
  # salinity
  
  glm_mean_sal <- glm(response ~  sal_mean, data = som_glm)
  mt_sum <- summary(glm_mean_sal)
  model_AIC[3,2] <- mt_sum$aic
  model_AIC[3,1] <- "Mean Salinity"
  model_AIC[3,3] <- cor(som_glm$sal_mean, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[3,4] <- coef(mt_sum)[2,4]
  
  glm_coeff_sal <- glm(response ~  sal_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_sal)
  model_AIC[12,2] <- ct_sum$aic
  model_AIC[12,1] <- "Coeff. Var. Salinity"
  model_AIC[12,3] <- cor(som_glm$sal_coeff, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[12,4] <- coef(ct_sum)[2,4]
  
  # NO3
  
  glm_mean_no3 <- glm(response ~  NO3_mean, data = som_glm)
  mn_sum <- summary(glm_mean_no3)
  model_AIC[4,2] <- mn_sum$aic
  model_AIC[4,1] <- "Mean NO3"
  model_AIC[4,3] <- cor(som_glm$NO3_mean, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[4,4] <- coef(mn_sum)[2,4]
  
  glm_coeff_no3 <- glm(response ~  NO3_coeff, data = som_glm)
  cn_sum <- summary(glm_coeff_no3)
  model_AIC[13,2] <- cn_sum$aic
  model_AIC[13,1] <- "Coeff. Var. NO3"
  model_AIC[13,3] <- cor(som_glm$NO3_coeff, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[13,4] <- coef(cn_sum)[2,4]
  
  # PO4
  
  glm_mean_po4 <- glm(response ~  PO4_mean, data = som_glm)
  mp_sum <- summary(glm_mean_po4)
  model_AIC[5,2] <- mp_sum$aic
  model_AIC[5,1] <- "Mean PO4"
  model_AIC[5,3] <- cor(som_glm$PO4_mean, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[5,4] <- coef(mp_sum)[2,4]
  
  glm_coeff_po4 <- glm(response ~  PO4_coeff, data = som_glm)
  cp_sum <- summary(glm_coeff_po4)
  model_AIC[14,2] <- cp_sum$aic
  model_AIC[14,1] <- "Coeff. Var. PO4"
  model_AIC[14,3] <- cor(som_glm$PO4_coeff, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[14,4] <- coef(cp_sum)[2,4]
  
  # SiO3
  
  glm_mean_sio3 <- glm(response ~  SiO3_mean, data = som_glm)
  ms_sum <- summary(glm_mean_sio3)
  model_AIC[6,2] <- ms_sum$aic
  model_AIC[6,1] <- "Mean SiO4"
  model_AIC[6,3] <- cor(som_glm$SiO3_mean, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[6,4] <- coef(ms_sum)[2,4]
  
  glm_coeff_sio3 <- glm(response ~  SiO3_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_sio3)
  model_AIC[15,2] <- cs_sum$aic
  model_AIC[15,1] <- "Coeff. Var. SiO4"
  model_AIC[15,3] <- cor(som_glm$SiO3_coeff, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[15,4] <- coef(cs_sum)[2,4]
  
  # MLD
  
  glm_mean_mld <- glm(response ~  MLD_mean, data = som_glm)
  ms_sum <- summary(glm_mean_mld)
  model_AIC[7,2] <- ms_sum$aic
  model_AIC[7,1] <- "Mean MLD"
  model_AIC[7,3] <- cor(som_glm$MLD_mean, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[7,4] <- coef(ms_sum)[2,4]
  
  glm_coeff_mld <- glm(response ~  MLD_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_mld)
  model_AIC[16,2] <- cs_sum$aic
  model_AIC[16,1] <- "Coeff. Var. MLD"
  model_AIC[16,3] <- cor(som_glm$MLD_coeff, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[16,4] <- coef(cs_sum)[2,4]
  
  # NC Depth
  
  glm_mean_nc <- glm(response ~  NC_mean, data = som_glm)
  ms_sum <- summary(glm_mean_nc)
  model_AIC[8,2] <- ms_sum$aic
  model_AIC[8,1] <- "Mean NCD"
  model_AIC[8,3] <- cor(som_glm$NC_mean, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[8,4] <- coef(ms_sum)[2,4]
  
  glm_coeff_nc <- glm(response ~  NC_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_nc)
  model_AIC[17,2] <- cs_sum$aic
  model_AIC[17,1] <- "Coeff. Var. NCD"
  model_AIC[17,3] <- cor(som_glm$NC_coeff, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[17,4] <- coef(cs_sum)[2,4]
  
  # Distance to Coast
  
  glm_mean_dc <- glm(response ~  Dist_mean, data = som_glm)
  ms_sum <- summary(glm_mean_dc)
  model_AIC[18,2] <- ms_sum$aic
  model_AIC[18,1] <- "Distance to Coast"
  model_AIC[18,3] <- cor(som_glm$Dist_mean, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[18,4] <- coef(ms_sum)[2,4]
  
  # Chlorophyll
  
  glm_mean_chl <- glm(response ~  Chl_mean, data = som_glm)
  ms_sum <- summary(glm_mean_chl)
  model_AIC[9,2] <- ms_sum$aic
  model_AIC[9,1] <- "Mean Chl-a"
  model_AIC[9,3] <- cor(som_glm$Chl_mean, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[9,4] <- coef(ms_sum)[2,4]
  
  glm_coeff_chl <- glm(response ~  Chl_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_chl)
  model_AIC[19,2] <- cs_sum$aic
  model_AIC[19,1] <- "Coeff. Var. Chl-a"
  model_AIC[19,3] <- cor(som_glm$Chl_coeff, som_glm$response, use = "pairwise.complete.obs")
  model_AIC[19,4] <- coef(cs_sum)[2,4]
  
  return(model_AIC)
  
}


full_aic_table_figure_diversity_sign <- function(in_group_list = c("cyano_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                                                                   "plastid_16s", "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                                   "hapto_18sv9", "metazoa_18sv9"),
                                                 in_group_names = c("Cyanobacteria", "Flavobacteriales","Rhodobacterales", "Sar Clade", "Archaea",
                                                                    "Eukaryotic\n Plastids", "Diatoms",
                                                                    "Dinoflagellates", "Syndiniales", "Haptophytes", "Metazoans"),
                                                 figure_name_2 = paste0("figures/group_aic_plot_logit_even_",Sys.Date(),".pdf"),
                                                 col = 27, width_plot = 15, tsize = 12){
  # load data
  
  map_list <- list()
  p_val_list <- list()
  
  for (i in 1:length(in_group_list)) {
    
    load(paste0("output/",in_group_list[i], "_map.Rdata"))
    AIC_table <- aic_table_func_diveristy_sign(som_maps = som_maps, col_num = col)
    AIC_table <- as.data.frame(AIC_table, stringsAsFactors = FALSE)
    colnames(AIC_table) <- c("Variables","AIC", "Slope", "P_val")
    AIC_table <- AIC_table[c(1,3,4,5,6,9,8,10,12,13,14,15,19,17,18),]
    AIC_table[,2] <- as.numeric(AIC_table[,2])
    AIC_table[,2] <- round(AIC_table[,2], digits = 2)
    AIC_table[,3] <- as.numeric(AIC_table[,3])
    AIC_table[,3] <- round(AIC_table[,3], digits = 3)
    
    AIC_table$sig <- "significant"
    AIC_table$sig[which(as.numeric(AIC_table$P_val) > 0.05)] <- "not significant"
    
    AIC_table$P_val <- NULL
    
    map_list[[i]] <- AIC_table[,1:3]
    p_val_list[[i]] <- AIC_table[,c(1,4)]
    
  }
  
  AIC_full <- full_join(map_list[[1]], map_list[[2]], by = "Variables")
  p_full <- full_join(p_val_list[[1]], p_val_list[[2]], by = "Variables")
  
  for (i in 3:length(map_list)) {
    
    AIC_full <- full_join(AIC_full, map_list[[i]], by = "Variables")
    p_full <- full_join(p_full, p_val_list[[i]], by = "Variables")
  }
  
  
  colnames(p_full) <- c("Variables", in_group_names)
  
  AIC_full <- AIC_full[,c(1,seq(2,ncol(AIC_full),by = 2),seq(3,ncol(AIC_full),by = 2))]
  
  colnames(AIC_full) <- c("Variables", in_group_names, paste0(in_group_names,"_slope"))
  
  AIC_scaled <- AIC_full
  
  for (i in 2:(length(in_group_names)+1)) {
    
    zero_one_scale <- 1-(AIC_scaled[,i]-min(AIC_scaled[,i], na.rm = TRUE))/
      abs(min(AIC_scaled[,i], na.rm = TRUE) - max(AIC_scaled[,i], na.rm = TRUE))
   
    
    AIC_scaled[,i] <- 50^zero_one_scale
    
  }
  
  plot_df <- melt(AIC_scaled[,1:(length(in_group_names)+1)])
  
  plot_slope <- melt(AIC_scaled[,c(1,(length(in_group_names)+2):ncol(AIC_scaled))])
  
  colnames(plot_slope)[3] <- "slope"
  plot_slope$variable <- sub("_.*","",plot_slope$variable)
  plot_df <- right_join(plot_df,plot_slope)
  
  p_df <- melt(p_full, id.vars = "Variables")
  colnames(p_df)[3] <- "p_val"
  plot_df <- right_join(plot_df,p_df)
  
  colnames(plot_df) <- c("Variables", "Group", "AIC","slope","p_val")
  
  plot_df$slope[which(plot_df$p_val == "not significant")] <- NA
  
  # # removing plastids from plot
  # plot_df <- plot_df[-which(plot_df$Group == "Eukaryotic\nPlastid\n AIC"),]
  
  plot_df$Variables <- as.factor(plot_df$Variables)
  plot_df$Variables <- factor(plot_df$Variables, levels = c("Distance to Coast",
                                                            "Coeff. Var. NCD","Coeff. Var. Chl-a",
                                                            "Coeff. Var. SiO4","Coeff. Var. PO4",
                                                            "Coeff. Var. NO3", "Coeff. Var. Salinity",
                                                            "Coeff. Var. Temp",
                                                            "Mean NCD", "Mean Chl-a",
                                                            "Mean SiO4", "Mean PO4",
                                                            "Mean NO3","Mean Salinity",
                                                            "Mean Temp" ))
  
  
  plot_df$Group <- as.factor(plot_df$Group)
  plot_df$Group <- factor(plot_df$Group, levels = in_group_names)
  
  pdf(figure_name_2, width = width_plot, height = 8)
  print(ggplot(data = plot_df, aes(x = Group, y = Variables, size = AIC, fill = slope)) + 
          geom_point(color = "black", alpha = 0.6, shape = 21) +
          scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "white",
                               midpoint = 0, limits = c(-1,1), na.value = "grey90") +
          labs(fill = "Correlation") + ylab("Variable") +
          theme(panel.background = element_blank(),
                panel.border = element_rect(color = "black", fill = NA),
                legend.position = "right",
                panel.grid.major.y = element_line(color = "grey", linetype = 2),
                plot.title = element_text(hjust = 0, size = tsize),
                axis.text.y = element_text(size = tsize),
                legend.title = element_text(size = tsize),
                legend.text = element_text(size = tsize),
                axis.text.x = element_text(size = tsize, angle = 0)) +
          scale_size_continuous(range = c(1,18), guide = FALSE) + xlab("") + ylab(""))
  dev.off()
  
  p <- ggplot(data = plot_df, aes(x = Group, y = Variables, size = AIC, fill = slope)) + 
    geom_point(color = "black", alpha = 0.6, shape = 21) +
    scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "white",
                         midpoint = 0, limits = c(-1,1), na.value = "grey90") +
    labs(fill = "Correlation") + ylab("Variable") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          legend.position = "right",
          panel.grid.major.y = element_line(color = "grey", linetype = 2),
          plot.title = element_text(hjust = 0, size = tsize),
          axis.text.y = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.text.x = element_text(size = tsize, angle = 0)) +
    scale_size_continuous(range = c(1,18), guide = FALSE) + xlab("") + ylab("")
  
  return(p)
  
}

##### Community Differences Plot #####

fig_commun_map_func <- function(in_all = "output/total_dissimilar.Rdata",
                       in_dat = "output/total_full_data.Rdata",
                       in_map = "output/total_map.Rdata",
                       community_diff_fig = "figures/figure_outline/fig_x.pdf",
                       tsize = 12, psize = 12, group = "Total ASVs"){
  
  load(in_all)
  load(in_dat)
  load(in_map)
  
  samps <- rownames(dissimilar)
  mat <- matrix(NA, nrow(dissimilar), 6)
  mat[,1] <- as.numeric(substr(samps,2,5))
  mat[,2] <- substr(samps,6,7)
  mat[,4] <- match(samps,full_dat$eco_name)
  
  mat[which(as.numeric(mat[,2]) < 3),3] <- "Winter" 
  mat[which(as.numeric(mat[,2]) == 4),3] <- "Spring" 
  mat[which(as.numeric(mat[,2]) > 5 & as.numeric(mat[,2]) < 9),3] <- "Summer"
  mat[which(as.numeric(mat[,2]) > 8),3] <- "Fall" 
  
  mat[which(as.numeric(mat[,1]) < 2017),5] <- "Early"
  mat[which(as.numeric(mat[,1]) > 2016 & as.numeric(mat[,1]) < 2019),5] <- "Late"
  
  depth_list <- strsplit(samps, "_")
  mat[,6] <- as.numeric(sapply(depth_list, "[[",4))
  
  all_diss <- dissimilar[which(mat[,5] == "Early" & as.numeric(mat[,6]) < 16),
                         which(mat[,5] == "Late" & as.numeric(mat[,6]) < 16)]
  
  all_comp <- as.data.frame(all_diss)
  all_comp$early_samps <- rownames(all_comp)
  
  all_long <- all_comp %>%
    pivot_longer(-early_samps,
                 names_to = "late_samps",
                 values_to = "b_c_vals")
  
  all_long <- all_long[which(substr(all_long$early_samps,9,19)==substr(all_long$late_samps,9,19)),]

  all_long$Sta_ID <- paste0(substr(all_long$early_samps,9,13)," ",substr(all_long$early_samps,15,19))
  
  all_mean <- all_long %>%
    group_by(Sta_ID) %>%
    summarise(mean_bc = mean(b_c_vals, na.rm = TRUE))
  
  som_maps$bc_vals <- 1-all_mean$mean_bc[match(som_maps$Sta_ID, all_mean$Sta_ID)]
  
  map <- map_data("world")

  example_plot <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = som_maps[complete.cases(som_maps$bc_vals),],
               aes_string(x = "long", y = "lat", fill = "bc_vals"),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Bray-Curtis\nSimilarity", low = "white", high = "red") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,1),
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle(group,paste0("Bray-Curtis Similarity\n(2014-2016) vs (2017-2018)"))
  
  
  pdf(file = community_diff_fig, width = 8, height = 8)
  print(example_plot)
  dev.off()
  
  
}

fig_commun_map_surf_deep_func <- function(in_all = "output/total_dissimilar.Rdata",
                                in_dat = "output/total_full_data.Rdata",
                                in_map = "output/total_map.Rdata",
                                community_diff_fig = "figures/figure_outline/fig_x2.pdf",
                                tsize = 12, psize = 12, group = "Total ASVs"){
  
  load(in_all)
  load(in_dat)
  load(in_map)
  
  samps <- rownames(dissimilar)
  mat <- matrix(NA,nrow(dissimilar), 6)
  mat[,1] <- as.numeric(substr(samps,2,5))
  mat[,2] <- substr(samps,6,7)
  mat[,4] <- match(samps,full_dat$eco_name)
  
  mat[which(as.numeric(mat[,2]) < 3),3] <- "Winter" 
  mat[which(as.numeric(mat[,2]) == 4),3] <- "Spring" 
  mat[which(as.numeric(mat[,2]) > 5 & as.numeric(mat[,2]) < 9),3] <- "Summer"
  mat[which(as.numeric(mat[,2]) > 8),3] <- "Fall" 
  
  mat[which(as.numeric(mat[,1]) < 2017),5] <- "Early"
  mat[which(as.numeric(mat[,1]) > 2016 & as.numeric(mat[,1]) < 2019),5] <- "Late"
  
  depth_list <- strsplit(samps, "_")
  mat[,6] <- as.numeric(sapply(depth_list, "[[",4))
  
  surf_diss <- dissimilar[which(mat[,5] == "Early" & as.numeric(mat[,6]) < 16),
                         which(mat[,5] == "Late" & as.numeric(mat[,6]) < 16)]
  
  deep_diss <- dissimilar[which(mat[,5] == "Early" & as.numeric(mat[,6]) > 15),
                          which(mat[,5] == "Late" & as.numeric(mat[,6]) > 15)]
  
  surf_warm <- dissimilar[which(mat[,5] == "Early" & as.numeric(mat[,6]) < 16),
                          which(mat[,5] == "Early" & as.numeric(mat[,6]) < 16)]
  
  deep_warm <- dissimilar[which(mat[,5] == "Early" & as.numeric(mat[,6]) > 15),
                          which(mat[,5] == "Early" & as.numeric(mat[,6]) > 15)]
  
  surf_cool <- dissimilar[which(mat[,5] == "Late" & as.numeric(mat[,6]) < 16),
                          which(mat[,5] == "Late" & as.numeric(mat[,6]) < 16)]
  
  deep_cool <- dissimilar[which(mat[,5] == "Late" & as.numeric(mat[,6]) > 15),
                          which(mat[,5] == "Late" & as.numeric(mat[,6]) > 15)]
  
  diag(surf_warm) <- NA
  diag(deep_warm) <- NA
  diag(surf_cool) <- NA
  diag(deep_cool) <- NA
  
  surf_comp <- as.data.frame(surf_diss)
  surf_comp$early_samps <- rownames(surf_comp)
  
  deep_comp <- as.data.frame(deep_diss)
  deep_comp$early_samps <- rownames(deep_comp)
  
  surf_warmin <- as.data.frame(surf_warm)
  surf_warmin$early_samps <- rownames(surf_warmin)
  
  deep_warmin <- as.data.frame(deep_warm)
  deep_warmin$early_samps <- rownames(deep_warmin)
  
  surf_coolin <- as.data.frame(surf_cool)
  surf_coolin$early_samps <- rownames(surf_coolin)
  
  deep_coolin <- as.data.frame(deep_cool)
  deep_coolin$early_samps <- rownames(deep_coolin)
  
  surf_long <- surf_comp %>%
    pivot_longer(-early_samps,
                 names_to = "late_samps",
                 values_to = "b_c_vals")
  
  deep_long <- deep_comp %>%
    pivot_longer(-early_samps,
                 names_to = "late_samps",
                 values_to = "b_c_vals")
  
  surf_warm_long <- surf_warmin %>%
    pivot_longer(-early_samps,
                 names_to = "late_samps",
                 values_to = "b_c_vals")
  
  deep_warm_long <- deep_warmin %>%
    pivot_longer(-early_samps,
                 names_to = "late_samps",
                 values_to = "b_c_vals")
  
  surf_cool_long <- surf_coolin %>%
    pivot_longer(-early_samps,
                 names_to = "late_samps",
                 values_to = "b_c_vals")
  
  deep_cool_long <- deep_coolin %>%
    pivot_longer(-early_samps,
                 names_to = "late_samps",
                 values_to = "b_c_vals")
  
  surf_long <- surf_long[which(substr(surf_long$early_samps,9,19)==substr(surf_long$late_samps,9,19)),]
  deep_long <- deep_long[which(substr(deep_long$early_samps,9,19)==substr(deep_long$late_samps,9,19)),]
  
  surf_warm_long <- surf_warm_long[which(substr(surf_warm_long$early_samps,9,19)==substr(surf_warm_long$late_samps,9,19)),]
  deep_warm_long <- deep_warm_long[which(substr(deep_warm_long$early_samps,9,19)==substr(deep_warm_long$late_samps,9,19)),]
  
  surf_cool_long <- surf_cool_long[which(substr(surf_cool_long$early_samps,9,19)==substr(surf_cool_long$late_samps,9,19)),]
  deep_cool_long <- deep_cool_long[which(substr(deep_cool_long$early_samps,9,19)==substr(deep_cool_long$late_samps,9,19)),]
  
  surf_long$Sta_ID <- paste0(substr(surf_long$early_samps,9,13)," ",substr(surf_long$early_samps,15,19))
  deep_long$Sta_ID <- paste0(substr(deep_long$early_samps,9,13)," ",substr(deep_long$early_samps,15,19))
  
  surf_warm_long$Sta_ID <- paste0(substr(surf_warm_long$early_samps,9,13)," ",substr(surf_warm_long$early_samps,15,19))
  deep_warm_long$Sta_ID <- paste0(substr(deep_warm_long$early_samps,9,13)," ",substr(deep_warm_long$early_samps,15,19))
  
  surf_cool_long$Sta_ID <- paste0(substr(surf_cool_long$early_samps,9,13)," ",substr(surf_cool_long$early_samps,15,19))
  deep_cool_long$Sta_ID <- paste0(substr(deep_cool_long$early_samps,9,13)," ",substr(deep_cool_long$early_samps,15,19))
  
  surf_mean <- surf_long %>%
    group_by(Sta_ID) %>%
    summarise(mean_bc = mean(b_c_vals, na.rm = TRUE))
  
  deep_mean <- deep_long %>%
    group_by(Sta_ID) %>%
    summarise(mean_bc = mean(b_c_vals, na.rm = TRUE))
  
  surf_warm_mean <- surf_warm_long %>%
    group_by(Sta_ID) %>%
    summarise(mean_bc = mean(b_c_vals, na.rm = TRUE))
  
  deep_warm_mean <- deep_warm_long %>%
    group_by(Sta_ID) %>%
    summarise(mean_bc = mean(b_c_vals, na.rm = TRUE))
  
  surf_cool_mean <- surf_cool_long %>%
    group_by(Sta_ID) %>%
    summarise(mean_bc = mean(b_c_vals, na.rm = TRUE))
  
  deep_cool_mean <- deep_cool_long %>%
    group_by(Sta_ID) %>%
    summarise(mean_bc = mean(b_c_vals, na.rm = TRUE))
  
  som_maps$surf_bc_vals <- 1-surf_mean$mean_bc[match(som_maps$Sta_ID, surf_mean$Sta_ID)]
  som_maps$deep_bc_vals <- 1-deep_mean$mean_bc[match(som_maps$Sta_ID, deep_mean$Sta_ID)]
  som_maps$surf_warm_bc_vals <- 1-surf_warm_mean$mean_bc[match(som_maps$Sta_ID, surf_warm_mean$Sta_ID)]
  som_maps$deep_warm_bc_vals <- 1-deep_warm_mean$mean_bc[match(som_maps$Sta_ID, deep_warm_mean$Sta_ID)]
  som_maps$surf_cool_bc_vals <- 1-surf_cool_mean$mean_bc[match(som_maps$Sta_ID, surf_cool_mean$Sta_ID)]
  som_maps$deep_cool_bc_vals <- 1-deep_cool_mean$mean_bc[match(som_maps$Sta_ID, deep_cool_mean$Sta_ID)]
  
  som_maps$bc_surf_bw_within <- som_maps$surf_bc_vals/som_maps$surf_warm_bc_vals
  som_maps$bc_deep_bw_within <- som_maps$deep_bc_vals/som_maps$deep_warm_bc_vals
  
  val_10 <- quantile(c(som_maps$surf_bc_vals, som_maps$deep_bc_vals), na.rm = TRUE, probs = seq(0,1,0.1))[2]
  val_90 <- quantile(c(som_maps$surf_bc_vals, som_maps$deep_bc_vals), na.rm = TRUE, probs = seq(0,1,0.1))[10]
  
  map <- map_data("world")
  
  plot_a <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = som_maps[complete.cases(som_maps$surf_bc_vals),],
               aes_string(x = "long", y = "lat", fill = "surf_bc_vals"),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Bray-Curtis\nSimilarity", low = "white", high = "red",
                        limits = c(val_10,val_90), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,1),
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle(group,paste0("Surface\n(2014-2016) vs (2017-2018)"))
  
  plot_b <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = som_maps[complete.cases(som_maps$deep_bc_vals),],
               aes_string(x = "long", y = "lat", fill = "deep_bc_vals"),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Bray-Curtis\nSimilarity", low = "white", high = "red",
                        limits = c(val_10,val_90), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,1),
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("DCM\n(2014-2016) vs (2017-2018)")
  
  example_plot <- plot_a + plot_b + plot_layout(guides = "collect")
  
  bc_mat <- matrix(NA,6*nrow(som_maps),7)
  bc_mat <- as.data.frame(bc_mat)
  colnames(bc_mat) <- c("bc_val", "Group", "Depth", "Group_Depth", "Distance", "NCDepth", "Phase")
  
  bc_mat$bc_val <- c(som_maps$surf_bc_vals, som_maps$deep_bc_vals,
                     som_maps$surf_warm_bc_vals, som_maps$deep_warm_bc_vals,
                     som_maps$surf_cool_bc_vals, som_maps$deep_cool_bc_vals)
  bc_mat$Group <- group
  bc_mat$Depth <- rep(c(rep("Surface", nrow(som_maps)),rep("DCM", nrow(som_maps))),3)
  bc_mat$Group_Depth <- rep(c(paste0(group, " ",rep("Surface", nrow(som_maps))),
                    paste0(group, " ", rep("DCM", nrow(som_maps)))),3)
  bc_mat$Distance <- rep(som_maps$Dist_mean, 6)
  bc_mat$NCDepth <- rep(som_maps$NC_mean, 6)
  bc_mat$Phase <- c(rep("Between", 2*nrow(som_maps)),
                    rep("Warm", 2*nrow(som_maps)),
                    rep("Cool", 2*nrow(som_maps)))
  
  bw_warm_surf <- wilcox.test(som_maps$surf_bc_vals, som_maps$surf_warm_bc_vals, paired = TRUE)
  bw_cool_surf <- wilcox.test(som_maps$surf_bc_vals, som_maps$surf_cool_bc_vals, paired = TRUE)
  bw_warm_deep <- wilcox.test(som_maps$deep_bc_vals, som_maps$deep_warm_bc_vals, paired = TRUE)
  bw_cool_deep <- wilcox.test(som_maps$deep_bc_vals, som_maps$deep_cool_bc_vals, paired = TRUE)
  
  p_val_mat <- as.data.frame(matrix(NA, 4, 5))
  colnames(p_val_mat) <- c("p_val", "Group", "Depth", "Group_Depth", "Phase")
  p_val_mat$p_val <- c(bw_warm_surf$p.value, bw_cool_surf$p.value, bw_warm_deep$p.value, bw_cool_deep$p.value)
  p_val_mat$Group <- group
  p_val_mat$Depth <- c("Surface", "Surface", "DCM", "DCM")
  p_val_mat$Group_Depth <- paste0(group, " ",p_val_mat$Depth)
  p_val_mat$Phase <- c("Warm", "Cool", "Warm", "Cool")
  
  pdf(file = community_diff_fig, width = 12, height = 8)
  print(example_plot)
  dev.off()
  
  ##### Second plot
  
  plot_surf_bw <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = som_maps[complete.cases(som_maps$bc_surf_bw_within),],
               aes_string(x = "long", y = "lat", fill = "bc_surf_bw_within"),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Bray-Curtis Similarity\nBetween / Within (cold)", low = "white", high = "red",
                        limits = c(0.7,1.1), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,1),
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle(group,paste0("Surface\n(2014-2016) vs (2017-2018)"))
  
  plot_deep_bw <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = som_maps[complete.cases(som_maps$bc_surf_bw_within),],
               aes_string(x = "long", y = "lat", fill = "bc_surf_bw_within"),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Bray-Curtis Similarity\nBetween / Within (cold)", low = "white", high = "red",
                        limits = c(0.7,1.1), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,1),
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("DCM\n(2014-2016) vs (2017-2018)")
  
  
  plot_surf_bw + plot_deep_bw + plot_layout(nrow = 2, guides = "collect")
  
  return(list(surf = plot_a, deep = plot_b, bc_mat = bc_mat, p_val_mat = p_val_mat))
  
}


##### Community Time Plots #####

community_comparison <- function(in_file = "output/euks_auto_18sv9_full_data.Rdata",
                                 similar_mat = "output/euks_auto_18sv9_dissimilar.Rdata",
                                 in_map = "output/euks_auto_18sv9_map.Rdata",
                                 out_diff_file = "output/euks_auto_18sv9_diffs.Rdata",
                                 title = "Eukaryotic Phytoplankton",
                                 upwelling_index = "output/upwelling_indicies.Rdata",
                                 index_plot = "figures/euks_auto_18sv9_index_plot.pdf"){
  
  load(in_file)
  load(similar_mat)
  load(in_map)
  
  # find centroids
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  wt_1 <- wt.centroid(x = centroid_df, p = 2)
  wt_2 <- wt.centroid(x = centroid_df, p = 3)
  
  clust1 <- which.max(c(wt_1@coords[1], wt_2@coords[1]))
  clust2 <- which.min(c(wt_1@coords[1], wt_2@coords[1]))
  
  if(clust1 == 1){nearshore_som <- 1}
  if(clust1 == 2){nearshore_som <- 2}
  
  full_dat$Year <- as.numeric(substr(full_dat$Cruise,1,4))
  
  full_dat <- full_dat[-which(as.numeric(substr(full_dat$Sta_ID,2,3)) < 76),]
  
  full_dat$ML_NC <- full_dat$MLD_Sigma - full_dat$NCDepth
  
  early_som_maps <- full_dat %>% 
    filter(Year < 2017 | Year == 2019) %>%
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
              richness = mean(richness, na.rm = TRUE), 
              ML_NC_mean = mean(ML_NC, na.rm = TRUE))
  
  
  centroid_df <- SpatialPointsDataFrame(coords = early_som_maps[,c(6,5)], data = early_som_maps)
  centroid1 <- wt.centroid(x =centroid_df , p = 2)
  centroid2 <- wt.centroid(x =centroid_df , p = 3)
  
  early_dist <- distHaversine(centroid1@coords, centroid2@coords)/100
  
  late_som_maps <- full_dat %>% 
    filter(Year > 2016 & Year < 2019) %>%
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
              richness = mean(richness, na.rm = TRUE), 
              ML_NC_mean = mean(ML_NC, na.rm = TRUE))
  
  centroid_df <- SpatialPointsDataFrame(coords = late_som_maps[,c(6,5)], data = late_som_maps)
  centroid1 <- wt.centroid(x =centroid_df , p = 2)
  centroid2 <- wt.centroid(x =centroid_df , p = 3)
  
  late_dist <- distHaversine(centroid1@coords, centroid2@coords)/1000     
  
  compare <- inner_join(early_som_maps, late_som_maps, by = "Sta_ID")
  
  compare$NC_diff <- compare$NC_mean.x - compare$NC_mean.y
  compare$Temp_diff <- compare$temp_mean.x - compare$temp_mean.y
  compare$diversity_diff <- compare$shannon.x - compare$shannon.y
  compare$even_diff <- compare$evenness.x - compare$evenness.y
  compare$rich_diff <- compare$richness.x - compare$richness.y
  compare$ML_NC_diff <- compare$ML_NC_mean.x - compare$ML_NC_mean.y
  
  
  map <- map_data("world")  
  
  diversity_diff <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = compare, aes(x = long.x, y = lat.x, fill = diversity_diff),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(paste(Delta," Shannon Diversity")),low = "blue", high = "red", mid = "white", midpoint = 0) +
    ggtitle(paste0(title,"\nDifference in Diversity Early-Late")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  print(diversity_diff)
  
  even_diff <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = compare, aes(x = long.x, y = lat.x, fill = even_diff),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(paste(Delta," Evenness")),low = "blue", high = "red", mid = "white", midpoint = 0) +
    ggtitle(paste0(title,"\nDifference in Evenness Early-Late")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  print(even_diff)
  
  som_cruise <- full_dat %>% 
    group_by(Cruise) %>%
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
              richness = mean(richness, na.rm = TRUE), 
              ML_NC_mean = mean(ML_NC, na.rm = TRUE), Date = mean(Date, na.rm = TRUE))
  
  # full_dat$dist_to_coast <- full_dat$dist_to_coast*1000
  
  results <-  full_dat %>% 
    group_by(Cruise) %>%
    do(model = lm(NCDepth ~ dist_to_coast, data = .)) %>%
    mutate(coef=coef(model)["dist_to_coast"])
  
  som_cruise$NC_slope <- results$coef
  
  som_cruise$phase <- c(rep("2014-2016",12),rep("2017-2018",8),rep("2019",4))
  som_cruise$season <- as.factor(rep(c("Winter", "Spring", "Summer", "Fall"),6))
  som_cruise$season <- factor(som_cruise$season, levels = c("Winter", "Spring", "Summer", "Fall"))
  
  gradient_plot <- ggplot(som_cruise, aes_string(x = "NC_slope", y = paste0("som_",nearshore_som))) +
    geom_point(size = 3, aes_string(color = "phase", shape = "season"), data = som_cruise) +
    scale_color_manual(values = c("red", "blue","gold3")) +
    scale_shape_manual(values = c(0,1,2,3)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Season", color = "Phase") + xlab("Nearshore-Offshore\nSlope in Nitracline") +
    ylab("Frequency of Nearshore Cluster") + ggtitle(title)
  
  
  which(!is.na(match(colnames(som_cruise), paste0("som_",nearshore_som))))
  
  formula1 <- as.formula(paste0("som_",nearshore_som, "~",
                   "NC_slope"))
  
  warm_lm <- summary(lm(formula = formula1,
                         data = som_cruise %>% filter(phase == "2014-2016")))
 
  warm_p <- warm_lm$coefficients[2,4]
  warm_rsq <- warm_lm$r.squared
  
  cool_lm <- summary(lm(formula = formula1,
                         data = som_cruise %>% filter(phase == "2017-2018")))
  
  cool_p <- cool_lm$coefficients[2,4]
  cool_rsq <- cool_lm$r.squared
  
  phase <- ggplot(som_cruise, aes_string(x = "NC_slope", y = paste0("som_",nearshore_som))) +
    geom_point(size = 3, aes_string(fill = "phase", color = "phase", shape = "season"), data = som_cruise) +
    stat_smooth(data = som_cruise %>% filter(phase != "2019"), 
                aes_string(x = "NC_slope", y = paste0("som_",nearshore_som), fill = "phase", color = "phase"),
                method="lm") +
    scale_fill_manual(values = c("red", "blue","gold3"), guide = FALSE) +
    scale_color_manual(values = c("red", "blue","gold3")) +
    scale_shape_manual(values = c(21,22,23,24)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Season", color = "Phase") + xlab("Nearshore-Offshore\nSlope in Nitracline") +
    ylab("Frequency of Nearshore Cluster") + ggtitle(title) +
    annotate(geom = "text", y = (max(som_cruise[,paste0("som_",nearshore_som)]) + 0.05), x = 0.15, 
             label = paste0("2014-2016\nR-Squared = ",
                            round(warm_rsq,3),
                            "\np-value = ",
                            round(warm_p, 3)),
             color = "red", size = 3) +
    annotate(geom = "text", y = (min(som_cruise[,paste0("som_",nearshore_som)]) + 0.05), x = 0.2, 
             label = paste0("2017-2018\nR-Squared = ",
                            round(cool_rsq,3),
                            "\np-value = ",
                            round(cool_p, 3)),
             color = "blue", size = 3)
  
  season <- ggplot(som_cruise, aes_string(x = "NC_slope", y = paste0("som_",nearshore_som))) +
    geom_point(size = 3, aes_string(color = "season", fill = "season", shape = "season"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_fill_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    stat_ellipse(aes_string(color = "season")) + 
    scale_shape_manual(values = c(21,22,23,24)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Season", color = "Season", fill = "Season") + xlab("Nearshore-Offshore\nSlope in Nitracline") +
    ylab("Frequency of Nearshore Cluster") + ggtitle(title)
  
  
  
  som_cruise$Year <- substr(som_cruise$Cruise, 1, 4)
  som_cruise$even_odd <- rep(c("even", "even","even","even",
                               "odd", "odd", "odd", "odd"),3)
  
  ts_plot <- ggplot(som_cruise, aes_string(x = "Date", y = paste0("som_",nearshore_som))) +
    geom_rect(aes(xmin = Date, xmax = dplyr::lead(Date), ymin = 0, ymax = 1, fill = even_odd), 
              alpha = 1, show.legend = FALSE) +
    scale_fill_manual(values = c("grey90", "white")) +
    geom_line(color = "black", size = 1) + 
    geom_point(size = 3, aes_string(color = "season"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("Date") +
    ylab("Frequency of Nearshore Cluster")  + ggtitle(title) +
    scale_y_continuous(expand = c(0,0)) 
  
  print(ts_plot)
  
  load(upwelling_index)
  
  index_mat <- matrix(NA,nrow = nrow(som_cruise), 3)
  
  for (i in 1:nrow(som_cruise)) {
    
    year <- substr(som_cruise$Cruise[i],1,4)
    month <- substr(som_cruise$Cruise[i],5,6) 
    
    index_yr <- substr(index_vals_mat$Date,1,4)
    index_month <- substr(index_vals_mat$Date,6,7)
    
    index_vals <- which(match(index_yr,year) & match(index_month,month))
    
    index_mat[i,] <- colMeans(index_vals_mat[index_vals,2:4], na.rm = TRUE)
    
  }
  
  som_cruise$CUTI <- index_mat[,1]
  som_cruise$BEUTI <- index_mat[,2]
  som_cruise$Nitrate <- index_mat[,3]
  
  som_cruise$log_BEUTI <- log(som_cruise$BEUTI)
  som_cruise$log_Nitrate <- log(som_cruise$Nitrate)
  
  cuti_plot <- ggplot(som_cruise, aes_string(x = "CUTI", y = paste0("som_",nearshore_som))) +
    geom_point(size = 3, aes_string(fill = "phase", color = "phase", shape = "season"), data = som_cruise) +
    scale_fill_manual(values = c("red", "blue","gold3"), guide = FALSE) +
    scale_color_manual(values = c("red", "blue","gold3")) +
    stat_ellipse(data = som_cruise %>% filter(phase != "2019"),aes_string(color = "phase")) + 
    scale_shape_manual(values = c(21,22,23,24)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Season", color = "Phase") + xlab("Coastal Upwelling Transport Index\n(CUTI)") +
    ylab("Frequency of Nearshore Cluster")  + ggtitle(title)
  
  
  formula1 <- as.formula(paste0("som_",nearshore_som, "~",
                                "BEUTI"))
  
  warm_lm <- summary(lm(formula = formula1,
                        data = som_cruise %>% filter(phase == "2014-2016")))
  
  warm_p <- warm_lm$coefficients[2,4]
  warm_rsq <- warm_lm$r.squared
  
  cool_lm <- summary(lm(formula = formula1,
                        data = som_cruise %>% filter(phase == "2017-2018")))
  
  cool_p <- cool_lm$coefficients[2,4]
  cool_rsq <- cool_lm$r.squared
  

  
  beuti_plot <- ggplot(som_cruise, aes_string(x = "BEUTI", y = paste0("som_",nearshore_som))) +
    geom_point(size = 3, aes_string(fill = "phase", color = "phase", shape = "season"), data = som_cruise) +
    stat_smooth(data = som_cruise %>% filter(phase != "2019"), 
                aes_string(x = "BEUTI", y = paste0("som_",nearshore_som), fill = "phase", color = "phase"),
                method="lm") +
    scale_fill_manual(values = c("red", "blue","gold3"),  guide = "none") +
    scale_color_manual(values = c("red", "blue","gold3")) +
    scale_shape_manual(values = c(21,22,23,24)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Season", color = "Phase") + xlab("Biologically Effective Upwelling Transport Index\n(BEUTI)") +
    ylab("Frequency of Nearshore Cluster") + ggtitle(title) +
    annotate(geom = "text", y = (max(som_cruise[,paste0("som_",nearshore_som)]) + 0.05), x = 0.15, 
             label = paste0("2014-2016\nR-Squared = ",
                            round(warm_rsq,3),
                            "\np-value = ",
                            round(warm_p, 3)),
             color = "red", size = 3) +
    annotate(geom = "text", y = (min(som_cruise[,paste0("som_",nearshore_som)]) + 0.05), x = 0.2, 
             label = paste0("2017-2018\nR-Squared = ",
                            round(cool_rsq,3),
                            "\np-value = ",
                            round(cool_p, 3)),
             color = "blue", size = 3)
  
  

  log_beuti_plot <- ggplot(som_cruise, aes_string(x = "log_BEUTI", y = paste0("som_",nearshore_som))) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") + 
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17,18)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") + 
    xlab("log(Biologically Effective Upwelling Transport Index)\n(BEUTI)") +
    ylab("") 
  
  reg_nitrate <- ggplot(som_cruise, aes_string(x = "Nitrate", y = paste0("som_",nearshore_som))) +
    geom_point(size = 3, aes_string(fill = "phase", color = "phase", shape = "season"), data = som_cruise) +
    scale_fill_manual(values = c("red", "blue","gold3"), guide = FALSE) +
    scale_color_manual(values = c("red", "blue","gold3")) +
    stat_ellipse(data = som_cruise %>% filter(phase != "2019"),aes_string(color = "phase")) + 
    scale_shape_manual(values = c(21,22,23,24)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Season", color = "Phase") + xlab("Regionally Availible Nitrate") +
    ylab("Frequency of Nearshore Cluster")  + ggtitle(title)
  
  log_reg_nitrate <- ggplot(som_cruise, aes_string(x = "log_Nitrate", y = paste0("som_",nearshore_som))) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") + 
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17,18)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("log(Regionally Availible Nitrate)") +
    ylab("") 
  
  
  title_plot <- ggdraw() + draw_label(title, fontface='bold')
  
  index_plots <- plot_grid(title_plot,
                           plot_grid(cuti_plot, beuti_plot, reg_nitrate, ncol = 3),
                           nrow = 2, rel_heights = c(0.1,1))
  
  pdf(index_plot, width = 12, height = 4)
  print(index_plots)
  dev.off()
  
  diss_melt <- melt(dissimilar)
  diss_melt <- diss_melt[-which(diss_melt$Var1 == diss_melt$Var2),]
  
  diss_melt$Var1_som <- full_dat$som_id[match(diss_melt$Var1,full_dat$eco_name)]
  diss_melt$Var2_som <- full_dat$som_id[match(diss_melt$Var2,full_dat$eco_name)]
  
  diss_melt$Var1_phase <- NA
  diss_melt$Var2_phase <- NA
  
  diss_melt$Var1_phase[which(substr(diss_melt$Var1,2,5) < 2017)] <- "Early"
  diss_melt$Var2_phase[which(substr(diss_melt$Var2,2,5) < 2017)] <- "Early"
  
  diss_melt$Var1_phase[which(substr(diss_melt$Var1,2,5) > 2016 & substr(diss_melt$Var1,2,5) < 2019)] <- "Late"
  diss_melt$Var2_phase[which(substr(diss_melt$Var2,2,5) > 2016 & substr(diss_melt$Var1,2,5) < 2019)] <- "Late"
  
  diss_melt$Var1_phase[which(substr(diss_melt$Var1,2,5) > 2018)] <- "2019"
  diss_melt$Var2_phase[which(substr(diss_melt$Var1,2,5) > 2018)] <- "2019"
  
  diss_melt$comp_phase <- paste0(diss_melt$Var1_phase,"-",diss_melt$Var2_phase)
  diss_melt$comp_som <- paste0(diss_melt$Var1_som,"-",diss_melt$Var2_som)
  
  diss_filt <- filter(diss_melt, comp_som == "1-1" | comp_som == "2-2")
  diss_filt <- filter(diss_filt, comp_phase == "Early-Early" | comp_phase == "Late-Late" | comp_phase == "2019-2019")
  
  diss_filt$comp_phase[diss_filt$comp_phase == "Early-Early"] = "2014-2016"
  diss_filt$comp_phase[diss_filt$comp_phase == "Late-Late"] = "2017-2018"
  diss_filt$comp_phase[diss_filt$comp_phase == "2019-2019"] = "2019"
  
  if(nearshore_som == 1){
    diss_filt$comp_som[diss_filt$comp_som == "1-1"] = "Nearshore"
    diss_filt$comp_som[diss_filt$comp_som == "2-2"] = "Offshore"
  }
  if(nearshore_som == 2){
    diss_filt$comp_som[diss_filt$comp_som == "1-1"] = "Offshore"
    diss_filt$comp_som[diss_filt$comp_som == "2-2"] = "Nearshore"
  }
  
  
  dissimilar_plot <- ggplot(diss_filt,
                            aes(y = 1-value, x = interaction(comp_som,comp_phase),fill = comp_phase)) +
    geom_boxplot() + ylab("Bray-Curtis Similarity") + xlab("Phase/Cluster") +
    labs(fill = "Phase") + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = c("Nearshore\n2014-2016", "Offshore\n2014-2016",
                                "Nearshore\n2017-2018", "Offshore\n2017-2018",
                                "Nearshore\n2019", "Offshore\n2019")) +
    ggtitle(title)
  
  print(dissimilar_plot)
  
  save(early_dist, late_dist, diversity_diff, even_diff, phase, season, ts_plot,
       cuti_plot, beuti_plot, reg_nitrate, dissimilar_plot, file = out_diff_file)
  
  
  
}

##### Line Specific Plots #####

community_comparison_line <- function(in_file = "output/euks_auto_18sv9_full_data.Rdata",
                                 similar_mat = "output/euks_auto_18sv9_dissimilar.Rdata",
                                 in_map = "output/euks_auto_18sv9_map.Rdata",
                                 out_diff_file = "output/euks_auto_18sv9_line_diffs.Rdata",
                                 title = "Eukaryotic Phytoplankton"){
  
  load(in_file)
  load(similar_mat)
  load(in_map)
  
  # find centroids
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  wt_1 <- wt.centroid(x = centroid_df, p = 2)
  wt_2 <- wt.centroid(x = centroid_df, p = 3)
  
  clust1 <- which.max(c(wt_1@coords[1], wt_2@coords[1]))
  clust2 <- which.min(c(wt_1@coords[1], wt_2@coords[1]))
  
  if(clust1 == 1){nearshore_som <- 1}
  if(clust1 == 2){nearshore_som <- 2}
  
  full_dat$Line <- substr(full_dat$Sta_ID,1,3)
  
  line_90 <- full_dat %>%
    filter(Line == "090")
  
  line_80 <- full_dat %>%
    filter(Line == "080")
  
  som_cruise_90 <- line_90 %>% 
    group_by(Cruise) %>%
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
              richness = mean(richness, na.rm = TRUE), Date = mean(Date, na.rm = TRUE))
  
  som_cruise_80 <- line_80 %>% 
    group_by(Cruise) %>%
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
              richness = mean(richness, na.rm = TRUE),  Date = mean(Date, na.rm = TRUE)) 
  
  results_90 <-  line_90 %>% 
    group_by(Cruise) %>%
    do(model = lm(NCDepth ~ dist_to_coast, data = .)) %>%
    mutate(coef=coef(model)["dist_to_coast"])
  
  results_80 <-  line_80 %>% 
    group_by(Cruise) %>%
    do(model = lm(NCDepth ~ dist_to_coast, data = .)) %>%
    mutate(coef=coef(model)["dist_to_coast"])
  
  som_cruise_90$NC_slope <- results_90$coef
  som_cruise_80$NC_slope <- results_80$coef
  
  som_cruise_90$phase <- c(rep("2014-2016",12),rep("2017-2018",8),rep("2019",4))
  som_cruise_90$season <- as.factor(rep(c("Winter", "Spring", "Summer", "Fall"),6))
  som_cruise_90$season <- factor(som_cruise_90$season, levels = c("Winter", "Spring", "Summer", "Fall"))
  
  som_cruise_80$phase <- c(rep("2014-2016",11),rep("2017-2018",8),rep("2019",4))
  som_cruise_80$season <- as.factor(c("Spring", "Summer", "Fall",rep(c("Winter", "Spring", "Summer", "Fall"),5)))
  som_cruise_80$season <- factor(som_cruise_80$season, levels = c("Winter", "Spring", "Summer", "Fall"))
  
  som_cruise_90$Line <- "90"
  som_cruise_80$Line <- "80"
  
  total_cruise <- bind_rows(som_cruise_80, som_cruise_90)
  
  
  ggplot(total_cruise, aes_string(x = "NC_slope", y = paste0("som_",nearshore_som))) +
    geom_point(size = 3, aes_string(fill = "Line", color = "Line", shape = "season"), data = total_cruise) 

  
  
  
  phase_90 <- ggplot(som_cruise_90, aes_string(x = "NC_slope", y = paste0("som_",nearshore_som))) +
    geom_point(size = 3, aes_string(fill = "phase", color = "phase", shape = "season"), data = som_cruise_90) +
    scale_fill_manual(values = c("red", "blue","gold3"), guide = FALSE) +
    scale_color_manual(values = c("red", "blue","gold3")) +
    stat_ellipse(data = som_cruise_90 %>% filter(phase != "2019"),aes_string(color = "phase")) + 
    scale_shape_manual(values = c(21,22,23,24)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Season", color = "Phase") + xlab("Nearshore-Offshore\nSlope in Nitracline") +
    ylab("Frequency of Nearshore Cluster") + ggtitle(title)
  

  

  
 
  save(phase_90, phase_80, file = out_diff_file)
  
  
  
}

##### Line Specific Plots #####

community_comparison_station <- function(in_file = "output/euks_auto_18sv9_full_data.Rdata",
                                      similar_mat = "output/euks_auto_18sv9_dissimilar.Rdata",
                                      in_map = "output/euks_auto_18sv9_map.Rdata",
                                      out_diff_file = "output/euks_auto_18sv9_station_diffs.Rdata",
                                      title = "Eukaryotic Phytoplankton"){
  
  load(in_file)
  load(similar_mat)
  load(in_map)
  
  # find centroids
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  wt_1 <- wt.centroid(x = centroid_df, p = 2)
  wt_2 <- wt.centroid(x = centroid_df, p = 3)
  
  clust1 <- which.max(c(wt_1@coords[1], wt_2@coords[1]))
  clust2 <- which.min(c(wt_1@coords[1], wt_2@coords[1]))
  
  if(clust1 == 1){nearshore_som <- 1}
  if(clust1 == 2){nearshore_som <- 2}
  
  full_dat$Year <- as.numeric(substr(full_dat$Cruise,1,4))
  
  som_cruise_warm <- full_dat %>% 
    filter(Year < 2017) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som_id == 1, na.rm = TRUE)/n(), som_2 = sum(som_id == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE),
              temp_mean = mean(T_degC, na.rm = TRUE), sal_mean = mean(Salnty, na.rm = TRUE),
              PO4_mean = mean(PO4ug, na.rm = TRUE), NO3_mean = mean(NO3ug, na.rm = TRUE),
              SiO3_mean = mean(SiO3ug, na.rm = TRUE), C14_mean = mean(IntC14, na.rm = TRUE),
              MLD_mean = mean(MLD_Sigma, na.rm = TRUE), NC_mean = mean(NCDepth, na.rm = TRUE),
              chl_mean = mean(ChlorA, na.rm = TRUE), chl_coeff = sd(ChlorA, na.rm = TRUE)/mean(ChlorA, na.rm = TRUE),
              temp_coeff = sd(T_degC, na.rm = TRUE)/mean(T_degC, na.rm = TRUE),
              sal_coeff = sd(Salnty, na.rm = TRUE)/mean(Salnty, na.rm = TRUE),
              PO4_coeff =  sd(PO4ug, na.rm = TRUE)/mean(PO4ug, na.rm = TRUE),
              NO3_coeff = sd(NO3ug, na.rm = TRUE)/mean(NO3ug, na.rm = TRUE),
              SiO3_coeff = sd(SiO3ug, na.rm = TRUE)/mean(SiO3ug, na.rm = TRUE),
              C14_coeff = sd(IntC14, na.rm = TRUE)/mean(IntC14, na.rm = TRUE),
              MLD_coeff = sd(MLD_Sigma, na.rm = TRUE)/mean(MLD_Sigma, na.rm = TRUE),
              NC_coeff = sd(NCDepth, na.rm = TRUE)/mean(NCDepth, na.rm = TRUE),
              evenness = mean(evenness, na.rm = TRUE), shannon = mean(shannon, na.rm = TRUE),
              richness = mean(richness, na.rm = TRUE), Date = mean(Date, na.rm = TRUE))
  
  som_cruise_cool <- full_dat %>% 
    filter(Year > 2016, Year < 2019) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som_id == 1, na.rm = TRUE)/n(), som_2 = sum(som_id == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE),
              temp_mean = mean(T_degC, na.rm = TRUE), sal_mean = mean(Salnty, na.rm = TRUE),
              PO4_mean = mean(PO4ug, na.rm = TRUE), NO3_mean = mean(NO3ug, na.rm = TRUE),
              SiO3_mean = mean(SiO3ug, na.rm = TRUE), C14_mean = mean(IntC14, na.rm = TRUE),
              MLD_mean = mean(MLD_Sigma, na.rm = TRUE), NC_mean = mean(NCDepth, na.rm = TRUE),
              chl_mean = mean(ChlorA, na.rm = TRUE), chl_coeff = sd(ChlorA, na.rm = TRUE)/mean(ChlorA, na.rm = TRUE),
              temp_coeff = sd(T_degC, na.rm = TRUE)/mean(T_degC, na.rm = TRUE),
              sal_coeff = sd(Salnty, na.rm = TRUE)/mean(Salnty, na.rm = TRUE),
              PO4_coeff =  sd(PO4ug, na.rm = TRUE)/mean(PO4ug, na.rm = TRUE),
              NO3_coeff = sd(NO3ug, na.rm = TRUE)/mean(NO3ug, na.rm = TRUE),
              SiO3_coeff = sd(SiO3ug, na.rm = TRUE)/mean(SiO3ug, na.rm = TRUE),
              C14_coeff = sd(IntC14, na.rm = TRUE)/mean(IntC14, na.rm = TRUE),
              MLD_coeff = sd(MLD_Sigma, na.rm = TRUE)/mean(MLD_Sigma, na.rm = TRUE),
              NC_coeff = sd(NCDepth, na.rm = TRUE)/mean(NCDepth, na.rm = TRUE),
              evenness = mean(evenness, na.rm = TRUE), shannon = mean(shannon, na.rm = TRUE),
              richness = mean(richness, na.rm = TRUE), Date = mean(Date, na.rm = TRUE))
  
  som_cruise_2019 <- full_dat %>% 
    filter(Year > 2018) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som_id == 1, na.rm = TRUE)/n(), som_2 = sum(som_id == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE),
              temp_mean = mean(T_degC, na.rm = TRUE), sal_mean = mean(Salnty, na.rm = TRUE),
              PO4_mean = mean(PO4ug, na.rm = TRUE), NO3_mean = mean(NO3ug, na.rm = TRUE),
              SiO3_mean = mean(SiO3ug, na.rm = TRUE), C14_mean = mean(IntC14, na.rm = TRUE),
              MLD_mean = mean(MLD_Sigma, na.rm = TRUE), NC_mean = mean(NCDepth, na.rm = TRUE),
              chl_mean = mean(ChlorA, na.rm = TRUE), chl_coeff = sd(ChlorA, na.rm = TRUE)/mean(ChlorA, na.rm = TRUE),
              temp_coeff = sd(T_degC, na.rm = TRUE)/mean(T_degC, na.rm = TRUE),
              sal_coeff = sd(Salnty, na.rm = TRUE)/mean(Salnty, na.rm = TRUE),
              PO4_coeff =  sd(PO4ug, na.rm = TRUE)/mean(PO4ug, na.rm = TRUE),
              NO3_coeff = sd(NO3ug, na.rm = TRUE)/mean(NO3ug, na.rm = TRUE),
              SiO3_coeff = sd(SiO3ug, na.rm = TRUE)/mean(SiO3ug, na.rm = TRUE),
              C14_coeff = sd(IntC14, na.rm = TRUE)/mean(IntC14, na.rm = TRUE),
              MLD_coeff = sd(MLD_Sigma, na.rm = TRUE)/mean(MLD_Sigma, na.rm = TRUE),
              NC_coeff = sd(NCDepth, na.rm = TRUE)/mean(NCDepth, na.rm = TRUE),
              evenness = mean(evenness, na.rm = TRUE), shannon = mean(shannon, na.rm = TRUE),
              richness = mean(richness, na.rm = TRUE), Date = mean(Date, na.rm = TRUE))

  som_cruise_warm$Phase <- "2014-2016"
  som_cruise_cool$Phase <- "2017-2018"
  som_cruise_2019$Phase <- "2019"
  
  all_comb <- bind_rows(som_cruise_warm, som_cruise_cool, som_cruise_2019)
  
  chl_plot <- ggplot(all_comb %>% filter(n_samps > 2, Phase != "2019"), aes_string(x = "chl_mean",
                              y = paste0("som_",nearshore_som),
                              color = "Phase", fill = "Phase")) +
    geom_point() + stat_smooth(method="glm", method.args = list(family = "binomial")) +
    scale_color_manual(values = c("red", "blue")) +  scale_fill_manual(values = c("red", "blue")) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black")) + xlab("Mean Chlorophyll") +
    ylab("Proportion of Samples\nIdentified as Nearshore") + ggtitle(title)
  
  
  
  save(chl_plot, file = out_diff_file)
  
  
  
}

community_comparison_station_2 <- function(in_file = "output/euks_auto_18sv9_full_data.Rdata",
                                         similar_mat = "output/euks_auto_18sv9_dissimilar.Rdata",
                                         in_map = "output/euks_auto_18sv9_map.Rdata",
                                         out_diff_file = "output/euks_auto_18sv9_station2_diffs.Rdata",
                                         title = "Eukaryotic Phytoplankton"){
  
  load(in_file)
  load(similar_mat)
  load(in_map)
  
  # find centroids
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  wt_1 <- wt.centroid(x = centroid_df, p = 2)
  wt_2 <- wt.centroid(x = centroid_df, p = 3)
  
  clust1 <- which.max(c(wt_1@coords[1], wt_2@coords[1]))
  clust2 <- which.min(c(wt_1@coords[1], wt_2@coords[1]))
  
  if(clust1 == 1){nearshore_som <- 1}
  if(clust1 == 2){nearshore_som <- 2}
  
  full_dat$Year <- as.numeric(substr(full_dat$Cruise,1,4))
  full_dat$Line <- as.numeric(substr(full_dat$Sta_ID,1,3))
  full_dat$Station <- substr(full_dat$Sta_ID,7,9)
  full_dat$Cruise <- as.character(full_dat$Cruise)
  
  full_dat <- full_dat %>% filter(Line > 75)
  
  soms <- full_dat$som_id
  
  soms[which(soms != nearshore_som)] <- "Offshore"
  soms[which(soms == nearshore_som)] <- "Nearshore"

  
  full_dat$som_name <- soms
  
  line_80_dat <- full_dat %>%
    filter(Line == 80)
  line_90_dat <- full_dat %>%
    filter(Line == 90)
  
  line_80 <- ggplot(line_80_dat, aes(x = Station, y = Cruise, fill = NCDepth)) +
    geom_tile(width = 1, height = 1) + 
    scale_fill_gradient(low = "aquamarine", high = "darkblue", na.value = "white") +
    geom_point(aes(color = som_name), size = 3) +
    scale_color_manual(values = c("red", "blue")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black")) +
    ggtitle(paste(title,"Line 80")) + 
    labs(color = "Community", fill = "Nitracline Depth (m)")
    
  
  line_90 <- ggplot(line_90_dat, aes(x = Station, y = Cruise, fill = NCDepth)) +
    geom_tile(width = 1, height = 1) + 
    scale_fill_gradient(low = "aquamarine", high = "darkblue", na.value = "white") +
    geom_point(aes(color = som_name), size = 3) +
    scale_color_manual(values = c("red", "blue")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.position = "none",
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) + 
    ggtitle(paste(title,"Line 90")) + 
    labs(color = "Community", fill = "Nitracline Depth (m)") 
  
  combo_plot <- line_80 + line_90 + plot_layout(guides = "collect")
  
  print(combo_plot)
  
  som_cruise_warm <- full_dat %>% 
    filter(Year < 2017) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som_id == 1, na.rm = TRUE)/n(), som_2 = sum(som_id == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE),
              temp_mean = mean(T_degC, na.rm = TRUE), sal_mean = mean(Salnty, na.rm = TRUE),
              PO4_mean = mean(PO4ug, na.rm = TRUE), NO3_mean = mean(NO3ug, na.rm = TRUE),
              SiO3_mean = mean(SiO3ug, na.rm = TRUE), C14_mean = mean(IntC14, na.rm = TRUE),
              MLD_mean = mean(MLD_Sigma, na.rm = TRUE), NC_mean = mean(NCDepth, na.rm = TRUE),
              chl_mean = mean(ChlorA, na.rm = TRUE), chl_coeff = sd(ChlorA, na.rm = TRUE)/mean(ChlorA, na.rm = TRUE),
              temp_coeff = sd(T_degC, na.rm = TRUE)/mean(T_degC, na.rm = TRUE),
              sal_coeff = sd(Salnty, na.rm = TRUE)/mean(Salnty, na.rm = TRUE),
              PO4_coeff =  sd(PO4ug, na.rm = TRUE)/mean(PO4ug, na.rm = TRUE),
              NO3_coeff = sd(NO3ug, na.rm = TRUE)/mean(NO3ug, na.rm = TRUE),
              SiO3_coeff = sd(SiO3ug, na.rm = TRUE)/mean(SiO3ug, na.rm = TRUE),
              C14_coeff = sd(IntC14, na.rm = TRUE)/mean(IntC14, na.rm = TRUE),
              MLD_coeff = sd(MLD_Sigma, na.rm = TRUE)/mean(MLD_Sigma, na.rm = TRUE),
              NC_coeff = sd(NCDepth, na.rm = TRUE)/mean(NCDepth, na.rm = TRUE),
              evenness = mean(evenness, na.rm = TRUE), shannon = mean(shannon, na.rm = TRUE),
              richness = mean(richness, na.rm = TRUE), Date = mean(Date, na.rm = TRUE))
  
  som_cruise_cool <- full_dat %>% 
    filter(Year > 2016, Year < 2019) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som_id == 1, na.rm = TRUE)/n(), som_2 = sum(som_id == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE),
              temp_mean = mean(T_degC, na.rm = TRUE), sal_mean = mean(Salnty, na.rm = TRUE),
              PO4_mean = mean(PO4ug, na.rm = TRUE), NO3_mean = mean(NO3ug, na.rm = TRUE),
              SiO3_mean = mean(SiO3ug, na.rm = TRUE), C14_mean = mean(IntC14, na.rm = TRUE),
              MLD_mean = mean(MLD_Sigma, na.rm = TRUE), NC_mean = mean(NCDepth, na.rm = TRUE),
              chl_mean = mean(ChlorA, na.rm = TRUE), chl_coeff = sd(ChlorA, na.rm = TRUE)/mean(ChlorA, na.rm = TRUE),
              temp_coeff = sd(T_degC, na.rm = TRUE)/mean(T_degC, na.rm = TRUE),
              sal_coeff = sd(Salnty, na.rm = TRUE)/mean(Salnty, na.rm = TRUE),
              PO4_coeff =  sd(PO4ug, na.rm = TRUE)/mean(PO4ug, na.rm = TRUE),
              NO3_coeff = sd(NO3ug, na.rm = TRUE)/mean(NO3ug, na.rm = TRUE),
              SiO3_coeff = sd(SiO3ug, na.rm = TRUE)/mean(SiO3ug, na.rm = TRUE),
              C14_coeff = sd(IntC14, na.rm = TRUE)/mean(IntC14, na.rm = TRUE),
              MLD_coeff = sd(MLD_Sigma, na.rm = TRUE)/mean(MLD_Sigma, na.rm = TRUE),
              NC_coeff = sd(NCDepth, na.rm = TRUE)/mean(NCDepth, na.rm = TRUE),
              evenness = mean(evenness, na.rm = TRUE), shannon = mean(shannon, na.rm = TRUE),
              richness = mean(richness, na.rm = TRUE), Date = mean(Date, na.rm = TRUE))
  
  som_cruise_2019 <- full_dat %>% 
    filter(Year > 2018) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som_id == 1, na.rm = TRUE)/n(), som_2 = sum(som_id == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE),
              temp_mean = mean(T_degC, na.rm = TRUE), sal_mean = mean(Salnty, na.rm = TRUE),
              PO4_mean = mean(PO4ug, na.rm = TRUE), NO3_mean = mean(NO3ug, na.rm = TRUE),
              SiO3_mean = mean(SiO3ug, na.rm = TRUE), C14_mean = mean(IntC14, na.rm = TRUE),
              MLD_mean = mean(MLD_Sigma, na.rm = TRUE), NC_mean = mean(NCDepth, na.rm = TRUE),
              chl_mean = mean(ChlorA, na.rm = TRUE), chl_coeff = sd(ChlorA, na.rm = TRUE)/mean(ChlorA, na.rm = TRUE),
              temp_coeff = sd(T_degC, na.rm = TRUE)/mean(T_degC, na.rm = TRUE),
              sal_coeff = sd(Salnty, na.rm = TRUE)/mean(Salnty, na.rm = TRUE),
              PO4_coeff =  sd(PO4ug, na.rm = TRUE)/mean(PO4ug, na.rm = TRUE),
              NO3_coeff = sd(NO3ug, na.rm = TRUE)/mean(NO3ug, na.rm = TRUE),
              SiO3_coeff = sd(SiO3ug, na.rm = TRUE)/mean(SiO3ug, na.rm = TRUE),
              C14_coeff = sd(IntC14, na.rm = TRUE)/mean(IntC14, na.rm = TRUE),
              MLD_coeff = sd(MLD_Sigma, na.rm = TRUE)/mean(MLD_Sigma, na.rm = TRUE),
              NC_coeff = sd(NCDepth, na.rm = TRUE)/mean(NCDepth, na.rm = TRUE),
              evenness = mean(evenness, na.rm = TRUE), shannon = mean(shannon, na.rm = TRUE),
              richness = mean(richness, na.rm = TRUE), Date = mean(Date, na.rm = TRUE))
  
  som_cruise_warm$Phase <- "2014-2016"
  som_cruise_cool$Phase <- "2017-2018"
  som_cruise_2019$Phase <- "2019"
  
  all_comb <- bind_rows(som_cruise_warm, som_cruise_cool, som_cruise_2019)
  
  chl_plot <- ggplot(all_comb %>% filter(n_samps > 2),
                             aes_string(x = "NC_mean",
                                        y = paste0("som_",nearshore_som),
                                        color = "Phase",
                                        shape = "Phase",
                                        fill = "chl_mean")) + 
    geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
    geom_point(size = 5, stroke = 2) + xlab("Mean Nitracline Depth (m)") +
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    scale_fill_gradient(low = "green", high = "darkgreen") +
    scale_shape_manual(values = c(21,22,24)) + 
    scale_color_manual(values = c("red", "blue","gold3")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black")) +
    ggtitle(title) + labs(fill = "Mean Chl-a")
  
  print(chl_plot)
  
  save(combo_plot, chl_plot, file = out_diff_file)
  
  
  
}


##### Diversity Time Plots #####

diversity_comparison <- function(in_file = "output/euks_auto_18sv9_full_data.Rdata",
                                 in_map = "output/euks_auto_18sv9_map.Rdata",
                                 in_raw = "data/18s_autotrophic_euks.Rdata",
                                 out_diff_file = "output/euks_auto_18sv9_diffs_div.Rdata",
                                 title = "Eukaryotic Phytoplankton",
                                 upwelling_index = "output/upwelling_indicies.Rdata",
                                 index_plot = "figures/euks_auto_18sv9_index_plot_div.pdf"){
  
  load(in_file)
  load(in_raw)
  
  full_dat$Year <- as.numeric(substr(full_dat$Cruise,1,4))
  
  full_dat <- full_dat[-which(as.numeric(substr(full_dat$Sta_ID,2,3)) < 76),]
  
  full_dat$ML_NC <- full_dat$MLD_Sigma - full_dat$NCDepth
  
  early_som_maps <- full_dat %>% 
    filter(Year < 2017 | Year == 2019) %>%
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
              richness = mean(richness, na.rm = TRUE), 
              ML_NC_mean = mean(ML_NC, na.rm = TRUE))
  
  
  centroid_df <- SpatialPointsDataFrame(coords = early_som_maps[,c(6,5)], data = early_som_maps)
  centroid1 <- wt.centroid(x =centroid_df , p = 2)
  centroid2 <- wt.centroid(x =centroid_df , p = 3)
  
  early_dist <- distHaversine(centroid1@coords, centroid2@coords)/100
  
  late_som_maps <- full_dat %>% 
    filter(Year > 2016 & Year < 2019) %>%
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
              richness = mean(richness, na.rm = TRUE), 
              ML_NC_mean = mean(ML_NC, na.rm = TRUE))
  
  centroid_df <- SpatialPointsDataFrame(coords = late_som_maps[,c(6,5)], data = late_som_maps)
  centroid1 <- wt.centroid(x =centroid_df , p = 2)
  centroid2 <- wt.centroid(x =centroid_df , p = 3)
  
  late_dist <- distHaversine(centroid1@coords, centroid2@coords)/1000     
  
  compare <- inner_join(early_som_maps, late_som_maps, by = "Sta_ID")
  
  compare$NC_diff <- compare$NC_mean.x - compare$NC_mean.y
  compare$Temp_diff <- compare$temp_mean.x - compare$temp_mean.y
  compare$diversity_diff <- compare$shannon.x - compare$shannon.y
  compare$even_diff <- compare$evenness.x - compare$evenness.y
  compare$rich_diff <- compare$richness.x - compare$richness.y
  compare$ML_NC_diff <- compare$ML_NC_mean.x - compare$ML_NC_mean.y
  
  
  map <- map_data("world")  
  
  diversity_diff <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = compare, aes(x = long.x, y = lat.x, fill = diversity_diff),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(paste(Delta," Shannon Diversity")),low = "blue", high = "red", mid = "white", midpoint = 0) +
    ggtitle(paste0(title,"\nDifference in Diversity Early-Late")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  print(diversity_diff)
  
  even_diff <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = compare, aes(x = long.x, y = lat.x, fill = even_diff),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(paste(Delta," Evenness")),low = "blue", high = "red", mid = "white", midpoint = 0) +
    ggtitle(paste0(title,"\nDifference in Evenness Early-Late")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  print(even_diff)
  
  som_cruise <- full_dat %>% 
    group_by(Cruise) %>%
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
              richness = mean(richness, na.rm = TRUE), 
              ML_NC_mean = mean(ML_NC, na.rm = TRUE), Date = mean(Date, na.rm = TRUE))
  
  # full_dat$dist_to_coast <- full_dat$dist_to_coast*1000
  
  results <-  full_dat %>% 
    group_by(Cruise) %>%
    do(model = lm(NCDepth ~ dist_to_coast, data = .)) %>%
    mutate(coef=coef(model)["dist_to_coast"])
  
  som_cruise$NC_slope <- results$coef
  
  som_cruise$phase <- c(rep("2014-2016",12),rep("2017-2018",8), rep("2019",4))
  som_cruise$season <- as.factor(rep(c("Winter", "Spring", "Summer", "Fall"),6))
  som_cruise$season <- factor(som_cruise$season, levels = c("Winter", "Spring", "Summer", "Fall"))
  
  # total diversity stats
  
  asv_table$cruise <- substr(rownames(asv_table),2,7)
  
  asv_cruise <- asv_table %>%
    group_by(cruise) %>% summarise(across(.cols = everything(), .fns = ~sum(.x,na.rm = TRUE)))
  
  asv_cruise$total_div <- diversity(asv_cruise[,-1], MARGIN = 1, index = "shannon")
  
  som_cruise$total_shannon <- asv_cruise$total_div[match(som_cruise$Cruise, asv_cruise$cruise)]
  
  formula1 <- as.formula(paste0("shannon", "~",
                                "NC_slope"))
  
  warm_lm <- summary(lm(formula = formula1,
                        data = som_cruise %>% filter(phase == "2014-2016")))
  
  warm_p <- warm_lm$coefficients[2,4]
  warm_rsq <- warm_lm$r.squared
  
  cool_lm <- summary(lm(formula = formula1,
                        data = som_cruise %>% filter(phase == "2017-2018")))
  
  cool_p <- cool_lm$coefficients[2,4]
  cool_rsq <- cool_lm$r.squared
  
  phase <- ggplot(som_cruise, aes_string(x = "NC_slope", y = "shannon")) +
    geom_point(size = 3, aes_string(fill = "phase", color = "phase", shape = "season"), data = som_cruise) +
    stat_smooth(data = som_cruise %>% filter(phase != "2019"), 
                aes_string(x = "NC_slope", y = "shannon", fill = "phase", color = "phase"),
                method="lm") +
    scale_fill_manual(values = c("red", "blue","gold3"), guide = FALSE) +
    scale_color_manual(values = c("red", "blue","gold3")) +
    scale_shape_manual(values = c(21,22,23,24)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Season", color = "Phase") + xlab("Nearshore-Offshore\nSlope in Nitracline") +
    ylab("Mean Alpha Diversity\nPer Cruise") + ggtitle(title) +
    annotate(geom = "text", y = (max(som_cruise[,"shannon"]) + 0.05), x = 0.16, 
             label = paste0("2014-2016\nR-Squared = ",
                            round(warm_rsq,3),
                            "\np-value = ",
                            round(warm_p, 3)),
             color = "red", size = 3) +
    annotate(geom = "text", y = (min(som_cruise[,"shannon"]) + 0.05), x = 0.15, 
             label = paste0("2017-2018\nR-Squared = ",
                            round(cool_rsq,3),
                            "\np-value = ",
                            round(cool_p, 3)),
             color = "blue", size = 3)
  
  formula2 <- as.formula(paste0("total_shannon", "~",
                                "NC_slope"))
  
  warm_lm <- summary(lm(formula = formula2,
                        data = som_cruise %>% filter(phase == "2014-2016")))
  
  warm_p <- warm_lm$coefficients[2,4]
  warm_rsq <- warm_lm$r.squared
  
  cool_lm <- summary(lm(formula = formula2,
                        data = som_cruise %>% filter(phase == "2017-2018")))
  
  cool_p <- cool_lm$coefficients[2,4]
  cool_rsq <- cool_lm$r.squared
  
  gradient_plot2 <- ggplot(som_cruise, aes_string(x = "NC_slope", y = "total_shannon")) +
    geom_point(size = 3, aes_string(fill = "phase", color = "phase", shape = "season"), data = som_cruise) +
    stat_smooth(data = som_cruise %>% filter(phase != "2019"), 
                aes_string(x = "NC_slope", y = "total_shannon", fill = "phase", color = "phase"),
                method="lm") +
    scale_fill_manual(values = c("red", "blue","gold3"), guide = FALSE) +
    scale_color_manual(values = c("red", "blue","gold3")) +
    scale_shape_manual(values = c(21,22,23,24)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Season", color = "Phase") + xlab("Nearshore-Offshore\nSlope in Nitracline") +
    ylab("Gamma Diversity\nPer Cruise") + ggtitle(title) +
    annotate(geom = "text", y = (max(som_cruise[,"total_shannon"]) + 0.05), x = 0.16, 
             label = paste0("2014-2016\nR-Squared = ",
                            round(warm_rsq,3),
                            "\np-value = ",
                            round(warm_p, 3)),
             color = "red", size = 3) +
    annotate(geom = "text", y = (min(som_cruise[,"total_shannon"]) + 0.05), x = 0.15, 
             label = paste0("2017-2018\nR-Squared = ",
                            round(cool_rsq,3),
                            "\np-value = ",
                            round(cool_p, 3)),
             color = "blue", size = 3)
  
  print(gradient_plot2)
  
  ts_plot <- ggplot(som_cruise, aes_string(x = "Date", y = "shannon")) +
    geom_line(color = "black", size = 1) +
    geom_point(size = 3, aes_string(color = "season"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("Date") +
    ylab("Mean Shannon Diversity\nPer Cruise")  + ggtitle(title)
  
  print(ts_plot)
  
  ts_plot2 <- ggplot(som_cruise, aes_string(x = "Date", y = "total_shannon")) +
    geom_line(color = "black", size = 1) +
    geom_point(size = 3, aes_string(color = "season"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("Date") +
    ylab("Mean Shannon Diversity\nPer Cruise")  + ggtitle(title)
  
  print(ts_plot2)
  
  load(upwelling_index)
  
  index_mat <- matrix(NA,nrow = nrow(som_cruise), 3)
  
  for (i in 1:nrow(som_cruise)) {
    
    year <- substr(som_cruise$Cruise[i],1,4)
    month <- substr(som_cruise$Cruise[i],5,6)
    
    index_yr <- substr(index_vals_mat$Date,1,4)
    index_month <- substr(index_vals_mat$Date,6,7)
    
    index_vals <- which(match(index_yr,year) & match(index_month,month))
    
    index_mat[i,] <- colMeans(index_vals_mat[index_vals,2:4], na.rm = TRUE)
    
  }
  
  som_cruise$CUTI <- index_mat[,1]
  som_cruise$BEUTI <- index_mat[,2]
  som_cruise$Nitrate <- index_mat[,3]
  
  som_cruise$log_BEUTI <- log(som_cruise$BEUTI)
  som_cruise$log_Nitrate <- log(som_cruise$Nitrate)
  
  cuti_plot <- ggplot(som_cruise, aes_string(x = "CUTI", y = "shannon")) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") +
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17,18)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("Coastal Upwelling Transport Index\n(CUTI)") +
    ylab("Mean Shannon Diversity\nPer Cruise")
  
  beuti_plot <- ggplot(som_cruise, aes_string(x = "BEUTI", y = "shannon")) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") +
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17,18)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") +
    xlab("Biologically Effective Upwelling Transport Index\n(BEUTI)") +
    ylab("")
  
  log_beuti_plot <- ggplot(som_cruise, aes_string(x = "log_BEUTI", y = "shannon")) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") +
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17,18)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") +
    xlab("log(Biologically Effective Upwelling Transport Index)\n(BEUTI)") +
    ylab("")
  
  reg_nitrate <- ggplot(som_cruise, aes_string(x = "Nitrate", y = "shannon")) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") +
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17,18)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("Regionally Availible Nitrate") +
    ylab("")
  
  log_reg_nitrate <- ggplot(som_cruise, aes_string(x = "log_Nitrate", y = "shannon")) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") +
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17,18)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("log(Regionally Availible Nitrate)") +
    ylab("")
  
  
  title_plot <- ggdraw() + draw_label(title, fontface='bold')
  
  index_plots <- plot_grid(title_plot,
                           plot_grid(cuti_plot, beuti_plot, reg_nitrate, ncol = 3),
                           nrow = 2, rel_heights = c(0.1,1))
  
  pdf(index_plot, width = 12, height = 4)
  print(index_plots)
  dev.off()
  
  save(early_dist, late_dist, diversity_diff, even_diff, 
       phase, season, ts_plot, gradient_plot2, ts_plot2,
       cuti_plot, beuti_plot, reg_nitrate, file = out_diff_file)
  
}


#### Extra ####

fig_x_func <- function(in_all = "output/total_dissimilar.Rdata",
                       in_dat = "output/total_full_data.Rdata",
                       in_map = "output/total_map.Rdata",
                       community_diff_fig = "figures/figure_outline/fig_x.pdf",
                       tsize = 12, psize = 12){
  
  load(in_all)
  load(in_dat)
  load(in_map)
  
  samps <- rownames(dissimilar)
  mat <- matrix(NA, 985, 5)
  mat[,1] <- as.numeric(substr(samps,2,5))
  mat[,2] <- substr(samps,6,7)
  mat[,4] <- match(samps,full_dat$eco_name)
  
  mat[which(as.numeric(mat[,2]) < 3),3] <- "Winter" 
  mat[which(as.numeric(mat[,2]) == 4),3] <- "Spring" 
  mat[which(as.numeric(mat[,2]) > 5 & as.numeric(mat[,2]) < 9),3] <- "Summer"
  mat[which(as.numeric(mat[,2]) > 8),3] <- "Fall" 
  
  mat[which(as.numeric(mat[,1]) < 2017),5] <- "Early"
  mat[which(as.numeric(mat[,1]) > 2016 & as.numeric(mat[,1]) < 2019),5] <- "Late"
  
  # wint_comp <- as.data.frame(dissimilar[which(mat[,1] < 2017 & mat[,3] == "Winter" & !is.na(mat[,4])),
  #                         which(mat[,1] > 2016 & mat[,1] < 2019 & mat[,3] == "Winter"& !is.na(mat[,4]))])
  # wint_comp$early_samps <- rownames(wint_comp)
  # 
  # spr_comp <- as.data.frame(dissimilar[which(mat[,1] < 2017 & mat[,3] == "Spring" & !is.na(mat[,4])),
  #                         which(mat[,1] > 2016 & mat[,1] < 2019 & mat[,3] == "Spring"& !is.na(mat[,4]))])
  # spr_comp$early_samps <- rownames(spr_comp)
  # 
  # sum_comp <- as.data.frame(dissimilar[which(mat[,1] < 2017 & mat[,3] == "Summer" & !is.na(mat[,4])),
  #                         which(mat[,1] > 2016 & mat[,1] < 2019 & mat[,3] == "Summer"& !is.na(mat[,4]))])
  # sum_comp$early_samps <- rownames(sum_comp)
  # 
  # fall_comp <- as.data.frame(dissimilar[which(mat[,1] < 2017 & mat[,3] == "Fall" & !is.na(mat[,4])),
  #                         which(mat[,1] > 2016 & mat[,1] < 2019 & mat[,3] == "Fall"& !is.na(mat[,4]))])
  # fall_comp$early_samps <- rownames(fall_comp)
  
  all_comp <- as.data.frame(dissimilar)
  all_comp$early_samps <- rownames(all_comp)
  
  wint_long <- wint_comp %>%
    pivot_longer(-early_samps,
                 names_to = "late_samps",
                 values_to = "b_c_vals")
  
  spr_long <- spr_comp %>%
    pivot_longer(-early_samps,
                 names_to = "late_samps",
                 values_to = "b_c_vals")
  
  sum_long <- sum_comp %>%
    pivot_longer(-early_samps,
                 names_to = "late_samps",
                 values_to = "b_c_vals")
  
  fall_long <- fall_comp %>%
    pivot_longer(-early_samps,
                 names_to = "late_samps",
                 values_to = "b_c_vals")
  
  all_long <- all_comp %>%
    pivot_longer(-early_samps,
                 names_to = "late_samps",
                 values_to = "b_c_vals")
  
  # wint_long <- wint_long[which(substr(wint_long$early_samps,9,19)==substr(wint_long$late_samps,9,19)),]
  # spr_long <- spr_long[which(substr(spr_long$early_samps,9,19)==substr(spr_long$late_samps,9,19)),]
  # sum_long <- sum_long[which(substr(sum_long$early_samps,9,19)==substr(sum_long$late_samps,9,19)),]
  # fall_long <- fall_long[which(substr(fall_long$early_samps,9,19)==substr(fall_long$late_samps,9,19)),]
  
  all_long <- all_long[which(substr(all_long$early_samps,9,19)==substr(all_long$late_samps,9,19)),]
  
  # wint_long$Sta_ID <- paste0(substr(wint_long$early_samps,9,13)," ",substr(wint_long$early_samps,15,19))
  # spr_long$Sta_ID <- paste0(substr(spr_long$early_samps,9,13)," ",substr(spr_long$early_samps,15,19))
  # sum_long$Sta_ID <- paste0(substr(sum_long$early_samps,9,13)," ",substr(sum_long$early_samps,15,19))
  # fall_long$Sta_ID <- paste0(substr(fall_long$early_samps,9,13)," ",substr(fall_long$early_samps,15,19))
  all_long$Year_ID <- paste0(substr(all_long$early_samps,2,5),"-",substr(all_long$late_samps,2,5))
  all_long$Sta_ID <- paste0(substr(all_long$early_samps,9,13)," ",substr(all_long$early_samps,15,19))
  all_long <- all_long[-which(is.na(match(all_long$Sta_ID, som_maps$Sta_ID[which(som_maps$n_samps > 30)]))),]
  
  all_long$Date_Diff <- full_dat$Date[match(all_long$early_samps, full_dat$eco_name)] -
    full_dat$Date[match(all_long$late_samps, full_dat$eco_name)]
  
  
  all_long <- all_long[which(as.numeric(substr(all_long$early_samps,2,7))==201404 | as.numeric(substr(all_long$early_samps,2,7))==201402),]
  
  all_long <- all_long[-which(as.numeric(all_long$Date_Diff) > 0),]
  
  all_long$Time_ID <- paste0(substr(all_long$early_samps,2,7),"-",substr(all_long$late_samps,2,7))
  
  all_long$Time_Sta <- substr(all_long$late_samps,2,19)
  
  all_long$Date_Diff <- abs(as.numeric(all_long$Date_Diff))
  all_long$Date_Diff[which(substr(all_long$early_samps,2,7) == "201404")] <- all_long$Date_Diff[which(substr(all_long$early_samps,2,7) == "201404")] + 62
  
  # wint_mean <- wint_long %>%
  #   group_by(Sta_ID) %>%
  #   summarise(mean_bc = mean(b_c_vals, na.rm = TRUE))
  # 
  # spr_mean <- spr_long %>%
  #   group_by(Sta_ID) %>%
  #   summarise(mean_bc = mean(b_c_vals, na.rm = TRUE))
  # 
  # sum_mean <- sum_long %>%
  #   group_by(Sta_ID) %>%
  #   summarise(mean_bc = mean(b_c_vals, na.rm = TRUE))
  # 
  # fall_mean <- fall_long %>%
  #   group_by(Sta_ID) %>%
  #   summarise(mean_bc = mean(b_c_vals, na.rm = TRUE))
  
  all_long <- all_long %>%
    group_by(Time_Sta) %>%
    summarise(mean_bc = mean(b_c_vals, na.rm = TRUE),
              Time_ID = first(Time_ID),
              Sta_ID = first(Sta_ID),
              Time_Diff = first(Date_Diff))
  
  # som_maps$wint_bc <- 1-wint_mean$mean_bc[match(som_maps$Sta_ID, wint_mean$Sta_ID)]
  # som_maps$spr_bc <- 1-spr_mean$mean_bc[match(som_maps$Sta_ID, spr_mean$Sta_ID)]
  # som_maps$sum_bc <- 1-sum_mean$mean_bc[match(som_maps$Sta_ID, sum_mean$Sta_ID)]
  # som_maps$fall_bc <- 1-fall_mean$mean_bc[match(som_maps$Sta_ID, fall_mean$Sta_ID)]
  all_long$b_c_vals <- 1-all_long$mean_bc
  
  start_date <- mdy("01-29-14")
  
  all_long$Date <- start_date + all_long$Time_Diff
  
  red_pal <- colorRampPalette(c("firebrick1","firebrick4"))
  blue_pal <- colorRampPalette(c("dodgerblue1","dodgerblue4"))
  
  col1 <- red_pal(4)
  col2 <- blue_pal(5)
  
  all_long$Line <- substr(all_long$Sta_ID,2,3)
  
  all_long$Sta_ID <- as.factor(all_long$Sta_ID)
  
  plot_90 <- ggplot(all_long %>% filter(Line == "90"), aes(x = Date, y = b_c_vals, color = Sta_ID)) + 
    geom_line(size = 1) + 
    scale_color_manual(values = c(col1, "forestgreen", col2), drop = FALSE) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.text.x = element_text(angle = -90, vjust = 0.5),
          plot.title = element_text(size = tsize),
          axis.text = element_text(size = tsize),
          axis.title = element_text(size = tsize)) +
    xlab("Time") + ylab("Bray-Curtis Similarity") + labs(color = "Station") + ggtitle("Line 90")
  
  plot_80 <- ggplot(all_long %>% filter(Line == "80"), aes(x = Date, y = b_c_vals, color = Sta_ID)) + 
    geom_line(size = 1) + 
    scale_color_manual(values = c(col1, "forestgreen", col2), drop = FALSE) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.text.x = element_text(angle = -90, vjust = 0.5),
          plot.title = element_text(size = tsize),
          axis.text = element_text(size = tsize),
          axis.title = element_text(size = tsize)) +
    xlab("Time") + ylab("Bray-Curtis Similarity") + labs(color = "Station") + ggtitle("Line 80")
  
  plot_81 <- ggplot(all_long %>% filter(Line == "81"), aes(x = Date, y = b_c_vals, color = Sta_ID)) + 
    geom_line(size = 1) + 
    scale_color_manual(values = c(col1, "forestgreen", col2), drop = FALSE) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.text.x = element_text(angle = -90, vjust = 0.5),
          plot.title = element_text(size = tsize),
          axis.text = element_text(size = tsize),
          axis.title = element_text(size = tsize)) +
    xlab("Time") + ylab("Bray-Curtis Similarity") + labs(color = "Station") + ggtitle("Line 81")
  
  example_plot <- plot_80 + plot_81 + plot_90 + guide_area() + plot_layout(guides = "collect")
  
  
  # map <- map_data("world")    
  # 
  # wint <- ggplot() + 
  #   geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  #   coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  #   xlab("Longitude") + ylab("Latitude") + 
  #   geom_point(data = som_maps[complete.cases(som_maps$wint_bc),], aes_string(x = "long", y = "lat", fill = "wint_bc"),
  #              color = "black", size =psize, stroke = 0.1, shape = 21) +
  #   scale_fill_gradient(name = "Bray-Curtis\nSimilarity", low = "white", high = "red", 
  #                       limits = c(0.2,0.4), oob = scales::squish) +
  #   theme(panel.background = element_blank(),
  #         panel.border = element_rect(fill = NA, colour = "black", linetype = "solid"),
  #         plot.title = element_text(), axis.line = element_blank(),
  #         legend.justification=c(1,1), 
  #         legend.position=c(0.97, 0.97),
  #         legend.background = element_rect(fill = "white", color = "black"),
  #         axis.text = element_text(size = tsize),
  #         legend.text = element_text(size = tsize),
  #         axis.title = element_text(size = tsize),
  #         legend.title = element_text(size = tsize)) +
  #   ggtitle("A. Winter Bray-Curtis Similarity")
  # 
  # spr <- ggplot() + 
  #   geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  #   coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  #   xlab("Longitude") + ylab("Latitude") + 
  #   geom_point(data = som_maps[complete.cases(som_maps$spr_bc),], aes_string(x = "long", y = "lat", fill = "spr_bc"),
  #              color = "black", size =psize, stroke = 0.1, shape = 21) +
  #   scale_fill_gradient(name = "Bray-Curtis\nSimilarity", low = "white", high = "red", 
  #                       limits = c(0.2,0.4), oob = scales::squish) +
  #   theme(panel.background = element_blank(),
  #         panel.border = element_rect(fill = NA, colour = "black", linetype = "solid"),
  #         plot.title = element_text(), axis.line = element_blank(),
  #         legend.justification=c(1,1), 
  #         legend.position=c(0.97, 0.97),
  #         legend.background = element_rect(fill = "white", color = "black"),
  #         axis.text = element_text(size = tsize),
  #         legend.text = element_text(size = tsize),
  #         axis.title = element_text(size = tsize),
  #         legend.title = element_text(size = tsize)) +
  #   ggtitle("B. Spring Bray-Curtis Similarity")
  # 
  # summ <- ggplot() + 
  #   geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  #   coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  #   xlab("Longitude") + ylab("Latitude") + 
  #   geom_point(data = som_maps[complete.cases(som_maps$sum_bc),], aes_string(x = "long", y = "lat", fill = "sum_bc"),
  #              color = "black", size =psize, stroke = 0.1, shape = 21) +
  #   scale_fill_gradient(name = "Bray-Curtis\nSimilarity", low = "white", high = "red",
  #                       limits = c(0.2,0.4), oob = scales::squish) +
  #   theme(panel.background = element_blank(),
  #         panel.border = element_rect(fill = NA, colour = "black", linetype = "solid"),
  #         plot.title = element_text(), axis.line = element_blank(),
  #         legend.justification=c(1,1), 
  #         legend.position=c(0.97, 0.97),
  #         legend.background = element_rect(fill = "white", color = "black"),
  #         axis.text = element_text(size = tsize),
  #         legend.text = element_text(size = tsize),
  #         axis.title = element_text(size = tsize),
  #         legend.title = element_text(size = tsize)) +
  #   ggtitle("C. Summer Bray-Curtis Similarity")
  # 
  # 
  # fall <- ggplot() + 
  #   geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  #   coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  #   xlab("Longitude") + ylab("Latitude") + 
  #   geom_point(data = som_maps[complete.cases(som_maps$fall_bc),], aes_string(x = "long", y = "lat", fill = "fall_bc"),
  #              color = "black", size =psize, stroke = 0.1, shape = 21) +
  #   scale_fill_gradient(name = "Bray-Curtis\nSimilarity", low = "white", high = "red",
  #                       limits = c(0.2,0.4), oob = scales::squish) +
  #   theme(panel.background = element_blank(),
  #         panel.border = element_rect(fill = NA, colour = "black", linetype = "solid"),
  #         plot.title = element_text(), axis.line = element_blank(),
  #         legend.justification=c(1,1), 
  #         legend.position=c(0.97, 0.97),
  #         legend.background = element_rect(fill = "white", color = "black"),
  #         axis.text = element_text(size = tsize),
  #         legend.text = element_text(size = tsize),
  #         axis.title = element_text(size = tsize),
  #         legend.title = element_text(size = tsize)) +
  #   ggtitle("D. Fall Bray-Curtis Similarity")
  # 
  # all <- ggplot() + 
  #   geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  #   coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  #   xlab("Longitude") + ylab("Latitude") + 
  #   geom_point(data = som_maps[complete.cases(som_maps$all_bc),], aes_string(x = "long", y = "lat", fill = "all_bc"),
  #              color = "black", size =psize, stroke = 0.1, shape = 21) +
  #   scale_fill_gradient(name = "Bray-Curtis\nSimilarity", low = "white", high = "red",
  #                       limits = c(0.2,0.4), oob = scales::squish) +
  #   theme(panel.background = element_blank(),
  #         panel.border = element_rect(fill = NA, colour = "black", linetype = "solid"),
  #         plot.title = element_text(), axis.line = element_blank(),
  #         legend.justification=c(1,1), 
  #         legend.position=c(0.97, 0.97),
  #         legend.background = element_rect(fill = "white", color = "black"),
  #         axis.text = element_text(size = tsize),
  #         legend.text = element_text(size = tsize),
  #         axis.title = element_text(size = tsize),
  #         legend.title = element_text(size = tsize)) +
  #   ggtitle("Bray-Curtis Similarity\n(2014-2016) vs (2017-2018)")
  # 
  # plots <- wint + spr + summ + fall + plot_layout(guides = "collect")
  
  pdf(file = community_diff_fig, width = 8, height = 8)
  print(example_plot)
  dev.off()
  
  
}

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
