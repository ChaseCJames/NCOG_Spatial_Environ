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
library(gganimate)
library(gifski)
library(chron)
library(purrr)
library(magick)
library(reshape2)
library(metR)
library(ncdf4)
library(patchwork)

#### SOM Figure ####

som_figure <- function(map_file = "output/bacteria_m_euks_16s_map.Rdata",
                       figure_name = paste0("figures/som_maps/bacteria_16s_som_",Sys.Date(),".pdf"),
                       main = "16s Bacteria", cluster1 = "Nearshore", cluster2 = "Offshore",
                       tsize = 12, psize = 6){
  

  
  map <- map_data("world")   

  load(map_file)
  
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
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",clust1)), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "darkred", limits = c(0,1)) +
    ggtitle(paste0(cluster1)) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank(),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize))

  if(clust1 == 1){p1 <- p1 + geom_point(aes(x = wt_1@coords[1], y = wt_1@coords[2]), color = "blue", size = 5, pch = 10)}
  if(clust1 == 2){p1 <- p1 + geom_point(aes(x = wt_2@coords[1], y = wt_2@coords[2]), color = "blue", size = 5, pch = 10)}
  
  p2 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",clust2)), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "darkblue", limits = c(0,1)) +
    ggtitle(paste0(cluster2)) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank(),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize))
  
  if(clust2 == 1){p2 <- p2 + geom_point(aes(x = wt_1@coords[1], y = wt_1@coords[2]), color = "red", size = 5, pch = 10)}
  if(clust2 == 2){p2 <- p2 + geom_point(aes(x = wt_2@coords[1], y = wt_2@coords[2]), color = "red", size = 5, pch = 10)}
  
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

  if(clust1 == 1){som_plots$cluster[which(som_plots$cluster == "som_1")] <- "Nearshore"}
  if(clust1 == 2){som_plots$cluster[which(som_plots$cluster == "som_2")] <- "Nearshore"}
  
  if(clust2 == 1){som_plots$cluster[which(som_plots$cluster == "som_1")] <- "Offshore"}
  if(clust2 == 2){som_plots$cluster[which(som_plots$cluster == "som_2")] <- "Offshore"}
  
  
  reg_plot <- ggplot(som_plots, aes_string(x = var, y = "freq", color = "cluster")) + geom_point() +
    stat_smooth(method="glm", method.args = list(family = "binomial")) + 
    theme(legend.title = element_blank(), 
          panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    xlab(var_name) + ylab("Frequency") + scale_color_manual(values = c("red", "blue")) + ggtitle(main)
  
  
  
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
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = mean(som_maps$shannon)) +
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
  
  asv_sums <- asv_sums[which(!is.na(asv_sums$station)),]
  
  station_sums <- asv_sums %>% 
    group_by(station) %>% 
    summarise_all(sum, na.rm = TRUE)
  
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
    geom_point(aes(y = shannon, color = "Alpha Diversity")) +
    stat_smooth(aes(y = shannon), method = "loess", level = 0.95, color = "black") +
    geom_point(aes(y = Gamma_Diversity, color = "Gamma Diversity")) +
    stat_smooth(aes(y = Gamma_Diversity), method = "loess", level = 0.95, color = "black") +
    scale_y_continuous(sec.axis = sec_axis(~., name = "Gamma Diversity")) +
    ylab("Alpha Diversity") + xlab("Distance to Coast (km)") +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.position = "bottom", plot.title = element_text(hjust =0.5)) +
    ggtitle(paste0(main,"\nAlpha Diversity vs\nGamma Diversity per Station")) +
    scale_color_manual(values = c("royalblue2", "seagreen3")) 
  
  pdf(file = figure_name, width = 8, height = 6)
  print(div_plot)
  dev.off()
  
  return(div_plot)
  
}


beta_diversity_figure <- function(full_data_file = "output/bacteria_m_euks_16s_full_data.Rdata",
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

### AIC Tables ####

aic_table_func <- function(som_maps = cyano_plots){  
  
  som_glm <- som_maps
  
  model_AIC <- matrix(NA,19,2)
  
  # temperature
  
  glm_mean_temp <- glm(som_1 ~  temp_mean, data = som_glm, family = binomial)
  mt_sum <- summary(glm_mean_temp)
  model_AIC[1,2] <- mt_sum$aic
  model_AIC[1,1] <- "Mean Temp"
  
  glm_coeff_temp <- glm(som_1 ~  temp_coeff, data = som_glm, family = binomial)
  ct_sum <- summary(glm_coeff_temp)
  model_AIC[10,2] <- ct_sum$aic
  model_AIC[10,1] <- "Coeff. Var. Temp"
  
  # sea surface temperature
  
  glm_mean_sst <- glm(som_1 ~  sst_mean, data = som_glm, family = binomial)
  mt_sum <- summary(glm_mean_sst)
  model_AIC[2,2] <- mt_sum$aic
  model_AIC[2,1] <- "Mean SST"
  
  glm_coeff_sst <- glm(som_1 ~  sst_coeff, data = som_glm, family = binomial)
  ct_sum <- summary(glm_coeff_sst)
  model_AIC[11,2] <- ct_sum$aic
  model_AIC[11,1] <- "Coeff. Var. SST"
  
  # salinity
  
  glm_mean_sal <- glm(som_1 ~  sal_mean, data = som_glm, family = binomial)
  mt_sum <- summary(glm_mean_sal)
  model_AIC[3,2] <- mt_sum$aic
  model_AIC[3,1] <- "Mean Salinity"
  
  glm_coeff_sal <- glm(som_1 ~  sal_coeff, data = som_glm, family = binomial)
  ct_sum <- summary(glm_coeff_sal)
  model_AIC[12,2] <- ct_sum$aic
  model_AIC[12,1] <- "Coeff. Var. Salinity"
  
  # NO3
  
  glm_mean_no3 <- glm(som_1 ~  NO3_mean, data = som_glm, family = binomial)
  mn_sum <- summary(glm_mean_no3)
  model_AIC[4,2] <- mn_sum$aic
  model_AIC[4,1] <- "Mean NO3"
  
  glm_coeff_no3 <- glm(som_1 ~  NO3_coeff, data = som_glm, family = binomial)
  cn_sum <- summary(glm_coeff_no3)
  model_AIC[13,2] <- cn_sum$aic
  model_AIC[13,1] <- "Coeff. Var. NO3"
  
  # PO4
  
  glm_mean_po4 <- glm(som_1 ~  PO4_mean, data = som_glm, family = binomial)
  mp_sum <- summary(glm_mean_po4)
  model_AIC[5,2] <- mp_sum$aic
  model_AIC[5,1] <- "Mean PO4"
  
  glm_coeff_po4 <- glm(som_1 ~  PO4_coeff, data = som_glm, family = binomial)
  cp_sum <- summary(glm_coeff_po4)
  model_AIC[14,2] <- cp_sum$aic
  model_AIC[14,1] <- "Coeff. Var. PO4"
  
  # SiO3
  
  glm_mean_sio3 <- glm(som_1 ~  SiO3_mean, data = som_glm, family = binomial)
  ms_sum <- summary(glm_mean_sio3)
  model_AIC[6,2] <- ms_sum$aic
  model_AIC[6,1] <- "Mean SiO3"
  
  glm_coeff_sio3 <- glm(som_1 ~  SiO3_coeff, data = som_glm, family = binomial)
  cs_sum <- summary(glm_coeff_sio3)
  model_AIC[15,2] <- cs_sum$aic
  model_AIC[15,1] <- "Coeff. Var. SiO3"
  
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
  
  glm_coeff_mld <- glm(som_1 ~  MLD_coeff, data = som_glm, family = binomial)
  cs_sum <- summary(glm_coeff_mld)
  model_AIC[16,2] <- cs_sum$aic
  model_AIC[16,1] <- "Coeff. Var. MLD"
  
  # NC Depth
  
  glm_mean_nc <- glm(som_1 ~  NC_mean, data = som_glm, family = binomial)
  ms_sum <- summary(glm_mean_nc)
  model_AIC[8,2] <- ms_sum$aic
  model_AIC[8,1] <- "Mean NCD"
  
  glm_coeff_nc <- glm(som_1 ~  NC_coeff, data = som_glm, family = binomial)
  cs_sum <- summary(glm_coeff_nc)
  model_AIC[17,2] <- cs_sum$aic
  model_AIC[17,1] <- "Coeff. Var. NCD"
  
  # Distance to Coast
  
  glm_mean_dc <- glm(som_1 ~  Dist_mean, data = som_glm, family = binomial)
  ms_sum <- summary(glm_mean_dc)
  model_AIC[18,2] <- ms_sum$aic
  model_AIC[18,1] <- "Distance to Coast"
  
  # Chlorophyll
  
  glm_mean_chl <- glm(som_1 ~  Chl_mean, data = som_glm, family = binomial)
  ms_sum <- summary(glm_mean_chl)
  model_AIC[9,2] <- ms_sum$aic
  model_AIC[9,1] <- "Mean Chl-a"
  
  glm_coeff_chl <- glm(som_1 ~  Chl_coeff, data = som_glm, family = binomial)
  cs_sum <- summary(glm_coeff_chl)
  model_AIC[19,2] <- cs_sum$aic
  model_AIC[19,1] <- "Coeff. Var. Chl-a"
  
  
  return(model_AIC)
  
}

aic_table_func_diveristy <- function(som_maps = cyano_plots, i = 2){  
  
  som_glm <- som_maps
  
  colnames(som_glm)[i] <- "response"
  
  model_AIC <- matrix(NA,19,2)
  model_p_val <- matrix(NA,19,2)
  
  
  # temperature
  
  glm_mean_temp <- glm(response ~  temp_mean, data = som_glm)
  mt_sum <- summary(glm_mean_temp)
  model_AIC[1,2] <- mt_sum$aic
  model_AIC[1,1] <- "Mean Temp"
  
  glm_coeff_temp <- glm(response ~  temp_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_temp)
  model_AIC[10,2] <- ct_sum$aic
  model_AIC[10,1] <- "Coeff. Var. Temp"
  
  # sea surface temperature
  
  glm_mean_sst <- glm(response ~  sst_mean, data = som_glm)
  mt_sum <- summary(glm_mean_sst)
  model_AIC[2,2] <- mt_sum$aic
  model_AIC[2,1] <- "Mean SST"
  
  glm_coeff_sst <- glm(response ~  sst_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_sst)
  model_AIC[11,2] <- ct_sum$aic
  model_AIC[11,1] <- "Coeff. Var. SST"
  
  # salinity
  
  glm_mean_sal <- glm(response ~  sal_mean, data = som_glm)
  mt_sum <- summary(glm_mean_sal)
  model_AIC[3,2] <- mt_sum$aic
  model_AIC[3,1] <- "Mean Salinity"
  
  glm_coeff_sal <- glm(response ~  sal_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_sal)
  model_AIC[12,2] <- ct_sum$aic
  model_AIC[12,1] <- "Coeff. Var. Salinity"
  
  # NO3
  
  glm_mean_no3 <- glm(response ~  NO3_mean, data = som_glm)
  mn_sum <- summary(glm_mean_no3)
  model_AIC[4,2] <- mn_sum$aic
  model_AIC[4,1] <- "Mean NO3"
  
  glm_coeff_no3 <- glm(response ~  NO3_coeff, data = som_glm)
  cn_sum <- summary(glm_coeff_no3)
  model_AIC[13,2] <- cn_sum$aic
  model_AIC[13,1] <- "Coeff. Var. NO3"
  
  # PO4
  
  glm_mean_po4 <- glm(response ~  PO4_mean, data = som_glm)
  mp_sum <- summary(glm_mean_po4)
  model_AIC[5,2] <- mp_sum$aic
  model_AIC[5,1] <- "Mean PO4"
  
  glm_coeff_po4 <- glm(response ~  PO4_coeff, data = som_glm)
  cp_sum <- summary(glm_coeff_po4)
  model_AIC[14,2] <- cp_sum$aic
  model_AIC[14,1] <- "Coeff. Var. PO4"
  
  # SiO3
  
  glm_mean_sio3 <- glm(response ~  SiO3_mean, data = som_glm)
  ms_sum <- summary(glm_mean_sio3)
  model_AIC[6,2] <- ms_sum$aic
  model_AIC[6,1] <- "Mean SiO3"
  
  glm_coeff_sio3 <- glm(response ~  SiO3_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_sio3)
  model_AIC[15,2] <- cs_sum$aic
  model_AIC[15,1] <- "Coeff. Var. SiO3"
  
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
  
  glm_coeff_mld <- glm(response ~  MLD_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_mld)
  model_AIC[16,2] <- cs_sum$aic
  model_AIC[16,1] <- "Coeff. Var. MLD"
  
  # NC Depth
  
  glm_mean_nc <- glm(response ~  NC_mean, data = som_glm)
  ms_sum <- summary(glm_mean_nc)
  model_AIC[8,2] <- ms_sum$aic
  model_AIC[8,1] <- "Mean NCD"
  
  glm_coeff_nc <- glm(response ~  NC_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_nc)
  model_AIC[17,2] <- cs_sum$aic
  model_AIC[17,1] <- "Coeff. Var. NCD"
  
  # Distance to Coast
  
  glm_mean_dc <- glm(response ~  Dist_mean, data = som_glm)
  ms_sum <- summary(glm_mean_dc)
  model_AIC[18,2] <- ms_sum$aic
  model_AIC[18,1] <- "Distance to Coast"
  
  # Chlorophyll
  
  glm_mean_chl <- glm(response ~  Chl_mean, data = som_glm)
  ms_sum <- summary(glm_mean_chl)
  model_AIC[9,2] <- ms_sum$aic
  model_AIC[9,1] <- "Mean Chl-a"
  
  glm_coeff_chl <- glm(response ~  Chl_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_chl)
  model_AIC[19,2] <- cs_sum$aic
  model_AIC[19,1] <- "Coeff. Var. Chl-a"
  
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
  
  for (i in 1:length(in_group_list)) {
    
    load(paste0("output/",in_group_list[i], "_map.Rdata"))
    som_maps2 <- som_maps[which(som_maps$n_samps > (minimum_tp-1)),]
    AIC_table <- aic_table_func(som_maps = som_maps2)
    AIC_table <- as.data.frame(AIC_table)
    colnames(AIC_table) <- c("Variables","AIC")
    AIC_table <- AIC_table[c(1,3,4,5,6,9,8,10,12,13,14,15,19,17,18),]
    AIC_table[,2] <- as.numeric(as.character(AIC_table[,2]))
    AIC_table[,2] <- round(AIC_table[,2], digits = 2)
    
    map_list[[i]] <- AIC_table

    }
  
  AIC_full <- full_join(map_list[[1]], map_list[[2]], by = "Variables")
  
  for (i in 3:length(map_list)) {
    
    AIC_full <- full_join(AIC_full, map_list[[i]], by = "Variables")
    
  }
  
  colnames(AIC_full) <- c("Variables", in_group_names)
  
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
  
  plot_df <- melt(AIC_scaled)
  
  
  
  colnames(plot_df) <- c("Variables", "Group", "AIC")
  
  # # removing plastids from plot
  # plot_df <- plot_df[-which(plot_df$Group == "Eukaryotic\nPlastid\n AIC"),]
  
  plot_df$Variables <- as.factor(plot_df$Variables)
  plot_df$Variables <- factor(plot_df$Variables, levels = c("Distance to Coast",
                                                            "Coeff. Var. NCD","Coeff. Var. Chl-a",
                                                            "Coeff. Var. SiO3","Coeff. Var. PO4",
                                                            "Coeff. Var. NO3", "Coeff. Var. Salinity",
                                                            "Coeff. Var. Temp",
                                                            "Mean NCD", "Mean Chl-a",
                                                            "Mean SiO3", "Mean PO4",
                                                            "Mean NO3","Mean Salinity",
                                                            "Mean Temp" ))
  
  pdf(figure_name_2, width = width_plot, height = 8)
  
  print(ggplot(data = plot_df, aes(x = Group, y = Variables, size = AIC)) + 
          geom_point(fill = "red", color = "black", alpha = 0.6, shape = 21) +
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
  
  a <- ggplot(data = plot_df, aes(x = Group, y = Variables, size = AIC)) + 
    geom_point(fill = "red", color = "black", alpha = 0.6, shape = 21) +
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
  
  for (i in 1:length(in_group_list)) {
    
    load(paste0("output/",in_group_list[i], "_map.Rdata"))
    AIC_table <- aic_table_func_diveristy(som_maps = som_maps, i = col)
    AIC_table <- as.data.frame(AIC_table, stringsAsFactors = FALSE)
    colnames(AIC_table) <- c("Variables","AIC")
    AIC_table <- AIC_table[c(1,3,4,5,6,9,8,10,12,13,14,15,19,17,18),]
    AIC_table[,2] <- as.numeric(AIC_table[,2])
    AIC_table[,2] <- round(AIC_table[,2], digits = 2)
    
    map_list[[i]] <- AIC_table
    
  }
  
  AIC_full <- full_join(map_list[[1]], map_list[[2]], by = "Variables")
  
  for (i in 3:length(map_list)) {
    
    AIC_full <- full_join(AIC_full, map_list[[i]], by = "Variables")
    
  }
  
  colnames(AIC_full) <- c("Variables", in_group_names)
  
  
  
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
  
  colnames(plot_df) <- c("Variables", "Group", "AIC")
  
  # # removing plastids from plot
  # plot_df <- plot_df[-which(plot_df$Group == "Eukaryotic\nPlastid\n AIC"),]
  
  plot_df$Variables <- as.factor(plot_df$Variables)
  plot_df$Variables <- factor(plot_df$Variables, levels = c("Distance to Coast",
                                                            "Coeff. Var. NCD","Coeff. Var. Chl-a",
                                                            "Coeff. Var. SiO3","Coeff. Var. PO4",
                                                            "Coeff. Var. NO3", "Coeff. Var. Salinity",
                                                            "Coeff. Var. Temp",
                                                            "Mean NCD", "Mean Chl-a",
                                                            "Mean SiO3", "Mean PO4",
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

aic_table_func_diveristy_sign <- function(som_maps = cyano_plots, i = 2){  
  
  som_glm <- som_maps
  
  colnames(som_glm)[i] <- "response"
  
  model_AIC <- matrix(NA,19,3)
  # temperature
  
  glm_mean_temp <- glm(response ~  temp_mean, data = som_glm)
  mt_sum <- summary(glm_mean_temp)
  model_AIC[1,2] <- mt_sum$aic
  model_AIC[1,1] <- "Mean Temp"
  model_AIC[1,3] <- cor(som_glm$temp_mean, som_glm$response, use = "pairwise.complete.obs")
  
  glm_coeff_temp <- glm(response ~  temp_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_temp)
  model_AIC[10,2] <- ct_sum$aic
  model_AIC[10,1] <- "Coeff. Var. Temp"
  model_AIC[10,3] <- cor(som_glm$temp_coeff, som_glm$response, use = "pairwise.complete.obs")
  
  # sea surface temperature
  
  glm_mean_sst <- glm(response ~  sst_mean, data = som_glm)
  mt_sum <- summary(glm_mean_sst)
  model_AIC[2,2] <- mt_sum$aic
  model_AIC[2,1] <- "Mean SST"
  model_AIC[2,3] <- cor(som_glm$sst_mean, som_glm$response, use = "pairwise.complete.obs")
  
  glm_coeff_sst <- glm(response ~  sst_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_sst)
  model_AIC[11,2] <- ct_sum$aic
  model_AIC[11,1] <- "Coeff. Var. SST"
  model_AIC[11,3] <- cor(som_glm$sst_coeff, som_glm$response, use = "pairwise.complete.obs")
  
  # salinity
  
  glm_mean_sal <- glm(response ~  sal_mean, data = som_glm)
  mt_sum <- summary(glm_mean_sal)
  model_AIC[3,2] <- mt_sum$aic
  model_AIC[3,1] <- "Mean Salinity"
  model_AIC[3,3] <- cor(som_glm$sal_mean, som_glm$response, use = "pairwise.complete.obs")
  
  glm_coeff_sal <- glm(response ~  sal_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_sal)
  model_AIC[12,2] <- ct_sum$aic
  model_AIC[12,1] <- "Coeff. Var. Salinity"
  model_AIC[12,3] <- cor(som_glm$sal_coeff, som_glm$response, use = "pairwise.complete.obs")
  
  # NO3
  
  glm_mean_no3 <- glm(response ~  NO3_mean, data = som_glm)
  mn_sum <- summary(glm_mean_no3)
  model_AIC[4,2] <- mn_sum$aic
  model_AIC[4,1] <- "Mean NO3"
  model_AIC[4,3] <- cor(som_glm$NO3_mean, som_glm$response, use = "pairwise.complete.obs")
  
  glm_coeff_no3 <- glm(response ~  NO3_coeff, data = som_glm)
  cn_sum <- summary(glm_coeff_no3)
  model_AIC[13,2] <- cn_sum$aic
  model_AIC[13,1] <- "Coeff. Var. NO3"
  model_AIC[13,3] <- cor(som_glm$NO3_coeff, som_glm$response, use = "pairwise.complete.obs")
  
  # PO4
  
  glm_mean_po4 <- glm(response ~  PO4_mean, data = som_glm)
  mp_sum <- summary(glm_mean_po4)
  model_AIC[5,2] <- mp_sum$aic
  model_AIC[5,1] <- "Mean PO4"
  model_AIC[5,3] <- cor(som_glm$PO4_mean, som_glm$response, use = "pairwise.complete.obs")
  
  glm_coeff_po4 <- glm(response ~  PO4_coeff, data = som_glm)
  cp_sum <- summary(glm_coeff_po4)
  model_AIC[14,2] <- cp_sum$aic
  model_AIC[14,1] <- "Coeff. Var. PO4"
  model_AIC[14,3] <- cor(som_glm$PO4_coeff, som_glm$response, use = "pairwise.complete.obs")
  
  # SiO3
  
  glm_mean_sio3 <- glm(response ~  SiO3_mean, data = som_glm)
  ms_sum <- summary(glm_mean_sio3)
  model_AIC[6,2] <- ms_sum$aic
  model_AIC[6,1] <- "Mean SiO3"
  model_AIC[6,3] <- cor(som_glm$SiO3_mean, som_glm$response, use = "pairwise.complete.obs")
  
  glm_coeff_sio3 <- glm(response ~  SiO3_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_sio3)
  model_AIC[15,2] <- cs_sum$aic
  model_AIC[15,1] <- "Coeff. Var. SiO3"
  model_AIC[15,3] <- cor(som_glm$SiO3_coeff, som_glm$response, use = "pairwise.complete.obs")
  
  # MLD
  
  glm_mean_mld <- glm(response ~  MLD_mean, data = som_glm)
  ms_sum <- summary(glm_mean_mld)
  model_AIC[7,2] <- ms_sum$aic
  model_AIC[7,1] <- "Mean MLD"
  model_AIC[7,3] <- cor(som_glm$MLD_mean, som_glm$response, use = "pairwise.complete.obs")
  
  glm_coeff_mld <- glm(response ~  MLD_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_mld)
  model_AIC[16,2] <- cs_sum$aic
  model_AIC[16,1] <- "Coeff. Var. MLD"
  model_AIC[16,3] <- cor(som_glm$MLD_coeff, som_glm$response, use = "pairwise.complete.obs")
  
  # NC Depth
  
  glm_mean_nc <- glm(response ~  NC_mean, data = som_glm)
  ms_sum <- summary(glm_mean_nc)
  model_AIC[8,2] <- ms_sum$aic
  model_AIC[8,1] <- "Mean NCD"
  model_AIC[8,3] <- cor(som_glm$NC_mean, som_glm$response, use = "pairwise.complete.obs")
  
  glm_coeff_nc <- glm(response ~  NC_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_nc)
  model_AIC[17,2] <- cs_sum$aic
  model_AIC[17,1] <- "Coeff. Var. NCD"
  model_AIC[17,3] <- cor(som_glm$NC_coeff, som_glm$response, use = "pairwise.complete.obs")
  
  # Distance to Coast
  
  glm_mean_dc <- glm(response ~  Dist_mean, data = som_glm)
  ms_sum <- summary(glm_mean_dc)
  model_AIC[18,2] <- ms_sum$aic
  model_AIC[18,1] <- "Distance to Coast"
  model_AIC[18,3] <- cor(som_glm$Dist_mean, som_glm$response, use = "pairwise.complete.obs")
  
  # Chlorophyll
  
  glm_mean_chl <- glm(response ~  Chl_mean, data = som_glm)
  ms_sum <- summary(glm_mean_chl)
  model_AIC[9,2] <- ms_sum$aic
  model_AIC[9,1] <- "Mean Chl-a"
  model_AIC[9,3] <- cor(som_glm$Chl_mean, som_glm$response, use = "pairwise.complete.obs")
  
  glm_coeff_chl <- glm(response ~  Chl_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_chl)
  model_AIC[19,2] <- cs_sum$aic
  model_AIC[19,1] <- "Coeff. Var. Chl-a"
  model_AIC[19,3] <- cor(som_glm$Chl_coeff, som_glm$response, use = "pairwise.complete.obs")
  
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
  
  for (i in 1:length(in_group_list)) {
    
    load(paste0("output/",in_group_list[i], "_map.Rdata"))
    AIC_table <- aic_table_func_diveristy_sign(som_maps = som_maps, i = col)
    AIC_table <- as.data.frame(AIC_table, stringsAsFactors = FALSE)
    colnames(AIC_table) <- c("Variables","AIC", "Slope")
    AIC_table <- AIC_table[c(1,3,4,5,6,9,8,10,12,13,14,15,19,17,18),]
    AIC_table[,2] <- as.numeric(AIC_table[,2])
    AIC_table[,2] <- round(AIC_table[,2], digits = 2)
    AIC_table[,3] <- as.numeric(AIC_table[,3])
    AIC_table[,3] <- round(AIC_table[,3], digits = 3)
    
    map_list[[i]] <- AIC_table
    
  }
  
  AIC_full <- full_join(map_list[[1]], map_list[[2]], by = "Variables")
  
  for (i in 3:length(map_list)) {
    
    AIC_full <- full_join(AIC_full, map_list[[i]], by = "Variables")
    
  }
  
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
  
  plot_df$slope <- plot_slope$value
  
  colnames(plot_df) <- c("Variables", "Group", "AIC","slope")
  
  # # removing plastids from plot
  # plot_df <- plot_df[-which(plot_df$Group == "Eukaryotic\nPlastid\n AIC"),]
  
  plot_df$Variables <- as.factor(plot_df$Variables)
  plot_df$Variables <- factor(plot_df$Variables, levels = c("Distance to Coast",
                                                            "Coeff. Var. NCD","Coeff. Var. Chl-a",
                                                            "Coeff. Var. SiO3","Coeff. Var. PO4",
                                                            "Coeff. Var. NO3", "Coeff. Var. Salinity",
                                                            "Coeff. Var. Temp",
                                                            "Mean NCD", "Mean Chl-a",
                                                            "Mean SiO3", "Mean PO4",
                                                            "Mean NO3","Mean Salinity",
                                                            "Mean Temp" ))
  
  pdf(figure_name_2, width = width_plot, height = 8)
  print(ggplot(data = plot_df, aes(x = Group, y = Variables, size = AIC, fill = slope)) + 
          geom_point(color = "black", alpha = 0.6, shape = 21) +
    scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "white", midpoint = 0) +
          labs(fill = "Correlation") + ylab("Variable") +
          theme(panel.background = element_blank(),
                panel.border = element_rect(color = "black", fill = NA),
                legend.position = "right",
                panel.grid.major.y = element_line(color = "grey", linetype = 2),
                plot.title = element_text(hjust = 0, size = tsize),
                axis.text.x = element_text(angle = 0)) +
          scale_size_continuous(range = c(1,18), guide = FALSE) + xlab("") + ylab(""))
  dev.off()
  
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
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") + 
    geom_point(size = 3, aes_string(color = "phase", shape = "season"), data = som_cruise) +
    scale_color_manual(values = c("red", "blue","gold3")) +
    scale_shape_manual(values = c(0,1,2,3)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Season", color = "Phase") + xlab("Nearshore-Offshore\nSlope in Nitracline") +
    ylab("Frequency of Nearshore Cluster") + ggtitle(title)
  
  print(gradient_plot)
  
  ts_plot <- ggplot(som_cruise, aes_string(x = "Date", y = paste0("som_",nearshore_som))) +
    geom_line(color = "black", size = 1) + 
    geom_point(size = 3, aes_string(color = "season"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("Date") +
    ylab("Frequency of Nearshore Cluster")  + ggtitle(title)
  
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
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") + 
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17,18)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("Coastal Upwelling Transport Index\n(CUTI)") +
    ylab("Frequency of Nearshore Cluster") 
  
  beuti_plot <- ggplot(som_cruise, aes_string(x = "BEUTI", y = paste0("som_",nearshore_som))) +
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
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") + 
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17,18)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("Regionally Availible Nitrate") +
    ylab("") 
  
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
  
  save(early_dist, late_dist, diversity_diff, even_diff, gradient_plot, ts_plot,
       cuti_plot, beuti_plot, reg_nitrate, dissimilar_plot, file = out_diff_file)
  
  
  
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
  
  asv_cruise <- asv_table %>% group_by(cruise) %>% summarise_all(list(~ sum(.,na.rm = TRUE)))
  
  asv_cruise$total_div <- diversity(asv_cruise[,-1], MARGIN = 1, index = "shannon")
  
  som_cruise$total_shannon <- asv_cruise$total_div[match(som_cruise$Cruise, asv_cruise$cruise)]
  
  
  gradient_plot <- ggplot(som_cruise, aes_string(x = "NC_slope", y = "shannon")) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") + 
    geom_point(size = 3, aes_string(color = "phase", shape = "season"), data = som_cruise) +
    scale_color_manual(values = c("red", "blue", "gold3")) +
    scale_shape_manual(values = c(0,1,2,3)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Season", color = "Phase") + xlab("Nearshore-Offshore\nSlope in Nitracline") +
    ylab("Mean Shannon Diversity\nPer Cruise") + ggtitle(title)
  
  print(gradient_plot)
  
  gradient_plot2 <- ggplot(som_cruise, aes_string(x = "NC_slope", y = "total_shannon")) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") + 
    geom_point(size = 3, aes_string(color = "phase", shape = "season"), data = som_cruise) +
    scale_color_manual(values = c("red", "blue", "gold3")) +
    scale_shape_manual(values = c(0,1,2,3)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Season", color = "Phase") + xlab("Nearshore-Offshore\nSlope in Nitracline") +
    ylab("Total Shannon Diversity\nPer Cruise") + ggtitle(title)
  
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
       gradient_plot, ts_plot, gradient_plot2, ts_plot2,
       cuti_plot, beuti_plot, reg_nitrate, file = out_diff_file)
  
}

###### Outputs #####

#### SOM Maps ####

# All

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9", "bacteria_m_euks_16s",
                  "plastid_16s", "cyano_16s", "euks_hetero_18sv9")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                   "Sar Clade", "Archaea","Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans", "Bacteria",
                   "Eukaryotic Phytoplankton (Plastids)", "Cyanobacteria",
                   "Eukaryotic Protists")


plot_list <- list()

for (i in 1:length(in_group_list)) {
  
  plot_list[[i]] <- som_figure(map_file = paste0("output/",in_group_list[i],"_map.Rdata"),
                               figure_name = paste0("figures/som_maps/", in_group_list[i],"_map_plot.pdf"),
                               main = in_group_names[i], cluster1 = "Nearshore", cluster2 = "Offshore")
  
  
}

#### Regression Plots #####

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9", "bacteria_m_euks_16s",
                  "plastid_16s", "cyano_16s", "euks_hetero_18sv9")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                   "Sar Clade", "Archaea","Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans", "Bacteria",
                   "Eukaryotic Phytoplankton (Plastids)", "Cyanobacteria",
                   "Eukaryotic Protists")

var_list = c("temp_mean", "sal_mean", "PO4_mean", "NO3_mean", "SiO3_mean", "NC_mean",
             "temp_coeff", "sal_coeff", "PO4_coeff", "NO3_coeff", "SiO3_coeff", "NC_coeff", "Dist_mean")

var_name_list = c("Mean Temperature (°C)", "Mean Salinity", "Mean PO4ug", "Mean NO3ug",
                  "Mean SiO3ug", "Mean Nitracline Depth (m)",
             "Coeff. Var. Temperature", "Coeff. Var. Salinity",
             "Coeff. Var. PO4", "Coeff. Var. NO3", "Coeff. Var. SiO3",
             "Coeff. Var. Nitracline Depth (m)", "Distance to Coast (km)")

for (i in 1:length(in_group_list)) {
  for (j in 1:length(var_list)) {
    regression_figure(glm_file = paste0("output/",in_group_list[i], "_glm.Rdata"),
                      map_file = paste0("output/",in_group_list[i],"_map.Rdata"),   
                      figure_name = paste0("figures/glm_plots/", in_group_list[i],"_"),
                      main = in_group_names[i], cluster1 = "Nearshore", cluster2 = "Offshore",
                      var = var_list[j], var_name = var_name_list[j])
  }
}


#### Diversity Figures ####

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9", "bacteria_m_euks_16s",
                  "plastid_16s", "cyano_16s", "euks_hetero_18sv9")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                   "Sar Clade", "Archaea","Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans", "Bacteria",
                   "Eukaryotic Phytoplankton (Plastids)", "Cyanobacteria",
                   "Eukaryotic Protists")

in_group_list_basic = c("16s_pro", "16s_syne","16s_flavo", "16s_rhodo", "16s_sar", "16s_archaea",
                  "18s_diatom","18s_dino", "18s_syndin",
                  "18s_hapto", "18s_chloro", "18s_metazoa", "16s_bacteria_m_euks",
                  "16s_plastids", "16s_cyanos", "18s_heterotrophic_euks")


for (i in 1:length(in_group_list)) {
  
  diveristy_figure(map_file = paste0("output/", in_group_list[i], "_map.Rdata"),
                   full_dat = paste0("output/", in_group_list[i], "_full_data.Rdata"),
                   figure_start = paste0("figures/diversity/", in_group_list[i], "_"),
                   main = in_group_names[i])
  
  alpha_versus_gamma_figure(full_data_file = paste0("output/", in_group_list[i], "_full_data.Rdata"),
                            raw_data_file = paste0("data/", in_group_list_basic[i], ".Rdata"),
                            map_file = paste0("output/", in_group_list[i], "_map.Rdata"), minimum_tp = 8,
                            figure_name = paste0("figures/diversity/", in_group_list[i], "_alpha_gamma.pdf"),
                            main = in_group_names[i])
  
  beta_diversity_figure(full_data_file = paste0("output/", in_group_list[i], "_full_data.Rdata"),
                        bc_data_file = paste0("output/", in_group_list[i], "_dissimilar.Rdata"),
                        map_file = paste0("output/", in_group_list[i], "_map.Rdata"), minimum_tp = 8,
                        figure_name = paste0("figures/diversity/", in_group_list[i],"_beta.pdf"),
                        main = in_group_names[i])
  
}

#### AIC Figures ####

full_aic_table_figure(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                        "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                        "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                      in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar Clade", 
                                         "Diatoms",
                                         "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"),
                      minimum_tp = 4, width_plot = 16,
                      figure_name_2 = paste0("figures/aic_figures/small_group_aic_plot_logit",".pdf"),
                      title_name = "Variable Importance")

full_aic_table_figure(in_group_list = c("archaea_16s","bacteria_m_euks_16s", "cyano_16s","plastid_16s",
                                        "euks_hetero_18sv9"),
                      in_group_names = c("Archaea","Bacteria", "Cyanobacteria", "Eukaryotic Phytoplankton\n(Plastids)",
                                         "Eukaryotic\n Protists"),
                      minimum_tp = 4, width_plot = 10,
                      figure_name_2 = paste0("figures/aic_figures/big_group_aic_plot_logit",".pdf"),
                      title_name = "Variable Importance")


# Diversity

full_aic_table_figure_diversity(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                                in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar Clade", 
                                                   "Diatoms",
                                                   "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"),
                                figure_name_2 = paste0("figures/aic_figures/small_group_aic_plot_logit_even",".pdf"),
                                title_name = "Variable Importance Evenness", # col 2 = even, col 3 = shan col 4 = rich
                                col = 2, color_fill = "purple")

full_aic_table_figure_diversity(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                                in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar Clade", 
                                                   "Diatoms",
                                                   "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"),
                                figure_name_2 = paste0("figures/aic_figures/small_group_aic_plot_logit_shannon",".pdf"),
                                title_name = "Variable Importance Shannon Diversity", # col 2 = even, col 3 = shan col 4 = rich
                                col = 3, color_fill = "red")

full_aic_table_figure_diversity(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                                in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar Clade", 
                                                   "Diatoms",
                                                   "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"),
                                figure_name_2 = paste0("figures/aic_figures/small_group_aic_plot_logit_rich",".pdf"),
                                title_name = "Variable Importance Richness", # col 2 = even, col 3 = shan col 4 = rich
                                col = 4, color_fill = "blue")

# Basic Groups

full_aic_table_figure_diversity(in_group_list = c("archaea_16s","bacteria_m_euks_16s", "cyano_16s","plastid_16s",
                                                   "euks_hetero_18sv9"),
                                in_group_names = c("Archaea","Bacteria", "Cyanobacteria", "Eukaryotic Phytoplankton\n(Plastids)",
                                                   "Eukaryotic\n Protists"),
                                figure_name_2 = paste0("figures/aic_figures/big_group_aic_plot_logit_even",".pdf"),
                                title_name = "Variable Importance Evenness", # col 2 = even, col 3 = shan col 4 = rich
                                col = 2, color_fill = "purple", width_plot = 8)

full_aic_table_figure_diversity(in_group_list = c("archaea_16s","bacteria_m_euks_16s", "cyano_16s","plastid_16s",
                                                   "euks_hetero_18sv9"),
                                in_group_names = c("Archaea","Bacteria", "Cyanobacteria", "Eukaryotic Phytoplankton\n(Plastids)",
                                                   "Eukaryotic\n Protists"),
                                figure_name_2 = paste0("figures/aic_figures/big_group_aic_plot_logit_shannon",".pdf"),
                                title_name = "Variable Importance\nShannon Diversity", # col 2 = even, col 3 = shan col 4 = rich
                                col = 3, color_fill = "red", width_plot = 8)

full_aic_table_figure_diversity(in_group_list = c("archaea_16s","bacteria_m_euks_16s", "cyano_16s","plastid_16s",
                                                  "euks_hetero_18sv9"),
                                in_group_names = c("Archaea","Bacteria", "Cyanobacteria", "Eukaryotic Phytoplankton\n(Plastids)",
                                                   "Eukaryotic\n Protists"),
                                figure_name_2 = paste0("figures/aic_figures/big_group_aic_plot_logit_rich",".pdf"),
                                title_name = "Variable Importance Richness", # col 2 = even, col 3 = shan col 4 = rich
                                col = 4, color_fill = "blue", width_plot = 8)

###### Community Time Plots ########

# All


in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9", "bacteria_m_euks_16s",
                  "plastid_16s", "cyano_16s", "euks_hetero_18sv9")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                   "Sar Clade", "Archaea","Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans", "Bacteria",
                   "Eukaryotic Phytoplankton (Plastids)", "Cyanobacteria",
                   "Eukaryotic Protists")

for (i in 1:length(in_group_list)) {
  
  community_comparison(in_file = paste0("output/",in_group_list[i],"_full_data.Rdata"),
                       in_map = paste0("output/",in_group_list[i],"_map.Rdata"),
                       similar_mat = paste0("output/",in_group_list[i],"_dissimilar.Rdata"),
                       out_diff_file = paste0("output/",in_group_list[i],"_diffs.Rdata"),
                       title = in_group_names[i],
                       upwelling_index = "output/upwelling_indicies.Rdata",
                       index_plot = paste0("figures/",in_group_list[i],"_index_plot.pdf"))
  
}

###### Diversity Time Plots #####

# All

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9", "bacteria_m_euks_16s",
                  "plastid_16s", "cyano_16s", "euks_hetero_18sv9")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                   "Sar Clade", "Archaea","Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans", "Bacteria",
                   "Eukaryotic Phytoplankton (Plastids)", "Cyanobacteria",
                   "Eukaryotic Protists")

in_group_list_basic = c("16s_pro", "16s_syne","16s_flavo", "16s_rhodo", "16s_sar", "16s_archaea",
                        "18s_diatom","18s_dino", "18s_syndin",
                        "18s_hapto", "18s_chloro", "18s_metazoa", "16s_bacteria_m_euks",
                        "16s_plastids", "16s_cyanos", "18s_heterotrophic_euks")

for (i in 1:length(in_group_list)) {
  
  diversity_comparison(in_file = paste0("output/",in_group_list[i],"_full_data.Rdata"),
                       in_map = paste0("output/",in_group_list[i],"_map.Rdata"),
                       in_raw = paste0("data/", in_group_list_basic[i],".Rdata"),
                       out_diff_file = paste0("output/",in_group_list[i],"_diffs_div.Rdata"),
                       title = in_group_names[i],
                       upwelling_index = "output/upwelling_indicies.Rdata",
                       index_plot = paste0("figures/",in_group_list[i],"_index_plot_div.pdf"))
  
}

###### Main Text Figures #####

#### Figure 1: Physical Conditions ####

fig_1_func <- function(in_vel = "output/uv_velocity_table.Rdata",
                         in_bath = "output/CALCOFI_bathymetry_table.Rdata",
                         in_temp = "output/CALCOFI_temp_tables.Rdata",
                         in_cyano = "output/cyano_16s_map.Rdata",
                         fig_name = "figures/figure_outline/fig_1.pdf",
                         tsize = 12, psize = 6){
  
  map <- map_data("world")    
  
  load(in_vel)
  load(in_bath)
  load(in_temp)
  load(in_cyano)
  
  
  stations <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "n_samps"),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "# of Samples", low = "white", high = "red") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("A.")
  
  
  reduced_uv_table <- uv_table[seq(1,1633,6),]
  elevation_table$value <- abs(elevation_table$value)
  
  vel_bath <- ggplot() + 
    geom_raster(data = elevation_table, aes(x = lon, y = lat, fill = value), interpolate = FALSE) +
    scale_fill_gradient(low = "darkblue", high = "cyan", name = "Depth (m)", trans = 'reverse') +
    metR::geom_vector(data = reduced_uv_table, aes(x = lon, y = lat, dx = Mean_U, dy = Mean_V), 
                      arrow.angle = 15, arrow.type = "open", arrow.length = unit(0.5, "inches"), 
                      pivot = 0,preserve.dir = TRUE, direction = "ccw",
                      min.mag = 0, show.legend = NA, color = "white") +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    labs(x = "Longitude", y = "Latitude", mag = "Speed (m/s)", color = "Depth (m)")  +
    theme(panel.background = element_blank(),
          legend.key = element_rect(fill = "grey"),
          panel.border = element_rect(fill = NA, colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.box = "horizontal",
          legend.title = element_text(size = tsize)) + 
    guides(fill = guide_legend(order = 2),mag = guide_legend(order = 1)) +
    ggtitle("B.")
  
  
  sst <- ggplot() + 
    geom_tile(data = coeff_table, aes(x = lon, y = lat, fill = coeff_var), width =0.26, height = 0.26) +
    scale_fill_gradient2(name = "Coeff. Var SST", low = "darkblue", mid = "white", high = "darkred", limits = c(0.09,0.12), oob = squish, midpoint = 0.1066851) +geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    ggtitle("Coeff. Var. SST") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  sst_mean <- ggplot() + 
    geom_tile(data = mean_table, aes(x = lon, y = lat, fill = Mean), width =0.26, height = 0.26) +
    scale_fill_gradient2(name = "SST Mean (°C)", low = "darkblue", mid = "white", high = "darkred", limits = c(15,18), oob = squish, midpoint = 16.5) +geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) + 
    ggtitle("D.")
  
  nc_depth <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "NC_mean"),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean Nitracline\nDepth (m)", low = "darkblue", high = "cyan", trans = 'reverse') +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("C.")
  
  
  stations <- stations + theme(axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank())
  
  vel_bath <- vel_bath + theme(axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),
                               axis.title.y=element_blank(),
                               axis.text.y=element_blank(),
                               axis.ticks.y=element_blank())
  
  sst_mean <- sst_mean + theme(axis.title.y=element_blank(),
                               axis.text.y=element_blank(),
                               axis.ticks.y=element_blank())
  
  patch <- stations + vel_bath + nc_depth + sst_mean
  
  pdf(fig_name, width = 12, height = 12)
  print(patch)
  dev.off()
  
}

#### Figure 2: SOMs ####

fig_2_func <- function(in_list = plot_list, file_name = "figures/figure_outline/fig_2.pdf", tsize  = 12){
  
  in_list[[6]]$bp <- in_list[[6]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("A. Archaea")
  in_list[[13]]$bp <- in_list[[13]]$bp + theme(axis.title.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize),
                                               legend.position = "none") + ggtitle("B. Bacteria")
  in_list[[15]]$bp <- in_list[[15]]$bp + theme(axis.title.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize),
                                               legend.position = "none") + ggtitle("C. Cyanobacteria")
  in_list[[14]]$bp <- in_list[[14]]$bp + theme(axis.title.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize),
                                               legend.position = "none") + ggtitle("D. Eukaryotic Phytoplankton")
  in_list[[16]]$bp <- in_list[[16]]$bp + theme(axis.title.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize)) + ggtitle("E. Eukaryotic Protists")
  
  in_list[[6]]$rp <- in_list[[6]]$rp + theme(legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[13]]$rp <- in_list[[13]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[15]]$rp <- in_list[[15]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[14]]$rp <- in_list[[14]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[16]]$rp <- in_list[[16]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0)) + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  
  
  
  pdf(file = file_name, width = 15, height = 6)
  print(plot_list[[6]]$bp + plot_list[[13]]$bp + plot_list[[15]]$bp +
          plot_list[[14]]$bp + plot_list[[16]]$bp +
          plot_list[[6]]$rp + plot_list[[13]]$rp + plot_list[[15]]$rp +
          plot_list[[14]]$rp + plot_list[[16]]$rp +
          plot_layout(ncol = 5))
  dev.off()
}

#### Figure 3: Reg and Var Import ####

fig_3_func <- function(file_name = "figures/figure_outline/fig_3.pdf", tsize = 12){

cyano <- regression_figure(glm_file = "output/cyano_16s_glm.Rdata",
                  map_file = "output/cyano_16s_map.Rdata",   
                  figure_name = "figures/glm_plots/cyano_16s_som_",
                  main = "16s Cyanobacteria", cluster1 = "Nearshore", cluster2 = "Offshore",
                  var = "NC_mean", var_name = "Nitracline Depth (m)")

cyano <- cyano + ggtitle("A.") + 
  theme(legend.justification=c(1,0.5), 
        legend.position=c(0.95, 0.5),
        legend.background = element_rect(fill = "white", color = "black"),
        axis.text = element_text(size = tsize),
        legend.text = element_text(size = tsize),
        axis.title = element_text(size = tsize),
        plot.title = element_text(hjust = 0, size = tsize))

aic_plot <- full_aic_table_figure(in_group_list = c("archaea_16s","bacteria_m_euks_16s", "cyano_16s","plastid_16s",
                                        "euks_hetero_18sv9"),
                      in_group_names = c("Archaea","Bacteria", "Cyanobacteria", "Eukaryotic Phytoplankton\n(Plastids)",
                                         "Eukaryotic\n Protists"),
                      minimum_tp = 4, width_plot = 10,
                      figure_name_2 = paste0("figures/aic_figures/big_group_aic_plot_logit",".pdf"),
                      title_name = "Variable Importance")

aic_plot <- aic_plot + ggtitle("B.") + 
  theme(axis.text = element_text(size = tsize),
        axis.title = element_text(size = tsize),
        plot.title = element_text(hjust = 0, size = tsize))

a <- cyano + aic_plot + plot_layout(widths = c(1,2))

pdf(file = file_name, width = 15, height = 8)
print(a)
dev.off()

}

#### Figure 4: Diversity Example ####

fig_4_func <- function(file_name = "figures/figure_outline/fig_4.pdf",
                       in_group = "diatom_18sv9", basic = "18s_diatom", name = "Diatoms", tsize = 12){
  
  
  map <- diveristy_figure(map_file = paste0("output/", in_group, "_map.Rdata"),
                   full_dat = paste0("output/", in_group, "_full_data.Rdata"),
                   figure_start = paste0("figures/diversity/", in_group, "_"),
                   main = in_group_names[i])
  
  alpha_gamma <- alpha_versus_gamma_figure(full_data_file = paste0("output/", in_group, "_full_data.Rdata"),
                            raw_data_file = paste0("data/", basic, ".Rdata"),
                            map_file = paste0("output/", in_group, "_map.Rdata"), minimum_tp = 8,
                            figure_name = paste0("figures/diversity/", in_group, "_alpha_gamma.pdf"),
                            main = in_group_names[i])
  
  beta <- beta_diversity_figure(full_data_file = paste0("output/", in_group, "_full_data.Rdata"),
                        bc_data_file = paste0("output/", in_group, "_dissimilar.Rdata"),
                        map_file = paste0("output/", in_group, "_map.Rdata"), minimum_tp = 8,
                        figure_name = paste0("figures/diversity/", in_group,"_beta.pdf"),
                        main = in_group_names[i])
  
  
  map <- map + ggtitle("A.") +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(0,1), 
          legend.position=c(0.03, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) + 
    labs(fill = "Mean Alpha\nShannon Diversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  
  alpha_gamma <- alpha_gamma + ggtitle("B.") +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(1,0.5), 
          legend.position=c(0.97, 0.5),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize)) + 
    labs(fill = "Mean Alpha\nShannon Diversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=5))
  
  beta <- beta + ggtitle("C.") +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position= "none",
          axis.text = element_text(size = tsize),
          axis.title = element_text(size = tsize)) + 
    geom_violin(fill = "grey90") + geom_boxplot(width = 0.05, fill = "white")
  
  
  plots <- map + alpha_gamma + beta
  
  pdf(file = file_name, width = 18, height = 7)
  print(plots)
  dev.off()
  
}

##### Figure 5: Diversity Importance #####

full_aic_table_figure_diversity_sign(in_group_list = c("archaea_16s","bacteria_m_euks_16s", "cyano_16s","plastid_16s",
                                                       "euks_hetero_18sv9"),
                                     in_group_names = c("Archaea","Bacteria", "Cyanobacteria", "Eukaryotic Phytoplankton\n(Plastids)",
                                                        "Eukaryotic\n Protists"),
                                     figure_name_2 = "figures/figure_outline/fig_5.pdf",
                                     col = 27, width_plot = 10, tsize = 12)

###### Figure 6: #####

##### Figure 7: Community vs Time #####

fig_7_func <- function(in_phyto = "output/plastid_16s_diffs.Rdata",
                            in_euks = "output/euks_hetero_18sv9_diffs.Rdata",
                            in_cyano = "output/cyano_16s_diffs.Rdata",
                            in_bact = "output/bacteria_m_euks_16s_diffs.Rdata",
                            gradient_plot_file = "figures/figure_outline/fig_7.pdf",
                       tsize = 12){
  
  
  load(in_phyto)
  phyto_gradient <- gradient_plot
  
  load(in_euks)
  euk_gradient <- gradient_plot
  
  load(in_cyano)
  cyano_gradient <- gradient_plot
 
  load(in_bact)
  bact_gradient <- gradient_plot
  
  phyto_gradient <- phyto_gradient + theme(legend.position = "none") + xlab("")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.title.y = element_text(size = tsize),
          axis.text.y = element_text(size = tsize)) + ggtitle("A. Eukaryotic Phytoplankton")
  
  euk_gradient <- euk_gradient + theme(legend.position = "none") + xlab("") + ylab("") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize)) + ggtitle("B. Eukaryotic Protists")
    
  cyano_gradient <- cyano_gradient + theme(legend.position = "none")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("C. Cyanobacteria")
  
  bact_gradient <- bact_gradient +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("D. Bacteria") +
    ylab("") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
  
  grad_plot <- phyto_gradient + euk_gradient + cyano_gradient + bact_gradient + plot_layout(guides = "collect")
  
  pdf(file = gradient_plot_file, width = 9, height = 7)
  print(grad_plot)
  dev.off()
  
  
}

##### Figure 8: Diversity vs Time #####

fig_8_func <- function(in_phyto = "output/plastid_16s_diffs_div.Rdata",
                            in_euks = "output/euks_hetero_18sv9_diffs_div.Rdata",
                            in_cyano = "output/cyano_16s_diffs_div.Rdata",
                            in_bact = "output/bacteria_m_euks_16s_diffs_div.Rdata",
                            gradient_plot_file = "figures/figure_outline/fig_8.pdf",
                            tsize = 12){
  
  load(in_phyto)
  phyto_gradient <- gradient_plot
  
  load(in_euks)
  euk_gradient <- gradient_plot
  
  load(in_cyano)
  cyano_gradient <- gradient_plot
  
  load(in_bact)
  bact_gradient <- gradient_plot
  
  phyto_gradient <- phyto_gradient + theme(legend.position = "none") + xlab("")  + 
    ylab("Mean Shannon Diversity\nPer Cruise") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.title.y = element_text(size = tsize),
          axis.text.y = element_text(size = tsize)) + ggtitle("A. Eukaryotic Phytoplankton")
  
  euk_gradient <- euk_gradient + theme(legend.position = "none") + xlab("") + ylab("") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize)) + ggtitle("B. Eukaryotic Protists")
  
  cyano_gradient <- cyano_gradient + theme(legend.position = "none")  + 
    ylab("Mean Shannon Diversity\nPer Cruise") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("C. Cyanobacteria")
  
  bact_gradient <- bact_gradient +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("D. Bacteria") +
    ylab("") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
  
  grad_plot <- phyto_gradient + euk_gradient + cyano_gradient + bact_gradient + plot_layout(guides = "collect")
  
  pdf(file = gradient_plot_file, width = 9, height = 7)
  print(grad_plot)
  dev.off()
  
  
}


##### Supplementary Figures #####


##### Suppl Figure 2: SOMs small groups #####

suppl_fig_2_func <- function(in_list = plot_list, file_name = "figures/figure_outline/supp_fig_2.pdf", tsize  = 12){
  
  in_list[[1]]$bp <- in_list[[1]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("A. Prochlorococcus")
  in_list[[2]]$bp <- in_list[[2]]$bp + theme(axis.title.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize),
                                               legend.position = "none") + ggtitle("B. Synechococcus")
  in_list[[3]]$bp <- in_list[[3]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("C. Flavobacteriales")
  in_list[[4]]$bp <- in_list[[4]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("D. Rhodobacterales")
  in_list[[5]]$bp <- in_list[[5]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize)) + ggtitle("E. Sar 11 Clade")
  in_list[[7]]$bp <- in_list[[7]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("F. Diatoms")
  in_list[[8]]$bp <- in_list[[8]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("G. Dinoflagellates")
  in_list[[9]]$bp <- in_list[[9]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("H. Syndiniales")
  in_list[[10]]$bp <- in_list[[10]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("I. Haptophytes")
  in_list[[11]]$bp <- in_list[[11]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("J. Chlorophytes")
  in_list[[12]]$bp <- in_list[[12]]$bp + theme(axis.title.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize),
                                               legend.position = "none") + ggtitle("K. Metazoans")
  
  in_list[[1]]$rp <- in_list[[1]]$rp + theme(legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[2]]$rp <- in_list[[2]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[3]]$rp <- in_list[[3]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[4]]$rp <- in_list[[4]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[5]]$rp <- in_list[[5]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0)) + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  
  in_list[[7]]$rp <- in_list[[7]]$rp + theme(legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[8]]$rp <- in_list[[8]]$rp + theme(axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[9]]$rp <- in_list[[9]]$rp + theme(axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[10]]$rp <- in_list[[9]]$rp + theme(axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[11]]$rp <- in_list[[10]]$rp + theme(axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[12]]$rp <- in_list[[11]]$rp + theme(axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))

  
  
  
  pdf(file = file_name, width = 25, height = 12)
  print(in_list[[1]]$bp + in_list[[2]]$bp + in_list[[3]]$bp + in_list[[4]]$bp + in_list[[5]]$bp + plot_spacer() +
  in_list[[1]]$rp + in_list[[2]]$rp + in_list[[3]]$rp + in_list[[4]]$rp + in_list[[5]]$rp + plot_spacer() +
  in_list[[7]]$bp + in_list[[8]]$bp + in_list[[9]]$bp + in_list[[10]]$bp + in_list[[11]]$bp + in_list[[12]]$bp + 
  in_list[[7]]$rp + in_list[[8]]$rp + in_list[[9]]$rp + in_list[[10]]$rp + in_list[[11]]$rp + in_list[[12]]$rp +
  plot_layout(ncol = 6, guides = "collect"))
  dev.off()
}

##### Suppl Figure 3: Variable AIC small groups #####

full_aic_table_figure(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                        "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                        "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                      in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar Clade", 
                                         "Diatoms",
                                         "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"),
                      minimum_tp = 4, width_plot = 16,
                      figure_name_2 = paste0("figures/figure_outline/supp_fig_3",".pdf"),
                      title_name = "", tsize = 12)

#### Suppl Figure 4: Diversity Maps #####

suppl_fig_4_func <- function(file_name = "figures/figure_outline/supp_fig_4.pdf",
                             in_group_list = c("archaea_16s","bacteria_m_euks_16s", "cyano_16s","plastid_16s",
                                               "euks_hetero_18sv9"),
                             in_group_names = c("Archaea","Bacteria", "Cyanobacteria", "Eukaryotic Phytoplankton\n(Plastids)",
                                                "Eukaryotic\n Protists"), tsize = 12){
  
  
}

##### Suppl Figure 5: Diversity Importance Small Groups #####

full_aic_table_figure_diversity_sign(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                                       "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                       "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                                     in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar Clade", 
                                                        "Diatoms",
                                                        "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"),
                                     figure_name_2 = "figures/figure_outline/supp_fig_5.pdf",
                                     col = 27, width_plot = 16, tsize = 12)

##### Suppl Figure 6: Total Diversity vs Time #####

suppl_fig_6_func <- function(in_phyto = "output/plastid_16s_diffs_div.Rdata",
                       in_euks = "output/euks_hetero_18sv9_diffs_div.Rdata",
                       in_cyano = "output/cyano_16s_diffs_div.Rdata",
                       in_bact = "output/bacteria_m_euks_16s_diffs_div.Rdata",
                       gradient_plot_file = "figures/figure_outline/supp_fig_6.pdf",
                       tsize = 12){
  
  load(in_phyto)
  phyto_gradient2 <- gradient_plot2
  
  load(in_euks)
  euk_gradient2 <- gradient_plot2

  load(in_cyano)
  cyano_gradient2 <- gradient_plot2
  
  load(in_bact)
  bact_gradient2 <- gradient_plot2

  phyto_gradient2 <- phyto_gradient2 + theme(legend.position = "none") + xlab("")  + 
    ylab("Total Shannon Diversity\nPer Cruise") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.title.y = element_text(size = tsize),
          axis.text.y = element_text(size = tsize)) + ggtitle("A. Eukaryotic Phytoplankton")
  
  euk_gradient2 <- euk_gradient2 + theme(legend.position = "none") + xlab("") + ylab("") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize)) + ggtitle("B. Eukaryotic Protists")
  
  cyano_gradient2 <- cyano_gradient2 + theme(legend.position = "none")  + 
    ylab("Total Shannon Diversity\nPer Cruise") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("C. Cyanobacteria")
  
  bact_gradient2 <- bact_gradient2 +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("D. Bacteria") +
    ylab("") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
  
  grad_plot <- phyto_gradient2 + euk_gradient2 + cyano_gradient2 + bact_gradient2 + plot_layout(guides = "collect")
  
  pdf(file = gradient_plot_file, width = 9, height = 7)
  print(grad_plot)
  dev.off()
  
  
}

