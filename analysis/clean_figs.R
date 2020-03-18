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

#### Physical Conditions Figure ####

physical_fig <- function(in_vel = "output/uv_velocity_table.Rdata",
                         in_bath = "output/CALCOFI_bathymetry_table.Rdata",
                         in_temp = "output/CALCOFI_temp_tables.Rdata",
                         in_cyano = "output/cyano_16s_map.Rdata",
                         fig_name = "figures/summary_physical_figure_clean.pdf"){
  
  map <- map_data("world")    
  
  load(in_vel)
  load(in_bath)
  load(in_temp)
  load(in_cyano)
  
  
  stations <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "n_samps"), color = "black", size =6, stroke = 0.1, shape = 21) +
    ggtitle("CalCOFI NCOG Stations") + scale_fill_gradient(name = "# of Samples", low = "white", high = "red") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  reduced_uv_table <- uv_table[seq(1,1633,6),]
  
  vel_bath <- ggplot() + 
    geom_raster(data = elevation_table, aes(x = lon, y = lat, fill = value), interpolate = FALSE) +
    scale_fill_gradient(low = "darkblue", high = "cyan", name = "Depth (m)",
                        limits = c(-5000,0), oob = scales::squish) +
    metR::geom_vector(data = reduced_uv_table, aes(x = lon, y = lat, dx = Mean_U, dy = Mean_V), 
                      arrow.angle = 15, arrow.type = "open", arrow.length = unit(0.5, "inches"), 
                      pivot = 0,preserve.dir = TRUE, direction = "ccw",
                      min.mag = 0, show.legend = NA, color = "white") +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    theme(legend.position = "right",
          panel.background = element_blank(),
          legend.key = element_rect(fill = "grey"),
          panel.border = element_rect(fill = NA, colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank()) +
    scale_mag(max = 0.1, name = "Speed (m/s)", max_size = 0.75) +
    labs(x = "Longitude", y = "Latitude", color = "Depth (m)") + ggtitle("Mean Geostrophic Current Velocity") 
  
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
    ggtitle("SST Mean (°C)") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  
  pdf(fig_name, width = 12, height = 12)
  print(plot_grid(stations, vel_bath,
                  sst_mean, sst, ncol = 2, nrow = 2, align = "hv", labels = c("A.", "B.", "C.", "D.")))
  dev.off()
  
}

#### SOM Figure ####

som_figure <- function(map_file = "output/bacteria_m_euks_16s_map.Rdata",
                       figure_name = paste0("figures/som_maps/bacteria_16s_som_",Sys.Date(),".pdf"),
                       main = "16s Bacteria", cluster1 = "Nearshore", cluster2 = "Offshore"){
  

  
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
    ggtitle(paste0("B. ",cluster1," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())

  if(clust1 == 1){p1 <- p1 + geom_point(aes(x = wt_1@coords[1], y = wt_1@coords[2]), color = "blue", size = 5, pch = 10)}
  if(clust1 == 2){p1 <- p1 + geom_point(aes(x = wt_2@coords[1], y = wt_2@coords[2]), color = "blue", size = 5, pch = 10)}
  
  p2 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",clust2)), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "darkblue", limits = c(0,1)) +
    ggtitle(paste0("A. ", cluster2," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  if(clust2 == 1){p2 <- p2 + geom_point(aes(x = wt_1@coords[1], y = wt_1@coords[2]), color = "red", size = 5, pch = 10)}
  if(clust2 == 2){p2 <- p2 + geom_point(aes(x = wt_2@coords[1], y = wt_2@coords[2]), color = "red", size = 5, pch = 10)}
  
  title <- ggdraw() + draw_label(main, fontface='bold')
  a <- plot_grid(title,plot_grid(p2,p1), nrow = 2,
            rel_heights = c(0.1,1))
  
  pdf(file = figure_name, width = 6, height = 4)
  print(a)
  dev.off()
  
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
  
  
  
  pdf(file = paste0(figure_name,var,"_",Sys.Date(),".pdf"), width = 5, height = 3)
  print(reg_plot)
  dev.off()
  
  
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
    geom_boxplot(fill = "white", width = 0.1) + labs(fill = "") +
    scale_fill_manual(values = c("red","blue")) + xlab("") +
    ylab("Bray-Curtis\nDissimilarity") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
    ggtitle(main)
  
  pdf(file = figure_name, width = 4, height = 5)
  print(bc_plot)
  dev.off()
  
  
}

### AIC Tables ####

aic_table_func <- function(som_maps = cyano_plots){  
  
  som_glm <- som_maps[,c(1:3,7:25,29:32)]
  
  model_AIC <- matrix(NA,21,2)
  
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
  
  glm_mean_sla <- glm(som_1 ~  sla_mean, data = som_glm, family = binomial)
  ms_sum <- summary(glm_mean_sla)
  model_AIC[7,2] <- ms_sum$aic
  model_AIC[7,1] <- "Mean SLA"
  
  glm_coeff_sla <- glm(som_1 ~  sla_coeff, data = som_glm, family = binomial)
  cs_sum <- summary(glm_coeff_sla)
  model_AIC[16,2] <- cs_sum$aic
  model_AIC[16,1] <- "Coeff. Var. SLA"
  
  # MLD
  
  glm_mean_mld <- glm(som_1 ~  MLD_mean, data = som_glm, family = binomial)
  ms_sum <- summary(glm_mean_mld)
  model_AIC[8,2] <- ms_sum$aic
  model_AIC[8,1] <- "Mean MLD"
  
  glm_coeff_mld <- glm(som_1 ~  MLD_coeff, data = som_glm, family = binomial)
  cs_sum <- summary(glm_coeff_mld)
  model_AIC[17,2] <- cs_sum$aic
  model_AIC[17,1] <- "Coeff. Var. MLD"
  
  # NC Depth
  
  glm_mean_nc <- glm(som_1 ~  NC_mean, data = som_glm, family = binomial)
  ms_sum <- summary(glm_mean_nc)
  model_AIC[9,2] <- ms_sum$aic
  model_AIC[9,1] <- "Mean NCD"
  
  glm_coeff_nc <- glm(som_1 ~  NC_coeff, data = som_glm, family = binomial)
  cs_sum <- summary(glm_coeff_nc)
  model_AIC[18,2] <- cs_sum$aic
  model_AIC[18,1] <- "Coeff. Var. NCD"
  
  # Distance to Coast
  
  glm_mean_dc <- glm(som_1 ~  Dist_mean, data = som_glm, family = binomial)
  ms_sum <- summary(glm_mean_dc)
  model_AIC[19,2] <- ms_sum$aic
  model_AIC[19,1] <- "Distance to Coast"
  
  # Chlorophyll
  
  glm_mean_chl <- glm(som_1 ~  Chl_mean, data = som_glm, family = binomial)
  ms_sum <- summary(glm_mean_chl)
  model_AIC[20,2] <- ms_sum$aic
  model_AIC[20,1] <- "Mean Chl-a"
  
  glm_coeff_chl <- glm(som_1 ~  Chl_coeff, data = som_glm, family = binomial)
  cs_sum <- summary(glm_coeff_chl)
  model_AIC[21,2] <- cs_sum$aic
  model_AIC[21,1] <- "Coeff. Var. Chl-a"
  
  
  return(model_AIC)
  
}

aic_table_func_diveristy <- function(som_maps = cyano_plots, i = 2){  
  
  som_glm <- som_maps[,c(1,26:28,7:25,29:32)]
  
  colnames(som_glm)[i] <- "response"
  
  model_AIC <- matrix(NA,21,2)
  model_p_val <- matrix(NA,21,2)
  
  
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
  
  glm_mean_sla <- glm(response ~  sla_mean, data = som_glm)
  ms_sum <- summary(glm_mean_sla)
  model_AIC[7,2] <- ms_sum$aic
  model_AIC[7,1] <- "Mean SLA"
  
  glm_coeff_sla <- glm(response ~  sla_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_sla)
  model_AIC[16,2] <- cs_sum$aic
  model_AIC[16,1] <- "Coeff. Var. SLA"
  
  # MLD
  
  glm_mean_mld <- glm(response ~  MLD_mean, data = som_glm)
  ms_sum <- summary(glm_mean_mld)
  model_AIC[8,2] <- ms_sum$aic
  model_AIC[8,1] <- "Mean MLD"
  
  glm_coeff_mld <- glm(response ~  MLD_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_mld)
  model_AIC[17,2] <- cs_sum$aic
  model_AIC[17,1] <- "Coeff. Var. MLD"
  
  # NC Depth
  
  glm_mean_nc <- glm(response ~  NC_mean, data = som_glm)
  ms_sum <- summary(glm_mean_nc)
  model_AIC[9,2] <- ms_sum$aic
  model_AIC[9,1] <- "Mean NCD"
  
  glm_coeff_nc <- glm(response ~  NC_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_nc)
  model_AIC[18,2] <- cs_sum$aic
  model_AIC[18,1] <- "Coeff. Var. NCD"
  
  # Distance to Coast
  
  glm_mean_dc <- glm(response ~  Dist_mean, data = som_glm)
  ms_sum <- summary(glm_mean_dc)
  model_AIC[19,2] <- ms_sum$aic
  model_AIC[19,1] <- "Distance to Coast"
  
  # Chlorophyll
  
  glm_mean_chl <- glm(response ~  Chl_mean, data = som_glm)
  ms_sum <- summary(glm_mean_chl)
  model_AIC[20,2] <- ms_sum$aic
  model_AIC[20,1] <- "Mean Chl-a"
  
  glm_coeff_chl <- glm(response ~  Chl_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_chl)
  model_AIC[21,2] <- cs_sum$aic
  model_AIC[21,1] <- "Coeff. Var. Chl-a"
  
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
                                  minimum_tp = 8, width_plot = 15,
                                  figure_name = paste0("figures/full_aic_table_logit_",Sys.Date(),".pdf"),
                                  figure_name_2 = paste0("figures/full_aic_plot_logit_",Sys.Date(),".pdf"),
                                  title_name = "Variable Importance"){
  # load data

  map_list <- list()
  
  for (i in 1:length(in_group_list)) {
    
    load(paste0("output/",in_group_list[i], "_map.Rdata"))
    som_maps <- som_maps[which(som_maps$n_samps > (minimum_tp-1)),]
    AIC_table <- aic_table_func(som_maps = som_maps)
    AIC_table <- as.data.frame(AIC_table)
    colnames(AIC_table) <- c("Variables","AIC")
    AIC_table <- AIC_table[c(1,3,4,5,6,20,9,10,12,13,14,15,21,18,19),]
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
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 0)) + ggtitle(title_name) +
          scale_size_continuous(range = c(1,18)) + xlab("") + ylab(""))
  dev.off()
  
  
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
    AIC_table <- AIC_table[c(1,3,4,5,6,20,9,10,12,13,14,15,21,18,19),]
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


###### Outputs #####

#### Physical Conditions Figure ####

physical_fig()

#### SOM Maps ####

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
  
  som_figure(map_file = paste0("output/",in_group_list[i],"_map.Rdata"),
             figure_name = paste0("figures/som_maps/", in_group_list[i],"_",Sys.Date(),".pdf"),
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

full_aic_table_figure(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                                        "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                        "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                      in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar Clade", "Archaea",
                                         "Diatoms",
                                         "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"),
                      minimum_tp = 8, width_plot = 18,
                      figure_name_2 = paste0("figures/aic_figures/small_group_aic_plot_logit",".pdf"),
                      title_name = "Variable Importance")

full_aic_table_figure(in_group_list = c("bacteria_m_euks_16s", "cyano_16s","plastid_16s",
                                        "euks_hetero_18sv9"),
                      in_group_names = c("Bacteria", "Cyanobacteria", "Eukaryotic Phytoplankton\n(Plastids)",
                                         "Eukaryotic\n Protists"),
                      minimum_tp = 8, width_plot = 8,
                      figure_name_2 = paste0("figures/aic_figures/big_group_aic_plot_logit",".pdf"),
                      title_name = "Variable Importance")



# Diversity

full_aic_table_figure_diversity(in_group_list = c("cyano_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                                                  "plastid_16s", "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                  "hapto_18sv9", "metazoa_18sv9"),
                                in_group_names = c("Cyanobacteria", "Flavobacteriales","Rhodobacterales", "Sar Clade", "Archaea",
                                                   "Eukaryotic\n Plastids", "Diatoms",
                                                   "Dinoflagellates", "Syndiniales", "Haptophytes", "Metazoans"),
                                figure_name_2 = paste0("figures/aic_figures/small_group_aic_plot_logit_even",".pdf"),
                                title_name = "Variable Importance Evenness", # col 2 = even, col 3 = shan col 4 = rich
                                col = 2, color_fill = "purple")

full_aic_table_figure_diversity(in_group_list = c("cyano_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                                                  "plastid_16s", "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                  "hapto_18sv9", "metazoa_18sv9"),
                                in_group_names = c("Cyanobacteria", "Flavobacteriales","Rhodobacterales", "Sar Clade", "Archaea",
                                                   "Eukaryotic\n Plastids", "Diatoms",
                                                   "Dinoflagellates", "Syndiniales", "Haptophytes", "Metazoans"),
                                figure_name_2 = paste0("figures/aic_figures/small_group_aic_plot_logit_shannon",".pdf"),
                                title_name = "Variable Importance Shannon", # col 2 = even, col 3 = shan col 4 = rich
                                col = 3, color_fill = "red")

full_aic_table_figure_diversity(in_group_list = c("cyano_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                                                  "plastid_16s", "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                  "hapto_18sv9", "metazoa_18sv9"),
                                in_group_names = c("Cyanobacteria", "Flavobacteriales","Rhodobacterales", "Sar Clade", "Archaea",
                                                   "Eukaryotic\n Plastids", "Diatoms",
                                                   "Dinoflagellates", "Syndiniales", "Haptophytes", "Metazoans"),
                                figure_name_2 = paste0("figures/aic_figures/small_group_aic_plot_logit_rich",".pdf"),
                                title_name = "Variable Importance Richness", # col 2 = even, col 3 = shan col 4 = rich
                                col = 4, color_fill = "blue")

# Basic Groups

full_aic_table_figure_diversity(in_group_list = c("bacteria_m_euks_16s", "cyano_16s","plastid_16s",
                                                   "euks_hetero_18sv9"),
                                in_group_names = c("Bacteria", "Cyanobacteria", "Eukaryotic Phytoplankton\n(Plastids)",
                                                   "Eukaryotic\n Protists"),
                                figure_name_2 = paste0("figures/aic_figures/big_group_aic_plot_logit_even",".pdf"),
                                title_name = "Variable Importance Evenness", # col 2 = even, col 3 = shan col 4 = rich
                                col = 2, color_fill = "purple", width_plot = 8)

full_aic_table_figure_diversity(in_group_list = c("bacteria_m_euks_16s", "cyano_16s","plastid_16s",
                                                   "euks_hetero_18sv9"),
                                in_group_names = c("Bacteria", "Cyanobacteria", "Eukaryotic Phytoplankton\n(Plastids)",
                                                   "Eukaryotic\n Protists"),
                                figure_name_2 = paste0("figures/aic_figures/big_group_aic_plot_logit_shannon",".pdf"),
                                title_name = "Variable Importance\nShannon Diversity", # col 2 = even, col 3 = shan col 4 = rich
                                col = 3, color_fill = "red", width_plot = 8)

full_aic_table_figure_diversity(in_group_list = c("bacteria_m_euks_16s", "cyano_16s","plastid_16s",
                                                  "euks_hetero_18sv9"),
                                in_group_names = c("Bacteria", "Cyanobacteria", "Eukaryotic Phytoplankton\n(Plastids)",
                                                   "Eukaryotic\n Protists"),
                                figure_name_2 = paste0("figures/aic_figures/big_group_aic_plot_logit_rich",".pdf"),
                                title_name = "Variable Importance Richness", # col 2 = even, col 3 = shan col 4 = rich
                                col = 4, color_fill = "blue", width_plot = 8)



