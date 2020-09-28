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
library(extrafont)

fig_1_func <- function(#in_vel = "output/uv_velocity_table.Rdata",
                       #in_bath = "output/CALCOFI_bathymetry_table.Rdata",
                       in_temp = "output/CALCOFI_temp_tables.Rdata",
                       in_cyano = "output/cyano_16s_map.Rdata",
                       fig_name = "figures/grant_fig.pdf",
                       tsize = 12, psize = 6, tsize2 = 16, ffamily = "Times New Roman"){
  
  map <- map_data("world")    
  
  # load(in_vel)
  # load(in_bath)
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
          plot.title = element_text(size = tsize2, family = ffamily), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize, family = ffamily),
          legend.text = element_text(size = tsize, family = ffamily),
          axis.title = element_text(size = tsize, family = ffamily),
          legend.title = element_text(size = tsize, family = ffamily)) +
    ggtitle("a. Samples per station")
  
  
  # reduced_uv_table <- uv_table[seq(1,1633,6),]
  # elevation_table$value <- abs(elevation_table$value)
  # 
  # vel_bath <- ggplot() + 
  #   geom_raster(data = elevation_table, aes(x = lon, y = lat, fill = value), interpolate = FALSE) +
  #   scale_fill_gradient(low = "darkblue", high = "cyan", name = "Depth (m)", trans = 'reverse') +
  #   metR::geom_vector(data = reduced_uv_table, aes(x = lon, y = lat, dx = Mean_U, dy = Mean_V), 
  #                     arrow.angle = 15, arrow.type = "open", arrow.length = unit(0.5, "inches"), 
  #                     pivot = 0,preserve.dir = TRUE, direction = "ccw",
  #                     min.mag = 0, show.legend = NA, color = "white") +
  #   geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
  #   coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  #   labs(x = "Longitude", y = "Latitude", mag = "Speed (m/s)", color = "Depth (m)")  +
  #   theme(panel.background = element_blank(),
  #         legend.key = element_rect(fill = "grey"),
  #         panel.border = element_rect(fill = NA, colour = "black", linetype = "solid"),
  #         plot.title = element_text(), axis.line = element_blank(),
  #         legend.justification=c(1,1), 
  #         legend.position=c(0.97, 0.97),
  #         legend.background = element_rect(fill = "white", color = "black"),
  #         axis.text = element_text(size = tsize),
  #         legend.text = element_text(size = tsize),
  #         axis.title = element_text(size = tsize),
  #         legend.box = "horizontal",
  #         legend.title = element_text(size = tsize)) + 
  #   guides(fill = guide_legend(order = 2),mag = guide_legend(order = 1)) +
  #   ggtitle("B. Mean Geostrophic Currents")
  
  
  sst <- ggplot() + 
    geom_tile(data = coeff_table, aes(x = lon, y = lat, fill = coeff_var), width =0.26, height = 0.26) +
    scale_fill_gradient2(name = "Coeff. Var SST", low = "darkblue", mid = "white", high = "darkred", limits = c(0.09,0.12), oob = squish, midpoint = 0.1066851) + geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    ggtitle("d. Coeff. var. SST") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(size = tsize2), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize))
  
  sst_mean <- ggplot() + 
    geom_tile(data = mean_table, aes(x = lon, y = lat, fill = Mean), width =0.26, height = 0.26) +
    scale_fill_gradient2(name = "SST Mean (°C)", low = "darkblue", mid = "white", high = "darkred", limits = c(15,18), oob = squish, midpoint = 16.5) +geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize2), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) + 
    ggtitle("c. Mean SST (°C)")
  
  nc_depth <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "NC_mean"),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean Nitracline\nDepth (m)", low = "darkblue", high = "cyan", trans = 'reverse') +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize2), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("b. Mean nitracline depth")
  
  
  stations <- stations + theme(axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank())
  
  nc_depth <- nc_depth + theme(axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),
                               axis.title.y=element_blank(),
                               axis.text.y=element_blank(),
                               axis.ticks.y=element_blank())
  
  sst <- sst + theme(axis.title.y=element_blank(),
                               axis.text.y=element_blank(),
                               axis.ticks.y=element_blank())
  
  patch <- stations  + nc_depth + sst_mean + sst
  
  pdf(fig_name, width = 12, height = 12)
  print(patch)
  dev.off()
  
}