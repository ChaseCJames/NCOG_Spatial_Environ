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
library(ggmap)
library(patchwork)
library(sp)
library(gstat)

rgb_map <- function(in_dat = "output/cyano_16s_dissimilar.Rdata",
                    in_map = "output/cyano_16s_map.Rdata",
                    in_full = "output/cyano_16s_full_data.Rdata",
                    in_bound = "output/boundary_grid.Rdata",
                    color_pal = "output/2dcolorpal.Rdata",
                    tsize = 12, out_file = "figures/cyano_pca_plots.pdf"){
  
  load(in_dat)
  load(in_map)
  load(in_full)
  load(in_bound)
  load(color_pal)
    
  splits <- strsplit(rownames(dissimilar), split = "_")
  depths <- vector()
  
  for (i in 1:length(splits)) {depths[i] <- splits[[i]][4]}
  
  depths <- as.numeric(depths)
  
  dissimilar <- dissimilar[which(depths < 15),which(depths < 15)]
  
  pca_out <- prcomp(dissimilar)
  # pca_out <- metaMDS(dissimilar,)
  
  pca_rgb <- pca_out$x[,1:3]
  # pca_rgb <- as.data.frame(pca_out$points)
  # pca_rgb$MDS2_1 <- pca_rgb$MDS2
  
  pca_rgb <- apply(pca_rgb,2,function(x) (round(((x-min(x))/(max(x)-min(x))*255),digits = 0)))
  
  # red blue pallete
  pca_rgb[,3] <- pca_rgb[,2]
  pca_rgb[,2] <- 0

  rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)
  
  # pca_hex <- apply(pca_rgb, 1, function(x) rgb2hex(r = x[1], g = x[2], b = x[3]))
  
  full_dat$red <- pca_rgb[match(paste0("X",full_dat$Sample.Name),rownames(pca_rgb)),1]
  full_dat$green <- pca_rgb[match(paste0("X",full_dat$Sample.Name),rownames(pca_rgb)),2]
  full_dat$blue <- pca_rgb[match(paste0("X",full_dat$Sample.Name),rownames(pca_rgb)),3]
  
  full_dat <- full_dat[which(!is.na(full_dat$red)),]
  full_dat <- full_dat[which(!is.na(full_dat$dist_to_coast)),]
  
  full_dat <- full_dat %>%
    filter(as.numeric(substr(Sta_ID,1,3)) >= 76)
  
  map <- map_data("world")   
  
  cruises <- unique(full_dat$Cruise)
  
  map_list <- list()
  
  # set up grid
  
  full_poly <- SpatialPolygons(list(Polygons(list(Polygon(coordinates(boundary_mat[,4:5]))), ID=1)))
  
  grd <- makegrid(full_poly, n = 1000)
  colnames(grd) <- c('x','y')
  
  grd_pts <- SpatialPoints(coords = grd, 
                           proj4string=CRS(proj4string(full_poly)))
  
  grd_pts_in <- grd_pts[full_poly, ]
  
  gdf <- as.data.frame(coordinates(grd_pts_in)) 
  
  coordinates(gdf) <- ~ x + y
  
  krig_dat <- full_dat
  
  coordinates(krig_dat) <- ~ Lon_Dec + Lat_Dec
  
  # red.vgm <- variogram(red~1, krig_dat) # calculates sample variogram values 
  # red.fit <- fit.variogram(red.vgm, model=vgm(c("Gau"))) # fit model
  # 
  # blue.vgm <- variogram(blue~1, krig_dat) # calculates sample variogram values 
  # blue.fit <- fit.variogram(blue.vgm, model=vgm(c("Gau"))) # fit model
  # 
  # green.vgm <- variogram(green~1, krig_dat) # calculates sample variogram values 
  # green.fit <- fit.variogram(green.vgm, model=vgm(c("Gau"))) # fit model
  
  for (i in 1:length(cruises)) {
    
    cruise_dat <- full_dat %>%
      filter(Cruise == cruises[i])
    
    cruise_krig <- cruise_dat
    coordinates(cruise_krig) <- ~ Lon_Dec + Lat_Dec
    
    # red.kriged <- krige(red ~ 1, cruise_krig, gdf, model=red.fit)
    # blue.kriged <- krige(blue ~ 1, cruise_krig, gdf, model=blue.fit)
    # green.kriged <- krige(green ~ 1, cruise_krig, gdf, model=green.fit)
    
    red.idw <- idw(formula=red ~ 1, locations=cruise_krig, newdata=gdf)
    red_df <- as.data.frame(red.idw)
    blue.idw <- idw(formula=blue ~ 1, locations=cruise_krig, newdata=gdf)
    blue_df <- as.data.frame(blue.idw)
    green.idw <- idw(formula=green ~ 1, locations=cruise_krig, newdata=gdf)
    green_df <- as.data.frame(green.idw)
    
    
    # red_df <- as.data.frame(red.kriged)
    # blue_df <- as.data.frame(blue.kriged)
    # green_df <- as.data.frame(green.kriged)
  
    red_df <- red_df[,1:3]
    colnames(red_df) <- c("Lon", "Lat", "Red")  
    red_df$Blue <- blue_df$var1.pred
    red_df$Green <- green_df$var1.pred
    red_df$Green <- 0
    
    pca_hex <- apply(red_df, 1, function(x) rgb2hex(r = x[3], g = x[5], b = x[4]))
    red_df$Hex <- pca_hex
    
    plot <- ggplot() + 
      geom_tile(data = red_df, aes(x = Lon, y = Lat, fill = Hex), show.legend = FALSE) +
      scale_fill_manual(values = red_df$Hex) +
      geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
      coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
      xlab("Longitude") + ylab("Latitude") + 
      geom_point(data = cruise_dat, aes(x = Lon_Dec, y = Lat_Dec), color = "black",fill = "white", shape = 21, size = 3) +
      ggtitle(paste0(substr(cruises[i],5,6),"-",substr(cruises[i],1,4))) +
      theme(legend.title = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
            plot.title = element_text(size = tsize), axis.line = element_blank(),
            axis.text = element_text(size = tsize),
            legend.text = element_text(size = tsize),
            axis.title = element_text(size = tsize))
    
    # print(plot)
    map_list[[i]] <- plot
    
  }
  
  out_plot <- map_list[[1]] + map_list[[2]] + map_list[[3]] + map_list[[4]] +
    map_list[[5]] + map_list[[6]] + map_list[[7]] + map_list[[8]] +
    map_list[[9]] + map_list[[10]] + map_list[[11]] + map_list[[12]] +
    map_list[[13]] + map_list[[14]] + map_list[[15]] + map_list[[16]] +
    map_list[[17]] + map_list[[18]] + map_list[[19]] + map_list[[20]] +
    map_list[[21]] + map_list[[22]] + map_list[[23]] + map_list[[24]] +
    plot_spacer() + out_pal + plot_spacer() +
    plot_layout(nrow =4, ncol = 7, byrow = FALSE)
    
  pdf(file = out_file, width = 27, height = 15)
  print(out_plot)
  dev.off()
  
  }


rgb_map(in_dat = "output/cyano_16s_dissimilar.Rdata",
        in_map = "output/cyano_16s_map.Rdata",
        in_full = "output/cyano_16s_full_data.Rdata",
        tsize = 12, out_file = "figures/cyano_pca_plots.pdf")

rgb_map(in_dat = "output/bacteria_m_euks_16s_dissimilar.Rdata",
        in_map = "output/bacteria_m_euks_16s_map.Rdata",
        in_full = "output/bacteria_m_euks_16s_full_data.Rdata",
        tsize = 12, out_file = "figures/bact_pca_plots.pdf")

rgb_map(in_dat = "output/euks_auto_18sv9_dissimilar.Rdata",
        in_map = "output/euks_auto_18sv9_map.Rdata",
        in_full = "output/euks_auto_18sv9_full_data.Rdata",
        tsize = 12, out_file = "figures/euks_auto_pca_plots.pdf")

rgb_map(in_dat = "output/euks_hetero_18sv9_dissimilar.Rdata",
        in_map = "output/euks_hetero_18sv9_map.Rdata",
        in_full = "output/euks_hetero_18sv9_full_data.Rdata",
        tsize = 12, out_file = "figures/euks_hetero_pca_plots.pdf")

rgb_map(in_dat = "output/archaea_16s_dissimilar.Rdata",
        in_map = "output/archaea_16s_map.Rdata",
        in_full = "output/archaea_16s_full_data.Rdata",
        tsize = 12, out_file = "figures/archaea_pca_plots.pdf")








