library(tidyverse)
library(ncdf4)
library(chron)
library(lattice)
library(RColorBrewer)
library(rworldmap)
library(ggplot2)
library(ggmap)
library(cowplot)
library(reshape2)
library(metR)
library(oce)

# GLOBAL OCEAN GRIDDED L4 SEA SURFACE HEIGHTS AND DERIVED VARIABLES REPROCESSED (1993-ONGOING)
# http://marine.copernicus.eu/services-portfolio/access-to-products/

  # set path and filename
  
  ncpath <- "data/"
  ncname <- "dataset-duacs-rep-global-merged-allsat-93-19" # use full dataset
  ncfname <- paste(ncpath, ncname, ".nc", sep="")
  dname_u <- "ugos"  # note: sea surface height relative to geoid
  dname_v <- "vgos"
  
  ncin <- nc_open(ncfname)
  
  # get lat and long
  
  lon <- ncin$var$ugos$dim[[1]]$vals
  lat <- ncin$var$ugos$dim[[2]]$vals
  
  # get time
  
  time <- ncvar_get(ncin,"time")
  tunits <- ncatt_get(ncin,"time","units")
  
  # Get temperature
  
  u_array <- ncvar_get(ncin,dname_u)
  d_u_units <- ncatt_get(ncin,dname_u,"units")
  
  v_array <- ncvar_get(ncin,dname_v)
  d_v_units <- ncatt_get(ncin,dname_v,"units")
  
  fillvalue_u <- ncatt_get(ncin,dname_u,"_FillValue")
  fillvalue_v <- ncatt_get(ncin,dname_v,"_FillValue")
  
  # Convert Time
  
  tustr <- strsplit(tunits$value, " ")
  tdstr <- strsplit(unlist(tustr)[3], "-")
  tmonth <- as.integer(unlist(tdstr)[2])
  tday <- as.integer(unlist(tdstr)[3])
  tyear <- as.integer(unlist(tdstr)[1])
  ts <- chron(time,origin=c(tmonth, tday, tyear))
  
  u_array[u_array==fillvalue_u$value] <- NA
  v_array[v_array==fillvalue_v$value] <- NA
  
  # getting dims for array
  
  min_lat = 27.5
  max_lat = 39.5
  min_lon = 360-127.5
  max_lon = 360-115.5
  
  lat_grid <- which(lat > min_lat & lat < max_lat)
  lon_grid <-which(lon > min_lon & lon < max_lon)
  
  u_slice <- u_array[lon_grid,lat_grid,]
  v_slice <- v_array[lon_grid,lat_grid,]

  
  # save the slice of data
save(u_slice, v_slice,ts,lat_grid,lon_grid,lon,lat, file = "output/CALCOFI_Copernicus_GSC_Slice.Rdata")

# Calculate Coeff. Var.

load("output/CALCOFI_Copernicus_GSC_Slice.Rdata")

# get coeff of variation u 
grid_mean_u <- matrix(NA,length(lon_grid),length(lat_grid))
grid_sd_u <- matrix(NA,length(lon_grid),length(lat_grid))
grid_coeff_u <- matrix(NA,length(lon_grid),length(lat_grid))

rownames(grid_coeff_u) <- lon[lon_grid]-360
rownames(grid_sd_u) <- lon[lon_grid]-360
rownames(grid_mean_u) <- lon[lon_grid]-360

colnames(grid_coeff_u) <- lat[lat_grid]
colnames(grid_sd_u) <- lat[lat_grid]
colnames(grid_mean_u) <- lat[lat_grid]

# get coeff of variation v
grid_mean_v <- matrix(NA,length(lon_grid),length(lat_grid))
grid_sd_v <- matrix(NA,length(lon_grid),length(lat_grid))
grid_coeff_v <- matrix(NA,length(lon_grid),length(lat_grid))

rownames(grid_coeff_v) <- lon[lon_grid]-360
rownames(grid_sd_v) <- lon[lon_grid]-360
rownames(grid_mean_v) <- lon[lon_grid]-360

colnames(grid_coeff_v) <- lat[lat_grid]
colnames(grid_sd_v) <- lat[lat_grid]
colnames(grid_mean_v) <- lat[lat_grid]

for(i in 1:length(lon_grid)){
  for (j in 1:length(lat_grid)) {
    
    vals_u <- list()
    vals_v <- list()
    
    for(k in 1:1826){
      
      vals_u[[k]] <- u_slice[i,j,k]
      vals_v[[k]] <- v_slice[i,j,k]
      
    }
    
    val_vect_u <- unlist(vals_u)
    val_vect_v <- unlist(vals_v)
    
    mean_val_u <- mean(val_vect_u, na.rm = TRUE)
    sd_val_u <- sd(val_vect_u, na.rm = TRUE)
    
    mean_val_v <- mean(val_vect_v, na.rm = TRUE)
    sd_val_v <- sd(val_vect_v, na.rm = TRUE)
    
    grid_mean_u[i,j] <- mean_val_u
    grid_sd_u[i,j] <- sd_val_u
    grid_coeff_u[i,j] <- sd_val_u/mean_val_u
    
    grid_mean_v[i,j] <- mean_val_v
    grid_sd_v[i,j] <- sd_val_v
    grid_coeff_v[i,j] <- sd_val_v/mean_val_v
    
    
  }
}

coeff_table_u <- melt(grid_coeff_u)
mean_table_u <- melt(grid_mean_u)
sd_table_u <- melt(grid_sd_u)

coeff_table_v <- melt(grid_coeff_v)
mean_table_v <- melt(grid_mean_v)
sd_table_v <- melt(grid_sd_v)

colnames(coeff_table_u) <- c("lon", "lat", "coeff_var")
colnames(mean_table_u) <- c("lon", "lat", "Mean")
colnames(sd_table_u) <- c("lon", "lat", "SD")

colnames(coeff_table_v) <- c("lon", "lat", "coeff_var")
colnames(mean_table_v) <- c("lon", "lat", "Mean")
colnames(sd_table_v) <- c("lon", "lat", "SD")

coeff_table_u <- coeff_table_u[which(!is.na(coeff_table_u$coeff_var)),]
mean_table_u <- mean_table_u[which(!is.na(mean_table_u$Mean)),]
sd_table_u <- mean_table_u[which(!is.na(mean_table_u$SD)),]

coeff_table_v <- coeff_table_v[which(!is.na(coeff_table_v$coeff_var)),]
mean_table_v <- mean_table_v[which(!is.na(mean_table_v$Mean)),]
sd_table_v <- mean_table_v[which(!is.na(mean_table_v$SD)),]

save(coeff_table_u = coeff_table_u, mean_table_u = mean_table_u, sd_table_u = sd_table_u,
     coeff_table_v = coeff_table_v, mean_table_v = mean_table_v, sd_table_v = sd_table_v,
     file = "output/CALCOFI_uv_current_tables.Rdata")

# Figures

colnames(mean_table_u)[3] <- "Mean_U"
colnames(mean_table_v)[3] <- "Mean_V"

map <- map_data("world")   

uv_table <- inner_join(mean_table_u, mean_table_v, by = c("lon", "lat"))

uv_table$vel = sqrt(uv_table$Mean_U^2+uv_table$Mean_V^2)

save(uv_table, file = "output/uv_velocity_table.Rdata")

# https://semba-blog.netlify.com/03/20/2019/plotting-streamlines-of-surface-current-with-ggplot2-and-metr-package/

current_plot <- ggplot() +
  geom_raster(data = uv_table, aes(x = lon, y = lat, fill = vel), interpolate = FALSE) + 
  metR::geom_vector(data = uv_table,aes(x = lon, y = lat, dx = Mean_U, dy = Mean_V), 
                    arrow.angle = 15, arrow.type = "open", arrow.length = unit(0.5, "inches"), 
                    pivot = 0,preserve.dir = TRUE, direction = "ccw",
                    min.mag = 0, show.legend = NA) +
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  scale_fill_gradientn(name = "Speed\n(m/s)",colours = oceColorsVelocity(120), 
                       limits = c(0,0.1), breaks =seq(0,0.1,.02), oob = scales::squish) +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.height = unit(1.4, "cm"), 
        legend.background = element_blank(),
        axis.text = element_text(size = 12, colour = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_mag(max = 0.1, name = "Speed", max_size = 0.75, guide = "none") +
  labs(x = "Longitude", y = "Latitude") + ggtitle("Mean Surface Current (1993-2018)")

pdf(file = paste0("figures/current_plot_",Sys.Date(),".pdf"), width = 8, height = 8)
current_plot
dev.off()
