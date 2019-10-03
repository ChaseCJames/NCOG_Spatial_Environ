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

# NOAA OI SST V2 High Resolution Dataset
# https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.highres.html 

  # set path and filename
  
  ncpath <- "data/"
  ncname <- "ssh_rel_geoid"
  ncfname <- paste(ncpath, ncname, ".nc", sep="")
  dname <- "sshg"  # note: sea surface height relative to geoid
  
  ncin <- nc_open(ncfname)
  
  # get lat and long
  
  lon <- ncvar_get(ncin,"lon")
  lat <- ncvar_get(ncin,"lat")
  
  # get time
  
  time <- ncvar_get(ncin,"time")
  tunits <- ncatt_get(ncin,"time","units")
  
  # Get temperature
  
  ssh_array <- ncvar_get(ncin,dname)
  dunits <- ncatt_get(ncin,dname,"units")
  fillvalue <- ncatt_get(ncin,dname,"_FillValue")
  
  
  # Convert Time
  
  tustr <- strsplit(tunits$value, " ")
  tdstr <- strsplit(unlist(tustr)[3], "-")
  tmonth <- as.integer(unlist(tdstr)[2])
  tday <- as.integer(unlist(tdstr)[3])
  tyear <- as.integer(unlist(tdstr)[1])
  ts <- chron(time,origin=c(tmonth, tday, tyear))
  
  ssh_array[ssh_array==fillvalue$value] <- NA
  
  # getting dims for array
  
  min_lat = 27.5
  max_lat = 39.5
  min_lon = 360-127.5
  max_lon = 360-115.5
  
  lat_grid <- which(lat > min_lat & lat < max_lat)
  lon_grid <-which(lon > min_lon & lon < max_lon)
  
  ssh_slice <- ssh_array[lon_grid,lat_grid,]

# save the slice of data
save(ssh_slice,file = "output/CALCOFI_GODAS_SSHG_Slice.Rdata")

# Calculate Coeff. Var.

load("output/CALCOFI_GODAS_SSHG_Slice.Rdata")

# get coeff of variation
grid_mean <- matrix(NA,length(lon_grid),length(lat_grid))
grid_sd <- matrix(NA,length(lon_grid),length(lat_grid))
grid_coeff <- matrix(NA,length(lon_grid),length(lat_grid))

rownames(grid_coeff) <- lon[lon_grid]-360
rownames(grid_sd) <- lon[lon_grid]-360
rownames(grid_mean) <- lon[lon_grid]-360

colnames(grid_coeff) <- lat[lat_grid]
colnames(grid_sd) <- lat[lat_grid]
colnames(grid_mean) <- lat[lat_grid]

for(i in 1:length(lon_grid)){
  for (j in 1:length(lat_grid)) {
    
    vals <- list()
    
    for(k in 1:(12*5)){
      
      vals[[k]] <- ssh_slice[i,j,k]
      
    }
    
    val_vect <- unlist(vals)
    
    mean_val <- mean(val_vect, na.rm = TRUE)
    sd_val <- sd(val_vect, na.rm = TRUE)
    
    grid_mean[i,j] <- mean_val
    grid_sd[i,j] <- sd_val
    grid_coeff[i,j] <- sd_val/mean_val
    
    
  }
}

coeff_table <- melt(grid_coeff)
mean_table <- melt(grid_mean)
sd_table <- melt(grid_sd)

colnames(coeff_table) <- c("lon", "lat", "coeff_var")
colnames(mean_table) <- c("lon", "lat", "Mean")
colnames(sd_table) <- c("lon", "lat", "SD")

map <- map_data("world")    

coeff_table <- coeff_table[which(!is.na(coeff_table$coeff_var)),]
mean_table <- mean_table[which(!is.na(mean_table$Mean)),]

save(coeff_table = coeff_table, mean_table = mean_table, sd_table = sd_table, file = "output/CALCOFI_sshg_tables.Rdata")

# Figures

sshg <- ggplot() + 
  geom_tile(data = coeff_table, aes(x = lon, y = lat, fill = coeff_var), width = 1, height = 1) +
  scale_fill_viridis(limits = c(0.12,1.2), oob = squish) +geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("A. SSHG Coeff. Var.") +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
        title = element_text(hjust = 0.5), axis.line = element_blank())

sshg_mean <- ggplot() + 
  geom_tile(data = mean_table, aes(x = lon, y = lat, fill = Mean), width =1, height = 1) +
  scale_fill_viridis(limits = c(0.02,0.28), oob = squish, option = "A", direction = -1) + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("B. SSHG Mean") +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
        title = element_text(hjust = 0.5), axis.line = element_blank())








