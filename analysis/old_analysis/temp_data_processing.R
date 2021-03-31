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

year_list <- 1991:2019

tmp_slice_list <- list()

ts_list <- vector()

for (i in 1:length(year_list)) {
  
  print(i)

# set path and filename

ncpath <- "data/"
ncname <- paste0("sst.day.mean.",year_list[i])
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "sst"  # note: tmp means temperature (not temporary)

ncin <- nc_open(ncfname)

# get lat and long

lon <- ncvar_get(ncin,"lon")
lat <- ncvar_get(ncin,"lat")

# get time

time <- ncvar_get(ncin,"time")
tunits <- ncatt_get(ncin,"time","units")

# Get temperature

tmp_array <- ncvar_get(ncin,dname)
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")


# Convert Time

tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth <- as.integer(unlist(tdstr)[2])
tday <- as.integer(unlist(tdstr)[3])
tyear <- as.integer(unlist(tdstr)[1])
ts <- chron(time,origin=c(tmonth, tday, tyear))

ts_list <- c(ts_list,ts)

tmp_array[tmp_array==fillvalue$value] <- NA

# getting dims for array

min_lat = 27.5
max_lat = 39.5
min_lon = 360-127.5
max_lon = 360-115.5

lat_grid <- which(lat > min_lat & lat < max_lat)
lon_grid <-which(lon > min_lon & lon < max_lon)

tmp_slice <- tmp_array[lon_grid,lat_grid,]

tmp_slice_list[[i]] <- tmp_slice

}

# save the slice of data
save(tmp_slice_list,ts_list,lat_grid,lon_grid,lon,lat, file = "output/CALCOFI_OI_SST_Slice.Rdata")

# Calculate Coeff. Var.

load("output/CALCOFI_OI_SST_Slice.Rdata")

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
    
   for(k in 1:length(tmp_slice_list)){
     
     vals[[k]] <- tmp_slice_list[[k]][i,j,]
     
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

save(coeff_table = coeff_table, mean_table = mean_table, file = "output/CALCOFI_temp_tables.Rdata")









