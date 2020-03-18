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
library(geosphere)
library(numDeriv)

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

ow_array <- array(NA, dim = c(47,47,dim(u_slice)[3]))
eddy_array <- array(NA, dim = c(47,47,dim(u_slice)[3]))

for (i in 1:dim(u_slice)[3]) {
  
  u_mat <- u_slice[,,i]
  v_mat <- v_slice[,,i]
  
  du_dx_mat <- matrix(NA,nrow = nrow(mat), ncol = (ncol(mat)-1))
  dv_dx_mat <- matrix(NA,nrow = nrow(mat), ncol = (ncol(mat)-1))
  dv_dy_mat <- matrix(NA,nrow = (nrow(mat)-1), ncol = ncol(mat))
  du_dy_mat <- matrix(NA,nrow = (nrow(mat)-1), ncol = ncol(mat))
  
  
  # normal gradients
  for (j in 1:length(lat_grid)) {
      du <- diff(u_mat[j,])
      dx_mat <- matrix(NA,(length(lon_grid)-1),4)
      dx_mat[,1] <- lon[lon_grid[1:(length(lon_grid)-1)]]
      dx_mat[,3] <- lon[lon_grid[2:length(lon_grid)]]
      dx_mat[,2] <- rep(lat[lat_grid[j]])
      dx_mat[,4] <- rep(lat[lat_grid[j]])
      dx <- distHaversine(dx_mat[,1:2], dx_mat[,3:4])/1000
      du_dx_mat[j,] <- du/dx
  }
  
  for (k in 1:length(lon_grid)) {
    dv <- diff(v_mat[,k])
    dy_mat <- matrix(NA,(length(lat_grid)-1),4)
    dy_mat[,2] <- lat[lat_grid[1:(length(lat_grid)-1)]]
    dy_mat[,4] <- lat[lat_grid[2:length(lat_grid)]]
    dy_mat[,1] <- rep(lon[lon_grid[j]])
    dy_mat[,3] <- rep(lon[lon_grid[j]])
    dy <- distHaversine(dx_mat[,1:2], dx_mat[,3:4])/1000
    dv_dy_mat[,k] <- dv/dy
  }
  
  # shear
  
  for (j in 1:length(lat_grid)) {
    dv <- diff(v_mat[j,])
    dx_mat <- matrix(NA,(length(lon_grid)-1),4)
    dx_mat[,1] <- lon[lon_grid[1:(length(lon_grid)-1)]]
    dx_mat[,3] <- lon[lon_grid[2:length(lon_grid)]]
    dx_mat[,2] <- rep(lat[lat_grid[j]])
    dx_mat[,4] <- rep(lat[lat_grid[j]])
    dx <- distHaversine(dx_mat[,1:2], dx_mat[,3:4])/1000
    dv_dx_mat[j,] <- dv/dx
  }
  
  for (k in 1:length(lon_grid)) {
    du <- diff(u_mat[,k])
    dy_mat <- matrix(NA,(length(lat_grid)-1),4)
    dy_mat[,2] <- lat[lat_grid[1:(length(lat_grid)-1)]]
    dy_mat[,4] <- lat[lat_grid[2:length(lat_grid)]]
    dy_mat[,1] <- rep(lon[lon_grid[j]])
    dy_mat[,3] <- rep(lon[lon_grid[j]])
    dy <- distHaversine(dx_mat[,1:2], dx_mat[,3:4])/1000
    du_dy_mat[,k] <- du/dy
  }
  
  ow_mat <- (du_dx_mat[-nrow(du_dx_mat),] - dv_dy_mat[,-ncol(dv_dy_mat)])^2 +
    (dv_dx_mat[-nrow(dv_dx_mat),] + du_dy_mat[,-ncol(du_dy_mat)])^2 -
    (dv_dx_mat[-nrow(dv_dx_mat),] - du_dy_mat[,-ncol(du_dy_mat)])^2
  
  eddy_mat <- ow_mat
  eddy_mat[eddy_mat >= (-2*10^-10)] = 0

  ow_array[,,i] <- ow_mat
  eddy_array[,,i] <- eddy_mat
  
  }

save(ow_array, eddy_array, lon, lat, lon_grid, lat_grid, file = "output/ow_eddy_slice.Rdata")

load("output/ow_eddy_slice.Rdata")

ow_mean <- matrix(NA,nrow(ow_array),ncol(ow_array))
ow_coeff <- matrix(NA,nrow(ow_array),ncol(ow_array))
eddy_sum <- matrix(NA,nrow(ow_array),ncol(ow_array))

rownames(ow_mean) <- lon[lon_grid[1:47]]-360
rownames(ow_coeff) <- lon[lon_grid[1:47]]-360
rownames(eddy_sum) <- lon[lon_grid[1:47]]-360

colnames(ow_mean) <- lat[lat_grid[1:47]]
colnames(ow_coeff) <- lat[lat_grid[1:47]]
colnames(eddy_sum) <- lat[lat_grid[1:47]]

for (i in 1:nrow(ow_array)) {
  for (j in 1:ncol(ow_array)) {
    
    ow_vect <- vector()
    eddy_vect <- vector()
    
    for (k in 1:dim(ow_array)[3]) {
      
      ow_vect[k] <- ow_array[i,j,k]
      eddy_vect[k] <- eddy_array[i,j,k]
      
    }
    
    ow_mean[i,j] <- mean(ow_vect)
    ow_coeff[i,j] <- sd(ow_vect)/mean(ow_vect)
    
    eddy_vect[eddy_vect != 0] <- 1
    
    eddy_sum[i,j] <- sum(eddy_vect)/dim(ow_array)[3]
   
    
  }
}


coeff_table <- melt(ow_coeff)
mean_table <- melt(ow_mean)
sum_table <- melt(eddy_sum)

colnames(coeff_table) <- c("lon", "lat", "coeff_var")
colnames(mean_table) <- c("lon", "lat", "Mean")
colnames(sum_table) <- c("lon", "lat", "Sum")

map <- map_data("world")    

coeff_table <- coeff_table[which(!is.na(coeff_table$coeff_var)),]
mean_table <- mean_table[which(!is.na(mean_table$Mean)),]
sum_table <- sum_table[which(!is.na(sum_table$Sum)),]

save(coeff_table = coeff_table, mean_table = mean_table, sum_table = sum_table,
     file = "output/CALCOFI_okubo_weiss_tables.Rdata")



mean_ow <- ggplot() + 
  geom_tile(data = mean_table, aes(x = lon, y = lat, fill = Mean), width =0.3, height = 0.3) +
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) + 
  theme(panel.background = element_blank()) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0.000017, limits = c(), oob = scales::squish) +
  labs(fill = "Mean\nOkubo-Weiss\nParameter") + xlab("Longitude") + ylab("Latitude")

coeff_ow <- ggplot() + 
  geom_tile(data = coeff_table, aes(x = lon, y = lat, fill = coeff_var), width =0.3, height = 0.3) +
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) + 
  theme(panel.background = element_blank()) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 3.2,
                       limits = c(2,6), oob = scales::squish) +
  labs(fill = "Coeff. Var.\nOkubo-Weiss\n Parameter") + xlab("Longitude") + ylab("Latitude")

eddy_freq <- ggplot() + 
  geom_tile(data = sum_table, aes(x = lon, y = lat, fill = Sum), width =0.3, height = 0.3) +
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) + 
  theme(panel.background = element_blank()) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.35,
                       limits = c(0.25,0.45), oob = scales::squish) +
  labs(fill = "Eddy\nFrequency") + xlab("Longitude") + ylab("Latitude")

pdf("figures/okubo_weiss_figure.pdf", height = 4, width = 12)
plot_grid(mean_ow,coeff_ow,eddy_freq, ncol = 3)
dev.off()
