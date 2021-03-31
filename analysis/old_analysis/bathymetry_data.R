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

# Gridded Bathymetry
# https://www.gebco.net/data_and_products/gridded_bathymetry_data/

# set path and filename

ncpath <- "data/"
ncname <- "GEBCO_2019_california"
ncfname <- paste(ncpath, ncname, ".nc", sep="")
dname <- "elevation"  # note: sea surface height relative to geoid

ncin <- nc_open(ncfname)

# get lat and long

lon <- ncin$var$elevation$dim[[1]]$vals
lat <- ncin$var$elevation$dim[[2]]$vals

# Get temperature

elevation_array <- ncvar_get(ncin,dname)
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")

# getting dims for array

min_lat = 27.5
max_lat = 39.5
min_lon = (360-127.5)-360
max_lon = (360-115.5)-360

lat_grid <- which(lat > min_lat & lat < max_lat)
lon_grid <-which(lon > min_lon & lon < max_lon)

elevation_slice <- elevation_array[lon_grid,lat_grid]

# save the slice of data
save(elevation_slice,file = "output/CALCOFI_GEBCO_bathymetry_Slice.Rdata")

colnames(elevation_slice) <- lat[lat_grid]
rownames(elevation_slice) <- lon[lon_grid]
elevation_table <- melt(elevation_slice)

colnames(elevation_table) <- c("lon", "lat", "value")


elevation_table <- elevation_table[which(!is.na(elevation_table$value)),]

save(elevation_table = elevation_table, file = "output/CALCOFI_bathymetry_table.Rdata")

elevation_table$value[elevation_table$value >= 0] <- NA

# Figures

map <- map_data("world")    

load("output/uv_velocity_table.Rdata")
load("output/CALCOFI_bathymetry_table.Rdata")

ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  geom_raster(data = elevation_table, aes(x = lon, y = lat, fill = value), interpolate = FALSE) + 
  scale_fill_viridis_c(name = "Depth (m)")


currents <- ggplot() + 
  geom_raster(data = elevation_table, aes(x = lon, y = lat, fill = value), interpolate = FALSE) + 
  scale_fill_viridis_c(name = "Depth (m)", option = "E") + 
  metR::geom_vector(data = uv_table,aes(x = lon, y = lat, dx = Mean_U, dy = Mean_V), 
                    arrow.angle = 15, arrow.type = "open", arrow.length = unit(0.5, "inches"), 
                    pivot = 0,preserve.dir = TRUE, direction = "ccw",
                    min.mag = 0, show.legend = NA, color = "white") +
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
  coord_fixed(xlim = c(-126.5, -116.5),ylim= c(28,37), 1.3) +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.height = unit(1.4, "cm"), 
        legend.background = element_blank(),
        axis.text = element_text(size = 12, colour = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_mag(max = 0.1, name = "Speed (m/s)", max_size = 0.75) +
  labs(x = "Longitude", y = "Latitude") + ggtitle("")

pdf("figures/currents_bathymetry.pdf", width = 8, height = 8)
currents
dev.off()

reduced_uv_table <- uv_table[seq(1,1633,3),]

black_vel <- ggplot() + 
  geom_raster(data = elevation_table, aes(x = lon, y = lat, fill = value), interpolate = FALSE) + 
  scale_fill_gradient(low = "darkblue", high = "cyan", name = "Depth (m)",
                      limits = c(-5000,0), oob = scales::squish) + 
  metR::geom_vector(data = reduced_uv_table,aes(x = lon, y = lat, dx = Mean_U, dy = Mean_V), 
                    arrow.angle = 15, arrow.type = "open", arrow.length = unit(0.5, "inches"), 
                    pivot = 0,preserve.dir = TRUE, direction = "ccw",
                    min.mag = 0, show.legend = NA, color = "black") +
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
  coord_fixed(xlim = c(-126.5, -116.5),ylim= c(28,37), 1.3) +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.height = unit(1.4, "cm"), 
        legend.background = element_blank(),
        axis.text = element_text(size = 12, colour = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_mag(max = 0.1, name = "Speed (m/s)", max_size = 0.75) +
  labs(x = "Longitude", y = "Latitude") + ggtitle("Mean Geostrophic Current Velocity (1993-Present)")

white_vel <- ggplot() + 
  geom_raster(data = elevation_table, aes(x = lon, y = lat, fill = value), interpolate = FALSE) + 
  scale_fill_gradient(low = "darkblue", high = "cyan", name = "Depth (m)",
                      limits = c(-5000,0), oob = scales::squish) + 
  metR::geom_vector(data = reduced_uv_table,aes(x = lon, y = lat, dx = Mean_U, dy = Mean_V), 
                    arrow.angle = 15, arrow.type = "open", arrow.length = unit(0.5, "inches"), 
                    pivot = 0,preserve.dir = TRUE, direction = "ccw",
                    min.mag = 0, show.legend = NA, color = "white") +
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
  coord_fixed(xlim = c(-126.5, -116.5),ylim= c(28,37), 1.3) +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.height = unit(1.4, "cm"), 
        legend.background = element_blank(),
        axis.text = element_text(size = 12, colour = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_mag(max = 0.1, name = "Speed (m/s)", max_size = 0.75, guide = "none") +
  labs(x = "Longitude", y = "Latitude") + ggtitle("Mean Geostrophic Current Velocity (1993-Present)")


pdf("figures/currents_bathymetry_white.pdf", width = 8, height = 8)
white_vel
dev.off()

pdf("figures/currents_bathymetry_black.pdf", width = 8, height = 8)
black_vel
dev.off()
