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
library(directlabels)

# set path and filename

ncpath <- "data/"
ncname_90 <- "total_z_90" 
ncname_80 <- "total_z_80" 
ncfname_90 <- paste(ncpath, ncname_90, ".nc", sep="")
ncfname_80 <- paste(ncpath, ncname_80, ".nc", sep="")
dname <- "potential_density"  # note: sea surface height anomaly

ncin_90 <- nc_open(ncfname_90)
ncin_80 <- nc_open(ncfname_80)

# get lat and long

dist_90 <- ncin_90$var$longitude$dim[[1]]$vals
dist_80 <- ncin_80$var$longitude$dim[[1]]$vals

depth_90 <- ncin_90$var$potential_density$dim[[1]]$vals
depth_80 <- ncin_80$var$potential_density$dim[[1]]$vals
# get time

time_90 <- ncvar_get(ncin_90,"time")
tunits_90 <- ncatt_get(ncin_90,"time","units")

time_80 <- ncvar_get(ncin_80,"time")
tunits_80 <- ncatt_get(ncin_80,"time","units")

# Get sla

pd_90_array <- ncvar_get(ncin_90,dname)
pd_80_array <- ncvar_get(ncin_80,dname)
dunits_90 <- ncatt_get(ncin_90,dname,"units")
dunits_90 <- ncatt_get(ncin_80,dname,"units")
fillvalue_90 <- ncatt_get(ncin_90,dname,"_FillValue")
fillvalue_80 <- ncatt_get(ncin_80,dname,"_FillValue")

# Convert Time

tustr <- strsplit(tunits_90$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth <- as.integer(unlist(tdstr)[2])
tday <- as.integer(strsplit(unlist(tdstr)[3],"T")[[1]][1])
tyear <- as.integer(unlist(tdstr)[1])
ts_90 <- as.Date(as.POSIXct(time_90,origin=paste0(tyear,"-",tmonth,"-",tday)))
ts_80 <- as.Date(as.POSIXct(time_80,origin=paste0(tyear,"-",tmonth,"-",tday)))

pd_90_array[pd_90_array==fillvalue_90$value] <- NA
pd_80_array[pd_80_array==fillvalue_80$value] <- NA

save(pd_80_array, pd_90_array, ts_80, ts_90, dist_80, dist_90, depth_80, depth_90,
     file = "spray_data_80_90.Rdata")

