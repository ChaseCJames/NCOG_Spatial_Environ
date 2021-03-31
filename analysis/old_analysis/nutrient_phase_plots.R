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
library(reshape2)
library(chron)


###### CalCOFI Nutrient Data #########

calcofi_dat <- read.csv("data/194903-201806_Bottle.csv", header = TRUE, stringsAsFactors = FALSE)
physical_dat <- read.csv("data/NCOG_sample_metadata.csv", header = TRUE, stringsAsFactors = FALSE)

calcofi_dat$year <- as.numeric(paste0(substr(calcofi_dat$Depth_ID,1,2),substr(calcofi_dat$Depth_ID,4,5)))

calcofi_dat <- calcofi_dat[which(calcofi_dat$year > 2013),]

# surface samples only

calcofi_dat_surf <- calcofi_dat[which(calcofi_dat$Depthm < 20),]
calcofi_dat_50 <- calcofi_dat[which(calcofi_dat$Depthm == 50),]

# deep samples only

# physical_dat <- physical_dat[which(physical_dat$CC_Depth > 20),]

# add lat and long

calcofi_dat_50$Lat_Dec <- physical_dat$Lat_Dec[match(calcofi_dat_50$Sta_ID, physical_dat$Sta_ID)]
calcofi_dat_50$Lon_Dec <- physical_dat$Lon_Dec[match(calcofi_dat_50$Sta_ID, physical_dat$Sta_ID)]

calcofi_dat_surf$Lat_Dec <- physical_dat$Lat_Dec[match(calcofi_dat_surf$Sta_ID, physical_dat$Sta_ID)]
calcofi_dat_surf$Lon_Dec <- physical_dat$Lon_Dec[match(calcofi_dat_surf$Sta_ID, physical_dat$Sta_ID)]

calcofi_dat_50 <- calcofi_dat_50[-which(is.na(calcofi_dat_50$Lat_Dec)),]
calcofi_dat_surf <- calcofi_dat_surf[-which(is.na(calcofi_dat_surf$Lat_Dec)),]
# check variables

data_50 <- calcofi_dat_50[colnames(calcofi_dat_50)[c(3,75,76,77,32,26)]]
data_surf <- calcofi_dat_surf[colnames(calcofi_dat_surf)[c(3,75,76,77,32,26)]]

data_50 <- na.omit(data_50)
data_surf <- na.omit(data_surf)

data <- data_surf

warm_dat <- data[which(data$year<2017),]
cool_dat <- data[which(data$year>2016),]

warm_maps <- warm_dat %>% 
  group_by(Sta_ID) %>%
  summarise(n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE),
            PO4_mean = mean(PO4uM, na.rm = TRUE), NO3_mean = mean(NO3uM, na.rm = TRUE),
            PO4_coeff =  sd(PO4uM, na.rm = TRUE)/mean(PO4uM, na.rm = TRUE),
            NO3_coeff = sd(NO3uM, na.rm = TRUE)/mean(NO3uM, na.rm = TRUE)
  )

cool_maps <- cool_dat %>% 
  group_by(Sta_ID) %>%
  summarise(n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE),
            PO4_mean = mean(PO4uM, na.rm = TRUE), NO3_mean = mean(NO3uM, na.rm = TRUE),
            PO4_coeff =  sd(PO4uM, na.rm = TRUE)/mean(PO4uM, na.rm = TRUE),
            NO3_coeff = sd(NO3uM, na.rm = TRUE)/mean(NO3uM, na.rm = TRUE)
  )

line <- as.numeric(substr(warm_maps$Sta_ID,2,3))
n_stations <- which(line < 76)
if(length(n_stations) > 0){warm_maps <- warm_maps[-n_stations,]}

line <- as.numeric(substr(cool_maps$Sta_ID,2,3))
n_stations <- which(line < 76)
if(length(n_stations) > 0){cool_maps <- cool_maps[-n_stations,]}

map <- map_data("world") 

warm_list <- list()
cool_list <- list()

min_max_matrix <- as.data.frame(matrix(data = NA, nrow = 9, ncol = 2))

for (i in 5:6) {
  
  min_max_matrix[i-4,1] <- min(c(min(warm_maps[,i],na.rm = TRUE), min(cool_maps[,i], na.rm=TRUE)), na.rm = TRUE)
  min_max_matrix[i-4,2] <- max(c(min(warm_maps[,i], na.rm = TRUE), max(cool_maps[,i], na.rm = TRUE)), na.rm = TRUE)
}

colnames(min_max_matrix) <- c("min","max")

title_names <- c("Mean PO4 ug/L", "Mean NO3 ug/L", "PO4 Coeff. Var.", "NO3 Coeff. Var.")

limit_diffs <- c(0.3,4,5,0,0.2,1,0.1,0,3)

difference_high_limits <- c(0.1,1,1,2,0.25,1,0.5,2,0)
difference_low_limits <- c(-0.05,-2,-1,-2,-0.25,-2,-0.5,-2,-3)

diff_list <- list()

for (i in 5:6) {
  
  p <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = warm_maps, aes_string(x = "long", y = "lat", fill = colnames(warm_maps)[i]), color = "black", size =3, stroke = 0.1, shape = 21) +
    scale_fill_gradient(limits = c(min_max_matrix$min[i-4], min_max_matrix$max[i-4]-limit_diffs[i-4]),
                        na.value = "grey",
                        oob = scales::squish,
                        low = "white", high = "red") +
    ggtitle(paste0(title_names[i-4]," 2014-2016")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          plot.title = element_text(hjust = 0.5))
  
  print(p)
  
  warm_list[[i-4]] <- p
  
  q <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = cool_maps, aes_string(x = "long", y = "lat", fill = colnames(cool_maps)[i]), color = "black", size =3, stroke = 0.1, shape = 21) +
    scale_fill_gradient(limits = c(min_max_matrix$min[i-4], min_max_matrix$max[i-4]-limit_diffs[i-4]),
                        na.value = "grey",
                        oob = scales::squish,
                        low = "white", high = "red") +
    ggtitle(paste0(title_names[i-4]," 2017-2018")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          plot.title = element_text(hjust = 0.5))
  
  print(q)
  
  cool_list[[i-4]] <- q
  
  diff_plots <- as.data.frame(matrix(NA,nrow(warm_maps),3))
  diff_plots[,1] <- warm_maps$long
  diff_plots[,2] <- warm_maps$lat
  diff_plots[,3] <- warm_maps[,i]-cool_maps[match(warm_maps$Sta_ID, cool_maps$Sta_ID),i]
  
  colnames(diff_plots) <- c("long","lat","Difference")
  
  r <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = diff_plots, aes_string(x = "long", y = "lat", fill = "Difference"), color = "black", size =3, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(na.value = "grey", low = "blue", mid = "white", high = "red", 
                         limits = c(difference_low_limits[i-4], difference_high_limits[i-4]),
                         oob = scales::squish) +
    ggtitle(paste0("Difference in ",title_names[i-4],"\n (2014-2016) - (2017-2018)")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          plot.title = element_text(hjust = 0.5))
  
  print(r)
  
  diff_list[[i-4]] <- r 
  
}

pdf(file = "figures/calcofi_hydrographic_data_means.pdf", width = 12, height = 6)
plot_grid(warm_list[[2]], cool_list[[2]], diff_list[[2]],
          warm_list[[1]], cool_list[[1]], diff_list[[1]],
          ncol = 3, labels = c("A","B","C",
                               "D", "E", "F"))
dev.off()

pdf(file = "figures/calcofi_hydrographic_data_coeff.pdf", width = 12, height = 11)
plot_grid(warm_list[[6]], cool_list[[6]], diff_list[[6]],
          warm_list[[5]], cool_list[[5]], diff_list[[5]],
          warm_list[[7]], cool_list[[7]], diff_list[[7]],
          ncol = 3, labels = c("A","B","C",
                               "D", "E", "F",
                               "G", "H", "I"))
dev.off()
