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

load("output/plastid_16s_glm.Rdata")

som_plots$cluster[which(som_plots$cluster == "som_1")] <- "Offshore"
som_plots$cluster[which(som_plots$cluster == "som_2")] <- "Nearshore"

ggplot(som_plots, aes(x = NC_mean, y = freq, color = cluster)) + geom_point() +
  stat_smooth(method="glm", method.args = list(family = "binomial")) + 
  theme(legend.title = element_blank(), legend.background = element_rect(fill = NA, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black")) +
  xlab("Nitracline Depth (m)") + ylab("Frequency") + scale_color_manual(values = c("red", "blue"))


map <- map_data("world") 

anom_temp <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "temp_anom_mean"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, aes(min = min(temp_anom_mean, na.rm = TRUE, max = max(temp_anom_mean, na.rm = TRUE)), limits = c(min,max))) +
  ggtitle(paste0("Mean Temperature Anomalies")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))

mean_temp <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes(x = long, y = lat, fill = temp_season_mean), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 15) +
  ggtitle(paste0("Mean Climatic Temperature")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))

diff_temp <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes(x = long, y = lat, fill = abs(temp_anom_mean)/temp_season_mean), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.075) +
  ggtitle(paste0("N'/N_bar Temperature")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))


anom_sal <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "sal_anom_mean"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, aes(min = min(sal_anom_mean, na.rm = TRUE, max = max(sal_anom_mean, na.rm = TRUE)), limits = c(min,max))) +
  ggtitle(paste0("Mean Salinity Anomalies")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))

mean_sal <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes(x = long, y = lat, fill = sal_season_mean), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 33.3) +
  ggtitle(paste0("Mean Climatic Salinity")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))

diff_sal <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes(x = long, y = lat, fill = abs(sal_anom_mean)/sal_season_mean), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.005) +
  ggtitle(paste0("N'/N_bar Salinity")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))

anom_no3 <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "NO3_anom_mean"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, aes(min = min(NO3_anom_mean, na.rm = TRUE, max = max(NO3_anom_mean, na.rm = TRUE)), limits = c(min,max))) +
  ggtitle(paste0("Mean NO3 Anomalies")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))

mean_no3 <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes(x = long, y = lat, fill = NO3_season_mean), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 3) +
  ggtitle(paste0("Mean Climatic NO3")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))

diff_no3 <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes(x = long, y = lat, fill = abs(NO3_anom_mean)/NO3_season_mean), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.55) +
  ggtitle(paste0("N'/N_bar NO3")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))

anom_PO4 <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "PO4_anom_mean"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, aes(min = min(PO4_anom_mean, na.rm = TRUE, max = max(PO4_anom_mean, na.rm = TRUE)), limits = c(min,max))) +
  ggtitle(paste0("Mean PO4 Anomalies")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))

mean_po4 <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes(x = long, y = lat, fill = PO4_season_mean), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5) +
  ggtitle(paste0("Mean Climatic PO4")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))

diff_po4 <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes(x = long, y = lat, fill = abs(PO4_anom_mean)/PO4_season_mean), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.3) +
  ggtitle(paste0("N'/N_bar PO4")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))

anom_SiO3 <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "SiO3_anom_mean"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, aes(min = min(SiO3_anom_mean, na.rm = TRUE, max = max(SiO3_anom_mean, na.rm = TRUE)), limits = c(min,max))) +
  ggtitle(paste0("Mean SiO3 Anomalies")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))

mean_sio3 <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes(x = long, y = lat, fill = SiO3_season_mean), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 4) +
  ggtitle(paste0("Mean Climatic SiO3")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))

diff_sio3 <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes(x = long, y = lat, fill = abs(SiO3_anom_mean)/SiO3_season_mean), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.3) +
  ggtitle(paste0("N'/N_bar SiO3")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5))

coeff_temp <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "temp_anom_coeff"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-2,2), oob = scales::squish) +
  ggtitle(paste0("Coeff. Temperature Anomalies")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        title = element_text(hjust = 0.5))

coeff_sal <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "sal_anom_coeff"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-4,4), oob = scales::squish) +
  ggtitle(paste0("Coeff. Salinity Anomalies")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        title = element_text(hjust = 0.5))

coeff_no3 <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "NO3_anom_coeff"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-4,4), oob = scales::squish) +
  ggtitle(paste0("Coeff. NO3 Anomalies")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        title = element_text(hjust = 0.5))

coeff_PO4 <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "PO4_anom_coeff"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-4,4), oob = scales::squish) +
  ggtitle(paste0("Coeff. PO4 Anomalies")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        title = element_text(hjust = 0.5))

coeff_SiO3 <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "SiO3_anom_coeff"), color = "black", size =5, stroke = 0.1, shape = 21) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,limits = c(-4,4), oob = scales::squish) +
  ggtitle(paste0("Coeff. SiO3 Anomalies")) +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        title = element_text(hjust = 0.5))


pdf(file = "figures/anom_spatial.pdf", width = 10, height = 20)
plot_grid(anom_temp, mean_temp, diff_temp,
          anom_sal, mean_sal, diff_sal,
          anom_no3, mean_no3, diff_no3,
          anom_PO4, mean_po4, diff_po4,
          anom_SiO3, mean_sio3, diff_sio3,
          ncol = 3, nrow = 5)
dev.off()

load("data/Station_anomaly.Rdata")

stations <- vector()
for (i in 1:length(Station_Anoms)) {stations[i] <- as.character(Station_Anoms[[i]]$Sta_ID)}



off <- Station_Anoms[[10]]$anom_mat
on <- Station_Anoms[[14]]$anom_mat


plot(on[,1], type = "l", ylim = c(12,20))
points(off[,1], type = "l", color = "red")


