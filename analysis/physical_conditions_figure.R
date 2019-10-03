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
library(scales)

map <- map_data("world")    

load("output/uv_velocity_table.Rdata")
load("output/CALCOFI_bathymetry_table.Rdata")
load("output/CALCOFI_temp_tables.Rdata")
load("output/cyano_16s_map.Rdata")


stations <-  ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "n_samps"), color = "black", size =4, stroke = 0.1, shape = 21) +
  ggtitle("CalCOFI NCOG Stations") + scale_fill_gradient(name = "# of Samples", low = "white", high = "red") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
        title = element_text(hjust = 0.5), axis.line = element_blank())

reduced_uv_table <- uv_table[seq(1,1633,6),]

vel_bath <- ggplot() + 
  geom_raster(data = elevation_table, aes(x = lon, y = lat, fill = value), interpolate = FALSE) + 
  scale_fill_gradient(low = "darkblue", high = "cyan", name = "Depth (m)",
                      limits = c(-5000,0), oob = scales::squish) + 
  metR::geom_vector(data = reduced_uv_table, aes(x = lon, y = lat, dx = Mean_U, dy = Mean_V), 
                    arrow.angle = 15, arrow.type = "open", arrow.length = unit(0.5, "inches"), 
                    pivot = 0,preserve.dir = TRUE, direction = "ccw",
                    min.mag = 0, show.legend = NA, color = "black") +
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  theme(legend.position = "right",
        legend.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linetype = "solid", size = 1),
        title = element_text(hjust = 0.5), axis.line = element_blank()) +
  scale_mag(max = 0.1, name = "Speed (m/s)", max_size = 0.75) +
  labs(x = "Longitude", y = "Latitude", color = "Depth (m)") + ggtitle("Mean Geostrophic Current Velocity")

sst <- ggplot() + 
  geom_tile(data = coeff_table, aes(x = lon, y = lat, fill = coeff_var), width =0.26, height = 0.26) +
  scale_fill_gradient2(name = "Coeff. Var SST", low = "darkblue", mid = "white", high = "darkred", limits = c(0.09,0.12), oob = squish, midpoint = 0.1066851) +geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Coeff. Var. SST") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
        title = element_text(hjust = 0.5), axis.line = element_blank())

sst_mean <- ggplot() + 
  geom_tile(data = mean_table, aes(x = lon, y = lat, fill = Mean), width =0.26, height = 0.26) +
  scale_fill_gradient2(name = "SST Mean (°C)", low = "darkblue", mid = "white", high = "darkred", limits = c(15,18), oob = squish, midpoint = 16.5) +geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("SST Mean (°C)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
        title = element_text(hjust = 0.5), axis.line = element_blank())

theme_set(theme_cowplot(font_size=12))
pdf("figures/summary_physical_figure.pdf", width = 12, height = 12)
plot_grid(stations, vel_bath,
          sst_mean, sst, ncol = 2, nrow = 2, align = "hv")
dev.off()


# bath line plot
vel_bath <- ggplot() + 
  geom_contour(data = elevation_table, aes(x = lon, y = lat, z = value,
                                           color = factor(..level.. == -2000,
                                                          levels = c(T),
                                                          labels = c("2000"))),
               breaks = -2000, size = 1) + 
  scale_color_manual(name = "Depth (m)", values = c("red")) +
  metR::geom_vector(data = reduced_uv_table, aes(x = lon, y = lat, dx = Mean_U, dy = Mean_V), 
                    arrow.angle = 15, arrow.type = "open", arrow.length = unit(0.5, "inches"), 
                    pivot = 0,preserve.dir = TRUE, direction = "ccw",
                    min.mag = 0, show.legend = NA, color = "black") +
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
  coord_fixed(xlim = c(-126.5, -116.5),ylim= c(28,37), 1.3) +
  theme(legend.position = "right",
        legend.key.height = unit(1.4, "cm"), 
        legend.background = element_blank(),
        axis.text = element_text(size = 12, colour = 1),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black", linetype = "solid", size = 1),
        title = element_text(hjust = 0.5), axis.line = element_blank()) +
  scale_mag(max = 0.1, name = "Speed (m/s)", max_size = 0.75) +
  labs(x = "Longitude", y = "Latitude", color = "Depth (m)") + ggtitle("Mean Geostrophic Current Velocity")

