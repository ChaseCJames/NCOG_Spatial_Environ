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


glider_dat_analysis <- function(in_data = "output/euks_auto_18sv9_full_data.Rdata",
                                in_glider = "output/spray_data_80_90.Rdata"){


load(in_data)
load(in_glider)  
  
# grab only glider lines

full_dat$line <- substr(full_dat$Sta_ID,2,3)

glider_dat <- full_dat %>%
  filter(line == 80 | line == 90)

full_som_cruise <- full_dat %>% 
  group_by(Cruise) %>%
  summarise(som_1 = sum(som_id == 1, na.rm = TRUE)/n(),
            som_2 = sum(som_id == 2, na.rm = TRUE)/n(),
            n_samps = n(), 
            evenness = mean(evenness, na.rm = TRUE),
            shannon = mean(shannon, na.rm = TRUE),
            richness = mean(richness, na.rm = TRUE), 
            Date = mean(Date, na.rm = TRUE))

glider_som_cruise <- glider_dat %>% 
  group_by(Cruise) %>%
  summarise(som_1 = sum(som_id == 1, na.rm = TRUE)/n(),
            som_2 = sum(som_id == 2, na.rm = TRUE)/n(),
            n_samps = n(), 
            evenness = mean(evenness, na.rm = TRUE),
            shannon = mean(shannon, na.rm = TRUE),
            richness = mean(richness, na.rm = TRUE), 
            Date = mean(Date, na.rm = TRUE))


total_80_matrix <- matrix(NA,length(depth_80),length(dist_80))
total_90_matrix <- matrix(NA,length(depth_90),length(dist_90))

for (j in 1:length(dist_80)) {
  for(k in 1:length(depth_80)){
    
    total_80_matrix[k,j] <- mean(pd_80_array[k,,j], na.rm = TRUE)
    
  }
}

for (j in 1:length(dist_90)) {
  for(k in 1:length(depth_90)){
    
    total_90_matrix[k,j] <- mean(pd_90_array[k,,j], na.rm = TRUE)
    
  }
} 

  
mean_80_array <- array(NA,dim = c(length(depth_80),length(dist_80),nrow(full_som_cruise)))
mean_90_array <- array(NA,dim = c(length(depth_90),length(dist_90),nrow(full_som_cruise)))

anomaly_80_array <- array(NA,dim = c(length(depth_80),length(dist_80),nrow(full_som_cruise)))
anomaly_90_array <- array(NA,dim = c(length(depth_90),length(dist_90),nrow(full_som_cruise)))

for (i in 1:nrow(full_som_cruise)) {
  
  first_date <- seq(full_som_cruise$Date[i], length = 2, by = "-1 months")[2]
  
  dates <- seq.Date(as.Date(first_date),as.Date(full_som_cruise$Date[i]), by = "days")
  
  line_80_vals <- which(!is.na(match(ts_80, dates)))
  line_90_vals <- which(!is.na(match(ts_90, dates)))
  
  for (j in 1:length(dist_80)) {
    for(k in 1:length(depth_80)){
      
      mean_80_array[k,j,i] <- mean(pd_80_array[k,line_80_vals,j])
      anomaly_80_array[k,j,i] <- mean(pd_80_array[k,line_80_vals,j])-total_80_matrix[k,j]
    }
  }
  
  for (j in 1:length(dist_90)) {
    for(k in 1:length(depth_90)){
      
      mean_90_array[k,j,i] <- mean(pd_90_array[k,line_90_vals,j])
      anomaly_90_array[k,j,i] <- mean(pd_90_array[k,line_90_vals,j])-total_90_matrix[k,j]
    }
  } 
  
  mat_80 <- mean_80_array[,,i]
  mat_90 <- mean_90_array[,,i]
  
  colnames(mat_80) <- dist_80
  rownames(mat_80) <- depth_80
  
  colnames(mat_90) <- dist_90
  rownames(mat_90) <- depth_90
  
  plot_df_80 <- melt(mat_80, varnames = c("Depth", "Distance"))
  plot_df_90 <- melt(mat_90, varnames = c("Depth", "Distance"))
  
  subset_dat <- glider_dat[which(!is.na(match(glider_dat$Cruise, full_som_cruise$Cruise[i]))),]
  
  subset_dat_90 <- subset_dat %>%
    filter(line == 90)
  subset_dat_80 <- subset_dat %>%
    filter(line == 80)
  
  dens_plot_80 <- ggplot(plot_df_80, aes(y = Depth, x = Distance, fill = value)) + 
    geom_tile() + 
    geom_contour(data = plot_df_80, aes(z = value,
                                     color = factor(..level.. == 26,
                                                    levels = c(T),
                                                    labels = c("26"))),
                 breaks = 24:28, size = 1) +
    scale_color_manual(values = "red", na.value = NA, na.translate = FALSE, guide = FALSE) +
    scale_fill_viridis_c(direction = -1, na.value = "grey") +
    geom_point(data = subset_dat_80, aes(x = dist_to_coast, y = Depthm, shape = as.factor(som_id)),
               stat = "identity", inherit.aes = FALSE) + 
    scale_shape_manual(values = c(17,19), guide = FALSE) +
    scale_x_reverse() + scale_y_reverse() +
    xlab("Distance to Coast (km)") + ylab("Depth (m)") +
    labs(fill = "Potential\nDensity", shape = "SOM ID") +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust= 0.5)) +
    ggtitle("Line 80")
  
  leg <- get_legend(dens_plot_80)
  
  dens_plot_90 <- ggplot(plot_df_90, aes(y = Depth, x = Distance, fill = value)) + 
    geom_tile() + 
    geom_contour(data = plot_df_90, aes(z = value,
                                        color = factor(..level.. == 26,
                                                       levels = c(T),
                                                       labels = c("26"))),
                 breaks = 24:28, size = 1) +
    scale_color_manual(values = "red", na.value = NA, na.translate = FALSE, guide = FALSE) +
    scale_fill_viridis_c(direction = -1, na.value = "grey") +
    geom_point(data = subset_dat_90, aes(x = dist_to_coast, y = Depthm, shape = as.factor(som_id)),
               stat = "identity", inherit.aes = FALSE) + 
    scale_shape_manual(values = c(17,19), guide = FALSE) +
    scale_x_reverse() + scale_y_reverse() +
    xlab("Distance to Coast (km)") + ylab("Depth (m)") +
    labs(fill = "Potential\nDensity", shape = "SOM ID") +
    theme(panel.background = element_blank(),
          plot.title = element_text(hjust= 0.5)) +
    ggtitle("Line 90")
  
  dens_plot_cont_80 <- direct.label(dens_plot_80, "top.pieces")
  dens_plot_cont_90 <- direct.label(dens_plot_90, "top.pieces")
  
  title <- ggdraw() + draw_label(full_som_cruise$Cruise[i], fontface='bold')
  full_plot <- plot_grid(title,
  plot_grid(dens_plot_cont_80, dens_plot_cont_90, leg, rel_widths = c(1,1,0.3), ncol = 3),
  nrow = 2, rel_heights = c(0.1,1))
  
  print(full_plot)
  
  # mat_80 <- anomaly_80_array[,,i]
  # mat_90 <- anomaly_90_array[,,i]
  # 
  # colnames(mat_80) <- dist_80
  # rownames(mat_80) <- depth_80
  # 
  # colnames(mat_90) <- dist_90
  # rownames(mat_90) <- depth_90
  # 
  # plot_df_80_anomaly <- melt(mat_80, varnames = c("Depth", "Distance"))
  # plot_df_90_anomaly <- melt(mat_90, varnames = c("Depth", "Distance"))
  # 
  # 
  # 
  # dens_ano_plot_80 <- ggplot(plot_df_80_anomaly, aes(y = Depth, x = Distance, fill = value)) + 
  #   geom_tile() + 
  #   scale_fill_gradient2(low = "blue", mid = "white", high = "red",
  #                        midpoint = 0, na.value = "grey") +
  #   geom_point(data = subset_dat_80, aes(x = dist_to_coast, y = Depthm, shape = as.factor(som_id)),
  #              stat = "identity", inherit.aes = FALSE) + 
  #   scale_shape_manual(values = c(17,19), guide = FALSE) +
  #   scale_x_reverse() + scale_y_reverse() +
  #   xlab("Distance to Coast (km)") + ylab("Depth (m)") +
  #   labs(fill = "Potential\nDensity Anomaly", shape = "SOM ID") +
  #   theme(panel.background = element_blank(),
  #         plot.title = element_text(hjust= 0.5)) +
  #   ggtitle("Line 80")
  # 
  # dens_ano_plot_90 <- ggplot(plot_df_90_anomaly, aes(y = Depth, x = Distance, fill = value)) + 
  #   geom_tile() + 
  #   scale_fill_gradient2(low = "blue", mid = "white", high = "red",
  #                        midpoint = 0, na.value = "grey") +
  #   geom_point(data = subset_dat_90, aes(x = dist_to_coast, y = Depthm, shape = as.factor(som_id)),
  #              stat = "identity", inherit.aes = FALSE) + 
  #   scale_shape_manual(values = c(17,19), guide = FALSE) +
  #   scale_x_reverse() + scale_y_reverse() +
  #   xlab("Distance to Coast (km)") + ylab("Depth (m)") +
  #   labs(fill = "Potential\nDensity Anomaly", shape = "SOM ID") +
  #   theme(panel.background = element_blank(),
  #         plot.title = element_text(hjust= 0.5)) +
  #   ggtitle("Line 90")
  # 
  # 
  # anomaly <- plot_grid(dens_ano_plot_80, dens_ano_plot_90, ncol = 2)
  
  # print(anomaly)
  
  }








}


