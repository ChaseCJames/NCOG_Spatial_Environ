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


###### Figure Code #####

som_figure <- function(SST_file = "output/CALCOFI_temp_tables.Rdata",
                       map_file = "output/bacteria_16s_map.Rdata",
                       reg_file = "output/bacteria_16s_glm.Rdata",
                       figure_name = paste0("figures/bacteria_16s_som_",Sys.Date(),".pdf"),
                       som_color_vector = c("darkred", "darkblue"),
                       som_name_vector = c("Nearshore", "Offshore"),
                       centroid_color = c("blue","red"), switch = TRUE,
                       main = "16s Bacteria", mean_position = c(15,0.5),
                       coeff_position = c(0.123,0.5)){
  
  
  load(SST_file)
  load(map_file)
  load(reg_file)
  
  map <- map_data("world")   
  
  sst <- ggplot() + 
    geom_tile(data = coeff_table, aes(x = lon, y = lat, fill = coeff_var), width =0.26, height = 0.26) +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", limits = c(0.09,0.12), oob = squish, midpoint = 0.1066851) +geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    ggtitle("A. SST Coeff. Var.") +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          title = element_text(hjust = 0.5), axis.line = element_blank())
  
  sst_mean <- ggplot() + 
    geom_tile(data = mean_table, aes(x = lon, y = lat, fill = Mean), width =0.26, height = 0.26) +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", limits = c(15,18), oob = squish, midpoint = 16.5) +geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    ggtitle("B. SST Mean (°C)") +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          title = element_text(hjust = 0.5), axis.line = element_blank())
  
  # find centroids
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  if(switch == TRUE){Letter_1 <- "D. "}else{Letter_1 <- "C. "}
  
  p1 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",1)), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = som_color_vector[1], limits = c(0,1)) +
    ggtitle(paste0(Letter_1,som_name_vector[1]," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          title = element_text(hjust = 0.5), axis.line = element_blank())
  
  centroid1 <- wt.centroid(x =centroid_df , p = 2)
  p1 <- p1 + geom_point(aes(x = centroid1@coords[1], y = centroid1@coords[2]), color = centroid_color[1], size = 5, pch = 10)
  
  
  if(switch == TRUE){Letter_2 <- "C. "}else{Letter_2 <- "D. "}
  
  p2 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",2)), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = som_color_vector[2], limits = c(0,1)) +
    ggtitle(paste0(Letter_2, som_name_vector[2]," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          title = element_text(hjust = 0.5), axis.line = element_blank())
  
  centroid2 <- wt.centroid(x =centroid_df , p = 3)
  p2 <- p2 + geom_point(aes(x = centroid2@coords[1], y = centroid2@coords[2]), color = centroid_color[2], size = 5, pch = 10)
  
  som_plots$cluster[som_plots$cluster == "som_1"] <- som_name_vector[1]
  som_plots$cluster[som_plots$cluster == "som_2"] <- som_name_vector[2]
  
  summary_fit <- summary(lm(som_maps$som_1~som_maps$sst_mean))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_mean <- summary_fit$r.squared
  
  summary_fit <- summary(lm(som_maps$som_1~som_maps$sst_coeff))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_coeff <- summary_fit$r.squared
  
  
  mean_plot <- ggplot(som_plots, aes(x = sst_mean, y = freq, color = cluster)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid")) +
    stat_smooth(method = "lm", level = 0.95) +
    scale_color_manual(values = c("red", "blue")) +
    xlab("Mean SST (°C)") + ylab("Frequency") +
    annotate("text", x = mean_position[1], y = mean_position[2], label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 3) + ggtitle("F. SOM Freq ~ Mean SST")
  
  coeff_plot <- ggplot(som_plots, aes(x = sst_coeff, y = freq, color = cluster)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none") +
    stat_smooth(method = "lm", level = 0.95) +
    scale_color_manual(values = c("red", "blue")) +
    xlab("SST Coeff. Var") + ylab("Frequency") +
    annotate("text", x = coeff_position[1], y = coeff_position[2], label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 3) + ggtitle("E. SOM Freq ~ Coeff. Var. SST")
  
  temp_som_plots <- plot_grid(coeff_plot, mean_plot, ncol = 2,rel_widths = c(0.85,1))
  
  title <- ggdraw() + draw_label(main, fontface='bold')
  theme_set(theme_cowplot(font_size=8))
  
  
  if(switch == TRUE){
    pdf(file = figure_name, width = 9, height = 12)
    a <- plot_grid(title, plot_grid(sst,sst_mean, ncol = 2), 
              plot_grid(p2,p1, ncol = 2),
              temp_som_plots, ncol = 1, nrow = 4, rel_heights = c(0.1,1,1,1))
    print(a)
    dev.off()
    }
  if(switch != TRUE){
    pdf(file = figure_name, width = 9, height = 12)
    a <- plot_grid(title, plot_grid(sst,sst_mean, ncol = 2), 
              plot_grid(p1,p2, ncol = 2),
              temp_som_plots, ncol = 1, nrow = 4, rel_heights = c(0.1,1,1,1))
    print(a)
    dev.off()
  }
  
}

aic_table <- function(map_file = "output/bacteria_16s_map.Rdata",
                      figure_name = paste0("figures/bacteria_16s_aic_",Sys.Date(),".pdf"),
                      minimum_tp = 8){
  
  load(map_file)
  
  som_maps <- som_maps[which(som_maps$n_samps > (minimum_tp-1)),]
  
  # trying to fit a glm
  
  som_glm <- som_maps[,c(1:9,16,10:12,17)]

  glm_mean_var <- glm(som_1 ~ PO4_mean + NO3_mean + SiO3_mean + temp_mean + 
                        PO4_coeff + NO3_coeff + SiO3_coeff + temp_coeff, data = som_glm)
  
  step_AIC_result <- stepAIC(glm_mean_var)
  
  model_AIC <- as.data.frame(matrix(NA,10,2))
  colnames(model_AIC) <- c("Model","AIC")
  
  # Temp
  
  glm_mean_temp <- glm(som_1 ~  temp_mean, data = som_glm)
  mt_sum <- summary(glm_mean_temp)
  model_AIC[1,2] <- mt_sum$aic
  model_AIC[1,1] <- "Mean Temp"
  
  glm_coeff_temp <- glm(som_1 ~  temp_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_temp)
  model_AIC[2,2] <- ct_sum$aic
  model_AIC[2,1] <- "Coeff Var. Temp"
  
  # NO3
  
  glm_mean_no3 <- glm(som_1 ~  NO3_mean, data = som_glm)
  mn_sum <- summary(glm_mean_no3)
  model_AIC[3,2] <- mn_sum$aic
  model_AIC[3,1] <- "Mean NO3"
  
  glm_coeff_no3 <- glm(som_1 ~  NO3_coeff, data = som_glm)
  cn_sum <- summary(glm_coeff_no3)
  model_AIC[4,2] <- cn_sum$aic
  model_AIC[4,1] <- "Coeff Var. NO3"
  
  # PO4
  
  glm_mean_po4 <- glm(som_1 ~  PO4_mean, data = som_glm)
  mp_sum <- summary(glm_mean_po4)
  model_AIC[5,2] <- mp_sum$aic
  model_AIC[5,1] <- "Mean PO4"
  
  glm_coeff_po4 <- glm(som_1 ~  PO4_coeff, data = som_glm)
  cp_sum <- summary(glm_coeff_po4)
  model_AIC[6,2] <- cp_sum$aic
  model_AIC[6,1] <- "Coeff Var. PO4"
  
  # SiO3
  
  glm_mean_sio3 <- glm(som_1 ~  SiO3_mean, data = som_glm)
  ms_sum <- summary(glm_mean_sio3)
  model_AIC[7,2] <- ms_sum$aic
  model_AIC[7,1] <- "Mean SiO3"
  
  glm_coeff_sio3 <- glm(som_1 ~  SiO3_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_sio3)
  model_AIC[8,2] <- cs_sum$aic
  model_AIC[8,1] <- "Coeff Var. SiO3"
  
  # Everything
  
  glm_mean_var <- glm(som_1 ~ PO4_mean + NO3_mean + SiO3_mean + temp_mean + 
                        PO4_coeff + NO3_coeff + SiO3_coeff + temp_coeff, data = som_glm)
  
  all_sum <- summary(glm_mean_var)
  model_AIC[9,2] <- all_sum$aic
  model_AIC[9,1] <- "Full Model"
  
  # best fit
  glm_simple <- glm(step_AIC_result$formula, data = som_glm)
  
  simple_sum <- summary(glm_simple)
  model_AIC[10,2] <- simple_sum$aic
  model_AIC[10,1] <- as.character(step_AIC_result$formula)[3]
  
  AIC_table <- model_AIC[order(model_AIC$AIC, decreasing = TRUE),]
  
  AIC_table$AIC <- round(AIC_table$AIC, 3)
  aic <- tableGrob(AIC_table, rows = NULL)
  
  pdf(file = figure_name, height = 5, width = 7)
  grid.draw(aic)
  dev.off()
  
  
}

random_forest_figure <- function(full_data_file = "output/bacteria_16s_full_data.Rdata",
                                 figure_name = paste0("figures/bacteria_16s_rf_",Sys.Date(),".pdf"),
                                 main = "16s Bacteria"){
  
  
  load(full_data_file)
  
  forest_dat <- full_dat[,-which(is.na(rowSums(full_dat[,c(32,33,37,38,39,59)])))]
  
  forest_dat$som_id <- as.factor(forest_dat$som_id)
  
  train <- sample(nrow(forest_dat), 0.7*nrow(forest_dat), replace = FALSE)
  
  forest_out <- randomForest(som_id ~  CC_Depth + T_degC + PO4ug + SiO3ug +
                               NO3ug, data = forest_dat, importance = TRUE,
                             na.action = na.omit, subset = train, mtry = 4)
  
  
  # Predicting on validation set
  predTrain <- predict(forest_out, forest_dat[train,], type = "class")
  # Checking classification accuracy
  table(predTrain, forest_dat$som_id[train]) 
  
  # Predicting on validation set
  predValid <- predict(forest_out, forest_dat[-train,], type = "class")
  # Checking classification accuracy
  table(predValid, forest_dat$som_id[-train])  
  
  plot <- varImpPlot(forest_out, type = 2, main = paste0(main," Random Forest"))
  
  pdf(file = figure_name, width = 5, height = 4)
  varImpPlot(forest_out, type = 2, main = paste0(main," Random Forest"))
  dev.off()
  
  
}
  
diversity_figure <- function(map_file = "output/bacteria_16s_map.Rdata",
                             full_data_file = "output/bacteria_16s_full_data.Rdata",
                             reg_file = "output/bacteria_16s_glm.Rdata",
                             figure_name = paste0("figures/bacteria_16s_diversity_",Sys.Date(),".pdf"),
                             minimum_tp = 8,
                             div_mean_position = c(15,4.73), div_coeff_position = c(0.125,4.78),
                             even_mean_position = c(15.2,0.73), even_coeff_position = c(0.128,0.73),
                             rich_mean_position = c(14.7,500), rich_coeff_position = c(0.123,550),
                             shannon_limits = c(3.8,5.4), even_limits = c(0.6,0.8), rich_limits = c(300,500),
                             main = "16s Bacteria"){
  
  load(map_file)
  load(full_data_file)
  load(reg_file)
  
  map <- map_data("world") 
  
  # shannon
  shannon <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "shannon"), color = "black", size =5, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "blue", aes(min = min(shannon, na.rm = TRUE, max = max(shannon, na.rm = TRUE)), limits = c(min,max))) +
    ggtitle(paste0("Shannon Diversity")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          title = element_text(hjust = 0.5))
  
  evenness <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "evenness"), color = "black", size =5, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "red", aes(min = min(evenness, na.rm = TRUE, max = max(evenness, na.rm = TRUE)), limits = c(min,max))) +
    ggtitle(paste0("Evenness")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          title = element_text(hjust = 0.5))
  
  richness <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "richness"), color = "black", size =5, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "darkgreen", aes(min = min(richness, na.rm = TRUE, max = max(richness, na.rm = TRUE)), limits = c(min,max))) +
    ggtitle(paste0("Species Richness")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          title = element_text(hjust = 0.5))
  
  # glms
  
  som_maps <- som_maps[which(som_maps$n_samps > (minimum_tp-1)),]
  
  summary_fit <- summary(lm(som_maps$shannon~som_maps$temp_mean))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_mean <- summary_fit$r.squared
  
  summary_fit <- summary(lm(som_maps$shannon~som_maps$temp_coeff))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_coeff <- summary_fit$r.squared
  
  mean_plot <- ggplot(som_plots, aes(x = temp_mean, y = shannon)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid")) +
    stat_smooth(method = "lm", level = 0.95) +
    xlab("Mean SST (°C)") + ylab("Shannon Index") +
    annotate("text", x = div_mean_position[1], y = div_mean_position[2], label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 4)
  
  coeff_plot <- ggplot(som_plots, aes(x = temp_coeff, y = shannon)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none") +
    stat_smooth(method = "lm", level = 0.95) +
    xlab("SST Coeff. Var") + ylab("Shannon Index") +
    annotate("text", x = div_coeff_position[1], y = div_coeff_position[2], label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 4)
  
  theme_set(theme_cowplot(font_size=8))
  temp_diversity_plots <- plot_grid(coeff_plot, mean_plot, ncol = 1,rel_widths = c(1,1))
  
  # Evenness
  
  summary_fit <- summary(lm(som_maps$evenness~som_maps$temp_mean))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_mean <- summary_fit$r.squared
  
  summary_fit <- summary(lm(som_maps$evenness~som_maps$temp_coeff))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_coeff <- summary_fit$r.squared
  
  mean_plot <- ggplot(som_plots, aes(x = temp_mean, y = evenness)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid")) +
    stat_smooth(method = "lm", level = 0.95) +
    xlab("Mean SST (°C)") + ylab("Evenness") +
    annotate("text", x = even_mean_position[1], y = even_mean_position[2], label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 4)
  
  coeff_plot <- ggplot(som_plots, aes(x = temp_coeff, y = evenness)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none") +
    stat_smooth(method = "lm", level = 0.95) +
    xlab("SST Coeff. Var") + ylab("Evenness") +
    annotate("text", x = even_coeff_position[1], y = even_coeff_position[2], label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 4)
  
  theme_set(theme_cowplot(font_size=8))
  temp_evenness_plots <- plot_grid(coeff_plot, mean_plot, ncol = 1,rel_widths = c(1,1))
  
  # Richness
  
  summary_fit <- summary(lm(som_maps$richness~som_maps$temp_mean))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_mean <- summary_fit$r.squared
  
  summary_fit <- summary(lm(som_maps$richness~som_maps$temp_coeff))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_coeff <- summary_fit$r.squared
  
  mean_plot <- ggplot(som_plots, aes(x = temp_mean, y = richness)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid")) +
    stat_smooth(method = "lm", level = 0.95) +
    xlab("Mean SST (°C)") + ylab("Richness") +
    annotate("text", x = rich_mean_position[1], y = rich_mean_position[2], label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 4)
  
  coeff_plot <- ggplot(som_plots, aes(x = temp_coeff, y = richness)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none") +
    stat_smooth(method = "lm", level = 0.95) +
    xlab("SST Coeff. Var") + ylab("Richness") +
    annotate("text", x = rich_coeff_position[1], y = rich_coeff_position[2], label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 4)
  
  theme_set(theme_cowplot(font_size=8))
  temp_richness_plots <- plot_grid(coeff_plot, mean_plot, ncol = 1,rel_widths = c(1,1))
  
  shannon_depth_dist <- ggplot(full_dat, aes(x = dist_to_coast, y = CC_Depth, fill = shannon)) + 
    geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
    scale_fill_viridis(name = "Shannon Diversity", option = "D", limits = shannon_limits, oob = scales::squish) +
    ylab("Depth (m)") + xlab("Distance to Coast (km)") 
  
  even_depth_dist <- ggplot(full_dat, aes(x = dist_to_coast, y = CC_Depth, fill = evenness)) + 
    geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
    scale_fill_viridis(name = "Evenness", option = "D", limits = even_limits, oob = scales::squish) +
    ylab("Depth (m)") + xlab("Distance to Coast (km)") 
  
  rich_depth_dist <- ggplot(full_dat, aes(x = dist_to_coast, y = CC_Depth, fill = richness)) + 
    geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
    scale_fill_viridis(name = "Richness", option = "D", limits = rich_limits, oob = scales::squish) +
    ylab("Depth (m)") + xlab("Distance to Coast (km)") 
  
  
  title1 <- ggdraw() + draw_label(paste0(main, " Diversity/Evenness/Richness"), fontface='bold')
  pdf(file = figure_name, width = 12, height = 12)
  a <- plot_grid(title1, plot_grid(shannon,evenness,richness, ncol = 3, labels = c("A","B","C")),
            plot_grid(temp_diversity_plots, temp_evenness_plots, temp_richness_plots, ncol = 3, labels = c("D","E","F")),
            plot_grid(shannon_depth_dist, even_depth_dist, rich_depth_dist, ncol = 3, labels = c("G","H", "I")),
            ncol = 1, nrow = 4, rel_heights = c(0.1,1,2,1))
  print(a)
  dev.off()
  
}

alpha_versus_gamma_figure <- function(full_data_file = "output/bacteria_16s_full_data.Rdata",
                                      raw_data_file = "data/16s_bacteria.Rdata",
                                      map_file = "output/bacteria_16s_map.Rdata", minimum_tp = 8,
                                      figure_name = paste0("figures/bacteria_16s_diversity_",Sys.Date(),".pdf"),
                                      main = "16s Bacteria", gamma_position = c(0.128,5.7), alpha_position = c(0.10,4.3)){
  
  load(full_data_file)
  load(raw_data_file)  
  load(map_file)
  
  som_maps <- som_maps[which(som_maps$n_samps > (minimum_tp-1)),]
  
  # Gamma Diversity
  
  asv_sums <- scaled_inputs
  
  asv_sums <- as.data.frame(asv_sums)
  
  asv_sums$station <- full_dat$Sta_ID[match(rownames(asv_sums), full_dat$eco_name)]
  
  station_sums <- asv_sums %>% 
    group_by(station) %>% 
    summarise_all(mean, na.rm = TRUE)
  
  station_sums <- as.data.frame(station_sums)
  
  rownames(station_sums) <- station_sums$station
  
  station_sums$station <- NULL
  
  station_plot_df <- as.data.frame(matrix(NA,NROW(station_sums),4))
  
  colnames(station_plot_df) <- c("Station", "Diversity", "Latitude", "Longitude")
  
  station_plot_df$Station <- rownames(station_sums)
  
  # diversity
  
  station_plot_df$Diversity <- diversity(station_sums, MARGIN = 1, index = "shannon")
  
  # station location
  
  station_position <- full_dat %>% 
    group_by(Sta_ID) %>% 
    summarise(Lat = mean(Lat_Dec, na.rm = TRUE), Lon = mean(Lon_Dec, na.rm = TRUE))
  
  station_plot_df$Latitude <- station_position$Lat[match(station_plot_df$Station, station_position$Sta_ID)]
  station_plot_df$Longitude <- station_position$Lon[match(station_plot_df$Station, station_position$Sta_ID)]
  
  # remove northern transects
  
  line <- as.numeric(substr(station_plot_df$Station,2,3))
  n_stations <- which(line < 76)
  if(length(n_stations) > 0){station_plot_df <- station_plot_df[-n_stations,]}
  
  # Alpha Gamma Plots
  
  summary_fit <- summary(lm(som_maps$shannon~som_maps$temp_coeff))
  
  alpha_coeff_p_val <- summary_fit$coefficients[2,4]
  alpha_r_sq_coeff <- summary_fit$r.squared
  
  som_maps$Gamma_Diversity <- station_plot_df$Diversity[match(som_maps$Sta_ID, station_plot_df$Station)]
  
  summary_fit <- summary(lm(som_maps$Gamma_Diversity~som_maps$temp_coeff))
  
  gamma_coeff_p_val <- summary_fit$coefficients[2,4]
  gamma_r_sq_coeff <- summary_fit$r.squared
  
  div_plot <- ggplot(som_maps, aes(x = temp_coeff)) + 
    geom_point(aes(y = shannon, color = "Alpha Diversity")) +
    stat_smooth(aes(y = shannon), method = "lm", level = 0.95, color = "black") +
    geom_point(aes(y = Gamma_Diversity, color = "Gamma Diversity")) +
    stat_smooth(aes(y = Gamma_Diversity), method = "lm", level = 0.95, color = "black") +
    scale_y_continuous(sec.axis = sec_axis(~., name = "Gamma Diversity")) +
    ylab("Alpha Diversity") + xlab("Coeff. Var. SST") +
    theme(legend.title = element_blank(),
          legend.box.background = element_rect(color = "black"),
          panel.background = element_blank()) +
    ggtitle(paste0(main," Mean Alpha Diversity vs Gamma Diversity per Station")) +
    scale_color_manual(values = c("royalblue2", "seagreen3")) +
    annotate("text", x = alpha_position[1], y = alpha_position[2], label = paste0("R Squared = ", round(alpha_r_sq_coeff, digits = 3)), size = 4) +
    annotate("text", x = gamma_position[1], y = gamma_position[2], label = paste0("R Squared = ", round(gamma_r_sq_coeff, digits = 3)), size = 4)
  
  theme_set(theme_cowplot(font_size=10))
  
  pdf(file = figure_name, width = 8, height = 6)
  print(div_plot)
  dev.off()
  
}

diversity_over_time_figure <- function(full_data_file = "output/bacteria_16s_full_data.Rdata",
                                       figure_name = paste0("figures/bacteria_16s_diversity_time_",Sys.Date(),".pdf"),
                                       shannon_limits = c(3.8,5.4), even_limits = c(0.6,0.8), rich_limits = c(300,500),
                                       main = "16s Bacteria"){
  
  load(full_data_file)
  
  shannon_timeplot <- ggplot(full_dat, aes(x = Date, y = dist_to_coast, fill = shannon)) + 
    geom_jitter(size = 4, pch = 21, color = "black",height = 0, width = 10) + 
    scale_fill_viridis(name = "Shannon", option = "D", limits = shannon_limits, oob = scales::squish) +
    ylab("Distance to Coast") + xlab("Time") 
  
  even_timeplot <- ggplot(full_dat, aes(x = Date, y = dist_to_coast, fill = evenness)) + 
    geom_jitter(size = 4, pch = 21, color = "black",height = 0, width = 10) + 
    scale_fill_viridis(name = "Evenness", option = "D", limits = even_limits, oob = scales::squish) +
    ylab("Distance to Coast") + xlab("Time") 
  
  rich_timeplot <- ggplot(full_dat, aes(x = Date, y = dist_to_coast, fill = richness)) + 
    geom_jitter(size = 4, pch = 21, color = "black",height = 0, width = 10) + 
    scale_fill_viridis(name = "Richness", option = "D", limits = rich_limits, oob = scales::squish) +
    ylab("Distance to Coast") + xlab("Time") 
  
  title2 <- ggdraw() + draw_label(paste0(main," Diversity/Evenness/Richness over Time"), fontface='bold')
  plot <- plot_grid(title2, shannon_timeplot, even_timeplot, rich_timeplot,
            ncol = 1, nrow = 4, rel_heights = c(0.1,1,1,1))
  
  pdf(file = figure_name, width = 8, height = 12)
  print(plot) 
  dev.off()
  
}

aic_table_func <- function(som_maps = cyano_plots){  
  
  som_glm <- som_maps[,c(1:3,7:22,26:29)]
  
  model_AIC <- matrix(NA,18,2)
  
  # temperature
  
  glm_mean_temp <- glm(som_1 ~  temp_mean, data = som_glm)
  mt_sum <- summary(glm_mean_temp)
  model_AIC[1,2] <- mt_sum$aic
  model_AIC[1,1] <- "Mean Temp"
  
  glm_coeff_temp <- glm(som_1 ~  temp_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_temp)
  model_AIC[10,2] <- ct_sum$aic
  model_AIC[10,1] <- "Coeff. Var. Temp"
  
  # sea surface temperature
  
  glm_mean_sst <- glm(som_1 ~  sst_mean, data = som_glm)
  mt_sum <- summary(glm_mean_sst)
  model_AIC[2,2] <- mt_sum$aic
  model_AIC[2,1] <- "Mean SST"
  
  glm_coeff_sst <- glm(som_1 ~  sst_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_sst)
  model_AIC[11,2] <- ct_sum$aic
  model_AIC[11,1] <- "Coeff. Var. SST"
  
  # salinity
  
  glm_mean_sal <- glm(som_1 ~  sal_mean, data = som_glm)
  mt_sum <- summary(glm_mean_sal)
  model_AIC[3,2] <- mt_sum$aic
  model_AIC[3,1] <- "Mean Salinity"
  
  glm_coeff_sal <- glm(som_1 ~  sal_coeff, data = som_glm)
  ct_sum <- summary(glm_coeff_sal)
  model_AIC[12,2] <- ct_sum$aic
  model_AIC[12,1] <- "Coeff. Var. Salinity"
  
  # NO3
  
  glm_mean_no3 <- glm(som_1 ~  NO3_mean, data = som_glm)
  mn_sum <- summary(glm_mean_no3)
  model_AIC[4,2] <- mn_sum$aic
  model_AIC[4,1] <- "Mean NO3"
  
  glm_coeff_no3 <- glm(som_1 ~  NO3_coeff, data = som_glm)
  cn_sum <- summary(glm_coeff_no3)
  model_AIC[13,2] <- cn_sum$aic
  model_AIC[13,1] <- "Coeff. Var. NO3"
  
  # PO4
  
  glm_mean_po4 <- glm(som_1 ~  PO4_mean, data = som_glm)
  mp_sum <- summary(glm_mean_po4)
  model_AIC[5,2] <- mp_sum$aic
  model_AIC[5,1] <- "Mean PO4"
  
  glm_coeff_po4 <- glm(som_1 ~  PO4_coeff, data = som_glm)
  cp_sum <- summary(glm_coeff_po4)
  model_AIC[14,2] <- cp_sum$aic
  model_AIC[14,1] <- "Coeff. Var. PO4"
  
  # SiO3
  
  glm_mean_sio3 <- glm(som_1 ~  SiO3_mean, data = som_glm)
  ms_sum <- summary(glm_mean_sio3)
  model_AIC[6,2] <- ms_sum$aic
  model_AIC[6,1] <- "Mean SiO3"
  
  glm_coeff_sio3 <- glm(som_1 ~  SiO3_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_sio3)
  model_AIC[15,2] <- cs_sum$aic
  model_AIC[15,1] <- "Coeff. Var. SiO3"
  
  # C14
  
  # glm_mean_c14 <- glm(som_1 ~  C14_mean, data = som_glm)
  # ms_sum <- summary(glm_mean_c14)
  # model_AIC[7,2] <- ms_sum$aic
  # model_AIC[7,1] <- "Mean C14"
  # 
  # glm_coeff_c14 <- glm(som_1 ~  C14_coeff, data = som_glm)
  # cs_sum <- summary(glm_coeff_c14)
  # model_AIC[17,2] <- cs_sum$aic
  # model_AIC[17,1] <- "Coeff. Var. C14"
  
  # SLA
  
  glm_mean_sla <- glm(som_1 ~  sla_mean, data = som_glm)
  ms_sum <- summary(glm_mean_sla)
  model_AIC[7,2] <- ms_sum$aic
  model_AIC[7,1] <- "Mean SLA"
  
  glm_coeff_sla <- glm(som_1 ~  sla_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_sla)
  model_AIC[16,2] <- cs_sum$aic
  model_AIC[16,1] <- "Coeff. Var. SLA"
  
  # MLD
  
  glm_mean_mld <- glm(som_1 ~  MLD_mean, data = som_glm)
  ms_sum <- summary(glm_mean_mld)
  model_AIC[8,2] <- ms_sum$aic
  model_AIC[8,1] <- "Mean MLD"
  
  glm_coeff_mld <- glm(som_1 ~  MLD_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_mld)
  model_AIC[17,2] <- cs_sum$aic
  model_AIC[17,1] <- "Coeff. Var. MLD"
  
  # NC Depth
  
  glm_mean_nc <- glm(som_1 ~  NC_mean, data = som_glm)
  ms_sum <- summary(glm_mean_nc)
  model_AIC[9,2] <- ms_sum$aic
  model_AIC[9,1] <- "Mean NCD"
  
  glm_coeff_nc <- glm(som_1 ~  NC_coeff, data = som_glm)
  cs_sum <- summary(glm_coeff_nc)
  model_AIC[18,2] <- cs_sum$aic
  model_AIC[18,1] <- "Coeff. Var. NCD"
  
  return(model_AIC)
  
}

full_aic_table_figure <- function(in_plastid = "output/plastid_16s_map.Rdata",
                                  in_cyano = "output/cyano_16s_map.Rdata",
                             in_bacteria = "output/bacteria_m_euks_16s_map.Rdata",
                             in_euks = "output/euks_hetero_18sv9_map.Rdata",
                             in_phyto = "output/euks_auto_18sv9_map.Rdata", minimum_tp = 8,
                             figure_name = paste0("figures/full_aic_table_",Sys.Date(),".pdf"),
                             figure_name_2 = paste0("figures/full_aic_plot_",Sys.Date(),".pdf"),
                             title_name = "Variable Importance"){
  # load an rename data
  
  load(in_plastid)
  som_maps <- som_maps[which(som_maps$n_samps > (minimum_tp-1)),]
  plastid_plots <- som_maps
  
  load(in_cyano)
  som_maps <- som_maps[which(som_maps$n_samps > (minimum_tp-1)),]
  cyano_plots <- som_maps
  
  load(in_bacteria)
  som_maps <- som_maps[which(som_maps$n_samps > (minimum_tp-1)),]
  bacteria_plots <- som_maps
  
  load(in_euks)
  som_maps <- som_maps[which(som_maps$n_samps > (minimum_tp-1)),]
  eukaryota_plots <- som_maps
  
  load(in_phyto) 
  som_maps <- som_maps[which(som_maps$n_samps > (minimum_tp-1)),]
  phyto_plots <- som_maps
  
  plastid_AIC <- aic_table_func(som_maps = plastid_plots)
  plastid_AIC <- as.data.frame(plastid_AIC, stringsAsFactors = FALSE)
  colnames(plastid_AIC) <- c("Variables","AIC")
  plastid_AIC$AIC <- as.numeric(plastid_AIC$AIC)
  plastid_AIC$AIC <- round(plastid_AIC$AIC, digits = 2)
  
  cyano_AIC <- aic_table_func(som_maps = cyano_plots)
  cyano_AIC <- as.data.frame(cyano_AIC, stringsAsFactors = FALSE)
  colnames(cyano_AIC) <- c("Variables","AIC")
  cyano_AIC$AIC <- as.numeric(cyano_AIC$AIC)
  cyano_AIC$AIC <- round(cyano_AIC$AIC, digits = 2)
  
  bacteria_AIC <- aic_table_func(som_maps = bacteria_plots)
  bacteria_AIC <- as.data.frame(bacteria_AIC, stringsAsFactors = FALSE)
  colnames(bacteria_AIC) <- c("Variables","AIC")
  bacteria_AIC$AIC <- as.numeric(bacteria_AIC$AIC)
  bacteria_AIC$AIC <- round(bacteria_AIC$AIC, digits = 2)
  
  eukaryota_AIC <- aic_table_func(som_maps = eukaryota_plots)
  eukaryota_AIC <- as.data.frame(eukaryota_AIC, stringsAsFactors = FALSE)
  colnames(eukaryota_AIC) <- c("Variables","AIC")
  eukaryota_AIC$AIC <- as.numeric(eukaryota_AIC$AIC)
  eukaryota_AIC$AIC <- round(eukaryota_AIC$AIC, digits = 2)
  
  phyto_AIC  <- aic_table_func(som_maps = phyto_plots)
  phyto_AIC <- as.data.frame(phyto_AIC, stringsAsFactors = FALSE)
  colnames(phyto_AIC) <- c("Variables","AIC")
  phyto_AIC$AIC <- as.numeric(phyto_AIC$AIC)
  phyto_AIC$AIC <- round(phyto_AIC$AIC, digits = 2)
  

  AIC_full <- full_join(plastid_AIC, cyano_AIC, by = "Variables")
  AIC_full <- full_join(AIC_full, bacteria_AIC, by = "Variables")
  AIC_full <- full_join(AIC_full, eukaryota_AIC, by = "Variables")
  AIC_full <- full_join(AIC_full, phyto_AIC, by = "Variables")
  
  colnames(AIC_full) <- c("Variables", "Eukaryotic Plastid\n AIC", "Cyanobacteria\n AIC",
                          "Bacteria/Archaea\n AIC", "Heterotrophic Eukaryotes\n AIC", "Eukaryotic Phytoplankton\n AIC")
  
  colfunc <- colorRampPalette(c("white", "red"))
  
  plastid_col <- round(seq(max(plastid_AIC$AIC, na.rm = TRUE),min(plastid_AIC$AIC, na.rm = TRUE), by = -0.01), digits = 2)
  scale <- colfunc(length(plastid_col))
  plastid_fill <- scale[match(round(plastid_AIC$AIC, digits = 2), plastid_col)]
  
  cyano_col <- round(seq(max(cyano_AIC$AIC, na.rm = TRUE),min(cyano_AIC$AIC, na.rm = TRUE), by = -0.01), digits = 2)
  scale <- colfunc(length(cyano_col))
  cyano_fill <- scale[match(round(cyano_AIC$AIC, digits = 2), cyano_col)]
  cyano_fill[5] <- scale[3903]
  
  bacteria_col <- round(seq(max(bacteria_AIC$AIC, na.rm = TRUE),min(bacteria_AIC$AIC, na.rm = TRUE), by = -0.01), digits = 2)
  scale <- colfunc(length(bacteria_col))
  bacteria_fill <- scale[match(round(bacteria_AIC$AIC, digits = 2), bacteria_col)]
  
  eukaryota_col <- round(seq(max(eukaryota_AIC$AIC, na.rm = TRUE),min(eukaryota_AIC$AIC, na.rm = TRUE), by = -0.01), digits = 2)
  scale <- colfunc(length(eukaryota_col))
  eukaryota_fill <- scale[match(round(eukaryota_AIC$AIC, digits = 2), eukaryota_col)]
  
  phyto_col <- round(seq(max(phyto_AIC$AIC, na.rm = TRUE),min(phyto_AIC$AIC, na.rm = TRUE), by = -0.01), digits = 2)
  scale <- colfunc(length(phyto_col))
  phyto_fill <- scale[match(round(phyto_AIC$AIC, digits = 2), phyto_col)]

  
  t0 <- tableGrob(AIC_full["Variables"], 
                  theme=ttheme_default(
                    core=list(bg_params = list(fill="grey90", col = "black"),
                              fg_params = list(fontface="bold")),
                    colhead = list(bg_params=list(fill="white", col="black"))), 
                  rows = NULL)
  t1 <- tableGrob(AIC_full["Eukaryotic Plastid\n AIC"], 
                  theme=ttheme_default(
                    core=list(bg_params = list(fill=plastid_fill, col = "black")),
                    colhead = list(bg_params=list(fill="white", col="black"))), 
                  rows = NULL)
  t2 <- tableGrob(AIC_full["Cyanobacteria\n AIC"],
                  theme=ttheme_default(
                    core=list(bg_params = list(fill=cyano_fill, col = "black")),
                    colhead = list(bg_params=list(fill="white", col="black"))),
                  rows = NULL)
  t3 <- tableGrob(AIC_full["Bacteria/Archaea\n AIC"],
                  theme=ttheme_default(
                    core=list(bg_params = list(fill=bacteria_fill, col = "black")),
                    colhead = list(bg_params=list(fill="white", col="black"))),
                  rows = NULL)
  t4 <- tableGrob(AIC_full["Heterotrophic Eukaryotes\n AIC"],
                  theme=ttheme_default(
                    core=list(bg_params = list(fill=eukaryota_fill, col = "black")),
                    colhead = list(bg_params=list(fill="white", col="black"))),
                  rows = NULL)
  t5 <- tableGrob(AIC_full["Eukaryotic Phytoplankton\n AIC"],
                  theme=ttheme_default(
                    core=list(bg_params = list(fill=phyto_fill, col = "black")),
                    colhead = list(bg_params=list(fill="white", col="black"))),
                  rows = NULL)
  
  # join tables
  tab <- gtable_combine(t0,t1,t2,t3,t4,t5)
  
  pdf(file = figure_name, width = 12, height = 6)
  grid.arrange(tab)
  dev.off()
  
  AIC_scaled <- AIC_full
  
  for (i in 2:ncol(AIC_scaled)) {
    zero_one_scale <- 1-(AIC_scaled[,i]-min(AIC_scaled[,i], na.rm = TRUE))/
      abs(min(AIC_scaled[,i], na.rm = TRUE) - max(AIC_scaled[,i], na.rm = TRUE))
    AIC_scaled[,i] <- 50^zero_one_scale
    
  }
  
  plot_df <- melt(AIC_scaled)
  
  colnames(plot_df) <- c("Variables", "Group", "AIC")
  plot_df$Variables <- as.factor(plot_df$Variables)
  plot_df$Variables <- factor(plot_df$Variables, levels = c("Coeff. Var. NCD","Coeff. Var. MLD",
                                                            "Coeff. Var. SLA",#"Coeff. Var. C14",
                                                            "Coeff. Var. SiO3","Coeff. Var. PO4",
                                                            "Coeff. Var. NO3", "Coeff. Var. Salinity",
                                                            "Coeff. Var. SST","Coeff. Var. Temp",
                                                            "Mean NCD", "Mean MLD",
                                                            "Mean SLA",#"Mean C14",
                                                            "Mean SiO3", "Mean PO4",
                                                            "Mean NO3","Mean Salinity",
                                                            "Mean SST", "Mean Temp" ))
  
  pdf(figure_name_2, width = 8, height = 8)
  print(ggplot(data = plot_df, aes(x = Group, y = Variables, size = AIC)) + 
    geom_point(fill = "red", color = "black", alpha = 0.6, shape = 21) +
    labs(size = "Variable\n Importance") + ylab("Variable") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          legend.position = "none",
          panel.grid.major.y = element_line(color = "grey", linetype = 2),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = -45)) + ggtitle(title_name) +
    scale_size_continuous(range = c(1,15)) + xlab("") + ylab(""))
  dev.off()
  
  
}


som_current_animation <- function(vel_table = "output/CALCOFI_Copernicus_GSC_Slice.Rdata",
                                  som_data = "output/bacteria_16s_full_data.Rdata", type = "bacteria_16s_vel",
                                  figure_name = paste0("figures/bacteria_16s_som_animation_",Sys.Date(),".gif")){
  
  load(vel_table)
  
  load(som_data)
  
  # split by cruise
  
  cruise_dates <- paste0("01-",substr(full_dat$Cruise,5,6),"-",substr(full_dat$Cruise,1,4))
  
  full_dat$Cruise_Date <- as.Date(cruise_dates, format = "%d-%m-%Y")
  
  # ts to date
  tmonth <- as.integer(1)
  tday <- as.integer(1)
  tyear <- as.integer(1950)
  ts <- chron(ts,origin=c(tmonth, tday, tyear))
  date_vect <- as.Date(ts, format = "%m/%d/%Y")
  
  unique_cruise_dates <- unique(full_dat$Cruise_Date)
  
  u_array_list <- array(NA,c(48,48,20))
  v_array_list <- array(NA,c(48,48,20))
  
  for (i in 1:length(unique_cruise_dates)) {
    
    # 3 month mean
    
    m_1 <- unique_cruise_dates[i] %m-% months(1)
    p_1 <- unique_cruise_dates[i] %m+% months(1)
    
    m_string <- c(m_1, unique_cruise_dates[i], p_1)
     
    z_vec <- which(!is.na(match(substr(date_vect,1,7),substr(m_string,1,7))))
    
  for (j in 1:48) {
      for (k in 1:48) {
        u_array_list[j,k,i] <- mean(u_slice[j,k,z_vec], na.rm = TRUE)
        v_array_list[j,k,i] <- mean(v_slice[j,k,z_vec], na.rm = TRUE)
      }
    }
  }
    
  u_mat <- u_array_list[,,1]
  colnames(u_mat) <- lon[lon_grid]
  rownames(u_mat) <- lat[lat_grid]
  
  mean_table_u <- melt(u_mat)
  
  mean_table_u$Cruise <- unique_cruise_dates[1]
  mean_table_u <- mean_table_u[seq(1,nrow(mean_table_u),3),]
  
  v_mat <- v_array_list[,,1]
  colnames(v_mat) <- lon[lon_grid]
  rownames(v_mat) <- lat[lat_grid]
  
  mean_table_v <- melt(v_mat)
  
  mean_table_v$Cruise <- unique_cruise_dates[1]
  mean_table_v <- mean_table_v[seq(1,nrow(mean_table_v),3),]
  
  for (i in 2:length(unique_cruise_dates)) {
    
    u_mat <- u_array_list[,,i]
    colnames(u_mat) <- lon[lon_grid]
    rownames(u_mat) <- lat[lat_grid]
    
    mean_table_u_i <- melt(u_mat)
    mean_table_u_i$Cruise <- unique_cruise_dates[i]
    mean_table_u_i <- mean_table_u_i[seq(1,nrow(mean_table_u_i),3),]
    
    v_mat <- v_array_list[,,i]
    colnames(v_mat) <- lon[lon_grid]
    rownames(v_mat) <- lat[lat_grid]
    
    mean_table_v_i <- melt(v_mat)
    mean_table_v_i$Cruise <- unique_cruise_dates[i]
    mean_table_v_i <- mean_table_v_i[seq(1,nrow(mean_table_v_i),3),]
    
    mean_table_u <- bind_rows(mean_table_u, mean_table_u_i)
    mean_table_v <- bind_rows(mean_table_v, mean_table_v_i)
  }
  
  colnames(mean_table_u) <- c("lat", "lon", "Mean_U", "Cruise")
  colnames(mean_table_v) <- c("lat", "lon", "Mean_V", "Cruise")
  
  uv_table <- inner_join(mean_table_u, mean_table_v, by = c("lon", "lat", "Cruise"))
  colnames(uv_table)[4] <- "Cruise_Date"
  uv_table$lon <- uv_table$lon - 360
  
  surf_dat <- full_dat[which(full_dat$CC_Depth < 15),]
  surf_dat <- surf_dat[-which(as.numeric(substr(surf_dat$Sta_ID,2,5)) < 76.7),]
  
  map <- map_data("world") 
  surf_dat$som_id <- as.character(surf_dat$som_id)
  
  
  for (i in 1:length(unique_cruise_dates)) {
    
    subset_surf <- surf_dat[which(surf_dat$Cruise_Date == unique_cruise_dates[i]),]
    subset_uv <- uv_table[which(uv_table$Cruise_Date == unique_cruise_dates[i]),]
    
    vel_som <- ggplot() + 
      geom_vector(data = subset_uv, aes(x = lon, y = lat, dx = Mean_U, dy = Mean_V),
                        arrow.angle = 15, arrow.type = "open", arrow.length = unit(0.2, "inches"),
                        pivot = 0,preserve.dir = TRUE, direction = "ccw",
                        min.mag = 0, show.legend = NA, color = "black") +
      geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
      coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
      geom_point(data = subset_surf, aes_string(x = "Lon_Dec", y = "Lat_Dec", fill = "som_id"), color = "black",
      size = 6, stroke = 0.1, shape = 21) +
      scale_fill_manual(name = "SOM", values = c("red","blue")) +
      theme(legend.position = "right",
            legend.background = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = "black", linetype = "solid", size = 1),
            title = element_text(hjust = 0.5), axis.line = element_blank()) +
      scale_mag(max = 0.1, name = "Speed (m/s)", max_size = 0.2) +
      labs(x = "Longitude", y = "Latitude") + ggtitle(subset_surf$Cruise_Date)
    if(i < 10){
      png(filename = paste0("figures/animations/",type,"_0",i,".png"), width = 6, height = 6, units = "in", res = 72)
      print(vel_som)
      dev.off()
      }
    if(i > 9){
      png(filename = paste0("figures/animations/",type,"_",i,".png"), width = 6, height = 6, units = "in", res = 72)
      print(vel_som)
      dev.off()
      }
    
    
    
      
    
  }
  
list.files(path = "figures/animations/", pattern = type, full.names = TRUE) %>%
  map(image_read) %>%
  image_join() %>%
  image_animate(fps = 0.5) %>%
  image_write(figure_name)
  

  
  
}

som_temp_animation <- function(temp_table = "output/CALCOFI_OI_SST_Slice.Rdata",
                                  som_data = "output/bacteria_16s_full_data.Rdata", type = "bacteria_16s_temp",
                                  figure_name = paste0("figures/bacteria_16s_som_animation_",Sys.Date(),".gif")){
  
  load(temp_table)
  
  load(som_data)
  
  # split by cruise
  
  cruise_dates <- paste0("01-",substr(full_dat$Cruise,5,6),"-",substr(full_dat$Cruise,1,4))
  
  full_dat$Cruise_Date <- as.Date(cruise_dates, format = "%d-%m-%Y")
  
  # ts to date
  tmonth <- as.integer(1)
  tday <- as.integer(1)
  tyear <- as.integer(1800)
  ts <- chron(ts_list,origin=c(tmonth, tday, tyear))
  date_vect <- as.Date(ts, format = "%m/%d/%Y")
  
  unique_cruise_dates <- unique(full_dat$Cruise_Date)
  
  t_array_list <- array(NA,c(48,48,20))
  
  # get length of each slice
  
  length_vector <- vector()
  seq_vector <- vector()
  
  for (i in 1:length(tmp_slice_list)) {
    
    length_vector[i] <- dim(tmp_slice_list[[i]])[3]  
    seq_vector <- c(seq_vector, rep(i,times = dim(tmp_slice_list[[i]])[3]))
  }
  
  for (i in 1:length(unique_cruise_dates)) {
    
    # 3 month mean
    
    m_1 <- unique_cruise_dates[i] %m-% months(1)
    p_1 <- unique_cruise_dates[i] %m+% months(1)
    
    m_string <- c(m_1, unique_cruise_dates[i], p_1)
    
    z_vec <- which(!is.na(match(substr(date_vect,1,7),substr(m_string,1,7))))
    
    list_subset <- seq_vector[z_vec]
    list_entry <- unique(list_subset)
    
    if(length(list_entry) == 1){
    
    days <- z_vec - sum(length_vector[1:(list_entry-1)])
    
    for (j in 1:48) {
      for (k in 1:48) {
         t_array_list[j,k,i] <- sd(tmp_slice_list[[list_entry]][j,k,days], na.rm = TRUE)/mean(tmp_slice_list[[list_entry]][j,k,days], na.rm = TRUE)
        }
      }
    }
    
    if(length(list_entry) == 2){
      
      
      days_1 <- z_vec[which(list_subset == list_entry[1])] - sum(length_vector[1:(list_entry[1]-1)]) 
      days_2 <- z_vec[which(list_subset == list_entry[2])] - sum(length_vector[1:(list_entry[2]-1)]) 
      
      for (j in 1:48) {
        for (k in 1:48) {
          t_array_list[j,k,i] <- sd(x = c(tmp_slice_list[[list_entry[1]]][j,k,days_1],tmp_slice_list[[list_entry[2]]][j,k,days_2]), na.rm = TRUE)/
            mean(x = c(tmp_slice_list[[list_entry[1]]][j,k,days_1],tmp_slice_list[[list_entry[2]]][j,k,days_2]), na.rm = TRUE)
        }
      }
    }
      
    }
    

    

  
  t_mat <- t_array_list[,,1]
  colnames(t_mat) <- lon[lon_grid]
  rownames(t_mat) <- lat[lat_grid]
  
  mean_table_t <- melt(t_mat)
  
  mean_table_t$Cruise <- unique_cruise_dates[1]
  
  for (i in 2:length(unique_cruise_dates)) {
    
    t_mat <- t_array_list[,,i]
    colnames(t_mat) <- lon[lon_grid]
    rownames(t_mat) <- lat[lat_grid]
    
    mean_table_t_i <- melt(t_mat)
    mean_table_t_i$Cruise <- unique_cruise_dates[i]
    
    mean_table_t <- bind_rows(mean_table_t, mean_table_t_i)

  }
  
  colnames(mean_table_t) <- c("lat", "lon", "Coeff", "Cruise_Date")

  mean_table_t$lon <- mean_table_t$lon - 360
  
  surf_dat <- full_dat[which(full_dat$CC_Depth < 15),]
  surf_dat <- surf_dat[-which(as.numeric(substr(surf_dat$Sta_ID,2,5)) < 76.7),]
  
  map <- map_data("world") 
  surf_dat$som_id <- as.character(surf_dat$som_id)
  
  for (i in 1:length(unique_cruise_dates)) {
    
    subset_surf <- surf_dat[which(surf_dat$Cruise_Date == unique_cruise_dates[i]),]
    subset_t <- mean_table_t[which(mean_table_t$Cruise_Date == unique_cruise_dates[i]),]
    
    t_som <- ggplot() + 
      geom_tile(data = subset_t, aes(x = lon, y = lat, fill = Coeff), width =0.26, height = 0.26) +
      scale_fill_gradient2(name = "Coeff. Var SST", low = "darkblue", mid = "white", high = "darkred",
                           limits = c(0.015,0.15), oob = squish, midpoint = 0.07) +
      geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
      coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
      geom_point(data = subset_surf, aes_string(x = "Lon_Dec", y = "Lat_Dec", color = "som_id"),
                 size = 4, stroke = 0.1, shape = 13) +
      scale_color_manual(name = "SOM Cluster", values = c("red","blue")) +
      theme(legend.position = "right",
            legend.background = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = "black", linetype = "solid", size = 1),
            title = element_text(hjust = 0.5), axis.line = element_blank()) +
      labs(x = "Longitude", y = "Latitude") + ggtitle(subset_surf$Cruise_Date)
    
    if(i < 10){
      png(filename = paste0("figures/animations/",type,"_0",i,".png"), width = 6, height = 6, units = "in", res = 72)
      print(t_som)
      dev.off()
      }
    if(i > 9){
      png(filename = paste0("figures/animations/",type,"_",i,".png"), width = 6, height = 6, units = "in", res = 72)
      print(t_som)
      dev.off()
      }
    
  
    
    
  }
  
  list.files(path = "figures/animations/", pattern = type, full.names = TRUE) %>%
    map(image_read) %>%
    image_join() %>%
    image_animate(fps = 0.5) %>%
    image_write(figure_name)
  
  
  
  
}

som_stacked_figure <- function(in_plastid = "output/plastid_16s_map.Rdata",
                               in_cyano = "output/cyano_16s_map.Rdata",
                               in_bacteria = "output/bacteria_m_euks_16s_map.Rdata",
                               in_euks = "output/euks_hetero_18sv9_map.Rdata",
                               in_phyto = "output/euks_auto_18sv9_map.Rdata",
                               figure_name = paste0("figures/som_summary_",Sys.Date(),".pdf")){
  
  
  map <- map_data("world")   
  
  # eukaryotic phytoplaknton
  
  load(in_phyto)
  
  som_color_vector = c("darkred", "darkblue")
  som_name_vector = c("Nearshore", "Offshore")
  centroid_color = c("blue","red")
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  p1 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",1)), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = som_color_vector[1], limits = c(0,1)) +
    ggtitle(paste0("B. ",som_name_vector[1]," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  centroid1 <- wt.centroid(x =centroid_df , p = 2)
  p1 <- p1 + geom_point(aes(x = centroid1@coords[1], y = centroid1@coords[2]), color = centroid_color[1], size = 5, pch = 10)
  
  p2 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",2)), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = som_color_vector[2], limits = c(0,1)) +
    ggtitle(paste0("A. ",som_name_vector[2]," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  centroid2 <- wt.centroid(x =centroid_df , p = 3)
  p2 <- p2 + geom_point(aes(x = centroid2@coords[1], y = centroid2@coords[2]), color = centroid_color[2], size = 5, pch = 10)
  
  title <- ggdraw() + draw_label("Eukaryotic\n Phytoplankton", fontface='bold')
  phyto <- plot_grid(title, plot_grid(p2,p1, ncol = 2), nrow = 1, rel_widths = c(0.2,1))
  
  # heterotrophic eukaryotes
  
  load(in_euks)
  
  som_color_vector = c("darkblue", "darkred")
  som_name_vector = c("Offshore", "Nearshore")
  centroid_color = c("red","blue")
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  p1 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",1)), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = som_color_vector[1], limits = c(0,1)) +
    ggtitle(paste0("C. ",som_name_vector[1]," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  centroid1 <- wt.centroid(x =centroid_df , p = 2)
  p1 <- p1 + geom_point(aes(x = centroid1@coords[1], y = centroid1@coords[2]), color = centroid_color[1], size = 5, pch = 10)
  
  p2 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",2)), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = som_color_vector[2], limits = c(0,1)) +
    ggtitle(paste0("D. ",som_name_vector[2]," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  centroid2 <- wt.centroid(x =centroid_df , p = 3)
  p2 <- p2 + geom_point(aes(x = centroid2@coords[1], y = centroid2@coords[2]), color = centroid_color[2], size = 5, pch = 10)
  
  title <- ggdraw() + draw_label("Heterotrophic\n Eukaryotes", fontface='bold')
  euks <- plot_grid(title, plot_grid(p1,p2, ncol = 2), nrow = 1, rel_widths = c(0.2,1))
  
  # cyanobacteria
  
  load(in_cyano)
  
  som_color_vector = c("darkblue", "darkred")
  som_name_vector = c("Offshore", "Nearshore")
  centroid_color = c("red","blue")
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  p1 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",1)), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = som_color_vector[1], limits = c(0,1)) +
    ggtitle(paste0("E. ",som_name_vector[1]," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  centroid1 <- wt.centroid(x =centroid_df , p = 2)
  p1 <- p1 + geom_point(aes(x = centroid1@coords[1], y = centroid1@coords[2]), color = centroid_color[1], size = 5, pch = 10)
  
  p2 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",2)), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = som_color_vector[2], limits = c(0,1)) +
    ggtitle(paste0("F. ",som_name_vector[2]," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  centroid2 <- wt.centroid(x =centroid_df , p = 3)
  p2 <- p2 + geom_point(aes(x = centroid2@coords[1], y = centroid2@coords[2]), color = centroid_color[2], size = 5, pch = 10)
  
  title <- ggdraw() + draw_label("Cyanobacteria", fontface='bold')
  cyano <- plot_grid(title, plot_grid(p1,p2, ncol = 2), nrow = 1, rel_widths = c(0.2,1))
  
  # heterotrophic bacteria archaea
  
  load(in_bacteria)
  
  som_color_vector = c("darkblue", "darkred")
  som_name_vector = c("Offshore", "Nearshore")
  centroid_color = c("red","blue")
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  p1 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",1)), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = som_color_vector[1], limits = c(0,1)) +
    ggtitle(paste0("H. ",som_name_vector[1]," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  centroid1 <- wt.centroid(x =centroid_df , p = 2)
  p1 <- p1 + geom_point(aes(x = centroid1@coords[1], y = centroid1@coords[2]), color = centroid_color[1], size = 5, pch = 10)
  
  p2 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",2)), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = som_color_vector[2], limits = c(0,1)) +
    ggtitle(paste0("I. ",som_name_vector[2]," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  centroid2 <- wt.centroid(x =centroid_df , p = 3)
  p2 <- p2 + geom_point(aes(x = centroid2@coords[1], y = centroid2@coords[2]), color = centroid_color[2], size = 5, pch = 10)
  
  title <- ggdraw() + draw_label("Heterotrophic\n Bacteria/Archaea", fontface='bold')
  bact <- plot_grid(title, plot_grid(p1,p2, ncol = 2), nrow = 1, rel_widths = c(0.2,1))
  
  
  all <- plot_grid(phyto, euks, cyano, bact, nrow = 4)
  
  pdf(figure_name, width = 10, height = 12)
  print(all)
  dev.off()
  
}

lm_stacked_figure <- function(in_plastid = "output/plastid_16s_glm.Rdata",
                               in_cyano = "output/cyano_16s_glm.Rdata",
                               in_bacteria = "output/bacteria_m_euks_16s_glm.Rdata",
                               in_euks = "output/euks_18sv9_glm.Rdata",
                               in_phyto = "output/euks_auto_18sv9_glm.Rdata",
                               figure_name = paste0("figures/glm_summary_",Sys.Date(),".pdf")){
  
  
  # Eukaryotic Phytoplankton
  
  load(in_phyto)
  
  summary_fit <- summary(lm(som_plots$freq[which(som_plots$cluster == "som_1")]~som_plots$temp_mean[which(som_plots$cluster == "som_1")]))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_mean <- summary_fit$r.squared
  
  summary_fit <- summary(lm(som_plots$freq[which(som_plots$cluster == "som_1")]~som_plots$temp_coeff[which(som_plots$cluster == "som_1")]))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_coeff <- summary_fit$r.squared
  
  
  mean_plot <- ggplot(som_plots, aes(x = temp_mean, y = freq, color = cluster)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", linetype = "solid"),
          panel.background = element_blank(),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    stat_smooth(method = "lm", level = 0.95) +
    scale_color_manual(values = c("red", "blue")) +
    xlab("Mean SST (°C)") + ylab("Frequency") + 
    annotate("text", x = 15.05, y = 0.125, label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 3) + ggtitle("A. SOM Freq ~ Mean SST")
  
  coeff_plot <- ggplot(som_plots, aes(x = temp_coeff, y = freq, color = cluster)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", linetype = "solid"),
          panel.background = element_blank(),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    stat_smooth(method = "lm", level = 0.95) +
    scale_color_manual(values = c("red", "blue")) +
    xlab("SST Coeff. Var") + ylab("Frequency") + 
    annotate("text", x = 0.1, y = 0.5, label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 3) + ggtitle("B. SOM Freq ~ Coeff. Var. SST")
  
  title <- ggdraw() + draw_label("Eukaryotic\n Phytoplankton", fontface='bold')
  phyto_som_plots <- plot_grid(title, mean_plot, coeff_plot, ncol = 3,rel_widths = c(0.4,1,1))
  
# heterotrophic euks
  
  load(in_euks)
  
  summary_fit <- summary(lm(som_plots$freq[which(som_plots$cluster == "som_1")]~som_plots$temp_mean[which(som_plots$cluster == "som_1")]))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_mean <- summary_fit$r.squared
  
  summary_fit <- summary(lm(som_plots$freq[which(som_plots$cluster == "som_1")]~som_plots$temp_coeff[which(som_plots$cluster == "som_1")]))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_coeff <- summary_fit$r.squared
  
  
  mean_plot <- ggplot(som_plots, aes(x = temp_mean, y = freq, color = cluster)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", linetype = "solid"),
          panel.background = element_blank(),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    stat_smooth(method = "lm", level = 0.95) +
    scale_color_manual(values = c("red", "blue")) +
    xlab("Mean SST (°C)") + ylab("Frequency") + 
    annotate("text", x = 16, y = 0.5, label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 3) + ggtitle("C. SOM Freq ~ Mean SST")
  
  coeff_plot <- ggplot(som_plots, aes(x = temp_coeff, y = freq, color = cluster)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", linetype = "solid"),
          panel.background = element_blank(),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    stat_smooth(method = "lm", level = 0.95) +
    scale_color_manual(values = c("red", "blue")) +
    xlab("SST Coeff. Var") + ylab("Frequency") + 
    annotate("text", x = 0.1, y = 0.5, label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 3) + ggtitle("D. SOM Freq ~ Coeff. Var. SST")
  
  title <- ggdraw() + draw_label("Heterotrophic\n Eukaryotes", fontface='bold')
  hetero_euk_som_plots <- plot_grid(title, mean_plot, coeff_plot, ncol = 3,rel_widths = c(0.4,1,1))
  
  # cyanobacteria
  
  load(in_cyano)
  
  summary_fit <- summary(lm(som_plots$freq[which(som_plots$cluster == "som_1")]~som_plots$temp_mean[which(som_plots$cluster == "som_1")]))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_mean <- summary_fit$r.squared
  
  summary_fit <- summary(lm(som_plots$freq[which(som_plots$cluster == "som_1")]~som_plots$temp_coeff[which(som_plots$cluster == "som_1")]))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_coeff <- summary_fit$r.squared
  
  
  mean_plot <- ggplot(som_plots, aes(x = temp_mean, y = freq, color = cluster)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", linetype = "solid"),
          panel.background = element_blank(),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    stat_smooth(method = "lm", level = 0.95) +
    scale_color_manual(values = c("red", "blue")) +
    xlab("Mean SST (°C)") + ylab("Frequency") + 
    annotate("text", x = 15.3, y = 0.88, label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 3) + ggtitle("E. SOM Freq ~ Mean SST")
  
  coeff_plot <- ggplot(som_plots, aes(x = temp_coeff, y = freq, color = cluster)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", linetype = "solid"),
          panel.background = element_blank(),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    stat_smooth(method = "lm", level = 0.95) +
    scale_color_manual(values = c("red", "blue")) +
    xlab("SST Coeff. Var") + ylab("Frequency") + 
    annotate("text", x = 0.095, y = 0.5, label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 3) + ggtitle("F. SOM Freq ~ Coeff. Var. SST")
  
  title <- ggdraw() + draw_label("Cyanobacteria", fontface='bold')
  cyano_som_plots <- plot_grid(title, mean_plot, coeff_plot, ncol = 3,rel_widths = c(0.4,1,1))
  
  # bacteria archaea
  
  load(in_bacteria)
  
  summary_fit <- summary(lm(som_plots$freq[which(som_plots$cluster == "som_1")]~som_plots$temp_mean[which(som_plots$cluster == "som_1")]))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_mean <- summary_fit$r.squared
  
  summary_fit <- summary(lm(som_plots$freq[which(som_plots$cluster == "som_1")]~som_plots$temp_coeff[which(som_plots$cluster == "som_1")]))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_coeff <- summary_fit$r.squared
  
  
  mean_plot <- ggplot(som_plots, aes(x = temp_mean, y = freq, color = cluster)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", linetype = "solid"),
          panel.background = element_blank(),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    stat_smooth(method = "lm", level = 0.95) +
    scale_color_manual(values = c("red", "blue")) +
    xlab("Mean SST (°C)") + ylab("Frequency") + 
    annotate("text", x = 15.3, y = 0.5, label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 3) + ggtitle("H. SOM Freq ~ Mean SST")
  
  coeff_plot <- ggplot(som_plots, aes(x = temp_coeff, y = freq, color = cluster)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", linetype = "solid"),
          panel.background = element_blank(),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    stat_smooth(method = "lm", level = 0.95) +
    scale_color_manual(values = c("red", "blue")) +
    xlab("SST Coeff. Var") + ylab("Frequency") + 
    annotate("text", x = 0.125, y = 0.5, label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 3) + ggtitle("I. SOM Freq ~ Coeff. Var. SST")
  
  title <- ggdraw() + draw_label("Heterotrophic\n Bacteria/Archaea", fontface='bold')
  bact_som_plots <- plot_grid(title, mean_plot, coeff_plot, ncol = 3,rel_widths = c(0.4,1,1))
  
  p <- ggplot(som_plots, aes(x = temp_coeff, y = freq, color = cluster)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", linetype = "solid"),
          panel.background = element_blank(),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          plot.title = element_text(hjust = 0.5)) +
    stat_smooth(method = "lm", level = 0.95) +
    scale_color_manual(labels = c("Nearshore", "Offshore"), values = c("red", "blue")) +
    xlab("SST Coeff. Var") + ylab("Frequency") + 
    annotate("text", x = 0.125, y = 0.5, label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 3) + ggtitle("I. SOM Freq ~ Coeff. Var. SST")
  
  legend_plot <- get_legend(p)
  
  full_plot <- plot_grid(plot_grid(phyto_som_plots, hetero_euk_som_plots,
                      cyano_som_plots, bact_som_plots, ncol = 1),
            legend_plot, rel_widths = c(1,0.2))
  
  pdf(file = figure_name, width = 12, height = 12)
  print(full_plot)
  dev.off()
  
  }

###### Running Figures #####

######## som figures #######

som_figure(SST_file = "output/CALCOFI_temp_tables.Rdata",
           map_file = "output/plastid_16s_map.Rdata",
           reg_file = "output/plastid_16s_glm.Rdata",
           figure_name = paste0("figures/plastid_16s_som_",Sys.Date(),".pdf"),
           som_color_vector = c("darkblue", "darkred"),
           som_name_vector = c("Offshore", "Nearshore"),
           centroid_color = c("red","blue"), switch = FALSE,
           main = "16s Plastid", mean_position = c(15,1),
           coeff_position = c(0.123,1))

som_figure(SST_file = "output/CALCOFI_temp_tables.Rdata",
           map_file = "output/bacteria_m_euks_16s_map.Rdata",
           reg_file = "output/bacteria_m_euks_16s_glm.Rdata",
           figure_name = paste0("figures/bacteria_m_euks_16s_som_",Sys.Date(),".pdf"),
           som_color_vector = c("darkblue", "darkred"),
           som_name_vector = c("Offshore", "Nearshore"),
           centroid_color = c("red","blue"), switch = FALSE,
           main = "16s Bacteria (no Euks)", mean_position = c(15,0.5),
           coeff_position = c(0.123,0.5))

som_figure(SST_file = "output/CALCOFI_temp_tables.Rdata",
           map_file = "output/cyano_16s_map.Rdata",
           reg_file = "output/cyano_16s_glm.Rdata",
           figure_name = paste0("figures/cyano_16s_som_",Sys.Date(),".pdf"),
           som_color_vector = c("darkblue", "darkred"),
           som_name_vector = c("Offshore", "Nearshore"),
           centroid_color = c("red","blue"), switch = FALSE,
           main = "16s Cyanobacteria", mean_position = c(15,1.2),
           coeff_position = c(0.123,1.2))

som_figure(SST_file = "output/CALCOFI_temp_tables.Rdata",
           map_file = "output/euks_auto_18sv9_map.Rdata",
           reg_file = "output/euks_auto_18sv9_glm.Rdata",
           figure_name = paste0("figures/euks_auto_18sv9_som_",Sys.Date(),".pdf"),
           som_color_vector = c("darkred", "darkblue"),
           som_name_vector = c("Nearshore", "Offshore"),
           centroid_color = c("blue","red"), switch = TRUE,
           main = "18sv9 Autotrophic Eukaryota", mean_position = c(15,1.2),
           coeff_position = c(0.123,1.2))

som_figure(SST_file = "output/CALCOFI_temp_tables.Rdata",
           map_file = "output/euks_hetero_18sv9_map.Rdata",
           reg_file = "output/euks_hetero_18sv9_glm.Rdata",
           figure_name = paste0("figures/euks_hetero_18sv9_som_",Sys.Date(),".pdf"),
           som_color_vector = c("darkblue", "darkred"),
           som_name_vector = c("Offshore", "Nearshore"),
           centroid_color = c("red","blue"), switch = FALSE,
           main = "18sv9 Heterotrophic Eukaryota", mean_position = c(15,1.2),
           coeff_position = c(0.123,1.2))

# Surface

som_figure(SST_file = "output/CALCOFI_temp_tables.Rdata",
           map_file = "output/plastid_16s_surf_map.Rdata",
           reg_file = "output/plastid_16s_surf_glm.Rdata",
           figure_name = paste0("figures/plastid_16s_surf_som_",Sys.Date(),".pdf"),
           som_color_vector = c("darkred", "darkblue"),
           som_name_vector = c("Nearshore", "Offshore"),
           centroid_color = c("blue","red"), switch = TRUE,
           main = "16s Plastid Surface", mean_position = c(16.2,1),
           coeff_position = c(0.123,1))

som_figure(SST_file = "output/CALCOFI_temp_tables.Rdata",
           map_file = "output/bacteria_m_euks_16s_surf_map.Rdata",
           reg_file = "output/bacteria_m_euks_16s_surf_glm.Rdata",
           figure_name = paste0("figures/bacteria_m_euks_16s_surf_som_",Sys.Date(),".pdf"),
           som_color_vector = c("darkred", "darkblue"),
           som_name_vector = c("Nearshore", "Offshore"),
           centroid_color = c("blue","red"), switch = TRUE,
           main = "16s Bacteria (no Euks) Surface", mean_position = c(16,1),
           coeff_position = c(0.1,1))

som_figure(SST_file = "output/CALCOFI_temp_tables.Rdata",
           map_file = "output/cyano_16s_surf_map.Rdata",
           reg_file = "output/cyano_16s_surf_glm.Rdata",
           figure_name = paste0("figures/cyano_16s_surf_som_",Sys.Date(),".pdf"),
           som_color_vector = c("darkblue", "darkred"),
           som_name_vector = c("Offshore", "Nearshore"),
           centroid_color = c("red","blue"), switch = FALSE,
           main = "16s Cyanobacteria Surface", mean_position = c(16.2,1.2),
           coeff_position = c(0.123,1.2))

som_figure(SST_file = "output/CALCOFI_temp_tables.Rdata",
           map_file = "output/euks_auto_18sv9_surf_map.Rdata",
           reg_file = "output/euks_auto_18sv9_surf_glm.Rdata",
           figure_name = paste0("figures/euks_auto_18sv9_surf_som_",Sys.Date(),".pdf"),
           som_color_vector = c("darkred", "darkblue"),
           som_name_vector = c("Nearshore", "Offshore"),
           centroid_color = c("blue","red"), switch = FALSE,
           main = "18sv9 Autotrophic Eukaryota Surface", mean_position = c(16,1),
           coeff_position = c(0.123,1.2))

som_figure(SST_file = "output/CALCOFI_temp_tables.Rdata",
           map_file = "output/euks_hetero_18sv9_surf_map.Rdata",
           reg_file = "output/euks_hetero_18sv9_surf_glm.Rdata",
           figure_name = paste0("figures/euks_hetero_18sv9_surf_som_",Sys.Date(),".pdf"),
           som_color_vector = c("darkblue", "darkred"),
           som_name_vector = c("Offshore", "Nearshore"),
           centroid_color = c("red","blue"), switch = FALSE,
           main = "18sv9 Heterotrophic Eukaryota Surface", mean_position = c(15.5,1),
           coeff_position = c(0.123,1))

# Depth

som_figure(SST_file = "output/CALCOFI_temp_tables.Rdata",
           map_file = "output/plastid_16s_deep_map.Rdata",
           reg_file = "output/plastid_16s_deep_glm.Rdata",
           figure_name = paste0("figures/plastid_16s_deep_som_",Sys.Date(),".pdf"),
           som_color_vector = c("darkred", "darkblue"),
           som_name_vector = c("Nearshore", "Offshore"),
           centroid_color = c("blue","red"), switch = TRUE,
           main = "16s Plastid Depth", mean_position = c(16.2,1),
           coeff_position = c(0.123,1))

som_figure(SST_file = "output/CALCOFI_temp_tables.Rdata",
           map_file = "output/bacteria_m_euks_16s_deep_map.Rdata",
           reg_file = "output/bacteria_m_euks_16s_deep_glm.Rdata",
           figure_name = paste0("figures/bacteria_m_euks_16s_deep_som_",Sys.Date(),".pdf"),
           som_color_vector = c("darkred", "darkblue"),
           som_name_vector = c("Nearshore", "Offshore"),
           centroid_color = c("blue","red"), switch = TRUE,
           main = "16s Bacteria (no Euks) Deptj", mean_position = c(16,1),
           coeff_position = c(0.1,1))

som_figure(SST_file = "output/CALCOFI_temp_tables.Rdata",
           map_file = "output/cyano_16s_deep_map.Rdata",
           reg_file = "output/cyano_16s_deep_glm.Rdata",
           figure_name = paste0("figures/cyano_16s_deep_som_",Sys.Date(),".pdf"),
           som_color_vector = c("darkblue", "darkred"),
           som_name_vector = c("Offshore", "Nearshore"),
           centroid_color = c("red","blue"), switch = FALSE,
           main = "16s Cyanobacteria Depth", mean_position = c(16.2,1.2),
           coeff_position = c(0.123,1.2))

som_figure(SST_file = "output/CALCOFI_temp_tables.Rdata",
           map_file = "output/euks_auto_18sv9_deep_map.Rdata",
           reg_file = "output/euks_auto_18sv9_deep_glm.Rdata",
           figure_name = paste0("figures/euks_auto_18sv9_deep_som_",Sys.Date(),".pdf"),
           som_color_vector = c("darkred", "darkblue"),
           som_name_vector = c("Nearshore", "Offshore"),
           centroid_color = c("blue","red"), switch = FALSE,
           main = "18sv9 Autotrophic Eukaryota Depth", mean_position = c(16,1),
           coeff_position = c(0.123,1.2))

som_figure(SST_file = "output/CALCOFI_temp_tables.Rdata",
           map_file = "output/euks_hetero_18sv9_deep_map.Rdata",
           reg_file = "output/euks_hetero_18sv9_deep_glm.Rdata",
           figure_name = paste0("figures/euks_hetero_18sv9_deep_som_",Sys.Date(),".pdf"),
           som_color_vector = c("darkblue", "darkred"),
           som_name_vector = c("Offshore", "Nearshore"),
           centroid_color = c("red","blue"), switch = FALSE,
           main = "18sv9 Heterotrophic Eukaryota Depth", mean_position = c(15.5,1),
           coeff_position = c(0.123,1))

###### aic figures #########

aic_table(map_file = "output/bacteria_16s_map.Rdata",
          figure_name = paste0("figures/bacteria_16s_aic_",Sys.Date(),".pdf"),
          minimum_tp = 8)

aic_table(map_file = "output/plastid_16s_map.Rdata",
          figure_name = paste0("figures/plastid_16s_aic_",Sys.Date(),".pdf"),
          minimum_tp = 8)

aic_table(map_file = "output/bacteria_m_euks_16s_map.Rdata",
          figure_name = paste0("figures/bacteria_m_euks_16s_aic_",Sys.Date(),".pdf"),
          minimum_tp = 8)

aic_table(map_file = "output/cyano_16s_map.Rdata",
          figure_name = paste0("figures/cyano_16s_aic_",Sys.Date(),".pdf"),
          minimum_tp = 8)

aic_table(map_file = "output/euks_18sv9_map.Rdata",
          figure_name = paste0("figures/euks_18sv9_aic_",Sys.Date(),".pdf"),
          minimum_tp = 8)

aic_table(map_file = "output/euks_auto_18sv9_map.Rdata",
          figure_name = paste0("figures/euks_auto_18sv9_aic_",Sys.Date(),".pdf"),
          minimum_tp = 8)

aic_table(map_file = "output/euks_hetero_18sv9_map.Rdata",
          figure_name = paste0("figures/euks_hetero_18sv9_aic_",Sys.Date(),".pdf"),
          minimum_tp = 8)

###### random forest figures #########

random_forest_figure(full_data_file = "output/bacteria_16s_full_data.Rdata",
                     figure_name = paste0("figures/bacteria_16s_rf_",Sys.Date(),".pdf"),
                    main = "16s Bacteria")

random_forest_figure(full_data_file = "output/plastid_16s_full_data.Rdata",
                     figure_name = paste0("figures/plastid_16s_rf_",Sys.Date(),".pdf"),
                     main = "16s Plastid")

random_forest_figure(full_data_file = "output/bacteria_m_euks_16s_full_data.Rdata",
                     figure_name = paste0("figures/bacteria_m_euks_16s_rf_",Sys.Date(),".pdf"),
                     main = "16s Bacteria (No Euks)")

random_forest_figure(full_data_file = "output/cyano_16s_full_data.Rdata",
                     figure_name = paste0("figures/cyano_16s_rf_",Sys.Date(),".pdf"),
                     main = "16s Cyanobacteria")

random_forest_figure(full_data_file = "output/euks_18sv9_full_data.Rdata",
                     figure_name = paste0("figures/euks_18sv9_rf_",Sys.Date(),".pdf"),
                     main = "18sv9 Eukaryota")

random_forest_figure(full_data_file = "output/euks_auto_18sv9_full_data.Rdata",
                     figure_name = paste0("figures/euks_auto_18sv9_rf_",Sys.Date(),".pdf"),
                     main = "18sv9 Auto Eukaryota")

random_forest_figure(full_data_file = "output/euks_hetero_18sv9_full_data.Rdata",
                     figure_name = paste0("figures/euks_hetero_18sv9_rf_",Sys.Date(),".pdf"),
                     main = "18sv9 Hetero Eukaryota")

###### diversity figures ######

diversity_figure(map_file = "output/bacteria_16s_map.Rdata",
                 full_data_file = "output/bacteria_16s_full_data.Rdata",
                 reg_file = "output/bacteria_16s_glm.Rdata",
                 figure_name = paste0("figures/bacteria_16s_diversity_",Sys.Date(),".pdf"),
                 minimum_tp = 8,
                 div_mean_position = c(15,4.73), div_coeff_position = c(0.125,4.78),
                 even_mean_position = c(15.2,0.73), even_coeff_position = c(0.128,0.73),
                 rich_mean_position = c(14.7,500), rich_coeff_position = c(0.123,550),
                 shannon_limits = c(3.8,5.4), even_limits = c(0.6,0.8), rich_limits = c(300,500),
                 main = "16s Bacteria")

diversity_figure(map_file = "output/plastid_16s_map.Rdata",
                 full_data_file = "output/plastid_16s_full_data.Rdata",
                 reg_file = "output/plastid_16s_glm.Rdata",
                 figure_name = paste0("figures/plastid_16s_diversity_",Sys.Date(),".pdf"),
                 minimum_tp = 8,
                 div_mean_position = c(15,3.8), div_coeff_position = c(0.121,3.8),
                 even_mean_position = c(15,0.77), even_coeff_position = c(0.128,0.77),
                 rich_mean_position = c(17,100), rich_coeff_position = c(0.122,140),
                 shannon_limits = c(3.3,3.8), even_limits = c(0.71,0.77), rich_limits = c(80,180),
                 main = "16s Plastid")

diversity_figure(map_file = "output/bacteria_m_euks_16s_map.Rdata",
                 full_data_file = "output/bacteria_m_euks_16s_full_data.Rdata",
                 reg_file = "output/bacteria_m_euks_16s_glm.Rdata",
                 figure_name = paste0("figures/bacteria_m_euks_16s_diversity_",Sys.Date(),".pdf"),
                 minimum_tp = 8,
                 div_mean_position = c(15,4.73), div_coeff_position = c(0.125,4.78),
                 even_mean_position = c(15.2,0.73), even_coeff_position = c(0.128,0.73),
                 rich_mean_position = c(14.7,500), rich_coeff_position = c(0.123,550),
                 shannon_limits = c(3.8,5.4), even_limits = c(0.6,0.8), rich_limits = c(300,500),
                 main = "16s Bacteria (No Euks)")

diversity_figure(map_file = "output/cyano_16s_map.Rdata",
                 full_data_file = "output/cyano_16s_full_data.Rdata",
                 reg_file = "output/cyano_16s_glm.Rdata",
                 figure_name = paste0("figures/cyano_16s_diversity_",Sys.Date(),".pdf"),
                 minimum_tp = 8,
                 div_mean_position = c(15,3.5), div_coeff_position = c(0.121,3.5),
                 even_mean_position = c(15,0.5), even_coeff_position = c(0.128,0.5),
                 rich_mean_position = c(15,100), rich_coeff_position = c(0.122,100),
                 shannon_limits = c(2.1,3.6), even_limits = c(0.50,0.77), rich_limits = c(55,100),
                 main = "16s Cyanobacteria")

diversity_figure(map_file = "output/euks_18sv9_map.Rdata",
                 full_data_file = "output/euks_18sv9_full_data.Rdata",
                 reg_file = "output/euks_18sv9_glm.Rdata",
                 figure_name = paste0("figures/euks_18sv9_diversity_",Sys.Date(),".pdf"),
                 minimum_tp = 8,
                 div_mean_position = c(17,3.45), div_coeff_position = c(0.125,5),
                 even_mean_position = c(17,0.6), even_coeff_position = c(0.128,0.8),
                 rich_mean_position = c(15,700), rich_coeff_position = c(0.123,750),
                 shannon_limits = c(3.8,5.4), even_limits = c(0.6,0.8), rich_limits = c(300,650),
                 main = "18sv9 Eukaryota")

diversity_figure(map_file = "output/euks_auto_18sv9_map.Rdata",
                 full_data_file = "output/euks_auto_18sv9_full_data.Rdata",
                 reg_file = "output/euks_auto_18sv9_glm.Rdata",
                 figure_name = paste0("figures/euks_auto_18sv9_diversity_",Sys.Date(),".pdf"),
                 minimum_tp = 8,
                 div_mean_position = c(17,3.3), div_coeff_position = c(0.1,3.4),
                 even_mean_position = c(17,0.7), even_coeff_position = c(0.128,0.7),
                 rich_mean_position = c(15,180), rich_coeff_position = c(0.123,180),
                 shannon_limits = c(3.5,4.2), even_limits = c(0.7,0.82), rich_limits = c(100,200),
                 main = "18sv9 Autotrophic Eukaryota")

diversity_figure(map_file = "output/euks_hetero_18sv9_map.Rdata",
                 full_data_file = "output/euks_hetero_18sv9_full_data.Rdata",
                 reg_file = "output/euks_hetero_18sv9_glm.Rdata",
                 figure_name = paste0("figures/euks_hetero_18sv9_diversity_",Sys.Date(),".pdf"),
                 minimum_tp = 8,
                 div_mean_position = c(16,3.2), div_coeff_position = c(0.1,3.35),
                 even_mean_position = c(16,0.55), even_coeff_position = c(0.128,0.73),
                 rich_mean_position = c(15,500), rich_coeff_position = c(0.123,600),
                 shannon_limits = c(3.2,4.5), even_limits = c(0.5,0.82), rich_limits = c(220,600),
                 main = "18sv9 Heterotrophic Eukaryota")


######## alpha versus gamma figures #######

alpha_versus_gamma_figure(full_data_file = "output/bacteria_16s_full_data.Rdata",
                          raw_data_file = "data/16s_bacteria.Rdata",
                          map_file = "output/bacteria_16s_map.Rdata", minimum_tp = 8,
                          figure_name = paste0("figures/bacteria_16s_alpha_gamma_",Sys.Date(),".pdf"),
                          main = "16s Bacteria", gamma_position = c(0.128,5.6), alpha_position = c(0.10,4.3))

alpha_versus_gamma_figure(full_data_file = "output/plastid_16s_full_data.Rdata",
                          raw_data_file = "data/16s_plastids.Rdata",
                          map_file = "output/plastid_16s_map.Rdata", minimum_tp = 8,
                          figure_name = paste0("figures/plastid_16s_alpha_gamma_",Sys.Date(),".pdf"),
                          main = "16s Plastid", gamma_position = c(0.128,4.4), alpha_position = c(0.10,3.2))

alpha_versus_gamma_figure(full_data_file = "output/bacteria_m_euks_16s_full_data.Rdata",
                          raw_data_file = "data/16s_bacteria_m_euks.Rdata",
                          map_file = "output/bacteria_m_euks_16s_map.Rdata", minimum_tp = 8,
                          figure_name = paste0("figures/bacteria_m_euks_16s_alpha_gamma_",Sys.Date(),".pdf"),
                          main = "16s Bacteria (No Euks)", gamma_position = c(0.128,5.6), alpha_position = c(0.10,4.3))

alpha_versus_gamma_figure(full_data_file = "output/cyano_16s_full_data.Rdata",
                          raw_data_file = "data/16s_cyanos.Rdata",
                          map_file = "output/cyano_16s_map.Rdata", minimum_tp = 8,
                          figure_name = paste0("figures/cyano_16s_alpha_gamma_",Sys.Date(),".pdf"),
                          main = "16s Cyanobacteria", gamma_position = c(0.128,4.4), alpha_position = c(0.10,2))

alpha_versus_gamma_figure(full_data_file = "output/euks_18sv9_full_data.Rdata",
                          raw_data_file = "data/18s_euks.Rdata",
                          map_file = "output/euks_18sv9_map.Rdata", minimum_tp = 8,
                          figure_name = paste0("figures/euks_18sv9_alpha_gamma_",Sys.Date(),".pdf"),
                          main = "18sv9 Eukaryota", gamma_position = c(0.128,6.2), alpha_position = c(0.10,3.9))

alpha_versus_gamma_figure(full_data_file = "output/euks_auto_18sv9_full_data.Rdata",
                          raw_data_file = "data/18s_autotrophic_euks.Rdata",
                          map_file = "output/euks_auto_18sv9_map.Rdata", minimum_tp = 8,
                          figure_name = paste0("figures/euks_auto_18sv9_alpha_gamma_",Sys.Date(),".pdf"),
                          main = "18sv9 Autotrophic Eukaryota", gamma_position = c(0.128,5.1), alpha_position = c(0.10,3.5))

alpha_versus_gamma_figure(full_data_file = "output/euks_hetero_18sv9_full_data.Rdata",
                          raw_data_file = "data/18s_heterotrophic_euks.Rdata",
                          map_file = "output/euks_hetero_18sv9_map.Rdata", minimum_tp = 8,
                          figure_name = paste0("figures/euks_hetero_18sv9_alpha_gamma_",Sys.Date(),".pdf"),
                          main = "18sv9 Heterotrophic Eukaryota", gamma_position = c(0.128,5.8), alpha_position = c(0.10,3.5))

######## Diversity over Time ########

diversity_over_time_figure(full_data_file = "output/bacteria_16s_full_data.Rdata",
                           figure_name = paste0("figures/bacteria_16s_diversity_time_",Sys.Date(),".pdf"),
                           shannon_limits = c(3.8,5.4), even_limits = c(0.6,0.8), rich_limits = c(300,500),
                           main = "16s Bacteria")

diversity_over_time_figure(full_data_file = "output/plastid_16s_full_data.Rdata",
                           figure_name = paste0("figures/plastid_16s_diversity_time_",Sys.Date(),".pdf"),
                           shannon_limits = c(3.3,3.8), even_limits = c(0.71,0.77), rich_limits = c(80,180),
                           main = "16s Plastids")

diversity_over_time_figure(full_data_file = "output/bacteria_m_euks_16s_full_data.Rdata",
                           figure_name = paste0("figures/bacteria_m_euks_16s_diversity_time_",Sys.Date(),".pdf"),
                           shannon_limits = c(3.8,5.4), even_limits = c(0.6,0.8), rich_limits = c(300,500),
                           main = "16s Bacteria (No Euks)")

diversity_over_time_figure(full_data_file = "output/cyano_16s_full_data.Rdata",
                           figure_name = paste0("figures/cyano_16s_diversity_time_",Sys.Date(),".pdf"),
                           shannon_limits = c(2.1,3.6), even_limits = c(0.50,0.77), rich_limits = c(55,100),
                           main = "16s Cyanobacteria")

diversity_over_time_figure(full_data_file = "output/euks_18sv9_full_data.Rdata",
                           figure_name = paste0("figures/euks_18sv9_diversity_time_",Sys.Date(),".pdf"),
                           shannon_limits = c(3.8,5.4), even_limits = c(0.6,0.8), rich_limits = c(300,650),
                           main = "18sv9 Eukaryota")

diversity_over_time_figure(full_data_file = "output/euks_auto_18sv9_full_data.Rdata",
                           figure_name = paste0("figures/euks_auto_18sv9_diversity_time_",Sys.Date(),".pdf"),
                           shannon_limits = c(3.5,4.2), even_limits = c(0.7,0.82), rich_limits = c(100,200),
                           main = "18sv9 Autotrophic Eukaryota")

diversity_over_time_figure(full_data_file = "output/euks_hetero_18sv9_full_data.Rdata",
                           figure_name = paste0("figures/euks_hetero_18sv9_diversity_time_",Sys.Date(),".pdf"),
                           shannon_limits = c(3.2,4.5), even_limits = c(0.5,0.82), rich_limits = c(220,600),
                           main = "18sv9 Heterotrophic Eukaryota")

####### Full AIC Table ######

full_aic_table_figure(in_plastid = "output/plastid_16s_map.Rdata", in_cyano = "output/cyano_16s_map.Rdata",
                       in_bacteria = "output/bacteria_m_euks_16s_map.Rdata", in_euks = "output/euks_hetero_18sv9_map.Rdata",
                       in_phyto = "output/euks_auto_18sv9_map.Rdata", minimum_tp = 8,
                       figure_name = paste0("figures/full_aic_table_",Sys.Date(),".pdf"),
                      figure_name_2 = paste0("figures/full_aic_plot_",Sys.Date(),".pdf"),
                      title_name = "Variable Importance All Samples")

# Surface

full_aic_table_figure(in_plastid = "output/plastid_16s_surf_map.Rdata", in_cyano = "output/cyano_16s_surf_map.Rdata",
                      in_bacteria = "output/bacteria_m_euks_16s_surf_map.Rdata", in_euks = "output/euks_hetero_18sv9_surf_map.Rdata",
                      in_phyto = "output/euks_auto_18sv9_surf_map.Rdata", minimum_tp = 4,
                      figure_name = paste0("figures/full_aic_table_surface_",Sys.Date(),".pdf"),
                      figure_name_2 = paste0("figures/full_aic_plot_surface_",Sys.Date(),".pdf"),
                      title_name = "Variable Importance Surface Samples")

# depth

full_aic_table_figure(in_plastid = "output/plastid_16s_deep_map.Rdata", in_cyano = "output/cyano_16s_deep_map.Rdata",
                      in_bacteria = "output/bacteria_m_euks_16s_deep_map.Rdata", in_euks = "output/euks_hetero_18sv9_deep_map.Rdata",
                      in_phyto = "output/euks_auto_18sv9_deep_map.Rdata", minimum_tp = 4,
                      figure_name = paste0("figures/full_aic_table_deep_",Sys.Date(),".pdf"),
                      figure_name_2 = paste0("figures/full_aic_plot_deep_",Sys.Date(),".pdf"),
                      title_name = "Variable Importance Deep Samples")


####### Animations ##########

# bacteria

som_current_animation(vel_table = "output/CALCOFI_Copernicus_GSC_Slice.Rdata",
                      som_data = "output/bacteria_16s_full_data.Rdata", type = "bacteria_16s_vel",
                      figure_name = paste0("figures/bacteria_16s_som_animation_vel_",Sys.Date(),".gif"))

som_temp_animation(temp_table = "output/CALCOFI_OI_SST_Slice.Rdata",
                   som_data = "output/bacteria_16s_full_data.Rdata", type = "bacteria_16s_temp",
                   figure_name = paste0("figures/bacteria_16s_som_animation_temp_",Sys.Date(),".gif"))

# plastid

som_current_animation(vel_table = "output/CALCOFI_Copernicus_GSC_Slice.Rdata",
                      som_data = "output/plastid_16s_full_data.Rdata", type = "plastid_16s_vel",
                      figure_name = paste0("figures/plastid_16s_som_animation_vel_",Sys.Date(),".gif"))

som_temp_animation(temp_table = "output/CALCOFI_OI_SST_Slice.Rdata",
                   som_data = "output/plastid_16s_full_data.Rdata", type = "plastid_16s_temp",
                   figure_name = paste0("figures/plastid_16s_som_animation_temp_",Sys.Date(),".gif"))

# cyano

som_current_animation(vel_table = "output/CALCOFI_Copernicus_GSC_Slice.Rdata",
                      som_data = "output/cyano_16s_full_data.Rdata", type = "cyano_16s_vel",
                      figure_name = paste0("figures/cyano_16s_som_animation_vel_",Sys.Date(),".gif"))

som_temp_animation(temp_table = "output/CALCOFI_OI_SST_Slice.Rdata",
                   som_data = "output/cyano_16s_full_data.Rdata", type = "cyano_16s_temp",
                   figure_name = paste0("figures/cyano_16s_som_animation_temp_",Sys.Date(),".gif"))

# euk phyto

som_current_animation(vel_table = "output/CALCOFI_Copernicus_GSC_Slice.Rdata",
                      som_data = "output/euks_auto_18sv9_full_data.Rdata", type = "phyto_18sv9_vel",
                      figure_name = paste0("figures/phyto_18sv9_som_animation_vel_",Sys.Date(),".gif"))

som_temp_animation(temp_table = "output/CALCOFI_OI_SST_Slice.Rdata",
                   som_data = "output/euks_auto_18sv9_full_data.Rdata", type = "phyto_18sv9_temp",
                   figure_name = paste0("figures/phyto_18sv9_som_animation_temp_",Sys.Date(),".gif"))

# euk hetero

som_current_animation(vel_table = "output/CALCOFI_Copernicus_GSC_Slice.Rdata",
                      som_data = "output/euks_18sv9_full_data.Rdata", type = "euks_18sv9_vel",
                      figure_name = paste0("figures/euks_18sv9_som_animation_vel_",Sys.Date(),".gif"))

som_temp_animation(temp_table = "output/CALCOFI_OI_SST_Slice.Rdata",
                   som_data = "output/euks_18sv9_full_data.Rdata", type = "euks_18sv9_temp",
                   figure_name = paste0("figures/euks_18sv9_som_animation_temp_",Sys.Date(),".gif"))

#### Stacked Figures #####

som_stacked_figure(in_plastid = "output/plastid_16s_map.Rdata",
                   in_cyano = "output/cyano_16s_map.Rdata",
                   in_bacteria = "output/bacteria_m_euks_16s_map.Rdata",
                   in_euks = "output/euks_hetero_18sv9_map.Rdata",
                   in_phyto = "output/euks_auto_18sv9_map.Rdata",
                   figure_name = paste0("figures/som_summary_",Sys.Date(),".pdf"))

som_stacked_figure(in_plastid = "output/plastid_16s_surf_map.Rdata",
                   in_cyano = "output/cyano_16s_surf_map.Rdata",
                   in_bacteria = "output/bacteria_m_euks_16s_surf_map.Rdata",
                   in_euks = "output/euks_hetero_18sv9_surf_map.Rdata",
                   in_phyto = "output/euks_auto_18sv9_surf_map.Rdata",
                   figure_name = paste0("figures/som_surface_summary_",Sys.Date(),".pdf"))

som_stacked_figure(in_plastid = "output/plastid_16s_deep_map.Rdata",
                   in_cyano = "output/cyano_16s_deep_map.Rdata",
                   in_bacteria = "output/bacteria_m_euks_16s_deep_map.Rdata",
                   in_euks = "output/euks_hetero_18sv9_deep_map.Rdata",
                   in_phyto = "output/euks_auto_18sv9_deep_map.Rdata",
                   figure_name = paste0("figures/som_depth_summary_",Sys.Date(),".pdf"))
