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

diversity_comparison <- function(in_file = "output/euks_auto_18sv9_full_data.Rdata",
                                   in_map = "output/euks_auto_18sv9_map.Rdata",
                                   out_diff_file = "output/euks_auto_18sv9_diffs_div.Rdata",
                                   title = "Eukaryotic Phytoplankton",
                                   upwelling_index = "output/upwelling_indicies.Rdata"){
  
  load(in_file)
  load(similar_mat)
  
  full_dat$Year <- as.numeric(substr(full_dat$Cruise,1,4))
  
  full_dat <- full_dat[-which(as.numeric(substr(full_dat$Sta_ID,2,3)) < 76),]
  
  full_dat$ML_NC <- full_dat$MLD_Sigma - full_dat$NCDepth
  
  early_som_maps <- full_dat %>% 
    filter(Year < 2017) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som_id == 1, na.rm = TRUE)/n(), som_2 = sum(som_id == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE),
              temp_mean = mean(T_degC, na.rm = TRUE), sal_mean = mean(Salnty, na.rm = TRUE),
              PO4_mean = mean(PO4ug, na.rm = TRUE), NO3_mean = mean(NO3ug, na.rm = TRUE),
              SiO3_mean = mean(SiO3ug, na.rm = TRUE), C14_mean = mean(IntC14, na.rm = TRUE),
              MLD_mean = mean(MLD_Sigma, na.rm = TRUE), NC_mean = mean(NCDepth, na.rm = TRUE),
              temp_coeff = sd(T_degC, na.rm = TRUE)/mean(T_degC, na.rm = TRUE),
              sal_coeff = sd(Salnty, na.rm = TRUE)/mean(Salnty, na.rm = TRUE),
              PO4_coeff =  sd(PO4ug, na.rm = TRUE)/mean(PO4ug, na.rm = TRUE),
              NO3_coeff = sd(NO3ug, na.rm = TRUE)/mean(NO3ug, na.rm = TRUE),
              SiO3_coeff = sd(SiO3ug, na.rm = TRUE)/mean(SiO3ug, na.rm = TRUE),
              C14_coeff = sd(IntC14, na.rm = TRUE)/mean(IntC14, na.rm = TRUE),
              MLD_coeff = sd(MLD_Sigma, na.rm = TRUE)/mean(MLD_Sigma, na.rm = TRUE),
              NC_coeff = sd(NCDepth, na.rm = TRUE)/mean(NCDepth, na.rm = TRUE),
              evenness = mean(evenness, na.rm = TRUE), shannon = mean(shannon, na.rm = TRUE),
              richness = mean(richness, na.rm = TRUE), 
              ML_NC_mean = mean(ML_NC, na.rm = TRUE))
  
  
  centroid_df <- SpatialPointsDataFrame(coords = early_som_maps[,c(6,5)], data = early_som_maps)
  centroid1 <- wt.centroid(x =centroid_df , p = 2)
  centroid2 <- wt.centroid(x =centroid_df , p = 3)
  
  early_dist <- distHaversine(centroid1@coords, centroid2@coords)/100
  
  late_som_maps <- full_dat %>% 
    filter(Year > 2016) %>%
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som_id == 1, na.rm = TRUE)/n(), som_2 = sum(som_id == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE),
              temp_mean = mean(T_degC, na.rm = TRUE), sal_mean = mean(Salnty, na.rm = TRUE),
              PO4_mean = mean(PO4ug, na.rm = TRUE), NO3_mean = mean(NO3ug, na.rm = TRUE),
              SiO3_mean = mean(SiO3ug, na.rm = TRUE), C14_mean = mean(IntC14, na.rm = TRUE),
              MLD_mean = mean(MLD_Sigma, na.rm = TRUE), NC_mean = mean(NCDepth, na.rm = TRUE),
              temp_coeff = sd(T_degC, na.rm = TRUE)/mean(T_degC, na.rm = TRUE),
              sal_coeff = sd(Salnty, na.rm = TRUE)/mean(Salnty, na.rm = TRUE),
              PO4_coeff =  sd(PO4ug, na.rm = TRUE)/mean(PO4ug, na.rm = TRUE),
              NO3_coeff = sd(NO3ug, na.rm = TRUE)/mean(NO3ug, na.rm = TRUE),
              SiO3_coeff = sd(SiO3ug, na.rm = TRUE)/mean(SiO3ug, na.rm = TRUE),
              C14_coeff = sd(IntC14, na.rm = TRUE)/mean(IntC14, na.rm = TRUE),
              MLD_coeff = sd(MLD_Sigma, na.rm = TRUE)/mean(MLD_Sigma, na.rm = TRUE),
              NC_coeff = sd(NCDepth, na.rm = TRUE)/mean(NCDepth, na.rm = TRUE),
              evenness = mean(evenness, na.rm = TRUE), shannon = mean(shannon, na.rm = TRUE),
              richness = mean(richness, na.rm = TRUE), 
              ML_NC_mean = mean(ML_NC, na.rm = TRUE))
  
  centroid_df <- SpatialPointsDataFrame(coords = late_som_maps[,c(6,5)], data = late_som_maps)
  centroid1 <- wt.centroid(x =centroid_df , p = 2)
  centroid2 <- wt.centroid(x =centroid_df , p = 3)
  
  late_dist <- distHaversine(centroid1@coords, centroid2@coords)/1000     
  
  compare <- inner_join(early_som_maps, late_som_maps, by = "Sta_ID")
  
  compare$NC_diff <- compare$NC_mean.x - compare$NC_mean.y
  compare$Temp_diff <- compare$temp_mean.x - compare$temp_mean.y
  compare$diversity_diff <- compare$shannon.x - compare$shannon.y
  compare$even_diff <- compare$evenness.x - compare$evenness.y
  compare$rich_diff <- compare$richness.x - compare$richness.y
  compare$ML_NC_diff <- compare$ML_NC_mean.x - compare$ML_NC_mean.y
  
  
  map <- map_data("world")  
  
  diversity_diff <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = compare, aes(x = long.x, y = lat.x, fill = diversity_diff),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(paste(Delta," Shannon Diversity")),low = "blue", high = "red", mid = "white", midpoint = 0) +
    ggtitle(paste0(title,"\nDifference in Diversity Early-Late")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  print(diversity_diff)
  
  even_diff <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = compare, aes(x = long.x, y = lat.x, fill = even_diff),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(paste(Delta," Evenness")),low = "blue", high = "red", mid = "white", midpoint = 0) +
    ggtitle(paste0(title,"\nDifference in Evenness Early-Late")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  print(even_diff)
  
  som_cruise <- full_dat %>% 
    group_by(Cruise) %>%
    summarise(som_1 = sum(som_id == 1, na.rm = TRUE)/n(), som_2 = sum(som_id == 2, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE),
              temp_mean = mean(T_degC, na.rm = TRUE), sal_mean = mean(Salnty, na.rm = TRUE),
              PO4_mean = mean(PO4ug, na.rm = TRUE), NO3_mean = mean(NO3ug, na.rm = TRUE),
              SiO3_mean = mean(SiO3ug, na.rm = TRUE), C14_mean = mean(IntC14, na.rm = TRUE),
              MLD_mean = mean(MLD_Sigma, na.rm = TRUE), NC_mean = mean(NCDepth, na.rm = TRUE),
              temp_coeff = sd(T_degC, na.rm = TRUE)/mean(T_degC, na.rm = TRUE),
              sal_coeff = sd(Salnty, na.rm = TRUE)/mean(Salnty, na.rm = TRUE),
              PO4_coeff =  sd(PO4ug, na.rm = TRUE)/mean(PO4ug, na.rm = TRUE),
              NO3_coeff = sd(NO3ug, na.rm = TRUE)/mean(NO3ug, na.rm = TRUE),
              SiO3_coeff = sd(SiO3ug, na.rm = TRUE)/mean(SiO3ug, na.rm = TRUE),
              C14_coeff = sd(IntC14, na.rm = TRUE)/mean(IntC14, na.rm = TRUE),
              MLD_coeff = sd(MLD_Sigma, na.rm = TRUE)/mean(MLD_Sigma, na.rm = TRUE),
              NC_coeff = sd(NCDepth, na.rm = TRUE)/mean(NCDepth, na.rm = TRUE),
              evenness = mean(evenness, na.rm = TRUE), shannon = mean(shannon, na.rm = TRUE),
              richness = mean(richness, na.rm = TRUE), 
              ML_NC_mean = mean(ML_NC, na.rm = TRUE), Date = mean(Date, na.rm = TRUE))
  
  # full_dat$dist_to_coast <- full_dat$dist_to_coast*1000
  
  results <-  full_dat %>% 
    group_by(Cruise) %>%
    do(model = lm(NCDepth ~ dist_to_coast, data = .)) %>%
    mutate(coef=coef(model)["dist_to_coast"])
  
  som_cruise$NC_slope <- results$coef
  
  som_cruise$phase <- c(rep("Early",12),rep("Late",8))
  som_cruise$season <- as.factor(rep(c("Winter", "Spring", "Summer", "Fall"),5))
  som_cruise$season <- factor(som_cruise$season, levels = c("Winter", "Spring", "Summer", "Fall"))
  
  gradient_plot <- ggplot(som_cruise, aes_string(x = "NC_slope", y = "shannon")) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") + 
    geom_point(size = 3, aes_string(color = "phase", shape = "season"), data = som_cruise) +
    scale_color_manual(values = c("red", "blue")) +
    scale_shape_manual(values = c(0,1,2,3)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("Nearshore-Offshore\nSlope in Nitracline") +
    ylab("Mean Shannon Diversity\nPer Cruise") + ggtitle(title)
  
  print(gradient_plot)

  ts_plot <- ggplot(som_cruise, aes_string(x = "Date", y = "shannon")) +
    geom_line(color = "black", size = 1) +
    geom_point(size = 3, aes_string(color = "season"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("Date") +
    ylab("Mean Shannon Diversity\nPer Cruise")  + ggtitle(title)

  print(ts_plot)

  load(upwelling_index)

  index_mat <- matrix(NA,nrow = nrow(som_cruise), 3)

  for (i in 1:nrow(som_cruise)) {

    year <- substr(som_cruise$Cruise[i],1,4)
    month <- substr(som_cruise$Cruise[i],5,6)

    index_yr <- substr(index_vals_mat$Date,1,4)
    index_month <- substr(index_vals_mat$Date,6,7)

    index_vals <- which(match(index_yr,year) & match(index_month,month))

    index_mat[i,] <- colMeans(index_vals_mat[index_vals,2:4], na.rm = TRUE)

  }

  som_cruise$CUTI <- index_mat[,1]
  som_cruise$BEUTI <- index_mat[,2]
  som_cruise$Nitrate <- index_mat[,3]

  som_cruise$log_BEUTI <- log(som_cruise$BEUTI)
  som_cruise$log_Nitrate <- log(som_cruise$Nitrate)

  cuti_plot <- ggplot(som_cruise, aes_string(x = "CUTI", y = "shannon")) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") +
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("Coastal Upwelling Transport Index\n(CUTI)") +
    ylab("Mean Shannon Diversity\nPer Cruise")

  beuti_plot <- ggplot(som_cruise, aes_string(x = "BEUTI", y = "shannon")) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") +
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") +
    xlab("Biologically Effective Upwelling Transport Index\n(BEUTI)") +
    ylab("")

  log_beuti_plot <- ggplot(som_cruise, aes_string(x = "log_BEUTI", y = "shannon")) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") +
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") +
    xlab("log(Biologically Effective Upwelling Transport Index)\n(BEUTI)") +
    ylab("")

  reg_nitrate <- ggplot(som_cruise, aes_string(x = "Nitrate", y = "shannon")) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") +
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("Regionally Availible Nitrate") +
    ylab("")

  log_reg_nitrate <- ggplot(som_cruise, aes_string(x = "log_Nitrate", y = "shannon")) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") +
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("log(Regionally Availible Nitrate)") +
    ylab("")


  title_plot <- ggdraw() + draw_label(title, fontface='bold')

  index_plots <- plot_grid(title_plot,
                           plot_grid(cuti_plot, beuti_plot, reg_nitrate, ncol = 3),
                           nrow = 2, rel_heights = c(0.1,1))

  pdf(index_plot, width = 12, height = 4)
  print(index_plots)
  dev.off()

  save(early_dist, late_dist, diversity_diff, even_diff, gradient_plot, ts_plot,
       cuti_plot, beuti_plot, reg_nitrate, file = out_diff_file)
  
}


in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9", "bacteria_m_euks_16s",
                  "plastid_16s", "cyano_16s", "euks_hetero_18sv9")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                   "Sar Clade", "Archaea","Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans", "Bacteria",
                   "Eukaryotic Phytoplankton (Plastids)", "Cyanobacteria",
                   "Eukaryotic Protists")

for (i in 1:length(in_group_list)) {
  
  diversity_comparison(in_file = paste0("output/",in_group_list[i],"_full_data.Rdata"),
                         in_map = paste0("output/",in_group_list[i],"_map.Rdata"),
                         out_diff_file = paste0("output/",in_group_list[i],"_diffs_div.Rdata"),
                         title = in_group_names[i],
                         upwelling_index = "output/upwelling_indicies.Rdata",
                         index_plot = paste0("figures/",in_group_list[i],"_index_plot.pdf"))
  
}

diversity_plots <- function(in_phyto = "output/plastid_16s_diffs_div.Rdata",
                              in_euks = "output/euks_hetero_18sv9_diffs_div.Rdata",
                              in_cyano = "output/cyano_16s_diffs_div.Rdata",
                              in_bact = "output/bacteria_m_euks_16s_diffs_div.Rdata",
                              in_phyto_full = "output/plastid_16s_map.Rdata",
                              in_euks_full = "output/euks_hetero_18sv9_map.Rdata",
                              in_cyano_full = "output/cyano_16s_map.Rdata",
                              in_bact_full = "output/bacteria_m_euks_16s_map.Rdata",
                              gradient_plot_file = "figures/gradient_plot_div.pdf",
                              time_series_plot = "figures/ts_plot_div.pdf"){

  dist <- as.data.frame(matrix(NA,4,3))
  colnames(dist) <- c("Early", "Late", "Diff")
  
  load(in_phyto)
  load(in_phyto_full)
  phyto_ts <- ts_plot
  phyto_gradient <- gradient_plot

  phyto_full <- som_maps
  
  
  load(in_euks)
  load(in_euks_full)
  euk_ts <- ts_plot
  euk_gradient <- gradient_plot
  euk_full <- som_maps
  
  load(in_cyano)
  load(in_cyano_full)
  cyano_ts <- ts_plot
  cyano_gradient <- gradient_plot
  cyano_full <- som_maps
  
  load(in_bact)
  load(in_bact_full)
  bact_ts <- ts_plot
  bact_gradient <- gradient_plot
  bact_full <- som_maps
  
  legs <- get_legend(phyto_gradient)
  phyto_gradient <- phyto_gradient + theme(legend.position = "none") + xlab("")  + 
    ylab("Mean Shannon Diversity\nPer Cruise")
  euk_gradient <- euk_gradient + theme(legend.position = "none") + xlab("") + ylab("")
  cyano_gradient <- cyano_gradient + theme(legend.position = "none")  + 
    ylab("Mean Shannon Diversity\nPer Cruise") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
  bact_gradient <- bact_gradient + theme(legend.position = "none") + ylab("") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
  
  grad_plot <- plot_grid(plot_grid(phyto_gradient, euk_gradient,
                                   cyano_gradient, bact_gradient),
                         legs, ncol = 2, rel_widths = c(1,0.2))
  
  pdf(file = gradient_plot_file, width = 9, height = 7)
  print(grad_plot)
  dev.off()
  
  ts_leg <- get_legend(phyto_ts)
  phyto_ts <- phyto_ts + theme(legend.position = "none") + xlab("")
  euk_ts <- euk_ts + theme(legend.position = "none") + xlab("")
  cyano_ts <- cyano_ts + theme(legend.position = "none") + xlab("")
  bact_ts <- bact_ts + theme(legend.position = "none")

  ts_plot <- plot_grid(plot_grid(phyto_ts, euk_ts, cyano_ts, bact_ts, nrow = 4), ts_leg,
                       ncol = 2, rel_widths = c(1,0.2))

  pdf(file = time_series_plot, width = 8, height = 12)
  print(ts_plot)
  dev.off()
  
}


