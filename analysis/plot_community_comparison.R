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

community_comparison <- function(in_file = "output/euks_auto_18sv9_full_data.Rdata",
                                   similar_mat = "output/euks_auto_18sv9_dissimilar.Rdata",
                                   in_map = "output/euks_auto_18sv9_map.Rdata",
                                   out_diff_file = "output/euks_auto_18sv9_diffs.Rdata",
                                  title = "Eukaryotic Phytoplankton",
                                   upwelling_index = "output/upwelling_indicies.Rdata",
                                   index_plot = "figures/euks_auto_18sv9_index_plot.pdf"){
  
  load(in_file)
  load(similar_mat)
  load(in_map)
  
  # find centroids
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  wt_1 <- wt.centroid(x = centroid_df, p = 2)
  wt_2 <- wt.centroid(x = centroid_df, p = 3)
  
  clust1 <- which.max(c(wt_1@coords[1], wt_2@coords[1]))
  clust2 <- which.min(c(wt_1@coords[1], wt_2@coords[1]))
  
  if(clust1 == 1){nearshore_som <- 1}
  if(clust1 == 2){nearshore_som <- 2}
  
  full_dat$Year <- as.numeric(substr(full_dat$Cruise,1,4))
  
  full_dat <- full_dat[-which(as.numeric(substr(full_dat$Sta_ID,2,3)) < 76),]
  
  full_dat$ML_NC <- full_dat$MLD_Sigma - full_dat$NCDepth
  
  early_som_maps <- full_dat %>% 
    filter(Year < 2017 | Year == 2019) %>%
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
    filter(Year > 2016 & Year < 2019) %>%
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
 
 som_cruise$phase <- c(rep("2014-2016",12),rep("2017-2018",8),rep("2019",4))
 som_cruise$season <- as.factor(rep(c("Winter", "Spring", "Summer", "Fall"),6))
 som_cruise$season <- factor(som_cruise$season, levels = c("Winter", "Spring", "Summer", "Fall"))
  
  gradient_plot <- ggplot(som_cruise, aes_string(x = "NC_slope", y = paste0("som_",nearshore_som))) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") + 
    geom_point(size = 3, aes_string(color = "phase", shape = "season"), data = som_cruise) +
    scale_color_manual(values = c("red", "blue","gold3")) +
    scale_shape_manual(values = c(0,1,2,3)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Season", color = "Phase") + xlab("Nearshore-Offshore\nSlope in Nitracline") +
    ylab("Frequency of Nearshore Cluster") + ggtitle(title)
  
  print(gradient_plot)
                          
  ts_plot <- ggplot(som_cruise, aes_string(x = "Date", y = paste0("som_",nearshore_som))) +
    geom_line(color = "black", size = 1) + 
    geom_point(size = 3, aes_string(color = "season"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("Date") +
    ylab("Frequency of Nearshore Cluster")  + ggtitle(title)
  
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
  
 cuti_plot <- ggplot(som_cruise, aes_string(x = "CUTI", y = paste0("som_",nearshore_som))) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") + 
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17,18)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("Coastal Upwelling Transport Index\n(CUTI)") +
    ylab("Frequency of Nearshore Cluster") 
  
  beuti_plot <- ggplot(som_cruise, aes_string(x = "BEUTI", y = paste0("som_",nearshore_som))) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") + 
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17,18)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") + 
    xlab("Biologically Effective Upwelling Transport Index\n(BEUTI)") +
    ylab("") 
  
  log_beuti_plot <- ggplot(som_cruise, aes_string(x = "log_BEUTI", y = paste0("som_",nearshore_som))) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") + 
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17,18)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") + 
    xlab("log(Biologically Effective Upwelling Transport Index)\n(BEUTI)") +
    ylab("") 
  
  reg_nitrate <- ggplot(som_cruise, aes_string(x = "Nitrate", y = paste0("som_",nearshore_som))) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") + 
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17,18)) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    labs(shape = "Phase", color = "Season") + xlab("Regionally Availible Nitrate") +
    ylab("") 
  
  log_reg_nitrate <- ggplot(som_cruise, aes_string(x = "log_Nitrate", y = paste0("som_",nearshore_som))) +
    # geom_smooth(method = 'glm', formula = y~x, se = FALSE, color = "black") + 
    geom_point(size = 3, aes_string(color = "season", shape = "phase"), data = som_cruise) +
    scale_color_manual(values = c("slategray3", "springgreen3", "gold3", "darkorange3")) +
    scale_shape_manual(values = c(15,17,18)) + 
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
  
  diss_melt <- melt(dissimilar)
  diss_melt <- diss_melt[-which(diss_melt$Var1 == diss_melt$Var2),]
  
  diss_melt$Var1_som <- full_dat$som_id[match(diss_melt$Var1,full_dat$eco_name)]
  diss_melt$Var2_som <- full_dat$som_id[match(diss_melt$Var2,full_dat$eco_name)]
  
  diss_melt$Var1_phase <- NA
  diss_melt$Var2_phase <- NA
  
  diss_melt$Var1_phase[which(substr(diss_melt$Var1,2,5) < 2017)] <- "Early"
  diss_melt$Var2_phase[which(substr(diss_melt$Var2,2,5) < 2017)] <- "Early"
  
  diss_melt$Var1_phase[which(substr(diss_melt$Var1,2,5) > 2016 & substr(diss_melt$Var1,2,5) < 2019)] <- "Late"
  diss_melt$Var2_phase[which(substr(diss_melt$Var2,2,5) > 2016 & substr(diss_melt$Var1,2,5) < 2019)] <- "Late"
  
  diss_melt$Var1_phase[which(substr(diss_melt$Var1,2,5) > 2018)] <- "2019"
  diss_melt$Var2_phase[which(substr(diss_melt$Var1,2,5) > 2018)] <- "2019"
  
  diss_melt$comp_phase <- paste0(diss_melt$Var1_phase,"-",diss_melt$Var2_phase)
  diss_melt$comp_som <- paste0(diss_melt$Var1_som,"-",diss_melt$Var2_som)
  
  diss_filt <- filter(diss_melt, comp_som == "1-1" | comp_som == "2-2")
  diss_filt <- filter(diss_filt, comp_phase == "Early-Early" | comp_phase == "Late-Late" | comp_phase == "2019-2019")
  
  diss_filt$comp_phase[diss_filt$comp_phase == "Early-Early"] = "2014-2016"
  diss_filt$comp_phase[diss_filt$comp_phase == "Late-Late"] = "2017-2018"
  diss_filt$comp_phase[diss_filt$comp_phase == "2019-2019"] = "2019"
  
  if(nearshore_som == 1){
    diss_filt$comp_som[diss_filt$comp_som == "1-1"] = "Nearshore"
    diss_filt$comp_som[diss_filt$comp_som == "2-2"] = "Offshore"
  }
  if(nearshore_som == 2){
    diss_filt$comp_som[diss_filt$comp_som == "1-1"] = "Offshore"
    diss_filt$comp_som[diss_filt$comp_som == "2-2"] = "Nearshore"
  }
 
  
  dissimilar_plot <- ggplot(diss_filt,
         aes(y = 1-value, x = interaction(comp_som,comp_phase),fill = comp_phase)) +
    geom_boxplot() + ylab("Bray-Curtis Similarity") + xlab("Phase/Cluster") +
    labs(fill = "Phase") + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = c("Nearshore\n2014-2016", "Offshore\n2014-2016",
                                "Nearshore\n2017-2018", "Offshore\n2017-2018",
                                "Nearshore\n2019", "Offshore\n2019")) +
    ggtitle(title)
  
  print(dissimilar_plot)
  
  save(early_dist, late_dist, diversity_diff, even_diff, gradient_plot, ts_plot,
       cuti_plot, beuti_plot, reg_nitrate, dissimilar_plot, file = out_diff_file)
  
  
  
  }


# 16s

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "bacteria_m_euks_16s", "plastid_16s", "cyano_16s")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                   "Sar Clade", "Archaea", "Bacteria",
                   "Eukaryotic Phytoplankton (Plastids)", "Cyanobacteria")

# All


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
  
  community_comparison(in_file = paste0("output/",in_group_list[i],"_full_data.Rdata"),
                         in_map = paste0("output/",in_group_list[i],"_map.Rdata"),
                         similar_mat = paste0("output/",in_group_list[i],"_dissimilar.Rdata"),
                         out_diff_file = paste0("output/",in_group_list[i],"_diffs.Rdata"),
                         title = in_group_names[i],
                         upwelling_index = "output/upwelling_indicies.Rdata",
                         index_plot = paste0("figures/",in_group_list[i],"_index_plot.pdf"))
  
}




community_plots <- function(in_phyto = "output/plastid_16s_diffs.Rdata",
                              in_euks = "output/euks_hetero_18sv9_diffs.Rdata",
                              in_cyano = "output/cyano_16s_diffs.Rdata",
                              in_bact = "output/bacteria_m_euks_16s_diffs.Rdata",
                              in_phyto_full = "output/plastid_16s_map.Rdata",
                              in_euks_full = "output/euks_hetero_18sv9_map.Rdata",
                              in_cyano_full = "output/cyano_16s_map.Rdata",
                              in_bact_full = "output/bacteria_m_euks_16s_map.Rdata",
                              gradient_plot_file = "figures/gradient_plot.pdf",
                              time_series_plot = "figures/ts_plot.pdf",
                              diff_plots = "figures/bray_curt_sim.pdf",
                              div_plots = "figures/diversity_diff.pdf",
                              even_plots = "figures/evenness_diff.pdf"){
  
  dist <- as.data.frame(matrix(NA,4,3))
  colnames(dist) <- c("Early", "Late", "Diff")
  
  load(in_phyto)
  load(in_phyto_full)
  dist[1,] <- c(early_dist, late_dist, early_dist - late_dist)
  phyto_div <- diversity_diff
  phyto_even <- even_diff
  phyto_ts <- ts_plot
  phyto_gradient <- gradient_plot
  phyto_diss <- dissimilar_plot
  phyto_full <- som_maps
  
  
  load(in_euks)
  load(in_euks_full)
  dist[2,] <- c(early_dist, late_dist, early_dist - late_dist)
  euk_div <- diversity_diff
  euk_even <- even_diff
  euk_ts <- ts_plot
  euk_gradient <- gradient_plot
  euk_diss <- dissimilar_plot
  euk_full <- som_maps
  
  load(in_cyano)
  load(in_cyano_full)
  dist[3,] <- c(early_dist, late_dist, early_dist - late_dist)
  cyano_div <- diversity_diff
  cyano_even <- even_diff
  cyano_ts <- ts_plot
  cyano_gradient <- gradient_plot
  cyano_diss <- dissimilar_plot
  cyano_full <- som_maps
  
  load(in_bact)
  load(in_bact_full)
  dist[4,] <- c(early_dist, late_dist, early_dist - late_dist)
  bact_div <- diversity_diff
  bact_even <- even_diff
  bact_ts <- ts_plot
  bact_gradient <- gradient_plot
  bact_diss <- dissimilar_plot
  bact_full <- som_maps
  
  legs <- get_legend(phyto_gradient)
  phyto_gradient <- phyto_gradient + theme(legend.position = "none") + xlab("")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore")
  euk_gradient <- euk_gradient + theme(legend.position = "none") + xlab("") + ylab("")
  cyano_gradient <- cyano_gradient + theme(legend.position = "none")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
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
  
  
  diss_leg <- get_legend(phyto_diss)
  phyto_diss <- phyto_diss + theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
    xlab("") + ggtitle("Eukaryotic Phytoplankton")
  euk_diss <- euk_diss +
    theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
    xlab("") + ggtitle("Heterotrophic Eukaryotes") + ylab("")
  cyano_diss <- cyano_diss + theme(legend.position = "none") +
    theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  ggtitle("Cyanobacteria")
  bact_diss <- bact_diss + theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
    ggtitle("Heterotrophic Bacteria/Archaea") + ylab("")
  
  diss_plot <- plot_grid(plot_grid(phyto_diss, euk_diss,
                                   cyano_diss, bact_diss),
                         diss_leg, ncol = 2, rel_widths = c(1,0.2))
  
  pdf(file = diff_plots, width = 8, height = 6)
  print(diss_plot)
  dev.off()
  
  map <- map_data("world")  
  
  phyto_tot <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = phyto_full, aes(x = long, y = lat, fill = shannon),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = "Mean Alpha\nDiversity", limits = c(3.25,3.75), oob = scales::squish,
                         low = "blue", high = "red", mid = "white", midpoint = 3.5) +
    ggtitle("Mean Alpha Diversity (2014-2018)") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  euk_tot <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = euk_full, aes(x = long, y = lat, fill = shannon),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = "Mean Alpha\nDiversity", limits = c(2.5,4.5), oob = scales::squish,
                         low = "blue", high = "red", mid = "white", midpoint = 3.5) +
    ggtitle("") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  cyano_tot <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = cyano_full, aes(x = long, y = lat, fill = shannon),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = "Mean Alpha\nDiversity", limits = c(2,3), oob = scales::squish,
                         low = "blue", high = "red", mid = "white", midpoint = 2.5) +
    ggtitle("") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  bact_tot <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = bact_full, aes(x = long, y = lat, fill = shannon),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = "Mean Alpha\nDiversity", limits = c(4.25,4.75), oob = scales::squish,
                         low = "blue", high = "red", mid = "white", midpoint = 4.5) +
    ggtitle("") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  
  div_leg <- get_legend(phyto_div)
  phyto_div <- phyto_div +
    theme(legend.position = "none") + ggtitle(expression(paste(Delta," Alpha Diversity (Warm - Cool)")))
  euk_div <- euk_div + 
    theme(legend.position = "none") + ggtitle("")
  cyano_div <- cyano_div + 
    theme(legend.position = "none") + ggtitle("")
  bact_div <- bact_div + 
    theme(legend.position = "none") + ggtitle("")
  
  phyto_title <- ggdraw() + draw_label("Eukaryotic\nPhytoplankton", fontface='bold')
  euk_title <- ggdraw() + draw_label("Heterotrophic\nEukaryotes", fontface='bold')
  cyano_title <- ggdraw() + draw_label("Cyanobacteria", fontface='bold')
  bact_title <- ggdraw() + draw_label("Heterotrophic\nBacteria/Archaea", fontface='bold')
  
  
  div_plot <- plot_grid(
    plot_grid(phyto_title, euk_title,
              cyano_title, bact_title, ncol = 1),
    plot_grid(phyto_tot, euk_tot,
              cyano_tot, bact_tot, ncol = 1),
    plot_grid(phyto_div, euk_div,
              cyano_div, bact_div, ncol = 1),
                         div_leg, ncol = 4, rel_widths = c(0.3,1,1,0.3))
  
  pdf(file = div_plots, width = 14, height = 14)
  print(div_plot)
  dev.off()
  
  
  phyto_et <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = phyto_full, aes(x = long, y = lat, fill = evenness),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = "Mean\nEvenness", limits = c(0.65,0.85), oob = scales::squish,
                         low = "blue", high = "red", mid = "white", midpoint = .75) +
    ggtitle("Mean Evenness (2014-2018)") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  euk_et <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = euk_full, aes(x = long, y = lat, fill = evenness),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = "Mean\nEvenness", limits = c(0.55,0.75), oob = scales::squish,
                         low = "blue", high = "red", mid = "white", midpoint = .65) +
    ggtitle("") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  cyano_et <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = cyano_full, aes(x = long, y = lat, fill = evenness),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = "Mean\nEvenness", limits = c(0.525,0.675), oob = scales::squish,
                         low = "blue", high = "red", mid = "white", midpoint = .6) +
    ggtitle("") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  bact_et <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = bact_full, aes(x = long, y = lat, fill = evenness),
               color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = "Mean\nEvenness", limits = c(0.725,0.775), oob = scales::squish,
                         low = "blue", high = "red", mid = "white", midpoint = .75) +
    ggtitle("") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  
  
  even_leg <- get_legend(phyto_even)
  phyto_even <- phyto_even + theme(legend.position = "none") +
    ggtitle(expression(paste(Delta," Evenness (Warm - Cool)"))) 
  euk_even <- euk_even +
    theme(legend.position = "none") + ggtitle("")
  cyano_even <- cyano_even + 
    theme(legend.position = "none") + ggtitle("")
  bact_even <- bact_even + 
    theme(legend.position = "none")+ ggtitle("")
  
  even_plot <- plot_grid(
    plot_grid(phyto_title, euk_title,
              cyano_title, bact_title, ncol = 1),
    plot_grid(phyto_et, euk_et,
              cyano_et, bact_et, ncol = 1),
    plot_grid(phyto_even, euk_even,
              cyano_even, bact_even, ncol = 1),
    even_leg, ncol = 4, rel_widths = c(0.3,1,1,0.3))
  
  pdf(file = even_plots, width = 14, height = 14)
  print(even_plot)
  dev.off()
  
  }


