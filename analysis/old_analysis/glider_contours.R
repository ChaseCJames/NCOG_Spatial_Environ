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
library(pracma)


glider_contours <- function(in_glider = "output/spray_data_80_90.Rdata",
                            in_data = "output/euks_auto_18sv9_full_data.Rdata",
                            isopycnal = 26){
  
  load(in_glider)
  load(in_data)
  
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
              richness = mean(richness, na.rm = TRUE), Date = mean(Date, na.rm = TRUE))
  
  
    line_90_contours <- matrix(NA, nrow = nrow(som_cruise), ncol = length(dist_90))
    line_80_contours <- matrix(NA, nrow = nrow(som_cruise), ncol = length(dist_80))
    line_90_mld_contour <- matrix(NA, nrow = nrow(som_cruise), ncol = length(dist_90))
    line_80_mld_contour <- matrix(NA, nrow = nrow(som_cruise), ncol = length(dist_80))
  
    for (i in 1:nrow(som_cruise)) {
      
      date1 <- som_cruise$Date[i]
      date2 <- seq(date1, length = 2, by = "-2 months")[2]
      possible_days <- seq.Date(date2,date1, by ="days")
      
      mat_80 <- matrix(NA,50,74)
      
      for (j in 1:length(depth_80)) {
        for (k in 1:length(dist_80)) {
          
          mat_80[j,k] <- mean(pd_80_array[j,which(!is.na(match(ts_80, possible_days))),k])
          
        }
      }
      
      mat_90 <- matrix(NA,50,107)
      
      for (j in 1:length(depth_90)) {
        for (k in 1:length(dist_90)) {
          
          mat_90[j,k] <- mean(pd_90_array[j,which(!is.na(match(ts_90, possible_days))),k])
          
        }
      }
    
      int_mat_80 <- matrix(NA, length(10:500),length(dist_80))
      rownames(int_mat_80) <- 10:500
      for (j in 1:ncol(mat_80)) {
      if(length(which(!is.na(mat_80[,j]))) > 10){
        int_mat_80[,j] <- interp1(x = as.numeric(depth_80), y = as.numeric(mat_80[,j]), xi = min(depth_80):max(depth_80))
      }else{}
      }
      
      int_mat_90 <- matrix(NA, length(10:500),length(dist_90))
      rownames(int_mat_90) <- 10:500
      for (j in 1:ncol(mat_90)) {
        if(length(which(!is.na(mat_90[,j]))) > 10){
          int_mat_90[,j] <- interp1(x = as.numeric(depth_90), y = as.numeric(mat_90[,j]), xi = min(depth_90):max(depth_90))
        }else{}
      }

      line_80_contours[i,] <- as.numeric(rownames(int_mat_80)[apply(int_mat_80, 2, function(x) min(which(x > isopycnal)))])
      line_80_mld_contour[i,] <- as.numeric(rownames(int_mat_80)[apply(int_mat_80, 2, function(x) min(which(x > (x[1]+0.03))))])
      
      line_90_contours[i,] <- as.numeric(rownames(int_mat_90)[apply(int_mat_90, 2, function(x) min(which(x > isopycnal)))])
      line_90_mld_contour[i,] <- as.numeric(rownames(int_mat_90)[apply(int_mat_90, 2, function(x) min(which(x > (x[1]+0.03))))])
      
      }
  
  colnames(line_80_contours) <- dist_80
  colnames(line_90_contours) <- dist_90
  colnames(line_80_mld_contour) <- dist_80
  colnames(line_90_mld_contour) <- dist_90
  
  rownames(line_80_contours) <- som_cruise$Cruise
  rownames(line_90_contours) <- som_cruise$Cruise
  rownames(line_80_mld_contour) <- som_cruise$Cruise
  rownames(line_90_mld_contour) <- som_cruise$Cruise
  
  contour_80_melt <- melt(line_80_contours)
  contour_90_melt <- melt(line_90_contours)
  contour_80_mld_melt <- melt(line_80_mld_contour)
  contour_90_mld_melt <- melt(line_90_mld_contour)
  
  colnames(contour_80_melt) <- c("Cruise", "Dist", "Depth")
  colnames(contour_90_melt) <- c("Cruise", "Dist", "Depth")
  colnames(contour_80_mld_melt ) <- c("Cruise", "Dist", "Depth")
  colnames(contour_90_mld_melt) <- c("Cruise", "Dist", "Depth")
  
  early <- unique(contour_80_melt$Cruise)[1:12]
  late <- unique(contour_80_melt$Cruise)[13:20]
  
  winter <- unique(contour_80_melt$Cruise)[seq(1,20,by = 4)]
  spring <- unique(contour_80_melt$Cruise)[seq(2,20,by = 4)]
  summer <- unique(contour_80_melt$Cruise)[seq(3,20,by = 4)]
  fall <- unique(contour_80_melt$Cruise)[seq(4,20,by = 4)]
  
  contour_80_melt$Phase <- rep(NA, nrow(contour_80_melt))
  contour_80_melt$Season <- rep(NA, nrow(contour_80_melt))
  
  contour_90_melt$Phase <- rep(NA, nrow(contour_90_melt))
  contour_90_melt$Season <- rep(NA, nrow(contour_90_melt))
  
  contour_80_mld_melt$Phase <- rep(NA, nrow(contour_80_melt))
  contour_80_mld_melt$Season <- rep(NA, nrow(contour_80_melt))
  
  contour_90_mld_melt$Phase <- rep(NA, nrow(contour_90_melt))
  contour_90_mld_melt$Season <- rep(NA, nrow(contour_90_melt))
  
  # phase
  
  contour_80_melt$Phase[which(!is.na(match(contour_80_melt$Cruise, early)))] <- "Early"
  contour_80_melt$Phase[which(is.na(contour_80_melt$Phase))] <- "Late"
  
  contour_90_melt$Phase[which(!is.na(match(contour_90_melt$Cruise, early)))] <- "Early"
  contour_90_melt$Phase[which(is.na(contour_90_melt$Phase))] <- "Late"
  
  contour_80_mld_melt$Phase[which(!is.na(match(contour_80_mld_melt$Cruise, early)))] <- "Early"
  contour_80_mld_melt$Phase[which(is.na(contour_80_mld_melt$Phase))] <- "Late"
  
  contour_90_mld_melt$Phase[which(!is.na(match(contour_90_mld_melt$Cruise, early)))] <- "Early"
  contour_90_mld_melt$Phase[which(is.na(contour_90_mld_melt$Phase))] <- "Late"
  
  # season
  
  contour_80_melt$Season[which(!is.na(match(contour_80_melt$Cruise, winter)))] <- "Winter"
  contour_80_melt$Season[which(!is.na(match(contour_80_melt$Cruise, spring)))] <- "Spring"
  contour_80_melt$Season[which(!is.na(match(contour_80_melt$Cruise, summer)))] <- "Summer"
  contour_80_melt$Season[which(!is.na(match(contour_80_melt$Cruise, fall)))] <- "Fall"
  
  contour_80_mld_melt$Season[which(!is.na(match(contour_80_mld_melt$Cruise, winter)))] <- "Winter"
  contour_80_mld_melt$Season[which(!is.na(match(contour_80_mld_melt$Cruise, spring)))] <- "Spring"
  contour_80_mld_melt$Season[which(!is.na(match(contour_80_mld_melt$Cruise, summer)))] <- "Summer"
  contour_80_mld_melt$Season[which(!is.na(match(contour_80_mld_melt$Cruise, fall)))] <- "Fall"
  
  
  line_80 <- ggplot(contour_80_melt, aes(x = Dist, y = Depth, group = interaction(Cruise, Phase),color = Phase)) +
    geom_line(size = 1) +
    theme(panel.background = element_blank(),
                        panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    scale_x_reverse() + scale_y_reverse() +
    xlab("Distance to Coast (km)") + ylab(paste0("Depth of ",isopycnal," Isopycnal")) + 
    ggtitle("Line 80")
  
  
  line_80_mld <- ggplot(contour_80_mld_melt, aes(x = Dist, y = Depth, group = interaction(Cruise, Phase),color = Phase)) +
    geom_line(size = 1) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    scale_x_reverse() + scale_y_reverse() +
    xlab("Distance to Coast (km)") + ylab("Depth of Mixed Layer") + 
    ggtitle("Line 80")
  
  line_90 <- ggplot(contour_90_melt, aes(x = Dist, y = Depth, group = interaction(Cruise, Phase),color = Phase)) +
    geom_line(size = 1) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    scale_x_reverse() + scale_y_reverse() +
    xlab("Distance to Coast (km)") + ylab(paste0("Depth of ",isopycnal," Isopycnal")) + 
    ggtitle("Line 90")
  
  line_90_mld <- ggplot(contour_90_mld_melt, aes(x = Dist, y = Depth, group = interaction(Cruise, Phase),color = Phase)) +
    geom_line(size = 1) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    scale_x_reverse() + scale_y_reverse() +
    xlab("Distance to Coast (km)") + ylab("Depth of Mixed Layer") + 
    ggtitle("Line 90")
  
  pdf("figures/contours_80_90_phase.pdf", width = 12, height = 4)
  plot_grid(line_80, line_90)
  dev.off()
  
  # isopycnal
  
  winter <- ggplot(filter(contour_80_melt, Season == "Winter"), aes(x = Dist, y = Depth, group = interaction(Cruise, Phase),color = Phase)) +
    geom_line(size = 1) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    scale_x_reverse() + scale_y_reverse(limits = c(200,0)) +
    xlab("Distance to Coast (km)") + ylab(paste0("Depth of ",isopycnal," Isopycnal")) + 
    ggtitle("Winter Cruises")
  
  spring <- ggplot(filter(contour_80_melt, Season == "Spring"), aes(x = Dist, y = Depth, group = interaction(Cruise, Phase),color = Phase)) +
    geom_line(size = 1) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    scale_x_reverse() + scale_y_reverse(limits = c(200,0)) +
    xlab("Distance to Coast (km)") + ylab(paste0("Depth of ",isopycnal," Isopycnal")) + 
    ggtitle("Spring Cruises")
  
  summer <- ggplot(filter(contour_80_melt, Season == "Summer"), aes(x = Dist, y = Depth, group = interaction(Cruise, Phase),color = Phase)) +
    geom_line(size = 1) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    scale_x_reverse() + scale_y_reverse(limits = c(200,0)) +
    xlab("Distance to Coast (km)") + ylab(paste0("Depth of ",isopycnal," Isopycnal")) + 
    ggtitle("Summer Cruises")
  
  fall <- ggplot(filter(contour_80_melt, Season == "Fall"), aes(x = Dist, y = Depth, group = interaction(Cruise, Phase),color = Phase)) +
    geom_line(size = 1) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    scale_x_reverse() + scale_y_reverse(limits = c(200,0)) +
    xlab("Distance to Coast (km)") + ylab(paste0("Depth of ",isopycnal," Isopycnal")) + 
    ggtitle("Fall Cruises")
  
  # MLD
  
  winter_mld <- ggplot(filter(contour_80_mld_melt, Season == "Winter"), aes(x = Dist, y = Depth, group = interaction(Cruise, Phase),color = Phase)) +
    geom_line(size = 1) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    scale_x_reverse() + scale_y_reverse(limits = c(100,0)) +
    xlab("Distance to Coast (km)") + ylab("Depth of Mixed Layer") + 
    ggtitle("Winter Cruises")
  
  spring_mld <- ggplot(filter(contour_80_mld_melt, Season == "Spring"), aes(x = Dist, y = Depth, group = interaction(Cruise, Phase),color = Phase)) +
    geom_line(size = 1) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    scale_x_reverse() + scale_y_reverse(limits = c(100,0)) +
    xlab("Distance to Coast (km)") + ylab("Depth of Mixed Layer") + 
    ggtitle("Spring Cruises")
  
  summer_mld <- ggplot(filter(contour_80_mld_melt, Season == "Summer"), aes(x = Dist, y = Depth, group = interaction(Cruise, Phase),color = Phase)) +
    geom_line(size = 1) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    scale_x_reverse() + scale_y_reverse(limits = c(100,0)) +
    xlab("Distance to Coast (km)") + ylab("Depth of Mixed Layer") + 
    ggtitle("Summer Cruises")
  
  fall_mld <- ggplot(filter(contour_80_mld_melt, Season == "Fall"), aes(x = Dist, y = Depth, group = interaction(Cruise, Phase),color = Phase)) +
    geom_line(size = 1) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(hjust=0.5)) +
    scale_x_reverse() + scale_y_reverse(limits = c(100,0)) +
    xlab("Distance to Coast (km)") + ylab("Depth of Mixed Layer") + 
    ggtitle("Fall Cruises")
  
  pdf("figures/line_80_contours_season.pdf", width = 12, height = 12)
  plot_grid(winter,spring,summer,fall, ncol = 2)
  dev.off()
  
  pdf("figures/line_80_contours_season_mld.pdf", width = 12, height = 12)
  plot_grid(winter_mld,spring_mld,summer_mld,fall_mld, ncol = 2)
  dev.off()
  
  contour_80_melt$Type <- "Isopycnal"
  contour_80_mld_melt$Type <- "MLD"
  
  all_80_contours <- bind_rows(contour_80_melt, contour_80_mld_melt)
  
  cruise_list <- list()
  
  for (i in 1:nrow(som_cruise)) {
    
    cruise_plot <- ggplot() +
      geom_line(data = filter(all_80_contours, Cruise == som_cruise$Cruise[i]),
                aes(x = Dist, y = Depth, linetype = Type)) +
      geom_point(data = filter(full_dat, Cruise == som_cruise$Cruise[i] & substr(full_dat$Sta_ID,2,3) == 80),
                 aes(x = dist_to_coast, y = Depthm, color = as.factor(som_id)), size = 4) +
      scale_color_manual(name = "SOM Cluster", values = c("red", "blue")) +
      theme(panel.background = element_blank(),
            panel.border = element_rect(fill = NA, color = "black"),
            plot.title = element_text(hjust=0.5)) +
      scale_x_reverse() + scale_y_reverse(limits = c(200,0)) +
      xlab("Distance to Coast (km)") + ylab("Depth (m)") +
      ggtitle(som_cruise$Cruise[i])
    
    print(cruise_plot)
    
    cruise_list[[i]] <- cruise_plot
    
  }
  
  pdf("figures/cruise_mld_iso_som.pdf", width = 20, height = 25)
  plot_grid(cruise_list[[1]],cruise_list[[2]], cruise_list[[3]], cruise_list[[4]],
            cruise_list[[5]],cruise_list[[6]], cruise_list[[7]], cruise_list[[8]],
            cruise_list[[9]],cruise_list[[10]], cruise_list[[11]], cruise_list[[12]],
            cruise_list[[13]],cruise_list[[14]], cruise_list[[15]], cruise_list[[16]],
            cruise_list[[17]],cruise_list[[18]], cruise_list[[19]], cruise_list[[20]],
            ncol = 4)
  dev.off()
  
  }