library(tidyverse)
library(ggmap)
library(patchwork)

# Suppl Fig X - Diversity between the two phases

in_group_names = c("Archaea","Bacteria", 
                   "Cyanobacteria",
                   "Photosynthetic Eukaryotic\nProtists",
                   "Heterotrophic\nEukaryotic Protists")

in_group_list = c( "archaea_16s", "bacteria_m_euks_16s",
                  "cyano_16s","euks_auto_18sv9","euks_hetero_18sv9")

fig_list <- list()

map <- map_data("world")  

for (i in 1:length(in_group_names)) {
  
  load(paste0("output/", in_group_list[i] ,"_full_data_S.Rdata"))
  
  early <- full_dat %>% filter(as.numeric(substr(Cruise,1,4)) < 2017) %>%
    group_by(Sta_ID) %>% summarise(mean_sha = mean(shannon, na.rm = TRUE))
  
  late <- full_dat %>% filter(as.numeric(substr(Cruise,1,4)) > 2016,
                              as.numeric(substr(Cruise,1,4)) < 2019) %>%
    group_by(Sta_ID) %>% summarise(mean_sha = mean(shannon, na.rm = TRUE))
  
  late <- late %>% filter(Sta_ID %in% early$Sta_ID)
  early <- early %>% filter(Sta_ID %in% late$Sta_ID)
  
  early$late_mean <- late$mean_sha[match(early$Sta_ID, late$Sta_ID)]
  
  early$diff <- early$mean_sha - early$late_mean
  
  early$lat <- full_dat$Lat_Dec[match(early$Sta_ID, full_dat$Sta_ID)]
  early$long <- full_dat$Lon_Dec[match(early$Sta_ID, full_dat$Sta_ID)]
  
  max_lim <- max(c(max(early$diff), abs(min(early$diff))))
  
  fig_list[[i]] <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early, aes(x = long, y = lat, fill = diff), color = "black", size = 3, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(low = "darkred", high = "darkblue", mid = "white", midpoint =  0, limits = c(-max_lim,max_lim)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          axis.line = element_blank(),
          axis.text = element_text(size = 8),
          legend.text = element_text(size = 8),
          axis.title = element_text(size = 8),
          legend.title = element_text(size = 8)) +
    ggtitle(in_group_names[i]) + labs(fill = expression(paste(Delta,H)))
  
}


out_plot <- fig_list[[1]] + fig_list[[2]] + 
  fig_list[[3]]  + fig_list[[4]]  +
  fig_list[[5]] + plot_layout(nrow = 1)

pdf("figures_S/supp_fig_X_S.pdf", width = 15, height = 3)
plot(out_plot)
dev.off()



