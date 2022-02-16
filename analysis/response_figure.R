library(ragg)
source("analysis/figure_functions_S.R")

in_group_list = c( "archaea_16s", "bacteria_m_euks_16s",
                  "cyano_16s", "euks_hetero_18sv9", "euks_auto_18sv9", "total")

in_group_names = c("Archaea", "Bacteria",
                   "Cyanobacteria","Eukaryotic Protists", "Photosynthetic Eukaryotic Protists", "All ASVs")

in_group_list_basic = c("16s_archaea", "16s_bacteria_m_euks",
                        "16s_cyanos", "18s_heterotrophic_euks", "18s_autotrophic_euks", "totals")

letters <- c("b", "c", "d", "e", "f","a")


diveristy_figure_response <- function(map_file = paste0("output/", in_group_list[i], "_map_S.Rdata"),
                             full_dat = paste0("output/", in_group_list[i], "_full_data_S.Rdata"),
                             figure_start = paste0("figures/diversity/", in_group_list[i], "_"),
                             main = in_group_names[i]){
  
  
  map <- map_data("world")   
  
  load(map_file)
  load(full_dat)
  

  surf <- full_dat %>% filter(Depthm <= 15) %>% group_by(Sta_ID) %>%
    summarise(lat = mean(Lat_Dec, na.rm = TRUE), long = mean(Lon_Dec, na.rm = TRUE), shannon = mean(shannon, na.rm = TRUE))
  
  dcm <- full_dat %>% filter(Depthm > 15) %>% group_by(Sta_ID) %>%
    summarise(lat = mean(Lat_Dec, na.rm = TRUE), long = mean(Lon_Dec, na.rm = TRUE), shannon = mean(shannon, na.rm = TRUE))
  
  max_v <- max(c(dcm$shannon, surf$shannon))
  min_v <- min(c(dcm$shannon, surf$shannon))
  
  shannon_surf <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = surf, aes(x = long, y = lat, fill = shannon), color = "black", size =4, stroke = 0.1, shape = 21) +
    scale_fill_viridis(limits = c(min_v,max_v), breaks = scales::pretty_breaks(3)) +
    ggtitle(paste0(main,"\nSurface")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 10),
          axis.line = element_blank(),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.direction = "horizontal",
          legend.position = c(0.975,0.025),
          legend.justification = c(1,0),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          axis.text.x = element_blank(),
          axis.title.x = element_blank()) + labs(fill = "Mean Shannon Diversity") +
    scale_x_continuous(breaks = c(-125,-122,-119)) + 
    annotate(x = -126.95, y = 36.95, label = letters[i], size = 7, geom = "text", fontface = "bold")
    
  
  shannon_dcm <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = dcm, aes(x = long, y = lat, fill = shannon), color = "black", size =4, stroke = 0.1, shape = 21) +
    scale_fill_viridis(limits = c(min_v,max_v), breaks = scales::pretty_breaks(3)) +
    ggtitle(paste0("Deep Chlorophyll Maximum")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 10),
          axis.line = element_blank(),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.direction = "horizontal",
          legend.position = c(0.975,0.025),
          legend.justification = c(1,0),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8)) + labs(fill = "Mean Shannon Diversity") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  shannon <- list(p1 = shannon_surf, p2 = shannon_dcm)
  
  return(shannon)
  
}

plot_list <- list()


for (i in 1:length(in_group_list)) {
  
  plot_list[[i]] <- diveristy_figure_response(map_file = paste0("output/", in_group_list[i], "_map_S.Rdata"),
                                              full_dat = paste0("output/", in_group_list[i], "_full_data_S.Rdata"),
                                              figure_start = paste0("figures/diversity/", in_group_list[i], "_"),
                                              main = in_group_names[i])
  
}


out <- plot_list[[6]]$p1 + plot_list[[1]]$p1 + plot_list[[2]]$p1 +
  plot_list[[6]]$p2 + plot_list[[1]]$p2 + plot_list[[2]]$p2 +
  plot_list[[3]]$p1 + plot_list[[4]]$p1 + plot_list[[5]]$p1 +  
  plot_list[[3]]$p2 + plot_list[[4]]$p2 + plot_list[[5]]$p2 + plot_layout(nrow = 4)

pdf("figures_S/supp_fig_7.pdf", width = 12, height = 15)
plot(out)
dev.off()


