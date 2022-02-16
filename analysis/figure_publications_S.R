source("analysis/figure_functions_S.R")

###### Main Text Figures #####
library(ggsn)
library(cmocean)
#### Figure 1 Revised: Physical Conditions and ASV compare ####

# CHECK VENN DIAGRAM

fig_1_rev_func <- function(in_temp = "output/CALCOFI_temp_tables.Rdata",
                       in_cyano = "output/cyano_16s_map.Rdata",
                       in_full = "output/cyano_16s_full_data.Rdata",
                       in_venn = "output/venn_diagram_S.Rdata",
                       fig_name = "figures_S/fig_1_S.pdf",
                       tsize = 16, psize = 6){
  
  map <- map_data("world")    
  
  load(in_temp)
  load(in_cyano)
  load(in_full)
  load(in_venn)
  
  som_maps <- som_maps %>% filter(substr(Sta_ID,2,3) > 75)
  full_dat <- full_dat %>% filter(substr(Sta_ID,2,3) > 75)
  
  som_maps$line <- substr(som_maps$Sta_ID,1,5) %>% as.numeric()
  
  som_line <- som_maps %>% group_by(line) %>%
    summarise(x_coord = min(long), y_coord = min(lat)) %>%
    filter(!line %in% c(81.8,93.4))
  
  som_maps$cardinal <- "Cardinal Stations"
  som_maps$cardinal[which(som_maps$n_samps < 35)] <- NA
  
  stations <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "n_samps"),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    geom_text(data = som_line, aes(x = x_coord - 0.43, y = y_coord - 0.05, label = line),
              hjust = 1, size = 6) +
    geom_point(data = som_maps %>% filter(!is.na(cardinal)), aes(x = long, y = lat, shape = cardinal),
               show.legend = TRUE, pch = 0 , size = 6) + 
    annotate(x = -126.95, y = 36.99, label = "a", size = 7, geom = "text", fontface = "bold") +
    scale_fill_gradient(name = "# of Samples", low = "white", high = "red") +
    scalebar(x.min = -127, x.max = -116, 
             y.min =  28, y.max = 37, dist = 200, dist_unit = "km",
             transform = TRUE, location = "bottomright", st.dist = -0.05, height = 0.025, st.size = 4,  st.color = "black") +
    scale_x_continuous(breaks = c(-125,-122,-119)) + 
    scale_shape_manual(na.value = NA) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.99, 0.99),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          title = element_text(size = tsize),
          plot.margin = ggplot2::margin(5,5,5,5)) 

  
  depths <- ggplot() +
    geom_point(data = full_dat, aes(x = dist_to_coast, y = Depthm), size = psize/3, pch = 21, alpha = 0.5, fill = "grey60") +
    scale_y_reverse() + scale_x_continuous(breaks = c(0,200,400)) +
    annotate(x = 0, y = 0, label = "b", size = 7, geom = "text", fontface = "bold") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.margin = ggplot2::margin(5,5,5,5)) +
    labs(x = "Distance to Coast (km)", y = "Depth (m)")
  
  sst_mean <- ggplot() + 
    geom_tile(data = mean_table, aes(x = lon, y = lat, fill = Mean), width =0.26, height = 0.26) +
    geom_point(data = som_maps, aes_string(x = "long", y = "lat"),
               color = "grey90", fill = NA,size =psize, stroke = 1.5, shape = 21) +
    scale_fill_gradient2(name = "SST Mean (°C)", low = "darkblue", mid = "white", high = "darkred",
                         limits = c(15,18), oob = squish, midpoint = 16.5,
                         breaks = c(15,16,17,18), labels = c("<15", "16", "17", ">18")) +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    annotate(x = -126.95, y = 36.99, label = "d", size = 7, geom = "label", fontface = "bold") +
    xlab("Longitude") + ylab("Latitude") +
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.99, 0.99),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-2),
          plot.margin = ggplot2::margin(5,5,5,5)) 
  
  nc_depth <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "c", size = 7, geom = "text", fontface = "bold") +
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "NC_mean"),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean\nNitracline\nDepth (m)", low = "darkblue", high = "cyan", trans = 'reverse') +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.99, 0.99),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-2),
          plot.margin = ggplot2::margin(5,5,5,5)) 
  
  
  stations <- stations + theme(axis.title.x=element_blank(),
                               axis.text.x=element_blank())

  
  sst_mean <- sst_mean + theme(axis.title.y=element_blank(),
                               axis.text.y=element_blank())
  
  
  
  layout <- "ABEE
             CDEE"
  
  grobTree(venn_plot)
  
  patch <- stations + depths +
    nc_depth + sst_mean +
    wrap_elements(grobTree(venn_plot)) +
    plot_layout(design = layout) 
  
  patch <- patch + plot_annotation(tag_levels = list(c("", "", "", "", "e"))) & 
    theme(plot.tag = element_text(size = 20, hjust = 0, vjust = 0, face = "bold"))
  
  pdf(fig_name, width = 22, height = 12)
  print(patch)
  dev.off()
  
}

#### Figure 2: SOMs ####

fig_2_func <- function(in_list = plot_list, file_name = "figures_S/fig_2_S.pdf", tsize  = 17){
  
  in_list[[6]]$bp <- in_list[[6]]$bp + theme(axis.title.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             plot.title = element_text(hjust = 0.5, size = tsize),
                                             legend.position = "none",
                                             plot.background = element_blank()) + ggtitle("Archaea") +
    annotate(x = -126.95, y = 36.95, label = "a", size = 6.5, geom = "text", fontface = "bold") +
    annotate(x = -126.99, y = 28.2, label = "Offshore", size = 6, geom = "text", hjust = 0)
  in_list[[13]]$bp <- in_list[[13]]$bp + theme(axis.title.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0.5, size = tsize),
                                               legend.position = "none",
                                               plot.background = element_blank()) + ggtitle(" Bacteria") +
    annotate(x = -126.95, y = 36.95, label = "b", size = 6.5, geom = "text", fontface = "bold")
  in_list[[14]]$bp <- in_list[[14]]$bp + theme(axis.title.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0.5, size = tsize),
                                               legend.position = "none",
                                               plot.background = element_blank()) + ggtitle("Cyanobacteria")  +
    annotate(x = -126.95, y = 36.95, label = "c", size = 6.5, geom = "text", fontface = "bold")
  in_list[[16]]$bp <- in_list[[16]]$bp + theme(axis.title.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0.5, size = tsize),
                                               legend.position = "none",
                                               plot.background = element_blank()) + ggtitle("Photosynthetic eukaryotic\n Protists") +
    annotate(x = -126.95, y = 36.97, label = "d", size = 6.5, geom = "text", fontface = "bold")
  in_list[[15]]$bp <- in_list[[15]]$bp + theme(axis.title.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0.5, size = tsize),
                                               plot.tag = element_text(angle = -90, hjust = 0.25, size = tsize),
                                               plot.tag.position = c(1.05,0.5),
                                               plot.background = element_blank(),
                                               plot.margin = ggplot2::margin(1,25,1,1)) + ggtitle("Heterotrophic\neukaryotic protists") +
    labs(tag = "Frequency") + annotate(x = -126.95, y = 36.95, label = "e", size = 6.5, geom = "text", fontface = "bold")
    
  in_list[[6]]$rp <- in_list[[6]]$rp + theme(legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) +
    annotate(x = -126.99, y = 28.2, label = "Nearshore", size = 6, geom = "text", hjust = 0)
  in_list[[13]]$rp <- in_list[[13]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))

  in_list[[14]]$rp <- in_list[[14]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[16]]$rp <- in_list[[16]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[15]]$rp <- in_list[[15]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               plot.tag = element_text(angle = -90, hjust = 0.5, size = tsize),
                                               plot.tag.position = c(1.05,0.5)) + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) +
    labs(tag = "Frequency")
  
  
  
  pdf(file = file_name, width = 13, height = 6.5)
  print(in_list[[6]]$bp + in_list[[13]]$bp + in_list[[14]]$bp +
          in_list[[16]]$bp + in_list[[15]]$bp +
          in_list[[6]]$rp + in_list[[13]]$rp + in_list[[14]]$rp +
          in_list[[16]]$rp + in_list[[15]]$rp +
          plot_layout(ncol = 5))
  dev.off()
}

#### Figure 3: Reg and Var Import ####

fig_3_func <- function(file_name = "figures_S/fig_3_S.pdf", tsize = 15){
  
  cyano <- regression_figure(glm_file = "output/cyano_16s_glm_S.Rdata",
                             map_file = "output/cyano_16s_map_S.Rdata",   
                             figure_name = "figures/glm_plots/cyano_16s_som_",
                             main = "16s Cyanobacteria", cluster1 = "Nearshore", cluster2 = "Offshore",
                             var = "NC_mean", var_name = "Nitracline Depth (m)")
  
  cyano <- cyano + ggtitle("") + 
    theme(legend.justification=c(1,0.5), 
          legend.position=c(0.90, 0.5),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          plot.title = element_text(hjust = 0, size = tsize),
          plot.margin = ggplot2::margin(1,1,1,1)) +
    annotate(x = 0, y = 1, geom = "text", label = "a", size = 6.5, fontface = "bold")
  
  aic_plot <- full_aic_table_figure(in_group_list = c("archaea_16s","bacteria_m_euks_16s", "cyano_16s","euks_auto_18sv9",
                                                      "euks_hetero_18sv9"),
                                    in_group_names = c("Archaea","Bacteria", "Cyanobacteria", "Photosynthetic Eukaryotic\nProtists",
                                                       "Heterotrophic\nEukaryotic Protists"),
                                    minimum_tp = 4, width_plot = 10,
                                    figure_name_2 = paste0("figures/aic_figures/big_group_aic_plot_logit",".pdf"),
                                    title_name = "Variable Importance")
  
  aic_plot <- aic_plot + ggtitle("") + 
    theme(axis.title = element_blank(),
          axis.text.x = element_text(size = (tsize - 2.5)),
          plot.title = element_text(hjust = 0, size = tsize),
          plot.margin = ggplot2::margin(1,1,1,1)) +
    annotate(xmin = 0.4, ymin = 13.25, xmax = 0.6, ymax = 14.75, geom = "rect", fill = "white", color= "white") +
    annotate(x = 0.5, y = 14, geom = "text", label = "b", size = 6.5, fontface = "bold") 
    
  
  a <- cyano + aic_plot + plot_layout(widths = c(1,2)) 
  
  pdf(file = file_name, width = 16, height = 8)
  print(a)
  dev.off()
  
}

#### Figure 4: Diversity Example ####

fig_4_func <- function(file_name = "figures_S/fig_4_S.pdf",
                       in_group = c("total","diatom_18sv9"), basic = c("totals","18s_diatom"),
                       name = c("All ASVs", "Diatoms"), tsize = 18){
  
  
  
  all_map <- diveristy_figure(map_file = paste0("output/", in_group[1], "_map_S.Rdata"),
                          full_dat = paste0("output/", in_group[1], "_full_data_S.Rdata"),
                          figure_start = paste0("figures/diversity/", in_group[1], "_"),
                          main = name[1])
  
  all_alpha_gamma <- alpha_versus_gamma_figure(full_data_file = paste0("output/", in_group[1], "_full_data_S.Rdata"),
                                           raw_data_file = paste0("data/", basic[1], "_S.Rdata"),
                                           map_file = paste0("output/", in_group[1], "_map_S.Rdata"), minimum_tp = 8,
                                           figure_name = paste0("figures/diversity/", in_group[1], "_alpha_gamma_S.pdf"),
                                           main = name[1])

  
  map <- diveristy_figure(map_file = paste0("output/", in_group[2], "_map_S.Rdata"),
                          full_dat = paste0("output/", in_group[2], "_full_data_S.Rdata"),
                          figure_start = paste0("figures/diversity/", in_group[2], "_"),
                          main = name[2])
  
  alpha_gamma <- alpha_versus_gamma_figure(full_data_file = paste0("output/", in_group[2], "_full_data_S.Rdata"),
                                           raw_data_file = paste0("data/", basic[2], "_S.Rdata"),
                                           map_file = paste0("output/", in_group[2], "_map_S.Rdata"), minimum_tp = 8,
                                           figure_name = paste0("figures/diversity/", in_group[2], "_alpha_gamma_S.pdf"),
                                           main = name[2])
  
  
  all_map <- all_map + ggtitle("All ASVs") +
    annotate(x = -126.95, y = 36.99, label = "a", size = 7, geom = "text", fontface = "bold") +
    theme(plot.title = element_text(hjust = 0, size = tsize+3),
          legend.justification=c(1,1), 
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text.y  = element_text(size = tsize+3),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(size = (tsize-3)),
          axis.title.y = element_text(size = tsize+3),
          legend.title = element_text(size = (tsize-3)),
          plot.margin = ggplot2::margin(5,1,1,1)) + 
    labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) +
    scale_fill_viridis(breaks = c(5.6,5.2,4.8))
  
  all_alpha_gamma <- all_alpha_gamma + ggtitle("All ASVs") +
    annotate(x = 0, y = 7, label = "b", size = 7, geom = "text", fontface = "bold") +
    theme(plot.title = element_text(hjust = 0, size = tsize+3),
          legend.justification=c(1,0), 
          legend.position=c(0.95, 0.05),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text.y = element_text(size = tsize+3),
          legend.text = element_text(size = tsize),
          axis.title.y = element_text(size = tsize+3),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = ggplot2::margin(5,1,1,5)) + 
    labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=5))
  
  
  map <- map + ggtitle("Diatoms") +
    annotate(x = -126.95, y = 36.99, label = "c", size = 7, geom = "text", fontface = "bold") +
    theme(plot.title = element_text(hjust = 0, size = tsize+3),
          legend.justification=c(1,1), 
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize+3),
          legend.text = element_text(size = (tsize-3)),
          axis.title = element_text(size = tsize+3),
          legend.title = element_text(size = (tsize-3)),
          plot.margin = ggplot2::margin(5,1,1,1)) + 
    labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) +
    scale_fill_viridis()
  
  alpha_gamma <- alpha_gamma + ggtitle("Diatoms") +
    annotate(x = 0, y = 4.35, label = "d", size = 7, geom = "text", fontface = "bold") +
    theme(plot.title = element_text(hjust = 0, size = tsize+3),
          legend.justification=c(1,0.5), 
          legend.position=c(0.97, 0.46),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize+3),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize+3),
          plot.margin = ggplot2::margin(5,1,1,5)) + 
    labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=5))
  
  out_plot <- full_aic_table_figure_diversity_sign(in_group_list = c("archaea_16s",
                                                                     "bacteria_m_euks_16s",
                                                                     "cyano_16s","euks_auto_18sv9",
                                                                     "euks_hetero_18sv9"),
                                                   in_group_names = c("Archaea","Bacteria", 
                                                                      "Cyanobacteria",
                                                                      "Photosynthetic Eukaryotic\nProtists",
                                                                      "Heterotrophic\nEukaryotic Protists"),
                                                   figure_name_2 = NULL,
                                                   col = 25, width_plot = 12, tsize = 17)
  
  
    out_plot <- out_plot +
    annotate(xmin = 0.4, ymin = 13.25, xmax = 0.6, ymax = 14.75, geom = "rect", fill = "white", color= "white") +
    annotate(x = 0.5, y = 14, geom = "text", label = "e", size = 7, fontface = "bold") 
  
  layout <- "ABEE
             CDEE"
  
  plots <- all_map + all_alpha_gamma + map +
    alpha_gamma + out_plot + 
    plot_layout(design = layout)
    
  pdf(file = file_name, width = 26, height = 13)
  print(plots)
  dev.off()
  
}


###### Figure 5: Time Nutrients ect. #####

fig_5_func <- function(in_upwell = "output/upwelling_plots.Rdata",
                       nutrient_dat = "output/mld_mean_profiles_S.Rdata",
                       time_plot = "figures_S/fig_5_S.pdf",
                       tsize = 15, psize = 6, in_list = fig_list){
  
  map <- map_data("world") 
  
  load(in_upwell)
  load(nutrient_dat)
  
  yearly_cuti <- yearly_cuti + theme(plot.title = element_text(size = tsize),
                                     axis.text = element_text(size = tsize),
                                     axis.title = element_text(size = tsize),
                                     legend.text = element_text(size = tsize),
                                     legend.title = element_text(size = tsize),
                                     legend.position = "none") + 
    ggtitle("Coastal Upwelling\nTransport Index (CUTI)") +
    xlab("Year Day") + ylab("CUTI") +
    annotate(x = 1, y = 1, geom = "text", size = 7, fontface = "bold", label = "a")
  
  yearly_beuti <- yearly_beuti + theme(plot.title = element_text(size = tsize),
                                       axis.text = element_text(size = tsize),
                                       axis.title = element_text(size = tsize),
                                       legend.text = element_text(size = tsize),
                                       legend.title = element_text(size = tsize),
                                       legend.position = "none") +
    ggtitle("Biologically Effective Upwelling\nTransport Index (BEUTI)") +
    xlab("Year Day") + ylab("BEUTI") +
    annotate(x = 1, y = 8.1, geom = "text", size = 7, fontface = "bold", label = "b")
  
  
  yearly_nitrate <- yearly_nitrate + theme(plot.title = element_text(size = tsize),
                                           axis.text = element_text(size = tsize),
                                           axis.title = element_text(size = tsize),
                                           legend.text = element_text(size = tsize),
                                           legend.title = element_text(size = tsize),
                                           legend.key = element_blank(),
                                           legend.position = c(0.99,0.99),
                                           legend.justification = c(1,1)) +
    ggtitle("Regionally available\n nitrate (BEUTI/CUTI)") +
    ylab("BEUTI/CUTI") +
    annotate(x = 1, y = 8.3, geom = "text", size = 7, fontface = "bold", label = "c")
  
  early_phase <- early_phase %>% filter(Sta_ID %in% late_phase$Sta_ID)
  
  early_phase$n_diff <- early_phase$mean_no3 - late_phase$mean_no3
  early_phase$p_diff <- early_phase$mean_po4 - late_phase$mean_po4
  
  
  early_n <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_no3),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(paste("Mean ",NO[3]," ",mu,M)), low = "white", high = "red",
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.position=c(0.97,0.03),
          legend.direction = "horizontal",
          legend.justification = c(1,0),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize-2),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-2),
          plot.margin = ggplot2::margin(2,2,0,2)) + 
    ggtitle("Mean nitrate μM (2014-2016)") +
    scale_x_continuous(breaks = c(-125,-122,-119)) + 
    annotate(x = -126.95, y = 36.99, label = "d", size = 7, geom = "text", fontface = "bold") 
  
  late_n <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_no3),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(paste("Mean ",NO[3]," ",mu,M)), low = "white", high = "red",
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.position=c(0.97,0.03),
          legend.direction = "horizontal",
          legend.justification = c(1,0),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize-2),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-2),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = ggplot2::margin(2,2,0,2)) +
    ggtitle("Mean nitrate μM (2017-2018)") +
    scale_x_continuous(breaks = c(-125,-122,-119)) + 
    annotate(x = -126.95, y = 36.99, label = "e", size = 7, geom = "text", fontface = "bold")
  
  diff_n <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = n_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(paste(Delta~NO[3]," ", mu,M)),
                        low = "blue", high = "white",
                        limits = c(-6,0), oob = scales::squish,
                        breaks = c(-6,-4,-2,0)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.position=c(0.97,0.03),
          legend.direction = "horizontal",
          legend.justification = c(1,0),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize-2),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-2),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = ggplot2::margin(2,2,0,2)) + 
    ggtitle("Changes in mean nitrate\n(2014-2016) - (2017-2018)") +
    scale_x_continuous(breaks = c(-125,-122,-119)) + 
    annotate(x = -126.95, y = 36.99, label = "f", size = 7, geom = "text", fontface = "bold")
  
  low <- 0.2
  high <- 0.6
  
  arch <- (in_list[[6]]$surf + ggtitle("Archaea") +
      scale_fill_gradient(limits = c(low, high),
                          low = "white", high = "red", oob = scales::squish) +
      labs(fill = "Bray-Curtis\nSimilarity")+
        theme(legend.position = "none",
              plot.margin = ggplot2::margin(1,2,2,2))) + 
    annotate(x = -126.95, y = 36.99, label = "g", size = 7, geom = "text", fontface = "bold") +
    scale_x_continuous(breaks = c(-125,-122,-119))
    
    bact <- (in_list[[13]]$surf + ggtitle("Bacteria") +
       scale_fill_gradient(limits = c(low,high),
                           low = "white", high = "red", oob = scales::squish) +
       labs(fill = "Bray-Curtis\nSimilarity") +
       theme(axis.title.y = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank())+
         theme(legend.position = "none",
               plot.margin = ggplot2::margin(1,2,2,2))) + 
      annotate(x = -126.95, y = 36.99, label = "h", size = 7, geom = "text", fontface = "bold") +
      scale_x_continuous(breaks = c(-125,-122,-119))
    
    cyano <- (in_list[[14]]$surf + ggtitle("Cyanobacteria") +
       scale_fill_gradient(limits = c(low,high),
                           low = "white", high = "red", oob = scales::squish) +
       labs(fill = "Bray-Curtis\nSimilarity") +
       theme(axis.title.y = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank())+
         theme(legend.position = "none",
               plot.margin = ggplot2::margin(1,2,2,2))) + 
      annotate(x = -126.95, y = 36.99, label = "i", size = 7, geom = "text", fontface = "bold") +
      scale_x_continuous(breaks = c(-125,-122,-119))
    
    photo <- (in_list[[16]]$surf + ggtitle("Photosynthetic Eukaryotic\nProtists") +
       scale_fill_gradient(limits = c(low,high),
                           low = "white", high = "red", oob = scales::squish) +
       labs(fill = "Bray-Curtis\nSimilarity") +
       theme(axis.title.y = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank())+
         theme(legend.position = "none",
               plot.margin = ggplot2::margin(1,2,2,2))) + 
      annotate(x = -126.95, y = 36.99, label = "j", size = 7, geom = "text", fontface = "bold") +
      scale_x_continuous(breaks = c(-125,-122,-119))
    
    hetero <- (in_list[[15]]$surf + ggtitle("Heterotrophic Eukaryotic\nProtists") +
       scale_fill_gradient(limits = c(low,high),
                           low = "white", high = "red", oob = scales::squish,
                           breaks = c(0.2,0.4,0.6)) +
       labs(fill = "Bray-Curtis\nSimilarity") +
       theme(axis.title.y = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),
             legend.position=c(0.97,0.03),
             legend.direction = "horizontal",
             legend.justification = c(1,0),
             legend.text = element_text(size = tsize -4),
             legend.title = element_text(size = tsize - 4),
             plot.margin = ggplot2::margin(1,2,2,2),
             legend.margin = ggplot2::margin(2,15,2,2)))+ 
      annotate(x = -126.95, y = 36.99, label = "k", size = 7, geom = "text", fontface = "bold") +
      scale_x_continuous(breaks = c(-125,-122,-119))
    
    fig_plot <- (arch + bact + cyano + photo + hetero + plot_layout(nrow = 1)) 
    
  layout = "ABC
            DEF
            GGG"
  
  plot1 <- yearly_cuti + yearly_beuti + yearly_nitrate  +
    early_n + late_n + diff_n  + fig_plot +
    plot_layout(design = layout, heights = c(1,1,0.8))

  pdf(file = time_plot, width = 14, height = 14)
  print(plot1)
  dev.off()

}


##### Figure 6: Community vs Time #####

fig_6_func <- function(in_phyto = "output/euks_auto_18sv9_diffs_S.Rdata",
                       in_euks = "output/euks_hetero_18sv9_diffs_S.Rdata",
                       in_cyano = "output/cyano_16s_diffs_S.Rdata",
                       in_bact = "output/bacteria_m_euks_16s_diffs_S.Rdata",
                       in_arch = "output/archaea_16s_diffs_S.Rdata",
                       ncd_data = "data/ncd_data.csv",
                       cyano_dat = "output/cyano_16s_full_data_S.Rdata",
                       gradient_plot_file = "figures_S/fig_6_S.pdf",
                       gradient_plot_file2 = "figures_S/fig_6_S_schematic.pdf",
                       tsize = 15){
  
  
  load(in_phyto)
  phyto_gradient <- phase
  
  load(in_euks)
  euk_gradient <- phase
  
  load(in_cyano)
  cyano_gradient <- phase
  
  load(in_bact)
  bact_gradient <- phase
  
  load(in_arch)
  arch_gradient <- phase
  
  load(cyano_dat)
  meta <- read.csv("data/NCOG_sample_log_DNA_stvx_meta_2014-2020_prim_prod.csv")
  full_dat$sample_match <- gsub("_S", "", full_dat$Sample.Name)
  
  full_dat$IntC14 <- meta$IntC14[match(full_dat$sample_match, meta$Sample.Name)]
  full_dat$IntC14_day <- full_dat$IntC14*2
  full_dat$ef_ratio <- ((0.5857-(0.0165*full_dat$T_degC))*full_dat$IntC14_day)/(51.7 + full_dat$IntC14_day)
  
  ncd_df <- read.csv(ncd_data)
  
  sum <- ncd_df %>% filter(CruiseAlias == 201907)
  wint <- ncd_df %>% filter(CruiseAlias == 201501)
  
  full_dat$som_id[which(full_dat$som_id == 1)] <- "Offshore"
  full_dat$som_id[which(full_dat$som_id == 2)] <- "Nearshore"
  
  sum_cyano <- full_dat %>% filter(Cruise == 201907)
  wint_cyano <- full_dat %>% filter(Cruise == 201501)
  
  sum_pie <- sum_cyano %>% group_by(som_id) %>%
    summarise(count = n()) %>%
    mutate(prop = count / sum(count) *100) %>%
    mutate(ypos = cumsum(prop)- 0.4*prop ) 
  
  wint_pie <- wint_cyano %>% group_by(som_id) %>%
    summarise(count = n()) %>%
    mutate(prop = count / sum(count) *100) 
  wint_pie$ypos <- c(95,45)
  
  mean(sum_cyano$ef_ratio, na.rm = TRUE)
  mean(wint_cyano$ef_ratio, na.rm = TRUE)

  sub_a <- ggplot(wint_pie, aes(x="", y=prop, fill=som_id)) +
    geom_bar(stat="identity", width=1, color=NA) +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="none") +
    geom_text(aes(y = ypos, label = count), color = "white", size=4, fontface = "bold") +
    scale_fill_manual(values = c("blue","red"))
  
  sub_b <- ggplot(sum_pie, aes(x="", y=prop, fill=som_id)) +
    geom_bar(stat="identity", width=1, color=NA) +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.position="none") +
    geom_text(aes(y = ypos, label = count), color = "white", size=4, fontface = "bold") +
    scale_fill_manual(values = c("blue","red"))
  
  gradient <- data.frame(x1 = rep(-10:380,times(5000)), y1 = rep(seq(0,150, length = 5000), each = 391))
  
  a <- ggplot() +
    # geom_line(data = gradient, aes(x = x1, y = y1, color = y1, group = y1), show.legend = FALSE) + 
    geom_smooth(data = wint, aes(x = abs(Distance.from.Shore), y = NCDepth),
                method = "glm", color = "black") +
    geom_point(data = wint_cyano, aes(x = abs(Distance), y = Depthm, fill = som_id),
               pch =21, size = 4) +
    scale_y_reverse(limits = c(140,0)) + scale_x_reverse(limits = c(369,0)) +
    # geom_hline(yintercept = 0) +
    labs(x = "Distance from Coast (km)", y = "Depth (m)") +
    scale_fill_manual(values = c("blue","red")) +
    # scale_color_gradient(high = "dodgerblue3", low = "steelblue1") +
    # annotate(geom = "text", x = 184.5, y = -6, label = "surface", size = 5) +
    annotate(geom = "text", x = 50, y = 135, label = "Winter 2015", size = 6) +
    annotate(geom = "text", x = 368, y = 0, label = "a", size = 6, fontface= "bold") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),
          legend.position = c(0.95,0.12),
          legend.justification = c(1,0),
          legend.key = element_blank(),
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.background = element_blank()) 
  
  a <- a + inset_element(sub_a, left = 0.35, right = 0.75, bottom = 0.05, top = 0.3)
  
  b <- ggplot() +
    # geom_line(data = gradient, aes(x = x1, y = y1, color = y1, group = y1), show.legend = FALSE) + 
    geom_smooth(data = sum, aes(x = abs(Distance.from.Shore), y = NCDepth),
                method = "glm", color = "black") +
    geom_point(data = sum_cyano, aes(x = abs(Distance), y = Depthm, fill = som_id),
               pch =21, size = 4) +
    scale_y_reverse(limits = c(140,0)) + scale_x_reverse(limits = c(369,0)) +
    # geom_hline(yintercept = 0) +
    labs(x = "Distance from Coast (km)", y = "Depth (m)") +
    scale_fill_manual(values = c("blue","red")) +
    # scale_color_gradient(high = "dodgerblue3", low = "steelblue1") +
    # annotate(geom = "text", x = 184.5, y = -6, label = "surface", size = 5) +
    annotate(geom = "text", x = 50, y = 135, label = "Summer 2019", size = 6) +
    annotate(geom = "text", x = 368, y = 0, label = "b", size = 6, fontface= "bold") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),
          legend.position = c(0.99,0.15),
          legend.justification = c(1,0),
          legend.key = element_blank(),
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.background = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) 

 b <- b + inset_element(sub_b, left = 0.35, right = 0.75, bottom = 0.05, top = 0.3)
 
 gradient <- data.frame(x1 = rep(1:10,times(10000)), y1 = rep(seq(0,1024, length = 10000), each = 10))
 floor <- data.frame(x1 = c(seq(1,10, by = 0.01),10.2,1), y1 = c(2^(seq(1,10, by = 0.01)),-10,-10))
 floor_line <- data.frame(x1 = seq(1,10, by = 0.01), y1 = 2^(seq(1,10, by = 0.01)))
 ceiling <- data.frame(x1 = seq(1,10, by = 0.01), y1 = rep(1024, times = 901))
 out_edge <- data.frame(x1 = rep(1, times = 2), y1 = c(0,1024))
 
 s1 <- data.frame(x1 = seq(1.4,8, length = 901), y1 = seq(500,1024,length = 901))
 a1 <- data.frame(x1 = seq(1,3.2, length = 200), y1 = rep(500,times = 200))
 diatoms <- data.frame(x1 = rep(8.5, times = 10), y1 = rep(850, times = 10))
 cyano <- data.frame(x1 = rep(3.8, times = 50), y1 = rep(880, times = 50))
 cyano2 <- data.frame(x1 = rep(5, times = 100), y1 = rep(945, times = 100))
 s2 <- data.frame(x1 = seq(1.4,9.895, length = 901), y1 = seq(750,950,length = 901))
 a2 <- data.frame(x1 = seq(1,3.2, length = 200), y1 = rep(750,times = 200))
 
 set.seed(765)
 
 a_sch <- ggplot() +
   geom_line(data = gradient, aes(x = x1, y = y1, color = y1, group = y1), show.legend = FALSE) + 
   geom_polygon(data = floor, aes(x = x1, y = y1), lwd = 1, fill = "white") +
   geom_line(data = floor_line, aes(x = x1, y = y1), lwd = 1) +
   geom_line(data = ceiling, aes(x = x1, y = y1), lwd = 1) +
   geom_line(data = out_edge, aes(x = x1, y = y1), lwd = 1) + 
   geom_line(data = s2, aes(x = x1, y = y1), lwd = 1) +
   geom_line(data = a2, aes(x = x1, y = y1), lty = 2, lwd = 1) +
   geom_curve(aes(yend = 776, y = 750, xend = 2.5, x = 2.52), lwd = 1) +
   geom_jitter(data = cyano2, aes(x = x1, y = y1), width = 2, height = 50,
               size = 2, pch = 21, fill = "chartreuse3", stroke = 0.1) +
   annotate(geom = "text", x = 4.2, y = 760, label = expression(alpha*" = slope"), size = 6) +
   annotate(geom = "text", x = 2, y = 880, label = "N-", size = 8) +
   annotate(geom = "text", x = 2, y = 220, label = "N+", size = 8) +
   annotate(geom = "text", x = 5, y = 1100, label = "Oligotrophic", size = 6) +
   annotate(geom = "text", x = 1, y = 1100, label = "a", size = 6, fontface = "bold") +
   scale_color_gradient(low = "dodgerblue3", high = "steelblue1") +
   theme_void()
 
 set.seed(65)
 
 b_sch <- ggplot() +
   geom_line(data = gradient, aes(x = x1, y = y1, color = y1, group = y1), show.legend = FALSE) + 
   geom_polygon(data = floor, aes(x = x1, y = y1), lwd = 1, fill = "white") +
   geom_line(data = floor_line, aes(x = x1, y = y1), lwd = 1) +
   geom_line(data = ceiling, aes(x = x1, y = y1), lwd = 1) +
   geom_line(data = out_edge, aes(x = x1, y = y1), lwd = 1) + 
   geom_line(data = s1, aes(x = x1, y = y1), lwd = 1) +
   geom_line(data = a1, aes(x = x1, y = y1), lty = 2, lwd = 1) +
   geom_curve(aes(yend = 587.3333, y = 500, xend = 2.5, x = 2.6), lwd = 1) +
   geom_jitter(data = diatoms, aes(x = x1, y = y1), width = 1, height = 100,
               size = 5, pch = 23, fill = "chartreuse4") +
   geom_jitter(data = cyano, aes(x = x1, y = y1), width = 1, height = 100,
               size = 2, pch = 21, fill = "chartreuse3", stroke = 0.1) +
   annotate(geom = "text", x = 4, y = 550, label = expression(alpha*" = slope"), size = 6) +
   annotate(geom = "text", x = 2, y = 880, label = "N-", size = 8) +
   annotate(geom = "text", x = 2, y = 220, label = "N+", size = 8) +
   annotate(geom = "text", x = 3.3, y = 1100, label = "Oligotrophic", size = 6) +
   annotate(geom = "text", x = 8.8, y = 1100, label = "Eutrophic", size = 6) +
   annotate(geom = "text", x = 1, y = 1100, label = "b", size = 6, fontface = "bold") +
   geom_point(aes(x = 4, y = 350), size = 5, pch = 23, fill = "chartreuse4") +
   geom_point(aes(x = 4, y = 250), size = 2, pch = 21, stroke = 0.1, fill = "chartreuse3") +
   annotate(geom = "text", x = 5.6, y = 350, label = "= Copiotroph", size = 6) +
   annotate(geom = "text", x = 5.5, y = 250, label = "= Oligotroph", size = 6) +
   scale_color_gradient(low = "dodgerblue3", high = "steelblue1") +
   theme_void()
  
  phyto_gradient <- phyto_gradient + theme(legend.position = "none") +  
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    theme(axis.text.x  = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.title.y = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          legend.position = "none") + ggtitle("Photosynthetic eukaryotic protists") +
    annotate(x = 0.17, y = 0.64, geom = "text", label = "c", fontface = "bold", size = 6)
  
  euk_gradient <- euk_gradient + theme(legend.position = "none") + 
    theme(axis.text.x  = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_blank(),
          legend.position = "none") + ggtitle("Heterotrophic eukaryotic protists") +
    annotate(x = 0.17, y = 0.77, geom = "text", label = "d", fontface = "bold", size = 6)
  
  cyano_gradient <- cyano_gradient + theme(legend.position = "none")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          legend.key = element_blank()) + ggtitle("Cyanobacteria") +
    annotate(x = 0.17, y = 0.675, geom = "text", label = "e", fontface = "bold", size = 6)
  
  bact_gradient <- bact_gradient +
    theme(axis.title.y = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          legend.key = element_blank()) + ggtitle("Bacteria") +
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    annotate(x = 0.17, y = 0.71, geom = "text", label = "f", fontface = "bold", size = 6)
  
  arch_gradient <- arch_gradient +
    theme(axis.title.y = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize),
          legend.position = "none") + ggtitle("Archaea") +
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    annotate(x = 0.17, y = 0.7, geom = "text", label = "g", fontface = "bold", size = 6)
  
  grad_plot <- (a + b) / (phyto_gradient + euk_gradient + guide_area() +
    cyano_gradient + bact_gradient + arch_gradient + plot_layout(guides = "collect")) +
    plot_layout(heights = c(0.5,1))
  
  grad_plot2 <- (a_sch + b_sch) / (phyto_gradient + euk_gradient + guide_area() +
                            cyano_gradient + bact_gradient + arch_gradient + plot_layout(guides = "collect")) +
    plot_layout(heights = c(0.5,1))
  
  pdf(file = gradient_plot_file, width = 14, height = 14)
  print(grad_plot)
  dev.off()
  
  pdf(file = gradient_plot_file2, width = 14, height = 14)
  print(grad_plot2)
  dev.off()
  
}

alt_fig_6_func <- function(in_phyto = "output/euks_auto_18sv9_diffs_S.Rdata",
                       in_euks = "output/euks_hetero_18sv9_diffs_S.Rdata",
                       in_cyano = "output/cyano_16s_diffs_S.Rdata",
                       in_bact = "output/bacteria_m_euks_16s_diffs_S.Rdata",
                       in_arch = "output/archaea_16s_diffs_S.Rdata",
                       gradient_plot_file = "figures_S/fig_6_S_alt.pdf",
                       tsize = 12){
  
  
  load(in_phyto)
  phyto_gradient <- beuti_plot
  
  load(in_euks)
  euk_gradient <- beuti_plot
  
  load(in_cyano)
  cyano_gradient <- beuti_plot
  
  load(in_bact)
  bact_gradient <- beuti_plot
  
  load(in_arch)
  arch_gradient <- beuti_plot
  
  phyto_gradient <- phyto_gradient + theme(legend.position = "none") + xlab("")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.title.y = element_text(size = tsize),
          axis.text.y = element_text(size = tsize)) + ggtitle("Photosynthetic eukaryotic protists")
  
  euk_gradient <- euk_gradient + theme(legend.position = "none") + xlab("") + ylab("") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize)) + ggtitle("Heterotrophic eukaryotic protists")
  
  cyano_gradient <- cyano_gradient + theme(legend.position = "none")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") + 
    xlab("Biologically Effective Upwelling Transport Index\n(BEUTI)") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("Cyanobacteria")
  
  bact_gradient <- bact_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("Bacteria") +
    ylab("") + xlab("Biologically Effective Upwelling Transport Index\n(BEUTI)")
  
  arch_gradient <- arch_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("Archaea") +
    ylab("") + xlab("Biologically Effective Upwelling Transport Index\n(BEUTI)")
  
  grad_plot <- phyto_gradient + euk_gradient + guide_area() +
    cyano_gradient + bact_gradient + arch_gradient + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 16))
  
  pdf(file = gradient_plot_file, width = 12, height = 7)
  print(grad_plot)
  dev.off()
  
  
}

##### Supplementary Figures #####

##### Suppl Figure 1: Physical Maps #####

suppl_fig_1_func <- function(in_dat = "output/mld_mean_profiles_S.Rdata",
                             out_plot = "figures_S/supp_fig_1_S.pdf",
                             tsize = 18, psize = 6){
  
  load(in_dat)
  
  map <- map_data("world")
  
  # Physical 
  
  t_m <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat, aes(x = lon, y = lat, fill = mean_temp),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = "°C", low = "blue", high = "red", mid = "white", midpoint = 16.5) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text.y = element_text(size = tsize),
          axis.text.x = element_blank(),
          legend.text = element_text(size = tsize-4),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-4),
          axis.title.x = element_blank()) +
    ggtitle("Mean temperature (°C)") + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "a", size = 7, geom = "text", fontface = "bold")
  
  t_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_temp), aes(x = lon, y = lat, fill = coeff_temp),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text.y = element_text(size = tsize),
          axis.text.x = element_blank(),
          legend.text = element_text(size = tsize-4),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-4),
          axis.title.x = element_blank()) +
    ggtitle("Coeff. var. temperature")  + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "e", size = 7, geom = "text", fontface = "bold")
  
  sal_m <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat, aes(x = lon, y = lat, fill = mean_sal),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = "PSU", low = "blue", high = "gold2", mid = "white", midpoint = 33.35) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_blank(),
          legend.text = element_text(size = tsize-4),
          axis.title = element_blank(),
          legend.title = element_text(size = tsize-4)) +
    ggtitle("Mean salinity")  + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "b", size = 7, geom = "text", fontface = "bold")
  
  sal_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_sal), aes(x = lon, y = lat, fill = coeff_sal),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_blank(),
          legend.text = element_text(size = tsize-4),
          axis.title = element_blank(),
          legend.title = element_text(size = tsize-4)) +
    ggtitle("Coeff. var. salinity")  + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "f", size = 7, geom = "text", fontface = "bold")
  
  # Nutrients
  
  no3_m <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat, aes(x = lon, y = lat, fill = mean_no3),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(paste(NO[3]," ",mu,M)), low = "white", high = "purple") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_blank(),
          legend.text = element_text(size = tsize-4),
          axis.title = element_blank(),
          legend.title = element_text(size = tsize-4)) +
    ggtitle(expression(paste("Mean ", NO[3]))) + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "j", size = 7, geom = "text", fontface = "bold")
  
  no3_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_no3), aes(x = lon, y = lat, fill = coeff_no3),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text.x = element_text(size = tsize),
          axis.text.y = element_blank(),
          legend.text = element_text(size = tsize-4),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = tsize),
          legend.title = element_text(size = tsize-4)) +
    ggtitle(expression(paste("Coeff. var. ", NO[3]))) + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "n", size = 7, geom = "text", fontface = "bold")
  
  po4_m <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat, aes(x = lon, y = lat, fill = mean_po4),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(paste(PO[4]," ",mu,M)), low = "white", high = "purple") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_blank(),
          legend.text = element_text(size = tsize-4),
          axis.title = element_blank(),
          legend.title = element_text(size = tsize-4)) +
    ggtitle(expression(paste("Mean ", PO[4]))) + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "k", size = 7, geom = "text", fontface = "bold")
  
  po4_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_po4), aes(x = lon, y = lat, fill = coeff_po4),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text.x = element_text(size = tsize),
          axis.text.y = element_blank(),
          legend.text = element_text(size = tsize-4),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = tsize),
          legend.title = element_text(size = tsize-4)) +
    ggtitle(expression(paste("Coeff. var. ", PO[4]))) + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "o", size = 7, geom = "text", fontface = "bold")
  
  sio4_m <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat, aes(x = lon, y = lat, fill = mean_sio4),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name =  expression(paste(SiO[4]," ",mu,M)), low = "white", high = "purple") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_blank(),
          legend.text = element_text(size = tsize-4),
          axis.title = element_blank(),
          legend.title = element_text(size = tsize-4)) +
    ggtitle(expression(paste("Mean ", SiO[4]))) + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "l", size = 7, geom = "text", fontface = "bold")
  
  sio4_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_sio4), aes(x = lon, y = lat, fill = coeff_sio4),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text.x = element_text(size = tsize),
          axis.text.y = element_blank(),
          legend.text = element_text(size = tsize-4),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = tsize),
          legend.title = element_text(size = tsize-4)) +
    ggtitle(expression(paste("Coeff. var. ", SiO[4]))) + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "p", size = 7, geom = "text", fontface = "bold")
  
  # Chlorophyll
  
  chl_m <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat, aes(x = lon, y = lat, fill = mean_chl),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(mu*g/L), low = "white", high = "darkgreen") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text.y = element_text(size = tsize),
          axis.text.x = element_blank(),
          legend.text = element_text(size = tsize-4),
          axis.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          legend.title = element_text(size = tsize-4)) +
    ggtitle("Mean chl-a") + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "i", size = 7, geom = "text", fontface = "bold")
  
  chl_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_chl), aes(x = lon, y = lat, fill = coeff_chl),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize-4),
          axis.title.y = element_text(size = tsize),
          axis.title.x = element_text(size = tsize),
          legend.title = element_text(size = tsize-4)) +
    ggtitle("Coeff. var. chl-a") + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "m", size = 7, geom = "text", fontface = "bold")
  
  # MLD and NCD
  
  mld_m <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat, aes(x = lon, y = lat, fill = mean_mld),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Depth (m)", low = "darkblue", high = "cyan", trans = 'reverse') +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_blank(),
          legend.text = element_text(size = tsize-4),
          axis.title = element_blank(),
          legend.title = element_text(size = tsize-4)) +
    ggtitle("Mean mixed layer depth") + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "c", size = 7, geom = "text", fontface = "bold")
  
  mld_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_mld), aes(x = lon, y = lat, fill = coeff_mld),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_blank(),
          legend.text = element_text(size = tsize-4),
          axis.title = element_blank(),
          legend.title = element_text(size = tsize-4)) +
    ggtitle("Coeff. var. mixed layer depth") + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "g", size = 7, geom = "text", fontface = "bold")
  
  ncd_m <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat, aes(x = lon, y = lat, fill = mean_ncd),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Depth (m)", low = "darkblue", high = "cyan", trans = 'reverse') +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_blank(),
          legend.text = element_text(size = tsize-4),
          axis.title = element_blank(),
          legend.title = element_text(size = tsize-4)) +
    ggtitle("Mean nitracline depth") + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "d", size = 7, geom = "text", fontface = "bold")
  
  ncd_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_ncd), aes(x = lon, y = lat, fill = coeff_ncd),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_blank(),
          legend.text = element_text(size = tsize-4),
          axis.title = element_blank(),
          legend.title = element_text(size = tsize-4)) +
    ggtitle("Coeff. var. nitracline depth") + 
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.99, label = "h", size = 7, geom = "text", fontface = "bold")
  
  
  plot1 <- t_m + sal_m + mld_m + ncd_m + 
    t_c + sal_c + mld_c + ncd_c + 
    chl_m + no3_m + po4_m + sio4_m +
    chl_c + no3_c + po4_c + sio4_c + plot_layout(ncol = 4)
  
  pdf(file = out_plot, height = 18, width = 18)
  print(plot1)
  dev.off()
  
  }


##### Suppl Figure 2: NCOG, TARA, Polar Endemics #######

suppl_fig_2_func <- function(in_ncog = "data/all_18Sv9_rare.Rdata",
                             in_tara = "data/18sv9_tara_oceans.Rdata",
                             in_polar = "data/18sv9_tara_polar.Rdata",
                             output = "figures_S/supp_fig_2_S.pdf",
                             output_table = "figures_S/supp_table_2_S.pdf"){
  
  # load data
  load(in_ncog)
  load(in_tara)
  load(in_polar)
  
  eight_rare <- eight_rare[-which(as.numeric(substr(rownames(eight_rare),10,11)) < 76),]
  
  calcofi <- names(which(colSums(eight_rare) != 0))
  tara <- names(which(colSums(tara_dat) != 0))
  polar <- names(which(colSums(polar_dat) != 0))
  
  # tax
  
  eight_tax_id <- eight_tax_id[,1:3]
  colnames(eight_tax_id) <- c("ASV", "Taxon", "Confidence")
  
  colnames(polar_tax) <- c("ASV", "Taxon", "Confidence")
  colnames(tara_tax) <- c("ASV", "Taxon", "Confidence")
  
  tax_full <- bind_rows(eight_tax_id, polar_tax, tara_tax)
  
  tax_full$CalCOFI <- 0
  tax_full$Tara <- 0
  tax_full$Polar <- 0
  
  tax_full <- tax_full[!duplicated(tax_full$ASV),]
  
  tax_full$CalCOFI[which(!is.na(match(tax_full$ASV,calcofi)))] <- 1
  tax_full$Tara[which(!is.na(match(tax_full$ASV,tara)))] <- 1
  tax_full$Polar[which(!is.na(match(tax_full$ASV,polar)))] <- 1
  
  split_taxa <- separate(tax_full, Taxon, sep = ";", into = c("A","B","C", "D", "E", "F", "G", "H", "I"))
  
  split_taxa$endemic <- NA
  
  split_taxa$endemic[which(split_taxa$CalCOFI == 1 & split_taxa$Tara == 0 & split_taxa$Polar == 0)] <- "NCOG"
  split_taxa$endemic[which(split_taxa$CalCOFI == 0 & split_taxa$Tara == 1 & split_taxa$Polar == 0)] <- "TARA"
  split_taxa$endemic[which(split_taxa$CalCOFI == 0 & split_taxa$Tara == 0 & split_taxa$Polar == 1)] <- "TARA Polar"
  
  split_taxa$C[which(is.na(split_taxa$C))] <- "Other"
  
  split_taxa$C[which(split_taxa$C == "Alveolata_X")] <- "Alveolata"
  split_taxa$C[which(split_taxa$C == "Eukaryota_XX")] <- "Eukaryota"
  split_taxa$C[which(split_taxa$C == "Opisthokonta_X")] <- "Opisthokonta"
  split_taxa$C[which(split_taxa$C == "Stramenopiles_X")] <- "Stramenopiles"
  
  mean_groups <- split_taxa %>%
    filter(!is.na(endemic), !is.na(C)) %>%
    group_by(endemic,C) %>% summarise(count = n()) %>% mutate(prop = count/sum(count))
  
  order <- split_taxa %>%
    filter(!is.na(endemic), !is.na(C)) %>%
    group_by(C) %>% summarise(count = n()) %>% mutate(prop = count/sum(count)) %>% arrange(prop)
  
  mean_groups$C <- factor(mean_groups$C, levels = order$C)
  
  proportion_plot <- ggplot(mean_groups, aes(x = endemic, y = prop, fill = C)) + 
    geom_bar(stat = "identity") + theme_classic() +
    labs(x = "", y = "Proportion", fill = "Division") +
    scale_fill_manual(values = c("#87bb37","#5362da","#4ec253","#aa60d9",
                                 "#4e8d2b","#d74ebb","#4dc484","#d73d84",
                                 "#4f985a","#a73e99","#c0b63a","#7655b2",
                                 "#dd9031","#6b83e5","#858c2f","#e473bf",
                                 "#92bb6f","#dc385a","#4cc8c6","#d04f2a",
                                 "#5fa2da","#a27f2c","#bc93e1","#2f7040",
                                 "#8c5396","#439f84","#d7685f","#5a6bac",
                                 "#cbad6b","#9c426a","#696d2e","#de95c6",
                                 "#9a5b2c","#b06b9c","#e6936d","#a5444b", "grey70")) +
    theme(panel.border = element_rect(fill = NA, color = "black"))
  
  pdf(output, width = 10, height = 8)
  print(proportion_plot)
  dev.off()

  split_taxa$endemic[which(split_taxa$CalCOFI == 1 & split_taxa$Tara == 1 & split_taxa$Polar == 0)] <- "NCOG + TARA"
  split_taxa$endemic[which(split_taxa$CalCOFI == 1 & split_taxa$Tara == 0 & split_taxa$Polar == 1)] <- "NCOG + TARA Polar"
  split_taxa$endemic[which(split_taxa$CalCOFI == 0 & split_taxa$Tara == 1 & split_taxa$Polar == 1)] <- "TARA + TARA Polar"
  split_taxa$endemic[which(split_taxa$CalCOFI == 1 & split_taxa$Tara == 1 & split_taxa$Polar == 1)] <- "All Datasets"
  
  mean_groups <- split_taxa %>%
    filter(!is.na(endemic), !is.na(C)) %>%
    group_by(endemic,C) %>% summarise(count = n()) %>% mutate(prop = count/sum(count))
  
 
  
  table <- pivot_wider(mean_groups[,1:3], names_from = "endemic", values_from = "count")  
  
  colnames(table)[1] <- "Division"
  
  table <- table[c(1:23,25:36,24),c(1,3,6,8,4,5,7,2)]
  
  table[is.na(table)] <- 0
  
  tot <- data.frame(matrix(c("Total",colSums(table[,2:8], na.rm = TRUE)),1,8))
  colnames(tot) <- colnames(table)
  
  table <- rbind(table, tot)
  
  pdf(output_table, width = 11.5, height = 11.5)
  plot.new()
  print(grid.table(table, rows = NULL))
  dev.off()
  
}


######## Suppl Figure 3: Relative Mean Abundance Nearshore-Offshore  #######

source("analysis/suppl_fig_3_S.R")

######## Suppl Figure 4: SOMs vs ef-ratio  #######

source("analysis/suppl_ef_ratio_figure.R")

##### Suppl Figure 5: SOMs small groups #####

suppl_fig_5_func <- function(in_list = plot_list, file_name = "figures_S/supp_fig_5_S.pdf", tsize  = 17){
  
  in_list[[1]]$bp <- in_list[[1]]$bp + theme(axis.title.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none",
                                             plot.margin = ggplot2::margin(5,5,0,5)) + ggtitle("Prochlorococcus") +
  annotate(x = -126.95, y = 36.85, label = "a", size = 7, geom = "text", fontface = "bold") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  in_list[[2]]$bp <- in_list[[2]]$bp + theme(axis.title.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none",
                                             plot.margin = ggplot2::margin(5,5,0,5)) + ggtitle("Synechococcus") +
    annotate(x = -126.95, y = 36.85, label = "b", size = 7, geom = "text", fontface = "bold") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  in_list[[3]]$bp <- in_list[[3]]$bp + theme(axis.title.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none",
                                             plot.margin = ggplot2::margin(5,5,0,5)) + ggtitle("Flavobacteriales") +
    annotate(x = -126.95, y = 36.85, label = "c", size = 7, geom = "text", fontface = "bold") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  in_list[[4]]$bp <- in_list[[4]]$bp + theme(axis.title.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none",
                                             plot.margin = ggplot2::margin(5,5,0,5)) + ggtitle("Rhodobacterales") +
    annotate(x = -126.95, y = 36.85, label = "d", size = 7, geom = "text", fontface = "bold") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  in_list[[5]]$bp <- in_list[[5]]$bp + theme(axis.title.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             plot.margin = ggplot2::margin(5,5,0,5),
                                             legend.title = element_text(size = tsize, vjust = 0.5)) + ggtitle("SAR 11 Clade") +
    annotate(x = -126.95, y = 36.85, label = "e", size = 7, geom = "text", fontface = "bold") +
    scale_x_continuous(breaks = c(-125,-122,-119)) + labs(fill = "Frequency")
  
  in_list[[7]]$bp <- in_list[[7]]$bp + theme(axis.title.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none",
                                             plot.margin = ggplot2::margin(5,5,0,5)) + ggtitle("Diatoms") +
    annotate(x = -126.95, y = 36.85, label = "f", size = 7, geom = "text", fontface = "bold") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  in_list[[8]]$bp <- in_list[[8]]$bp + theme(axis.title.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none",
                                             plot.margin = ggplot2::margin(5,5,0,5)) + ggtitle("Dinoflagellates") +
    annotate(x = -126.95, y = 36.85, label = "g", size = 7, geom = "text", fontface = "bold") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  in_list[[9]]$bp <- in_list[[9]]$bp + theme(axis.title.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none",
                                             plot.margin = ggplot2::margin(5,5,0,5)) + ggtitle("Syndiniales") +
    annotate(x = -126.95, y = 36.85, label = "h", size = 7, geom = "text", fontface = "bold") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  in_list[[10]]$bp <- in_list[[10]]$bp + theme(axis.title.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize),
                                               legend.position = "none",
                                               plot.margin = ggplot2::margin(5,5,0,5)) + ggtitle("Haptophytes") +
    annotate(x = -126.95, y = 36.85, label = "i", size = 7, geom = "text", fontface = "bold") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  in_list[[11]]$bp <- in_list[[11]]$bp + theme(axis.title.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize),
                                               legend.position = "none",
                                               plot.margin = ggplot2::margin(5,5,0,5)) + ggtitle("Chlorophytes") +
    annotate(x = -126.95, y = 36.85, label = "j", size = 7, geom = "text", fontface = "bold") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  in_list[[12]]$bp <- in_list[[12]]$bp + theme(axis.title.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize),
                                               legend.position = "none",
                                               plot.margin = ggplot2::margin(5,5,0,5)) + ggtitle("Metazoans") +
    annotate(x = -126.95, y = 36.85, label = "k", size = 7, geom = "text", fontface = "bold") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  in_list[[1]]$rp <- in_list[[1]]$rp + theme(legend.position = "none",
                                             plot.margin = ggplot2::margin(0,5,5,5)) + 
    scale_x_continuous(breaks = c(-125,-122,-119))
  in_list[[2]]$rp <- in_list[[2]]$rp + theme(axis.title.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             legend.position = "none",
                                             plot.margin = ggplot2::margin(0,5,5,5)) + 
    scale_x_continuous(breaks = c(-125,-122,-119))
  in_list[[3]]$rp <- in_list[[3]]$rp + theme(axis.title.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             legend.position = "none",
                                             plot.margin = ggplot2::margin(0,5,5,5)) + 
    scale_x_continuous(breaks = c(-125,-122,-119))
  in_list[[4]]$rp <- in_list[[4]]$rp + theme(axis.title.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             legend.position = "none",
                                             plot.margin = ggplot2::margin(0,5,5,5)) + 
    scale_x_continuous(breaks = c(-125,-122,-119))
  in_list[[5]]$rp <- in_list[[5]]$rp + theme(axis.title.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             plot.margin = ggplot2::margin(0,5,5,5),
                                             legend.title = element_text(size = tsize, vjust = 0.5)) + 
    scale_x_continuous(breaks = c(-125,-122,-119)) + labs(fill = "Frequency")
  
  in_list[[7]]$rp <- in_list[[7]]$rp + theme(legend.position = "none") + 
    scale_x_continuous(breaks = c(-125,-122,-119))
  in_list[[8]]$rp <- in_list[[8]]$rp + theme(axis.title.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             legend.position = "none",
                                             plot.margin = ggplot2::margin(0,5,5,5)) + 
    scale_x_continuous(breaks = c(-125,-122,-119))
  in_list[[9]]$rp <- in_list[[9]]$rp + theme(axis.title.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             legend.position = "none",
                                             plot.margin = ggplot2::margin(0,5,5,5)) + 
    scale_x_continuous(breaks = c(-125,-122,-119))
  in_list[[10]]$rp <- in_list[[10]]$rp + theme(axis.title.y = element_blank(),
                                              axis.text.y = element_blank(),
                                              plot.title = element_text(hjust = 0),
                                              legend.position = "none",
                                              plot.margin = ggplot2::margin(0,5,5,5)) + 
    scale_x_continuous(breaks = c(-125,-122,-119))
  in_list[[11]]$rp <- in_list[[11]]$rp + theme(axis.title.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none",
                                               plot.margin = ggplot2::margin(0,5,5,5)) + 
    scale_x_continuous(breaks = c(-125,-122,-119))
  in_list[[12]]$rp <- in_list[[12]]$rp + theme(axis.title.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none",
                                               plot.margin = ggplot2::margin(0,5,5,5)) + 
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  layout <- "ABCDEK
             FGHIJK
             LMNOPQ
             RSTUVW"
  
  
  pdf(file = file_name, width = 16, height = 12)
  print(in_list[[1]]$bp + in_list[[2]]$bp + in_list[[3]]$bp + in_list[[4]]$bp + in_list[[5]]$bp + 
          in_list[[1]]$rp + in_list[[2]]$rp + in_list[[3]]$rp + in_list[[4]]$rp + in_list[[5]]$rp + guide_area() +
          in_list[[7]]$bp + in_list[[8]]$bp + in_list[[9]]$bp + in_list[[10]]$bp + in_list[[11]]$bp + in_list[[12]]$bp + 
          in_list[[7]]$rp + in_list[[8]]$rp + in_list[[9]]$rp + in_list[[10]]$rp + in_list[[11]]$rp + in_list[[12]]$rp +
          plot_layout(guides = "collect", design = layout))
  dev.off()
}

##### Suppl Figure 6: Variable AIC small groups #####

full_aic_table_figure(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                        "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                        "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                      in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "SAR 11 Clade", 
                                         "Diatoms",
                                         "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"),
                      minimum_tp = 4, width_plot = 17.5,
                      figure_name_2 = paste0("figures_S/supp_fig_6_S",".pdf"),
                      title_name = "", tsize = 15)

##### Suppl Figure 7: Surface and DCM Diversity Maps #####

source("analysis/response_figure.R")

#### Suppl Figure 8: Diversity Maps #####

suppl_fig_8_func <- function(file_name = "figures_S/supp_fig_8_S.pdf",
                               in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                                 "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                 "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                               in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar Clade", 
                                                  "Diatoms",
                                                  "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"), tsize = 20){
  
  div_map <- list()
  
  for (i in 1:length(in_group_list)) {
    
    div_map[[i]] <- diveristy_figure(map_file = paste0("output/", in_group_list[i], "_map_S.Rdata"),
                                     full_dat = paste0("output/", in_group_list[i], "_full_data_S.Rdata"),
                                     figure_start = paste0("figures/diversity/", in_group_list[i], "_"),
                                     main = in_group_names[i])
    
  }
  
  pro <- div_map[[1]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+2),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize-7),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-7),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    ggtitle("Prochlorococcus") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.9, label = "a", size = 7, geom = "text", fontface = "bold") 
  
  syn <- div_map[[2]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+2),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize-7),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-7),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Synecococcus") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.9, label = "b", size = 7, geom = "text", fontface = "bold") 
  
  flavo <- div_map[[3]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+2),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize-7),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-7),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Flavobacteriales") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.9, label = "c", size = 7, geom = "text", fontface = "bold")  
  
  rho <- div_map[[4]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+2),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize-7),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-7),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Rhodobacterales") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.9, label = "d", size = 7, geom = "text", fontface = "bold") 
  
  sar <- div_map[[5]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+2),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize-7),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-7),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    ggtitle("SAR 11 Clade") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.9, label = "e", size = 7, geom = "text", fontface = "bold") 
  
  diatom <- div_map[[6]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+2),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize-7),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-7),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Diatoms") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.9, label = "f", size = 7, geom = "text", fontface = "bold") 
  
  dino <- div_map[[7]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+2),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize-7),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-7),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Dinoflagellates") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.9, label = "g", size = 7, geom = "text", fontface = "bold") 
  
  sindin <- div_map[[8]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+2),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize-7),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-7),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = tsize),
          axis.title.x = element_text(size = tsize)) +
    ggtitle("Syndiniales") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.9, label = "h", size = 7, geom = "text", fontface = "bold")  
  
  hapto <- div_map[[9]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+2),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize-7),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-7)) +
    ggtitle("Haptophytes") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.9, label = "i", size = 7, geom = "text", fontface = "bold")  
  
  chloro <- div_map[[10]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+2),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize-7),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-7),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Chlorophytes") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.9, label = "j", size = 7, geom = "text", fontface = "bold")  
  
  meta <- div_map[[11]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+2),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize-7),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize-7),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Metazoans") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks = c(-125,-122,-119)) +
    annotate(x = -126.95, y = 36.9, label = "k", size = 7, geom = "text", fontface = "bold")  
  
  plot_div <-  pro + syn + flavo + rho +
    sar + diatom + dino + sindin +
    hapto + chloro + meta + plot_spacer() + plot_layout(ncol = 4)
  
  
  pdf(file = file_name, width = 19, height = 17)
  print(plot_div)
  dev.off()
  
}


##### Suppl Figure 9: Diversity Importance Small Groups #####

full_aic_table_figure_diversity_sign(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                                       "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                       "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                                     in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar 11 Clade", 
                                                        "Diatoms",
                                                        "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"),
                                     figure_name_2 = "figures_S/supp_fig_9_S.pdf",
                                     col = 27, width_plot = 17.5, tsize = 14)

####### Suppl Figure 10: Diversity - Productivity ########


suppl_fig_10_func <- function(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                               "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                               "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                             in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales",
                                                "Rhodobacterales", "SAR 11 Clade", "Diatoms",
                                                "Dinoflagellates", "Syndiniales", "Haptophytes",
                                                "Chlorophytes","Metazoans"), tsize = 15,
                             out_fig_name = "figures_S/supp_fig_10_S.pdf"){
  
  
  out_plot_list <- list()
  
  #pull IntC14 Data
  metadata <- read.csv("data/NCOG_sample_log_DNA_meta_2014-2020.csv")
  lab_vect <- c("a","b","c","d","e",
                "f","g","h","i","j", "k")
  
  for (i in 1:length(in_group_list)) {
    
    load(paste0("output/", in_group_list[i],"_full_data_S.RData"))
    
    full_dat$samp <- gsub("_S", "", full_dat$Sample.Name)
    
    full_dat$IntC14 <- metadata$IntC14[match(full_dat$samp, metadata$Sample.Name)]
    
    prodo_dat <- full_dat %>%
      filter(!is.na(IntC14))
    
    prodo_dat$phase <- substr(prodo_dat$Cruise,1,4)
    
    prodo_dat$phase[which(prodo_dat$phase == "2014" | prodo_dat$phase == "2015" | prodo_dat$phase == "2016")] <- "2014-2016"
    prodo_dat$phase[which(prodo_dat$phase == "2017" | prodo_dat$phase == "2018")] <- "2017-2018"
    prodo_dat$phase[which(prodo_dat$phase == "2019" | prodo_dat$phase == "2020")] <- "2019-2020"
    
    prodo_dat$phase <- as.factor(prodo_dat$phase)
    prodo_dat$phase <- factor(prodo_dat$phase, levels = c("2014-2016", "2017-2018", "2019-2020"))
    
    gam_dat <- prodo_dat 
    
    gam_dat$IntC14_day <- gam_dat$IntC14*2
    
    warm_lm <- summary(gam(richness ~ IntC14_day,data = gam_dat %>% filter(phase == "2014-2016")))
    
    warm_p <- warm_lm$p.pv[2]
    
    cool_lm <- summary(gam(richness ~ IntC14_day,data = gam_dat %>% filter(phase == "2017-2018")))
    
    cool_p <- cool_lm$p.pv[2]
    
    new_lm <- summary(gam(richness ~ IntC14_day,data = gam_dat %>% filter(phase == "2019-2020")))
    
    new_p <- new_lm$p.pv[2]
    
    fit_mat <- matrix(c("2014-2016", "2017-2018", "2019-2020",
                        warm_p, cool_p, new_p),3,2, byrow = FALSE) %>% as.data.frame()
    
    fit_mat$V2 <- as.numeric(fit_mat$V2)
    
    gam_dat$sig <- "Y"
    gam_dat$sig[which(!is.na(match(gam_dat$phase, fit_mat$V1[which(fit_mat$V2 > 0.05)])))] <- NA
    
    
    round_list <- c(round(warm_p, 3),round(cool_p, 3),round(new_p, 3))
    round_list[which(round_list == 0)] <- "< 0.01"
    round_list[which(round_list > 0)] <- paste("= ", round_list[which(round_list  > 0)])
    
    p_window <- dist(range(gam_dat$richness)) * 0.1
    
    out_plot_list[[i]] <- ggplot() +
      geom_point(data = gam_dat, aes(x = IntC14*2, y = richness, color = phase)) + scale_x_log10() + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
      stat_smooth(data = gam_dat %>% filter(sig == "Y"), aes(x = IntC14*2, y = richness, color = phase, fill = phase),
                  method = "gam", se = 0.99) + labs(color = "Phase") +
      scale_color_manual(values = c("red", "blue","gold3")) +
      scale_fill_manual(values = c("red", "blue","gold3")) + ylab("Richness") +
      xlab("mgC/m3 per one light day") +
      theme(panel.background = element_blank(),
            panel.border = element_rect(fill = NA, color = "black"),
            legend.text = element_text(size = tsize),
            legend.title = element_text(size = tsize),
            axis.text = element_text(size = tsize),
            axis.title = element_text(size = tsize),
            plot.title = element_text(size = tsize),
            legend.key = element_blank()) +
      ggtitle(in_group_names[i]) + guides(fill = FALSE) +
      annotate(geom = "text", x = 100, y = max(prodo_dat$richness, na.rm = TRUE)*0.95, label = lab_vect[i],
               fontface = "bold", size = 7) +
      annotate(geom = "rect", ymax = -2,
               ymin = -Inf, xmin = -Inf, xmax = Inf,
               fill = "white", color = "white") +
      annotate(geom = "text", y = -p_window*0.75, x = 260, 
               label = paste0("p ",
                              round_list[1]),
               color = "red", size = 2.5, fontface = "bold") +
      annotate(geom = "text", y = -p_window*0.75, x = 650, 
               label = paste0("p ",
                              round_list[2]),
               color = "blue", size = 2.5, fontface = "bold") +
      annotate(geom = "text", y = -p_window*0.75, x = 1600, 
               label = paste0("p ",
                              round_list[3]),
               color = "gold3", size = 2.5, fontface = "bold") +
      coord_cartesian(ylim = c(-p_window,max(gam_dat$richness)))
    
    
  }
  
  
  out_plot <- (out_plot_list[[1]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank())) +
    (out_plot_list[[2]] + theme(axis.text.x = element_blank(), axis.title = element_blank())) + 
    (out_plot_list[[3]] + theme(axis.text.x = element_blank(), axis.title = element_blank())) +
    (out_plot_list[[4]] + theme(axis.text.x = element_blank(), axis.title = element_blank())) +
    (out_plot_list[[5]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank())) +
    (out_plot_list[[6]] + theme(axis.text.x = element_blank(), axis.title = element_blank())) +
    (out_plot_list[[7]] + theme(axis.text.x = element_blank(), axis.title = element_blank())) +
    (out_plot_list[[8]] + theme(axis.title.y = element_blank())) +
    (out_plot_list[[9]] + theme()) +
    (out_plot_list[[10]] + theme(axis.title.y = element_blank())) +
    (out_plot_list[[11]] + theme(axis.title.y = element_blank())) +
    guide_area()  + plot_layout(ncol = 4, guides = "collect") 
  
  pdf(out_fig_name, width = 12, height = 9)
  print(out_plot)
  dev.off()
  
}







##### Suppl Figure 11: Environmental Differences #####

suppl_fig_11_func <- function(nutrient_dat = "output/mld_mean_profiles_S.Rdata",
                             suppl_plot = "figures_S/supp_fig_11_S.pdf",
                             tsize = 12, psize = 4){
  
  map <- map_data("world")
  
  load(nutrient_dat)
  
  early_phase <- early_phase %>% filter(Sta_ID %in% late_phase$Sta_ID)
  
  early_phase$n_diff <- early_phase$mean_no3 - late_phase$mean_no3
  early_phase$p_diff <- early_phase$mean_po4 - late_phase$mean_po4
  early_phase$si_diff <- early_phase$mean_sio4 - late_phase$mean_sio4
  early_phase$nc_diff <- early_phase$mean_ncd - late_phase$mean_ncd
  early_phase$ml_diff <- early_phase$mean_mld - late_phase$mean_mld
  early_phase$t_diff <- early_phase$mean_temp - late_phase$mean_temp
  early_phase$sa_diff <- early_phase$mean_sal - late_phase$mean_sal
  early_phase$chl_diff <- early_phase$mean_chl - late_phase$mean_chl
  
  colors <- c("black", "black", "black", "black",
              "black", "black", "black", "black")
  
  early_t <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_temp),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "°C", low = "white", high = colors[1],
                        limits = c(13,20), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(1.1,1.08),
          plot.tag = element_text(size = tsize + 5),
          plot.margin = ggplot2::margin(40,5,5,5)) + 
    ggtitle("2014-2016") +
    annotate(x = -126.95, y = 36.85, label = "a", size = 7, geom = "text", fontface = "bold") +
    labs(tag = "Mean Temperature (°C)") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  late_t <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_temp),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "°C", low = "white", high = colors[1],
                        limits = c(13,20), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.background = element_blank()) +
    ggtitle("2017-2018") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  diff_t <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = t_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(Delta~" °C"), low = "blue", high = "red",
                         mid = "white", midpoint = 0,
                         limits = c(-0.5,2), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.tag.position = c(0,1.08),
          plot.tag = element_text(size = tsize + 5, hjust = 0),
          plot.margin = ggplot2::margin(40,5,5,5)) + 
    ggtitle("(2014-2016) - (2017-2018)") +
    labs(tag = expression(Delta~Temperature)) +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  early_sa <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_sal),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "PSU", low = "white", high = colors[2],
                        limits = c(33,34), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize + 3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(1.1,1.08),
          plot.tag = element_text(size = tsize + 5),
          plot.margin = ggplot2::margin(40,5,5,5)) + 
    ggtitle("2014-2016") +
    annotate(x = -126.95, y = 36.85, label = "c", size = 7, geom = "text", fontface = "bold") +
    labs(tag = "Mean Salinty (PSU)") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  late_sa <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_sal),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "PSU", low = "white", high = colors[2],
                        limits = c(33,34), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.background = element_blank()) +
    ggtitle("2017-2018") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  diff_sa <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = sa_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(Delta~" PSU"), low = "blue", high = "red",
                         mid = "white", midpoint = 0,
                         limits = c(-0.5,0), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.tag.position = c(0,1.08),
          plot.tag = element_text(size = tsize + 5, hjust = 0),
          plot.margin = ggplot2::margin(40,5,5,5)) + 
    ggtitle("(2014-2016) - (2017-2018)") +
    labs(tag = expression(Delta~Salinity)) +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  early_ml <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_mld),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Depth (m)", low = "white", high = colors[3],
                        limits = c(10,50), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(1.1,1.08),
          plot.tag = element_text(size = tsize + 5),
          plot.margin = ggplot2::margin(40,5,5,5)) + 
    ggtitle("2014-2016") +
    annotate(x = -126.95, y = 36.85, label = "e", size = 7, geom = "text", fontface = "bold") +
    labs(tag = "Mean Mixed Layer Depth (m)") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  late_ml <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_mld),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Depth (m)", low = "white", high = colors[3],
                        limits = c(10,50), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.background = element_blank()) +
    ggtitle("2017-2018")  +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  diff_ml <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = ml_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(Delta~"m"), low = "blue", high = "red",
                         mid = "white", midpoint = 0,
                         limits = c(-20,10), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.tag.position = c(0,1.08),
          plot.tag = element_text(size = tsize + 5, hjust = 0),
          plot.margin = ggplot2::margin(40,5,5,5)) + 
    ggtitle("(2014-2016) - (2017-2018)") +
    labs(tag = expression(Delta~`Mixed Layer Depth`)) +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  early_nc <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_ncd),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Depth (m)", low = "white", high = colors[4],
                        limits = c(0,120), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(1.1,1.08),
          plot.tag = element_text(size = tsize + 5),
          plot.margin = ggplot2::margin(40,5,5,5)) + 
    ggtitle("2014-2016") +
    annotate(x = -126.95, y = 36.85, label = "g", size = 7, geom = "text", fontface = "bold") +
    labs(tag = "Mean Nitracline Depth (m)") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  late_nc <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_ncd),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Depth (m)", low = "white", high = colors[4],
                        limits = c(0,120), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.background = element_blank()) +
    ggtitle("2017-2018") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  diff_nc <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = nc_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(Delta~"m"), low = "blue", high = "red",
                         mid = "white", midpoint = 0,
                         limits = c(-20,20), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.tag.position = c(0,1.08),
          plot.tag = element_text(size = tsize + 5, hjust = 0),
          plot.margin = ggplot2::margin(40,5,5,5)) + 
    ggtitle("(2014-2016) - (2017-2018)") +
    labs(tag = expression(Delta~`Nitracline Depth`)) +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  
  early_chl <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_chl),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "μg/L", low = "white", high = colors[5],
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(1.1,1.08),
          plot.tag = element_text(size = tsize + 5),
          plot.margin = ggplot2::margin(40,5,5,5)) + 
    ggtitle("2014-2016") +
    annotate(x = -126.95, y = 36.85, label = "b", size = 7, geom = "text", fontface = "bold") +
    labs(tag = "Mean Chlorophyll-a (μg/L)") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  late_chl <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_chl),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "μg/L", low = "white", high = colors[5],
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.background = element_blank()) +
    ggtitle("2017-2018") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  diff_chl <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = chl_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(Delta~"μg/L"), low = "blue", high = "red",
                         mid = "white", midpoint = 0,
                         limits = c(-2,1), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.tag.position = c(0,1.08),
          plot.tag = element_text(size = tsize + 5, hjust = 0),
          plot.margin = ggplot2::margin(40,5,5,5)) + 
    ggtitle("(2014-2016) - (2017-2018)") +
    labs(tag = expression(Delta~`Chlorophyll-a`)) +
    scale_x_continuous(breaks = c(-125,-122,-119))

  early_n <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_no3),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(paste(NO[3]," ", mu,M)), low = "white", high = colors[6],
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(1.1,1.08),
          plot.tag = element_text(size = tsize + 5),
          plot.margin = ggplot2::margin(40,5,5,5)) + 
    ggtitle("2014-2016") +
    annotate(x = -126.95, y = 36.85, label = "d", size = 7, geom = "text", fontface = "bold") +
    labs(tag = "Mean Nitrate (μM)") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  late_n <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_no3),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(paste(NO[3]," ", mu,M)), low = "white", high = colors[6],
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.background = element_blank()) +
    ggtitle("2017-2018") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  diff_n <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = n_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(Delta~"NO3 μM"),
                        low = "blue", high = "white",
                        limits = c(-6,0), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.tag.position = c(0,1.08),
          plot.tag = element_text(size = tsize + 5, hjust = 0),
          plot.margin = ggplot2::margin(40,5,5,5)) + 
    ggtitle("(2014-2016) - (2017-2018)")  +
    labs(tag = expression(Delta~Nitrate)) +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  early_p <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_po4),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(paste(PO[4]," ", mu,M)), low = "white", high = colors[7],
                        limits = c(0,1), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0),
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(1.1,1.08),
          plot.tag = element_text(size = tsize + 5),
          plot.margin = ggplot2::margin(40,5,5,5)) +
    ggtitle("2014-2016") +
    annotate(x = -126.95, y = 36.85, label = "f", size = 7, geom = "text", fontface = "bold") +
    labs(tag = "Mean Phosphate (μM)") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  late_p <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_po4),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(paste(PO[4]," ", mu,M)), low = "white", high = colors[7],
                        limits = c(0,1), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0),
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.background = element_blank()) +
    ggtitle("2017-2018") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  diff_p <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = p_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(Delta~"PO4 μM"), low = "blue", high = "red",
                         mid = "white", midpoint = 0,
                        limits = c(-0.25,0.25), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0), axis.line = element_blank(),
          legend.justification=c(1,0),
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.tag.position = c(0,1.08),
          plot.tag = element_text(size = tsize + 5, hjust = 0),
          plot.margin = ggplot2::margin(40,5,5,5)) + ggtitle("(2014-2016) - (2017-2018)")  +
    labs(tag = expression(Delta~Phosphate)) +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  early_si <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_sio4),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(paste(SiO[4]," ", mu,M)), low = "white", high = colors[8],
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0),
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(1.1,1.08),
          plot.tag = element_text(size = tsize + 5),
          plot.margin = ggplot2::margin(40,5,5,5)) +
    ggtitle("2014-2016") +
    annotate(x = -126.95, y = 36.85, label = "h", size = 7, geom = "text", fontface = "bold") +
    labs(tag = "Mean Silicate (μM)") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  late_si <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_sio4),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(paste(SiO[4]," ", mu,M)), low = "white", high = colors[8],
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0.5), axis.line = element_blank(),
          legend.justification=c(1,0),
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.background = element_blank()) +
    ggtitle("2017-2018") +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
  diff_si <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = si_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(Delta~"SiO4 μM"), low = "blue", high = "red",
                         mid = "white", midpoint = 0,
                         limits = c(-1.5,1.5), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize+3, hjust = 0), axis.line = element_blank(),
          legend.justification=c(1,0),
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          plot.tag.position = c(0,1.08),
          plot.tag = element_text(size = tsize + 5, hjust = 0),
          plot.margin = ggplot2::margin(40,5,5,5)) + ggtitle("(2014-2016) - (2017-2018)")  +
    labs(tag = expression(Delta~Silicate)) +
    scale_x_continuous(breaks = c(-125,-122,-119))
  
 plot_comb <- (early_t + late_t + diff_t + plot_layout(ncol = 3, guides = "collect"))/
   (early_sa + late_sa + diff_sa + plot_layout(ncol = 3, guides = "collect"))/
   (early_ml + late_ml + diff_ml + plot_layout(ncol = 3, guides = "collect"))/
   (early_nc + late_nc + diff_nc + plot_layout(ncol = 3, guides = "collect"))|
   (early_chl + late_chl + diff_chl + plot_layout(ncol = 3, guides = "collect"))/
   (early_n + late_n + diff_n + plot_layout(ncol = 3, guides = "collect"))/
   (early_p + late_p + diff_p + plot_layout(ncol = 3, guides = "collect"))/
   (early_si + late_si + diff_si + plot_layout(ncol = 3, guides = "collect")) 
 
 pdf(file = suppl_plot, width = 22, height = 18)
 print(plot_comb)
 dev.off()
 
  }

##### Suppl Figure 12: Community Differences ######

suppl_fig_12_func <- function(in_list = fig_list, file_name = "figures_S/supp_fig_12_S.pdf", tsize  = 12){
  
    low <- 0.125
    high <- 0.6
    
    all_plots <- (in_list[[6]]$deep + ggtitle("Archaea") +
                    scale_fill_gradient(limits = c(low,high),
                                        low = "white", high = "red", oob = scales::squish) +
                    labs(fill = "Bray-Curtis\nSimilarity") +
                    theme(legend.position = "none",
                          plot.background = element_blank()) +
                    scale_x_continuous(breaks = c(-125,-122,-119)) +
                    annotate(x = -126.95, y = 36.85, label = "a", size = 7, geom = "text", fontface = "bold")) +
      (in_list[[13]]$deep + ggtitle("Bacteria") +
         scale_fill_gradient(limits = c(low,high),
                             low = "white", high = "red", oob = scales::squish) +
         labs(fill = "Bray-Curtis\nSimilarity") +
         theme(axis.title.y = element_blank(),
               axis.text.y = element_blank())+
         theme(legend.position = "none",
               plot.background = element_blank()) +
         scale_x_continuous(breaks = c(-125,-122,-119)) +
         annotate(x = -126.95, y = 36.85, label = "b", size = 7, geom = "text", fontface = "bold")) +
      (in_list[[14]]$deep + ggtitle("Cyanobacteria") +
         scale_fill_gradient(limits = c(low,high),
                             low = "white", high = "red", oob = scales::squish) +
         labs(fill = "Bray-Curtis\nSimilarity") +
         theme(axis.title.y = element_blank(),
               axis.text.y = element_blank())+
         theme(legend.position = "none",
               plot.background = element_blank()) +
         scale_x_continuous(breaks = c(-125,-122,-119)) +
         annotate(x = -126.95, y = 36.85, label = "c", size = 7, geom = "text", fontface = "bold")) +
      (in_list[[15]]$deep + ggtitle("Photosynthetic eukaryotic\nprotists") +
         scale_fill_gradient(limits = c(low,high),
                             low = "white", high = "red", oob = scales::squish) +
         labs(fill = "Bray-Curtis\nSimilarity") +
         theme(axis.title.y = element_blank(),
               axis.text.y = element_blank())+
         theme(legend.position = "none",
               plot.background = element_blank()) +
         scale_x_continuous(breaks = c(-125,-122,-119)) +
         annotate(x = -126.95, y = 36.85, label = "d", size = 7, geom = "text", fontface = "bold")) +
      (in_list[[16]]$deep + ggtitle("Heterotrophic eukaryotic\nprotists") +
         scale_fill_gradient(limits = c(low,high),
                             low = "white", high = "red", oob = scales::squish) +
         labs(fill = "Bray-Curtis\nSimilarity") +
         theme(axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               plot.background = element_blank()) +
         scale_x_continuous(breaks = c(-125,-122,-119)) +
         annotate(x = -126.95, y = 36.85, label = "e", size = 7, geom = "text", fontface = "bold")) +
      plot_layout(ncol = 5, nrow = 1, byrow = TRUE, guides = "collect") +
      plot_annotation(title = "DCM (2014-2016) vs (2017-2018)",
                      theme = theme(plot.title = element_text(size = 18, hjust = 0.5))) 
    
    pdf(file = file_name, width = 15, height = 4)
    print(all_plots)
    dev.off()
  
}

##### Suppl Figure 13: Difference Plots ######

suppl_fig_13_func <- function(in_list = fig_list, file_name = "figures_S/supp_fig_13_S.pdf", tsize = 12,
                              in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales",
                                                 "Rhodobacterales", "SAR 11 Clade", "Diatoms",
                                                 "Dinoflagellates", "Syndiniales", "Haptophytes",
                                                 "Chlorophytes","Metazoans")){
  
  groups <- c(1:5,7:12)
  
  low = 0.1
  high = 0.8
  
  lab_list <- c("a","b","c","d","e","f","g","h", "i", "j", "k")
  
  plot_list <- list()
  for (i in 1:length(groups)) {
    
    plot_list[[i]] <- in_list[[groups[i]]]$surf + 
       scale_fill_gradient(limits = c(low,high),low = "white", high = "red", oob = scales::squish,
                           breaks = scales::pretty_breaks(n = 3)) +
       labs(fill = "Bray-Curtis\nSimilarity", subtitle = NULL, title = in_group_names[i])+
       theme(plot.margin = ggplot2::margin(5,5,5,5),
             legend.direction = "vertical",
             legend.position = c(0.025,0.025), 
             legend.justification = c(0,0),
             legend.margin = ggplot2::margin(4,15,4,4)) + 
      annotate(x = -126.95, y = 36.99, label = lab_list[i], size = 7, geom = "text", fontface = "bold") +
      scale_x_continuous(breaks = c(-125,-122,-119))
    
  }

  fig_plot <- (plot_list[[1]] + theme(legend.position = "none",
                          axis.text.x = element_blank(), axis.title.x = element_blank())) +
    (plot_list[[2]] + theme(legend.position = "none",
                            axis.text.x = element_blank(), axis.title.x = element_blank(),
                            axis.text.y = element_blank(), axis.title.y = element_blank())) +
    (plot_list[[3]] + theme(legend.position = "none",
                            axis.text.x = element_blank(), axis.title.x = element_blank(),
                            axis.text.y = element_blank(), axis.title.y = element_blank())) +
    (plot_list[[4]] + theme(legend.position = "none",
                            axis.text.x = element_blank(), axis.title.x = element_blank(),
                            axis.text.y = element_blank(), axis.title.y = element_blank())) +
    (plot_list[[5]] + theme(legend.position = "none",
                            axis.text.x = element_blank(), axis.title.x = element_blank())) +
    (plot_list[[6]] + theme(legend.position = "none",
                            axis.text.x = element_blank(), axis.title.x = element_blank(),
                            axis.text.y = element_blank(), axis.title.y = element_blank())) +
    (plot_list[[7]] + theme(legend.position = "none",
                            axis.text.x = element_blank(), axis.title.x = element_blank(),
                            axis.text.y = element_blank(), axis.title.y = element_blank())) +
    (plot_list[[8]] + theme(legend.position = "none",
                            axis.text.y = element_blank(), axis.title.y = element_blank())) +
    (plot_list[[9]] + theme(legend.position = "none")) +
    (plot_list[[10]] + theme(legend.position = "none",
                            axis.text.y = element_blank(), axis.title.y = element_blank())) +
    (plot_list[[11]] + theme(axis.text.y = element_blank(), axis.title.y = element_blank())) + guide_area() +
    plot_layout(guides = "collect") +
    plot_annotation(title = "Surface\n(2014-2016) vs (2017-2018)",
                    theme = theme(plot.title = element_text(size = 14, hjust = 0.5)))
  
  
  
  pdf(file = file_name, width = 13, height = 10)
  print(fig_plot)
  dev.off()
  
}

###### Suppl Figure 14: Nitracline slope #####

source("analysis/CUTI_BEUTI_Plots.R")

###### Suppl Figure 15: Prop Nearshore over Time #####

suppl_fig_15_func <- function(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                                 "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                 "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                               in_group_names = c("a Prochlorococcus", "b Synecococcus", "c Flavobacteriales",
                                                  "d Rhodobacterales", "e SAR 11 Clade", "f Diatoms",
                                                  "g Dinoflagellates", "h Syndiniales", "i Haptophytes",
                                                  "j Chlorophytes","k Metazoans"),
                               gradient_plot_file = "figures_S/supp_fig_15_S.pdf",
                               tsize = 12){
  
  
  plot_list <- list()
  
  for (i in  1:length(in_group_list)) {
    
    load(paste0("output/",in_group_list[i],"_diffs_S.Rdata"))
    
    plot_list[[i]] <- ts_plot
    
  }
  
  
  pro <- plot_list[[1]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    ggtitle("Prochlorococcus") + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    annotate(geom = "text", x = dmy("1-1-2014"), y = 0.95, 
             label = "a", fontface = "bold", size = 6) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
  
  syn <- plot_list[[2]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Synecococcus")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    annotate(geom = "text", x = dmy("1-1-2014"), y = 0.95, 
             label = "b", fontface = "bold", size = 6) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
  
  flavo <- plot_list[[3]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Flavobacteriales")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    annotate(geom = "text", x = dmy("1-1-2014"), y = 0.95, 
             label = "c", fontface = "bold", size = 6) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
  
  rho <- plot_list[[4]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Rhodobacterales")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    annotate(geom = "text", x = dmy("1-1-2014"), y = 0.95, 
             label = "d", fontface = "bold", size = 6) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
  
  sar <- plot_list[[5]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    ggtitle("SAR 11 Clade") + 
    ylab("Proportion of Samples\nIdentified as Nearshore") + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    annotate(geom = "text", x = dmy("1-1-2014"), y = 0.95, 
             label = "e", fontface = "bold", size = 6) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
  
  diatom <- plot_list[[6]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Diatoms")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    annotate(geom = "text", x = dmy("1-1-2014"), y = 0.95, 
             label = "f", fontface = "bold", size = 6) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
  
  dino <- plot_list[[7]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Dinoflagellates")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    annotate(geom = "text", x = dmy("1-1-2014"), y = 0.95, 
             label = "g", fontface = "bold", size = 6) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
  
  sindin <- plot_list[[8]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Syndiniales") + 
    xlab("Date")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    annotate(geom = "text", x = dmy("1-1-2014"), y = 0.95, 
             label = "h", fontface = "bold", size = 6) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
  
  hapto <- plot_list[[9]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Haptophytes") + 
    xlab("Date") + 
    ylab("Proportion of Samples\nIdentified as Nearshore") + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    annotate(geom = "text", x = dmy("1-1-2014"), y = 0.95, 
             label = "i", fontface = "bold", size = 6) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
  
  chloro <- plot_list[[10]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Chlorophytes") + 
    xlab("Date") + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    annotate(geom = "text", x = dmy("1-1-2014"), y = 0.95, 
             label = "j", fontface = "bold", size = 6) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
  
  meta <- plot_list[[11]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Metazoans") + 
    xlab("Date") + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    annotate(geom = "text", x = dmy("1-1-2014"), y = 0.95, 
             label = "k", fontface = "bold", size = 6) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
  
  plot_grad <-  pro + syn + flavo + rho +
    sar + diatom + dino + sindin +
    hapto + chloro + meta + guide_area() + plot_layout(ncol = 4, guides = "collect")
  
  
  pdf(file = gradient_plot_file, width = 14, height = 10)
  print(plot_grad)
  dev.off()
  
  
}

##### Suppl Figure 16: Diversity vs Slope #####

suppl_fig_16_func <- function(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                                 "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                 "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                               in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales",
                                                  "Rhodobacterales", "SAR 11 Clade", "Diatoms",
                                                  "Dinoflagellates", "Syndiniales", "Haptophytes",
                                                  "Chlorophytes","Metazoans"),
                             gradient_plot_file = "figures_S/supp_fig_16_S.pdf",
                             tsize = 15){
  
  
  plot_list <- list()
  
  for (i in  1:length(in_group_list)) {
    
    load(paste0("output/",in_group_list[i],"_diffs_div_S.Rdata"))
    
    plot_list[[i]] <- phase
    
  }
  
  pro <- plot_list[[1]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(0.19,0.88),
          plot.tag = element_text(size = 14, face = "bold")) +
    ggtitle("Prochlorococcus") + 
    ylab("Mean Alpha Diversity\nPer Cruise") +
    labs(tag = "a")
  
  syn <- plot_list[[2]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(0.1,0.88),
          plot.tag = element_text(size = 14, face = "bold"),
          axis.title.y = element_blank()) +
    ggtitle("Synecococcus")  +
    labs(tag = "b")
  
  flavo <- plot_list[[3]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(0.1,0.88),
          plot.tag = element_text(size = 14, face = "bold"),
          axis.title.y = element_blank()) +
    ggtitle("Flavobacteriales") +
    labs(tag = "c")
  
  rho <- plot_list[[4]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(0.1,0.88),
          plot.tag = element_text(size = 14, face = "bold"),
          axis.title.y = element_blank()) +
    ggtitle("Rhodobacterales") +
    labs(tag = "d") 
  
  sar <- plot_list[[5]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(0.19,0.88),
          plot.tag = element_text(size = 14, face = "bold")) +
    ggtitle("SAR 11 Clade") + 
    ylab("Mean Alpha Diversity\nPer Cruise") +
    labs(tag = "e")
  
  diatom <- plot_list[[6]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(0.1,0.88),
          plot.tag = element_text(size = 14, face = "bold"),
          axis.title.y = element_blank()) +
    ggtitle("Diatoms") +
    labs(tag = "f") 
  
  dino <- plot_list[[7]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(0.1,0.88),
          plot.tag = element_text(size = 14, face = "bold"),
          axis.title.y = element_blank()) +
    ggtitle("Dinoflagellates") +
    labs(tag = "g") 
  
  sindin <- plot_list[[8]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(0.1,0.88),
          plot.tag = element_text(size = 14, face = "bold"),
          axis.title.y = element_blank()) +
    ggtitle("Syndiniales") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    labs(tag = "h")
  
  hapto <- plot_list[[9]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(0.19,0.88),
          plot.tag = element_text(size = 14, face = "bold")) +
    ggtitle("Haptophytes") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") + 
    ylab("Mean Alpha Diversity\nPer Cruise") +
    labs(tag = "i")
  
  chloro <- plot_list[[10]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(0.1,0.88),
          plot.tag = element_text(size = 14, face = "bold"),
          axis.title.y = element_blank()) +
    ggtitle("Chlorophytes") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    labs(tag = "j")
  
  meta <- plot_list[[11]] +
    theme(plot.title = element_text(hjust = 0, size = tsize+4),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          plot.tag.position = c(0.1,0.88),
          plot.tag = element_text(size = 14, face = "bold"),
          axis.title.y = element_blank()) +
    ggtitle("Metazoans") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    labs(tag = "k")
 
  plot_grad <-  pro + syn + flavo + rho +
    sar + diatom + dino + sindin +
    hapto + chloro + meta + guide_area() + plot_layout(ncol = 4, guides = "collect") 
  
  
  pdf(file = gradient_plot_file, width = 19, height = 13)
  print(plot_grad)
  dev.off()
  
  
}

###### Suppl Figure 17: Mock Communities ########

source("analysis/suppl_fig_17_S.R")

##### Suppl Figure XX: Violin Plots ######

suppl_fig_11x_func <- function(in_list = fig_list, file_name = "figures_S/supp_fig_XX_S.pdf", tsize = 12){
  
  groups <- c(1:5,7:12)
  
  total_mat <- bind_rows(in_list[[1]]$bc_mat, in_list[[2]]$bc_mat)
  p_val_comb <- bind_rows(in_list[[1]]$p_val_mat[1:2,], in_list[[2]]$p_val_mat[1:2,])
  
  for (i in 3:11) {
    total_mat <- bind_rows(total_mat, in_list[[groups[i]]]$bc_mat)
    p_val_comb <- bind_rows(p_val_comb, in_list[[groups[i]]]$p_val_mat[1:2,])
  }
  
  total_mat$Group <-  as.factor(total_mat$Group)
  
  total_mat$Group <- factor(total_mat$Group,
                            levels =  c("Prochlorococcus", "Synecococcus", "Flavobacteriales",
                                        "Rhodobacterales", "SAR 11 Clade", "Diatoms",
                                        "Dinoflagellates", "Syndiniales", "Haptophytes",
                                        "Chlorophytes", "Metazoans"))
  
  total_mat$bc_val[is.nan(total_mat$bc_val)] <- NA
  
  p_val_comb$Group <-  as.factor(p_val_comb$Group)
  
  p_val_comb$Group <- factor(p_val_comb$Group,
                             levels =  c("Prochlorococcus", "Synecococcus", "Flavobacteriales",
                                         "Rhodobacterales", "SAR 11 Clade", "Diatoms",
                                         "Dinoflagellates", "Syndiniales", "Haptophytes",
                                         "Chlorophytes", "Metazoans"))
  
  within_means <- total_mat %>%
    filter(Depth == "Surface", Phase != "Between") %>%
    group_by(Group, Phase) %>%
    summarise(mean_bc = mean(bc_val, na.rm = TRUE))
  
  within <- total_mat %>%
    filter(Depth == "Surface", Phase != "Between")
  
  between <- total_mat %>%
    filter(Depth == "Surface", Phase == "Between")
  
  p_val_comb <- p_val_comb %>% filter(p_val < 0.05)
  
  group_order <-  c("Prochlorococcus", "Synecococcus", "Flavobacteriales",
                    "Rhodobacterales", "SAR 11 Clade", "Diatoms",
                    "Dinoflagellates", "Syndiniales", "Haptophytes",
                    "Chlorophytes", "Metazoans")
  
  p_val_comb$x_id <- match(p_val_comb$Group, group_order)
  
  p_val_comb$X_val <- NA
  p_val_comb$Y_val <- NA
  
  median_lines <- as.data.frame(matrix(NA,33,6))
  
  grab_plot1 <- ggplot() +
    geom_boxplot(data = within, aes(x = Group, y = bc_val, fill = Phase),
                 width = 1.5, alpha = 1)
  grab_plot2 <- ggplot() +
    geom_boxplot(data = between, aes(x = Group, y = bc_val, fill = Phase),
                 width = 1.5, alpha = 1)
  grab_plot3 <- ggplot() +
    geom_split_violin(data = within, aes(x = Group, y = bc_val, fill = Phase),
                      width = 1.5, alpha = 1)
  grab_plot4 <- ggplot() +
    geom_violin(data = between, aes(x = Group, y = bc_val, fill = Phase),
                width = 0.7, alpha = 1)
  
  dat_w <- ggplot_build(grab_plot1)$data[[1]]
  dat_w_vio <- ggplot_build(grab_plot3)$data[[1]]
  dat_b <- ggplot_build(grab_plot2)$data[[1]]
  dat_b_vio <- ggplot_build(grab_plot4)$data[[1]]
  
  median_lines$V1 <- c(dat_w$middle, dat_b$middle)
  median_lines$V2 <- c(dat_w$middle, dat_b$middle)
  median_lines$V3 <-  c(rep(1:11, each = 2),1:11)
  median_lines$V4 <-  c(rep(1:11, each = 2),1:11)
  median_lines$V5 <- c(rep("Within",22),rep("Between",11))
  median_lines$V6 <- c(rep(c("#F8766D", "#00BFC4"),11),rep("#F8766D",11))
  
  for (i in 1:nrow(median_lines)) {
    
    if(i < 23){
      subset <- dat_w_vio %>% filter(x == median_lines$V3[i], fill == median_lines$V6[i])
      width <- subset$violinwidth[which(abs(subset$y-median_lines$V1[i])==
                                          min(abs(subset$y-median_lines$V1[i])))] * 0.75
      if(median_lines$V6[i] == "#F8766D"){median_lines$V4[i] <- median_lines$V3[i] - width}
      if(median_lines$V6[i] == "#00BFC4"){median_lines$V4[i] <- median_lines$V3[i] + width}
    }
    
    if(i > 22){
      
      subset <- dat_b_vio %>% filter(x == median_lines$V3[i], fill == median_lines$V6[i])  
      width <- subset$violinwidth[which(abs(subset$y-median_lines$V1[i])==
                                          min(abs(subset$y-median_lines$V1[i])))] * 0.35
      median_lines$V3[i] <- median_lines$V3[i] - width
      median_lines$V4[i] <- median_lines$V4[i] + width
      
    }
    
  }
  
  p_val_comb$color <- NA
  p_val_comb$color[p_val_comb$Phase == "Cool"] <- "#F8766D"
  p_val_comb$color[p_val_comb$Phase == "Warm"] <- "#00BFC4"
  
  for (i in 1:nrow(p_val_comb)) {
    
    y.1 <- median_lines$V1[which(median_lines$V6 == p_val_comb$color[i] &
                                   median_lines$V3 == p_val_comb$x_id[i] &
                                   median_lines$V5 == "Within")]
    
    p_val_comb$Y_val[i] <- y.1
    
    x.1 <- median_lines$V4[which(median_lines$V6 == p_val_comb$color[i] &
                                   median_lines$V3 == p_val_comb$x_id[i] &
                                   median_lines$V5 == "Within")]
    
    subset <- dat_b_vio %>% filter(x == p_val_comb$x_id[i])  
    width <- subset$violinwidth[which(abs(subset$y-y.1)==
                                        min(abs(subset$y-y.1)))] * 0.35
    
    if(p_val_comb$Phase[i] == "Cool"){
      x.2 <- p_val_comb$x_id[i] - width
      p_val_comb$X_val[i] <- mean(c(x.1,x.2))
    }
    
    if(p_val_comb$Phase[i] == "Warm"){
      x.2 <- p_val_comb$x_id[i] + width
      p_val_comb$X_val[i] <- mean(c(x.1,x.2))
    }
    
  }
  
  with_med <- median_lines %>% filter(V5 == "Within")
  b_med <- median_lines %>% filter(V5 == "Between")
  
  fig_plot <- ggplot() +
    geom_split_violin(data = within, aes(x = Group, y = bc_val, fill = Phase),
                      width = 1.5, alpha = 1) +
    annotate("segment", x = with_med$V3, xend = with_med$V4,
             y = with_med$V2, yend = with_med$V1,
             colour = "black") +
    geom_violin(data = between, aes(x = Group, y = bc_val, fill = Phase),
                width = 0.7, alpha = 1) +
    annotate("segment", x = b_med$V3, xend = b_med$V4,
             y = b_med$V2, yend = b_med$V1,
             colour = "black") +
    annotate("point", x = p_val_comb$X_val, y = p_val_comb$Y_val, size = 3) +
    scale_fill_manual(values = c("white", "dodgerblue1", "firebrick1")) +
    ylab("Bray-Curtis Similarity") + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.title = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          axis.text.y = element_text(size = tsize+2),
          axis.text.x = element_text(size = tsize),
          axis.title.x = element_blank(),
          legend.text = element_text(size = tsize +2),
          legend.title = element_text(size = tsize +2)) 
  
  pdf(file = file_name, width = 16.75, height = 6)
  print(fig_plot)
  dev.off()
  
}

##### Suppl Figure XX: Total Diversity vs Slope #####

suppl_fig_15x_func <- function(in_phyto = "output/euks_auto_18sv9_diffs_div_S.Rdata",
                             in_euks = "output/euks_hetero_18sv9_diffs_div_S.Rdata",
                             in_cyano = "output/cyano_16s_diffs_div_S.Rdata",
                             in_bact = "output/bacteria_m_euks_16s_diffs_div_S.Rdata",
                             in_arch = "output/archaea_16s_diffs_div_S.Rdata",
                             gradient_plot_file = "figures_S/supp_fig_XX_S.pdf",
                             tsize = 12){
  
  load(in_phyto)
  phyto_gradient2 <- gradient_plot2
  
  load(in_euks)
  euk_gradient2 <- gradient_plot2
  
  load(in_cyano)
  cyano_gradient2 <- gradient_plot2
  
  load(in_bact)
  bact_gradient2 <- gradient_plot2
  
  load(in_arch)
  arch_gradient2 <- gradient_plot2
  
  phyto_gradient2 <- phyto_gradient2 + theme(legend.position = "none") + xlab("")  + 
    ylab("Gamma Diversity\nPer Cruise") +
    theme(plot.title = element_text(hjust=0, size = tsize + 4),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.key = element_blank(),
          legend.position = "none") + ggtitle("Photosynthetic eukaryotic protists") +
    annotate(geom = "text", x =0.17, y = 5.2, label = "a", fontface ="bold", size = 5)
  
  euk_gradient2 <- euk_gradient2 + theme(legend.position = "none") + xlab("") + ylab("") +
    theme(axis.text  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize + 4),
          axis.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.key = element_blank(),
          legend.position = "none") + ggtitle("Heterotrophic eukaryotic protists") +
    annotate(geom = "text", x =0.17, y = 6.8, label = "b", fontface ="bold", size = 5)
  
  cyano_gradient2 <- cyano_gradient2 + theme(legend.position = "none")  + 
    ylab("Gamma Diversity\nPer Cruise") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    theme(plot.title = element_text(hjust=0, size = tsize + 4),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize),
          legend.key = element_blank(),
          legend.position = "none") + ggtitle("Cyanobacteria") +
    annotate(geom = "text", x =0.17, y = 2.37, label = "c", fontface ="bold", size = 5)
  
  bact_gradient2 <- bact_gradient2 +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize + 4),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize),
          legend.key = element_blank(),
          legend.text = element_text(size = tsize),
          legend.title = element_text(size = tsize)) + ggtitle("Bacteria") +
    ylab("") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    annotate(geom = "text", x =0.17, y = 6.14, label = "d", fontface ="bold", size = 5)
  
  arch_gradient2 <- arch_gradient2 +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize + 4),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize),
          legend.key = element_blank(),
          legend.position = "none") + ggtitle("Archaea") +
    ylab("") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    annotate(geom = "text", x =0.17, y = 4.25, label = "e", fontface ="bold", size = 5)
  
  
  grad_plot <- phyto_gradient2 + euk_gradient2 + guide_area() +
    cyano_gradient2 + bact_gradient2 + arch_gradient2 + plot_layout(guides = "collect")
  
  pdf(file = gradient_plot_file, width = 13, height = 8)
  print(grad_plot)
  dev.off()
  
  
}

###### Suppl Figure XX: Filter QC ########

# source("analysis/suppl_fig_15.R")

##### Suppl Figure XX: Community vs Time Season #####


suppl_fig_xx_func <- function(in_phyto = "output/euks_auto_18sv9_diffs.Rdata",
                       in_euks = "output/euks_hetero_18sv9_diffs.Rdata",
                       in_cyano = "output/cyano_16s_diffs.Rdata",
                       in_bact = "output/bacteria_m_euks_16s_diffs.Rdata",
                       in_arch = "output/archaea_16s_diffs.Rdata",
                       gradient_plot_file = "figures/figure_outline/supp_fig_11.pdf",
                       tsize = 12){
  
  
  load(in_phyto)
  phyto_gradient <- season
  
  load(in_euks)
  euk_gradient <- season
  
  load(in_cyano)
  cyano_gradient <- season
  
  load(in_bact)
  bact_gradient <- season
  
  load(in_arch)
  arch_gradient <- season
  
  phyto_gradient <- phyto_gradient + theme(legend.position = "none") + xlab("")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.title.y = element_text(size = tsize),
          axis.text.y = element_text(size = tsize)) + ggtitle("A. Photosynthetic Eukaryotic Protists")
  
  euk_gradient <- euk_gradient + theme(legend.position = "none") + xlab("") + ylab("") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize)) + ggtitle("B. Heterotrophic Eukaryotic Protists")
  
  cyano_gradient <- cyano_gradient + theme(legend.position = "none")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("C. Cyanobacteria")
  
  bact_gradient <- bact_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("D. Bacteria") +
    ylab("") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
 
  arch_gradient <- arch_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("E. Archaea") +
    ylab("") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
  
  grad_plot <- phyto_gradient + euk_gradient + guide_area() +
    cyano_gradient + bact_gradient + arch_gradient + plot_layout(guides = "collect")
  
  pdf(file = gradient_plot_file, width = 12, height = 7)
  print(grad_plot)
  dev.off()
  
  
}





#### Suppl Figure XX: Diversity v Season #####

suppl_fig_xx_func <- function(in_phyto = "output/euks_auto_18sv9_diffs_div.Rdata",
                       in_euks = "output/euks_hetero_18sv9_diffs_div.Rdata",
                       in_cyano = "output/cyano_16s_diffs_div.Rdata",
                       in_bact = "output/bacteria_m_euks_16s_diffs_div.Rdata",
                       in_arch = "output/archaea_16s_diffs_div.Rdata",
                       gradient_plot_file = "figures/figure_outline/supp_fig_12.pdf",
                       tsize = 12){
  
  load(in_phyto)
  phyto_gradient <- season
  
  load(in_euks)
  euk_gradient <- season
  
  load(in_cyano)
  cyano_gradient <- season
  
  load(in_bact)
  bact_gradient <- season
  
  load(in_arch)
  arch_gradient <- season
  
  phyto_gradient <- phyto_gradient + theme(legend.position = "none") + xlab("")  + 
    ylab("Mean Shannon Diversity\nPer Cruise") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.title.y = element_text(size = tsize),
          axis.text.y = element_text(size = tsize)) + ggtitle("A. Photosynthetic Eukaryotic Protists")
  
  euk_gradient <- euk_gradient + theme(legend.position = "none") + xlab("") + ylab("") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize)) + ggtitle("B. Heterotrophic Eukaryotic Protists")
  
  cyano_gradient <- cyano_gradient + theme(legend.position = "none")  + 
    ylab("Mean Shannon Diversity\nPer Cruise") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("C. Cyanobacteria")
  
  bact_gradient <- bact_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("D. Bacteria") +
    ylab("") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
  
  arch_gradient <- arch_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("E. Archaea") +
    ylab("") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
  
  grad_plot <- phyto_gradient + euk_gradient + guide_area() +
    cyano_gradient + bact_gradient + arch_gradient + plot_layout(guides = "collect")
  
  pdf(file = gradient_plot_file, width = 12, height = 7)
  print(grad_plot)
  dev.off()
  
  
}


########## Additional Figs ######

fig_7_func_2 <- function(in_list = fig_list, file_name = "figures/figure_outline/fig_7_2.pdf", tsize  = 12){
  
  
  all_plots <- (in_list[[18]]$surf + ggtitle("A. All ASVs") +
                  scale_fill_gradientn(colors = c("white", "red", "black"),
                                       values=rescale(c(0.05,0.2,0.5)),
                                       limits = c(0.05,0.5),
                                       oob = scales::squish) +
                  labs(fill = "Bray-Curtis\nSimilarity")) +
    (in_list[[18]]$deep +
       scale_fill_gradientn(colors = c("white", "red", "black"),
                            values=rescale(c(0.05,0.2,0.5)),
                            limits = c(0.05,0.5),
                            oob = scales::squish) +
       labs(fill = "Bray-Curtis\nSimilarity")) +
    (in_list[[7]]$surf + ggtitle("B. Diatoms") +
       scale_fill_gradientn(colors = c("white", "red", "black"),
                            values=rescale(c(0.05,0.2,0.5)),
                            limits = c(0.05,0.5),
                            oob = scales::squish) +
       labs(fill = "Bray-Curtis\nSimilarity")) +
    (in_list[[7]]$deep +
       scale_fill_gradientn(colors = c("white", "red", "black"),
                            values=rescale(c(0.05,0.2,0.5)),
                            limits = c(0.05,0.5),
                            oob = scales::squish) +
       labs(fill = "Bray-Curtis\nSimilarity")) +
    plot_layout(ncol = 2, nrow = 2, byrow = FALSE, guides = "collect")
  
  pdf(file = file_name, width = 10, height = 10)
  print(all_plots)
  dev.off()
  
}

cuti_func <- function(in_phyto = "output/euks_auto_18sv9_diffs.Rdata",
                       in_euks = "output/euks_hetero_18sv9_diffs.Rdata",
                       in_cyano = "output/cyano_16s_diffs.Rdata",
                       in_bact = "output/bacteria_m_euks_16s_diffs.Rdata",
                       in_arch = "output/archaea_16s_diffs.Rdata",
                       gradient_plot_file = "figures/figure_outline/fig_x_cuti_plot.pdf",
                       tsize = 12){
  
  
  load(in_phyto)
  phyto_gradient <- cuti_plot
  
  load(in_euks)
  euk_gradient <- cuti_plot
  
  load(in_cyano)
  cyano_gradient <- cuti_plot
  
  load(in_bact)
  bact_gradient <- cuti_plot
  
  load(in_arch)
  arch_gradient <- cuti_plot
  
  phyto_gradient <- phyto_gradient + theme(legend.position = "none") + xlab("")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.title.y = element_text(size = tsize),
          axis.text.y = element_text(size = tsize)) + ggtitle("A. Photosynthetic Eukaryotic Protists")
  
  euk_gradient <- euk_gradient + theme(legend.position = "none") + xlab("") + ylab("") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize)) + ggtitle("B. Heterotrophic Eukaryotic Protists")
  
  cyano_gradient <- cyano_gradient + theme(legend.position = "none")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") + 
    xlab("Coastal Upwelling Transport Index\n(CUTI)") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("C. Cyanobacteria")
  
  bact_gradient <- bact_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("D. Bacteria") +
    ylab("") + xlab("Coastal Upwelling Transport Index\n(CUTI)")
  
  arch_gradient <- arch_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("E. Archaea") +
    ylab("") + xlab("Coastal Upwelling Transport Index\n(CUTI)")
  
  grad_plot <- phyto_gradient + euk_gradient + guide_area() +
    cyano_gradient + bact_gradient + arch_gradient + plot_layout(guides = "collect")
  
  pdf(file = gradient_plot_file, width = 12, height = 7)
  print(grad_plot)
  dev.off()
  
  
}

beuti_func <- function(in_phyto = "output/euks_auto_18sv9_diffs.Rdata",
                      in_euks = "output/euks_hetero_18sv9_diffs.Rdata",
                      in_cyano = "output/cyano_16s_diffs.Rdata",
                      in_bact = "output/bacteria_m_euks_16s_diffs.Rdata",
                      in_arch = "output/archaea_16s_diffs.Rdata",
                      gradient_plot_file = "figures/figure_outline/fig_x_beuti_plot.pdf",
                      tsize = 12){
  
  
  load(in_phyto)
  phyto_gradient <- beuti_plot
  
  load(in_euks)
  euk_gradient <- beuti_plot
  
  load(in_cyano)
  cyano_gradient <- beuti_plot
  
  load(in_bact)
  bact_gradient <- beuti_plot
  
  load(in_arch)
  arch_gradient <- beuti_plot
  
  phyto_gradient <- phyto_gradient + theme(legend.position = "none") + xlab("")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.title.y = element_text(size = tsize),
          axis.text.y = element_text(size = tsize)) + ggtitle("A. Photosynthetic Eukaryotic Protists")
  
  euk_gradient <- euk_gradient + theme(legend.position = "none") + xlab("") + ylab("") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize)) + ggtitle("B. Heterotrophic Eukaryotic Protists")
  
  cyano_gradient <- cyano_gradient + theme(legend.position = "none")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") + 
    xlab("Biologically Effective Upwelling Transport Index\n(BEUTI)") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("C. Cyanobacteria")
  
  bact_gradient <- bact_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("D. Bacteria") +
    ylab("") + xlab("Biologically Effective Upwelling Transport Index\n(BEUTI)")
  
  arch_gradient <- arch_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("E. Archaea") +
    ylab("") + xlab("Biologically Effective Upwelling Transport Index\n(BEUTI)")
  
  grad_plot <- phyto_gradient + euk_gradient + guide_area() +
    cyano_gradient + bact_gradient + arch_gradient + plot_layout(guides = "collect")
  
  pdf(file = gradient_plot_file, width = 12, height = 7)
  print(grad_plot)
  dev.off()
  
  
}

nitrate_func <- function(in_phyto = "output/euks_auto_18sv9_diffs.Rdata",
                       in_euks = "output/euks_hetero_18sv9_diffs.Rdata",
                       in_cyano = "output/cyano_16s_diffs.Rdata",
                       in_bact = "output/bacteria_m_euks_16s_diffs.Rdata",
                       in_arch = "output/archaea_16s_diffs.Rdata",
                       gradient_plot_file = "figures/figure_outline/fig_x_nitrate_plot.pdf",
                       tsize = 12){
  
  
  load(in_phyto)
  phyto_gradient <- reg_nitrate
  
  load(in_euks)
  euk_gradient <- reg_nitrate
  
  load(in_cyano)
  cyano_gradient <- reg_nitrate
  
  load(in_bact)
  bact_gradient <- reg_nitrate
  
  load(in_arch)
  arch_gradient <- reg_nitrate
  
  phyto_gradient <- phyto_gradient + theme(legend.position = "none") + xlab("")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.title.y = element_text(size = tsize),
          axis.text.y = element_text(size = tsize)) + ggtitle("A. Photosynthetic Eukaryotic Protists")
  
  euk_gradient <- euk_gradient + theme(legend.position = "none") + xlab("") + ylab("") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize)) + ggtitle("B. Heterotrophic Eukaryotic Protists")
  
  cyano_gradient <- cyano_gradient + theme(legend.position = "none")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") + 
    xlab("Regionally Availible Nitrate") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("C. Cyanobacteria")
  
  bact_gradient <- bact_gradient + theme(legend.position = "none")  + 
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("D. Bacteria") +
    ylab("") + xlab("Regionally Availible Nitrate")
  
  arch_gradient <- arch_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("E. Archaea") +
    ylab("") + xlab("Regionally Availible Nitrate")
  
  grad_plot <- phyto_gradient + euk_gradient + guide_area() +
    cyano_gradient + bact_gradient + arch_gradient + plot_layout(guides = "collect")
  
  pdf(file = gradient_plot_file, width = 12, height = 7)
  print(grad_plot)
  dev.off()
  
  
}


line_spec_func <- function(in_phyto = "output/euks_auto_18sv9_line_diffs.Rdata",
                         in_euks = "output/euks_hetero_18sv9_line_diffs.Rdata",
                         in_cyano = "output/cyano_16s_line_diffs.Rdata",
                         in_bact = "output/bacteria_m_euks_16s_line_diffs.Rdata",
                         in_arch = "output/archaea_16s_line_diffs.Rdata",
                         gradient_plot_file = "figures/figure_outline/fig_x_line_plot.pdf",
                         tsize = 12){
  
  
  load(in_phyto)
  phyto_gradient1 <- phase_80
  phyto_gradient2 <- phase_90
  
  load(in_euks)
  euk_gradient1 <- phase_80
  euk_gradient2 <- phase_90
  
  load(in_cyano)
  cyano_gradient1 <- phase_80
  cyano_gradient2 <- phase_90
  
  load(in_bact)
  bact_gradient1 <- phase_80
  bact_gradient2 <- phase_90
  
  load(in_arch)
  arch_gradient1 <- phase_80
  arch_gradient2 <- phase_90
  
  phyto_gradient1 <- phyto_gradient1 + theme(legend.position = "none") 
  phyto_gradient2 <- phyto_gradient2 + theme(legend.position = "none") 
  
  euk_gradient1 <- euk_gradient1 + theme(legend.position = "none") 
  euk_gradient2 <- euk_gradient2 + theme(legend.position = "none") 
  
  cyano_gradient1 <- cyano_gradient1 + theme(legend.position = "none")  
  cyano_gradient2 <- cyano_gradient2 + theme(legend.position = "none")  
  
  bact_gradient1 <- bact_gradient1 + theme(legend.position = "none")  
  bact_gradient2 <- bact_gradient2 + theme(legend.position = "none")
  
  arch_gradient1 <- arch_gradient1 + theme(legend.position = "none")
  arch_gradient2 <- arch_gradient2 
  
  grad_plot <- phyto_gradient1 + euk_gradient1 +
    cyano_gradient1 + bact_gradient1 + arch_gradient1 +
    phyto_gradient2 + euk_gradient2 + cyano_gradient2 +
    bact_gradient2 + arch_gradient2 + plot_layout(ncol = 5, guides = "collect")
  
  pdf(file = gradient_plot_file, width = 20, height = 7)
  print(grad_plot)
  dev.off()
  
  
}

beuti_euks_func <- function(in_phyto = "output/euks_auto_18sv9_diffs.Rdata",
                       in_euks = "output/dino_18sv9_diffs.Rdata",
                       in_cyano = "output/diatom_18sv9_diffs.Rdata",
                       in_bact = "output/hapto_18sv9_diffs.Rdata",
                       in_arch = "output/chloro_18sv9_diffs.Rdata",
                       gradient_plot_file = "figures/figure_outline/fig_x_beuti_euks_plot.pdf",
                       tsize = 12){
  
  
  load(in_phyto)
  phyto_gradient <- beuti_plot
  
  load(in_euks)
  euk_gradient <- beuti_plot
  
  load(in_cyano)
  cyano_gradient <- beuti_plot
  
  load(in_bact)
  bact_gradient <- beuti_plot
  
  load(in_arch)
  arch_gradient <- beuti_plot
  
  phyto_gradient <- phyto_gradient + theme(legend.position = "none") + xlab("")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.title.y = element_text(size = tsize),
          axis.text.y = element_text(size = tsize)) + ggtitle("A. Photosynthetic Eukaryotic Protists")
  
  euk_gradient <- euk_gradient + theme(legend.position = "none") + xlab("") + ylab("") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize)) + ggtitle("B. Dinoflagellates")
  
  cyano_gradient <- cyano_gradient + theme(legend.position = "none")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") + 
    xlab("Biologically Effective Upwelling Transport Index\n(BEUTI)") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("C. Diatoms")
  
  bact_gradient <- bact_gradient + theme(legend.position = "none")  +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("D. Haptophytes") +
    ylab("") + xlab("Biologically Effective Upwelling Transport Index\n(BEUTI)")
  
  arch_gradient <- arch_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("E. Chlorophytes") +
    ylab("") + xlab("Biologically Effective Upwelling Transport Index\n(BEUTI)")
  
  grad_plot <- phyto_gradient + euk_gradient + guide_area() +
    cyano_gradient + bact_gradient + arch_gradient + plot_layout(guides = "collect")
  
  pdf(file = gradient_plot_file, width = 12, height = 7)
  print(grad_plot)
  dev.off()
  
  
}

station_func <- function(in_phyto = "output/euks_auto_18sv9_station_diffs.Rdata",
                       in_euks = "output/euks_hetero_18sv9_station_diffs.Rdata",
                       in_cyano = "output/cyano_16s_station_diffs.Rdata",
                       in_bact = "output/bacteria_m_euks_16s_station_diffs.Rdata",
                       in_arch = "output/archaea_16s_station_diffs.Rdata",
                       gradient_plot_file = "figures/figure_outline/fig_x_station_plot.pdf",
                       tsize = 12){
  
  
  load(in_phyto)
  phyto_gradient <- chl_plot
  
  load(in_euks)
  euk_gradient <- chl_plot
  
  load(in_cyano)
  cyano_gradient <- chl_plot
  
  load(in_bact)
  bact_gradient <- chl_plot
  
  load(in_arch)
  arch_gradient <- chl_plot
  
  phyto_gradient <- phyto_gradient + theme(legend.position = "none") + xlab("")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.title.y = element_text(size = tsize),
          axis.text.y = element_text(size = tsize)) + ggtitle("A. Photosynthetic Eukaryotic Protists")
  
  euk_gradient <- euk_gradient + theme(legend.position = "none") + xlab("") + ylab("") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize)) + ggtitle("B. Heterotrophic Eukaryotic Protists")
  
  cyano_gradient <- cyano_gradient + theme(legend.position = "none")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") + 
    xlab("Mean Chlorophyll") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("C. Cyanobacteria")
  
  bact_gradient <- bact_gradient +
    theme(legend.position = "none",
          axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("D. Bacteria") +
    ylab("") + xlab("Mean Chlorophyll")
  
  arch_gradient <- arch_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("E. Archaea") +
    ylab("") + xlab("Mean Chlorophyll")
  
  grad_plot <- phyto_gradient + euk_gradient + guide_area() +
    cyano_gradient + bact_gradient + arch_gradient + plot_layout(guides = "collect")
  
  pdf(file = gradient_plot_file, width = 12, height = 7)
  print(grad_plot)
  dev.off()
  
  
}

station_2_func <- function(in_phyto = "output/euks_auto_18sv9_station2_diffs.Rdata",
                         in_euks = "output/euks_hetero_18sv9_station2_diffs.Rdata",
                         in_cyano = "output/cyano_16s_station2_diffs.Rdata",
                         in_bact = "output/bacteria_m_euks_16s_station2_diffs.Rdata",
                         in_arch = "output/archaea_16s_station2_diffs.Rdata",
                         gradient_plot_file = "figures/figure_outline/fig_x_station2_plot.pdf",
                         tsize = 12){
  
  
  load(in_phyto)
  phyto_gradient <- chl_plot
  
  load(in_euks)
  euk_gradient <- chl_plot
  
  load(in_cyano)
  cyano_gradient <- chl_plot
  
  load(in_bact)
  bact_gradient <- chl_plot
  
  load(in_arch)
  arch_gradient <- chl_plot
  
  phyto_gradient <- phyto_gradient + theme(legend.position = "none") + xlab("")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.title.y = element_text(size = tsize),
          axis.text.y = element_text(size = tsize)) + ggtitle("A. Photosynthetic Eukaryotic Protists")
  
  euk_gradient <- euk_gradient + theme(legend.position = "none") + xlab("") + ylab("") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize)) + ggtitle("B. Heterotrophic Eukaryotic Protists")
  
  cyano_gradient <- cyano_gradient + theme(legend.position = "none")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") + 
    xlab("Mean Nitracline Depth (m)") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("C. Cyanobacteria")
  
  bact_gradient <- bact_gradient +
    theme(legend.position = "none",
          axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("D. Bacteria") +
    ylab("") + xlab("Mean Nitracline Depth (m)")
  
  arch_gradient <- arch_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("E. Archaea") +
    ylab("") + xlab("Mean Nitracline Depth (m)")
  
  grad_plot <- phyto_gradient + euk_gradient + guide_area() +
    cyano_gradient + bact_gradient + arch_gradient + plot_layout(guides = "collect")
  
  pdf(file = gradient_plot_file, width = 12, height = 7)
  print(grad_plot)
  dev.off()
  
  
}

station_2_euks_func <- function(in_phyto = "output/euks_auto_18sv9_station2_diffs.Rdata",
                           in_euks = "output/dino_18sv9_station2_diffs.Rdata",
                           in_cyano = "output/diatom_18sv9_station2_diffs.Rdata",
                           in_bact = "output/hapto_18sv9_station2_diffs.Rdata",
                           in_arch = "output/chloro_18sv9_station2_diffs.Rdata",
                           gradient_plot_file = "figures/figure_outline/fig_x_station2_euks_plot.pdf",
                           tsize = 12){
  
  
  load(in_phyto)
  phyto_gradient <- chl_plot
  
  load(in_euks)
  euk_gradient <- chl_plot
  
  load(in_cyano)
  cyano_gradient <- chl_plot
  
  load(in_bact)
  bact_gradient <- chl_plot
  
  load(in_arch)
  arch_gradient <- chl_plot
  
  phyto_gradient <- phyto_gradient + theme(legend.position = "none") + xlab("")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.title.y = element_text(size = tsize),
          axis.text.y = element_text(size = tsize)) + ggtitle("A. Photosynthetic Eukaryotic Protists")
  
  euk_gradient <- euk_gradient + theme(legend.position = "none") + xlab("") + ylab("") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize)) + ggtitle("B. Dinoflagellates")
  
  cyano_gradient <- cyano_gradient + theme(legend.position = "none")  + 
    ylab("Proportion of Samples\nIdentified as Nearshore") + 
    xlab("Mean Nitracline Depth (m)") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("C. Diatoms")
  
  bact_gradient <- bact_gradient +
    theme(legend.position = "none",
          axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("D. Haptophytes") +
    ylab("") + xlab("Mean Nitracline Depth (m)")
  
  arch_gradient <- arch_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("E. Chlorophytes") +
    ylab("") + xlab("Mean Nitracline Depth (m)")
  
  grad_plot <- phyto_gradient + euk_gradient + guide_area() +
    cyano_gradient + bact_gradient + arch_gradient + plot_layout(guides = "collect")
  
  pdf(file = gradient_plot_file, width = 12, height = 7)
  print(grad_plot)
  dev.off()
  
  
}
