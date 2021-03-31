source("analysis/figure_functions.R")

###### Main Text Figures #####

#### Figure 1: Physical Conditions and ASV compare ####

fig_1_func <- function(in_vel = "output/uv_velocity_table.Rdata",
                       in_bath = "output/CALCOFI_bathymetry_table.Rdata",
                       in_temp = "output/CALCOFI_temp_tables.Rdata",
                       in_cyano = "output/cyano_16s_map.Rdata",
                       in_venn = "output/venn_diagram.Rdata",
                       fig_name = "figures/figure_outline/fig_1.pdf",
                       tsize = 12, psize = 6){
  
  map <- map_data("world")    
  
  load(in_vel)
  load(in_bath)
  load(in_temp)
  load(in_cyano)
  load(in_venn)
  
  
  stations <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "n_samps"),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "# of Samples", low = "white", high = "red") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize))
  
  
  reduced_uv_table <- uv_table[seq(1,1633,6),]
  elevation_table$value <- abs(elevation_table$value)
  
  vel_bath <- ggplot() + 
    geom_raster(data = elevation_table, aes(x = lon, y = lat, fill = value), interpolate = FALSE) +
    scale_fill_gradient(low = "darkblue", high = "cyan", name = "Depth (m)", trans = 'reverse') +
    metR::geom_vector(data = reduced_uv_table, aes(x = lon, y = lat, dx = Mean_U, dy = Mean_V), 
                      arrow.angle = 15, arrow.type = "open", arrow.length = unit(0.2, "inches"), 
                      pivot = 0,preserve.dir = TRUE, direction = "ccw",
                      min.mag = 0, show.legend = NA, color = "white") +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    labs(x = "Longitude", y = "Latitude", mag = "Speed (m/s)", color = "Depth (m)")  +
    theme(panel.background = element_blank(),
          legend.key = element_rect(fill = "grey"),
          panel.border = element_rect(fill = NA, colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.box = "horizontal",
          legend.title = element_text(size = tsize)) + 
    guides(fill = guide_legend(order = 2),mag = guide_legend(order = 1)) 
  
  sst <- ggplot() + 
    geom_tile(data = coeff_table, aes(x = lon, y = lat, fill = coeff_var), width =0.26, height = 0.26) +
    scale_fill_gradient2(name = "Coeff. Var SST", low = "darkblue", mid = "white", high = "darkred", limits = c(0.09,0.12), oob = squish, midpoint = 0.1066851) + geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank())
  
  sst_mean <- ggplot() + 
    geom_tile(data = mean_table, aes(x = lon, y = lat, fill = Mean), width =0.26, height = 0.26) +
    scale_fill_gradient2(name = "SST Mean (°C)", low = "darkblue", mid = "white", high = "darkred", limits = c(15,18), oob = squish, midpoint = 16.5) +geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) 
  
  nc_depth <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "NC_mean"),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean Nitracline\nDepth (m)", low = "darkblue", high = "cyan", trans = 'reverse') +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,1), 
          legend.position=c(0.97, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) 
  
  
  stations <- stations + theme(axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank())
  
  vel_bath <- vel_bath + theme(axis.title.x=element_blank(),
                               axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),
                               axis.title.y=element_blank(),
                               axis.text.y=element_blank(),
                               axis.ticks.y=element_blank())
  
  sst_mean <- sst_mean + theme(axis.title.y=element_blank(),
                               axis.text.y=element_blank(),
                               axis.ticks.y=element_blank())
  

  
  layout <- "ABEE
             CDEE"
  
  patch <- stations + vel_bath +
    nc_depth + sst_mean +
    wrap_elements(grobTree(venn_plot)) +
    plot_layout(design = layout) 
  
  patch <- patch + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 16))
  
  pdf(fig_name, width = 22, height = 12)
  print(patch)
  dev.off()
  
}

#### Figure 2: SOMs ####

fig_2_func <- function(in_list = plot_list, file_name = "figures/figure_outline/fig_2.pdf", tsize  = 12){
  
  in_list[[6]]$bp <- in_list[[6]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle(expression(paste(bold("a"), " Archaea")))
  in_list[[13]]$bp <- in_list[[13]]$bp + theme(axis.title.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize),
                                               legend.position = "none") + ggtitle(expression(paste(bold("b"), " Bacteria")))
  in_list[[15]]$bp <- in_list[[15]]$bp + theme(axis.title.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize),
                                               legend.position = "none") + ggtitle(expression(paste(bold("c"), " Cyanobacteria")))
  in_list[[17]]$bp <- in_list[[17]]$bp + theme(axis.title.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize),
                                               legend.position = "none") + ggtitle(expression(atop(paste(bold("d"), " Photosynthetic eukaryotic"), "protists")))
  in_list[[16]]$bp <- in_list[[16]]$bp + theme(axis.title.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize)) + ggtitle(expression(atop(paste(bold("e"), " Heterotrophic eukaryotic"),"protists")))
  
  in_list[[6]]$rp <- in_list[[6]]$rp + theme(legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[13]]$rp <- in_list[[13]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[15]]$rp <- in_list[[15]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[17]]$rp <- in_list[[17]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[16]]$rp <- in_list[[16]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0)) + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  
  
  
  pdf(file = file_name, width = 15, height = 6)
  print(in_list[[6]]$bp + in_list[[13]]$bp + in_list[[15]]$bp +
          in_list[[17]]$bp + in_list[[16]]$bp +
          in_list[[6]]$rp + in_list[[13]]$rp + in_list[[15]]$rp +
          in_list[[17]]$rp + in_list[[16]]$rp +
          plot_layout(ncol = 5))
  dev.off()
}

#### Figure 3: Reg and Var Import ####

fig_3_func <- function(file_name = "figures/figure_outline/fig_3.pdf", tsize = 12){
  
  cyano <- regression_figure(glm_file = "output/cyano_16s_glm.Rdata",
                             map_file = "output/cyano_16s_map.Rdata",   
                             figure_name = "figures/glm_plots/cyano_16s_som_",
                             main = "16s Cyanobacteria", cluster1 = "Offshore", cluster2 = "Nearshore",
                             var = "NC_mean", var_name = "Nitracline Depth (m)")
  
  cyano <- cyano + ggtitle("") + 
    theme(legend.justification=c(1,0.5), 
          legend.position=c(0.95, 0.5),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          plot.title = element_text(hjust = 0, size = tsize))
  
  aic_plot <- full_aic_table_figure(in_group_list = c("archaea_16s","bacteria_m_euks_16s", "cyano_16s","euks_auto_18sv9",
                                                      "euks_hetero_18sv9"),
                                    in_group_names = c("Archaea","Bacteria", "Cyanobacteria", "Photosynthetic Eukaryotic\nProtists",
                                                       "Heterotrophic\nEukaryotic Protists"),
                                    minimum_tp = 4, width_plot = 10,
                                    figure_name_2 = paste0("figures/aic_figures/big_group_aic_plot_logit",".pdf"),
                                    title_name = "Variable Importance")
  
  aic_plot <- aic_plot + ggtitle("") + 
    theme(axis.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          plot.title = element_text(hjust = 0, size = tsize))
  
  a <- cyano + aic_plot + plot_layout(widths = c(1,2)) +
    plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 14))
  
  pdf(file = file_name, width = 16, height = 8)
  print(a)
  dev.off()
  
}

#### Figure 4: Diversity Example ####

fig_4_func <- function(file_name = "figures/figure_outline/fig_4.pdf",
                       in_group = c("total","diatom_18sv9"), basic = c("totals","18s_diatom"),
                       name = c("All ASVs", "Diatoms"), tsize = 12){
  
  
  
  all_map <- diveristy_figure(map_file = paste0("output/", in_group[1], "_map.Rdata"),
                          full_dat = paste0("output/", in_group[1], "_full_data.Rdata"),
                          figure_start = paste0("figures/diversity/", in_group[1], "_"),
                          main = name[1])
  
  all_alpha_gamma <- alpha_versus_gamma_figure(full_data_file = paste0("output/", in_group[1], "_full_data.Rdata"),
                                           raw_data_file = paste0("data/", basic[1], ".Rdata"),
                                           map_file = paste0("output/", in_group[1], "_map.Rdata"), minimum_tp = 8,
                                           figure_name = paste0("figures/diversity/", in_group[1], "_alpha_gamma.pdf"),
                                           main = name[1])

  
  map <- diveristy_figure(map_file = paste0("output/", in_group[2], "_map.Rdata"),
                          full_dat = paste0("output/", in_group[2], "_full_data.Rdata"),
                          figure_start = paste0("figures/diversity/", in_group[2], "_"),
                          main = name[2])
  
  alpha_gamma <- alpha_versus_gamma_figure(full_data_file = paste0("output/", in_group[2], "_full_data.Rdata"),
                                           raw_data_file = paste0("data/", basic[2], ".Rdata"),
                                           map_file = paste0("output/", in_group[2], "_map.Rdata"), minimum_tp = 8,
                                           figure_name = paste0("figures/diversity/", in_group[2], "_alpha_gamma.pdf"),
                                           main = name[2])
  
  
  all_map <- all_map + ggtitle("All ASVs") +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(0,1), 
          legend.position=c(0.03, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) + 
    labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) +
    scale_fill_viridis()
  
  all_alpha_gamma <- all_alpha_gamma + ggtitle("All ASVs") +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(1,0.5), 
          legend.position=c(0.97, 0.5),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize)) + 
    labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=5))
  
  
  map <- map + ggtitle("Diatoms") +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(0,1), 
          legend.position=c(0.03, 0.97),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) + 
    labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) +
    scale_fill_viridis()
  
  alpha_gamma <- alpha_gamma + ggtitle("Diatoms") +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(1,0.5), 
          legend.position=c(0.97, 0.5),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize)) + 
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
                                                   col = 27, width_plot = 12, tsize = 16)
  
  
  layout <- "ABEE
             CDEE"
  
  plots <- all_map + all_alpha_gamma + map +
    alpha_gamma + out_plot + 
    plot_layout(design = layout)
    
  plots <- plots + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 16))
  
  
  pdf(file = file_name, width = 28, height = 15)
  print(plots)
  dev.off()
  
}


###### Figure 5: Time Nutrients ect. #####

fig_5_func <- function(in_upwell = "output/upwelling_plots.Rdata",
                       nutrient_dat = "output/mld_mean_profiles.Rdata",
                       time_plot = "figures/figure_outline/fig_5.pdf",
                       tsize = 12, psize = 6,
                       in_list = fig_list){
  
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
    xlab("Year Day") + ylab("CUTI")
  
  yearly_beuti <- yearly_beuti + theme(plot.title = element_text(size = tsize),
                                       axis.text = element_text(size = tsize),
                                       axis.title = element_text(size = tsize),
                                       legend.text = element_text(size = tsize),
                                       legend.title = element_text(size = tsize),
                                       legend.position = "none") +
    ggtitle("Biologically Effective Upwelling\nTransport Index (BEUTI)") +
    xlab("Year Day") + ylab("BEUTI")
  
  yearly_nitrate <- yearly_nitrate + theme(plot.title = element_text(size = tsize),
                                           axis.text = element_text(size = tsize),
                                           axis.title = element_text(size = tsize),
                                           legend.text = element_text(size = tsize),
                                           legend.title = element_text(size = tsize)) +
    ggtitle("Regionally available\n nitrate (BEUTI/CUTI)") +
    ylab("BEUTI/CUTI")
  
  early_phase$n_diff <- early_phase$mean_no3 - late_phase$mean_no3
  early_phase$p_diff <- early_phase$mean_po4 - late_phase$mean_po4
  
  early_n <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_no3),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean NO3 mM", low = "white", high = "red",
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) + 
    ggtitle("Mean nitrate mM (2014-2016)")
  
  late_n <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_no3),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean NO3 mM", low = "white", high = "red",
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Mean nitrate mM (2017-2018)")
  
  diff_n <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = n_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(Delta~"NO3 mM"),
                        low = "blue", high = "white",
                        limits = c(-6,0), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) + 
    ggtitle("Changes in mean nitrate\n(2014-2016) - (2017-2018)")
  
  low <- 0.125
  high <- 0.6
  
  arch <- (in_list[[6]]$surf + ggtitle("Archaea") +
      scale_fill_gradient(limits = c(low,high),
                          low = "white", high = "red", oob = scales::squish) +
      labs(fill = "Bray-Curtis\nSimilarity") +
      theme(legend.position = "none")) 
    
    bact <- (in_list[[13]]$surf + ggtitle("Bacteria") +
       scale_fill_gradient(limits = c(low,high),
                           low = "white", high = "red", oob = scales::squish) +
       labs(fill = "Bray-Curtis\nSimilarity") +
       theme(axis.title.y = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank())+
       theme(legend.position = "none")) 
    
    cyano <- (in_list[[15]]$surf + ggtitle("Cyanobacteria") +
       scale_fill_gradient(limits = c(low,high),
                           low = "white", high = "red", oob = scales::squish) +
       labs(fill = "Bray-Curtis\nSimilarity") +
       theme(axis.title.y = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank())+
       theme(legend.position = "none")) 
    
    photo <- (in_list[[16]]$surf + ggtitle("Photosynthetic Eukaryotic\nProtists") +
       scale_fill_gradient(limits = c(low,high),
                           low = "white", high = "red", oob = scales::squish) +
       labs(fill = "Bray-Curtis\nSimilarity") +
       theme(axis.title.y = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank())+
       theme(legend.position = "none")) 
    
    hetero <- (in_list[[17]]$surf + ggtitle("Heterotrophic Eukaryotic\nProtists") +
       scale_fill_gradient(limits = c(low,high),
                           low = "white", high = "red", oob = scales::squish) +
       labs(fill = "Bray-Curtis\nSimilarity") +
       theme(axis.title.y = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank()))
    
    fig_plot <- arch + bact + cyano + photo + hetero + plot_layout(nrow = 1, guides = "collect")
  
  layout = "ABC
            DEF
            GGG"
  
  plot1 <- yearly_cuti + yearly_beuti + yearly_nitrate  +
    early_n + late_n + diff_n  + fig_plot +
    plot_layout(design = layout)
  
  plot1 <- plot1 + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 16))

  pdf(file = time_plot, width = 17, height = 17)
  print(plot1)
  dev.off()

}


##### Figure 6: Community vs Time #####

fig_6_func <- function(in_phyto = "output/euks_auto_18sv9_diffs.Rdata",
                       in_euks = "output/euks_hetero_18sv9_diffs.Rdata",
                       in_cyano = "output/cyano_16s_diffs.Rdata",
                       in_bact = "output/bacteria_m_euks_16s_diffs.Rdata",
                       in_arch = "output/archaea_16s_diffs.Rdata",
                       gradient_plot_file = "figures/figure_outline/fig_6.pdf",
                       tsize = 12){
  
  
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
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("Cyanobacteria")
  
  bact_gradient <- bact_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("Bacteria") +
    ylab("") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
  
  arch_gradient <- arch_gradient +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("Archaea") +
    ylab("") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
  
  grad_plot <- phyto_gradient + euk_gradient + guide_area() +
    cyano_gradient + bact_gradient + arch_gradient + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 16))
  
  pdf(file = gradient_plot_file, width = 12, height = 7)
  print(grad_plot)
  dev.off()
  
  
}

##### Supplementary Figures #####

##### Suppl Figure 1: Physical Maps #####

suppl_fig_1_func <- function(in_dat = "output/mld_mean_profiles.Rdata",
                             out_plot = "figures/figure_outline/supp_fig_1.pdf",
                             tsize = 12, psize = 6){
  
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
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Mean Temperature (°C)")
  
  t_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_temp), aes(x = lon, y = lat, fill = coeff_temp),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Coeff. Var. Temperature")
  
  sal_m <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat, aes(x = lon, y = lat, fill = mean_sal),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = "PSU", low = "blue", high = "gold2", mid = "white", midpoint = 33.35) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Mean Salinity")
  
  sal_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_sal), aes(x = lon, y = lat, fill = coeff_sal),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Coeff. Var. Salinity")
  
  # Nutrients
  
  no3_m <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat, aes(x = lon, y = lat, fill = mean_no3),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "NO3 mM", low = "white", high = "purple") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Mean NO3")
  
  no3_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_no3), aes(x = lon, y = lat, fill = coeff_no3),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Coeff. Var. NO3")
  
  po4_m <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat, aes(x = lon, y = lat, fill = mean_po4),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "PO4 mM", low = "white", high = "purple") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Mean PO4")
  
  po4_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_po4), aes(x = lon, y = lat, fill = coeff_po4),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Coeff. Var. PO4")
  
  sio4_m <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat, aes(x = lon, y = lat, fill = mean_sio4),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "SiO4 mM", low = "white", high = "purple") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Mean SiO4")
  
  sio4_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_sio4), aes(x = lon, y = lat, fill = coeff_sio4),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Coeff. Var. SiO4")
  
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
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Mean Chl-a")
  
  chl_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_chl), aes(x = lon, y = lat, fill = coeff_chl),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Coeff. Var. Chl-a")
  
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
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Mean Mixed Layer Depth")
  
  mld_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_mld), aes(x = lon, y = lat, fill = coeff_mld),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Coeff. Var. Mixed Layer Depth")
  
  ncd_m <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat, aes(x = lon, y = lat, fill = mean_ncd),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Depth (m)", low = "darkblue", high = "cyan", trans = 'reverse') +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Mean Nitracline Depth")
  
  ncd_c <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = all_dat %>% drop_na(coeff_ncd), aes(x = lon, y = lat, fill = coeff_ncd),
               color = "black", size = psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(sigma/mu), low = "white", high = "black", na.value = "none") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Coeff. Var. Nitracline Depth")
  
  
  plot1 <- t_m + sal_m + mld_m + ncd_m + 
    t_c + sal_c + mld_c + ncd_c + 
    chl_m + no3_m + po4_m + sio4_m +
    chl_c + no3_c + po4_c + sio4_c + plot_layout(ncol = 4)
  
  pdf(file = out_plot, height = 18, width = 18)
  print(plot1)
  dev.off()
  
  }


##### Suppl Figure 2: NCOG, TARA, Polar Endemics #######

suppl_fig_2_func <- function(in_ncog = "data/18sv9_all.Rdata",
                             in_tara = "data/18sv9_tara_oceans.Rdata",
                             in_polar = "data/18sv9_tara_polar.Rdata",
                             output = "figures/figure_outline/supp_fig_2.pdf",
                             output_table = "figures/figure_outline/supp_table_2.pdf"){
  
  # load data
  load(in_ncog)
  load(in_tara)
  load(in_polar)
  
  calcofi <- names(which(colSums(eighteen_s) != 0))
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
  
  mean_groups <- split_taxa %>%
    filter(!is.na(endemic), !is.na(C)) %>%
    group_by(endemic,C) %>% summarise(count = n()) %>% mutate(prop = count/sum(count))
  
  uniquec <- unique(mean_groups$C)
  
  mean_groups$C <- as.factor(mean_groups$C)
  mean_groups$C <- factor(mean_groups$C, levels = uniquec[c(1:24,26:36,25)])
  
  proportion_plot <- ggplot(mean_groups, aes(x = endemic, y = prop, fill = C)) + 
    geom_bar(stat = "identity") + theme_classic() +
    labs(x = "", y = "Proportion", fill = "Division") +
    scale_fill_manual(values = c("#dc7e83","#41c65d","#9f5dd2","#53a024","#c74eaf","#81c954","#5b6ed9",
                                 "#aebe32","#8358a0","#3d933f","#da4e8d","#51c586","#cf3e5e","#54c1ae",
                                 "#d34239","#4ab4dd","#df7032","#6584c7","#db9835","#cf8fd0","#809f39",
                                 "#9c4969","#6ea66a","#ac5030","#30866c","#dc996b","#34753e","#986039",
                                 "#52701e","#c3ab43","#576426","#abb871","#906d24","#7d8028","#807a45", "grey70")) +
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

source("analysis/suppl_fig_3.R")

##### Suppl Figure 4: SOMs small groups #####

suppl_fig_4_func <- function(in_list = plot_list, file_name = "figures/figure_outline/supp_fig_4.pdf", tsize  = 12){
  
  in_list[[1]]$bp <- in_list[[1]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("a Prochlorococcus")
  in_list[[2]]$bp <- in_list[[2]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("b Synechococcus")
  in_list[[3]]$bp <- in_list[[3]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("c Flavobacteriales")
  in_list[[4]]$bp <- in_list[[4]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("d Rhodobacterales")
  in_list[[5]]$bp <- in_list[[5]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize)) + ggtitle("e SAR 11 Clade")
  in_list[[7]]$bp <- in_list[[7]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("f Diatoms")
  in_list[[8]]$bp <- in_list[[8]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("g Dinoflagellates")
  in_list[[9]]$bp <- in_list[[9]]$bp + theme(axis.title.x = element_blank(),
                                             axis.ticks.x = element_blank(),
                                             axis.text.x = element_blank(),
                                             axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0, size = tsize),
                                             legend.position = "none") + ggtitle("h Syndiniales")
  in_list[[10]]$bp <- in_list[[10]]$bp + theme(axis.title.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize),
                                               legend.position = "none") + ggtitle("i Haptophytes")
  in_list[[11]]$bp <- in_list[[11]]$bp + theme(axis.title.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize),
                                               legend.position = "none") + ggtitle("j Chlorophytes")
  in_list[[12]]$bp <- in_list[[12]]$bp + theme(axis.title.x = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0, size = tsize),
                                               legend.position = "none") + ggtitle("k Metazoans")
  
  in_list[[1]]$rp <- in_list[[1]]$rp + theme(legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[2]]$rp <- in_list[[2]]$rp + theme(axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[3]]$rp <- in_list[[3]]$rp + theme(axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[4]]$rp <- in_list[[4]]$rp + theme(axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[5]]$rp <- in_list[[5]]$rp + theme(axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0)) + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  
  in_list[[7]]$rp <- in_list[[7]]$rp + theme(legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[8]]$rp <- in_list[[8]]$rp + theme(axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[9]]$rp <- in_list[[9]]$rp + theme(axis.title.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.text.y = element_blank(),
                                             plot.title = element_text(hjust = 0),
                                             legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[10]]$rp <- in_list[[10]]$rp + theme(axis.title.y = element_blank(),
                                              axis.ticks.y = element_blank(),
                                              axis.text.y = element_blank(),
                                              plot.title = element_text(hjust = 0),
                                              legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[11]]$rp <- in_list[[11]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  in_list[[12]]$rp <- in_list[[12]]$rp + theme(axis.title.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.text.y = element_blank(),
                                               plot.title = element_text(hjust = 0),
                                               legend.position = "none") + 
    scale_x_continuous(breaks= scales::pretty_breaks(n=3))
  
  
  
  
  pdf(file = file_name, width = 16, height = 12)
  print(in_list[[1]]$bp + in_list[[2]]$bp + in_list[[3]]$bp + in_list[[4]]$bp + in_list[[5]]$bp + plot_spacer() +
          in_list[[1]]$rp + in_list[[2]]$rp + in_list[[3]]$rp + in_list[[4]]$rp + in_list[[5]]$rp + guide_area() +
          in_list[[7]]$bp + in_list[[8]]$bp + in_list[[9]]$bp + in_list[[10]]$bp + in_list[[11]]$bp + in_list[[12]]$bp + 
          in_list[[7]]$rp + in_list[[8]]$rp + in_list[[9]]$rp + in_list[[10]]$rp + in_list[[11]]$rp + in_list[[12]]$rp +
          plot_layout(ncol = 6, guides = "collect"))
  dev.off()
}

##### Suppl Figure 5: Variable AIC small groups #####

full_aic_table_figure(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                        "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                        "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                      in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar Clade", 
                                         "Diatoms",
                                         "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"),
                      minimum_tp = 4, width_plot = 16,
                      figure_name_2 = paste0("figures/figure_outline/supp_fig_5",".pdf"),
                      title_name = "", tsize = 12)

#### Suppl Figure 6: Diversity Maps #####

suppl_fig_6_func <- function(file_name = "figures/figure_outline/supp_fig_6.pdf",
                               in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                                 "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                 "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                               in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar Clade", 
                                                  "Diatoms",
                                                  "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"), tsize = 12){
  
  div_map <- list()
  
  for (i in 1:length(in_group_list)) {
    
    div_map[[i]] <- diveristy_figure(map_file = paste0("output/", in_group_list[i], "_map.Rdata"),
                                     full_dat = paste0("output/", in_group_list[i], "_full_data.Rdata"),
                                     figure_start = paste0("figures/diversity/", in_group_list[i], "_"),
                                     main = in_group_names[i])
    
  }
  
  pro <- div_map[[1]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("a Prochlorococcus") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) 
  
  syn <- div_map[[2]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("b Synecococcus") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) 
  
  flavo <- div_map[[3]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("c Flavobacteriales") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) 
  
  rho <- div_map[[4]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("d Rhodobacterales") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) 
  
  sar <- div_map[[5]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("e SAR 11 Clade") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) 
  
  diatom <- div_map[[6]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("f Diatoms") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) 
  
  dino <- div_map[[7]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("g Dinoflagellates") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) 
  
  sindin <- div_map[[8]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("h Syndiniales") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) 
  
  hapto <- div_map[[9]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("i Haptophytes") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) 
  
  chloro <- div_map[[10]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("j Chlorophytes") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) 
  
  meta <- div_map[[11]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("k Metazoans") + labs(fill = "Mean Alpha\nDiversity") +
    scale_x_continuous(breaks= scales::pretty_breaks(n=3)) 
  
  plot_div <-  pro + syn + flavo + rho +
    sar + diatom + dino + sindin +
    hapto + chloro + meta + plot_spacer() + plot_layout(ncol = 4)
  
  
  pdf(file = file_name, width = 19, height = 17)
  print(plot_div)
  dev.off()
  
}


##### Suppl Figure 7: Diversity Importance Small Groups #####

full_aic_table_figure_diversity_sign(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                                       "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                       "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                                     in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar 11 Clade", 
                                                        "Diatoms",
                                                        "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"),
                                     figure_name_2 = "figures/figure_outline/supp_fig_7.pdf",
                                     col = 27, width_plot = 16, tsize = 12)

####### Suppl Figure 8: Diversity - Productivity ########


suppl_fig_8_func <- function(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                               "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                               "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                             in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales",
                                                "Rhodobacterales", "SAR 11 Clade", "Diatoms",
                                                "Dinoflagellates", "Syndiniales", "Haptophytes",
                                                "Chlorophytes","Metazoans"), tsize = 12,
                             out_fig_name = "figures/figure_outline/supp_fig_8.pdf"){
  
  
  out_plot_list <- list()
  
  for (i in 1:length(in_group_list)) {
    
    load(paste0("output/", in_group_list[i],"_full_data.RData"))
    
    prodo_dat <- full_dat %>%
      filter(!is.na(IntC14))
    
    prodo_dat$phase <- substr(prodo_dat$Cruise,1,4)
    
    prodo_dat$phase[which(prodo_dat$phase == "2014" | prodo_dat$phase == "2015" | prodo_dat$phase == "2016")] <- "Warm"
    prodo_dat$phase[which(prodo_dat$phase == "2017" | prodo_dat$phase == "2018")] <- "Cool"
    prodo_dat$phase[which(prodo_dat$phase == "2019")] <- "2019"
    
    prodo_dat$phase <- as.factor(prodo_dat$phase)
    prodo_dat$phase <- factor(prodo_dat$phase, levels = c("Warm", "Cool", "2019"))
    
    gam_dat <- prodo_dat %>%
      filter(phase != "2019")
    
    out_plot_list[[i]] <- ggplot() +
      geom_point(data = prodo_dat, aes(x = IntC14*2, y = richness, color = phase)) + scale_x_log10() + 
      stat_smooth(data = gam_dat, aes(x = IntC14*2, y = richness, color = phase, fill = phase),
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
            plot.title = element_text(size = tsize)) +
      ggtitle(in_group_names[i]) + guides(fill = FALSE)
    
    
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
    guide_area()  + plot_layout(ncol = 4, guides = "collect") + plot_annotation(tag_levels = "a")
  
  pdf(out_fig_name, width = 12, height = 9)
  print(out_plot)
  dev.off()
  
}







##### Suppl Figure 9: Environmental Differences #####

suppl_fig_9_func <- function(nutrient_dat = "output/mld_mean_profiles.Rdata",
                             suppl_plot = "figures/figure_outline/supp_fig_9.pdf",
                             tsize = 12, psize = 6){
  
  map <- map_data("world")
  
  load(nutrient_dat)
  
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
    scale_fill_gradient(name = "Mean Temperature (°C)", low = "white", high = colors[1],
                        limits = c(13,20), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) + 
    ggtitle("Mean temperature (°C) (2014-2016)")
  
  late_t <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_temp),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean Temperature (°C)", low = "white", high = colors[1],
                        limits = c(13,20), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Mean temperature (°C) (2017-2018)")
  
  diff_t <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = t_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(Delta~"Temperature (°C)"), low = "blue", high = "red",
                         mid = "white", midpoint = 0,
                         limits = c(-0.5,2), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) + 
    ggtitle("Changes in temperature\n(2014-2016) - (2017-2018)")
  
  early_sa <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_sal),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean Salinity (PSU)", low = "white", high = colors[2],
                        limits = c(33,34), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) + 
    ggtitle("Mean salinity (PSU) (2014-2016)")
  
  late_sa <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_sal),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean Salinity (PSU)", low = "white", high = colors[2],
                        limits = c(33,34), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Mean salinity (PSU) (2017-2018)")
  
  diff_sa <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = sa_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(Delta~"Salinity (PSU)"), low = "blue", high = "red",
                         mid = "white", midpoint = 0,
                         limits = c(-0.5,0), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) + 
    ggtitle("Changes in salinity\n(2014-2016) - (2017-2018)")
  
  early_ml <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_mld),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean Mixed Layer\nDepth(m)", low = "white", high = colors[3],
                        limits = c(10,50), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) + 
    ggtitle("Mean mixed layer depth (m) (2014-2016)")
  
  late_ml <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_mld),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean Mixed Layer\nDepth(m)", low = "white", high = colors[3],
                        limits = c(10,50), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Mean mixed layer depth (m) (2017-2018)")
  
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
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) + 
    ggtitle("Changes in mean mixed layer depth\n(2014-2016) - (2017-2018)")
  
  early_nc <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_ncd),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean Nitracline\nDepth(m)", low = "white", high = colors[4],
                        limits = c(0,120), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) + 
    ggtitle("Mean nitracline depth (m) (2014-2016)")
  
  late_nc <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_ncd),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean Nitracline\nDepth(m)", low = "white", high = colors[4],
                        limits = c(0,120), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Mean nitracline depth (m) (2017-2018)")
  
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
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) + 
    ggtitle("Changes in mean nitracline depth\n(2014-2016) - (2017-2018)")
  
  
  early_chl <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_chl),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean Chl-a ug/L", low = "white", high = colors[5],
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) + 
    ggtitle("Mean chl-a ug/L (2014-2016)")
  
  late_chl <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_chl),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean Chl-a ug/L", low = "white", high = colors[5],
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Mean chl-a ug/L (2017-2018)")
  
  diff_chl <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = chl_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(Delta~"Chl-a ug/L"), low = "blue", high = "red",
                         mid = "white", midpoint = 0,
                         limits = c(-2,1), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) + 
    ggtitle("Changes in mean chl-a\n(2014-2016) - (2017-2018)")

  early_n <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_no3),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean NO3 uM/L", low = "white", high = colors[6],
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) + 
    ggtitle("Mean nitrate uM/L (2014-2016)")
  
  late_n <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_no3),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean NO3 uM/L", low = "white", high = colors[6],
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Mean nitrate uM/L (2017-2018)")
  
  diff_n <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = n_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = expression(Delta~"NO3 uM/L"),
                        low = "blue", high = "white",
                        limits = c(-6,0), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0), 
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) + 
    ggtitle("Changes in mean nitrate\n(2014-2016) - (2017-2018)")
  
  early_p <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_po4),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean PO4 uM/L", low = "white", high = colors[7],
                        limits = c(0,1), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0),
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Mean phosphate uM/L (2014-2016)")
  
  late_p <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_po4),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean PO4 uM/L", low = "white", high = colors[7],
                        limits = c(0,1), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0),
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Mean Phosphate uM/L (2017-2018)")
  
  diff_p <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = p_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(Delta~"PO4 uM/L"), low = "blue", high = "red",
                         mid = "white", midpoint = 0,
                        limits = c(-0.25,0.25), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0),
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) + ggtitle("Changes in mean phosphate\n(2014-2016) - (2017-2018)")
  
  early_si <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = mean_sio4),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean SiO4 uM/L", low = "white", high = colors[8],
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0),
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Mean silicate uM/L (2014-2016)")
  
  late_si <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = late_phase, aes(x = lon, y = lat, fill = mean_sio4),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient(name = "Mean SiO4 uM/L", low = "white", high = colors[8],
                        limits = c(0,8), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0),
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) +
    ggtitle("Mean silicate uM/L (2017-2018)")
  
  diff_si <- ggplot() +
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") +
    geom_point(data = early_phase, aes(x = lon, y = lat, fill = si_diff),
               color = "black", size =psize, stroke = 0.1, shape = 21) +
    scale_fill_gradient2(name = expression(Delta~"SiO4 uM/L"), low = "blue", high = "red",
                         mid = "white", midpoint = 0,
                         limits = c(-1.5,1.5), oob = scales::squish) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid"),
          plot.title = element_text(size = tsize), axis.line = element_blank(),
          legend.justification=c(1,0),
          legend.position=c(0.97, 0.03),
          legend.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) + ggtitle("Changes in mean silicate\n(2014-2016) - (2017-2018)")
  
  
  plot_comb <- (early_t + late_t + diff_t + plot_layout(ncol = 3, guides = "collect"))/
    (early_sa + late_sa + diff_sa + plot_layout(ncol = 3, guides = "collect"))/
    (early_ml + late_ml + diff_ml + plot_layout(ncol = 3, guides = "collect"))/
    (early_nc + late_nc + diff_nc + plot_layout(ncol = 3, guides = "collect"))|
    (early_chl + late_chl + diff_chl + plot_layout(ncol = 3, guides = "collect"))/
    (early_n + late_n + diff_n + plot_layout(ncol = 3, guides = "collect"))/
    (early_p + late_p + diff_p + plot_layout(ncol = 3, guides = "collect"))/
    (early_si + late_si + diff_si + plot_layout(ncol = 3, guides = "collect"))
  
  pdf(file = suppl_plot, width = 28, height = 16)
  print(plot_comb)
  dev.off()
  
  }

##### Suppl Figure 10: Community Differences ######

suppl_fig_10_func <- function(in_list = fig_list, file_name = "figures/figure_outline/supp_fig_10.pdf", tsize  = 12){
  
    low <- 0.125
    high <- 0.6
    
    all_plots <- (in_list[[6]]$deep + ggtitle("Archaea\nDCM (2014-2016) vs (2017-2018)") +
                    scale_fill_gradient(limits = c(low,high),
                                        low = "white", high = "red", oob = scales::squish) +
                    labs(fill = "Bray-Curtis\nSimilarity") +
                    theme(legend.position = "none")) +
      (in_list[[13]]$deep + ggtitle("Bacteria\nDCM (2014-2016) vs (2017-2018)") +
         scale_fill_gradient(limits = c(low,high),
                             low = "white", high = "red", oob = scales::squish) +
         labs(fill = "Bray-Curtis\nSimilarity") +
         theme(axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())+
         theme(legend.position = "none")) +
      (in_list[[15]]$deep + ggtitle("Cyanobacteria\nDCM (2014-2016) vs (2017-2018)") +
         scale_fill_gradient(limits = c(low,high),
                             low = "white", high = "red", oob = scales::squish) +
         labs(fill = "Bray-Curtis\nSimilarity") +
         theme(axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())+
         theme(legend.position = "none")) +
      (in_list[[16]]$deep + ggtitle("Photosynthetic eukaryotic\nprotists\nDCM (2014-2016) vs (2017-2018)") +
         scale_fill_gradient(limits = c(low,high),
                             low = "white", high = "red", oob = scales::squish) +
         labs(fill = "Bray-Curtis\nSimilarity") +
         theme(axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())+
         theme(legend.position = "none")) +
      (in_list[[17]]$deep + ggtitle("Heterotrophic eukaryotic\nprotists\nDCM (2014-2016) vs (2017-2018)") +
         scale_fill_gradient(limits = c(low,high),
                             low = "white", high = "red", oob = scales::squish) +
         labs(fill = "Bray-Curtis\nSimilarity") +
         theme(axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())) +
      plot_layout(ncol = 5, nrow = 1, byrow = TRUE, guides = "collect") + 
      plot_annotation(tag_levels = "a")
    
    pdf(file = file_name, width = 17, height = 4)
    print(all_plots)
    dev.off()
  
}

##### Suppl Figure 11: Violin Plots ######

suppl_fig_11_func <- function(in_list = fig_list, file_name = "figures/figure_outline/supp_fig_11.pdf", tsize = 12){

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
                                      "Rhodobacterales", "Sar Clade", "Diatoms",
                                      "Dinoflagellates", "Syndiniales", "Haptophytes",
                                      "Chlorophytes", "Metazoans"))

total_mat$bc_val[is.nan(total_mat$bc_val)] <- NA

p_val_comb$Group <-  as.factor(p_val_comb$Group)

p_val_comb$Group <- factor(p_val_comb$Group,
                           levels =  c("Prochlorococcus", "Synecococcus", "Flavobacteriales",
                                       "Rhodobacterales", "Sar Clade", "Diatoms",
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
                  "Rhodobacterales", "Sar Clade", "Diatoms",
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
  annotate("text", x = 10.5, y = 0.95,
           label = "Wilcoxon Signed Rank Test, p < 0.05") +
  scale_fill_manual(values = c("white", "dodgerblue1", "firebrick1")) +
  ylab("Bray-Curtis Similarity") + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        plot.title = element_text(size = tsize),
        axis.title = element_text(size = tsize),
        axis.text.y = element_text(size = tsize),
        axis.text.x = element_text(size = tsize)) 

pdf(file = file_name, width = 17, height = 6)
print(fig_plot)
dev.off()

}

###### Suppl Figure 12: Prop Nearshore over Time #####



suppl_fig_12_func <- function(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                                 "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                 "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                               in_group_names = c("a Prochlorococcus", "b Synecococcus", "c Flavobacteriales",
                                                  "d Rhodobacterales", "e SAR 11 Clade", "f Diatoms",
                                                  "g Dinoflagellates", "h Syndiniales", "i Haptophytes",
                                                  "j Chlorophytes","k Metazoans"),
                               gradient_plot_file = "figures/figure_outline/supp_fig_12.pdf",
                               tsize = 12){
  
  
  plot_list <- list()
  
  for (i in  1:length(in_group_list)) {
    
    load(paste0("output/",in_group_list[i],"_diffs_div.Rdata"))
    
    plot_list[[i]] <- ts_plot
    
  }
  
  
  pro <- plot_list[[1]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("A. Prochlorococcus") + 
    ylab("Proportion of Samples\nIdentified as Nearshore")
  
  syn <- plot_list[[2]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("B. Synecococcus") 
  
  flavo <- plot_list[[3]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("C. Flavobacteriales") 
  
  rho <- plot_list[[4]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("D. Rhodobacterales") 
  
  sar <- plot_list[[5]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("E. Sar Clade") + 
    ylab("Proportion of Samples\nIdentified as Nearshore")
  
  diatom <- plot_list[[6]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("F. Diatoms") 
  
  dino <- plot_list[[7]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("G. Dinoflagellates") 
  
  sindin <- plot_list[[8]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("H. Syndiniales") + 
    xlab("Date") 
  
  hapto <- plot_list[[9]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("I. Haptophytes") + 
    xlab("Date") + 
    ylab("Proportion of Samples\nIdentified as Nearshore")
  
  chloro <- plot_list[[10]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("J. Chlorophytes") + 
    xlab("Date")
  
  meta <- plot_list[[11]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("K. Metazoans") + 
    xlab("Date")
  
  plot_grad <-  pro + syn + flavo + rho +
    sar + diatom + dino + sindin +
    hapto + chloro + meta + guide_area() + plot_layout(ncol = 4, guides = "collect")
  
  
  pdf(file = gradient_plot_file, width = 14, height = 10)
  print(plot_grad)
  dev.off()
  
  
}

##### Suppl Figure 13: Diversity vs Slope #####

suppl_fig_13_func <- function(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                                                 "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                                                 "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
                               in_group_names = c("a Prochlorococcus", "b Synecococcus", "c Flavobacteriales",
                                                  "d Rhodobacterales", "e Sar Clade", "f Diatoms",
                                                  "g Dinoflagellates", "h Syndiniales", "i Haptophytes",
                                                  "j Chlorophytes","k Metazoans"),
                             gradient_plot_file = "figures/figure_outline/supp_fig_13.pdf",
                             tsize = 12){
  
  
  plot_list <- list()
  
  for (i in  1:length(in_group_list)) {
    
    load(paste0("output/",in_group_list[i],"_diffs_div.Rdata"))
    
    plot_list[[i]] <- phase
    
  }
  
  
  pro <- plot_list[[1]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("Prochlorococcus") + 
    ylab("Mean Alpha Diversity\nPer Cruise")
  
  syn <- plot_list[[2]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("Synecococcus") 
  
  flavo <- plot_list[[3]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("Flavobacteriales") 
  
  rho <- plot_list[[4]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("Rhodobacterales") 
  
  sar <- plot_list[[5]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    ggtitle("SAR 11 Clade") + 
    ylab("Mean Alpha Diversity\nPer Cruise")
  
  diatom <- plot_list[[6]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("Diatoms") 
  
  dino <- plot_list[[7]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("Dinoflagellates") 
  
  sindin <- plot_list[[8]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("Syndiniales") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
  
  hapto <- plot_list[[9]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize)) +
    ggtitle("Haptophytes") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") + 
    ylab("Mean Alpha Diversity\nPer Cruise")
  
  chloro <- plot_list[[10]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          legend.position = "none",
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("Chlorophytes") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
  
  meta <- plot_list[[11]] +
    theme(plot.title = element_text(hjust = 0, size = tsize),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ggtitle("Metazoans") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
 
  plot_grad <-  pro + syn + flavo + rho +
    sar + diatom + dino + sindin +
    hapto + chloro + meta + guide_area() + plot_layout(ncol = 4, guides = "collect") +
    plot_annotation(tag_levels = "a")
  
  
  pdf(file = gradient_plot_file, width = 14, height = 10)
  print(plot_grad)
  dev.off()
  
  
}

##### Suppl Figure 14: Total Diversity vs Slope #####

suppl_fig_14_func <- function(in_phyto = "output/euks_auto_18sv9_diffs_div.Rdata",
                             in_euks = "output/euks_hetero_18sv9_diffs_div.Rdata",
                             in_cyano = "output/cyano_16s_diffs_div.Rdata",
                             in_bact = "output/bacteria_m_euks_16s_diffs_div.Rdata",
                             in_arch = "output/archaea_16s_diffs_div.Rdata",
                             gradient_plot_file = "figures/figure_outline/supp_fig_14.pdf",
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
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.title.y = element_text(size = tsize),
          axis.text.y = element_text(size = tsize)) + ggtitle("Photosynthetic eukaryotic protists")
  
  euk_gradient2 <- euk_gradient2 + theme(legend.position = "none") + xlab("") + ylab("") +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize)) + ggtitle("Heterotrophic eukaryotic protists")
  
  cyano_gradient2 <- cyano_gradient2 + theme(legend.position = "none")  + 
    ylab("Gamma Diversity\nPer Cruise") + 
    xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)") +
    theme(plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("Cyanobacteria")
  
  bact_gradient2 <- bact_gradient2 +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("Bacteria") +
    ylab("") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
  
  arch_gradient2 <- arch_gradient2 +
    theme(axis.text.y  = element_text(size = tsize),
          plot.title = element_text(hjust=0, size = tsize),
          axis.title = element_text(size = tsize),
          axis.text = element_text(size = tsize)) + ggtitle("Archaea") +
    ylab("") + xlab("Nearshore-Offshore\nSlope in Nitracline (m/km)")
  
  grad_plot <- phyto_gradient2 + euk_gradient2 + guide_area() +
    cyano_gradient2 + bact_gradient2 + arch_gradient2 + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a")
  
  pdf(file = gradient_plot_file, width = 12, height = 7)
  print(grad_plot)
  dev.off()
  
  
}

###### Suppl Figure 15: Filter QC ########

source("analysis/suppl_fig_15.R")

###### Suppl Figure 15: Mock Communities ########

source("analysis/suppl_fig_16.R")

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
