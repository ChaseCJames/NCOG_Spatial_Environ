library(mvpart)
library(MVPARTwrap)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(labdsv)
library(RColorBrewer)

multivariate_RT_analysis <- function(raw_data = "data/16s_cyanos.Rdata",
                                     env_data = "output/cyano_16s_full_data.Rdata",
                                     mvrt_output = "output/cyano_16s_mvrt.Rdata",
                                     mvrt_som_plot = "figures/cyano_16s_mvrt_som_compare.pdf",
                                     mvrt_som_plot2 = "figures/cyano_16s_mvrt_som_map.pdf",
                                     som_order = c(2,1), mvrt_order = c(1,2),
                                     title_name = "Cyanobacteria All Samples"){
  
  # load in datases for mvrt analysis
  load(raw_data)
  load(env_data)
  
  if(length(which(is.na(match(rownames(scaled_inputs), full_dat$eco_name)))) > 0){
  scaled_inputs <- scaled_inputs[-which(is.na(match(rownames(scaled_inputs), full_dat$eco_name))),]
  }
  # reorder env dataset
  ordered_dat <- full_dat[match(rownames(scaled_inputs),full_dat$eco_name),]
  
  # only pull explanatory variables
  env_dat <- ordered_dat[,c(33:34,37:39,47:48)]
  
  par(mfrow = c(1,1))
  mvrt_result <- mvpart(scaled_inputs ~.,
                        env_dat,
                        size = 2,
                        margin = 0.8,
                        xv = "none",
                        xval = 10,
                        xvmult = 1,
                        cp = 0)

  
  if(length(mvrt_result$na.action) > 0){final_dat <- ordered_dat[-mvrt_result$na.action,]}else{final_dat <- ordered_dat}
  
  final_dat$mvrt_id <- mvrt_result$where
  
  if(length(table(final_dat$mvrt_id)) > 2){
    
  level <- names(table(final_dat$mvrt_id)[order(table(final_dat$mvrt_id))][1])
  
  final_dat <- final_dat[-which(final_dat$mvrt_id == level),]
  
  
    
  }
  
  tax_dat <- final_dat
  final_dat <- final_dat[-which(as.numeric(substr(final_dat$Sta_ID,2,3)) < 76),]
    
  mvrt_plot <- ggplot(final_dat, aes(x = T_degC, y = NCDepth, color = as.factor(mvrt_id))) + geom_point() +
    scale_color_manual(labels = c("Cluster 1","Cluster 2"), values = c("red","blue")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA)) + 
    labs(color = "Multivariate\nRegression Tree\nClusters") + xlab("Temperature (°C)") +
    ylab("Nutricline Depth (m)") + scale_y_reverse()
  
  
  som_plot <- ggplot(final_dat, aes(x = T_degC, y = NCDepth, color = as.factor(som_id))) + geom_point() +
    scale_color_manual(labels = c("Cluster 1","Cluster 2"), values = c("orange","purple")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA)) + 
    labs(color = "Self-Organizing\nMap Clusters") + xlab("Temperature (°C)") + 
    ylab("Nutricline Depth (m)") + scale_y_reverse()
  
  # ts diagram 
  
  mvrt_plot_ts <- ggplot(final_dat, aes(x = T_degC, y = Salnty, color = as.factor(mvrt_id))) + geom_point() +
    scale_color_manual(labels = c("Cluster 1","Cluster 2"), values = c("red","blue")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA)) + 
    labs(color = "Multivariate\nRegression Tree\nClusters") + xlab("Temperature (°C)") +
    ylab("Salinity") 
  
  
  som_plot_ts <- ggplot(final_dat, aes(x = T_degC, y = Salnty, color = as.factor(som_id))) + geom_point() +
    scale_color_manual(labels = c("Cluster 1","Cluster 2"), values = c("orange","purple")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA)) + 
    labs(color = "Self-Organizing\nMap Clusters") + xlab("Temperature (°C)") + 
    ylab("Salinity") 
  
  # map
  
  map <- map_data("world")
  
  summary_maps <- final_dat %>% 
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som_id == 1, na.rm = TRUE)/n(), som_2 = sum(som_id == 2, na.rm = TRUE)/n(),
              mvrt_1 = sum(mvrt_id == 2, na.rm = TRUE)/n(), mvrt_2 = sum(mvrt_id == 3, na.rm = TRUE)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE))
  
  som_color_vector = c("darkred", "darkblue")
  som_name_vector = c("Nearshore", "Offshore")
  
  mvrt_map_1 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = summary_maps, aes_string(x = "long", y = "lat", fill = paste0("mvrt_",mvrt_order[1])), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = som_color_vector[1], limits = c(0,1)) +
    ggtitle(paste0("MVRT ",som_name_vector[1]," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          title = element_text(hjust = 0.5), axis.line = element_blank())
  
  som_map_1 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = summary_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",som_order[1])), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = som_color_vector[1], limits = c(0,1)) +
    ggtitle(paste0("SOM ",som_name_vector[1]," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          title = element_text(hjust = 0.5), axis.line = element_blank())
  
  mvrt_map_2 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = summary_maps, aes_string(x = "long", y = "lat", fill = paste0("mvrt_",mvrt_order[2])), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = som_color_vector[2], limits = c(0,1)) +
    ggtitle(paste0("MVRT ",som_name_vector[2]," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          title = element_text(hjust = 0.5), axis.line = element_blank())
  
  som_map_2 <-  ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = summary_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",som_order[2])), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = som_color_vector[2], limits = c(0,1)) +
    ggtitle(paste0("SOM ",som_name_vector[2]," Cluster")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          title = element_text(hjust = 0.5), axis.line = element_blank())
  
  
  title <- ggdraw() + draw_label(title_name, fontface='bold')
  
  pdf(mvrt_som_plot, width = 10, height = 6)
  print(plot_grid(title,
                  plot_grid(mvrt_plot,mvrt_plot_ts,som_plot,som_plot_ts, ncol = 2),
        ncol = 1, rel_heights = c(0.1,1)))
  
  dev.off()
  
  pdf(mvrt_som_plot2, width = 12, height = 12)
  print(plot_grid(title,
                  plot_grid(mvrt_map_1, mvrt_map_2, som_map_1, som_map_2, ncol = 2),
                  ncol = 1, rel_heights = c(0.1,1)))
  dev.off()

  
  if(length(mvrt_result$na.action) > 0){scaled_inputs <- scaled_inputs[-mvrt_result$na.action,]}
  if(length(which(colSums(scaled_inputs) == 0)) > 0){scaled_inputs <- scaled_inputs[,-which(colSums(scaled_inputs) == 0)]}
  
  
  ind_val_out <- indval(scaled_inputs, mvrt_result$where)
  pval.adj3 <- p.adjust(ind_val_out$pval)
  
  indicator_spp <- ind_val_out$maxcls[which(ind_val_out$pval <= 0.05)]
  
  indicator_spp <- indicator_spp[-which(indicator_spp == level)]
  
  save(mvrt_result, indicator_spp, final_dat, tax_dat, file = mvrt_output)
  
}

indicator_spp_identifier <- function(mvrt_dat = "output/euks_auto_18sv9_surf_mvrt.Rdata",
                                     taxa_id = "data/18s_autotrophic_euks.Rdata",
                                     amplicon = "eight", database = "PR2", level = 4,
                                     ind_plot = "figures/euks_auto_surf_ind.pdf",
                                     title_name = "Eukaryotic Phytoplankton\nSurface Samples",
                                     plt_height = 6){
  
  load(mvrt_dat)
  load(taxa_id)  
  
  names(indicator_spp)[which(substr(names(indicator_spp),1,1) == "X")] <- substr(names(indicator_spp[which(substr(names(indicator_spp),1,1) == "X")]),2,500) 
  
  if(amplicon == "eight"){
    if(database == "PR2"){taxa <- eight_tax_id$Taxon_PR2[match(names(indicator_spp),eight_tax_id$Feature.ID)]
      }else{taxa <- eight_tax_id$Taxon_Silva[match(names(indicator_spp),eight_tax_id$Feature.ID)]}
  }
 
  
  if(amplicon == "six"){
    if(database == "PR2"){taxa <- six_tax_id$Phytoref_Taxon[match(names(indicator_spp),six_tax_id$Feature.ID)]
    }else{taxa <- six_tax_id$Silva_Taxon[match(names(indicator_spp),six_tax_id$Feature.ID)]}
  }
  
  split_tax <- strsplit(taxa, ";")
  
  taxa_matrix <- as.data.frame(matrix(NA,length(indicator_spp),2))
  
  colnames(taxa_matrix) <- c("Group", "ID")
  
  taxa_matrix$ID <- indicator_spp
  
  tax_vect <- vector()
  for (i in 1:length(split_tax)){tax_vect[i] <- split_tax[[i]][level]}
  
  taxa_matrix$Group <- tax_vect
  
  taxa_table <- as.data.frame(table(taxa_matrix))
  
  colourCount = length(unique(taxa_table$Group))
  getPalette = colorRampPalette(brewer.pal(9, "Set2"))
  
  p1 <- ggplot(taxa_table, aes(x = ID, y = Freq ,fill = Group)) + geom_bar(stat = "identity") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA)) + xlab("MVRT Cluster") +
    ylab("Frequency") + scale_x_discrete(labels = c("Cluster 1","Cluster 2")) +
    scale_fill_manual(values = getPalette(colourCount))
  
  p2 <- ggplot(final_dat, aes(x = T_degC, y = NCDepth, color = as.factor(mvrt_id))) + 
    geom_point() + scale_y_reverse() + 
    scale_color_discrete(labels = c("Cluster 1","Cluster 2")) +
    labs(color = "MVRT Clusters") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA)) +
    xlab("Temperature (°C)") + ylab("Nitricline Depth")
  
  p3 <- ggplot(final_dat, aes(x = T_degC, y = MLD_Sigma-NCDepth, color = as.factor(mvrt_id))) + 
    geom_point() +  
    scale_color_discrete(labels = c("Cluster 1","Cluster 2")) +
    labs(color = "MVRT Clusters") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA)) +
    xlab("Temperature (°C)") + ylab("Mixed Layer Depth -\nNitricline Depth") +
    geom_hline(yintercept = 0, lty = 2)
  
  print(p3)
  
  title <- ggdraw() + draw_label(title_name, fontface='bold')
  
  pdf(ind_plot, width = 10, height = plt_height)
  print(plot_grid(title, plot_grid(p1,p2,ncol = 2), ncol = 1, rel_heights = c(0.1,1)))
  dev.off()
  
  }




multivariate_RT_analysis(raw_data = "data/16s_cyanos.Rdata",
                         env_data = "output/cyano_16s_surf_full_data.Rdata",
                         mvrt_output = "output/cyano_16s_surf_mvrt.Rdata",
                         mvrt_som_plot = "figures/cyano_16s_surf_mvrt_som_compare.pdf",
                         mvrt_som_plot2 = "figures/cyano_16s_surf_mvrt_som_map.pdf",
                         som_order = c(2,1), mvrt_order = c(1,2),
                         title_name = "Cyanobacteria Surface Samples")

multivariate_RT_analysis(raw_data = "data/16s_cyanos.Rdata",
                         env_data = "output/cyano_16s_deep_full_data.Rdata",
                         mvrt_output = "output/cyano_16s_deep_mvrt.Rdata",
                         mvrt_som_plot = "figures/cyano_16s_deep_mvrt_som_compare.pdf",
                         mvrt_som_plot2 = "figures/cyano_16s_deep_mvrt_som_map.pdf",
                         som_order = c(2,1), mvrt_order = c(1,2),
                         title_name = "Cyanobacteria Deep Samples")

multivariate_RT_analysis(raw_data = "data/16s_bacteria_m_euks.Rdata",
                         env_data = "output/bacteria_m_euks_16s_surf_full_data.Rdata",
                         mvrt_output = "output/bacteria_m_euks_16s_surf_mvrt.Rdata",
                         mvrt_som_plot = "figures/bacteria_m_euks_16s_surf_mvrt_som_compare.pdf",
                         mvrt_som_plot2 = "figures/bacteria_m_euks_16s_surf_mvrt_som_map.pdf",
                         som_order = c(1,2), mvrt_order = c(1,2),
                         title_name = "Bacteria/Archaea Surface Samples")

multivariate_RT_analysis(raw_data = "data/16s_bacteria_m_euks.Rdata",
                         env_data = "output/bacteria_m_euks_16s_deep_full_data.Rdata",
                         mvrt_output = "output/bacteria_m_euks_16s_deep_mvrt.Rdata",
                         mvrt_som_plot = "figures/bacteria_m_euks_16s_deep_mvrt_som_compare.pdf",
                         mvrt_som_plot2 = "figures/bacteria_m_euks_16s_deep_mvrt_som_map.pdf",
                         som_order = c(1,2), mvrt_order = c(2,1),
                         title_name = "Bacteria/Archaea Deep Samples")

multivariate_RT_analysis(raw_data = "data/18s_autotrophic_euks.Rdata",
                         env_data = "output/euks_auto_18sv9_surf_full_data.Rdata",
                         mvrt_output = "output/euks_auto_18sv9_surf_mvrt.Rdata",
                         mvrt_som_plot = "figures/euks_auto_18sv9_surf_mvrt_som_compare.pdf",
                         mvrt_som_plot2 = "figures/euks_auto_18sv9_surf_mvrt_som_map.pdf",
                         som_order = c(2,1), mvrt_order = c(1,2),
                         title_name = "Autotrophic Eukaryotes Surface Samples")

multivariate_RT_analysis(raw_data = "data/18s_autotrophic_euks.Rdata",
                         env_data = "output/euks_auto_18sv9_deep_full_data.Rdata",
                         mvrt_output = "output/euks_auto_18sv9_deep_mvrt.Rdata",
                         mvrt_som_plot = "figures/euks_auto_18sv9_deep_mvrt_som_compare.pdf",
                         mvrt_som_plot2 = "figures/euks_auto_18sv9_deep_mvrt_som_map.pdf",
                         som_order = c(2,1), mvrt_order = c(1,2),
                         title_name = "Autotrophic Eukaryotes Deep Samples")

multivariate_RT_analysis(raw_data = "data/18s_heterotrophic_euks.Rdata",
                         env_data = "output/euks_hetero_18sv9_surf_full_data.Rdata",
                         mvrt_output = "output/euks_hetero_18sv9_surf_mvrt.Rdata",
                         mvrt_som_plot = "figures/euks_hetero_18sv9_surf_mvrt_som_compare.pdf",
                         mvrt_som_plot2 = "figures/euks_hetero_18sv9_surf_mvrt_som_map.pdf",
                         som_order = c(2,1), mvrt_order = c(1,2),
                         title_name = "Heterotrophic Eukaryotes Surface Samples")

multivariate_RT_analysis(raw_data = "data/18s_heterotrophic_euks.Rdata",
                         env_data = "output/euks_hetero_18sv9_deep_full_data.Rdata",
                         mvrt_output = "output/euks_hetero_18sv9_deep_mvrt.Rdata",
                         mvrt_som_plot = "figures/euks_hetero_18sv9_deep_mvrt_som_compare.pdf",
                         mvrt_som_plot2 = "figures/euks_hetero_18sv9_deep_mvrt_som_map.pdf",
                         som_order = c(2,1), mvrt_order = c(2,1),
                         title_name = "Heterotrophic Eukaryotes Deep Samples")

multivariate_RT_analysis(raw_data = "data/18s_autotrophic_euks.Rdata",
                         env_data = "output/euks_auto_18sv9_full_data.Rdata",
                         mvrt_output = "output/euks_auto_18sv9_mvrt.Rdata",
                         mvrt_som_plot = "figures/euks_auto_18sv9_mvrt_som_compare.pdf",
                         mvrt_som_plot2 = "figures/euks_auto_18sv9_mvrt_som_map.pdf",
                         som_order = c(2,1), mvrt_order = c(1,2),
                         title_name = "Autotrophic Eukaryotes All Samples")


###### indicator species plots ######


indicator_spp_identifier(mvrt_dat = "output/euks_auto_18sv9_surf_mvrt.Rdata",
                         taxa_id = "data/18s_autotrophic_euks.Rdata",
                         amplicon = "eight", database = "PR2", level = 4,
                         ind_plot = "figures/euks_auto_surf_ind.pdf",
                         title_name = "Eukaryotic Phytoplankton\nSurface Samples",
                         plt_height = 6)

indicator_spp_identifier(mvrt_dat = "output/euks_auto_18sv9_deep_mvrt.Rdata",
                         taxa_id = "data/18s_autotrophic_euks.Rdata",
                         amplicon = "eight", database = "PR2", level = 4,
                         ind_plot = "figures/euks_auto_deep_ind.pdf",
                         title_name = "Eukaryotic Phytoplankton\nDeep Samples",
                         plt_height = 6)

indicator_spp_identifier(mvrt_dat = "output/cyano_16s_surf_mvrt.Rdata",
                         taxa_id = "data/16s_cyanos.Rdata",
                         amplicon = "six", database = "Silva", level = 6,
                         ind_plot = "figures/cyano_surf_ind.pdf",
                         title_name = "Cyanobacteria\nSurface Samples",
                         plt_height = 6)

indicator_spp_identifier(mvrt_dat = "output/cyano_16s_deep_mvrt.Rdata",
                         taxa_id = "data/16s_cyanos.Rdata",
                         amplicon = "six", database = "Silva", level = 6,
                         ind_plot = "figures/cyano_deep_ind.pdf",
                         title_name = "Cyanobacteria\nDeep Samples",
                         plt_height = 6)

indicator_spp_identifier(mvrt_dat = "output/euks_auto_18sv9_mvrt.Rdata",
                         taxa_id = "data/18s_autotrophic_euks.Rdata",
                         amplicon = "eight", database = "PR2", level = 4,
                         ind_plot = "figures/euks_auto_ind.pdf",
                         title_name = "Eukaryotic Phytoplankton\nAll Samples",
                         plt_height = 6)
