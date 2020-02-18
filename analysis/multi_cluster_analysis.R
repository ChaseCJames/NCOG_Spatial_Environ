library(mvpart)
library(MVPARTwrap)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(labdsv)
library(RColorBrewer)
library(pryr)
library(SOMbrero)
library(vegan)

multicluster_analysis <- function(raw_data = "data/18s_autotrophic_euks.Rdata",
                                     env_data = "output/euks_auto_18sv9_surf_full_data.Rdata",
                                     mvrt_n_clust = 7, som_n_clust = 3, n_clust = 4,
                                     multicluster_output = "output/euks_auto_18sv9_surf_multi_cluster.Rdata",
                                     mvrt_multi_plot = "figures/euks_auto_18sv9_surf_mvrt_mutlipanel.pdf",
                                     som_multi_plot = "figures/euks_auto_18sv9_surf_som_mutlipanel.pdf",
                                     title_name = "Autotrophic Eukaryotes Surface Samples"){
  
  # load in datases for mvrt analysis
  load(raw_data)
  load(env_data)
  
  full_dat <- full_dat[complete.cases(full_dat[,c(33:34,37:39,47:48)]),]
  
  # removing smaples with too few reads
  
  scaled_inputs <- scaled_inputs[-which(!is.na(match(rownames(scaled_inputs),rownames(asv_table)[which(rowSums(asv_table) < 2000)]))),] 
  
  if(length(which(!is.na(match(full_dat$eco_name,rownames(asv_table)[which(rowSums(asv_table) < 2000)]))))>0){
    full_dat <- full_dat[-which(!is.na(match(full_dat$eco_name,rownames(asv_table)[which(rowSums(asv_table) < 2000)]))),] 
  }
  
  
  if(length(which(is.na(match(rownames(scaled_inputs), full_dat$eco_name)))) > 0){
    scaled_inputs <- scaled_inputs[-which(is.na(match(rownames(scaled_inputs), full_dat$eco_name))),]
  }
  # reorder env dataset
  ordered_dat <- full_dat[match(rownames(scaled_inputs),full_dat$eco_name),]
  
  # only pull explanatory variables
  env_dat <- ordered_dat[,c(33:34,37:39,47:48)]
  
  eco.som <- trainSOM(x.data = scaled_inputs, dimension = c(5, 5), nb.save = 10, maxit = 2000, 
                      scaling = "none")
  
  mvrt_list <- list()
  mvrt_robust_clusters <- vector()
  som_list <- list()
  som_robust_clusters <- vector()
  
  for(i in 2:15){
  par(mfrow = c(1,1))
  mvrt_result <- mvpart(scaled_inputs ~.,
                        env_dat,
                        size = i,
                        margin = 0.8,
                        xv = "none",
                        xval = 10,
                        xvmult = 1,
                        cp = 0,
                        plot.add = FALSE)
  
  eco.clust <- superClass(eco.som, k = i)
  
  mvrt_list[[i-1]] <- mvrt_result
  
  clusters <- eco.clust$cluster
  ids <- eco.clust$som$clustering
  som_ids <- clusters[ids]
  som_list[[i-1]] <- som_ids
  
  mvrt_robust_clusters[i-1] <- length(which(table(mvrt_result$where)/length(mvrt_result$where) > 0.05))
  som_robust_clusters[i-1] <- length(which(table(som_ids)/length(som_ids) > 0.05))
  
  }
  
  print(plot(2:15, mvrt_robust_clusters, type = "l"))
  print(plot(2:15, som_robust_clusters, type = "l"))
  
  mvrt_result <- mvrt_list[[mvrt_n_clust]]
  som_result <- som_list[[som_n_clust]]
  
  mvrt_robust_clust <- names(table(mvrt_result$where)[which(table(mvrt_result$where)/length(mvrt_result$where) > 0.05)])
  som_robust_clust <- names(table(som_result)[which(table(som_result)/length(som_result) > 0.05)])
  
  ordered_dat$som_id <- som_result
  
  if(length(mvrt_result$na.action) > 0){final_dat <- ordered_dat[-mvrt_result$na.action,]}else{final_dat <- ordered_dat}
  
  final_dat$mvrt_id <- mvrt_result$where
  
  if(length(which(is.na(match(final_dat$mvrt_id, as.numeric(mvrt_robust_clust))))) > 0){
    final_dat <-  final_dat[which(!is.na(match(final_dat$mvrt_id, as.numeric(mvrt_robust_clust)))),]
  }
  
  if(length(which(is.na(match(final_dat$som_id, as.numeric(som_robust_clust))))) > 0){
    final_dat <-  final_dat[which(!is.na(match(final_dat$som_id, as.numeric(som_robust_clust)))),]
  }
  
  final_dat <- final_dat[-which(as.numeric(substr(final_dat$Sta_ID,2,3)) < 76),]
  
  MVRTcolorCount = length(unique(final_dat$mvrt_id))
  SOMcolorCount = length(unique(final_dat$som_id))
  getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
  
  mvrt_plot <- ggplot(final_dat, aes(x = T_degC, y = NCDepth, color = as.factor(mvrt_id))) + geom_point(size = 2) +
    scale_color_manual(labels = paste0("Cluster ", 1:length(table(final_dat$mvrt_id))), values = getPalette(MVRTcolorCount)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          legend.position = "none") + 
    labs(color = "Multivariate\nRegression Tree\nClusters") + xlab("Temperature (°C)") +
    ylab("Nutricline Depth (m)") + scale_y_reverse() + stat_ellipse()
  
  som_plot <- ggplot(final_dat, aes(x = T_degC, y = NCDepth, color = as.factor(som_id))) + geom_point(size = 2) +
    scale_color_manual(labels = paste0("Cluster ", 1:length(table(final_dat$mvrt_id))), values = getPalette(SOMcolorCount)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          legend.position = "none") + 
    labs(color = "Self-Organizing\nMap Clusters") + xlab("Temperature (°C)") +
    ylab("Nutricline Depth (m)") + scale_y_reverse() + stat_ellipse()
  
  # ts diagram 
  
  mvrt_plot_ts <- ggplot(final_dat, aes(x = T_degC, y = Salnty, color = as.factor(mvrt_id))) + geom_point(size = 2) +
    scale_color_manual(labels = paste0("Cluster ", 1:length(table(final_dat$mvrt_id))), values = getPalette(MVRTcolorCount)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          legend.position = "none") + 
    labs(color = "Multivariate\nRegression Tree\nClusters") + xlab("Temperature (°C)") +
    ylab("Salinity") + stat_ellipse()
  
  som_plot_ts <- ggplot(final_dat, aes(x = T_degC, y = Salnty, color = as.factor(som_id))) + geom_point(size = 2) +
    scale_color_manual(labels = paste0("Cluster ", 1:length(table(final_dat$mvrt_id))), values = getPalette(SOMcolorCount)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          legend.position = "none") + 
    labs(color = "Self-Organizing\nMap Clusters") + xlab("Temperature (°C)") +
    ylab("Salinity") + stat_ellipse()
  
  # map
  
  map <- map_data("world")
  
  mvrt_map_list <- list()
  som_map_list <- list()
  
  for (i in 1:length(table(final_dat$mvrt_id))) {
    
    summary_maps <- final_dat %>% 
      group_by(Sta_ID) %>%
      summarise(mvrt_1 = sum(mvrt_id == as.numeric(names(table(final_dat$mvrt_id))[i]), na.rm = TRUE)/n(), 
                som_1 = sum(som_id == as.numeric(names(table(final_dat$som_id))[i]), na.rm = TRUE)/n(),
                n_samps = n(), lat = mean(Lat_Dec, na.rm= TRUE), long = mean(Lon_Dec, na.rm = TRUE))
    
    mvrt_map <-  ggplot() + 
      geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
      coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
      xlab("Longitude") + ylab("Latitude") + 
      geom_point(data = summary_maps, aes(x = long, y = lat, fill = mvrt_1), color = "black", size =6, stroke = 0.1, shape = 21) +
      scale_fill_gradient(low = "white", high = getPalette(MVRTcolorCount)[i], limits = c(0,1)) +
      ggtitle(paste0("MVRT Cluster ", i)) +
      theme(legend.title = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
            title = element_text(hjust = 0.5), axis.line = element_blank())
    
    som_map <-  ggplot() + 
      geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
      coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
      xlab("Longitude") + ylab("Latitude") + 
      geom_point(data = summary_maps, aes(x = long, y = lat, fill = som_1), color = "black", size =6, stroke = 0.1, shape = 21) +
      scale_fill_gradient(low = "white", high = getPalette(SOMcolorCount)[i], limits = c(0,1)) +
      ggtitle(paste0("SOM Cluster ", i)) +
      theme(legend.title = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
            title = element_text(hjust = 0.5), axis.line = element_blank())
    
    print(mvrt_map)
    print(som_map)
    
    mvrt_map_list[[i]] <- mvrt_map
    som_map_list[[i]] <- som_map
  }
  
  mvrt_ts_matrix <- as.data.frame(matrix(NA,length(unique(final_dat$Cruise)),(length(table(final_dat$mvrt_id)))+1))
  colnames(mvrt_ts_matrix)[1] <- "Date"
  
  som_ts_matrix <- as.data.frame(matrix(NA,length(unique(final_dat$Cruise)),(length(table(final_dat$som_id)))+1))
  colnames(som_ts_matrix)[1] <- "Date"
  
  for (i in 1:length(table(final_dat$mvrt_id))) {
    
    summary_ts <- final_dat %>% 
      group_by(Cruise) %>%
      summarise(mvrt_1 = sum(mvrt_id == as.numeric(names(table(final_dat$mvrt_id))[i]), na.rm = TRUE)/n(), 
                som_1 = sum(som_id == as.numeric(names(table(final_dat$som_id))[i]), na.rm = TRUE)/n(),
                n_samps = n(), Date = mean(Date, na.rm = TRUE))
    
    if(i == 1){
      mvrt_ts_matrix[,i] <- summary_ts$Date
      som_ts_matrix[,i] <- summary_ts$Date
    }
    mvrt_ts_matrix[,i+1] <- summary_ts$mvrt_1
    som_ts_matrix[,i+1] <- summary_ts$som_1
    
    colnames(mvrt_ts_matrix)[i+1] <- paste0("Cluster ",i)
    colnames(som_ts_matrix)[i+1] <- paste0("Cluster ",i)
  }
  
mvrt_ts_plot_df <- gather(mvrt_ts_matrix, "Cluster", "Freq", -Date)
som_ts_plot_df <- gather(som_ts_matrix, "Cluster", "Freq", -Date)

mvrt_time_series_plot <- ggplot(mvrt_ts_plot_df , aes(x = Date, y = Freq, fill = Cluster)) + 
                    geom_area(position = "stack") +
                    scale_fill_manual(labels = paste0("Cluster ",
                    1:length(table(final_dat$mvrt_id))),
                    values = getPalette(MVRTcolorCount)) + 
                    theme(panel.background = element_blank(),
                    panel.border = element_rect(color= "black", fill = NA),
                    legend.position = "none") +
                    ylab("Frequency") +
                    labs(color = "Multivariate\nRegression Tree\nClusters")

som_time_series_plot <- ggplot(som_ts_plot_df, aes(x = Date, y = Freq, fill = Cluster)) + 
  geom_area(position = "stack", color = "black") +
  scale_fill_manual(labels = paste0("Cluster ",
                                    1:length(table(final_dat$mvrt_id))),
                    values = getPalette(SOMcolorCount)) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color= "black", fill = NA),
        legend.position = "none") +
  ylab("Frequency") +
  labs(color = "Self-Organizing\nMap Clusters")

months <- substr(final_dat$Cruise,5,6)

months[as.numeric(months) < 3] <- "Winter"
months[as.numeric(months) < 5] <- "Spring"
months[as.numeric(months) < 9] <- "Summer"
months[as.numeric(months) < 15] <- "Fall"

final_dat$Season <- months

mvrt_season_matrix <- as.data.frame(matrix(NA,length(unique(final_dat$Season)),(length(table(final_dat$mvrt_id)))+1))
colnames(mvrt_season_matrix)[1] <- "Season"

som_season_matrix <- as.data.frame(matrix(NA,length(unique(final_dat$Season)),(length(table(final_dat$som_id)))+1))
colnames(som_season_matrix)[1] <- "Season"

for (i in 1:length(table(final_dat$mvrt_id))) {
  
  summary_season <- final_dat %>% 
    group_by(Season) %>%
    summarise(mvrt_1 = sum(mvrt_id == as.numeric(names(table(final_dat$mvrt_id))[i]), na.rm = TRUE)/n(),
              som_1 = sum(som_id == as.numeric(names(table(final_dat$som_id))[i]), na.rm = TRUE)/n(),
              n_samps = n())
  
  if(i == 1){mvrt_season_matrix[,i] <- summary_season$Season}
  if(i == 1){som_season_matrix[,i] <- summary_season$Season}
  
  mvrt_season_matrix[,i+1] <- summary_season$mvrt_1
  som_season_matrix[,i+1] <- summary_season$som_1
  
  colnames(mvrt_season_matrix)[i+1] <- paste0("Cluster ",i)
  colnames(som_season_matrix)[i+1] <- paste0("Cluster ",i)
  
}

mvrt_season_plot_df <- gather(mvrt_season_matrix, "Cluster", "Freq", -Season)
som_season_plot_df <- gather(som_season_matrix, "Cluster", "Freq", -Season)

mvrt_season_plot <- ggplot(mvrt_season_plot_df, aes(x = Season, y = Freq, fill = Cluster)) + 
               geom_bar(stat = "identity") +
               scale_fill_manual(labels = paste0("Cluster ",
               1:length(table(final_dat$mvrt_id))),
               values = getPalette(MVRTcolorCount)) + 
               theme(panel.background = element_blank(),
               panel.border = element_rect(color= "black", fill = NA),
               legend.position = "bottom") +
               ylab("Frequency") +
               labs(fill = "Multivariate Regression Tree Clusters")

som_season_plot <- ggplot(som_season_plot_df, aes(x = Season, y = Freq, fill = Cluster)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(labels = paste0("Cluster ",
                                    1:length(table(final_dat$mvrt_id))),
                    values = getPalette(SOMcolorCount)) + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(color= "black", fill = NA),
        legend.position = "bottom") +
  ylab("Frequency") +
  labs(fill = "Self-Organizing Clusters")

# Region plot

# region_matrix <- as.data.frame(matrix(NA,length(unique(final_dat$Code_CCE)),(length(table(final_dat$mvrt_id)))+1))
# colnames(region_matrix)[1] <- "Region"
# 
# for (i in 1:length(table(final_dat$mvrt_id))) {
#   
#   summary_region <- final_dat %>% 
#     group_by(Code_CCE) %>%
#     summarise(mvrt_1 = sum(mvrt_id == as.numeric(names(table(final_dat$mvrt_id))[i]), na.rm = TRUE)/n(), 
#               n_samps = n())
#   
#   if(i == 1){region_matrix[,i] <- summary_region$Code_CCE}
#   region_matrix[,i+1] <- summary_region$mvrt_1
#   colnames(region_matrix)[i+1] <- paste0("Cluster ",i)
#   
# }
# 
# region_plot_df <- gather(region_matrix, "Cluster", "Freq", -Region)
# 
# region_plot_df <- region_plot_df[-which(region_plot_df$Region == ""),]
# 
# region_plot <- ggplot(region_plot_df, aes(x = Region, y = Freq, fill = Cluster)) + 
#   geom_bar(stat = "identity") +
#   scale_fill_manual(labels = paste0("Cluster ",
#                                     1:length(table(final_dat$mvrt_id))),
#                     values = getPalette(colourCount)) + 
#   theme(panel.background = element_blank(),
#         panel.border = element_rect(color= "black", fill = NA),
#         legend.position = "bottom") +
#   ylab("Frequency") +
#   labs(fill = "Multivariate Regression Tree Clusters")
# 
# print(region_plot)

mvrt_cluster_leg <- get_legend(mvrt_season_plot)
som_cluster_leg <- get_legend(som_season_plot)

mvrt_season_plot <- mvrt_season_plot + theme(legend.position = "none")
som_season_plot <- som_season_plot + theme(legend.position = "none")

mvrt_title <- ggdraw() + draw_label(paste0("MVRT Analysis ",title_name), fontface='bold')
som_title <- ggdraw() + draw_label(paste0("SOM Analysis ",title_name), fontface='bold')

# p1_tree %<a-% {
#   plot(mvrt_result, margin = 0.1)
#   text(mvrt_result)
# }

# get maps
mvrt_a <- plot_grid(mvrt_map_list[[1]],
                    mvrt_map_list[[2]],
                    mvrt_map_list[[3]],
                    mvrt_map_list[[4]], 
                    # mvrt_map_list[[5]],
                    ncol = n_clust)

som_a <- plot_grid(som_map_list[[1]],
                   som_map_list[[2]],
                   som_map_list[[3]],
                   som_map_list[[4]], 
                   # som_map_list[[5]],
                    ncol = n_clust)

# b <- ~p1_tree

mvrt_large_plot <- plot_grid(mvrt_title,
                        plot_grid(mvrt_plot, mvrt_plot_ts, ncol = 2),
                        mvrt_a,
                        plot_grid(mvrt_time_series_plot, mvrt_season_plot, rel_widths = c(1,0.6)),
                        mvrt_cluster_leg, nrow = 5, rel_heights = c(0.1,1,1,1,0.1))

som_large_plot <- plot_grid(som_title,
                             plot_grid(som_plot, som_plot_ts, ncol = 2),
                             som_a,
                             plot_grid(som_time_series_plot, som_season_plot, rel_widths = c(1,0.6)),
                             som_cluster_leg, nrow = 5, rel_heights = c(0.1,1,1,1,0.1))

  
pdf(mvrt_multi_plot, width = 15, height = 15)
print(mvrt_large_plot)
dev.off()

pdf(som_multi_plot, width = 15, height = 15)
print(som_large_plot)
dev.off()


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
  
  pdf(ind_plot, width = 15, height = plt_height)
  print(plot_grid(title, plot_grid(p1,p2,ncol = 2), ncol = 1, rel_heights = c(0.1,1)))
  dev.off()
  
}



# Multicluster Analysis

multicluster_analysis(raw_data = "data/18s_autotrophic_euks.Rdata",
                      env_data = "output/euks_auto_18sv9_surf_full_data.Rdata",
                      mvrt_n_clust = 7, som_n_clust = 3, n_clust = 4,
                      multicluster_output = "output/euks_auto_18sv9_surf_multi_cluster.Rdata",
                      mvrt_multi_plot = "figures/euks_auto_18sv9_surf_mvrt_mutlipanel.pdf",
                      som_multi_plot = "figures/euks_auto_18sv9_surf_som_mutlipanel.pdf",
                      title_name = "Autotrophic Eukaryotes Surface Samples")

multicluster_analysis(raw_data = "data/18s_autotrophic_euks.Rdata",
                      env_data = "output/euks_auto_18sv9_deep_full_data.Rdata",
                      mvrt_n_clust = 7, som_n_clust = 3, n_clust = 4,
                      multicluster_output = "output/euks_auto_18sv9_deep_multi_cluster.Rdata",
                      mvrt_multi_plot = "figures/euks_auto_18sv9_deep_mvrt_mutlipanel.pdf",
                      som_multi_plot = "figures/euks_auto_18sv9_deep_som_mutlipanel.pdf",
                      title_name = "Autotrophic Eukaryotes Deep Samples")

multicluster_analysis(raw_data = "data/18s_heterotrophic_euks.Rdata",
                      env_data = "output/euks_hetero_18sv9_surf_full_data.Rdata",
                      mvrt_n_clust = 3, som_n_clust = 3, n_clust = 4,
                      multicluster_output = "output/euks_hetero_18sv9_surf_multi_cluster.Rdata",
                      mvrt_multi_plot = "figures/euks_hetero_18sv9_surf_mvrt_mutlipanel.pdf",
                      som_multi_plot = "figures/euks_hetero_18sv9_surf_som_mutlipanel.pdf",
                      title_name = "Heterotrophic Eukaryotes Surface Samples")

multicluster_analysis(raw_data = "data/18s_heterotrophic_euks.Rdata",
                      env_data = "output/euks_hetero_18sv9_deep_full_data.Rdata",
                      mvrt_n_clust = 7, som_n_clust = 3, n_clust = 4,
                      multicluster_output = "output/euks_hetero_18sv9_deep_multi_cluster.Rdata",
                      mvrt_multi_plot = "figures/euks_hetero_18sv9_deep_mvrt_mutlipanel.pdf",
                      som_multi_plot = "figures/euks_hetero_18sv9_deep_som_mutlipanel.pdf",
                      title_name = "Heterotrophic Eukaryotes Deep Samples")

multicluster_analysis(raw_data = "data/16s_cyanos.Rdata",
                      env_data = "output/cyano_16s_surf_full_data.Rdata",
                      mvrt_n_clust = 3, som_n_clust = 3, n_clust = 4,
                      multicluster_output = "output/cyano_16s_surf_multi_cluster.Rdata",
                      mvrt_multi_plot = "figures/cyano_16s_surf_mvrt_mutlipanel.pdf",
                      som_multi_plot = "figures/cyano_16s_surf_som_mutlipanel.pdf",
                      title_name = "Cyanobacteria Surface Samples")

multicluster_analysis(raw_data = "data/16s_cyanos.Rdata",
                      env_data = "output/cyano_16s_deep_full_data.Rdata",
                      mvrt_n_clust = 3, som_n_clust = 3, n_clust = 4,
                      multicluster_output = "output/cyano_16s_deep_multi_cluster.Rdata",
                      mvrt_multi_plot = "figures/cyano_16s_deep_mvrt_mutlipanel.pdf",
                      som_multi_plot = "figures/cyano_16s_deep_som_mutlipanel.pdf",
                      title_name = "Cyanobacteria Deep Samples")

multicluster_analysis(raw_data = "data/16s_bacteria_m_euks.Rdata",
                      env_data = "output/bacteria_m_euks_16s_surf_full_data.Rdata",
                      mvrt_n_clust = 3, som_n_clust = 3, n_clust = 4,
                      multicluster_output = "output/bacteria_m_euks_16s_surf_multi_cluster.Rdata",
                      mvrt_multi_plot = "figures/bacteria_m_euks_16s_surf_mvrt_mutlipanel.pdf",
                      som_multi_plot = "figures/bacteria_m_euks_16s_surf_som_mutlipanel.pdf",
                      title_name = "Bacteria/Archaea Surface Samples")

multicluster_analysis(raw_data = "data/16s_bacteria_m_euks.Rdata",
                      env_data = "output/bacteria_m_euks_16s_deep_full_data.Rdata",
                      mvrt_n_clust = 3, som_n_clust = 3, n_clust = 4,
                      multicluster_output = "output/bacteria_m_euks_16s_deep_multi_cluster.Rdata",
                      mvrt_multi_plot = "figures/bacteria_m_euks_16s_deep_mvrt_mutlipanel.pdf",
                      som_multi_plot = "figures/bacteria_m_euks_16s_deep_som_mutlipanel.pdf",
                      title_name = "Bacteria/Archaea Deep Samples")

multicluster_analysis(raw_data = "data/16s_plastids.Rdata",
                      env_data = "output/plastid_16s_surf_full_data.Rdata",
                      mvrt_n_clust = 3, som_n_clust = 3, n_clust = 4,
                      multicluster_output = "output/plastid_16s_surf_multi_cluster.Rdata",
                      mvrt_multi_plot = "figures/plastids_16s_surf_mvrt_mutlipanel.pdf",
                      som_multi_plot = "figures/plastid_16s_surf_som_mutlipanel.pdf",
                      title_name = "Eukaryotic Plastid Surface Samples")

multicluster_analysis(raw_data = "data/16s_plastids.Rdata",
                      env_data = "output/plastid_16s_deep_full_data.Rdata",
                      mvrt_n_clust = 3, som_n_clust = 3, n_clust = 4,
                      multicluster_output = "output/plastid_16s_deep_multi_cluster.Rdata",
                      mvrt_multi_plot = "figures/plastid_16s_deep_mvrt_mutlipanel.pdf",
                      som_multi_plot = "figures/plastid_16s_deep_som_mutlipanel.pdf",
                      title_name = "Eukaryotic Plastid Deep Samples")

