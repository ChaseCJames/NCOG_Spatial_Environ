library(tidyverse)
library(spatialEco)
library(sp)
library(patchwork)
library(pscl)


meta <- read.csv("data/NCOG_sample_log_DNA_stvx_meta_2014-2020_prim_prod.csv")

in_group_list = c("archaea_16s", "bacteria_m_euks_16s","cyano_16s",
                  "euks_hetero_18sv9","euks_auto_18sv9")

in_group_names = c("Archaea","Bacteria","Cyanobacteria",
                   "Heterotrophic Eukaryotic Protists","Photosynthetic Eukaryotic Protists")

meta$IntC14_day <- meta$IntC14*2

meta$ef_ratio <- 0.04756*(0.78 - ((0.43*meta$T_degC)/30))*(meta$IntC14_day^0.307)

meta$ef_ratio <- ((0.5857-(0.0165*meta$T_degC))*meta$IntC14_day)/(51.7 + meta$IntC14_day)

meta_sta <- meta %>% group_by(Sta_ID) %>%
  summarise(mean_ef = mean(ef_ratio, na.rm= TRUE)) %>%
  filter(!is.nan(mean_ef))

plot_list <- list()

for (i in 1:length(in_group_list)) {
  
  load(paste0("output/",in_group_list[i],"_map_S.Rdata"))
  
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  wt_1 <- wt.centroid(x = centroid_df, p = 2)
  wt_2 <- wt.centroid(x = centroid_df, p = 3)
  
  clust1 <- which.max(c(wt_1@coords[1], wt_2@coords[1]))
  clust2 <- which.min(c(wt_1@coords[1], wt_2@coords[1]))
  
  som_maps$mean_ef <- meta_sta$mean_ef[match(som_maps$Sta_ID,meta_sta$Sta_ID)]
  
  som_piv <- som_maps[,c(1:3,5:6,30)] %>% pivot_longer(-c(Sta_ID,lat, long, mean_ef), names_to = "som", values_to = "freq")
  
  if(clust1 == 1){som_piv$som[which(som_piv$som == "som_1")] <- "Nearshore"}
  if(clust1 == 2){som_piv$som[which(som_piv$som == "som_2")] <- "Nearshore"}
  
  if(clust2 == 1){som_piv$som[which(som_piv$som == "som_1")] <- "Offshore"}
  if(clust2 == 2){som_piv$som[which(som_piv$som == "som_2")] <- "Offshore"}
  
  som_fit <- glm(som_1 ~ mean_ef, family = binomial, data = som_maps) 
  
  p_r2 <- round(pR2(som_fit)[6],digits = 3)
  p_val <- round(summary(som_fit)$coefficients[2,4],digits = 3)
  
  if(p_val == 0){p_val <- "< 0.001"}else{p_val <- paste0("= ",p_val)}
  

  plot_list[[i]] <- ggplot(som_piv, aes(x = mean_ef, y = freq, color = som)) +
    geom_point() +
    scale_y_continuous(limits = c(0,1.05), breaks = c(0,0.25,0.5,0.75,1)) +
    scale_x_continuous(limits = c(0.218,0.379)) +
    stat_smooth(method="glm", method.args = list(family = "binomial")) +
    scale_color_manual(values = c("blue", "red")) +
    labs(x = "Mean ef-ratio", y = "Frequency", color = "") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.key = element_blank(),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12)) + 
    ggtitle(in_group_names[i]) +
    annotate(geom = "text", x = 0.24, y = 0.55, label = paste0("R[pseudo]^2 == ", p_r2), parse = TRUE, size = 3) +
    annotate(geom = "text", x = 0.24, y = 0.45, label = paste0("p-value ", p_val), parse = FALSE, size = 3)
  
}


  out_plot <- (plot_list[[1]] + theme(axis.text.x = element_blank(),
                          axis.title.x = element_blank()) +
      annotate(geom = "text", x = 0.22, y = 1.042, label = "a", size = 6, fontface = "bold")) +
  (plot_list[[2]] + theme(axis.text = element_blank(),
                          axis.title = element_blank()) +
     annotate(geom = "text", x = 0.22, y = 1.042, label = "b", size = 6, fontface = "bold")) +
  guide_area() +
  (plot_list[[3]] +
     annotate(geom = "text", x = 0.22, y = 1.042, label = "c", size = 6, fontface = "bold")) +
  (plot_list[[5]] + theme(axis.text.y = element_blank(),
                          axis.title.y = element_blank()) +
     annotate(geom = "text", x = 0.22, y = 1.042, label = "d", size = 6, fontface = "bold")) +
  (plot_list[[4]] + theme(axis.text.y = element_blank(),
                          axis.title.y = element_blank()) +
     annotate(geom = "text", x = 0.22, y = 1.042, label = "e", size = 6, fontface = "bold")) +
  plot_layout(guides = "collect")

  pdf("figures_S/suppl_fig_4_S.pdf", width = 12, height = 7)
  plot(out_plot)
  dev.off()



