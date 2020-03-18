
library(tidyverse)
library(rworldmap)
library(ggplot2)
library(ggmap)
library(cowplot)
library(scales)
library(gridExtra)
library(lubridate)
library(vegan)
library(spatialEco)
library(reshape2)
library(rsq)


aic_table_func_bt <- function(som_maps = cyano_plots, i = 1){  
  
  som_glm <- som_maps[,c(61:63,33:34,37:39,41,48,59)]
  
  colnames(som_glm)[i] <- "response"
  
  model_AIC <- matrix(NA,8,2)
  model_p_val <- matrix(NA,8,2)
  
  
  # temperature
  
  glm_mean_temp <- glm(response ~  T_degC, data = som_glm)
  mt_sum <- rsq(glm_mean_temp)
  model_AIC[1,2] <- mt_sum
  model_AIC[1,1] <- "Temperature"
  
  # salinity
  
  glm_mean_sal <- glm(response ~  Salnty, data = som_glm)
  mt_sum <- rsq(glm_mean_sal)
  model_AIC[2,2] <- mt_sum
  model_AIC[2,1] <- "Salinity"
  
  # NO3
  
  glm_mean_no3 <- glm(response ~  NO3ug, data = som_glm)
  mn_sum <- rsq(glm_mean_no3)
  model_AIC[3,2] <- mn_sum
  model_AIC[3,1] <- "NO3"
  
  # PO4
  
  glm_mean_po4 <- glm(response ~  PO4ug, data = som_glm)
  mp_sum <- rsq(glm_mean_po4)
  model_AIC[4,2] <- mp_sum
  model_AIC[4,1] <- "PO4"
  
  # SiO3
  
  glm_mean_sio3 <- glm(response ~  SiO3ug, data = som_glm)
  ms_sum <- rsq(glm_mean_sio3)
  model_AIC[5,2] <- ms_sum
  model_AIC[5,1] <- "SiO3"
  
  # NC Depth
  
  glm_mean_nc <- glm(response ~  NCDepth, data = som_glm)
  ms_sum <- rsq(glm_mean_nc)
  model_AIC[6,2] <- ms_sum
  model_AIC[6,1] <- "Nutricline Depth (m)"
  
  # Distance to Coast
  
  glm_mean_dc <- glm(response ~  dist_to_coast, data = som_glm)
  ms_sum <- rsq(glm_mean_dc)
  model_AIC[7,2] <- ms_sum
  model_AIC[7,1] <- "Distance to Coast"
  
  # Chlorophyll
  
  glm_mean_chl <- glm(response ~  ChlorA, data = som_glm)
  ms_sum <- rsq(glm_mean_chl)
  model_AIC[8,2] <- ms_sum
  model_AIC[8,1] <- "Chl-a"
  
  return(model_AIC)
  
}

full_aic_table_figure_diversity <- function(
  in_cyano = "output/cyano_16s_full_data.Rdata",
  in_bacteria = "output/bacteria_m_euks_16s_full_data.Rdata",
  in_euks = "output/euks_hetero_18sv9_full_data.Rdata",
  in_phyto = "output/euks_auto_18sv9_full_data.Rdata", 
  figure_name = paste0("figures/full_aic_table_ns_rich_",Sys.Date(),".pdf"),
  figure_name_2 = paste0("figures/full_aic_plot_ns_rich_",Sys.Date(),".pdf"),
  title_name = "R-Squared Nearshore Richness", # col 1 = shan, col 2 = even col 3 = rich
  type = "Nearshore", color_fill = "purple", col = 1){
  # load an rename data
  
  # load(in_plastid)
  # som_maps <- som_maps[which(som_maps$n_samps > (minimum_tp-1)),]
  # plastid_plots <- som_maps
  
  # som 2 nearshore
  load(in_cyano)
  cyano_plots <- full_dat
  cyano_plots$som_id[which(cyano_plots$som_id == 2)] = "Nearshore"
  cyano_plots$som_id[which(cyano_plots$som_id == 1)] = "Offshore"
  
  load(in_bacteria)
  bacteria_plots <- full_dat
  bacteria_plots$som_id[which(bacteria_plots$som_id == 2)] = "Nearshore"
  bacteria_plots$som_id[which(bacteria_plots$som_id == 1)] = "Offshore"
  
  load(in_euks)
  eukaryota_plots <- full_dat
  eukaryota_plots$som_id[which(eukaryota_plots$som_id == 2)] = "Nearshore"
  eukaryota_plots$som_id[which(eukaryota_plots$som_id == 1)] = "Offshore"
  
  load(in_phyto) 
  phyto_plots <- full_dat
  phyto_plots$som_id[which(phyto_plots$som_id == 1)] = "Nearshore"
  phyto_plots$som_id[which(phyto_plots$som_id == 2)] = "Offshore"
    
  cyano_plots <- cyano_plots %>% filter(som_id == type)
  bacteria_plots <- bacteria_plots %>% filter(som_id == type)
  eukaryota_plots <- eukaryota_plots %>% filter(som_id == type)
  phyto_plots <- phyto_plots %>% filter(som_id == type)

  
  cyano_AIC <- aic_table_func_bt(som_maps = cyano_plots, i = col)
  cyano_AIC <- as.data.frame(cyano_AIC, stringsAsFactors = FALSE)
  colnames(cyano_AIC) <- c("Variables","AIC")
  cyano_AIC$AIC <- as.numeric(cyano_AIC$AIC)
  cyano_AIC$AIC <- round(cyano_AIC$AIC, digits = 2)
  
  bacteria_AIC <- aic_table_func_bt(som_maps = bacteria_plots, i = col)
  bacteria_AIC <- as.data.frame(bacteria_AIC, stringsAsFactors = FALSE)
  colnames(bacteria_AIC) <- c("Variables","AIC")
  bacteria_AIC$AIC <- as.numeric(bacteria_AIC$AIC)
  bacteria_AIC$AIC <- round(bacteria_AIC$AIC, digits = 2)
  
  eukaryota_AIC <- aic_table_func_bt(som_maps = eukaryota_plots, i = col)
  eukaryota_AIC <- as.data.frame(eukaryota_AIC, stringsAsFactors = FALSE)
  colnames(eukaryota_AIC) <- c("Variables","AIC")
  eukaryota_AIC$AIC <- as.numeric(eukaryota_AIC$AIC)
  eukaryota_AIC$AIC <- round(eukaryota_AIC$AIC, digits = 2)
  
  phyto_AIC  <- aic_table_func_bt(som_maps = phyto_plots, i = col)
  phyto_AIC <- as.data.frame(phyto_AIC, stringsAsFactors = FALSE)
  colnames(phyto_AIC) <- c("Variables","AIC")
  phyto_AIC$AIC <- as.numeric(phyto_AIC$AIC)
  phyto_AIC$AIC <- round(phyto_AIC$AIC, digits = 2)
  
  AIC_full <- full_join(cyano_AIC, bacteria_AIC, by = "Variables")
  AIC_full <- full_join(AIC_full, eukaryota_AIC, by = "Variables")
  AIC_full <- full_join(AIC_full, phyto_AIC, by = "Variables")
  
  colnames(AIC_full) <- c("Variables", "Cyanobacteria\n AIC",
                          "Bacteria/Archaea\n AIC", "Heterotrophic\nEukaryotes\n AIC", "Eukaryotic\nPhytoplankton\n AIC")
  
  colfunc <- colorRampPalette(c("white", color_fill))
  
  # cyano_col <- round(seq(max(cyano_AIC$AIC, na.rm = TRUE),min(cyano_AIC$AIC, na.rm = TRUE), by = -0.01), digits = 2)
  cyano_col <- seq(0,1,by = 0.01)
  scale <- colfunc(length(cyano_col))
  cyano_fill <- scale[match(cyano_AIC$AIC, cyano_col)]
  # cyano_fill[5] <- scale[3903]
  
  # bacteria_col <- round(seq(max(bacteria_AIC$AIC, na.rm = TRUE),min(bacteria_AIC$AIC, na.rm = TRUE), by = -0.01), digits = 2)
  bacteria_col <- seq(0,1,by = 0.01)
  scale <- colfunc(length(bacteria_col))
  bacteria_fill <- scale[match(bacteria_AIC$AIC, bacteria_col)]
  
  # eukaryota_col <- round(seq(max(eukaryota_AIC$AIC, na.rm = TRUE),min(eukaryota_AIC$AIC, na.rm = TRUE), by = -0.01), digits = 2)
  eukaryota_col <- seq(0,1,by = 0.01)
  scale <- colfunc(length(eukaryota_col))
  eukaryota_fill <- scale[match(round(eukaryota_AIC$AIC, digits = 2), eukaryota_col)]
  
  # phyto_col <- round(seq(max(phyto_AIC$AIC, na.rm = TRUE),min(phyto_AIC$AIC, na.rm = TRUE), by = -0.01), digits = 2)
  phyto_col <- seq(0,1,by = 0.01)
  scale <- colfunc(length(phyto_col))
  phyto_fill <- scale[match(round(phyto_AIC$AIC, digits = 2), phyto_col)]
  
  
  t0 <- tableGrob(AIC_full["Variables"], 
                  theme=ttheme_default(
                    core=list(bg_params = list(fill="grey90", col = "black"),
                              fg_params = list(fontface="bold")),
                    colhead = list(bg_params=list(fill="white", col="black"))), 
                  rows = NULL)
  # t1 <- tableGrob(AIC_full["Eukaryotic\nPlastid\n AIC"], 
  #                 theme=ttheme_default(
  #                   core=list(bg_params = list(fill=plastid_fill, col = "black")),
  #                   colhead = list(bg_params=list(fill="white", col="black"))), 
  #                 rows = NULL)
  t1 <- tableGrob(AIC_full["Cyanobacteria\n AIC"],
                  theme=ttheme_default(
                    core=list(bg_params = list(fill=cyano_fill, col = "black")),
                    colhead = list(bg_params=list(fill="white", col="black"))),
                  rows = NULL)
  t2 <- tableGrob(AIC_full["Bacteria/Archaea\n AIC"],
                  theme=ttheme_default(
                    core=list(bg_params = list(fill=bacteria_fill, col = "black")),
                    colhead = list(bg_params=list(fill="white", col="black"))),
                  rows = NULL)
  t3 <- tableGrob(AIC_full["Heterotrophic\nEukaryotes\n AIC"],
                  theme=ttheme_default(
                    core=list(bg_params = list(fill=eukaryota_fill, col = "black")),
                    colhead = list(bg_params=list(fill="white", col="black"))),
                  rows = NULL)
  t4 <- tableGrob(AIC_full["Eukaryotic\nPhytoplankton\n AIC"],
                  theme=ttheme_default(
                    core=list(bg_params = list(fill=phyto_fill, col = "black")),
                    colhead = list(bg_params=list(fill="white", col="black"))),
                  rows = NULL)
  
  # join tables
  tab <- gtable_combine(t0,t1,t2,t3,t4)
  
  pdf(file = figure_name, width = 10, height = 6)
  grid.arrange(tab)
  dev.off()
  
  AIC_scaled <- AIC_full
  
  for (i in 2:ncol(AIC_scaled)) {
    zero_one_scale <- AIC_scaled[,i]
    AIC_scaled[,i] <- 25*zero_one_scale
    
  }
  
  plot_df <- melt(AIC_scaled)
  
  
  
  colnames(plot_df) <- c("Variables", "Group", "AIC")
  
  # # removing plastids from plot
  # plot_df <- plot_df[-which(plot_df$Group == "Eukaryotic\nPlastid\n AIC"),]
  
  plot_df$Variables <- as.factor(plot_df$Variables)
  plot_df$Variables <- factor(plot_df$Variables, levels = c("Chl-a",
                                                            "Distance to Coast", "Nutricline Depth (m)",
                                                            "SiO3", "PO4",
                                                            "NO3","Salinity",
                                                            "Temperature" ))
  
  pdf(figure_name_2, width = 8, height = 8)
  print(ggplot(data = plot_df, aes(x = Group, y = Variables, size = AIC)) + 
          geom_point(fill = color_fill, color = "black", alpha = 0.6, shape = 21) +
          labs(size = "Variable\n Importance") + ylab("Variable") +
          theme(panel.background = element_blank(),
                panel.border = element_rect(color = "black", fill = NA),
                legend.position = "none",
                panel.grid.major.y = element_line(color = "grey", linetype = 2),
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 0)) + ggtitle(title_name) +
          scale_size_continuous(range = c(1,15)) + xlab("") + ylab(""))
  dev.off()
  
  
}

corr_table_func_bt <- function(som_maps = cyano_plots, i = 1){  
  
  som_glm <- som_maps[,c(61:63,33:34,37:39,41,48,59,64:75)]
  
  colnames(som_glm)[i] <- "response"
  
  cor_mat <- cor(som_glm, use = "pairwise.complete.obs")
  
  tab <- cor_mat[-c(1:3),i]
  
  model_AIC <- matrix(NA,20,2)
  
  model_AIC[,1] <- names(tab)
  model_AIC[,2] <- tab
  
  return(model_AIC)
  
}

full_aic_table_figure_diff <- function(
  in_cyano = "output/cyano_16s_full_data.Rdata",
  in_bacteria = "output/bacteria_m_euks_16s_full_data.Rdata",
  in_euks = "output/euks_hetero_18sv9_full_data.Rdata",
  in_phyto = "output/euks_auto_18sv9_full_data.Rdata", 
  figure_name = paste0("figures/full_table_diff_shan_",Sys.Date(),".pdf"),
  figure_name_2 = paste0("figures/full_table_plot_diff_shan_",Sys.Date(),".pdf"),
  title_name = "Correlation Nearshore Richness", # col 1 = shan, col 2 = even col 3 = rich
  col = 1){
  # load an rename data
  
  # load(in_plastid)
  # som_maps <- som_maps[which(som_maps$n_samps > (minimum_tp-1)),]
  # plastid_plots <- som_maps
  
  # som 2 nearshore
  load(in_cyano)
  cyano_plots <- full_dat
  cyano_plots$som_id[which(cyano_plots$som_id == 2)] = "Nearshore"
  cyano_plots$som_id[which(cyano_plots$som_id == 1)] = "Offshore"
  
  load(in_bacteria)
  bacteria_plots <- full_dat
  bacteria_plots$som_id[which(bacteria_plots$som_id == 2)] = "Nearshore"
  bacteria_plots$som_id[which(bacteria_plots$som_id == 1)] = "Offshore"
  
  load(in_euks)
  eukaryota_plots <- full_dat
  eukaryota_plots$som_id[which(eukaryota_plots$som_id == 2)] = "Nearshore"
  eukaryota_plots$som_id[which(eukaryota_plots$som_id == 1)] = "Offshore"
  
  load(in_phyto) 
  phyto_plots <- full_dat
  phyto_plots$som_id[which(phyto_plots$som_id == 1)] = "Nearshore"
  phyto_plots$som_id[which(phyto_plots$som_id == 2)] = "Offshore"
  
  # richness
  
  cyano_plots$cyano_rich <- cyano_plots$richness
  cyano_plots$bact_rich <- bacteria_plots$richness
  cyano_plots$eukary_rich <- eukaryota_plots$richness
  cyano_plots$phyto_rich <- phyto_plots$richness
  
  bacteria_plots$cyano_rich <- cyano_plots$richness
  bacteria_plots$bact_rich <- bacteria_plots$richness
  bacteria_plots$eukary_rich <- eukaryota_plots$richness
  bacteria_plots$phyto_rich <- phyto_plots$richness
  
  eukaryota_plots$cyano_rich <- cyano_plots$richness
  eukaryota_plots$bact_rich <- bacteria_plots$richness
  eukaryota_plots$eukary_rich <- eukaryota_plots$richness
  eukaryota_plots$phyto_rich <- phyto_plots$richness
  
  phyto_plots$cyano_rich <- cyano_plots$richness
  phyto_plots$bact_rich <- bacteria_plots$richness
  phyto_plots$eukary_rich <- eukaryota_plots$richness
  phyto_plots$phyto_rich <- phyto_plots$richness
  
  # evenness

  cyano_plots$cyano_even <- cyano_plots$evenness
  cyano_plots$bact_even <- bacteria_plots$evenness
  cyano_plots$eukary_even <- eukaryota_plots$evenness
  cyano_plots$phyto_even <- phyto_plots$evenness
  
  bacteria_plots$cyano_even <- cyano_plots$evenness
  bacteria_plots$bact_even <- bacteria_plots$evenness
  bacteria_plots$eukary_even <- eukaryota_plots$evenness
  bacteria_plots$phyto_even <- phyto_plots$evenness
  
  eukaryota_plots$cyano_even <- cyano_plots$evenness
  eukaryota_plots$bact_even <- bacteria_plots$evenness
  eukaryota_plots$eukary_even <- eukaryota_plots$evenness
  eukaryota_plots$phyto_even <- phyto_plots$evenness
  
  phyto_plots$cyano_even <- cyano_plots$evenness
  phyto_plots$bact_even <- bacteria_plots$evenness
  phyto_plots$eukary_even <- eukaryota_plots$evenness
  phyto_plots$phyto_even <- phyto_plots$evenness
  
  # shannon
  
  cyano_plots$cyano_shannon <- cyano_plots$shannon
  cyano_plots$bact_shannon <- bacteria_plots$shannon
  cyano_plots$eukary_shannon <- eukaryota_plots$shannon
  cyano_plots$phyto_shannon <- phyto_plots$shannon
  
  bacteria_plots$cyano_shannon <- cyano_plots$shannon
  bacteria_plots$bact_shannon <- bacteria_plots$shannon
  bacteria_plots$eukary_shannon <- eukaryota_plots$shannon
  bacteria_plots$phyto_shannon <- phyto_plots$shannon
  
  eukaryota_plots$cyano_shannon <- cyano_plots$shannon
  eukaryota_plots$bact_shannon <- bacteria_plots$shannon
  eukaryota_plots$eukary_shannon <- eukaryota_plots$shannon
  eukaryota_plots$phyto_shannon <- phyto_plots$shannon
  
  phyto_plots$cyano_shannon <- cyano_plots$shannon
  phyto_plots$bact_shannon <- bacteria_plots$shannon
  phyto_plots$eukary_shannon <- eukaryota_plots$shannon
  phyto_plots$phyto_shannon <- phyto_plots$shannon
  
  cyano_ns <- cyano_plots %>% filter(som_id == "Nearshore")
  bacteria_ns <- bacteria_plots %>% filter(som_id == "Nearshore")
  eukaryota_ns <- eukaryota_plots %>% filter(som_id == "Nearshore")
  phyto_ns <- phyto_plots %>% filter(som_id == "Nearshore")
  
  
  cyano_AIC <- corr_table_func_bt(som_maps = cyano_ns, i = col)
  cyano_AIC <- as.data.frame(cyano_AIC, stringsAsFactors = FALSE)
  colnames(cyano_AIC) <- c("Variables","AIC")
  cyano_AIC$AIC <- as.numeric(cyano_AIC$AIC)
  cyano_AIC$AIC <- round(cyano_AIC$AIC, digits = 2)
  
  bacteria_AIC <- corr_table_func_bt(som_maps = bacteria_ns, i = col)
  bacteria_AIC <- as.data.frame(bacteria_AIC, stringsAsFactors = FALSE)
  colnames(bacteria_AIC) <- c("Variables","AIC")
  bacteria_AIC$AIC <- as.numeric(bacteria_AIC$AIC)
  bacteria_AIC$AIC <- round(bacteria_AIC$AIC, digits = 2)
  
  eukaryota_AIC <- corr_table_func_bt(som_maps = eukaryota_ns, i = col)
  eukaryota_AIC <- as.data.frame(eukaryota_AIC, stringsAsFactors = FALSE)
  colnames(eukaryota_AIC) <- c("Variables","AIC")
  eukaryota_AIC$AIC <- as.numeric(eukaryota_AIC$AIC)
  eukaryota_AIC$AIC <- round(eukaryota_AIC$AIC, digits = 2)
  
  phyto_AIC  <- corr_table_func_bt(som_maps = phyto_ns, i = col)
  phyto_AIC <- as.data.frame(phyto_AIC, stringsAsFactors = FALSE)
  colnames(phyto_AIC) <- c("Variables","AIC")
  phyto_AIC$AIC <- as.numeric(phyto_AIC$AIC)
  phyto_AIC$AIC <- round(phyto_AIC$AIC, digits = 2)
  
  AIC_full <- full_join(cyano_AIC, bacteria_AIC, by = "Variables")
  AIC_full <- full_join(AIC_full, eukaryota_AIC, by = "Variables")
  NS_full <- full_join(AIC_full, phyto_AIC, by = "Variables")
  
  NS_full$Variables <- c("Temperature", "Salinity", "PO4ug",
                         "SiO3ug", "NO3ug","Chl-a", "Nutricline Depth",
                         "Distance to Coast", "Cyanobacteria Richness",
                         "Bacteria Richness", "Heterotrophic Eukaryote Richness",
                         "Eukaryotic Phytoplankton Richness", "Cyanobacteria Evenness",
                         "Bacteria Evenness", "Heterotrophic Eukaryote Evenness",
                         "Eukaryotic Phytoplankton Evenness", "Cyanobacteria Shannon",
                         "Bacteria Shannon", "Heterotrophic Eukaryote Shannon",
                         "Eukaryotic Phytoplankton Shannon")
  
  cyano_os <- cyano_plots %>% filter(som_id == "Offshore")
  bacteria_os <- bacteria_plots %>% filter(som_id == "Offshore")
  eukaryota_os <- eukaryota_plots %>% filter(som_id == "Offshore")
  phyto_os <- phyto_plots %>% filter(som_id == "Offshore")
  
  
  cyano_AIC <- corr_table_func_bt(som_maps = cyano_os, i = col)
  cyano_AIC <- as.data.frame(cyano_AIC, stringsAsFactors = FALSE)
  colnames(cyano_AIC) <- c("Variables","AIC")
  cyano_AIC$AIC <- as.numeric(cyano_AIC$AIC)
  cyano_AIC$AIC <- round(cyano_AIC$AIC, digits = 2)
  
  bacteria_AIC <- corr_table_func_bt(som_maps = bacteria_os, i = col)
  bacteria_AIC <- as.data.frame(bacteria_AIC, stringsAsFactors = FALSE)
  colnames(bacteria_AIC) <- c("Variables","AIC")
  bacteria_AIC$AIC <- as.numeric(bacteria_AIC$AIC)
  bacteria_AIC$AIC <- round(bacteria_AIC$AIC, digits = 2)
  
  eukaryota_AIC <- corr_table_func_bt(som_maps = eukaryota_os, i = col)
  eukaryota_AIC <- as.data.frame(eukaryota_AIC, stringsAsFactors = FALSE)
  colnames(eukaryota_AIC) <- c("Variables","AIC")
  eukaryota_AIC$AIC <- as.numeric(eukaryota_AIC$AIC)
  eukaryota_AIC$AIC <- round(eukaryota_AIC$AIC, digits = 2)
  
  phyto_AIC  <- corr_table_func_bt(som_maps = phyto_os, i = col)
  phyto_AIC <- as.data.frame(phyto_AIC, stringsAsFactors = FALSE)
  colnames(phyto_AIC) <- c("Variables","AIC")
  phyto_AIC$AIC <- as.numeric(phyto_AIC$AIC)
  phyto_AIC$AIC <- round(phyto_AIC$AIC, digits = 2)
  
  AIC_full <- full_join(cyano_AIC, bacteria_AIC, by = "Variables")
  AIC_full <- full_join(AIC_full, eukaryota_AIC, by = "Variables")
  OS_full <- full_join(AIC_full, phyto_AIC, by = "Variables")
  
  OS_full$Variables <- c("Temperature", "Salinity", "PO4ug",
                         "SiO3ug", "NO3ug","Chl-a", "Nutricline Depth",
                         "Distance to Coast", "Cyanobacteria Richness",
                         "Bacteria Richness", "Heterotrophic Eukaryote Richness",
                         "Eukaryotic Phytoplankton Richness", "Cyanobacteria Evenness",
                         "Bacteria Evenness", "Heterotrophic Eukaryote Evenness",
                         "Eukaryotic Phytoplankton Evenness", "Cyanobacteria Shannon",
                         "Bacteria Shannon", "Heterotrophic Eukaryote Shannon",
                         "Eukaryotic Phytoplankton Shannon")
  
  AIC_full <- NS_full
  AIC_full[,2:5] <- AIC_full[,2:5] - OS_full[,2:5]
  
  if(col == 1){
    NS_full[c(17),2] <- NA
    NS_full[c(18),3] <- NA
    NS_full[c(19),4] <- NA
    NS_full[c(20),5] <- NA
    
    OS_full[c(17),2] <- NA
    OS_full[c(18),3] <- NA
    OS_full[c(19),4] <- NA
    OS_full[c(20),5] <- NA
    
    AIC_full[c(17),2] <- NA
    AIC_full[c(18),3] <- NA
    AIC_full[c(19),4] <- NA
    AIC_full[c(20),5] <- NA
  }
  if(col == 2){  
    NS_full[c(13),2] <- NA
    NS_full[c(14),3] <- NA
    NS_full[c(15),4] <- NA
    NS_full[c(16),5] <- NA
    
    OS_full[c(13),2] <- NA
    OS_full[c(14),3] <- NA
    OS_full[c(15),4] <- NA
    OS_full[c(16),5] <- NA
    
    AIC_full[c(13),2] <- NA
    AIC_full[c(14),3] <- NA
    AIC_full[c(15),4] <- NA
    AIC_full[c(16),5] <- NA
  }
  if(col == 3){
    NS_full[c(9),2] <- NA
    NS_full[c(10),3] <- NA
    NS_full[c(11),4] <- NA
    NS_full[c(12),5] <- NA
    
    OS_full[c(9),2] <- NA
    OS_full[c(10),3] <- NA
    OS_full[c(11),4] <- NA
    OS_full[c(12),5] <- NA
    
    AIC_full[c(9),2] <- NA
    AIC_full[c(10),3] <- NA
    AIC_full[c(11),4] <- NA
    AIC_full[c(12),5] <- NA
  }
 
  
  colnames(AIC_full) <- c("Variables", "Cyanobacteria\n AIC",
                          "Bacteria/Archaea\n AIC", "Heterotrophic\nEukaryotes\n AIC", "Eukaryotic\nPhytoplankton\n AIC")
  
  colnames(NS_full) <- c("Variables", "Cyanobacteria\n AIC",
                          "Bacteria/Archaea\n AIC", "Heterotrophic\nEukaryotes\n AIC", "Eukaryotic\nPhytoplankton\n AIC")
  
  colnames(OS_full) <- c("Variables", "Cyanobacteria\n AIC",
                          "Bacteria/Archaea\n AIC", "Heterotrophic\nEukaryotes\n AIC", "Eukaryotic\nPhytoplankton\n AIC")
  
  colfunc <- colorRampPalette(c("blue","white","red"))
  
  cyano_col <- as.numeric(seq(-1,1,by = 0.01))
  scale <- colfunc(length(cyano_col))
  cyano_fill <- scale[match(round(AIC_full$`Cyanobacteria
 AIC`, digits = 2), round(cyano_col, digits = 2))]
  cyano_fill[which(is.na(cyano_fill))] <- "black"

  bacteria_col <- seq(-1,1,by = 0.01)
  scale <- colfunc(length(bacteria_col))
  bacteria_fill <- scale[match(round(AIC_full$`Bacteria/Archaea
 AIC`, digits = 2), round(cyano_col, digits = 2))]
  bacteria_fill[which(is.na(bacteria_fill))] <- "black"

  eukaryota_col <- seq(-1,1,by = 0.01)
  scale <- colfunc(length(eukaryota_col))
  eukaryota_fill <- scale[match(round(AIC_full$`Heterotrophic
Eukaryotes
 AIC`, digits = 2), round(cyano_col, digits = 2))]
  eukaryota_fill[which(is.na(eukaryota_fill))] <- "black"
  
  phyto_col <- seq(-1,1,by = 0.01)
  scale <- colfunc(length(phyto_col))
  phyto_fill <- scale[match(round(AIC_full$`Eukaryotic
Phytoplankton
 AIC`, digits = 2), round(cyano_col, digits = 2))]
  phyto_fill[which(is.na(phyto_fill))] <- "black"
  
  t0 <- tableGrob(AIC_full["Variables"], 
                  theme=ttheme_default(
                    core=list(bg_params = list(fill="grey90", col = "black"),
                              fg_params = list(fontface="bold")),
                    colhead = list(bg_params=list(fill="white", col="black"))), 
                  rows = NULL)
  t1 <- tableGrob(round(AIC_full["Cyanobacteria\n AIC"], digits = 2),
                  theme=ttheme_default(
                    core=list(bg_params = list(fill=cyano_fill, col = "black")),
                    colhead = list(bg_params=list(fill="white", col="black"))),
                  rows = NULL)
  t2 <- tableGrob(round(AIC_full["Bacteria/Archaea\n AIC"], digits = 2),
                  theme=ttheme_default(
                    core=list(bg_params = list(fill=bacteria_fill, col = "black")),
                    colhead = list(bg_params=list(fill="white", col="black"))),
                  rows = NULL)
  t3 <- tableGrob(round(AIC_full["Heterotrophic\nEukaryotes\n AIC"], digits = 2),
                  theme=ttheme_default(
                    core=list(bg_params = list(fill=eukaryota_fill, col = "black")),
                    colhead = list(bg_params=list(fill="white", col="black"))),
                  rows = NULL)
  t4 <- tableGrob(round(AIC_full["Eukaryotic\nPhytoplankton\n AIC"], digits = 2),
                  theme=ttheme_default(
                    core=list(bg_params = list(fill=phyto_fill, col = "black")),
                    colhead = list(bg_params=list(fill="white", col="black"))),
                  rows = NULL)
  
  # join tables
  tab <- gtable_combine(t0,t1,t2,t3,t4)
  
  pdf(file = figure_name, width = 9, height = 9)
  grid.arrange(tab)
  dev.off()
  
  AIC_scaled <- NS_full
  
  for (i in 2:ncol(AIC_scaled)) {
    zero_one_scale <- AIC_scaled[,i]
    AIC_scaled[,i] <- abs(zero_one_scale)
    
  }
  
  AIC_color <- AIC_full
  
  plot_df <- melt(AIC_scaled)
  color_df <- melt(AIC_color, id.vars = "Variables", value.name = "color")
  plot_df$color <- color_df$color
  
  
  
  colnames(plot_df) <- c("Variables", "Group", "Correlation", "Color")
  
  # # removing plastids from plot
  # plot_df <- plot_df[-which(plot_df$Group == "Eukaryotic\nPlastid\n AIC"),]
  
  plot_df$Variables <- as.factor(plot_df$Variables)
  plot_df$Variables <- factor(plot_df$Variables, levels = c("Distance to Coast","Eukaryotic Phytoplankton Evenness",
                                                            "Eukaryotic Phytoplankton Richness","Eukaryotic Phytoplankton Shannon",
                                                            "Heterotrophic Eukaryote Evenness","Heterotrophic Eukaryote Richness",
                                                            "Heterotrophic Eukaryote Shannon","Bacteria Evenness",
                                                            "Bacteria Richness","Bacteria Shannon",
                                                            "Cyanobacteria Evenness","Cyanobacteria Richness",
                                                            "Cyanobacteria Shannon","Nutricline Depth",
                                                            "Chl-a","SiO3ug",
                                                            "PO4ug","NO3ug",
                                                            "Salinity","Temperature"))
  
  pdf(figure_name_2, width = 8, height = 8)
  print(ggplot(data = plot_df, aes(x = Group, y = Variables, size = Correlation, fill = Color)) + 
          geom_point(color = "black", alpha = 0.6, shape = 21) +
          scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                               limits = c(-0.25,0.25), oob = scales::squish) +
          labs(size = "Correlation\n in Nearshore", fill = "Difference in\n Correlation") +
          ylab("Variable") +
          theme(panel.background = element_blank(),
                panel.border = element_rect(color = "black", fill = NA),
                legend.position = "right",
                panel.grid.major.y = element_line(color = "grey", linetype = 2),
                plot.title = element_text(hjust = 0.5),
                axis.text.x = element_text(angle = 0)) + ggtitle(title_name) +
          scale_size_continuous(range = c(1,10)) + xlab("") + ylab(""))
  dev.off()
  
  
}


# Nearshore

full_aic_table_figure_diversity(
  in_cyano = "output/cyano_16s_full_data.Rdata",
  in_bacteria = "output/bacteria_m_euks_16s_full_data.Rdata",
  in_euks = "output/euks_hetero_18sv9_full_data.Rdata",
  in_phyto = "output/euks_auto_18sv9_full_data.Rdata", 
  figure_name = paste0("figures/full_aic_table_ns_shan_",Sys.Date(),".pdf"),
  figure_name_2 = paste0("figures/full_aic_plot_ns_shan_",Sys.Date(),".pdf"),
  title_name = "R-Squared Nearshore Shannon", # col 1 = shan, col 2 = even col 3 = rich
  type = "Nearshore", color_fill = "red", col = 1)

full_aic_table_figure_diversity(
  in_cyano = "output/cyano_16s_full_data.Rdata",
  in_bacteria = "output/bacteria_m_euks_16s_full_data.Rdata",
  in_euks = "output/euks_hetero_18sv9_full_data.Rdata",
  in_phyto = "output/euks_auto_18sv9_full_data.Rdata", 
  figure_name = paste0("figures/full_aic_table_ns_even_",Sys.Date(),".pdf"),
  figure_name_2 = paste0("figures/full_aic_plot_ns_even_",Sys.Date(),".pdf"),
  title_name = "R-Squared Nearshore Evenness", # col 1 = shan, col 2 = even col 3 = rich
  type = "Nearshore", color_fill = "purple", col = 2)


full_aic_table_figure_diversity(
  in_cyano = "output/cyano_16s_full_data.Rdata",
  in_bacteria = "output/bacteria_m_euks_16s_full_data.Rdata",
  in_euks = "output/euks_hetero_18sv9_full_data.Rdata",
  in_phyto = "output/euks_auto_18sv9_full_data.Rdata", 
  figure_name = paste0("figures/full_aic_table_ns_rich_",Sys.Date(),".pdf"),
  figure_name_2 = paste0("figures/full_aic_plot_ns_rich_",Sys.Date(),".pdf"),
  title_name = "R-Squared Nearshore Richness", # col 1 = shan, col 2 = even col 3 = rich
  type = "Nearshore", color_fill = "blue", col = 3)

# Offshore

full_aic_table_figure_diversity(
  in_cyano = "output/cyano_16s_full_data.Rdata",
  in_bacteria = "output/bacteria_m_euks_16s_full_data.Rdata",
  in_euks = "output/euks_hetero_18sv9_full_data.Rdata",
  in_phyto = "output/euks_auto_18sv9_full_data.Rdata", 
  figure_name = paste0("figures/full_aic_table_os_shan_",Sys.Date(),".pdf"),
  figure_name_2 = paste0("figures/full_aic_plot_os_shan_",Sys.Date(),".pdf"),
  title_name = "R-Squared Offshore Shannon", # col 1 = shan, col 2 = even col 3 = rich
  type = "Offshore", color_fill = "red", col = 1)

full_aic_table_figure_diversity(
  in_cyano = "output/cyano_16s_full_data.Rdata",
  in_bacteria = "output/bacteria_m_euks_16s_full_data.Rdata",
  in_euks = "output/euks_hetero_18sv9_full_data.Rdata",
  in_phyto = "output/euks_auto_18sv9_full_data.Rdata", 
  figure_name = paste0("figures/full_aic_table_os_even_",Sys.Date(),".pdf"),
  figure_name_2 = paste0("figures/full_aic_plot_os_even_",Sys.Date(),".pdf"),
  title_name = "R-Squared Offshore Evenness", # col 1 = shan, col 2 = even col 3 = rich
  type = "Offshore", color_fill = "purple", col = 2)


full_aic_table_figure_diversity(
  in_cyano = "output/cyano_16s_full_data.Rdata",
  in_bacteria = "output/bacteria_m_euks_16s_full_data.Rdata",
  in_euks = "output/euks_hetero_18sv9_full_data.Rdata",
  in_phyto = "output/euks_auto_18sv9_full_data.Rdata", 
  figure_name = paste0("figures/full_aic_table_os_rich_",Sys.Date(),".pdf"),
  figure_name_2 = paste0("figures/full_aic_plot_os_rich_",Sys.Date(),".pdf"),
  title_name = "R-Squared Offshore Richness", # col 1 = shan, col 2 = even col 3 = rich
  type = "Offshore", color_fill = "blue", col = 3)


# Difference

full_aic_table_figure_diff(
  in_cyano = "output/cyano_16s_full_data.Rdata",
  in_bacteria = "output/bacteria_m_euks_16s_full_data.Rdata",
  in_euks = "output/euks_hetero_18sv9_full_data.Rdata",
  in_phyto = "output/euks_auto_18sv9_full_data.Rdata", 
  figure_name = paste0("figures/full_table_diff_shan_",Sys.Date(),".pdf"),
  figure_name_2 = paste0("figures/full_table_plot_diff_shan_",Sys.Date(),".pdf"),
  title_name = "Correlation Shannon Diversity", # col 1 = shan, col 2 = even col 3 = rich
  col = 1)

full_aic_table_figure_diff(
  in_cyano = "output/cyano_16s_full_data.Rdata",
  in_bacteria = "output/bacteria_m_euks_16s_full_data.Rdata",
  in_euks = "output/euks_hetero_18sv9_full_data.Rdata",
  in_phyto = "output/euks_auto_18sv9_full_data.Rdata", 
  figure_name = paste0("figures/full_table_diff_rich_",Sys.Date(),".pdf"),
  figure_name_2 = paste0("figures/full_table_plot_diff_rich_",Sys.Date(),".pdf"),
  title_name = "Correlation Richness", # col 1 = shan, col 2 = even col 3 = rich
  col = 3)

full_aic_table_figure_diff(
  in_cyano = "output/cyano_16s_full_data.Rdata",
  in_bacteria = "output/bacteria_m_euks_16s_full_data.Rdata",
  in_euks = "output/euks_hetero_18sv9_full_data.Rdata",
  in_phyto = "output/euks_auto_18sv9_full_data.Rdata", 
  figure_name = paste0("figures/full_table_diff_even_",Sys.Date(),".pdf"),
  figure_name_2 = paste0("figures/full_table_plot_diff_even_",Sys.Date(),".pdf"),
  title_name = "Correlation Evenness", # col 1 = shan, col 2 = even col 3 = rich
  col = 2)

############## Additional Plots #####

r1 <- ggplot(cyano_plots, aes(x = bact_rich, y = richness, color = som_id)) + geom_point() + geom_smooth(method = "glm") + xlab("Bacteria Richness") +
  ylab("Cyanobacteria Richness") + labs(color = "Cluster")

r2 <- ggplot(cyano_plots, aes(x = eukary_rich, y = richness, color = som_id)) + geom_point() + geom_smooth(method = "glm") + xlab("Heterotrophic Eukaryote Richness") +
  ylab("Cyanobacteria Richness") + labs(color = "Cluster")

r3 <- ggplot(cyano_plots, aes(x = phyto_rich, y = richness, color = som_id)) + geom_point() + geom_smooth(method = "glm") + xlab("Eukaryotic Phytoplankton Richness") +
  ylab("Cyanobacteria Richness") + labs(color = "Cluster")

nc <- ggplot(cyano_plots, aes(x = NCDepth, y = richness, color = som_id)) + geom_point() + geom_smooth(method = "glm") + xlab("Nutricline Depth (m)") +
  ylab("Cyanobacteria Richness") + labs(color = "Cluster")

plot_grid(r1,r2,r3,nc)