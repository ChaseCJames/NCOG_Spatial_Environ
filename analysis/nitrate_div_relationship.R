library(tidyverse)
library(patchwork)

in_group_list = c("bacteria_m_euks_16s", "cyano_16s", "archaea_16s",
                  "euks_hetero_18sv9", "euks_auto_18sv9")

in_group_names = c( "Bacteria", "Cyanobacteria","Archaea", 
                   "Heterotrophic eukaryotic protists", "Photosynthetic eukaryotic protists")

in_group_list_basic = c("16s_bacteria_m_euks", "16s_cyanos","16s_archaea",
                        "18s_heterotrophic_euks", "18s_autotrophic_euks")


load("output/pro_16s_full_data.Rdata")
load("output/upwelling_indicies.Rdata")



cruise <- full_dat %>%
  group_by(Cruise) %>%
  summarise(mean_nitrate = mean(NO3ug, na.rm = TRUE))

cruise$Date <- paste0("01-",substr(cruise$Cruise,5,6),"-",substr(cruise$Cruise,1,4))
cruise$Date <- as.Date(cruise$Date, format = "%d-%m-%Y")
cruise$Phase <- c(rep("Warm",12), rep("Cool",8), rep("2019",4))
cruise$Season <- rep(c("Winter", "Spring", "Summer", "Fall"),6)

index_mat <- matrix(NA,nrow = nrow(cruise), 3)

for (i in 1:nrow(cruise)) {
  
  year <- substr(cruise$Cruise[i],1,4)
  month <- substr(cruise$Cruise[i],5,6) 
  
  index_yr <- substr(index_vals_mat$Date,1,4)
  index_month <- substr(index_vals_mat$Date,6,7)
  
  index_vals <- which(match(index_yr,year) & match(index_month,month))
  
  index_mat[i,] <- colMeans(index_vals_mat[index_vals,2:4], na.rm = TRUE)
  
}

cruise$CUTI <- index_mat[,1]
cruise$BEUTI <- index_mat[,2]
cruise$Nitrate <- index_mat[,3]

n_plot <- ggplot() +
  geom_line(data = cruise, aes(x = Date, y = BEUTI), color = "black", lwd = 1, lty = 2) + 
  geom_point(data = cruise, aes(x = Date, y = BEUTI, fill = Phase, color = Phase, shape = Season), size = 3) + 
  ggtitle("Biologically Effective Upwelling Index (BEUTI)") + 
  ylab("BEUTI") +
  scale_fill_manual(values = c("gold3", "blue", "red")) +
  scale_color_manual(values = c("gold3", "blue", "red")) +
  scale_shape_manual(values = c(21,22,23,24)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"))

plot_list <- list()

cruise_mat <- as.data.frame(matrix(NA,24,6))

for (i in 1:length(in_group_list)) {
  
  load(paste0("output/", in_group_list[i], "_full_data.Rdata"))
  
  dat <- full_dat %>%
    group_by(Cruise) %>%
    summarise(mean_shannon = mean(shannon, na.rm = TRUE))
  
  dat$Date <- paste0("01-",substr(dat$Cruise,5,6),"-",substr(dat$Cruise,1,4))
  dat$Date <- as.Date(dat$Date, format = "%d-%m-%Y")
  dat$Phase <- c(rep("Warm",12), rep("Cool",8), rep("2019",4))
  dat$Season <- rep(c("Winter", "Spring", "Summer", "Fall"),6)
  
  plot_list[[i]] <- ggplot() +
    geom_line(data = dat, aes(x = Date, y = mean_shannon), color = "black", lwd = 1, lty = 2) + 
    geom_point(data = dat, aes(x = Date, y = mean_shannon, fill = Phase, color = Phase, shape = Season), size = 3) + 
    ggtitle(in_group_names[i]) + 
    ylab("Mean Alpha Diversity") +
    scale_fill_manual(values = c("gold3", "blue", "red")) +
    scale_color_manual(values = c("gold3", "blue", "red")) +
    scale_shape_manual(values = c(21,22,23,24)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"))
  
  cruise_mat[,i] <- dat$mean_shannon
  
  }

cruise_mat$V6 <- cruise$Date
colnames(cruise_mat) <- c("Bacteria","Cyanobacteria","Archaea", "Heterotrophic_eukaryotic_protists", "Photosynthetic_eukaryotic_protists", "Date")
cruise_mat$Phase <- cruise$Phase
cruise_mat$Season <- cruise$Season
cruise_mat$mean_nitrate <- cruise$mean_nitrate
cruise_mat$BEUTI <- cruise$BEUTI
cruise_mat$CUTI <- cruise$CUTI
cruise_mat$Nitrate <- cruise$Nitrate

long_plot <- cruise_mat %>%
  pivot_longer(-c(Date,Phase,Season, mean_nitrate, BEUTI, CUTI, Nitrate), names_to = "Group", values_to = "mean_alpha")

long_plot$Group <- gsub(long_plot$Group, pattern = "_", replacement = " ")

names <- c("Bacteria","Cyanobacteria","Archaea", "Heterotrophic eukaryotic protists", "Photosynthetic eukaryotic protists")

p_val <- vector()

for (i in 1:length(names)) {
  
  dat <- long_plot %>% filter(Group == names[i])
  
  glm_out <- summary(glm(mean_alpha ~ BEUTI, data = dat))
  
  p_val[i] <- glm_out$coefficients[2,4]
  
}

sig <- which(p_val < 0.05)

side_plot <- ggplot() +
  geom_point(data = long_plot, aes(x = BEUTI, y = mean_alpha, color = Group, fill = Group, shape = Season), size = 3) +
  stat_smooth(data = long_plot %>% filter(Group %in% names[sig]),
              aes(x = BEUTI, y = mean_alpha, color = Group, fill = Group), method = "glm", show.legend = FALSE) +
  scale_color_manual(values = c("#ba583b","#50b47b","#7f63b8","#9da140","#b84c7d")) +
  scale_fill_manual(values = c("#ba583b","#50b47b","#7f63b8","#9da140","#b84c7d")) + xlab("Biologically Effective Upwelling Index\nPer Cruise") +
  ylab("Mean Alpha Diversity\n Per Cruise") + scale_shape_manual(values = c(21,22,23,24)) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"))


out_plot <- (n_plot + plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] + plot_list[[5]] +
  plot_layout(nrow = 6, guides = "collect")) | side_plot 

pdf(file = "figures/cce_plot_BEUTI.pdf", width = 18, height = 9)
print(out_plot)
dev.off()


# all groups

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                   "Sar Clade", "Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans", "Bacteria",
                   "Eukaryotic Phytoplankton (Plastids)", "Cyanobacteria",
                   "Eukaryotic Protists", "Photosynthetic Eukaryotic Protists")

in_group_list_basic = c("16s_pro", "16s_syne","16s_flavo", "16s_rhodo", "16s_sar", 
                        "18s_diatom","18s_dino", "18s_syndin",
                        "18s_hapto", "18s_chloro", "18s_metazoa")


load("output/upwelling_indicies.Rdata")
load("output/pro_16s_full_data.Rdata")

cruise <- full_dat %>%
  group_by(Cruise) %>%
  summarise(mean_nitrate = mean(NO3ug, na.rm = TRUE))

cruise$Date <- paste0("01-",substr(cruise$Cruise,5,6),"-",substr(cruise$Cruise,1,4))
cruise$Date <- as.Date(cruise$Date, format = "%d-%m-%Y")
cruise$Phase <- c(rep("Warm",12), rep("Cool",8), rep("2019",4))
cruise$Season <- rep(c("Winter", "Spring", "Summer", "Fall"),6)

index_mat <- matrix(NA,nrow = nrow(cruise), 3)

for (i in 1:nrow(cruise)) {
  
  year <- substr(cruise$Cruise[i],1,4)
  month <- substr(cruise$Cruise[i],5,6) 
  
  index_yr <- substr(index_vals_mat$Date,1,4)
  index_month <- substr(index_vals_mat$Date,6,7)
  
  index_vals <- which(match(index_yr,year) & match(index_month,month))
  
  index_mat[i,] <- colMeans(index_vals_mat[index_vals,2:4], na.rm = TRUE)
  
}

cruise$CUTI <- index_mat[,1]
cruise$BEUTI <- index_mat[,2]
cruise$Nitrate <- index_mat[,3]

plot_list <- list()

cruise_mat <- as.data.frame(matrix(NA,24,12))

for (i in 1:length(in_group_list)) {
  
  load(paste0("output/", in_group_list[i], "_full_data.Rdata"))
  
  dat <- full_dat %>%
    group_by(Cruise) %>%
    summarise(mean_shannon = mean(shannon, na.rm = TRUE))
  
  dat$Date <- paste0("01-",substr(dat$Cruise,5,6),"-",substr(dat$Cruise,1,4))
  dat$Date <- as.Date(dat$Date, format = "%d-%m-%Y")
  dat$Phase <- c(rep("Warm",12), rep("Cool",8), rep("2019",4))
  dat$Season <- rep(c("Winter", "Spring", "Summer", "Fall"),6)
  
  plot_list[[i]] <- ggplot() +
    geom_line(data = dat, aes(x = Date, y = mean_shannon), color = "black", lwd = 1, lty = 2) + 
    geom_point(data = dat, aes(x = Date, y = mean_shannon, fill = Phase, color = Phase, shape = Season), size = 3) + 
    ggtitle(in_group_names[i]) + 
    ylab("Mean Alpha Diversity") +
    scale_fill_manual(values = c("gold3", "blue", "red")) +
    scale_color_manual(values = c("gold3", "blue", "red")) +
    scale_shape_manual(values = c(21,22,23,24)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"))
  
  cruise_mat[,i] <- dat$mean_shannon
  
}

cruise_mat$V12 <- cruise$Date
colnames(cruise_mat) <- c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                          "SAR_11_Clade", "Diatoms", "Dinoflagellates", "Syndiniales",
                          "Haptophytes", "Chlorophytes","Metazoans", "Date")
cruise_mat$Phase <- cruise$Phase
cruise_mat$Season <- cruise$Season
cruise_mat$mean_nitrate <- cruise$mean_nitrate
cruise_mat$BEUTI <- cruise$BEUTI
cruise_mat$CUTI <- cruise$CUTI
cruise_mat$Nitrate <- cruise$Nitrate

long_plot <- cruise_mat %>%
  pivot_longer(-c(Date,Phase,Season, mean_nitrate, BEUTI, CUTI, Nitrate), names_to = "Group", values_to = "mean_alpha")

long_plot$Group <- gsub(long_plot$Group, pattern = "_", replacement = " ")

names <- c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
           "SAR 11 Clade", "Diatoms", "Dinoflagellates", "Syndiniales",
           "Haptophytes", "Chlorophytes","Metazoans")

p_val <- vector()

for (i in 1:length(names)) {
  
  dat <- long_plot %>% filter(Group == names[i])
  
  gam_out <- summary(gam(mean_alpha ~ BEUTI, data = dat))
  
  p_val[i] <- gam_out$pTerms.pv
  
}

sig <- which(p_val < 0.05)

side_plot <- ggplot() +
  geom_point(data = long_plot, aes(x = BEUTI, y = mean_alpha, color = Group, fill = Group, shape = Season), size = 3) +
  stat_smooth(data = long_plot %>% filter(Group %in% names[sig]),
              aes(x = BEUTI, y = mean_alpha, color = Group, fill = Group), method = "gam", show.legend = FALSE) +
  scale_color_manual(values = c("#9a0035","#169d36","#3c2c9a","#ffca51","#01b9fd","#c13820",
                                "#510042","#664600","#ff6492","#f09583","#e98aaa")) +
  scale_fill_manual(values = c("#9a0035","#169d36","#3c2c9a","#ffca51","#01b9fd","#c13820",
                               "#510042","#664600","#ff6492","#f09583","#e98aaa")) + xlab("Biologically Effective Upwelling Index\nPer Cruise") +
  ylab("Mean Alpha Diversity\n Per Cruise") + scale_shape_manual(values = c(21,22,23,24)) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"))


out_plot <- side_plot 

pdf(file = "figures/cce_plot_taxa_BEUTI.pdf", width = 10, height = 8)
print(out_plot)
dev.off()








