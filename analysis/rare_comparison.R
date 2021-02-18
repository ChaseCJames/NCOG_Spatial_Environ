library(tidyverse)
library(patchwork)

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "bacteria_m_euks_16s", "plastid_16s", "cyano_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9", "euks_hetero_18sv9", "euks_auto_18sv9")

in_group_names = c("A. Prochlorococcus", "B. Synecococcus", "C. Flavobacteriales", "D. Rhodobacterales",
                   "E. Sar Clade", "F. Archaea", "G. Bacteria", "H. Eukaryotic Phytoplankton (Plastids)",
                   "I. Cyanobacteria", "J. Diatoms", "K. Dinoflagellates", "L. Syndiniales",
                   "M. Haptophytes", "N. Chlorophytes","O. Metazoans", "P. Eukaryotic Protists", "Q. Photosynthetic Eukaryotic Protists")

plot_list <- list()

for (i in 1:length(in_group_list)) {
  
  
  load(paste0("output/",in_group_list[i],"_full_data.Rdata"))
  
  full_dat <- full_dat %>%
    filter(!is.na(rare_richness))
  
  sha <- ggplot(full_dat, aes(x = shannon, y = rare_shannon)) +
    geom_point() + stat_smooth(method = "lm", color = "red", fill = "red") +
    xlab("Shannon Diversity") + ylab("Rarified Shannon Diversity") +
    theme(aspect.ratio = 1,
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black")) +
    xlim(min(full_dat$shannon, na.rm = TRUE),max(full_dat$shannon)) +
    ylim(min(full_dat$shannon, na.rm = TRUE),max(full_dat$shannon)) +
    ggtitle(in_group_names[i])
  
  rich <- ggplot(full_dat, aes(x = richness, y = rare_richness)) +
    geom_point() + stat_smooth(method = "lm", color = "green", fill = "green") +
    xlab("Richness") + ylab("Rarified Richness") +
    theme(aspect.ratio = 1,
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black")) +
    xlim(min(full_dat$richness, na.rm = TRUE),max(full_dat$richness)) +
    ylim(min(full_dat$richness, na.rm = TRUE),max(full_dat$richness))
  
  even <- ggplot(full_dat, aes(x = evenness, y = rare_evenness)) +
    geom_point() + stat_smooth(method = "lm", color = "blue", fill = "blue") +
    xlab("Evenness") + ylab("Rarified Evenness") +
    theme(aspect.ratio = 1,
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black")) +
    xlim(min(full_dat$evenness, na.rm = TRUE),max(full_dat$evenness)) +
    ylim(min(full_dat$evenness, na.rm = TRUE),max(full_dat$evenness))
  
  plot_list[[i]] <- list(plot1 = sha, plot2 = rich, plot3 = even)
  
}

out <- plot_list[[1]]$plot1 + plot_list[[1]]$plot2 + plot_list[[1]]$plot3 +
  plot_list[[2]]$plot1 + plot_list[[2]]$plot2 + plot_list[[2]]$plot3 +
  plot_list[[3]]$plot1 + plot_list[[3]]$plot2 + plot_list[[3]]$plot3 +
  plot_list[[4]]$plot1 + plot_list[[4]]$plot2 + plot_list[[4]]$plot3 +
  plot_list[[5]]$plot1 + plot_list[[5]]$plot2 + plot_list[[5]]$plot3 +
  plot_list[[6]]$plot1 + plot_list[[6]]$plot2 + plot_list[[6]]$plot3 +
  plot_list[[7]]$plot1 + plot_list[[7]]$plot2 + plot_list[[7]]$plot3 +
  plot_list[[8]]$plot1 + plot_list[[8]]$plot2 + plot_list[[8]]$plot3 +
  plot_list[[9]]$plot1 + plot_list[[9]]$plot2 + plot_list[[9]]$plot3 +
  plot_list[[10]]$plot1 + plot_list[[10]]$plot2 + plot_list[[10]]$plot3 +
  plot_list[[11]]$plot1 + plot_list[[11]]$plot2 + plot_list[[11]]$plot3 +
  plot_list[[12]]$plot1 + plot_list[[12]]$plot2 + plot_list[[12]]$plot3 +
  plot_list[[13]]$plot1 + plot_list[[13]]$plot2 + plot_list[[13]]$plot3 +
  plot_list[[14]]$plot1 + plot_list[[14]]$plot2 + plot_list[[14]]$plot3 +
  plot_list[[15]]$plot1 + plot_list[[15]]$plot2 + plot_list[[15]]$plot3 +
  plot_list[[16]]$plot1 + plot_list[[16]]$plot2 + plot_list[[16]]$plot3 +
  plot_list[[17]]$plot1 + plot_list[[17]]$plot2 + plot_list[[17]]$plot3 +
  plot_layout(ncol = 3)


pdf(file = "figures/rare_compare.pdf", width = 10, height = 45)  
  print(out)
  dev.off()
  
  
  

