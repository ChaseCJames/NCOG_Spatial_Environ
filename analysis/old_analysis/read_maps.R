library(tidyverse)
library(ggmap)

map <- map_data("world")   
tsize = 12

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9", "bacteria_m_euks_16s",
                  "plastid_16s", "cyano_16s", "euks_hetero_18sv9", "euks_auto_18sv9")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                   "Sar 11 Clade", "Archaea","Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans", "Bacteria",
                   "Eukaryotic Phytoplankton (Plastids)", "Cyanobacteria",
                   "Eukaryotic Protists", "Photosynthetic Eukaryotic Protists")

in_group_list_basic = c("16s_pro", "16s_syne","16s_flavo", "16s_rhodo", "16s_sar", "16s_archaea",
                        "18s_diatom","18s_dino", "18s_syndin",
                        "18s_hapto", "18s_chloro", "18s_metazoa", "16s_bacteria_m_euks",
                        "16s_plastids", "16s_cyanos", "18s_heterotrophic_euks", "18s_autotrophic_euks")


for (i in 1:length(in_group_list)) {
  
  load(paste0("data/",in_group_list_basic[i],".Rdata"))
  load(paste0("output/",in_group_list[i],"_map.Rdata"))
  
  total_df <- as.data.frame(matrix(NA,nrow = nrow(asv_table), ncol = 2))
  colnames(total_df) <- c("sample", "reads")
  
  total_df$reads <- rowSums(asv_table, na.rm = TRUE)
  
  total_df$sample <- gsub("_"," ",substr(rownames(asv_table),9,19))
  
  station_reads <- total_df %>% group_by(sample) %>% summarise(mean_reads = mean(reads, na.rm = TRUE))
  
  som_maps$reads <- station_reads$mean_reads[match(som_maps$Sta_ID, station_reads$sample)]
  
  print(ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes(x = long, y = lat, fill = reads), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "darkred") +
    ggtitle(in_group_names[i]) + labs(fill = "Mean Reads") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          axis.line = element_blank(),
          axis.text = element_text(size = tsize),
          legend.text = element_text(size = tsize),
          axis.title = element_text(size = tsize)))
  
  
}
