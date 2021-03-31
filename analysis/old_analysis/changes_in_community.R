library(tidyverse)
library(sp)
library(spatialEco)
library(RColorBrewer)

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9",
                  "bacteria_m_euks_16s","cyano_16s","archaea_16s", "euks_hetero_18sv9", "euks_auto_18sv9")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales", "Rhodobacterales",
                   "Sar Clade", "Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans",
                   "Bacteria", "Cyanobacteria", "Archaea", "Eukaryotic Protists",
                   "Photosynthetic Eukaryotic Protists")

in_group_list_basic = c("16s_pro", "16s_syne","16s_flavo", "16s_rhodo", "16s_sar", 
                        "18s_diatom","18s_dino", "18s_syndin",
                        "18s_hapto", "18s_chloro", "18s_metazoa", 
                        "16s_bacteria_m_euks", "16s_cyanos", "16s_archaea",
                        "18s_heterotrophic_euks", "18s_autotrophic_euks")

type_list = c("16s", "16s","16s", "16s", "16s", 
                        "18s","18s", "18s",
                        "18s", "18s", "18s", 
                        "16s", "16s", "16s",
                        "18s", "18s")
# big groups

in_group_list = c("bacteria_m_euks_16s","cyano_16s","archaea_16s", "euks_hetero_18sv9", "euks_auto_18sv9")

in_group_names = c("Bacteria", "Cyanobacteria", "Archaea", "Heterotrophic Eukaryotic Protists",
                   "Photosynthetic Eukaryotic Protists")

in_group_list_basic = c("16s_bacteria_m_euks", "16s_cyanos", "16s_archaea",
                        "18s_heterotrophic_euks", "18s_autotrophic_euks")

type_list = c("16s", "16s", "16s",
              "18s", "18s")

depth_list <- c(3,6,4,3,3)
vjust_list <- c(-12, -8, -5, -7, -5)
hjust_list <- c(0, 0, 1, 0, 0)
vjust_list2 <- c(1.75, 0, 0, 0, 0)

out_plot_list <- list()

for (i in 1:length(in_group_list)) {
  
  load(paste0("output/", in_group_list[i],"_full_data.Rdata"))
  load(paste0("output/", in_group_list[i],"_map.Rdata"))
  load(paste0("data/", in_group_list_basic[i],".Rdata"))
  
  centroid_df <- SpatialPointsDataFrame(coords = som_maps[,c(6,5)], data = som_maps)
  
  wt_1 <- wt.centroid(x = centroid_df, p = 2)
  wt_2 <- wt.centroid(x = centroid_df, p = 3)
  
  clust1 <- which.max(c(wt_1@coords[1], wt_2@coords[1]))
  clust2 <- which.min(c(wt_1@coords[1], wt_2@coords[1]))
  
  if(clust1 == 1){
    full_dat$som_id[which(full_dat$som_id == 1)] <- "Nearshore"
    full_dat$som_id[which(full_dat$som_id == 2)] <- "Offshore"
  }
  
  if(clust1 == 2){
    full_dat$som_id[which(full_dat$som_id == 1)] <- "Offshore"
    full_dat$som_id[which(full_dat$som_id == 2)] <- "Nearshore"
  }
  
  
  scaled_inputs <- as.data.frame(scaled_inputs)
  
  scaled_inputs <- scaled_inputs[which(!is.na(match(rownames(scaled_inputs), paste0("X",full_dat$Sample.Name)))),]
  
  scaled_inputs$som_id <- full_dat$som_id[match(rownames(scaled_inputs), paste0("X",full_dat$Sample.Name))]
  scaled_inputs$sample <- full_dat$Sample.Name[match(rownames(scaled_inputs), paste0("X",full_dat$Sample.Name))]
  
  asv_long <- scaled_inputs %>%
    pivot_longer(-c(som_id,sample), names_to = "ASV", values_to = "reads")
  
  # Offshore - Nearshore
  
  asv_diff <- asv_long %>%
    group_by(ASV, som_id) %>%
    summarise(mean = mean(reads, na.rm = TRUE)) %>%
    mutate(diff=c(NA,diff(mean))) %>%
    filter(!is.na(diff)) %>%
    arrange(desc(diff))
  
  asv_diff$preference <- "Offshore"
  asv_diff$preference[which(asv_diff$diff < 0)] <- "Nearshore"
  asv_diff$preference[which(asv_diff$diff == 0)] <- "Equal"
  
  asv_diff$preference <- as.factor(asv_diff$preference)
  asv_diff$preference <- factor(asv_diff$preference, levels = c("Nearshore", "Offshore", "Equal"))
  
  quants <- quantile(asv_diff$diff, probs = c(0.01,0.99))
  
  asv_most_diff <- asv_diff %>%
    filter(diff < quants[1] | diff > quants[2])
  
  if(type_list[i] == "16s"){
    six_tax_id$taxon <- six_tax_id$Silva_Taxon
    tax_id <- six_tax_id
    taxa <- tax_id$taxon[match(asv_most_diff$ASV,tax_id$Feature.ID)]
    taxa <- strsplit(taxa, ";")
    taxa <- sapply(taxa, "[", depth_list[i])
    taxa <- sub(paste0("D_",depth_list[i]-1,"__"), "", taxa)   
  }
  
  if(type_list[i] == "18s"){
    eight_tax_id$taxon <- eight_tax_id$PR2_Taxon
    tax_id <- eight_tax_id
    taxa <- tax_id$taxon[match(asv_most_diff$ASV,tax_id$Feature.ID)]
    taxa <- strsplit(taxa, ";")
    taxa <- sapply(taxa, "[", depth_list[i])
  }
  
  asv_most_diff$taxa <- taxa
  
  taxa_table <- asv_most_diff %>%
    group_by(preference, taxa) %>%
    tally()
  
  taxa_table$taxa[which(is.na(taxa_table$taxa))] <- "unclassified"
  
  taxa_table$taxa <- as.factor(taxa_table$taxa)

  # Define the number of colors you want
  nb.cols <- length(unique(taxa_table$taxa))
  mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
  
  
  p1 <- ggplot(asv_most_diff, aes(x = reorder(ASV,-diff), y = diff, fill = taxa)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
              fill = "red", alpha = 0.15) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = -Inf,
              fill = "blue", alpha = 0.15) +
    geom_bar(stat = "identity", color = "black", show.legend = FALSE, lwd = 0.25) +
    xlab(paste0(in_group_names[i]," ASVs")) + 
    scale_fill_manual(values = mycolors, drop = FALSE) + 
    scale_color_manual(values = mycolors, drop = FALSE) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black")) +
    ylab(expression(Delta*" Relative Abundance")) +
    labs(fill = "Groups", color = "Groups") +
    geom_hline(yintercept = 0, lty = 2, lwd = 1) +
    ggtitle(in_group_names[i])
  
  os <- ggplot(taxa_table %>% filter(preference == "Offshore"), aes(x="", y=n, fill=taxa)) +
    geom_bar(stat="identity", width=1, show.legend = FALSE, color = "black") +
    scale_fill_manual(values = mycolors, drop = FALSE) +
    coord_polar("y", start=0) + theme_void() +
    ggtitle("Offshore") + theme(plot.title = element_text(vjust = vjust_list[i], hjust = 0.5))
   
  ns <- ggplot(taxa_table %>% filter(preference == "Nearshore"), aes(x="", y=n, fill=taxa)) +
    geom_bar(stat="identity", width=1, color = "black") +
    scale_fill_manual(values = mycolors, drop = FALSE) +
    coord_polar("y", start=0) + theme_void() + labs(fill = "Groups") +
    ggtitle("Nearshore") + theme(plot.title = element_text(vjust = vjust_list[i], hjust = 0.5))
    
  layout <- c(
    area(t = 1, l = 1, b = 5, r = 6),
    area(t = (1 + vjust_list2[i]), l = (3.5 + hjust_list[i]), b = (2 + vjust_list2[i]), r = (4.5 + hjust_list[i])),
    area(t = 4, l = (3.5), b = 5, r = (4.5))
    
  )
  
  out_plot_list[[i]] <- p1 + os + ns + 
    plot_layout(design = layout, guides = "collect")

  pdf(file = paste0("figures/", in_group_list[i],"_som_comm_difs.pdf"), width = 10, height = 8)
  print(out_plot_list[[i]])
  dev.off()
    
}                                                                            


