library(tidyverse)
library(sp)
library(spatialEco)
library(RColorBrewer)
library(patchwork)

# big groups

in_group_list = c("bacteria_m_euks_16s","cyano_16s","archaea_16s", "euks_hetero_18sv9", "euks_auto_18sv9")

in_group_names = c("Bacteria", "Cyanobacteria", "Archaea", "Heterotrophic Eukaryotic Protists",
                   "Photosynthetic Eukaryotic Protists")

in_group_list_basic = c("16s_bacteria_m_euks", "16s_cyanos", "16s_archaea",
                        "18s_heterotrophic_euks", "18s_autotrophic_euks")


out_plot <- list()
split_plot <- list()
tsize = 22

###### bacteria #######

i = 1
  
  load(paste0("output/", in_group_list[i],"_full_data_S.Rdata"))
  load(paste0("output/", in_group_list[i],"_map_S.Rdata"))
  load(paste0("data/", in_group_list_basic[i],"_S.Rdata"))
  
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
  
  # sorting out ASVs
  
  six_tax_id$taxon <- six_tax_id$silva_Taxon
  tax_id <- six_tax_id
  taxa <- tax_id$taxon[match(asv_long$ASV,tax_id$Feature.ID)]
  taxa <- strsplit(taxa, ";")
  taxa_sha <- sapply(taxa, "[", 3)
  taxa_sha <- sub(paste0(" c__"), "", taxa_sha)   
  taxa_deep <- sapply(taxa, "[", 4)
  taxa_deep <- sub(paste0(" o__"), "", taxa_deep)  
  
  taxa_sha[which(taxa_deep == "Flavobacteriales")] <- "Flavobacteriales"
  taxa_sha[which(taxa_deep == "Rhodobacterales")] <- "Rhodobacterales"
  taxa_sha[which(taxa_deep == "SAR11_clade")] <- "SAR 11 clade"
  
  unique_tax <- unique(taxa_sha)
  
  others <- unique_tax[c(4:5,17:18,20,24,27,36:40,
                         42,44,50:51,53,55,57,59:60,
                         67,69,71,74:76,79,82:85,88,91,
                         98:99,104)]
  
  taxa_sha[which(taxa_sha %in% c(others, NA))] <- "Other Bacteria"
  
  asv_long$Group <- taxa_sha
  
  asv_mean_group <- asv_long %>%
    group_by(sample, Group, som_id) %>%
    summarise(grp_reads = sum(reads, na.rm = TRUE), n = n()) %>%
    group_by(som_id, Group) %>%
    summarise(mean_grp_reads = mean(grp_reads, na.rm = TRUE), mean_n = mean(n, na.rm = TRUE))
  
  asv_diff <- asv_long %>%
    group_by(sample, Group, som_id) %>%
    summarise(grp_reads = sum(reads, na.rm = TRUE), n = n()) %>%
    group_by(Group, som_id) %>%
    summarise(mean = mean(grp_reads, na.rm = TRUE)) %>%
    mutate(diff=c(NA,diff(mean))) %>%
    filter(!is.na(diff)) %>%
    arrange(desc(diff))

  order <- asv_diff$Group
  
  asv_mean_group$Group[which(is.na(asv_mean_group$Group))] <- "unclassified"
  
  asv_mean_group$Group <- factor(asv_mean_group$Group, levels = order[length(order):1])
  
  asv_mean_group$size <- (scale(asv_mean_group$mean_n)+1)
  
  offshore <- asv_mean_group %>% filter(som_id == "Offshore")
  offshore <- offshore[match(order, offshore$Group),]
  offshore$rank <- 1:nrow(offshore)
  offshore$diff <- asv_diff$diff[match(offshore$Group, asv_diff$Group)]
  
  nearshore <- asv_mean_group %>% filter(som_id == "Nearshore")
  nearshore <- nearshore[match(order, nearshore$Group),]
  nearshore$rank <- 1:nrow(nearshore)
  nearshore$diff <- asv_diff$diff[match(nearshore$Group, asv_diff$Group)]
  
  
  os <- ggplot(offshore, aes(x=mean_grp_reads, y=Group, color = diff, size = mean_n)) +
    geom_point(show.legend = TRUE) +
    geom_segment(aes(y=Group, 
                     yend=Group, 
                     x=0, 
                     xend=mean_grp_reads),size = 0.1, show.legend = FALSE) +
    scale_color_gradient2(low = "blue", high = "red", mid = "grey") +
    scale_size_continuous(range = c(0, 8)) +
    ylab("Rank") + xlab("Mean Relative Abundance") +
    ggtitle("Offshore") + theme(plot.title = element_text(vjust = 0, hjust = 0.5)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) +
    xlim(c(0,max(c(offshore$mean_grp_reads,nearshore$mean_grp_reads))))
  
  ns <- ggplot(nearshore, aes(x=mean_grp_reads, y= Group, color=diff, size = mean_n)) +
    geom_point(show.legend = FALSE) +
    geom_segment(aes(y=Group, 
                     yend=Group, 
                     x=0, 
                     xend=mean_grp_reads),size = 0.1, show.legend = FALSE) +
    scale_color_gradient2(low = "blue", high = "red", mid = "grey") +
    scale_size_continuous(range = c(0, 8)) + scale_x_reverse() +
    ylab("Rank") + xlab("Mean Relative Abundance") +
    ggtitle("Nearshore") + theme(plot.title = element_text(vjust = 0, hjust = 0.5)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.key = element_rect(fill = NA, color = NA)) +
    labs(color = expression(Delta*" Relative Abundance"),
         size = "# of ASVs") +
    xlim(c(max(c(offshore$mean_grp_reads,nearshore$mean_grp_reads)),0))
  
  out_plot[[i]] <-  ns + os + plot_layout(guides = "collect")
  
  pdf(file = paste0("figures/", in_group_list[i],"_som_comm_S.pdf"), width = 18, height = 8)
  print(out_plot[[i]])
  dev.off()
  
  split_plot[[i]] <- list(ns = ns, os = os)
  
  ###### cyanobacteria #######
  
  i = 2
  
  load(paste0("output/", in_group_list[i],"_full_data_S.Rdata"))
  load(paste0("output/", in_group_list[i],"_map_S.Rdata"))
  load(paste0("data/", in_group_list_basic[i],"_S.Rdata"))
  
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
  
  # sorting out ASVs
  
  six_tax_id$taxon <- six_tax_id$silva_Taxon
  tax_id <- six_tax_id
  taxa <- tax_id$taxon[match(asv_long$ASV,tax_id$Feature.ID)]
  taxa <- strsplit(taxa, ";")
  taxa_sha <- sapply(taxa, "[", 5)
  taxa_sha <- sub(paste0(" f__"), "", taxa_sha)   
  taxa_deep <- sapply(taxa, "[", 6)
  taxa_deep <- sub(paste0(" g__"), "", taxa_deep)  
  
  taxa_sha[which(taxa_deep == "Prochlorococcus_MIT9313")] <- "Prochlorococcus MIT9313"
  taxa_sha[which(taxa_deep == "Synechococcus_CC9902")] <- "Synechococcus CC9902"
  taxa_sha[which(taxa_deep == "Trichodesmium_IMS101")] <- "Trichodesmium IMS101"
  taxa_sha[which(taxa_deep == "Cyanobium_PCC-6307")] <- "Cyanobium PCC-6307"
  taxa_sha[which(taxa_deep == "Richelia_HH01")] <- "Richelia HH01"
  taxa_sha[which(taxa_deep == "Atelocyanobacterium_(UCYN-A)")] <- "Atelocyanobacterium (UCYN-A)"
  taxa_sha[which(taxa_deep == "Tychonema_CCAP_1459-11B")] <- "Tychonema CCAP 1459-11B"
  taxa_sha[which(taxa_deep == "Phormidesmis_ANT.L52.6")] <- "Phormidesmis ANT.L52.6"
  taxa_sha[which(taxa_deep == "Annamia_HOs24")] <- "Annamia HOs24"
  taxa_sha[which(taxa_deep == "Aliterella")] <- "Aliterella"
  taxa_sha[which(taxa_deep == "Candidatus_Obscuribacter")] <- "Candidatus Obscuribacter"
  
  unique_tax <- unique(taxa_sha)
  
  asv_long$Group <- taxa_sha
  
  asv_mean_group <- asv_long %>%
    group_by(sample, Group, som_id) %>%
    summarise(grp_reads = sum(reads, na.rm = TRUE), n = n()) %>%
    group_by(som_id, Group) %>%
    summarise(mean_grp_reads = mean(grp_reads, na.rm = TRUE), mean_n = mean(n, na.rm = TRUE))
  
  asv_diff <- asv_long %>%
    group_by(sample, Group, som_id) %>%
    summarise(grp_reads = sum(reads, na.rm = TRUE), n = n()) %>%
    group_by(Group, som_id) %>%
    summarise(mean = mean(grp_reads, na.rm = TRUE)) %>%
    mutate(diff=c(NA,diff(mean))) %>%
    filter(!is.na(diff)) %>%
    arrange(desc(diff))
  
  asv_mean_group$Group[which(is.na(asv_mean_group$Group))] <- "unclassified"
  asv_diff$Group[which(is.na(asv_diff$Group))] <- "unclassified"
  
  order <- asv_diff$Group
  
  asv_mean_group$Group <- factor(asv_mean_group$Group, levels = order[length(order):1])
  
  asv_mean_group$size <- (scale(asv_mean_group$mean_n)+1)
  
  offshore <- asv_mean_group %>% filter(som_id == "Offshore")
  offshore <- offshore[match(order, offshore$Group),]
  offshore$rank <- 1:nrow(offshore)
  offshore$diff <- asv_diff$diff[match(offshore$Group, asv_diff$Group)]
  
  nearshore <- asv_mean_group %>% filter(som_id == "Nearshore")
  nearshore <- nearshore[match(order, nearshore$Group),]
  nearshore$rank <- 1:nrow(nearshore)
  nearshore$diff <- asv_diff$diff[match(nearshore$Group, asv_diff$Group)]
  
  
  os <- ggplot(offshore, aes(x=mean_grp_reads, y=Group, color = diff, size = mean_n)) +
    geom_point(show.legend = TRUE) +
    geom_segment(aes(y=Group, 
                     yend=Group, 
                     x=0, 
                     xend=mean_grp_reads),size = 0.1, show.legend = FALSE) +
    scale_color_gradient2(low = "blue", high = "red", mid = "grey") +
    scale_size_continuous(range = c(0, 8)) +
    ylab("Rank") + xlab("Mean Relative Abundance") +
    ggtitle("Offshore") + theme(plot.title = element_text(vjust = 0, hjust = 0.5)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) +
    xlim(c(0,max(c(offshore$mean_grp_reads,nearshore$mean_grp_reads))))
  
  ns <- ggplot(nearshore, aes(x=mean_grp_reads, y= Group, color=diff, size = mean_n)) +
    geom_point(show.legend = FALSE) +
    geom_segment(aes(y=Group, 
                     yend=Group, 
                     x=0, 
                     xend=mean_grp_reads),size = 0.1, show.legend = FALSE) +
    scale_color_gradient2(low = "blue", high = "red", mid = "grey") +
    scale_size_continuous(range = c(0, 8)) + scale_x_reverse() +
    ylab("Rank") + xlab("Mean Relative Abundance") +
    ggtitle("Nearshore") + theme(plot.title = element_text(vjust = 0, hjust = 0.5)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.key = element_rect(fill = NA, color = NA)) +
    labs(color = expression(Delta*" Relative Abundance"),
         size = "# of ASVs") +
    xlim(c(max(c(offshore$mean_grp_reads,nearshore$mean_grp_reads)),0))
  
  out_plot[[i]] <-  ns + os + plot_layout(guides = "collect")
  
  pdf(file = paste0("figures/", in_group_list[i],"_som_comm_S.pdf"), width = 18, height = 8)
  print(out_plot[[i]])
  dev.off()
  
  split_plot[[i]] <- list(ns = ns, os = os)
  
  ####### Eukaryotic Phytoplankton ##############
  
  i = 5
  
  load(paste0("output/", in_group_list[i],"_full_data_S.Rdata"))
  load(paste0("output/", in_group_list[i],"_map_S.Rdata"))
  load(paste0("data/", in_group_list_basic[i],"_S.Rdata"))
  
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
  
  # sorting out ASVs
  
  eight_tax_id$taxon <- eight_tax_id$pr2_Taxon
  tax_id <- eight_tax_id
  taxa <- tax_id$taxon[match(asv_long$ASV,tax_id$Feature.ID)]
  taxa <- strsplit(taxa, ";")
  taxa_sha <- sapply(taxa, "[", 3)
  taxa_deep <- sapply(taxa, "[", 4)
  
  taxa_sha[which(taxa_deep == "Mamiellophyceae")] <- "Mamiellophyceae"
  taxa_sha[which(taxa_deep == "Bacillariophyta")] <- "Bacillariophyta"
  
  unique_tax <- unique(taxa_sha)
  
  asv_long$Group <- taxa_sha
  
  asv_mean_group <- asv_long %>%
    group_by(sample, Group, som_id) %>%
    summarise(grp_reads = sum(reads, na.rm = TRUE), n = n()) %>%
    group_by(som_id, Group) %>%
    summarise(mean_grp_reads = mean(grp_reads, na.rm = TRUE), mean_n = mean(n, na.rm = TRUE))
  
  asv_diff <- asv_long %>%
    group_by(sample, Group, som_id) %>%
    summarise(grp_reads = sum(reads, na.rm = TRUE), n = n()) %>%
    group_by(Group, som_id) %>%
    summarise(mean = mean(grp_reads, na.rm = TRUE)) %>%
    mutate(diff=c(NA,diff(mean))) %>%
    filter(!is.na(diff)) %>%
    arrange(desc(diff))
  
  order <- asv_diff$Group
  
  asv_mean_group$Group[which(is.na(asv_mean_group$Group))] <- "unclassified"
  
  asv_mean_group$Group <- as.factor(asv_mean_group$Group)
  
  asv_mean_group$Group <- factor(asv_mean_group$Group, levels = order[length(order):1])
  
  asv_mean_group$size <- (scale(asv_mean_group$mean_n)+1)
  
  offshore <- asv_mean_group %>% filter(som_id == "Offshore")
  offshore <- offshore[match(order, offshore$Group),]
  offshore$rank <- 1:nrow(offshore)
  offshore$diff <- asv_diff$diff[match(offshore$Group, asv_diff$Group)]
  
  nearshore <- asv_mean_group %>% filter(som_id == "Nearshore")
  nearshore <- nearshore[match(order, nearshore$Group),]
  nearshore$rank <- 1:nrow(nearshore)
  nearshore$diff <- asv_diff$diff[match(nearshore$Group, asv_diff$Group)]
  
  
  os <- ggplot(offshore, aes(x=mean_grp_reads, y=Group, color = diff, size = mean_n)) +
    geom_point(show.legend = TRUE) +
    geom_segment(aes(y=Group, 
                     yend=Group, 
                     x=0, 
                     xend=mean_grp_reads),size = 0.1, show.legend = FALSE) +
    scale_color_gradient2(low = "blue", high = "red", mid = "grey") +
    scale_size_continuous(range = c(0, 8)) +
    ylab("Rank") + xlab("Mean Relative Abundance") +
    ggtitle("Offshore") + theme(plot.title = element_text(vjust = 0, hjust = 0.5)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) +
    xlim(c(0,max(c(offshore$mean_grp_reads,nearshore$mean_grp_reads))))
  
  ns <- ggplot(nearshore, aes(x=mean_grp_reads, y= Group, color=diff, size = mean_n)) +
    geom_point(show.legend = FALSE) +
    geom_segment(aes(y=Group, 
                     yend=Group, 
                     x=0, 
                     xend=mean_grp_reads),size = 0.1, show.legend = FALSE) +
    scale_color_gradient2(low = "blue", high = "red", mid = "grey") +
    scale_size_continuous(range = c(0, 8)) + scale_x_reverse() +
    ylab("Rank") + xlab("Mean Relative Abundance") +
    ggtitle("Nearshore") + theme(plot.title = element_text(vjust = 0, hjust = 0.5)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.key = element_rect(fill = NA, color = NA)) +
    labs(color = expression(Delta*" Relative Abundance"),
         size = "# of ASVs") +
    xlim(c(max(c(offshore$mean_grp_reads,nearshore$mean_grp_reads)),0))
  
  out_plot[[i]] <-  ns + os + plot_layout(guides = "collect")
  
  pdf(file = paste0("figures/", in_group_list[i],"_som_comm_S.pdf"), width = 18, height = 8)
  print(out_plot[[i]])
  dev.off()
  
  split_plot[[i]] <- list(ns = ns, os = os)
  
  ####### Archaea ######  
  
  i = 3
  
  load(paste0("output/", in_group_list[i],"_full_data_S.Rdata"))
  load(paste0("output/", in_group_list[i],"_map_S.Rdata"))
  load(paste0("data/", in_group_list_basic[i],"_S.Rdata"))
  
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
  
  # sorting out ASVs
  
  six_tax_id$taxon <- six_tax_id$silva_Taxon
  tax_id <- six_tax_id
  taxa <- tax_id$taxon[match(asv_long$ASV,tax_id$Feature.ID)]
  taxa <- strsplit(taxa, ";")
  taxa_sha <- sapply(taxa, "[", 2)
  taxa_sha <- sub(paste0(" p__"), "", taxa_sha)   
  taxa_deep <- sapply(taxa, "[", 3)
  taxa_deep <- sub(paste0(" c__"), "", taxa_deep)  
  
  taxa_sha[which(taxa_deep == "Thermoplasmata")] <- "Thermoplasmata"
  taxa_sha[which(taxa_deep == "Nitrososphaeria")] <- "Nitrososphaeria"
  taxa_sha[which(taxa_deep == "Halobacteria")] <- "Halobacteria"
  taxa_sha[which(taxa_deep == "Altiarchaeia")] <- "Altiarchaeia"

  
  unique_tax <- unique(taxa_sha)
  
  others <- unique_tax[c(4)]
  
  taxa_sha[which(taxa_sha %in% others)] <- "Unclassified Archaea"
  
  asv_long$Group <- taxa_sha
  
  asv_mean_group <- asv_long %>%
    group_by(sample, Group, som_id) %>%
    summarise(grp_reads = sum(reads, na.rm = TRUE), n = n()) %>%
    group_by(som_id, Group) %>%
    summarise(mean_grp_reads = mean(grp_reads, na.rm = TRUE), mean_n = mean(n, na.rm = TRUE))
  
  asv_diff <- asv_long %>%
    group_by(sample, Group, som_id) %>%
    summarise(grp_reads = sum(reads, na.rm = TRUE), n = n()) %>%
    group_by(Group, som_id) %>%
    summarise(mean = mean(grp_reads, na.rm = TRUE)) %>%
    mutate(diff=c(NA,diff(mean))) %>%
    filter(!is.na(diff)) %>%
    arrange(desc(diff))
  
  order <- asv_diff$Group
  
  asv_mean_group$Group <- factor(asv_mean_group$Group, levels = order[55:1])
  
  asv_mean_group$size <- (scale(asv_mean_group$mean_n)+1)
  
  offshore <- asv_mean_group %>% filter(som_id == "Offshore")
  offshore <- offshore[match(order, offshore$Group),]
  offshore$rank <- 1:nrow(offshore)
  offshore$diff <- asv_diff$diff[match(offshore$Group, asv_diff$Group)]
  
  nearshore <- asv_mean_group %>% filter(som_id == "Nearshore")
  nearshore <- nearshore[match(order, nearshore$Group),]
  nearshore$rank <- 1:nrow(nearshore)
  nearshore$diff <- asv_diff$diff[match(nearshore$Group, asv_diff$Group)]
  
  
  os <- ggplot(offshore, aes(x=mean_grp_reads, y=Group, color = diff, size = mean_n)) +
    geom_point(show.legend = TRUE) +
    geom_segment(aes(y=Group, 
                     yend=Group, 
                     x=0, 
                     xend=mean_grp_reads),size = 0.1, show.legend = FALSE) +
    scale_color_gradient2(low = "blue", high = "red", mid = "grey") +
    scale_size_continuous(range = c(0, 8)) +
    ylab("Rank") + xlab("Mean Relative Abundance") +
    ggtitle("Offshore") + theme(plot.title = element_text(vjust = 0, hjust = 0.5)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) +
    xlim(c(0,max(c(offshore$mean_grp_reads,nearshore$mean_grp_reads))))
  
  ns <- ggplot(nearshore, aes(x=mean_grp_reads, y= Group, color=diff, size = mean_n)) +
    geom_point(show.legend = FALSE) +
    geom_segment(aes(y=Group, 
                     yend=Group, 
                     x=0, 
                     xend=mean_grp_reads),size = 0.1, show.legend = FALSE) +
    scale_color_gradient2(low = "blue", high = "red", mid = "grey") +
    scale_size_continuous(range = c(0, 8)) + scale_x_reverse() +
    ylab("Rank") + xlab("Mean Relative Abundance") +
    ggtitle("Nearshore") + theme(plot.title = element_text(vjust = 0, hjust = 0.5)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.key = element_rect(fill = NA, color = NA)) +
    labs(color = expression(Delta*" Relative Abundance"),
         size = "# of ASVs") +
    xlim(c(max(c(offshore$mean_grp_reads,nearshore$mean_grp_reads)),0))
  
  out_plot[[i]] <-  ns + os + plot_layout(guides = "collect")
  
  pdf(file = paste0("figures/", in_group_list[i],"_som_comm_S.pdf"), width = 18, height = 8)
  print(out_plot[[i]])
  dev.off()
  
  split_plot[[i]] <- list(ns = ns, os = os)
  
  ###### Heterotrophic Euks ######
  
  i = 4
  
  load(paste0("output/", in_group_list[i],"_full_data_S.Rdata"))
  load(paste0("output/", in_group_list[i],"_map_S.Rdata"))
  load(paste0("data/", in_group_list_basic[i],"_S.Rdata"))
  
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
  
  # sorting out ASVs
  
  eight_tax_id$taxon <- eight_tax_id$pr2_Taxon
  tax_id <- eight_tax_id
  taxa <- tax_id$taxon[match(asv_long$ASV,tax_id$Feature.ID)]
  taxa <- strsplit(taxa, ";")
  taxa_sha <- sapply(taxa, "[", 4)
  
  unique_tax <- unique(taxa_sha)
  
  others <- unique_tax[c(1,71)]
  
  taxa_sha[which(taxa_sha %in% others)] <- "Other Eukaryotic Protists"
  
  asv_long$Group <- taxa_sha
  
  asv_mean_group <- asv_long %>%
    group_by(sample, Group, som_id) %>%
    summarise(grp_reads = sum(reads, na.rm = TRUE), n = n()) %>%
    group_by(som_id, Group) %>%
    summarise(mean_grp_reads = mean(grp_reads, na.rm = TRUE), mean_n = mean(n, na.rm = TRUE))
  
  asv_diff <- asv_long %>%
    group_by(sample, Group, som_id) %>%
    summarise(grp_reads = sum(reads, na.rm = TRUE), n = n()) %>%
    group_by(Group, som_id) %>%
    summarise(mean = mean(grp_reads, na.rm = TRUE)) %>%
    mutate(diff=c(NA,diff(mean))) %>%
    filter(!is.na(diff)) %>%
    arrange(desc(diff))
  
  order <- asv_diff$Group
  
  asv_mean_group$Group <- factor(asv_mean_group$Group, levels = order[length(order):1])
  
  asv_mean_group$size <- (scale(asv_mean_group$mean_n)+1)
  
  offshore <- asv_mean_group %>% filter(som_id == "Offshore")
  offshore <- offshore[match(order, offshore$Group),]
  offshore$rank <- 1:nrow(offshore)
  offshore$diff <- asv_diff$diff[match(offshore$Group, asv_diff$Group)]

  
  nearshore <- asv_mean_group %>% filter(som_id == "Nearshore")
  nearshore <- nearshore[match(order, nearshore$Group),]
  nearshore$rank <- 1:nrow(nearshore)
  nearshore$diff <- asv_diff$diff[match(nearshore$Group, asv_diff$Group)]

  
  
  os <- ggplot(offshore, aes(x=mean_grp_reads, y=Group, color = diff, size = mean_n)) +
    geom_point(show.legend = TRUE) +
    geom_segment(aes(y=Group, 
                     yend=Group, 
                     x=0, 
                     xend=mean_grp_reads),size = 0.1, show.legend = FALSE) +
    scale_color_gradient2(low = "blue", high = "red", mid = "grey") +
    scale_size_continuous(range = c(0, 8)) +
    ylab("Rank") + xlab("Mean Relative Abundance") +
    ggtitle("Offshore") + theme(plot.title = element_text(vjust = 0, hjust = 0.5)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank()) +
    xlim(c(0,max(c(offshore$mean_grp_reads,nearshore$mean_grp_reads))))
  
  ns <- ggplot(nearshore, aes(x=mean_grp_reads, y= Group, color=diff, size = mean_n)) +
    geom_point(show.legend = FALSE) +
    geom_segment(aes(y=Group, 
                     yend=Group, 
                     x=0, 
                     xend=mean_grp_reads),size = 0.1, show.legend = FALSE) +
    scale_color_gradient2(low = "blue", high = "red", mid = "grey") +
    scale_size_continuous(range = c(0, 8)) + scale_x_reverse() +
    ylab("Rank") + xlab("Mean Relative Abundance") +
    ggtitle("Nearshore") + theme(plot.title = element_text(vjust = 0, hjust = 0.5)) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.key = element_rect(fill = NA, color = NA)) +
    labs(color = expression(Delta*" Relative Abundance"),
         size = "# of ASVs") +
    xlim(c(max(c(offshore$mean_grp_reads,nearshore$mean_grp_reads)),0))
  
  out_plot[[i]] <-  ns + os + plot_layout(guides = "collect")
  
  pdf(file = paste0("figures/", in_group_list[i],"_som_comm.pdf"), width = 18, height = 8)
  print(out_plot[[i]])
  dev.off()
  
  split_plot[[i]] <- list(ns = ns, os = os)
  
  layout <- "ABCD
             ABCD
             ABCD
             EFGH
             IJKK"
  
  plot_final <- (split_plot[[1]]$ns + theme(text = element_text(size = tsize))) +
    (split_plot[[1]]$os + labs(size = "# of ASVs",color = expression(Delta*" Relative Abundance")) + 
       theme(text = element_text(size = tsize))) +
    (split_plot[[4]]$ns + theme(text = element_text(size = tsize))) +
    (split_plot[[4]]$os + labs(size = "# of ASVs",color = expression(Delta*" Relative Abundance")) + 
       theme(text = element_text(size = tsize))) +
    (split_plot[[2]]$ns + theme(text = element_text(size = tsize))) +
    (split_plot[[2]]$os + labs(size = "# of ASVs",color = expression(Delta*" Relative Abundance")) + 
       theme(text = element_text(size = tsize))) +
    (split_plot[[3]]$ns + theme(text = element_text(size = tsize))) +
    (split_plot[[3]]$os + labs(size = "# of ASVs",color = expression(Delta*" Relative Abundance")) + 
       theme(text = element_text(size = tsize))) +
    (split_plot[[5]]$ns + theme(text = element_text(size = tsize))) +
    (split_plot[[5]]$os + labs(size = "# of ASVs",color = expression(Delta*" Relative Abundance")) + 
       theme(text = element_text(size = tsize))) +
    plot_spacer() + plot_layout(design = layout)
  
  
  pdf("figures_S/supp_fig_3_S.pdf", width = 25, height = 28)
  print(plot_final)
  dev.off()
  
  