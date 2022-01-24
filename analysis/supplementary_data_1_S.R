library(tidyverse)
library(sp)
library(spatialEco)
library(RColorBrewer)
library(gridExtra)
library(grid)
# 
# in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
#                   "diatom_18sv9","dino_18sv9", "syndin_18sv9",
#                   "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9",
#                   "bacteria_m_euks_16s","cyano_16s","archaea_16s", "euks_hetero_18sv9", "euks_auto_18sv9")
# 
# in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales", "Rhodobacterales",
#                    "Sar Clade", "Diatoms", "Dinoflagellates", "Syndiniales",
#                    "Haptophytes", "Chlorophytes","Metazoans",
#                    "Bacteria", "Cyanobacteria", "Archaea", "Eukaryotic Protists",
#                    "Photosynthetic Eukaryotic Protists")
# 
# in_group_list_basic = c("16s_pro", "16s_syne","16s_flavo", "16s_rhodo", "16s_sar", 
#                         "18s_diatom","18s_dino", "18s_syndin",
#                         "18s_hapto", "18s_chloro", "18s_metazoa", 
#                         "16s_bacteria_m_euks", "16s_cyanos", "16s_archaea",
#                         "18s_heterotrophic_euks", "18s_autotrophic_euks")
# 
# type_list = c("16s", "16s","16s", "16s", "16s", 
#                         "18s","18s", "18s",
#                         "18s", "18s", "18s", 
#                         "16s", "16s", "16s",
#                         "18s", "18s")
# big groups

in_group_list = c("bacteria_m_euks_16s","cyano_16s","archaea_16s", "euks_hetero_18sv9", "euks_auto_18sv9")

in_group_names = c("Bacteria", "Cyanobacteria", "Archaea", "Heterotrophic Eukaryotic Protists",
                   "Photosynthetic Eukaryotic Protists")

in_group_list_basic = c("16s_bacteria_m_euks", "16s_cyanos", "16s_archaea",
                        "18s_heterotrophic_euks", "18s_autotrophic_euks")

type_list = c("16s", "16s", "16s",
              "18s", "18s")

out_list <- list()

for (i in 1:length(in_group_list)) {
  
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
    six_tax_id$taxon <- six_tax_id$silva_Taxon
    tax_id <- six_tax_id
    taxa <- tax_id$taxon[match(asv_most_diff$ASV,tax_id$Feature.ID)]
  }
  
  if(type_list[i] == "18s"){
    eight_tax_id$taxon <- eight_tax_id$pr2_Taxon
    tax_id <- eight_tax_id
    taxa <- tax_id$taxon[match(asv_most_diff$ASV,tax_id$Feature.ID)]
  }
  
  asv_most_diff$taxa <- taxa
  asv_most_diff$Group <- in_group_names[i]
  
  asv_most_diff <- asv_most_diff[,4:7]

  out_list[[i]] <- asv_most_diff
    
}                                                                            

all_groups <- bind_rows(out_list[[1]], out_list[[2]])

all_groups <- bind_rows(all_groups, out_list[[3]])
all_groups <- bind_rows(all_groups, out_list[[4]])
all_groups <- bind_rows(all_groups, out_list[[5]])

colnames(all_groups) <- c("Delta Relative Abundance", "Dominant Cluster", "Taxonomy", "Major Group")

write.csv(x = all_groups, file = "output/supplementary_data_1_diff_abun_S.csv")


