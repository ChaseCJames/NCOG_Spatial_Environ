library(tidyverse)
library(rworldmap)
library(ggplot2)
library(ggmap)
library(cowplot)
library(scales)

station_cosmo_endemics <- function(in_data = "data/16s_cyanos.Rdata",
                             output = "output/cyano_16s_cosmo_end.Rdata",
                             type = "six", type2 = NULL){

  load(in_data)
  
splits <- strsplit(rownames(asv_table), "_")

stations <- vector()

for (i in 1:length(splits)) {
  stations[i] <- paste0(splits[[i]][2]," ",splits[[i]][3])
}

asv_table$station <- stations

asv_gather <- gather(asv_table, "ASV", "Count", -station)

station_sums <- asv_gather %>%
  group_by(station, ASV) %>%
  summarise(total_count = sum(Count,na.rm = TRUE))

scaled_inputs_df <- as.data.frame(scaled_inputs)
scaled_inputs_df$station <- stations

scale_gather <- gather(scaled_inputs_df, "ASV", "Count", -station)

station_rel <- scale_gather %>%
  group_by(station, ASV) %>%
  summarise(mean_rel = mean(Count,na.rm = TRUE))


wide_station <- pivot_wider(station_sums, names_from = ASV, values_from = total_count)

# Cosmopolitans

n_zeros <- vector()
for (i in 2:ncol(wide_station)) {
  n_zeros[i-1] <-  length(which(wide_station[,i] == 0))
}

names(n_zeros) <- colnames(wide_station)[-1]

cosmos <- which(n_zeros == 0)

mean_rel_abun_cosmos <- colMeans(scaled_inputs)[cosmos]

if(type == "six"){
six_tax_id$Feature.ID[which(!is.na(match(six_tax_id$Feature.ID, names(cosmos))))]
cosmo_names <- six_tax_id$Silva_Taxon[which(!is.na(match(six_tax_id$Feature.ID, names(cosmos))))]
}

if(type == "eight"){
  eight_tax_id$Feature.ID[which(!is.na(match(eight_tax_id$Feature.ID, names(cosmos))))]
  if(type2 == "auto"){
    cosmo_names <- eight_tax_id$Taxon_PR2[which(!is.na(match(eight_tax_id$Feature.ID, names(cosmos))))]
  }
  if(type2 == "hetero"){
    cosmo_names <- eight_tax_id$Taxon_Silva[which(!is.na(match(eight_tax_id$Feature.ID, names(cosmos))))]
  }
}

station_cosmos <- station_rel[which(!is.na(match(station_rel$ASV, names(cosmos)))),]

# endemics

endemics <- which(n_zeros == 76)

if(type == "six"){
  six_tax_id$Feature.ID[which(!is.na(match(six_tax_id$Feature.ID, names(endemics))))]
  end_names <- six_tax_id$Silva_Taxon[which(!is.na(match(six_tax_id$Feature.ID, names(endemics))))]
}

if(type == "eight"){
  eight_tax_id$Feature.ID[which(!is.na(match(eight_tax_id$Feature.ID, names(endemics))))]
  if(type2 == "auto"){
    end_names <- eight_tax_id$Taxon_PR2[which(!is.na(match(eight_tax_id$Feature.ID, names(endemics))))]
  }
  if(type2 == "hetero"){
    end_names <- eight_tax_id$Taxon_Silva[which(!is.na(match(eight_tax_id$Feature.ID, names(endemics))))]
  }
}

stations <- vector()

for (i in 1:length(endemics)) {
  
 stations[i] <- wide_station$station[which(wide_station[,(endemics[i]+1)] > 0)]
  
}

station_table <- table(stations)[order(table(stations), decreasing = TRUE)]

save(cosmos, cosmo_names, mean_rel_abun_cosmos, station_cosmos,
     endemics, end_names, station_table, file = output)

}


station_cosmo_endemics(in_data = "data/16s_cyanos.Rdata",
                       output = "output/cyano_16s_cosmo_end.Rdata",
                       type = "six", type2 = NULL)

station_cosmo_endemics(in_data = "data/16s_bacteria_m_euks.Rdata",
                       output = "output/bacter_m_euks_16s_cosmo_end.Rdata",
                       type = "six", type2 = NULL)

station_cosmo_endemics(in_data = "data/18s_autotrophic_euks.Rdata",
                       output = "output/euks_auto_18sv9_cosmo_end.Rdata",
                       type = "eight", type2 = "auto")

station_cosmo_endemics(in_data = "data/18s_heterotrophic_euks.Rdata",
                       output = "output/euks_hetero_18sv9_cosmo_end.Rdata",
                       type = "eight", type2 = "hetero")



cosmo_endemic_figs <- function(in_cosmo = "output/cyano_16s_cosmo_end.Rdata",
                               in_map = "output/cyano_16s_map.Rdata", 
                               name = "Cyanobacteria"){
  
  load(in_cosmo)
  load(in_map)
  
  sta_cosmo <- station_cosmos[-which(is.na(match(station_cosmos$station, som_maps$Sta_ID))),]
  
  sta_cosmo$Lat <- som_maps$lat[match(sta_cosmo$station, som_maps$Sta_ID)]
  sta_cosmo$Lon <- som_maps$long[match(sta_cosmo$station, som_maps$Sta_ID)]
  
  map <- map_data("world")   
  
  uni_cosmo <- unique(sta_cosmo$ASV)
  
  cosmo_vect <- list()
  
  for (i in 1:length(uni_cosmo)) {
    
    sub <- sta_cosmo %>% filter(ASV == uni_cosmo[i])
    
    max_cosmo <- max(sub$mean_rel)
    
    plots <- ggplot() + 
      geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
      coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
      xlab("Longitude") + ylab("Latitude") + 
      geom_point(data = sta_cosmo %>% filter(ASV == uni_cosmo[i]), aes(x = Lon, y = Lat, fill = mean_rel), color = "black", size =6, stroke = 0.1, shape = 21) +
      scale_fill_gradient(low = "pink", high = "red", limits = c(0,max_cosmo), oob = scales::squish) +
      theme(legend.title = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
            plot.title = element_text(hjust = 0.5), axis.line = element_blank()) +
      ggtitle(paste0(name,"\nCosmopolitan ASV #",i)) 
    
    print(plots)
    
    cosmo_vect[[i]] <- plots
    
  }

  som_maps$Endemics <- as.vector(station_table[match(som_maps$Sta_ID, names(station_table))])
  som_maps$Endemics[is.na(som_maps$Endemics)] <- 0
  
  som_maps$Endemic_Samp <- som_maps$Endemics/som_maps$n_samps
  
  endemic_plot <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes(x = long, y = lat, fill = Endemics), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "pink", high = "red", limits = c(0,(max(station_table))), oob = scales::squish) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank()) +
    ggtitle(paste0(name, "\n# of Endemics by Station")) 
  
  endemic_samp_plot <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes(x = long, y = lat, fill = Endemic_Samp), color = "black", size =6, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "red") +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA,colour = "black", linetype = "solid", size = 1),
          plot.title = element_text(hjust = 0.5), axis.line = element_blank()) +
    ggtitle(paste0(name, "\n# of Endemics by Station / # Samples")) 
  
  
  print(endemic_plot)
  print(endemic_samp_plot)
  
  return(list(cosmo_vect = cosmo_vect, endemic_plot = endemic_plot, sample_plot = endemic_samp_plot))
  
}

cyano_out <- cosmo_endemic_figs(in_cosmo = "output/cyano_16s_cosmo_end.Rdata",
                   in_map = "output/cyano_16s_map.Rdata", name = "Cyanobacteria")

bact_out <- cosmo_endemic_figs(in_cosmo = "output/bacter_m_euks_16s_cosmo_end.Rdata",
                   in_map = "output/bacteria_m_euks_16s_map.Rdata",
                   name = "Bacteria/Archaea")

auto_out <- cosmo_endemic_figs(in_cosmo = "output/euks_auto_18sv9_cosmo_end.Rdata",
                   in_map = "output/euks_auto_18sv9_map.Rdata",
                   name = "Eukaryotic Phytoplankton")

hetero_out <- cosmo_endemic_figs(in_cosmo = "output/euks_hetero_18sv9_cosmo_end.Rdata",                 
                   in_map = "output/euks_hetero_18sv9_map.Rdata",
                   name = "Heterotrophic Eukaryotes")


pdf(file = "figures/endemics_fig.pdf", width = 8, height = 8)
plot_grid(auto_out$endemic_plot, hetero_out$endemic_plot,
          cyano_out$endemic_plot, bact_out$endemic_plot,
          ncol = 2 , nrow = 2)
dev.off()

pdf(file = "figures/endemics_samp_fig.pdf", width = 8, height = 8)
plot_grid(auto_out$sample_plot, hetero_out$sample_plot,
          cyano_out$sample_plot, bact_out$sample_plot,
          ncol = 2 , nrow = 2)
dev.off()

load("output/cyano_16s_cosmo_end.Rdata")
cosmo_df <- as.data.frame(matrix(nrow = 1,ncol = 2))
cosmo_df$V1 <- cosmo_names
cosmo_df$V2 <- mean_rel_abun_cosmos
write.csv(cosmo_df, file = "output/cyano_16s_cosmo_names.csv")

load("output/bacter_m_euks_16s_cosmo_end.Rdata")
cosmo_df <- as.data.frame(matrix(nrow = 33,ncol = 2))
cosmo_df$V1 <- cosmo_names
cosmo_df$V2 <- mean_rel_abun_cosmos
write.csv(cosmo_df, file = "output/bacter_m_euks_16s_cosmo_names.csv")

load("output/euks_auto_18sv9_cosmo_end.Rdata")
cosmo_df <- as.data.frame(matrix(nrow = 10,ncol = 2))
cosmo_df$V1 <- cosmo_names
cosmo_df$V2 <- mean_rel_abun_cosmos
write.csv(cosmo_df, file = "output/euks_auto_18sv9_cosmo_names.csv")

load("output/euks_hetero_18sv9_cosmo_end.Rdata")
cosmo_df <- as.data.frame(matrix(nrow = 11,ncol = 2))
cosmo_df$V1 <- cosmo_names
cosmo_df$V2 <- mean_rel_abun_cosmos
write.csv(cosmo_df, file = "output/euks_hetero_18sv9_cosmo_names.csv")







