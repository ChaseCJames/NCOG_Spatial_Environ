library(tidyverse)
library(ggmap)


load("data/18s_dino.Rdata")

eight_tax_id$spp <- sapply(strsplit(eight_tax_id$PR2_Taxon, ";"),"[",8)

eight_tax_id <- eight_tax_id %>% filter(spp %in% c("Noctiluca_scintillans","Lingulodinium_polyedra"))
eight_tax_id$sp_id <- make.unique(eight_tax_id$spp,sep = "_")

load("output/dino_18sv9_map.Rdata")

scaled_inputs <- scaled_inputs[,which(!is.na(match(colnames(scaled_inputs), eight_tax_id$Feature.ID)))]
scaled_inputs <- as.data.frame(scaled_inputs)
scaled_inputs$station_name <- gsub("_"," ",substr(rownames(scaled_inputs),9,19))

scaled_long <- scaled_inputs %>% pivot_longer(-station_name, names_to = "ASV", values_to = "rel_abun")

scaled_long$spp <- eight_tax_id$sp_id[match(scaled_long$ASV, eight_tax_id$Feature.ID)]

station_means <- scaled_long %>% group_by(station_name, spp) %>%
  summarise(mean_rel = mean(rel_abun, na.rm = TRUE)) %>% filter(station_name %in% unique(som_maps$Sta_ID))

station_means$lat <- som_maps$lat[match(station_means$station_name, som_maps$Sta_ID)]
station_means$long <- som_maps$long[match(station_means$station_name, som_maps$Sta_ID)]

map <- map_data("world")  

station_means$spp <- as.factor(station_means$spp)
station_means$spp <- factor(station_means$spp, levels = unique(station_means$spp))

out_plot <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  xlab("Longitude") + ylab("Latitude") + geom_point(data = station_means, aes(x = long, y = lat, fill = mean_rel), pch = 21, size = 4) +
  facet_wrap(~spp) + labs(fill = "Mean Rel.\n Abundance")  +
  scale_fill_gradient(low = "white", high = "red") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA, color = "black"),
        strip.background = element_rect(fill = NA, color = "black"))

pdf("figures/HAB_rel_abun.pdf", width = 16, height = 14)
print(out_plot)
dev.off()
