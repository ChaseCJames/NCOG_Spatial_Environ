offshore <- ggplot() +
  metR::geom_vector(data = uv_table,aes(x = lon, y = lat, dx = Mean_U, dy = Mean_V), 
                    arrow.angle = 15, arrow.type = "open", arrow.length = unit(0.5, "inches"), 
                    pivot = 0,preserve.dir = TRUE, direction = "ccw",
                    min.mag = 0.05, show.legend = NA) +
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  geom_point(data = som_maps, aes_string(x = "long", y = "lat",
                                         fill = paste0("som_",1)), color = "black", size =8,
                                         stroke = 0.1, shape = 21) +
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0,1)) +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.height = unit(1.4, "cm"), 
        legend.background = element_blank(),
        axis.text = element_text(size = 12, colour = 1),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()) +
  scale_mag(max = 0.1, name = "Speed", max_size = 0.75, guide = "none") +
  labs(x = "Longitude", y = "Latitude") + ggtitle("16s Cyanobacteria Offshore Cluster")


nearshore <- ggplot() +
  metR::geom_vector(data = uv_table,aes(x = lon, y = lat, dx = Mean_U, dy = Mean_V), 
                    arrow.angle = 15, arrow.type = "open", arrow.length = unit(0.5, "inches"), 
                    pivot = 0,preserve.dir = TRUE, direction = "ccw",
                    min.mag = 0.05, show.legend = NA) +
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") +
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  geom_point(data = som_maps, aes_string(x = "long", y = "lat",
                                         fill = paste0("som_",2)), color = "black", size =8,
             stroke = 0.1, shape = 21) +
  scale_fill_gradient(low = "white", high = "darkred", limits = c(0,1)) +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.height = unit(1.4, "cm"), 
        legend.background = element_blank(),
        axis.text = element_text(size = 12, colour = 1),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()) +
  scale_mag(max = 0.1, name = "Speed", max_size = 0.75, guide = "none") +
  labs(x = "Longitude", y = "Latitude") + ggtitle("16s Cyanobacteria Nearshore Cluster")

pdf(file = "figures/current_som_cyano_test.pdf", height = 6, width = 14)
plot_grid(offshore, nearshore, ncol = 2)
dev.off()
