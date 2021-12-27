library(tidyverse)
library(lubridate)
library(patchwork)
library(ragg)

metadata <- read.csv("data/NCOG_sample_log_DNA_stvx_meta_2014-2020.csv")

# date by cruise
metadata$Date <- mdy(metadata$Date)

cruise_date <- metadata %>% group_by(Cruise) %>% summarise(mean_date = mean(Date, na.rm = TRUE))

#### Archaea ####

load("data/16s_archaea_S.Rdata")

arch_tax <- six_tax_id %>% filter(Feature.ID %in% colnames(asv_table))

arch_tax$taxon <- arch_tax$silva_Taxon
taxa <- strsplit(arch_tax$taxon, ";")
taxa_sha <- sapply(taxa, "[", 2)
taxa_sha <- sub(paste0(" p__"), "", taxa_sha)   
taxa_deep <- sapply(taxa, "[", 3)
taxa_deep <- sub(paste0(" c__"), "", taxa_deep)  

taxa_deeper <- sapply(taxa, "[", 4)
taxa_deeper <- sub(paste0(" c__"), "", taxa_deeper) 

taxa_sha[which(taxa_deep == "Thermoplasmata")] <- "Thermoplasmata"
taxa_sha[which(taxa_deep == "Nitrososphaeria")] <- "Nitrososphaeria"
taxa_sha[which(taxa_deep == "Halobacteria")] <- "Halobacteria"
taxa_sha[which(taxa_deep == "Altiarchaeia")] <- "Altiarchaeia"

taxa_sha[which(is.na(taxa_sha))] <- "Other Archaea"

arch_tax$Group <- taxa_sha

asv_table <- asv_table %>% as.data.frame()

asv_table$sample <- rownames(asv_table)

asv_long <- asv_table %>% pivot_longer(-sample, names_to = "hash", values_to = "abun")

asv_long$Group <- arch_tax$Group[match(asv_long$hash, arch_tax$Feature.ID)]

asv_long$cruise <- metadata$Cruise[match(asv_long$sample, paste0("X",metadata$Sample.Name))]

asv_sum <- asv_long %>% group_by(Group, cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE)) %>% 
  group_by(cruise) %>% mutate(sum_cruise = sum(sum_abun,na.rm=TRUE))

asv_rel_abun <- asv_long %>% group_by(cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE),
                                                            samps = n_distinct(sample)) 

asv_rel_abun$samps <- asv_rel_abun$samps*17000

asv_rel_abun$rel_abun_all <- asv_rel_abun$sum_abun/asv_rel_abun$samps

asv_div <- asv_long %>% filter(abun != 0) %>%
  group_by(Group, cruise) %>% summarise(div = n_distinct(hash)) %>% ungroup() %>%
  complete(Group, cruise, fill = list(div = 0))

asv_sum$rel_abun_cruise <- asv_sum$sum_abun/asv_sum$sum_cruise

asv_sum$mean_date <- cruise_date$mean_date[match(asv_sum$cruise,cruise_date$Cruise)]

asv_ratio <- asv_sum[,-c(3,4,6)] %>% filter(Group == "Thermoplasmata" | Group == "Nitrososphaeria") %>% 
  pivot_wider(names_from = Group, values_from = rel_abun_cruise)

asv_ratio$ratio <- asv_ratio$Thermoplasmata-asv_ratio$Nitrososphaeria

asv_ratio$rel_abun <- asv_rel_abun$rel_abun_all[match(asv_ratio$cruise, asv_rel_abun$cruise)]

asv_sum$samples <- asv_rel_abun$samps[match(asv_sum$cruise,asv_rel_abun$cruise)]

asv_sum$rel_16 <- asv_sum$sum_abun/asv_sum$samples

possible <- ggplot(asv_ratio, aes(x = rel_abun, y = ratio, fill = ratio)) +
  geom_point(pch = 21, size = 6) + geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Total Archaeal Relative Abundance\nPer Cruise",
                      y = "Thermoplasmata - Nitrososphaeria",
       fill = "Difference in\nwithin Archaea\nRel. Abun.") +
  scale_fill_gradient2(high = "#58a865", low = "#c65c8a", mid = "white", midpoint = 0) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
  
agg_png("figures_S/archaea_community_shift.png", width = 6, height = 6, res = 400, units = "in")
plot(possible)
dev.off()

asv_sum$Group <- factor(asv_sum$Group,
                        levels = c("Thermoplasmata", "Nitrososphaeria", 
                                   "Altiarchaeia", "Euryarchaeota", 
                                   "Halobacteria", "Nanoarchaeota",
                                   "Other Archaea"))

asv_sum$mean_date <- cruise_date$mean_date[match(asv_sum$cruise,cruise_date$Cruise)]
asv_div$mean_date <- cruise_date$mean_date[match(asv_div$cruise,cruise_date$Cruise)]

asv_div$Group <- factor(asv_div$Group,
                        levels = c("Thermoplasmata", "Nitrososphaeria", 
                                   "Altiarchaeia", "Euryarchaeota", 
                                   "Halobacteria", "Nanoarchaeota",
                                   "Other Archaea"))


a <- ggplot(asv_sum, aes(x = mean_date, y = rel_abun_cruise, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#58a865","#c65c8a", "#9da13c",
                               "#a361c7","#c28442","#648ace",
                               "#cc5242")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("Archaea")

 ggplot(asv_sum, aes(x = mean_date, y = rel_16, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#58a865","#c65c8a", "#9da13c",
                               "#a361c7","#c28442","#648ace",
                               "#cc5242")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("Archaea")


b <- ggplot(asv_div, aes(x = mean_date, y = div, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#58a865","#c65c8a", "#9da13c",
                               "#a361c7","#c28442","#648ace",
                               "#cc5242")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "Date", y = "Richness per Cruise") 

out <- a/b + plot_layout(guides = "collect")

agg_png(filename = "figures_S/archaea_changes.png", width = 10, height = 8, units = "in", res =400)
plot(out)
dev.off()




#### Cyanobacteria ####

load("data/16s_cyanos_S.Rdata")

cyano_tax <- six_tax_id %>% filter(Feature.ID %in% colnames(asv_table))

cyano_tax$taxon <- cyano_tax$silva_Taxon
taxa <- strsplit(cyano_tax$taxon, ";")
taxa_sha <- sapply(taxa, "[", 6)

taxa_sha[which(substr(taxa_sha,5,19) == "Prochlorococcus")] <- "Prochlorococcus"
taxa_sha[which(substr(taxa_sha,5,17) == "Synechococcus")] <- "Synechococcus"


taxa_sha[which(taxa_sha != "Prochlorococcus" & taxa_sha != "Synechococcus")] <- "Other Cyanobacteria"

taxa_sha[which(is.na(taxa_sha))] <- "Other Cyanobacteria"

cyano_tax$Group <- taxa_sha

asv_table <- asv_table %>% as.data.frame()

asv_table$sample <- rownames(asv_table)

asv_long <- asv_table %>% pivot_longer(-sample, names_to = "hash", values_to = "abun")

asv_long$Group <- cyano_tax$Group[match(asv_long$hash, cyano_tax$Feature.ID)]

asv_long$cruise <- metadata$Cruise[match(asv_long$sample, paste0("X",metadata$Sample.Name))]

asv_sum <- asv_long %>% group_by(Group, cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE)) %>% 
  group_by(cruise) %>% mutate(sum_cruise = sum(sum_abun,na.rm=TRUE))

asv_div <- asv_long %>% filter(abun != 0) %>%
  group_by(Group, cruise) %>% summarise(div = n_distinct(hash)) %>% ungroup() %>%
  complete(Group, cruise, fill = list(div = 0))

asv_sum$rel_abun_cruise <- asv_sum$sum_abun/asv_sum$sum_cruise

asv_sum$mean_date <- cruise_date$mean_date[match(asv_sum$cruise,cruise_date$Cruise)]

asv_sum$Group <- factor(asv_sum$Group,
                        levels = c("Prochlorococcus", "Synechococcus", 
                                   "Other Cyanobacteria"))


asv_div$mean_date <- cruise_date$mean_date[match(asv_div$cruise,cruise_date$Cruise)]

asv_div$Group <- factor(asv_div$Group,
                        levels = c("Prochlorococcus", "Synechococcus", 

                                   "Other Cyanobacteria"))


a <- ggplot(asv_sum, aes(x = mean_date, y = rel_abun_cruise, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#cb6751",
                               "#7aa457",
                               "#9e6ebd")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("Cyanobacteria")


b <- ggplot(asv_div, aes(x = mean_date, y = div, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#cb6751",
                               "#7aa457",
                               "#9e6ebd")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "Date", y = "Richness per Cruise") 

out <- a/b + plot_layout(guides = "collect")

agg_png(filename = "figures_S/cyano_changes.png", width = 10, height = 8, units = "in", res =400)
plot(out)
dev.off()

#### Diatoms ####

load("data/18s_diatom_S.Rdata")

diatom_tax <- eight_tax_id %>% filter(Feature.ID %in% colnames(asv_table))

diatom_tax$taxon <- diatom_tax$pr2_Taxon
taxa <- strsplit(diatom_tax$taxon, ";")
taxa_sha <- sapply(taxa, "[", 6)
taxa_deep <- sapply(taxa, "[", 7)


taxa_sha[which(is.na(taxa_sha))] <- "Other Diatoms"

diatom_tax$Group <- taxa_sha

asv_table <- asv_table %>% as.data.frame()

asv_table$sample <- rownames(asv_table)

asv_long <- asv_table %>% pivot_longer(-sample, names_to = "hash", values_to = "abun")

asv_long$Group <- diatom_tax$Group[match(asv_long$hash, diatom_tax$Feature.ID)]

asv_long$cruise <- metadata$Cruise[match(asv_long$sample, paste0("X",metadata$Sample.Name))]

asv_sum <- asv_long %>% group_by(Group, cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE)) %>% 
  group_by(cruise) %>% mutate(sum_cruise = sum(sum_abun,na.rm=TRUE))

asv_div <- asv_long %>% filter(abun != 0) %>%
  group_by(Group, cruise) %>% summarise(div = n_distinct(hash)) %>% ungroup() %>%
  complete(Group, cruise, fill = list(div = 0))

asv_sum$rel_abun_cruise <- asv_sum$sum_abun/asv_sum$sum_cruise

asv_sum$mean_date <- cruise_date$mean_date[match(asv_sum$cruise,cruise_date$Cruise)]

asv_div$mean_date <- cruise_date$mean_date[match(asv_div$cruise,cruise_date$Cruise)]

order <- asv_sum %>% group_by(Group) %>% summarise(mean = mean(rel_abun_cruise)) %>%
  arrange(mean)

asv_sum$Group <- factor(asv_sum$Group, levels = order$Group)
asv_div$Group <- factor(asv_div$Group, levels = order$Group)

a <- ggplot(asv_sum, aes(x = mean_date, y = rel_abun_cruise, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#c75a93",
                               "#5ba966",
                               "#8176cc",
                               "#ac973e",
                               "#cc5f42")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("Diatoms")


b <- ggplot(asv_div, aes(x = mean_date, y = div, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#c75a93",
                               "#5ba966",
                               "#8176cc",
                               "#ac973e",
                               "#cc5f42")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "Date", y = "Richness per Cruise") 

out <- a/b + plot_layout(guides = "collect")

agg_png(filename = "figures_S/diatom_changes.png", width = 10, height = 8, units = "in", res =400)
plot(out)
dev.off()

#### Syndiniales #####

load("data/18s_syndin_S.Rdata")

syndin_tax <- eight_tax_id %>% filter(Feature.ID %in% colnames(asv_table))

syndin_tax$taxon <- syndin_tax$pr2_Taxon
taxa <- strsplit(syndin_tax$taxon, ";")
taxa_sha <- sapply(taxa, "[", 5)

taxa_sha[which(taxa_sha == "Syndiniales_X")] <- "Other Syndiniales"

taxa_sha[which(is.na(taxa_sha))] <- "Other Syndiniales"

syndin_tax$Group <- taxa_sha

asv_table <- asv_table %>% as.data.frame()

asv_table$sample <- rownames(asv_table)

asv_long <- asv_table %>% pivot_longer(-sample, names_to = "hash", values_to = "abun")

asv_long$Group <- syndin_tax$Group[match(asv_long$hash, syndin_tax$Feature.ID)]

asv_long$cruise <- metadata$Cruise[match(asv_long$sample, paste0("X",metadata$Sample.Name))]

asv_sum <- asv_long %>% group_by(Group, cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE)) %>% 
  group_by(cruise) %>% mutate(sum_cruise = sum(sum_abun,na.rm=TRUE))

asv_div <- asv_long %>% filter(abun != 0) %>%
  group_by(Group, cruise) %>% summarise(div = n_distinct(hash)) %>% ungroup() %>%
  complete(Group, cruise, fill = list(div = 0))

asv_sum$rel_abun_cruise <- asv_sum$sum_abun/asv_sum$sum_cruise

asv_sum$mean_date <- cruise_date$mean_date[match(asv_sum$cruise,cruise_date$Cruise)]

asv_div$mean_date <- cruise_date$mean_date[match(asv_div$cruise,cruise_date$Cruise)]

order <- asv_sum %>% group_by(Group) %>% summarise(mean = mean(rel_abun_cruise)) %>%
  arrange(mean)

asv_sum$Group <- factor(asv_sum$Group, levels = order$Group)
asv_div$Group <- factor(asv_div$Group, levels = order$Group)

a <- ggplot(asv_sum, aes(x = mean_date, y = rel_abun_cruise, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#aa973d",
                               "#9e65c4",
                               "#5ba85f",
                               "#c95680",
                               "#6295cd",
                               "#cc6441")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("Syndiniales")


b <- ggplot(asv_div, aes(x = mean_date, y = div, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#aa973d",
                               "#9e65c4",
                               "#5ba85f",
                               "#c95680",
                               "#6295cd",
                               "#cc6441")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "Date", y = "Richness per Cruise") 

out <- a/b + plot_layout(guides = "collect")

agg_png(filename = "figures_S/syndin_changes.png", width = 10, height = 8, units = "in", res =400)
plot(out)
dev.off()
