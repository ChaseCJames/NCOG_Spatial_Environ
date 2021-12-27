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

taxa_deeper <- sapply(taxa, "[", 6)
taxa_deeper <- sub(paste0(" g__"), "", taxa_deeper) 

taxa_sha[which(taxa_deep == "Thermoplasmata")] <- "Thermoplasmata"
taxa_sha[which(taxa_deep == "Nitrososphaeria")] <- "Nitrososphaeria"
taxa_sha[which(taxa_deep == "Halobacteria")] <- "Halobacteria"
taxa_sha[which(taxa_deep == "Altiarchaeia")] <- "Altiarchaeia"
taxa_sha[which(taxa_deeper == "Marine_Group_II")] <- "Marine Group II"

taxa_sha[which(is.na(taxa_sha))] <- "Other Archaea"

arch_tax$Group <- taxa_sha

asv_table <- asv_table %>% as.data.frame()

asv_table$sample <- rownames(asv_table)

asv_long <- asv_table %>% pivot_longer(-sample, names_to = "hash", values_to = "abun")

asv_long$cruise <- metadata$Cruise[match(asv_long$sample, paste0("X",metadata$Sample.Name))]

asv_long$Group <- arch_tax$Group[match(asv_long$hash, arch_tax$Feature.ID)]

asv_sum <- asv_long %>% group_by(Group, cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE)) %>% 
  group_by(cruise) %>% mutate(sum_cruise = sum(sum_abun,na.rm=TRUE))

asv_rel_abun <- asv_long %>% group_by(cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE),
                                                            samps = n_distinct(sample)) 

asv_rel_abun$samps <- asv_rel_abun$samps*17000

asv_rel_abun$rel_abun_all <- asv_rel_abun$sum_abun/asv_rel_abun$samps

asv_sum$rel_abun_cruise <- asv_sum$sum_abun/asv_sum$sum_cruise

asv_sum$mean_date <- cruise_date$mean_date[match(asv_sum$cruise,cruise_date$Cruise)]

asv_sum$samples <- asv_rel_abun$samps[match(asv_sum$cruise,asv_rel_abun$cruise)]

asv_sum$rel_16 <- asv_sum$sum_abun/asv_sum$samples

asv_sum$Group <- factor(asv_sum$Group,
                        levels = c("Marine Group II",
                                   "Thermoplasmata", "Nitrososphaeria", 
                                   "Altiarchaeia", "Euryarchaeota", 
                                   "Halobacteria", "Nanoarchaeota",
                                   "Other Archaea"))

asv_sum$mean_date <- cruise_date$mean_date[match(asv_sum$cruise,cruise_date$Cruise)]

archaea <- ggplot(asv_sum, aes(x = mean_date, y = rel_abun_cruise, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#688dcd","#66b046","#8462cc","#a6953f",
                               "#c361aa","#55a47b","#cc566a","#cd6c39")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("Archaea")

 arch_rel <- ggplot(asv_sum, aes(x = mean_date, y = rel_16, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#688dcd","#66b046","#8462cc","#a6953f",
                               "#c361aa","#55a47b","#cc566a","#cd6c39")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("Archaea")


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

asv_rel_abun <- asv_long %>% group_by(cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE),
                                                            samps = n_distinct(sample)) 

asv_rel_abun$samps <- asv_rel_abun$samps*17000

asv_rel_abun$rel_abun_all <- asv_rel_abun$sum_abun/asv_rel_abun$samps

asv_sum$samples <- asv_rel_abun$samps[match(asv_sum$cruise,asv_rel_abun$cruise)]

asv_sum$rel_16 <- asv_sum$sum_abun/asv_sum$samples

cyano <- ggplot(asv_sum, aes(x = mean_date, y = rel_abun_cruise, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#cb6751",
                               "#7aa457",
                               "#9e6ebd")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("Cyanobacteria")

cyano_rel <- ggplot(asv_sum, aes(x = mean_date, y = rel_16, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#cb6751",
                               "#7aa457",
                               "#9e6ebd")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("Cyanobacteria")


##### Bacteria #### 

load("data/16s_bacteria_m_euks_S.Rdata")

bact_tax <- six_tax_id %>% filter(Feature.ID %in% colnames(asv_table))

bact_tax$taxon <- bact_tax$silva_Taxon
taxa <- strsplit(bact_tax$taxon, ";")
taxa_sha <- sapply(taxa, "[", 3)
taxa_sha <- sub(paste0(" c__"), "", taxa_sha)   
taxa_deep <- sapply(taxa, "[", 4)
taxa_deep <- sub(paste0(" o__"), "", taxa_deep)  

taxa_sha[which(taxa_deep == "Flavobacteriales")] <- "Flavobacteriales"
taxa_sha[which(taxa_deep == "Rhodobacterales")] <- "Rhodobacterales"
taxa_sha[which(taxa_deep == "SAR11_clade")] <- "SAR 11 clade"

in_groups <- c("SAR 11 clade", "Flavobacteriales", "Gammaproteobacteria",
               "Alphaproteobacteria", "Rhodobacterales", "Acidimicrobiia",
               "Verrucomicrobiae")

taxa_sha[which(!taxa_sha %in% in_groups)] <- "Other Bacteria"

bact_tax$Group <- taxa_sha

asv_table <- asv_table %>% as.data.frame()

asv_table$sample <- rownames(asv_table)

asv_long <- asv_table %>% pivot_longer(-sample, names_to = "hash", values_to = "abun")

asv_long$Group <- bact_tax$Group[match(asv_long$hash, bact_tax$Feature.ID)]

asv_long$cruise <- metadata$Cruise[match(asv_long$sample, paste0("X",metadata$Sample.Name))]

asv_sum <- asv_long %>% group_by(Group, cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE)) %>% 
  group_by(cruise) %>% mutate(sum_cruise = sum(sum_abun,na.rm=TRUE))

asv_sum$rel_abun_cruise <- asv_sum$sum_abun/asv_sum$sum_cruise

asv_sum$mean_date <- cruise_date$mean_date[match(asv_sum$cruise,cruise_date$Cruise)]

asv_sum$Group <- factor(asv_sum$Group,
                        levels = c("SAR 11 clade", "Flavobacteriales", "Gammaproteobacteria",
                                   "Alphaproteobacteria", "Rhodobacterales", "Acidimicrobiia",
                                   "Verrucomicrobiae", 
                                   "Other Bacteria"))

asv_rel_abun <- asv_long %>% group_by(cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE),
                                                            samps = n_distinct(sample)) 

asv_rel_abun$samps <- asv_rel_abun$samps*17000

asv_rel_abun$rel_abun_all <- asv_rel_abun$sum_abun/asv_rel_abun$samps

asv_sum$samples <- asv_rel_abun$samps[match(asv_sum$cruise,asv_rel_abun$cruise)]

asv_sum$rel_16 <- asv_sum$sum_abun/asv_sum$samples


bact <- ggplot(asv_sum, aes(x = mean_date, y = rel_abun_cruise, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#c55d93",
                               "#64ac48",
                               "#b05cc6",
                               "#a39440",
                               "#6e7fcc",
                               "#cf723a",
                               "#4aac8d",
                               "#cb5358")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("Bacteria")

bact_rel <- ggplot(asv_sum, aes(x = mean_date, y = rel_16, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#c55d93",
                               "#64ac48",
                               "#b05cc6",
                               "#a39440",
                               "#6e7fcc",
                               "#cf723a",
                               "#4aac8d",
                               "#cb5358")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("Bacteria")

##### Photosynthetic Eukaryotes ######

load("data/18s_autotrophic_euks_S.Rdata")

photo_tax <- eight_tax_id %>% filter(Feature.ID %in% colnames(asv_table))

photo_tax$taxon <- photo_tax$pr2_Taxon
taxa <- strsplit(photo_tax$taxon, ";")
taxa_sha <- sapply(taxa, "[", 3)
taxa_deep <- sapply(taxa, "[", 4)

taxa_sha[which(taxa_deep == "Mamiellophyceae")] <- "Mamiellophyceae"
taxa_sha[which(taxa_deep == "Bacillariophyta")] <- "Bacillariophyta"

photo_tax$Group <- taxa_sha

asv_table <- asv_table %>% as.data.frame()

asv_table$sample <- rownames(asv_table)

asv_long <- asv_table %>% pivot_longer(-sample, names_to = "hash", values_to = "abun")

asv_long$Group <- photo_tax$Group[match(asv_long$hash, photo_tax$Feature.ID)]

asv_long$cruise <- metadata$Cruise[match(asv_long$sample, paste0("X",metadata$Sample.Name))]

asv_sum <- asv_long %>% group_by(Group, cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE)) %>% 
  group_by(cruise) %>% mutate(sum_cruise = sum(sum_abun,na.rm=TRUE))

asv_sum$rel_abun_cruise <- asv_sum$sum_abun/asv_sum$sum_cruise

asv_sum$mean_date <- cruise_date$mean_date[match(asv_sum$cruise,cruise_date$Cruise)]

asv_order <- asv_sum %>% group_by(Group) %>% 
  summarise(mean = mean(rel_abun_cruise,na.rm = TRUE)) %>%
  arrange(mean)

asv_sum$Group <- factor(asv_sum$Group,
                        levels = asv_order$Group)

asv_rel_abun <- asv_long %>% group_by(cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE),
                                                            samps = n_distinct(sample)) 

asv_rel_abun$samps <- asv_rel_abun$samps*17000

asv_rel_abun$rel_abun_all <- asv_rel_abun$sum_abun/asv_rel_abun$samps

asv_sum$samples <- asv_rel_abun$samps[match(asv_sum$cruise,asv_rel_abun$cruise)]

asv_sum$rel_18 <- asv_sum$sum_abun/asv_sum$samples

photo <- ggplot(asv_sum, aes(x = mean_date, y = rel_abun_cruise, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#c55d93",
                               "#64ac48",
                               "#b05cc6",
                               "#a39440",
                               "#6e7fcc",
                               "#cf723a",
                               "#4aac8d",
                               "#cb5358")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("Photosynthetic Eukaryotic Protists")


photo_rel <- ggplot(asv_sum, aes(x = mean_date, y = rel_18, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#c55d93",
                               "#64ac48",
                               "#b05cc6",
                               "#a39440",
                               "#6e7fcc",
                               "#cf723a",
                               "#4aac8d",
                               "#cb5358")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("Photosynthetic Eukaryotic Protists")

#### Heterotrophic Eukaryotic Protists ####

load("data/18s_heterotrophic_euks_S.Rdata")

hetero_tax <- eight_tax_id %>% filter(Feature.ID %in% colnames(asv_table))

hetero_tax$taxon <- hetero_tax$pr2_Taxon
taxa <- strsplit(hetero_tax$taxon, ";")
taxa_sha <- sapply(taxa, "[", 4)

in_groups <- c("Syndiniales", "Polycystinea", "Spirotrichea",
               "Ascomycota", "RAD-B", "MAST-3", "Telonemia_X",
               "Acantharea", "MAST-1", "MAST-4")

taxa_sha[which(!taxa_sha %in% in_groups)] <- "Other Heterotrophic Eukaryotes"

taxa_sha <- gsub("_X","",taxa_sha)

hetero_tax$Group <- taxa_sha

asv_table <- asv_table %>% as.data.frame()

asv_table$sample <- rownames(asv_table)

asv_long <- asv_table %>% pivot_longer(-sample, names_to = "hash", values_to = "abun")

asv_long$Group <- hetero_tax$Group[match(asv_long$hash, hetero_tax$Feature.ID)]

asv_long$cruise <- metadata$Cruise[match(asv_long$sample, paste0("X",metadata$Sample.Name))]

asv_sum <- asv_long %>% group_by(Group, cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE)) %>% 
  group_by(cruise) %>% mutate(sum_cruise = sum(sum_abun,na.rm=TRUE))

asv_sum$rel_abun_cruise <- asv_sum$sum_abun/asv_sum$sum_cruise

asv_sum$mean_date <- cruise_date$mean_date[match(asv_sum$cruise,cruise_date$Cruise)]

asv_order <- asv_sum %>% group_by(Group) %>% 
  summarise(mean = mean(rel_abun_cruise,na.rm = TRUE)) %>%
  arrange(mean)

asv_sum$Group <- factor(asv_sum$Group,
                        levels = asv_order$Group)

asv_rel_abun <- asv_long %>% group_by(cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE),
                                                            samps = n_distinct(sample)) 

asv_rel_abun$samps <- asv_rel_abun$samps*17000

asv_rel_abun$rel_abun_all <- asv_rel_abun$sum_abun/asv_rel_abun$samps

asv_sum$samples <- asv_rel_abun$samps[match(asv_sum$cruise,asv_rel_abun$cruise)]

asv_sum$rel_18 <- asv_sum$sum_abun/asv_sum$samples

hetero <- ggplot(asv_sum, aes(x = mean_date, y = rel_abun_cruise, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#c75d9c",
                               "#65b74f",
                               "#8f62ca",
                               "#b7b43e",
                               "#6b8bcd",
                               "#e09042",
                               "#51bfa5",
                               "#cf5054",
                               "#428248",
                               "#ae673d",
                               "#86863d")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("Heterotrophic Eukaryotic Protists")

hetero_rel <- ggplot(asv_sum, aes(x = mean_date, y = rel_18, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#c75d9c",
                               "#65b74f",
                               "#8f62ca",
                               "#b7b43e",
                               "#6b8bcd",
                               "#e09042",
                               "#51bfa5",
                               "#cf5054",
                               "#428248",
                               "#ae673d",
                               "#86863d")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("Heterotrophic Eukaryotic Protists")



out <- archaea + bact + cyano + 
  photo + hetero + plot_layout(nrow = 1) + plot_annotation(tag_levels = "a")

pdf("figures_S/suppl_fig_XX.pdf", width = 30, height = 5) 
plot(out)
dev.off()

out2 <- arch_rel + bact_rel + cyano_rel + 
  photo_rel + hetero_rel + plot_layout(nrow = 1) + plot_annotation(tag_levels = "a")

pdf("figures_S/rel_abun_plots.pdf", width = 30, height = 5) 
plot(out2)
dev.off()