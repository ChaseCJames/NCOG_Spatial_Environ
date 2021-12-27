# Archaea GF/F vs Sterivex

metadata <- read.csv("data/NCOG_sample_log_DNA_stvx_meta_2014-2020.csv")

# date by cruise
metadata$Date <- mdy(metadata$Date)

cruise_date <- metadata %>% group_by(Cruise) %>% summarise(mean_date = mean(Date, na.rm = TRUE))

#### Archaea Sterivex ####

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

asv_long$Group <- arch_tax$Group[match(asv_long$hash, arch_tax$Feature.ID)]

asv_long$cruise <- metadata$Cruise[match(asv_long$sample, paste0("X",metadata$Sample.Name))]

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
                        levels = c("Thermoplasmata", "Nitrososphaeria", 
                                   "Altiarchaeia", "Euryarchaeota", 
                                   "Halobacteria", "Nanoarchaeota",
                                   "Other Archaea", "Marine Group II"))

asv_sum$mean_date <- cruise_date$mean_date[match(asv_sum$cruise,cruise_date$Cruise)]

sterivex <- ggplot(asv_sum, aes(x = mean_date, y = rel_abun_cruise, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#58a865","#c65c8a", "#9da13c",
                               "#a361c7","#c28442","#648ace",
                               "#cc5242")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("2014-2016 Sterivex")

p1_a <- ggplot(asv_sum, aes(x = mean_date, y = rel_16, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#58a865","#c65c8a", "#9da13c",
                               "#a361c7","#c28442","#648ace",
                               "#cc5242", "grey60")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("2014-2016 Sterivex Relative Abundance")

mg_ster <- asv_long %>% filter(Group == "Marine Group II")

mg_ster$type <- "Old Seq"
mg_ster$type[which(grepl("_S", mg_ster$sample, fixed = TRUE) == TRUE)] <- "New Seq"

asv_odd <- mg_ster %>% group_by(type,hash) %>%
  summarise(sum_new = sum(abun[which(type == "New Seq")]),
            sum_old = sum(abun[which(type == "Old Seq")])) %>%
  filter(sum_new == 0, sum_old > 50)

asv_short_sum <- asv_long %>% filter(hash %in% asv_odd$hash) %>% group_by(hash,cruise) %>%
  summarise(sum_abun = sum(abun, na.rm = TRUE))

asv_short_sum$samples <- asv_rel_abun$samps[match(asv_short_sum$cruise,asv_rel_abun$cruise)]
asv_short_sum$rel_16 <- asv_short_sum$sum_abun/asv_short_sum$samples
asv_short_sum$mean_date <- cruise_date$mean_date[match(asv_short_sum$cruise,cruise_date$Cruise)]

colors <- c("#755ea5","#53c24c","#a03cb8","#a6ba33","#905ed8","#81b94c","#505dcb",
  "#e6a62f","#5d85ea","#cbb63d","#d670dd","#398e2d","#d143aa","#3ac37d",
  "#e33e89","#5bc9aa","#e43b63","#41c5d5","#c03320","#53a0d9","#e37e21",
  "#aa83e4","#a18b23","#8f47a0","#7a882a","#b62a7c","#73b97b","#a84585",
  "#367c43","#e96da2","#4e6e1b","#d487d2","#69904b","#b53860","#389c84",
  "#e3602d","#4369a7","#e69b45","#9c9cde","#bb7123","#277257","#ed6054",
  "#4d662b","#df98c5","#cda046","#aa6899","#a9b36d","#bc3742","#6a6728",
  "#94486a","#cda96c","#8b3f4c","#8e6b24","#e28395","#947346","#dc766a",
  "#955330","#ee9f81","#a84a22","#ae595a","#eb8d5b","#be7e4c")

order <- asv_odd %>% arrange(sum_old)

asv_short_sum$hash <- factor(asv_short_sum$hash,
                             levels = order$hash)

p2_a <- ggplot(asv_short_sum, aes(x = mean_date, y = rel_16, fill = hash)) +
  geom_area(show.legend = FALSE) +
  scale_fill_manual(values = colors[length(colors):1]) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.07)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("2014-2016 Sterivex Relative Abundance")



##### GF/F #####

load("data/16s_archaea.Rdata")

arch_tax <- six_tax_id %>% filter(Feature.ID %in% colnames(asv_table))

arch_tax$taxon <- arch_tax$silva_Taxon
taxa <- strsplit(arch_tax$taxon, ";")
taxa_sha <- sapply(taxa, "[", 2)
taxa_sha <- sub(paste0(" p__"), "", taxa_sha)   
taxa_deep <- sapply(taxa, "[", 3)
taxa_deep <- sub(paste0(" c__"), "", taxa_deep)  

taxa_sha[which(taxa_deep == "Thermoplasmata")] <- "Thermoplasmata"
taxa_sha[which(taxa_deep == "Nitrososphaeria")] <- "Nitrososphaeria"
taxa_sha[which(taxa_deep == "Halobacteria")] <- "Halobacteria"
taxa_sha[which(taxa_deep == "Altiarchaeia")] <- "Altiarchaeia"

in_groups <- c("Thermoplasmata", "Nitrososphaeria", "Halobacteria",
               "Altiarchaeia", "Nanoarchaeota", "Halobacteria")

taxa_sha[which(!taxa_sha %in% in_groups)] <- "Other Archaea"

arch_tax$Group <- taxa_sha

asv_table <- asv_table %>% as.data.frame()

asv_table$sample <- rownames(asv_table)

asv_long <- asv_table %>% pivot_longer(-sample, names_to = "hash", values_to = "abun")

asv_long$Group <- arch_tax$Group[match(asv_long$hash, arch_tax$Feature.ID)]

asv_long$cruise <- metadata$Cruise[match(asv_long$sample, paste0("X",metadata$Sample.Name))]

asv_sum <- asv_long %>% group_by(Group, cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE)) %>% 
  group_by(cruise) %>% mutate(sum_cruise = sum(sum_abun,na.rm=TRUE))

asv_div <- asv_long %>% filter(abun != 0) %>%
  group_by(Group, cruise) %>% summarise(div = n_distinct(hash)) %>% ungroup() %>%
  complete(Group, cruise, fill = list(div = 0))

asv_sum$rel_abun_cruise <- asv_sum$sum_abun/asv_sum$sum_cruise

asv_sum$mean_date <- cruise_date$mean_date[match(asv_sum$cruise,cruise_date$Cruise)]

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


gf_f <- ggplot(asv_sum, aes(x = mean_date, y = rel_abun_cruise, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#58a865","#c65c8a", "#9da13c",
                               "#a361c7","#c28442","#648ace",
                               "#cc5242"), drop = FALSE) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("2014-2016 GF/F")

out <- sterivex / gf_f + plot_layout(guides = "collect")

agg_png("figures_S/sterviex_gf_f.png", width = 10, height = 6, res = 400, units = "in")
plot(out)
dev.off()


##### Non rareified data GF/F ####

six_s <- six_s[-which(grepl("_S",rownames(six_s)) == TRUE),]

taxa <- strsplit(six_tax_id$silva_Taxon, ";")

taxa_sha <- sapply(taxa, "[", 1)

six_tax_id$arch <- taxa_sha

arch_tax_gf <- six_tax_id %>% filter(arch == "d__Archaea")

asv_table_gf <- six_s

asv_table_gf <- asv_table_gf[,which(!is.na(match(colnames(asv_table_gf),arch_tax_gf$Feature.ID)))]

arch_tax_gf <- six_tax_id %>% filter(Feature.ID %in% colnames(asv_table_gf))

arch_tax_gf$taxon <- arch_tax_gf$silva_Taxon
taxa <- strsplit(arch_tax_gf$taxon, ";")
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

in_groups <- c("Thermoplasmata", "Nitrososphaeria", "Halobacteria",
               "Altiarchaeia", "Nanoarchaeota", "Halobacteria", "Marine Group II")

taxa_sha[which(!taxa_sha %in% in_groups)] <- "Other Archaea"

arch_tax_gf$Group <- taxa_sha

asv_table_gf <- asv_table_gf %>% as.data.frame()

asv_table_gf$sample <- rownames(asv_table_gf)
six_s$sample <- rownames(six_s)

asv_long_gf <- asv_table_gf %>% pivot_longer(-sample, names_to = "hash", values_to = "abun")
asv_long_gf$Group <- arch_tax_gf$Group[match(asv_long_gf$hash, arch_tax_gf$Feature.ID)]
asv_long_gf$cruise <- metadata$Cruise[match(asv_long_gf$sample, paste0("X",metadata$Sample.Name))]

asv_long_gf_six <- six_s %>% pivot_longer(-sample, names_to = "hash", values_to = "abun")
asv_long_gf_six$cruise <- metadata$Cruise[match(asv_long_gf_six$sample, paste0("X",metadata$Sample.Name))]

asv_sum_gf <- asv_long_gf %>% group_by(Group, cruise) %>% summarise(sum_abun = sum(abun, na.rm = TRUE)) %>% 
  group_by(cruise) %>% mutate(sum_cruise = sum(sum_abun,na.rm=TRUE))

total_sum <- asv_long_gf_six %>% group_by(cruise) %>% summarise(sum_total = sum(abun, na.rm = TRUE))

asv_sum_gf$rel_abun_six <- asv_sum_gf$sum_abun/total_sum$sum_total[match(asv_sum_gf$cruise,total_sum$cruise)]

asv_sum_gf$rel_abun_cruise <- asv_sum_gf$sum_abun/asv_sum_gf$sum_cruise

asv_sum_gf$mean_date <- cruise_date$mean_date[match(asv_sum_gf$cruise,cruise_date$Cruise)]

asv_sum_gf$Group <- factor(asv_sum_gf$Group,
                        levels = c("Thermoplasmata", "Nitrososphaeria", 
                                   "Altiarchaeia", "Euryarchaeota", 
                                   "Halobacteria", "Nanoarchaeota",
                                   "Other Archaea", "Marine Group II"))

asv_sum_gf$mean_date <- cruise_date$mean_date[match(asv_sum_gf$cruise,cruise_date$Cruise)]


gf_f_raw <- ggplot(asv_sum_gf, aes(x = mean_date, y = rel_abun_cruise, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#58a865","#c65c8a", "#9da13c",
                               "#a361c7","#c28442","#648ace",
                               "#cc5242"), drop = FALSE) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("2014-2016 GF/F")

p1_b <- ggplot(asv_sum_gf, aes(x = mean_date, y = rel_abun_six, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#58a865","#c65c8a", "#9da13c",
                               "#a361c7","#c28442","#648ace",
                               "#cc5242", "grey60"), drop = FALSE) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("2014-2016 GF/F Relative Abundance")

asv_short_sum_gf <- asv_long_gf %>% filter(hash %in% asv_odd$hash) %>% group_by(hash,cruise) %>%
  summarise(sum_abun = sum(abun, na.rm = TRUE))

asv_short_sum_gf$rel_abun_six <- asv_short_sum_gf$sum_abun/total_sum$sum_total[match(asv_short_sum_gf$cruise,total_sum$cruise)]
asv_short_sum_gf$mean_date <- cruise_date$mean_date[match(asv_short_sum_gf$cruise,cruise_date$Cruise)]

asv_short_sum_gf$hash <- factor(asv_short_sum_gf$hash,
                             levels = order$hash)

p2_b <- ggplot(asv_short_sum_gf, aes(x = mean_date, y = rel_abun_six, fill = hash)) +
  geom_area(show.legend = FALSE) +
  scale_fill_manual(values = colors[length(colors):1]) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.07)) +
  scale_x_date(expand = c(0,0), breaks = "1 year", date_labels = "%Y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance per Cruise") +
  ggtitle("2014-2016 GF/F Relative Abundance")


plot_1 <- p1_a / p1_b + plot_layout(guides = "collect")

plot_2 <- p2_a / p2_b 

agg_png("figures_S/archaea_sterviex_v_gff.png", width = 12, height = 8, units = "in", res = 400)
plot(plot_1)
dev.off()

agg_png("figures_S/mgii_asvs_sterviex_v_gff.png", width = 12, height = 8, units = "in", res = 400)
plot(plot_2)
dev.off()
