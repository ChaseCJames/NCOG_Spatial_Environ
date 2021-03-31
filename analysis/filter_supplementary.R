library(tidyverse)
library(patchwork)
library(ecodist)
library(vegan)
library(viridis)


###### 16S ######

six <- read.csv("data/16_asv_count_tax_final_19.csv", stringsAsFactors = FALSE)

six_ASV_ids <- six[,c(1,1033:1036)]

six_ASV_ids$Silva_Taxon <- gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__|D_7__|D_8__", "", six_ASV_ids$Silva_Taxon)
split_six_taxa <- separate(six_ASV_ids, Silva_Taxon, sep = ";", into = c("A","B","C", "D", "E", "F", "G", "H", "I"))

samps <- c("201611_093.3_026.7_24", "201611_093.3_026.7_24_T", "201611_093.3_045.0_33", "201611_093.3_045.0_33_T",
"201611_093.3_070.0_40", "201611_093.3_070.0_40_T", "201611_086.7_070.0_10","201611_086.7_070.0_10_T",
"201611_076.7_070.0_8","201611_076.7_070.0_8_T")

cruise <- colnames(six)[which(substr(colnames(six),2,7) == "201611")]

samps <- paste0("X", samps)

six <- six[,match(cruise, colnames(six))]

samples <- colnames(six)

six <- as.data.frame(t(six))

rownames(six) <- samples
colnames(six) <- six_ASV_ids$Feature.ID
six$sample <- samples

six_long <- six %>% pivot_longer(-sample)

six_long$occur <- 0
six_long$occur[which(six_long$value > 0)] <- 1

###### prok analysis #####

prok <- split_six_taxa$Feature.ID[which(split_six_taxa$A == "Bacteria" | split_six_taxa$A == "Archaea" & split_six_taxa$D != "Chloroplast")]
notprok <- split_six_taxa$Feature.ID[which(split_six_taxa$A != "Bacteria" & split_six_taxa$A != "Archaea" | split_six_taxa$D == "Chloroplast")]

length(prok)/ncol(six)
six$sample <- NULL

prok_sequences <- sum(six[,which(colnames(six) %in% prok)])
tot_sequences <- sum(six)

prok_sequences/tot_sequences

six_summary <- six_long %>%
  filter(name %in% prok) %>%
  group_by(sample) %>%
  summarise(total_reads = sum(value, na.rm = TRUE), richness = sum(occur, na.rm = TRUE))

six_summary$type <- "GF/F (201611 Cruise)"

six_summary$type[match(samps, six_summary$sample)] <- rep(c("GF/F", "Swinnex"), times = 5)

six_summary$sample_num <- NA

six_summary$sample_num[match(samps, six_summary$sample)] <- rep(1:5, each = 2)

wide_six <- six_summary %>% filter(!is.na(sample_num)) %>%
  pivot_wider(-sample, names_from = type, values_from = c(total_reads, richness))

total_reads_six <- t.test(wide_six$`total_reads_GF/F`, wide_six$total_reads_Swinnex, paired = TRUE)
richness_six <- t.test(wide_six$`richness_GF/F`, wide_six$richness_Swinnex, paired = TRUE)

##### prop 16 #####

six_long$group <- "prok"
six_long$group[which(!is.na(match(six_long$name, notprok)))] <- "notprok"

six_summary1 <- six_long %>%
  group_by(sample, group) %>%
  summarise(total_reads = sum(value, na.rm = TRUE), richness = sum(occur, na.rm = TRUE))

six_summary1 <- six_summary1 %>% pivot_wider(names_from = group, values_from = c(total_reads, richness))

six_summary1$type <- "GF/F (201611 Cruise)"

six_summary1$type[match(samps, six_summary1$sample)] <- rep(c("GF/F", "Swinnex"), times = 5)

six_summary1$sample_num <- NA

six_summary1$sample_num[match(samps, six_summary1$sample)] <- rep(1:5, each = 2)

six_summary1$prop_prok_reads <- six_summary1$total_reads_prok / (six_summary1$total_reads_notprok + six_summary1$total_reads_prok)
six_summary1$prop_prok_rich <- six_summary1$richness_prok / (six_summary1$richness_notprok + six_summary1$richness_prok)
  
wide_six1 <- six_summary1[,c(1,6:9)] %>% filter(!is.na(sample_num)) %>%
  pivot_wider(-sample, names_from = type, values_from = c(prop_prok_reads, prop_prok_rich))

prop_reads_six <- t.test(wide_six1$`prop_prok_reads_GF/F`, wide_six1$prop_prok_reads_Swinnex, paired = TRUE)
prop_richness_six <- t.test(wide_six1$`prop_prok_rich_GF/F`, wide_six1$prop_prok_rich_Swinnex, paired = TRUE)


####### 18S ##########

eight <- read.csv("data/18Sv9_asv_count_tax_final_19.csv", stringsAsFactors = FALSE)

eight_ASV_ids <- eight[,c(1,1030:1033)]

eight_ASV_ids$Silva_Taxon <- gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__|D_7__|D_8__", "", eight_ASV_ids$Silva_Taxon)
eight_ASV_ids$PR2_Taxon <- gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__|D_7__|D_8__", "", eight_ASV_ids$PR2_Taxon)

eight_ASV_ids$tax <- NA
eight_ASV_ids$tax[which(eight_ASV_ids$PR2_Confidence >= eight_ASV_ids$Silva_Confidence)] <- eight_ASV_ids$PR2_Taxon[which(eight_ASV_ids$PR2_Confidence >= eight_ASV_ids$Silva_Confidence)]
eight_ASV_ids$tax[which(eight_ASV_ids$PR2_Confidence < eight_ASV_ids$Silva_Confidence)] <- eight_ASV_ids$Silva_Taxon[which(eight_ASV_ids$PR2_Confidence < eight_ASV_ids$Silva_Confidence)]

split_eight_taxa <- separate(eight_ASV_ids, tax, sep = ";", into = c("A","B","C", "D", "E", "F", "G", "H", "I"))

samps <- c("201611_093.3_026.7_24", "201611_093.3_026.7_24_T", "201611_093.3_045.0_33", "201611_093.3_045.0_33_T",
           "201611_093.3_070.0_40", "201611_093.3_070.0_40_T", "201611_086.7_070.0_10","201611_086.7_070.0_10_T",
           "201611_076.7_070.0_8","201611_076.7_070.0_8_T")

cruise <- colnames(eight)[which(substr(colnames(eight),2,7) == "201611")]

samps <- paste0("X", samps)

eight <- eight[,match(cruise, colnames(eight))]

samples <- colnames(eight)

eight <- as.data.frame(t(eight))

rownames(eight) <- samples
colnames(eight) <- eight_ASV_ids$Feature.ID
eight$sample <- samples

eight_long <- eight %>% pivot_longer(-sample)

eight_long$occur <- 0
eight_long$occur[which(eight_long$value > 0)] <- 1

euk <- split_eight_taxa$Feature.ID[which(split_eight_taxa$A == "Eukaryota" )]
noteuk <- split_eight_taxa$Feature.ID[which(split_eight_taxa$A != "Eukaryota" )]

length(euk)/ncol(eight)
eight$sample <- NULL

euk_sequences <- sum(eight[,which(colnames(eight) %in% euk)])
tot_sequences <- sum(eight)

euk_sequences/tot_sequences


eight_summary <- eight_long %>%
  group_by(sample) %>%
  summarise(total_reads = sum(value, na.rm = TRUE), richness = sum(occur, na.rm = TRUE))

eight_summary$type <- "GF/F (201611 Cruise)"

eight_summary$type[match(samps, eight_summary$sample)] <- rep(c("GF/F", "Swinnex"), times = 5)

eight_summary$sample_num <- NA

eight_summary$sample_num[match(samps, eight_summary$sample)] <- rep(1:5, each = 2)

wide_eight <- eight_summary %>% filter(!is.na(sample_num)) %>%
  pivot_wider(-sample, names_from = type, values_from = c(total_reads, richness))

total_reads_eight <- t.test(wide_eight$`total_reads_GF/F`, wide_eight$total_reads_Swinnex, paired = TRUE)
richness_eight <- t.test(wide_eight$`total_reads_GF/F`, wide_eight$richness_Swinnex, paired = TRUE)

##### prop 18 #####

eight_long$group <- "euk"
eight_long$group[which(!is.na(match(eight_long$name, noteuk)))] <- "noteuk"

eight_summary1 <- eight_long %>%
  group_by(sample, group) %>%
  summarise(total_reads = sum(value, na.rm = TRUE), richness = sum(occur, na.rm = TRUE))

eight_summary1 <- eight_summary1 %>% pivot_wider(names_from = group, values_from = c(total_reads, richness))

eight_summary1$type <- "GF/F (201611 Cruise)"

eight_summary1$type[match(samps, eight_summary1$sample)] <- rep(c("GF/F", "Swinnex"), times = 5)

eight_summary1$sample_num <- NA

eight_summary1$sample_num[match(samps, eight_summary1$sample)] <- rep(1:5, each = 2)

eight_summary1$prop_euk_reads <- eight_summary1$total_reads_euk / (eight_summary1$total_reads_noteuk + eight_summary1$total_reads_euk)
eight_summary1$prop_euk_rich <- eight_summary1$richness_euk / (eight_summary1$richness_noteuk + eight_summary1$richness_euk)

wide_eight1 <- eight_summary1[,c(1,6:9)] %>% filter(!is.na(sample_num)) %>%
  pivot_wider(-sample, names_from = type, values_from = c(prop_euk_reads, prop_euk_rich))

prop_reads_eight <- t.test(wide_eight1$`prop_euk_reads_GF/F`, wide_eight1$prop_euk_reads_Swinnex, paired = TRUE)
prop_richness_eight <- t.test(wide_eight1$`prop_euk_rich_GF/F`, wide_eight1$prop_euk_rich_Swinnex, paired = TRUE)

##### combine for plot #####

six_summary$region <- "16SV4-V5"
eight_summary$region <- "18SV9"

total_summary <- bind_rows(six_summary, eight_summary)

total_summary$type <- as.factor(total_summary$type)
total_summary$type <- factor(total_summary$type, 
                             levels = c("GF/F", "Swinnex", 
                                        "GF/F (201611 Cruise)"))

richness <- ggplot(total_summary, aes(x = region, y = richness, fill = type)) +
  geom_boxplot() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "", y = "Richness", fill = "Filter Type") +
  annotate("text", x = 0.9, y = 900, label = paste0("p-value = ", round(richness_six$p.value,3))) +
  annotate("text", x = 1.9, y = 500, label = paste0("p-value = ", round(richness_eight$p.value,3)))

total_reads <- ggplot(total_summary, aes(x = region, y = total_reads, fill = type)) +
  geom_boxplot() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "", y = "# of Amplicons", fill = "Filter Type") +
  annotate("text", x = 0.9, y = 130000, label = paste0("p-value = ", round(total_reads_six$p.value,3))) +
  annotate("text", x = 1.9, y = 27000, label = paste0("p-value = ", round(total_reads_eight$p.value,3)))

out <- total_reads + richness + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect")

pdf(file = "figures/figure_outline/supp_fig_XX.pdf", width = 8.5, height = 5)
print(out)
dev.off()

####### Sar & prochloro #####

# SAR 11

sar_tax <- split_six_taxa$Feature.ID[which(split_six_taxa$D == "SAR11 clade")]

sar_summary <- six_long %>%
  filter(name %in% sar_tax) %>%
  group_by(sample) %>%
  summarise(total_reads = sum(value, na.rm = TRUE), richness = sum(occur, na.rm = TRUE))

sar_summary$type <- "GF/F (201611 Cruise)"

sar_summary$type[match(samps, sar_summary$sample)] <- rep(c("GF/F", "Swinnex"), times = 5)

sar_summary$sample_num <- NA

sar_summary$sample_num[match(samps, sar_summary$sample)] <- rep(1:5, each = 2)

sar_summary$group <- "Sar 11 Clade"

wide_sar <- sar_summary %>% filter(!is.na(sample_num)) %>%
  pivot_wider(-sample, names_from = type, values_from = c(total_reads, richness))

total_reads_sar <- t.test(wide_sar$`total_reads_GF/F`, wide_sar$total_reads_Swinnex, paired = TRUE)
richness_sar <- t.test(wide_sar$`richness_GF/F`, wide_sar$richness_Swinnex, paired = TRUE)

# pro

pro_tax <- split_six_taxa$Feature.ID[which(substr(split_six_taxa$F, 1, 9) == "Prochloro")]

pro_summary <- six_long %>%
  filter(name %in% pro_tax) %>%
  group_by(sample) %>%
  summarise(total_reads = sum(value, na.rm = TRUE), richness = sum(occur, na.rm = TRUE))

pro_summary$type <- "GF/F (201611 Cruise)"

pro_summary$type[match(samps, pro_summary$sample)] <- rep(c("GF/F", "Swinnex"), times = 5)

pro_summary$sample_num <- NA

pro_summary$sample_num[match(samps, pro_summary$sample)] <- rep(1:5, each = 2)

pro_summary$group <- "Prochlorococcus"

wide_pro <- pro_summary %>% filter(!is.na(sample_num)) %>%
  pivot_wider(-sample, names_from = type, values_from = c(total_reads, richness))

total_reads_pro <- t.test(wide_pro$`total_reads_GF/F`, wide_pro$total_reads_Swinnex, paired = TRUE)
richness_pro <- t.test(wide_pro$`richness_GF/F`, wide_pro$richness_Swinnex, paired = TRUE)

sar_pro_summary <- bind_rows(sar_summary, pro_summary)

sar_pro_summary$prop_16s <- sar_pro_summary$total_reads/six_summary$total_reads[match(sar_pro_summary$sample, six_summary$sample)]

sar_pro_summary$type <- as.factor(sar_pro_summary$type)
sar_pro_summary$type <- factor(sar_pro_summary$type, 
                             levels = c("GF/F", "Swinnex", 
                                        "GF/F (201611 Cruise)"))


richness1 <- ggplot(sar_pro_summary, aes(x = group, y = richness, fill = type)) +
  geom_boxplot() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "", y = "Richness", fill = "Filter Type") +
  annotate("text", x = 0.9, y = 25, label = paste0("p-value = ", round(richness_pro$p.value,3))) +
  annotate("text", x = 1.9, y = 55, label = paste0("p-value = ", round(richness_sar$p.value,3))) +
  scale_fill_manual(values = c("#9673c5","#ba6d75","#a09342"))

total_reads1 <- ggplot(sar_pro_summary, aes(x = group, y = total_reads, fill = type)) +
  geom_boxplot() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "", y = "# of Amplicons", fill = "Filter Type") +
  annotate("text", x = 0.9, y = 15000, label = paste0("p-value = ", round(total_reads_pro$p.value,3))) +
  annotate("text", x = 1.9, y = 20000, label = paste0("p-value = ", round(total_reads_sar$p.value,3))) +
  scale_fill_manual(values = c("#9673c5","#ba6d75","#a09342"))

prop_reads1 <- ggplot(sar_pro_summary, aes(x = group, y = prop_16s*100, fill = type)) +
  geom_boxplot() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "", y = "% of 16S V4-V5 Amplicons", fill = "Filter Type") +
  annotate("text", x = 0.9, y = 15, label = paste0("p-value = ", round(total_reads_pro$p.value,3))) +
  annotate("text", x = 1.9, y = 30, label = paste0("p-value = ", round(total_reads_sar$p.value,3))) +
  scale_fill_manual(values = c("#9673c5","#ba6d75","#a09342"))

out1 <- total_reads1 + richness1 + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect")

# pdf(file = "figures/figure_outline/supp_fig_pro_sar_XX.pdf", width = 8.5, height = 5)
# print(out1)
# dev.off()

#### PCOA ######

load("output/cyano_16s_full_data.Rdata")

six$sample <- NULL

bc_six <- vegdist(six, method = "bray")

six_pcoa <- pco(bc_six, negvals = "zero", dround = 0)

six_pcoa_ggplot <- as.data.frame(six_pcoa$vectors[,1:2])
six_pcoa_ggplot$sample <- samples
six_pcoa_ggplot$type <- "GF/F (201611 Cruise)"
six_pcoa_ggplot$type[match(samps, six_pcoa_ggplot$sample)] <- rep(c("GF/F", "Swinnex"), times = 5)

six_pcoa_ggplot$type <- as.factor(six_pcoa_ggplot$type)
six_pcoa_ggplot$type <- factor(six_pcoa_ggplot$type, 
                             levels = c("GF/F", "Swinnex", 
                                        "GF/F (201611 Cruise)"))

six_pcoa_ggplot$sample_num <- NA

six_pcoa_ggplot$sample_num[match(samps, six_pcoa_ggplot$sample)] <- rep(1:5, each = 2)

six_pcoa_ggplot$dist <- full_dat$dist_to_coast[match(gsub("_T", "",six_pcoa_ggplot$sample),full_dat$eco_name)]

six_plot <- ggplot() +
  geom_point(data = six_pcoa_ggplot,
             aes(x = V1, y = V2, shape = type, color = dist), size = 2) +
  geom_text(data = six_pcoa_ggplot %>% filter(type != "GF/F (201611 Cruise)"),
            aes(x = V1+0.009, y = V2+0.015, label = sample_num, color = dist), size = 2.5) +
  scale_color_viridis() + scale_shape_manual(values = c(3,4,19)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "PCoA1", y = "PCoA2", color = "Distance to Coast (km)") +
  ggtitle("16S V4-V5")

# eight

eight$sample <- NULL

bc_eight <- vegdist(eight, method = "bray")

eight_pcoa <- pco(bc_eight, negvals = "zero", dround = 0)

eight_pcoa_ggplot <- as.data.frame(eight_pcoa$vectors[,1:2])
eight_pcoa_ggplot$sample <- samples
eight_pcoa_ggplot$type <- "GF/F (201611 Cruise)"
eight_pcoa_ggplot$type[match(samps, eight_pcoa_ggplot$sample)] <- rep(c("GF/F", "Swinnex"), times = 5)

eight_pcoa_ggplot$type <- as.factor(eight_pcoa_ggplot$type)
eight_pcoa_ggplot$type <- factor(eight_pcoa_ggplot$type, 
                               levels = c("GF/F", "Swinnex", 
                                          "GF/F (201611 Cruise)"))

eight_pcoa_ggplot$sample_num <- NA

eight_pcoa_ggplot$sample_num[match(samps, eight_pcoa_ggplot$sample)] <- rep(1:5, each = 2)

eight_pcoa_ggplot$dist <- full_dat$dist_to_coast[match(gsub("_T", "",eight_pcoa_ggplot$sample),full_dat$eco_name)]

eight_plot <- ggplot() +
  geom_point(data = eight_pcoa_ggplot,
             aes(x = V1, y = V2, shape = type, color = dist), size = 2) +
  geom_text(data = eight_pcoa_ggplot %>% filter(type != "GF/F (201611 Cruise)"),
            aes(x = V1+0.009, y = V2+0.015, label = sample_num, color = dist), size = 2.5) +
  scale_color_viridis() + scale_shape_manual(values = c(3,4,19)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "PCoA1", y = "PCoA2", color = "Distance to Coast (km)") +
  ggtitle("18Sv9")


###### Fixed Plots #####

total_summary$prop <- total_summary$total_reads/sum(total_summary$total_reads)

richness_new <- ggplot(total_summary, aes(x = region, y = richness, fill = type)) +
  geom_boxplot() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "", y = "Richness", fill = "Filter Type") +
  annotate("text", x = 0.9, y = 900, label = paste0("p-value = ", round(richness_six$p.value,3))) +
  annotate("text", x = 1.9, y = 500, label = paste0("p-value = ", round(richness_eight$p.value,3))) +
  scale_fill_manual(values = c("#9673c5","#ba6d75","#a09342"))

perc_total <- ggplot(total_summary, aes(x = region, y = total_reads, fill = type)) +
  geom_boxplot() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "", y = "# of amplicons", fill = "Filter Type") +
  annotate("text", x = 0.9, y = 145000, label = paste0("p-value = ", round(total_reads_six$p.value,3))) +
  annotate("text", x = 1.9, y = 50000, label = paste0("p-value = ", round(total_reads_eight$p.value,3))) +
  scale_fill_manual(values = c("#9673c5","#ba6d75","#a09342"))

out_full <- perc_total + richness_new + six_plot +
  prop_reads1 + richness1 + eight_plot + plot_layout(widths = c(1,1,1))

# pdf(file = "figures/figure_outline/supp_fig_XX.pdf", width = 17, height = 8)
# print(out_full)
# dev.off()


###### proportion plots #####

##### combine for plot #####

six_summary1$region <- "16SV4-V5"
eight_summary1$region <- "18SV9"

six_summary1 <- six_summary1[,c(1,6:10)]
eight_summary1 <- eight_summary1[,c(1,6,7,8,9,10)]

colnames(six_summary1)[4:5] <- c("prop_reads", "prop_rich")
colnames(eight_summary1)[4:5] <- c("prop_reads", "prop_rich")

total_summary1 <- bind_rows(six_summary1, eight_summary1)

total_summary1$type <- as.factor(total_summary1$type)
total_summary1$type <- factor(total_summary1$type, 
                             levels = c("GF/F", "Swinnex", 
                                        "GF/F (201611 Cruise)"))

richness_prop <- ggplot(total_summary1, aes(x = region, y = prop_rich*100, fill = type)) +
  geom_boxplot() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "", y = "% of ASVs\nprokaryotic (16S) or eukaryotic (18S)", fill = "Filter Type") +
  annotate("text", x = 0.9, y = 95, label = paste0("p-value = ", round(prop_richness_six$p.value,3))) +
  annotate("text", x = 1.9, y = 90, label = paste0("p-value = ", round(prop_richness_eight$p.value,3))) +
  scale_fill_manual(values = c("#9673c5","#ba6d75","#a09342"))

total_prop <- ggplot(total_summary1, aes(x = region, y = prop_reads*100, fill = type)) +
  geom_boxplot() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "", y = "% of amplicons\nprokaryotic (16S) or eukaryotic (18S)", fill = "Filter Type") +
  annotate("text", x = 0.9, y = 100, label = paste0("p-value = ", round(prop_reads_six$p.value,3))) +
  annotate("text", x = 1.9, y = 100, label = paste0("p-value = ", round(prop_reads_eight$p.value,3))) +
  scale_fill_manual(values = c("#9673c5","#ba6d75","#a09342"))


out_full <- total_prop + richness_prop + six_plot +
  prop_reads1 + richness1 + eight_plot + plot_layout(widths = c(1,1,1))

# pdf(file = "figures/figure_outline/supp_fig_XX2.pdf", width = 17, height = 8)
# print(out_full)
# dev.off()

###### All Groups #####

in_group_list_basic = c("16s_pro", "16s_syne","16s_flavo", "16s_rhodo", "16s_sar", 
                        "18s_diatom","18s_dino", "18s_syndin",
                        "18s_hapto", "18s_chloro", "18s_metazoa")
in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar 11 Clade", 
                   "Diatoms", "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans")

regions <- c("six", "six", "six", "six", "six", "eight","eight","eight","eight","eight","eight")


data_list <- list()

samps <- c("201611_093.3_026.7_24", "201611_093.3_026.7_24_T", "201611_093.3_045.0_33", "201611_093.3_045.0_33_T",
           "201611_093.3_070.0_40", "201611_093.3_070.0_40_T", "201611_086.7_070.0_10","201611_086.7_070.0_10_T",
           "201611_076.7_070.0_8","201611_076.7_070.0_8_T")

samps <- paste0("X", samps)

for (i in 1:length(in_group_list_basic)) {
  
  load(paste0("data/",in_group_list_basic[i],".Rdata"))
  
  if(regions[i] == "six"){group_long <- six_long[which(!is.na(match(six_long$name, colnames(asv_table)))),]}
  if(regions[i] == "eight"){group_long <- eight_long[which(!is.na(match(eight_long$name, colnames(asv_table)))),]}
  
  group_summary <- group_long %>%
    group_by(sample) %>%
    summarise(total_reads = sum(value, na.rm = TRUE), richness = sum(occur, na.rm = TRUE))
  
  group_summary$type <- "GF/F (201611 Cruise)"
  
  group_summary$type[match(samps, group_summary$sample)] <- rep(c("GF/F", "Swinnex"), times = 5)
  
  group_summary$sample_num <- NA
  
  group_summary$sample_num[match(samps, group_summary$sample)] <- rep(1:5, each = 2)
  
  group_summary$group <- in_group_names[i]
  
  if(regions[i] == "six"){
    group_summary$prop_reads <- group_summary$total_reads/six_summary$total_reads[match(group_summary$sample, six_summary$sample)]
    group_summary$prop_richness <- group_summary$richness/six_summary$richness[match(group_summary$sample, six_summary$sample)]
    }
  
  if(regions[i] == "eight"){
    group_summary$prop_reads <- group_summary$total_reads/eight_summary$total_reads[match(group_summary$sample, eight_summary$sample)]
    group_summary$prop_richness <- group_summary$richness/eight_summary$richness[match(group_summary$sample, eight_summary$sample)]
    }
  
  wide_group <- group_summary[,-c(2:3)] %>% filter(!is.na(sample_num)) %>%
    pivot_wider(-sample, names_from = type, values_from = c(prop_reads, prop_richness))
  
  prop_reads_group <- t.test(wide_group$`prop_reads_GF/F`, wide_group$prop_reads_Swinnex, paired = TRUE)
  prop_richness_group <- t.test(wide_group$`prop_richness_GF/F`, wide_group$prop_richness_Swinnex, paired = TRUE)
  
  data_list[[i]] <- list(group_df = group_summary, p_reads = prop_reads_group$p.value, p_richness = prop_richness_group$p.value)
  
  print(i)
  
}

total_df <- bind_rows(data_list[[1]]$group_df, data_list[[2]]$group_df)
total_df <- bind_rows(total_df, data_list[[3]]$group_df)
total_df <- bind_rows(total_df, data_list[[4]]$group_df)
total_df <- bind_rows(total_df, data_list[[5]]$group_df)
total_df <- bind_rows(total_df, data_list[[6]]$group_df)
total_df <- bind_rows(total_df, data_list[[7]]$group_df)
total_df <- bind_rows(total_df, data_list[[8]]$group_df)
total_df <- bind_rows(total_df, data_list[[9]]$group_df)
total_df <- bind_rows(total_df, data_list[[10]]$group_df)
total_df <- bind_rows(total_df, data_list[[11]]$group_df)

total_df$group <- as.factor(total_df$group)
total_df$group <- factor(total_df$group, levels = in_group_names)

total_df$type <- as.factor(total_df$type)
total_df$type <- factor(total_df$type, 
                             levels = c("GF/F", "Swinnex", 
                                        "GF/F (201611 Cruise)"))

rel_group_rich <- ggplot(total_df, aes(x = group, y = prop_richness*100, fill = type)) +
  geom_boxplot() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "", y = "% of 16S V4-V5/18SV9 ASVs\n(relative group richness)", fill = "Filter Type") +
  annotate("text", x = 0.9, y = 10, label = paste0("p-value = ", round(data_list[[1]]$p_richness,3))) +
  annotate("text", x = 1.9, y = 10, label = paste0("p-value = ", round(data_list[[2]]$p_richness,3))) +
  annotate("text", x = 2.9, y = 30, label = paste0("p-value = ", round(data_list[[3]]$p_richness,3))) +
  annotate("text", x = 3.9, y = 10, label = paste0("p-value = ", round(data_list[[4]]$p_richness,3))) +
  annotate("text", x = 4.9, y = 20, label = paste0("p-value = ", round(data_list[[5]]$p_richness,3))) +
  annotate("text", x = 5.9, y = 20, label = paste0("p-value = ", round(data_list[[6]]$p_richness,3))) +
  annotate("text", x = 6.9, y = 60, label = paste0("p-value = ", round(data_list[[7]]$p_richness,3))) +
  annotate("text", x = 7.9, y = 35, label = paste0("p-value = ", round(data_list[[8]]$p_richness,3))) +
  annotate("text", x = 8.9, y = 10, label = paste0("p-value = ", round(data_list[[9]]$p_richness,3))) +
  annotate("text", x = 9.9, y = 10, label = paste0("p-value = ", round(data_list[[10]]$p_richness,3))) +
  annotate("text", x = 10.9, y = 15, label = paste0("p-value = ", round(data_list[[11]]$p_richness,3))) +
  scale_fill_manual(values = c("#9673c5","#ba6d75","#a09342"))

rel_group_abun <- ggplot(total_df, aes(x = group, y = prop_reads*100, fill = type)) +
  geom_boxplot() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  labs(x = "", y = "% of 16S V4-V5/18SV9 amplicons\n(relative group abundance)", fill = "Filter Type") +
  annotate("text", x = 0.9, y = 45, label = paste0("p-value = ", round(data_list[[1]]$p_reads,3))) +
  annotate("text", x = 1.9, y = 20, label = paste0("p-value = ", round(data_list[[2]]$p_reads,3))) +
  annotate("text", x = 2.9, y = 30, label = paste0("p-value = ", round(data_list[[3]]$p_reads,3))) +
  annotate("text", x = 3.9, y = 15, label = paste0("p-value = ", round(data_list[[4]]$p_reads,3))) +
  annotate("text", x = 4.9, y = 40, label = paste0("p-value = ", round(data_list[[5]]$p_reads,3))) +
  annotate("text", x = 5.9, y = 25, label = paste0("p-value = ", round(data_list[[6]]$p_reads,3))) +
  annotate("text", x = 6.9, y = 75, label = paste0("p-value = ", round(data_list[[7]]$p_reads,3))) +
  annotate("text", x = 7.9, y = 30, label = paste0("p-value = ", round(data_list[[8]]$p_reads,3))) +
  annotate("text", x = 8.9, y = 10, label = paste0("p-value = ", round(data_list[[9]]$p_reads,3))) +
  annotate("text", x = 9.9, y = 10, label = paste0("p-value = ", round(data_list[[10]]$p_reads,3))) +
  annotate("text", x = 10.9, y = 70, label = paste0("p-value = ", round(data_list[[11]]$p_reads,3))) +
  scale_fill_manual(values = c("#9673c5","#ba6d75","#a09342"))


layout <- "AABB
           CCCC
           DDDD
           EEFF"


out_total_group <- (total_prop + theme(legend.position = "none")) +
  (richness_prop + theme(legend.position = "none")) +
  rel_group_abun +
  (rel_group_rich + theme(legend.position = "none")) + 
  (six_plot + theme(legend.position = "none")) + eight_plot +
  plot_layout(design = layout) + plot_annotation(tag_levels = "a")





pdf("figures/figure_outline/supp_fig_15.pdf", width = 14, height = 18)
print(out_total_group)
dev.off()

