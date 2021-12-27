library(tidyverse)
library(lubridate)
library(ragg)

# sixteen s

load("data/all_16S_rare.Rdata")
metadata <- read.csv("data/NCOG_sample_log_DNA_stvx_meta_2014-2020.csv")

metadata$Date <- mdy(metadata$Date)
dates <- metadata %>% group_by(Cruise) %>% summarise(mean_date = mean(Date, na.rm = TRUE))

split_taxa_six <- separate(six_tax_id, silva_Taxon, sep = "; ", into = c("A","B","C", "D", "E", "F", "G"))

flavobacteriales <- which(split_taxa_six$D == "o__Flavobacteriales")
rhodobacterales <- which(split_taxa_six$D == "o__Rhodobacterales")
sar_11 <- which(split_taxa_six$D == "o__SAR11_clade")
prochlorococcus <- which(split_taxa_six$B == "p__Cyanobacteria" & substr(split_taxa_six$F,1,12) == "g__Prochloro")
synechococcus <- which(split_taxa_six$B == "p__Cyanobacteria" & substr(split_taxa_six$F,1,10) == "g__Synecho")
cyanobacteria <- which(split_taxa_six$B == "p__Cyanobacteria")
archaea <- which(split_taxa_six$A == "d__Archaea")

six_rare$sample <- rownames(six_rare)

six_pivot <- six_rare %>% pivot_longer(-sample, names_to = "Hash", values_to = "reads")

six_pivot$Group <- "Other Bacteria"

six_pivot$Group[which(six_pivot$Hash %in% split_taxa_six$Feature.ID[flavobacteriales])] <- "Flavobacteriales"
six_pivot$Group[which(six_pivot$Hash %in% split_taxa_six$Feature.ID[rhodobacterales])] <- "Rhodobacterales"
six_pivot$Group[which(six_pivot$Hash %in% split_taxa_six$Feature.ID[sar_11])] <- "Sar 11 Clade"
six_pivot$Group[which(six_pivot$Hash %in% split_taxa_six$Feature.ID[prochlorococcus])] <- "Prochlorococcus"
six_pivot$Group[which(six_pivot$Hash %in% split_taxa_six$Feature.ID[synechococcus])] <- "Synechococcus"
six_pivot$Group[which(six_pivot$Hash %in% split_taxa_six$Feature.ID[archaea])] <- "Archaea"

six_pivot$cruise  <- metadata$Cruise[match(six_pivot$sample, paste0("X", metadata$Sample.Name))]

six_rich <- six_pivot %>% filter(reads != 0)

# richness
six_ts_rich <- six_rich %>% group_by(Group, cruise) %>% summarise(grp_rich = n_distinct(Hash))

six_ts <- six_pivot %>% group_by(Group, cruise) %>% summarise(grp_reads = sum(reads,na.rm = TRUE))
six_reads <- six_pivot %>% group_by(cruise) %>% summarise(grp_reads = sum(reads,na.rm = TRUE))

six_ts$rel_abun <- six_ts$grp_reads/six_reads$grp_reads[match(six_ts$cruise, six_reads$cruise)]

six_ts$Group <- factor(six_ts$Group, levels = c("Archaea", "Flavobacteriales",
                                                "Prochlorococcus", "Rhodobacterales",
                                                "Sar 11 Clade", "Synechococcus",
                                                "Other Bacteria"))

six_ts_rich$Group <- factor(six_ts_rich$Group, levels = c("Archaea", "Flavobacteriales",
                                                "Prochlorococcus", "Rhodobacterales",
                                                "Sar 11 Clade", "Synechococcus",
                                                "Other Bacteria"))

six_ts$Date <- dates$mean_date[match(six_ts$cruise, dates$Cruise)]
six_ts_rich$Date <- dates$mean_date[match(six_ts_rich$cruise, dates$Cruise)]

six_plot <- ggplot(six_ts, aes(x = Date, y = rel_abun, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#969e3d", "#b260bd", "#50ac71",
                               "#c95574","#6783d0","#c9713e", "grey60")) +
  scale_x_date(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  theme(panel.border = element_rect(fill=NA, color = "black")) +
  labs(y = "Group Relative Abundance") + ggtitle("16S")

six_rich_plot <- ggplot(six_ts_rich, aes(x = Date, y = grp_rich, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#969e3d", "#b260bd", "#50ac71",
                               "#c95574","#6783d0","#c9713e", "grey60")) +
  scale_x_date(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  theme(panel.border = element_rect(fill=NA, color = "black"),
        panel.background = element_blank()) +
  labs(y = "Total Richness") + ggtitle("16S")

six_rich_plot <- ggplot(six_ts_rich %>% filter(Group != "Other Bacteria"), aes(x = Date, y = grp_rich, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#969e3d", "#b260bd", "#50ac71",
                               "#c95574","#6783d0","#c9713e", "grey60")) +
  scale_x_date(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  theme(panel.border = element_rect(fill=NA, color = "black"),
        panel.background = element_blank()) +
  labs(y = "Total Richness") + ggtitle("16S")

agg_png("figures_S/group_changes_16S_S.png", width = 10, height = 7, units = "in", res = 400)
plot(six_plot)
dev.off()

# eighteen s

load("data/all_18Sv9_rare.Rdata")

split_taxa_eight <- separate(eight_tax_id, pr2_Taxon, sep = ";", into = c("A","B","C", "D", "E", "F", "G"))

diatoms <- which(split_taxa_eight$D == "Bacillariophyta")
dinos_minus_syn <- which(split_taxa_eight$C == "Dinoflagellata" &
                           split_taxa_eight$D != "Syndiniales" & !is.na(split_taxa_eight$D))
sindins <- which(split_taxa_eight$D == "Syndiniales")
haptophytes <- which(split_taxa_eight$C == "Haptophyta")
chlorophytes <- which(split_taxa_eight$C == "Chlorophyta")
metazoa <- which(split_taxa_eight$C == "Metazoa")

eight_rare$sample <- rownames(eight_rare)

eight_pivot <- eight_rare %>% pivot_longer(-sample, names_to = "Hash", values_to = "reads")

eight_pivot$Group <- "Other Eukaryotes"

eight_pivot$Group[which(eight_pivot$Hash %in% split_taxa_eight$Feature.ID[diatoms])] <- "Diatoms"
eight_pivot$Group[which(eight_pivot$Hash %in% split_taxa_eight$Feature.ID[dinos_minus_syn])] <- "Dinoflagellates"
eight_pivot$Group[which(eight_pivot$Hash %in% split_taxa_eight$Feature.ID[sindins])] <- "Syndiniales"
eight_pivot$Group[which(eight_pivot$Hash %in% split_taxa_eight$Feature.ID[haptophytes])] <- "Haptophytes"
eight_pivot$Group[which(eight_pivot$Hash %in% split_taxa_eight$Feature.ID[chlorophytes])] <- "Chlorophytes"
eight_pivot$Group[which(eight_pivot$Hash %in% split_taxa_eight$Feature.ID[metazoa])] <- "Metazoa"

eight_pivot$cruise  <- metadata$Cruise[match(eight_pivot$sample, paste0("X", metadata$Sample.Name))]

eight_ts <- eight_pivot %>% group_by(Group, cruise) %>% summarise(grp_reads = sum(reads,na.rm = TRUE))

eight_rich <- eight_pivot %>% filter(reads != 0)
eight_ts_rich <- eight_rich %>% group_by(Group, cruise) %>% summarise(grp_rich = n_distinct(Hash))

eight_reads <- eight_pivot %>% group_by(cruise) %>% summarise(grp_reads = sum(reads,na.rm = TRUE))

eight_ts$rel_abun <- eight_ts$grp_reads/eight_reads$grp_reads[match(eight_ts$cruise, eight_reads$cruise)]

eight_ts$Group <- factor(eight_ts$Group, levels = c("Diatoms", "Dinoflagellates",
                                                "Chlorophytes", "Haptophytes",
                                                "Metazoa", "Syndiniales",
                                                "Other Eukaryotes"))


eight_ts_rich$Group <- factor(eight_ts_rich$Group, levels = c("Diatoms", "Dinoflagellates",
                                                    "Chlorophytes", "Haptophytes",
                                                    "Metazoa", "Syndiniales",
                                                    "Other Eukaryotes"))

eight_ts$Date <- dates$mean_date[match(eight_ts$cruise, dates$Cruise)]
eight_ts_rich$Date <- dates$mean_date[match(eight_ts_rich$cruise, dates$Cruise)]

eight_plot <- ggplot(eight_ts, aes(x = Date, y = rel_abun, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#969e3d", "#b260bd", "#50ac71",
                               "#c95574","#6783d0","#c9713e", "grey60")) +
  scale_x_date(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  theme(panel.border = element_rect(fill=NA, color = "black")) +
  labs(y = "Group Relative Abundance") + ggtitle("18Sv9")

eight_rich_plot <- ggplot(eight_ts_rich, aes(x = Date, y = grp_rich, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#969e3d", "#b260bd", "#50ac71",
                               "#c95574","#6783d0","#c9713e", "grey60")) +
  scale_x_date(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  theme(panel.border = element_rect(fill=NA, color = "black"),
        panel.background = element_blank()) +
  labs(y = "Group Relative Abundance") + ggtitle("18Sv9")

eight_rich_plot <- ggplot(eight_ts_rich %>% filter(Group != "Other Eukaryotes"), aes(x = Date, y = grp_rich, fill = Group)) +
  geom_area() +
  scale_fill_manual(values = c("#969e3d", "#b260bd", "#50ac71",
                               "#c95574","#6783d0","#c9713e", "grey60")) +
  scale_x_date(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
  theme(panel.border = element_rect(fill=NA, color = "black"),
        panel.background = element_blank()) +
  labs(y = "Group Relative Abundance") + ggtitle("18Sv9")

agg_png("figures_S/group_changes_18Sv9_S.png", width = 10, height = 7, units = "in", res = 400)
plot(eight_plot)
dev.off()
