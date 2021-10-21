library(tidyverse)
library(vegan)
library(ragg)
library(patchwork)

load("data/18sv9_all.Rdata")
metadata <- read.csv("data/NCOG_sample_log_DNA_stvx_meta_2014-2020.csv")

eighteen_s$sample <- gsub(pattern = "X","",rownames(eighteen_s))

metadata <- metadata %>% filter(as.numeric(substr(sample_num,3,6)) > 478)

eighteen_s <- eighteen_s %>% filter(sample %in% metadata$Sample.Name)
metadata <- metadata %>% filter(Sample.Name %in% eighteen_s$sample)

metadata$library_size <- rowSums(eighteen_s[,-ncol(eighteen_s)])[match(metadata$Sample.Name, eighteen_s$sample)]

metadata$sequencing <- "Previous Sequencing"
metadata$sequencing[which(as.numeric(substr(metadata$sample_num,3,6)) > 1150)] <- "NextSeq"

lib_size <- ggplot(metadata, aes(x = as.numeric(substr(sample_num,3,6)), y = library_size, color = sequencing)) +
  geom_point() + scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "Sample Number", y = "Library Size", color = "") +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"))

agg_png("figures/library_sizes_18Sv9.png", width = 8, height = 6, units = "in", res = 400)
plot(lib_size)
dev.off()

eight_rare <- rrarefy(eighteen_s[,-ncol(eighteen_s)], 17000)

raw_rich <- apply(eighteen_s[,-ncol(eighteen_s)],1,function(x) length(!is.na(which(x!=0))))
rare_rich <- apply(eight_rare,1,function(x) length(!is.na(which(x!=0))))

metadata$raw_rich <- raw_rich[match(paste0("X",metadata$Sample.Name), names(raw_rich))]
metadata$rare_rich <- rare_rich[match(paste0("X",metadata$Sample.Name), names(rare_rich))]


no_rare <- ggplot(metadata, aes(x = NCDepth, y = raw_rich, color = sequencing)) +
  geom_point() + scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "Nitracline Depth", y = "Richness", color = "") +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black")) +
  stat_smooth(method = "lm") + ggtitle("No Rarefaction")


rare <- ggplot(metadata, aes(x = NCDepth, y = rare_rich, color = sequencing)) +
  geom_point() + scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "Nitracline Depth", y = "Richness", color = "") +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black")) +
  stat_smooth(method = "lm") + ggtitle("Rarefy to 17,000 reads")

out <- no_rare / rare

agg_png("figures/rare_richness_18Sv9.png", width = 8, height = 10, units = "in", res = 400)
plot(out)
dev.off()

# by group

split_taxa <- separate(eight_tax_id, pr2_Taxon, sep = ";", into = c("A","B","C", "D", "E", "F", "G", "H", "I"))

split_taxa$group <- "Other"

split_taxa$group[which(split_taxa$C == "Chlorophyta")] <- "Chlorophytes"
split_taxa$group[which(split_taxa$C == "Dinoflagellata" &
                         split_taxa$D != "Syndiniales" | is.na(split_taxa$D))] <- "Dinoflagellates"
split_taxa$group[which(split_taxa$C == "Cryptophyta")] <- "Cryptophytes"
split_taxa$group[which(split_taxa$C == "Haptophyta")] <- "Haptophytes"
split_taxa$group[which(split_taxa$D == "Syndiniales")] <- "Syndiniales"
split_taxa$group[which(split_taxa$D == "Bacillariophyta")] <- "Diatoms"
split_taxa$group[which(split_taxa$C == "Metazoa")] <- "Metazoa"

eight_long <- eighteen_s %>% pivot_longer(-sample, names_to = "Hash", values_to = "reads")

eight_rare <- as.data.frame(eight_rare)
eight_rare$sample <- rownames(eight_rare)

eight_long_rare <- eight_rare %>% pivot_longer(-sample, names_to = "Hash", values_to = "reads")

eight_long$group <- split_taxa$group[match(eight_long$Hash, split_taxa$Feature.ID)]
eight_long_rare$group <- split_taxa$group[match(eight_long_rare$Hash, split_taxa$Feature.ID)]

eight_long_rich <- eight_long %>% group_by(sample, group) %>% summarise(rich = length(!is.na(which(reads !=  0))))
eight_long_rare_rich <- eight_long_rare %>% group_by(sample, group) %>% summarise(rich = length(!is.na(which(reads !=  0))))

eight_long_rich$NCDepth <- metadata$NCDepth[match(eight_long_rich$sample, metadata$Sample.Name)]
eight_long_rare_rich$NCDepth <- metadata$NCDepth[match(eight_long_rare_rich$sample, paste0("X",metadata$Sample.Name))]

eight_long_rich$seq <- metadata$sequencing[match(eight_long_rich$sample, metadata$Sample.Name)]
eight_long_rare_rich$seq <- metadata$sequencing[match(eight_long_rare_rich$sample, paste0("X",metadata$Sample.Name))]

eight_long_rich$group <- factor(eight_long_rich$group,
                                     levels = c("Chlorophytes", "Cryptophytes",
                                                "Diatoms", "Dinoflagellates",
                                                "Haptophytes", "Metazoa", "Syndiniales", "Other"))

eight_long_rare_rich$group <- factor(eight_long_rare_rich$group,
                                   levels = c("Chlorophytes", "Cryptophytes",
                                              "Diatoms", "Dinoflagellates",
                                              "Haptophytes", "Metazoa", "Syndiniales", "Other"))


no_rare <- ggplot(eight_long_rich, aes(x = NCDepth, y = rich, color = seq)) +
  geom_point() + facet_wrap(~group, scales = "free_y") +
  geom_smooth(method = "lm") + ggtitle("No Rarefaction") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        legend.position=c(.85,.15),
        legend.key = element_blank()) +
  labs(x = "Nitracline Depth (m)", y = "Richness", color = "")


rare <- ggplot(eight_long_rare_rich, aes(x = NCDepth, y = rich, color = seq)) +
  geom_point() + facet_wrap(~group, scales = "free_y") +
  geom_smooth(method = "lm") + ggtitle("Rareified to 17,000 reads") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        legend.position=c(.85,.15),
        legend.key = element_blank()) +
  labs(x = "Nitracline Depth (m)", y = "Richness", color = "")


agg_png("figures/richness_18Sv9_groups.png", width = 12, height = 10, units = "in", res = 400)
plot(no_rare)
dev.off()

agg_png("figures/richness_18Sv9_groups_rare.png", width = 12, height = 10, units = "in", res = 400)
plot(rare)
dev.off()


### Sterivex vs GF/F ######

load("data/18sv9_all.Rdata")
metadata <- read.csv("data/NCOG_sample_log_DNA_stvx_meta_2014-2020.csv")

eighteen_s$sample <- gsub(pattern = "X","",rownames(eighteen_s))


metadata <- metadata %>% filter(as.numeric(substr(sample_num,3,6)) < 478 |
                                  as.numeric(substr(sample_num,3,6)) > 1150)


eighteen_s <- eighteen_s %>% filter(sample %in% metadata$Sample.Name)
metadata <- metadata %>% filter(Sample.Name %in% eighteen_s$sample)

metadata$library_size <- rowSums(eighteen_s[,-ncol(eighteen_s)])[match(metadata$Sample.Name, eighteen_s$sample)]

metadata$sequencing <- "Previous Sequencing GF/F"
metadata$sequencing[which(as.numeric(substr(metadata$sample_num,3,6)) < 1150 &
                            as.numeric(substr(metadata$sample_num,3,6)) > 478)] <- "Previous Sequencing Sterivex"
metadata$sequencing[which(as.numeric(substr(metadata$sample_num,3,6)) > 1150)] <- "NextSeq Sterivex"

lib_size <- ggplot(metadata, aes(x = as.numeric(substr(sample_num,3,6)), y = library_size, color = sequencing)) +
  geom_point() + scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "Sample Number", y = "Library Size", color = "") +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black")) +
  ggtitle("GF/F vs Sterivex")

agg_png("figures/library_sizes_18Sv9_all_compare.png", width = 8, height = 6, units = "in", res = 400)
plot(lib_size)
dev.off()

eight_rare <- rrarefy(eighteen_s[,-ncol(eighteen_s)], 7100)

raw_rich <- apply(eighteen_s[,-ncol(eighteen_s)],1,function(x) length(!is.na(which(x!=0))))
rare_rich <- apply(eight_rare,1,function(x) length(!is.na(which(x!=0))))

eight_rare <- eight_rare %>% as.data.frame()
eight_rare$sample <- eighteen_s$sample

eighteen_s$cruise <- metadata$Cruise[match(eighteen_s$sample,metadata$Sample.Name)]
eight_rare$cruise <- metadata$Cruise[match(eight_rare$sample,metadata$Sample.Name)]

eight_no_rare <- eighteen_s %>% pivot_longer(-c(cruise,sample)) 
eight_rare <- eight_rare %>% pivot_longer(-c(cruise,sample))

eight_no_rare$seq <- metadata$sequencing[match(eight_no_rare$sample, metadata$Sample.Name)]
eight_rare$seq <- metadata$sequencing[match(eight_rare$sample, metadata$Sample.Name)]

eight_no_rare <- eight_no_rare %>% filter(value != 0) %>%
  group_by(cruise,seq) %>% summarise(rich = n_distinct(name))

eight_rare <- eight_rare %>% filter(value != 0) %>%
  group_by(cruise,seq) %>% summarise(rich = n_distinct(name))

eight_no_rare$rare <- "Raw Richness"
eight_rare$rare <- "Rarefied Richness"

metadata$raw_rich <- raw_rich[match(paste0("X",metadata$Sample.Name), names(raw_rich))]
metadata$rare_rich <- rare_rich[match(paste0("X",metadata$Sample.Name), names(rare_rich))]



no_rare <- ggplot(metadata, aes(x = NCDepth, y = raw_rich, color = sequencing)) +
  geom_point() + scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "Nitracline Depth", y = "Richness", color = "") +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black")) +
  stat_smooth(method = "lm") + ggtitle("No Rarefaction")


rare <- ggplot(metadata, aes(x = NCDepth, y = rare_rich, color = sequencing)) +
  geom_point() + scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "Nitracline Depth", y = "Richness", color = "") +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black")) +
  stat_smooth(method = "lm") + ggtitle("Rarefy to 17,000 reads")

ggplot(eight_no_rare, aes(x = cruise, y = rich, color = seq)) +
  geom_point()

ggplot(eight_rare, aes(x = cruise, y = rich, color = seq)) +
  geom_point()




out <- no_rare / rare

agg_png("figures/rare_richness_18Sv9.png", width = 8, height = 10, units = "in", res = 400)
plot(out)
dev.off()

# by group

split_taxa <- separate(eight_tax_id, pr2_Taxon, sep = ";", into = c("A","B","C", "D", "E", "F", "G", "H", "I"))

split_taxa$group <- "Other"

split_taxa$group[which(split_taxa$C == "Chlorophyta")] <- "Chlorophytes"
split_taxa$group[which(split_taxa$C == "Dinoflagellata" &
                         split_taxa$D != "Syndiniales" | is.na(split_taxa$D))] <- "Dinoflagellates"
split_taxa$group[which(split_taxa$C == "Cryptophyta")] <- "Cryptophytes"
split_taxa$group[which(split_taxa$C == "Haptophyta")] <- "Haptophytes"
split_taxa$group[which(split_taxa$D == "Syndiniales")] <- "Syndiniales"
split_taxa$group[which(split_taxa$D == "Bacillariophyta")] <- "Diatoms"
split_taxa$group[which(split_taxa$C == "Metazoa")] <- "Metazoa"

eight_long <- eighteen_s %>% pivot_longer(-sample, names_to = "Hash", values_to = "reads")

eight_rare <- as.data.frame(eight_rare)
eight_rare$sample <- rownames(eight_rare)

eight_long_rare <- eight_rare %>% pivot_longer(-sample, names_to = "Hash", values_to = "reads")

eight_long$group <- split_taxa$group[match(eight_long$Hash, split_taxa$Feature.ID)]
eight_long_rare$group <- split_taxa$group[match(eight_long_rare$Hash, split_taxa$Feature.ID)]

eight_long_rich <- eight_long %>% group_by(sample, group) %>% summarise(rich = length(!is.na(which(reads !=  0))))
eight_long_rare_rich <- eight_long_rare %>% group_by(sample, group) %>% summarise(rich = length(!is.na(which(reads !=  0))))

eight_long_rich$NCDepth <- metadata$NCDepth[match(eight_long_rich$sample, metadata$Sample.Name)]
eight_long_rare_rich$NCDepth <- metadata$NCDepth[match(eight_long_rare_rich$sample, paste0("X",metadata$Sample.Name))]

eight_long_rich$seq <- metadata$sequencing[match(eight_long_rich$sample, metadata$Sample.Name)]
eight_long_rare_rich$seq <- metadata$sequencing[match(eight_long_rare_rich$sample, paste0("X",metadata$Sample.Name))]

eight_long_rich$group <- factor(eight_long_rich$group,
                                levels = c("Chlorophytes", "Cryptophytes",
                                           "Diatoms", "Dinoflagellates",
                                           "Haptophytes", "Metazoa", "Syndiniales", "Other"))

eight_long_rare_rich$group <- factor(eight_long_rare_rich$group,
                                     levels = c("Chlorophytes", "Cryptophytes",
                                                "Diatoms", "Dinoflagellates",
                                                "Haptophytes", "Metazoa", "Syndiniales", "Other"))


no_rare <- ggplot(eight_long_rich, aes(x = NCDepth, y = rich, color = seq)) +
  geom_point() + facet_wrap(~group, scales = "free_y") +
  geom_smooth(method = "lm") + ggtitle("No Rarefaction") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        legend.position=c(.85,.15),
        legend.key = element_blank()) +
  labs(x = "Nitracline Depth (m)", y = "Richness", color = "")


rare <- ggplot(eight_long_rare_rich, aes(x = NCDepth, y = rich, color = seq)) +
  geom_point() + facet_wrap(~group, scales = "free_y") +
  geom_smooth(method = "lm") + ggtitle("Rareified to 17,000 reads") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        legend.position=c(.85,.15),
        legend.key = element_blank()) +
  labs(x = "Nitracline Depth (m)", y = "Richness", color = "")


agg_png("figures/richness_18Sv9_groups.png", width = 12, height = 10, units = "in", res = 400)
plot(no_rare)
dev.off()

agg_png("figures/richness_18Sv9_groups_rare.png", width = 12, height = 10, units = "in", res = 400)
plot(rare)
dev.off()








