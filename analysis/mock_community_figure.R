library(tidyverse)
library(patchwork)
library(stringdist)
library(data.table)


# read in expected distributions found in (https://doi.org/10.1101/866731)
# tab_16 <- read_tsv("data/210223_BGT_example_mocks_all-16S-seqs.with-tax.tsv")
# tab_16 <- tab_16[,c(1:3,9)]
# colnames(tab_16) <- c("Hash", "Even", "Staggered", "Taxonomy")
# tab_16 <- tab_16[-which(rowSums(tab_16[,2:3]) == 0),]
# 
# tab_18 <- read_tsv("data/210223_BGT_example_mocks_all-18S-seqs.with-SILVA132-tax.tsv")
# tab_18 <- tab_18[,c(1:3,9)]
# colnames(tab_18) <- c("Hash", "Even", "Staggered", "Taxonomy")
# tab_18 <- tab_18[-which(rowSums(tab_18[,2:3]) == 0),]

# read in 16s (2014-2019)
sixteen_s <- read.csv("data/16_asv_count_tax_final_19.csv", stringsAsFactors = FALSE)

six_id_names <- sixteen_s$Feature.ID

six_tax_id <- sixteen_s[,c(1,(ncol(sixteen_s)-3):ncol(sixteen_s))]

sixteen_s <- apply(sixteen_s, 2, as.numeric)

sixteen_s <- t(sixteen_s)

sixteen_s <- as.data.frame(sixteen_s)

six_tp <- rownames(sixteen_s)

colnames(sixteen_s) <- six_id_names

six_mock <- sixteen_s[c(656:661,665:670,863:868,875:878),]

# eighteen s

eighteen_s <- read.csv("data/18Sv9_asv_count_tax_final_19.csv", stringsAsFactors = FALSE)

eight_id_names <- eighteen_s$Feature.ID

eight_tax_id <- eighteen_s[,c(1,(ncol(eighteen_s)-3):ncol(eighteen_s))]

eighteen_s <- eighteen_s[,-c(1,(ncol(eighteen_s)-3):ncol(eighteen_s))]

eighteen_s <- apply(eighteen_s, 2, as.numeric)

eighteen_s <- t(eighteen_s)

eighteen_s <- as.data.frame(eighteen_s)

eight_tp <- rownames(eighteen_s)

colnames(eighteen_s) <- eight_id_names

eight_mock <- eighteen_s[c(479:484,667:672,871:874),]

#### sort mock communities

six_mock <- six_mock[,-which(colSums(six_mock) == 0)]
eight_mock <- eight_mock[,-which(colSums(eight_mock) == 0)]

# taxonomy for ASVs

colnames(six_mock) <- six_tax_id$Silva_Taxon[match(colnames(six_mock), six_tax_id$Feature.ID)]

colnames(eight_mock) <- eight_tax_id$Silva_Taxon[match(colnames(eight_mock), eight_tax_id$Feature.ID)]

six_mock$sample <- rownames(six_mock)
eight_mock$sample <- rownames(eight_mock)

# rename samples and make factor

six_mock$sample 

eight_mock$sample

six_long <- six_mock %>% pivot_longer(-sample, names_to = "ASV", values_to = "read") %>%
  group_by(sample) %>% mutate(rel_abun = read/sum(read, na.rm = TRUE))

eight_long <- eight_mock %>% pivot_longer(-sample, names_to = "ASV", values_to = "read") %>%
  group_by(sample) %>% mutate(rel_abun = read/sum(read, na.rm = TRUE))

eight_long <- eight_long %>% filter(sample != "X2017_mock_even2")

six_long$Hash <-  six_tax_id$Feature.ID[match(six_long$ASV,six_tax_id$Silva_Taxon)]
eight_long$Hash <-  eight_tax_id$Feature.ID[match(eight_long$ASV,eight_tax_id$Silva_Taxon)]

# filter

adj_six_long <- six_long %>% filter(rel_abun > 0.0005) %>%
  group_by(sample) %>% mutate(rel_abun_adj = read/sum(read, na.rm = TRUE))

adj_eight_long <- eight_long %>% filter(rel_abun > 0.0005) %>%
  group_by(sample) %>% mutate(rel_abun_adj = read/sum(read, na.rm = TRUE))

#expected samples
adj_six_long$ASV <- gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__|D_7__|D_8__", "", adj_six_long$ASV)
adj_eight_long$ASV <- gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__|D_7__|D_8__", "", adj_eight_long$ASV)

adj_six_long <- adj_six_long %>% group_by(sample,Hash) %>% summarise(rel_abun_adj = sum(rel_abun_adj, na.rm = TRUE), ASV = ASV)
adj_eight_long <- adj_eight_long %>% group_by(sample,Hash) %>% summarise(rel_abun_adj = sum(rel_abun_adj, na.rm = TRUE), ASV = ASV)


# remove duplicates

adj_six_long <- adj_six_long[!duplicated(adj_six_long),]
adj_eight_long <- adj_eight_long[!duplicated(adj_eight_long),]
adj_eight_long$PR2 <- eight_tax_id$PR2_Taxon[match(adj_eight_long$Hash, eight_tax_id$Feature.ID)]

# pivot wide

six_pivot <- adj_six_long[,c(1,3,4)] %>% pivot_wider(names_from = ASV, values_from = rel_abun_adj) %>% t()
eight_pivot <- adj_eight_long[,c(1,2,3)] %>% pivot_wider(names_from = c(Hash), values_from = rel_abun_adj) %>% t()

eight_pivot <- as.data.frame(eight_pivot)
eight_pivot$PR2_Taxon <- adj_eight_long$PR2[match(rownames(eight_pivot), adj_eight_long$Hash)]
eight_pivot$Silva_Taxon <- adj_eight_long$ASV[match(rownames(eight_pivot), adj_eight_long$Hash)]

write.csv(six_pivot, file = "output/16S_mock_community_low_filter.csv")
write.csv(eight_pivot, file = "output/18S_mock_community_low_filter.csv")

# expected distributions
#16s
expected_six <- read.csv("data/16S_expected.csv", stringsAsFactors = FALSE)

removes <- c(colnames(expected_six)[21:24], "Mock.even3", "Mock.stag3")

expected_six <- expected_six[,c(1,3:4)]
colnames(expected_six)[1] <- "ASV" 
expected_six <- expected_six %>%
  pivot_longer(-ASV, names_to = "sample", values_to = "rel_abun_adj")
expected_six$Hash <- adj_six_long$Hash[match(expected_six$ASV, adj_six_long$ASV)]
expected_six$Hash[which(is.na(expected_six$Hash))] <- "No Associated Hash"

adj_six_long <- bind_rows(adj_six_long, expected_six)

adj_six_long <- adj_six_long %>% filter(!sample %in% removes)

adj_six_long <- adj_six_long[complete.cases(adj_six_long),]
#18sv9

expected_eight <- read.csv("data/18S_expected.csv", stringsAsFactors = FALSE)
colnames(expected_eight) <- c("Hash", "PR2", "ASV", "taxa", "stag_expected", "even_expected")
expected_eight$taxa <- NULL
expected_eight <- expected_eight %>%
  pivot_longer(-c(ASV,PR2, Hash), names_to = "sample", values_to = "rel_abun_adj")

adj_eight_long <- bind_rows(adj_eight_long, expected_eight)

adj_eight_long <- adj_eight_long %>% filter(!sample %in% removes)

adj_eight_long <- adj_eight_long[complete.cases(adj_eight_long),]

# plots

six_samps <- unique(adj_six_long$sample)
eight_samps <- unique(adj_eight_long$sample)

adj_six_long$sample <- as.factor(adj_six_long$sample)
adj_six_long$sample <- factor(adj_six_long$sample, levels = six_samps[c(18,1:6,13,14,17,7:12,15,16)])

adj_eight_long$sample <- as.factor(adj_eight_long$sample)
adj_eight_long$sample <- factor(adj_eight_long$sample, levels = eight_samps[c(17,1:2,6:8,12:13,16,3:5,9:11,14:15)])

six_asv_order <- adj_six_long %>% group_by(ASV) %>%
  summarise(mean_rel = mean(rel_abun_adj, na.rm = TRUE)) %>%
  arrange(desc(mean_rel))

adj_six_long$ASV <- as.factor(adj_six_long$ASV)
adj_six_long$ASV <- factor(adj_six_long$ASV, levels = c(six_asv_order$ASV))

eight_asv_order <- adj_eight_long %>% group_by(ASV) %>%
  summarise(mean_rel = mean(rel_abun_adj, na.rm = TRUE)) %>%
  arrange(desc(mean_rel))

adj_eight_long$ASV <- as.factor(adj_eight_long$ASV)
adj_eight_long$ASV <- factor(adj_eight_long$ASV, levels = c(eight_asv_order$ASV))

six <- ggplot(adj_six_long, aes(x = sample , y = rel_abun_adj, fill = ASV)) + 
  geom_bar(stat = "identity", show.legend = TRUE) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = -45),
        legend.position = "bottom") + 
  guides(fill = guide_legend(ncol = 2)) +
  labs(x = "", y = "Relative Abundance", title = "16S V4-V5 Mock Communities", fill = "Taxonomy") + 
  scale_fill_manual(values = c("#5d6ad9","#5fbd47","#be6ee1","#a8b635", "#8845ac","#4b8833",
                               "#d24daf","#67bf79","#e33f73","#4dc0b1","#d54e23","#51a4d9",
                               "#d97f2e","#5269ac","#d2a138","#9a93df","#767420","#be4080",
                               "#37835d","#c5373e","#aeaf65","#96568e","#6a7139","#da8bc9",
                               "#975a28","#e17f93","#d49a67","#a14957","#d6664a","#e28274"))

eight <- ggplot(adj_eight_long, aes(x = sample , y = rel_abun_adj, fill = ASV)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = -45),
        legend.position = "bottom") + 
  guides(fill = guide_legend(ncol = 2)) +
  labs(x = "", y = "Relative Abundance", title = "18sv9 Mock Communities", fill = "Taxonomy") + 
  scale_fill_manual(values = c("#ab5061","#54bb4c","#c05ac7","#9cb835","#7065d8","#bfa939",
                               "#5a6ab4","#da9334","#4f9ed5","#d0512c","#3ebabd","#d64459",
                               "#58c89e","#d5438e","#4da664","#954d8d","#5e8a2a","#ab94e0",
                               "#92bc71","#dd86b5","#357b4f","#e38c71","#6d6e28","#9c5d2b",
                               "#bba365")) 


out <- six / eight

pdf(file = "figures/figure_outline/supp_fig_16.pdf", width = 30, height = 30)
print(out)
dev.off()

# 
# pdf(file = "figures/figure_outline/supp_fig_16_six.pdf", width = 30, height = 15)
# print(six)
# dev.off()









