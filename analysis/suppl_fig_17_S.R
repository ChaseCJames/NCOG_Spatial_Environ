library(tidyverse)
library(patchwork)
library(data.table)

# read in 16s (2014-2019)

sixteen_s <- read.csv("data/NCOG_21_16S_redo2_asv_count_tax.csv", header = TRUE)

six_id_names <- sixteen_s$Feature.ID

six_tax_id <- sixteen_s[,c(1,(ncol(sixteen_s)-1):ncol(sixteen_s))]

sixteen_s <- sixteen_s[,-c(1,(ncol(sixteen_s)-1):ncol(sixteen_s))]

sixteen_s <- t(sixteen_s)

sixteen_s <- as.data.frame(sixteen_s)

six_tp <- rownames(sixteen_s)

colnames(sixteen_s) <- six_id_names

mock_vals <- rownames(sixteen_s) %>% tolower() %>%
  grepl("mock",.) %>% which(. == TRUE)

six_mock <- sixteen_s[mock_vals,]

# eighteen s

eighteen_s <- read.csv("data/NCOG_18sV9_asv_count_tax_S.csv", header = TRUE)

# pull out mock samples

eight_id_names <- eighteen_s$Feature.ID

eight_tax_id <- eighteen_s[,c(1,(ncol(eighteen_s)-5):ncol(eighteen_s))]

eighteen_s <- eighteen_s[,-c(1,(ncol(eighteen_s)-5):ncol(eighteen_s))]

eighteen_s <- t(eighteen_s)

eighteen_s <- as.data.frame(eighteen_s)

eight_tp <- rownames(eighteen_s)

colnames(eighteen_s) <- eight_id_names

mock_vals <- rownames(eighteen_s) %>% tolower() %>%
  grepl("mock",.) %>% which(. == TRUE)

eight_mock <- eighteen_s[mock_vals,]

#### sort mock communities

six_mock <- six_mock[,-which(colSums(six_mock) == 0)]
eight_mock <- eight_mock[,-which(colSums(eight_mock) == 0)]

# taxonomy for ASVs

colnames(six_mock) <- six_tax_id$silva_Taxon[match(colnames(six_mock), six_tax_id$Feature.ID)]

colnames(eight_mock) <- eight_tax_id$silva_Taxon[match(colnames(eight_mock), eight_tax_id$Feature.ID)]

six_mock$sample <- rownames(six_mock)
eight_mock$sample <- rownames(eight_mock)

# rename samples and make factor

six_long <- six_mock %>% pivot_longer(-sample, names_to = "ASV", values_to = "read") %>%
  group_by(sample) %>% mutate(rel_abun = read/sum(read, na.rm = TRUE))

eight_long <- eight_mock %>% pivot_longer(-sample, names_to = "ASV", values_to = "read") %>%
  group_by(sample) %>% mutate(rel_abun = read/sum(read, na.rm = TRUE))

eight_long <- eight_long %>% filter(sample != "X2017_mock_even2")

six_long$Hash <-  six_tax_id$Feature.ID[match(six_long$ASV,six_tax_id$silva_Taxon)]
eight_long$Hash <-  eight_tax_id$Feature.ID[match(eight_long$ASV,eight_tax_id$silva_Taxon)]

# filter

adj_six_long <- six_long %>% filter(rel_abun > 0.001) %>%
  group_by(sample) %>% mutate(rel_abun_adj = read/sum(read, na.rm = TRUE))

adj_eight_long <- eight_long %>% filter(rel_abun > 0.0001) %>%
  group_by(sample) %>% mutate(rel_abun_adj = read/sum(read, na.rm = TRUE))

#expected samples
adj_six_long$ASV <- gsub("d__|p__|c__|o__|f__|g__|s__", "", adj_six_long$ASV)
adj_eight_long$ASV <- gsub("d__|p__|c__|o__|f__|g__|s__", "", adj_eight_long$ASV)

adj_six_long <- adj_six_long %>% group_by(sample,Hash) %>% summarise(rel_abun_adj = sum(rel_abun_adj, na.rm = TRUE), ASV = ASV)
adj_eight_long <- adj_eight_long %>% group_by(sample,Hash) %>% summarise(rel_abun_adj = sum(rel_abun_adj, na.rm = TRUE), ASV = ASV)


# remove duplicates

adj_six_long <- adj_six_long[!duplicated(adj_six_long),]
adj_eight_long <- adj_eight_long[!duplicated(adj_eight_long),]
adj_eight_long$PR2 <- eight_tax_id$pr2_Taxon[match(adj_eight_long$Hash, eight_tax_id$Feature.ID)]

# Remove space near semi-colon, add space instead of underscore

adj_six_long$ASV <- gsub("; ",";", adj_six_long$ASV)
adj_six_long$ASV <- gsub("_"," ", adj_six_long$ASV)

# filter to correct mocks
mocks_16 <- unique(adj_six_long$sample)

adj_six_long <- adj_six_long %>% filter(sample %in% mocks_16[c(1:10,17:34)])

# pivot wide

six_pivot <- adj_six_long[,c(1,3,2)] %>% pivot_wider(names_from = Hash, values_from = rel_abun_adj) %>% t()
eight_pivot <- adj_eight_long[,c(1,2,3)] %>% pivot_wider(names_from = c(Hash), values_from = rel_abun_adj) %>% t()

eight_pivot <- as.data.frame(eight_pivot)
six_pivot <- as.data.frame(six_pivot)
six_pivot$ASV <- adj_six_long$ASV[match(rownames(six_pivot), adj_six_long$Hash)]
eight_pivot$PR2_Taxon <- adj_eight_long$PR2[match(rownames(eight_pivot), adj_eight_long$Hash)]
eight_pivot$Silva_Taxon <- adj_eight_long$ASV[match(rownames(eight_pivot), adj_eight_long$Hash)]

# rearrange

six_rank <- adj_six_long %>% group_by(Hash) %>% summarise(mean = mean(rel_abun_adj)) %>% arrange(desc(mean))
eight_rank <- adj_eight_long %>% group_by(Hash) %>% summarise(mean = mean(rel_abun_adj)) %>% arrange(desc(mean))

six_pivot <- six_pivot[c(1,((match(rownames(six_pivot), six_rank$Hash))[2:nrow(six_pivot)]+1)),]
eight_pivot <- eight_pivot[c(1,((match(rownames(eight_pivot), eight_rank$Hash))[2:nrow(eight_pivot)]+1)),]

six_pivot <- six_pivot[,c(29,1:10,21:22,26:27,11:19,23:25,28)]
eight_pivot <- eight_pivot[,c(20:21,1:2,6:8,12:13,16:17,3:5,9:11,14:15,18:19)]

write.csv(six_pivot, file = "output/16S_mock_community_low_filter.csv")
write.csv(eight_pivot, file = "output/18S_mock_community_low_filter.csv")

# expected distributions
#16s
expected_six <- read.csv("data/16S_mock_S.csv", stringsAsFactors = FALSE)
expected_six <- expected_six[,c(2:5)]

colnames(expected_six)[1] <- "ASV" 
expected_six <- expected_six %>%
  pivot_longer(-c(ASV,taxa), names_to = "sample", values_to = "rel_abun_adj")

adj_six_long$taxa <- expected_six$taxa[match(adj_six_long$ASV,expected_six$ASV)]

expected_six$ASV <- NA
expected_six$Hash <- NA

expected_six <- expected_six %>% distinct(taxa,sample, .keep_all = TRUE)

adj_six_long <- bind_rows(adj_six_long, expected_six)

# adj_six_long <- adj_six_long[complete.cases(adj_six_long),]
#18sv9

expected_eight <- read.csv("data/18S_mock_S.csv", stringsAsFactors = FALSE)
colnames(expected_eight)[1:6] <- c("Hash", "PR2", "ASV",  "taxa", "stag_expected", "even_expected")

expected_eight <- expected_eight[,c(1:6)]

adj_eight_long$taxa <- expected_eight$taxa[match(adj_eight_long$ASV,expected_eight$ASV)]

expected_eight$ASV <- NA
expected_eight$Hash <- NA
expected_eight$PR2 <- NA
expected_eight <- expected_eight %>%
  pivot_longer(-c(ASV,taxa,PR2, Hash), names_to = "sample", values_to = "rel_abun_adj")

expected_eight <- expected_eight %>% distinct(taxa,sample, .keep_all = TRUE)

adj_eight_long <- bind_rows(adj_eight_long, expected_eight)

# plots

six_samps <- unique(adj_six_long$sample)
eight_samps <- unique(adj_eight_long$sample)

adj_six_long$sample <- factor(adj_six_long$sample, levels = six_samps[c(30,7:10,1:6,21:22,25:26,29,17:20,
                                                                        11:16,23:24,27:28)])

adj_eight_long$sample <- factor(adj_eight_long$sample, levels = eight_samps[c(21,1:2,6:8,12:13,16:17,
                                                                              20,3:5,9:11,14:15,18:19)])

adj_six_long$rel_abun_adj[which(is.na(adj_six_long$rel_abun_adj))] <- 0

six_asv_order <- adj_six_long %>% group_by(taxa) %>%
  summarise(mean_rel = mean(rel_abun_adj, na.rm = TRUE)) %>%
  arrange(desc(mean_rel))

adj_six_long$taxa <- factor(adj_six_long$taxa, levels = c(six_asv_order$taxa))

eight_asv_order <- adj_eight_long %>% group_by(taxa) %>%
  summarise(mean_rel = mean(rel_abun_adj, na.rm = TRUE)) %>%
  arrange(desc(mean_rel))

adj_eight_long$taxa <- factor(adj_eight_long$taxa, levels = c(eight_asv_order$taxa))

six <- ggplot(adj_six_long, aes(x = sample , y = rel_abun_adj, fill = taxa)) + 
  geom_bar(stat = "identity", show.legend = TRUE) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = -45, hjust = 0),
        legend.position = "bottom",
        plot.margin = ggplot2::margin(5,60,5,5),
        axis.text = element_text(size = 11),
        axis.title.y = element_text(size = 11)) + 
  guides(fill = guide_legend(ncol = 6)) +
  labs(x = "", y = "Relative Abundance", title = "16S V4-V5 Mock Communities", fill = "Taxonomy") + 
  scale_fill_manual(values = c("#55c1a1","#a056c9","#57b64d","#e26fd0","#a7b639",
                               "#6364cf","#d5a23c","#5b7ac1","#d74c2b","#4cadd8",
                               "#d1752c","#b392df","#627527","#b53994","#3e8757",
                               "#dd4080","#a3b167","#965089","#927232","#dc86b3",
                               "#a15a38","#de4757","#e49375","#b05167","#b24439", "grey80"))

eight <- ggplot(adj_eight_long, aes(x = sample , y = rel_abun_adj, fill = taxa)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = -45,hjust = 0),
        legend.position = "bottom",
        plot.margin = ggplot2::margin(5,60,5,5),
        axis.text = element_text(size = 11),
        axis.title.y = element_text(size = 11)) + 
  guides(fill = guide_legend(ncol = 5)) +
  labs(x = "", y = "Relative Abundance", title = "18sv9 Mock Communities", fill = "Taxonomy") + 
  scale_fill_manual(values = c("#acb039","#a859c9","#5db645","#5e6ed9","#de8f34",
                               "#5f97d3","#d35636","#4cb9b0","#d03c56","#5eb779",
                               "#d04798","#577936","#cb8ccf","#c4a15f","#785a9c",
                               "#9e6431","grey80")) 


out <- six / eight

pdf(file = "figures_S/supp_fig_17_S.pdf", width = 10, height = 10)
print(out)
dev.off()

# 
# pdf(file = "figures_S/supp_fig_16s_S.pdf", width = 30, height = 15)
# print(six)
# dev.off()









