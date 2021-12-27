library(tidyverse)
library(patchwork)
library(data.table)

# read in 16s (2014-2019)

sixteen_s <- read.csv("data/NCOG_16S_asv_count_tax_S.csv", header = TRUE)

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

adj_six_long <- six_long %>% filter(rel_abun > 0.0005) %>%
  group_by(sample) %>% mutate(rel_abun_adj = read/sum(read, na.rm = TRUE))

adj_eight_long <- eight_long %>% filter(rel_abun > 0.0005) %>%
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

adj_six_long <- adj_six_long %>% filter(sample %in% mocks_16[c(1:7,14:28)])

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
expected_six <- read.csv("data/16S_expected_S.csv", stringsAsFactors = FALSE)

expected_six <- expected_six[,c(1,3:4)]
colnames(expected_six)[1] <- "ASV" 
expected_six <- expected_six %>%
  pivot_longer(-ASV, names_to = "sample", values_to = "rel_abun_adj")
expected_six$Hash <- adj_six_long$Hash[match(expected_six$ASV, adj_six_long$ASV)]
expected_six$Hash[which(is.na(expected_six$Hash))] <- "No Associated Hash"

adj_six_long <- bind_rows(adj_six_long, expected_six)

adj_six_long <- adj_six_long[complete.cases(adj_six_long),]
#18sv9

expected_eight <- read.csv("data/18S_expected_S.csv", stringsAsFactors = FALSE)
colnames(expected_eight)[1:6] <- c("Hash", "ASV", "PR2",  "taxa", "stag_expected", "even_expected")

expected_eight <- expected_eight[,c(1:6)]

expected_eight$taxa <- NULL
expected_eight <- expected_eight %>%
  pivot_longer(-c(ASV,PR2, Hash), names_to = "sample", values_to = "rel_abun_adj")

adj_eight_long <- bind_rows(adj_eight_long, expected_eight)

adj_eight_long <- adj_eight_long[complete.cases(adj_eight_long),]

# plots

six_samps <- unique(adj_six_long$sample)
eight_samps <- unique(adj_eight_long$sample)

adj_six_long$sample <- factor(adj_six_long$sample, levels = six_samps[c(24,1:6,15:16,19:20,7,23,
                                                                        8:13,17:18,21:22,14)])

adj_eight_long$sample <- factor(adj_eight_long$sample, levels = eight_samps[c(21,1:2,6:8,12:13,16:17,
                                                                              20,3:5,9:11,14:15,18:19)])

six_asv_order <- adj_six_long %>% group_by(ASV) %>%
  summarise(mean_rel = mean(rel_abun_adj, na.rm = TRUE)) %>%
  arrange(desc(mean_rel))

adj_six_long$ASV <- factor(adj_six_long$ASV, levels = c(six_asv_order$ASV))

eight_asv_order <- adj_eight_long %>% group_by(ASV) %>%
  summarise(mean_rel = mean(rel_abun_adj, na.rm = TRUE)) %>%
  arrange(desc(mean_rel))

adj_eight_long$ASV <- factor(adj_eight_long$ASV, levels = c(eight_asv_order$ASV))

six <- ggplot(adj_six_long, aes(x = sample , y = rel_abun_adj, fill = ASV)) + 
  geom_bar(stat = "identity", show.legend = TRUE) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = -45),
        legend.position = "bottom") + 
  guides(fill = guide_legend(ncol = 2)) +
  labs(x = "", y = "Relative Abundance", title = "16S V4-V5 Mock Communities", fill = "Taxonomy") + 
  scale_fill_manual(values = c("#d97dd8","#60bf47","#9a5ad0","#acb836","#5b70dd",
                               "#cea339","#c042a4","#46c787","#d74989","#3a8c3b",
                               "#d63d56","#43c4c6","#d74a2c","#56a6db","#e08930",
                               "#4973b6","#688f2c","#7956a3","#80bc6c","#9a4972",
                               "#59af8a","#b65831","#a498e0","#73741d","#cf85b8",
                               "#32784e","#e38087","#6a7139","#a44750","#b0af66",
                               "#93662d","#de9a6c"))

eight <- ggplot(adj_eight_long, aes(x = sample , y = rel_abun_adj, fill = ASV)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = -45),
        legend.position = "bottom") + 
  guides(fill = guide_legend(ncol = 2)) +
  labs(x = "", y = "Relative Abundance", title = "18sv9 Mock Communities", fill = "Taxonomy") + 
  scale_fill_manual(values = c("#95528d","#62bf4b","#bf57c5","#b3b538","#7465ce",
                               "#4e8e36","#d4478f","#58c289","#d33e55","#49bfcb",
                               "#d0502b","#6587ca","#dc9030","#da8fca","#61702a",
                               "#ad5062","#388661","#e28471","#9eb46b","#9c5b2d",
                               "#d5a768","#937b2b")) 


out <- six / eight

pdf(file = "figures_S/supp_fig_16_S.pdf", width = 30, height = 30)
print(out)
dev.off()

# 
# pdf(file = "figures/figure_outline/supp_fig_16_six.pdf", width = 30, height = 15)
# print(six)
# dev.off()









