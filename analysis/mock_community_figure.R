library(tidyverse)
library(patchwork)
library(stringdist)
library(data.table)

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

# filter out rel abun < 0.5%

adj_six_long <- six_long %>% filter(rel_abun > 0.00025) %>%
  group_by(sample) %>% mutate(rel_abun_adj = read/sum(read, na.rm = TRUE))

adj_eight_long <- eight_long %>% filter(rel_abun > 0.00025) %>%
  group_by(sample) %>% mutate(rel_abun_adj = read/sum(read, na.rm = TRUE))

#expected samples
adj_six_long$ASV <- gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__|D_7__|D_8__", "", adj_six_long$ASV)
adj_eight_long$ASV <- gsub("D_0__|D_1__|D_2__|D_3__|D_4__|D_5__|D_6__|D_7__|D_8__", "", adj_eight_long$ASV)

adj_six_long <- adj_six_long %>% group_by(sample,ASV) %>%summarise(rel_abun_adj = sum(rel_abun_adj, na.rm = TRUE))
adj_eight_long <- adj_eight_long %>% group_by(sample,ASV) %>%summarise(rel_abun_adj = sum(rel_abun_adj, na.rm = TRUE))

# make expected distributions
six_even_exp <- adj_six_long %>% filter(sample == "even_mock_a1", rel_abun_adj > 0.01)
six_even_exp$rel_abun_adj <- 1/11

six_stag_exp <- adj_six_long %>% filter(sample == "stag_mock_a1")
six_stag_exp <- six_stag_exp[order(six_stag_exp$rel_abun_adj, decreasing = TRUE),]
six_even_exp$rel_abun_adj <- c(0.315,0.09,0.09,0.158,0.068,
                               0.009,0.018,0.018,0.018, 0.045,
                               0.009)


expected_six <- expected_six %>% 
  pivot_longer(-ASV, names_to = "sample", values_to = "rel_abun_adj")

expected_eight <- read.csv("data/eight_expected.csv")
expected_eight$ASV <- gsub("; ", ";", expected_eight$ASV)

match(substr(expected_eight$ASV,1,12), substr(adj_eight_long$ASV,1,12))


six_samps <- unique(adj_six_long$sample)
eight_samps <- unique(adj_eight_long$sample)

adj_six_long$sample <- as.factor(adj_six_long$sample)
adj_six_long$sample <- factor(adj_six_long$sample, levels = six_samps[c(1:6,13:15,19:20,7:12,16:18,21:22)])

adj_eight_long$sample <- as.factor(adj_eight_long$sample)
adj_eight_long$sample <- factor(adj_eight_long$sample, levels = eight_samps[c(1:2,6:8,12:13,3:5,9:11,14:15)])


six <- ggplot(adj_six_long, aes(x = sample , y = rel_abun_adj, fill = ASV)) + 
  geom_bar(stat = "identity", show.legend = TRUE) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = -45),
        legend.position = "bottom") + 
  guides(fill = guide_legend(ncol = 2)) +
  labs(x = "", y = "Relative Abundance", title = "16sv4 Mock Communities", fill = "Taxonomy") + 
  scale_fill_manual(values = c("#dd86b6",
                               "#4bbc50",
                               "#9958cb",
                               "#bcb72f",
                               "#5c6bd4",
                               "#80af3f",
                               "#d157ba",
                               "#428542",
                               "#d9417f",
                               "#53c5a1",
                               "#cf3a45",
                               "#4eb0db",
                               "#dd5c31",
                               "#5a78bc",
                               "#da9334",
                               "#b18dd9",
                               "#b9a64e",
                               "#974c87",
                               "#93b46f",
                               "#a14658",
                               "#308b75",
                               "#dd7979",
                               "#796d1d",
                               "#a4542a",
                               "#726f39",
                               "#d49465"))

eight <- ggplot(adj_eight_long, aes(x = sample , y = rel_abun_adj, fill = ASV)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text.x = element_text(angle = -45),
        legend.position = "bottom") + 
  guides(fill = guide_legend(ncol = 2)) +
  labs(x = "", y = "Relative Abundance", title = "18sv9 Mock Communities", fill = "Taxonomy") + 
  scale_fill_manual(values = c("#dd86b6",
                               "#4bbc50",
                               "#9958cb",
                               "#bcb72f",
                               "#5c6bd4",
                               "#80af3f",
                               "#d157ba",
                               "#428542",
                               "#d9417f",
                               "#53c5a1",
                               "#cf3a45",
                               "#4eb0db",
                               "#dd5c31",
                               "#5a78bc",
                               "#da9334",
                               "#b18dd9",
                               "#b9a64e",
                               "#974c87",
                               "#93b46f",
                               "#a14658",
                               "#308b75",
                               "#dd7979",
                               "#796d1d",
                               "#a4542a",
                               "#726f39",
                               "#d49465")) 


out <- six / eight

pdf(file = "figures/figure_outline/supp_fig_15.pdf", width = 30, height = 30)
print(out)
dev.off()





