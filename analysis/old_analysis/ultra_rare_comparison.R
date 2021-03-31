library(tidyverse)
library(patchwork)

load("data/18sv9_all.Rdata")
load("data/18sv9_tara_oceans.Rdata")
load("data/18sv9_tara_polar.Rdata")

# sp
asv_table <- eighteen_s

splits <- strsplit(rownames(asv_table), "_")

stations <- vector()

for (i in 1:length(splits)) {
  stations[i] <- paste0(splits[[i]][2]," ",splits[[i]][3])
}

asv_table <- asv_table[which(as.numeric(substr(stations,1,3)) > 75),]
stations <- stations[which(as.numeric(substr(stations,1,3)) > 75)]

asv_table <- asv_table[,-which(colSums(asv_table, na.rm = TRUE) == 0)]

rank_abun <- asv_table
rank_abun[rank_abun > 0] <- 1
ultra_rare <- names(which(colSums(rank_abun) == 1))

tara_binary <- as.matrix(tara_dat)
polar_binary <- as.matrix(polar_dat)

tara_binary[tara_binary > 0] <- 1
polar_binary[polar_binary > 0] <- 1

tara_rel <- as.matrix(tara_dat)
polar_rel <- as.matrix(polar_dat)

for (i in 1:nrow(tara_rel)) {
  if(i/10 == round(i/10)){print(i)}
  tara_rel[i,] <- tara_rel[i,]/sum(tara_rel[i,], na.rm = TRUE)
}

for (i in 1:nrow(polar_rel)) {
  if(i/10 == round(i/10)){print(i)}
  polar_rel[i,] <- polar_rel[i,]/sum(polar_rel[i,], na.rm = TRUE)
}

tara_summary <- as.data.frame(matrix(NA,nrow = ncol(tara_dat), ncol =5))
polar_summary <- as.data.frame(matrix(NA,nrow = ncol(polar_dat), ncol = 5))

colnames(tara_summary) <- c("ASV", "Mean_Rel", "Sum_Read", "Occurance", "NCOG")
colnames(polar_summary) <- c("ASV", "Mean_Rel", "Sum_Read", "Occurance", "NCOG")

# TARA

NCOG <- colnames(asv_table)

tara_summary$ASV <- colnames(tara_dat)
tara_summary$Mean_Rel <- colMeans(tara_rel, na.rm = TRUE)
tara_summary$Occurance <- colSums(tara_binary)
tara_summary$NCOG <- "Everything Else in NCOG"
tara_summary$Sum_Read <- colSums(tara_dat)

tara_summary <- tara_summary[tara_summary$ASV %in% NCOG,]
tara_summary$NCOG_Occur <- colSums(rank_abun[,match(tara_summary$ASV, NCOG)], na.rm = TRUE)
tara_summary$NCOG_Sum_Read <- colSums(asv_table[,match(tara_summary$ASV, NCOG)], na.rm = TRUE)

polar_summary$ASV <- colnames(polar_dat)
polar_summary$Mean_Rel <- colMeans(polar_rel, na.rm = TRUE)
polar_summary$Occurance <- colSums(polar_binary)
polar_summary$NCOG <- "Everything Else in NCOG"
polar_summary$Sum_Read <- colSums(polar_dat)

polar_summary <- polar_summary[polar_summary$ASV %in% NCOG,]
polar_summary$NCOG_Occur <- colSums(rank_abun[,match(polar_summary$ASV, NCOG)], na.rm = TRUE)
polar_summary$NCOG_Sum_Read <- colSums(asv_table[,match(polar_summary$ASV, NCOG)], na.rm = TRUE)

tara_summary$NCOG[which(!is.na(match(tara_summary$ASV, ultra_rare)))] <- "Ultra-Rare in NCOG"
polar_summary$NCOG[which(!is.na(match(polar_summary$ASV, ultra_rare)))] <- "Ultra-Rare in NCOG"


a <- ggplot(tara_summary, aes(x = NCOG, y = Occurance, fill = NCOG)) + 
  geom_violin(draw_quantiles = c(0.5), lwd = 1) + scale_y_log10() + theme_bw() +
  xlab("") + labs(fill = "") +
  ggtitle("NCOG ASVs in Tara Oceans")

  
b <- ggplot(tara_summary, aes(x = NCOG, y = Sum_Read, fill = NCOG)) + 
  geom_violin(draw_quantiles = c(0.5), lwd = 1) + theme_bw() + scale_y_log10() +
  xlab("") + labs(fill = "") + ylab("Total Reads") +
  ggtitle("NCOG ASVs in Tara Oceans")

c <- ggplot(polar_summary, aes(x = NCOG, y = Occurance, fill = NCOG)) + 
  geom_violin(draw_quantiles = c(0.5), lwd = 1) + scale_y_log10() + theme_bw() +
  xlab("") + labs(fill = "") + ggtitle("NCOG ASVs in Tara Polar")

d <- ggplot(polar_summary, aes(x = NCOG, y = Sum_Read, fill = NCOG)) + 
  geom_violin(draw_quantiles = c(0.5), lwd = 1) + theme_bw() + scale_y_log10() +
  xlab("") + labs(fill = "") + ylab("Total Reads") + 
  ggtitle("NCOG ASVs in Tara Polar")

# a <- ggplot() +
#   geom_hex(data = tara_summary, aes(x = Occurance, y = Sum_Read)) +
#   scale_y_log10() + scale_fill_gradient(low = "grey90", high = "black") +
#   geom_point(data = tara_summary %>% filter(NCOG == "Ultra-Rare in NCOG"),
#              aes(x = Occurance, y = Sum_Read, color = NCOG)) +
#   theme(panel.background = element_blank(),
#         panel.border = element_rect(fill = NA, color = "black")) +
#   ylab("Total Reads") + labs(color = "", fill = "# of ASVs") +
#   ggtitle("NCOG Ultra-Rare in Tara Oceans")
# 
# b <- ggplot() +
#   geom_hex(data = polar_summary, aes(x = Occurance, y = Sum_Read)) +
#   scale_y_log10() + scale_fill_gradient(low = "grey90", high = "black") +
#   geom_point(data = polar_summary %>% filter(NCOG == "Ultra-Rare in NCOG"),
#              aes(x = Occurance, y = Sum_Read, color = NCOG)) +
#   theme(panel.background = element_blank(),
#         panel.border = element_rect(fill = NA, color = "black")) +
#   ylab("Total Reads") + labs(color = "", fill = "# of ASVs") +
#   ggtitle("NCOG Ultra-Rare in Tara Polar")

out_plot <- a + b + c + d + plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect")

pdf(file = "figures/tara_ultra_dominance.pdf", width = 13, height = 10)
print(out_plot)
dev.off()


##### Test similarity of groups

taxonomy <- read.csv("data/18Sv9_taxa_summary.csv")

taxonomy$PR2_Taxon[which(!is.na(match(taxonomy$Feature.ID, ultra_rare)))]



taxa <- taxonomy$PR2_Taxon[match(colnames(asv_table),taxonomy$Feature.ID)]
taxa <- strsplit(taxa, ";")
taxa_sha <- sapply(taxa, "[", 3)
taxa_deep <- sapply(taxa, "[", 4)

taxa_sha[which(taxa_deep == "Mamiellophyceae")] <- "Mamiellophyceae"
taxa_sha[which(taxa_deep == "Bacillariophyta")] <- "Bacillariophyta"
taxa_sha[which(taxa_deep == "Chrysophyceae")] <- "Chrysophyceae"
taxa_sha[which(taxa_deep == "Arthropoda")] <- "Arthropoda"
taxa_sha[which(taxa_deep == "Syndiniales")] <- "Syndiniales"

unique_tax <- unique(taxa_sha)

others <- unique_tax[c(8)]

taxa_sha[which(taxa_sha %in% others)] <- "Other Eukaryotic Phytoplankton"

taxa_sha <- sapply(taxa, "[", 4)

unique_tax <- unique(taxa_sha)

others <- unique_tax[c(4,25)]

taxa_sha[which(taxa_sha %in% others)] <- "Other Eukaryotic Protists"

taxa_sha[which(taxa_sha == "Opisthokonta_X")] <- "Opisthokonta"
taxa_sha[which(taxa_sha == "Alveolata_X")] <- "Alveolata"
taxa_sha[which(taxa_sha == "Stramenopiles_X")] <- "Stramenopiles"

taxa_df <- as.data.frame(matrix(NA,nrow = length(taxa), ncol = 3))
colnames(taxa_df) <- c("ID", "Taxonomy", "Ultra_Rare")

taxa_df$ID <- colnames(asv_table)
taxa_df$Taxonomy <- taxa_sha
taxa_df$Ultra_Rare <- "N"

taxa_df$Ultra_Rare[which(!is.na(match(taxa_df$ID, ultra_rare)))] <- "Y"

taxa_summary <- taxa_df %>%
  group_by(Ultra_Rare, Taxonomy) %>%
  summarise(n = n())

taxa_summary <- taxa_summary[order(taxa_summary$n,decreasing = TRUE),]
taxa_summary$Taxonomy <- as.factor(taxa_summary$Taxonomy)
taxa_summary$Taxonomy <- factor(taxa_summary$Taxonomy, levels = unique(taxa_summary$Taxonomy)[38:1])

colors <- c("#d82774","#ff8199","#6a0020","#ff616e","#590c00",
            "#ff956f","#952900","#e56631","#8e5300","#9e8200",
            "#d9da63","#bbcf40","#4b6c00","#92e973","#93e88b",
            "#017623","#008e4e","#02b97d","#33eebc","#01a8fb",
            "#0186f3","#015bae","#85a9ff","#003472","#033aa8",
            "#9d89d0","#35115b","#d59fff","#812d9f","#ff96fd",
            "#945289","#52004a","#a92191","#ff95d1","#b8107f",
            "#75004c","#78244b","#d07496")

everything_else <- taxa_summary %>% filter(Ultra_Rare == "N")
everything_else$prop <- everything_else$n/sum(everything_else$n, na.rm = TRUE)

u_rare <- taxa_summary %>% filter(Ultra_Rare == "Y")
u_rare$prop <- u_rare$n/sum(u_rare$n, na.rm = TRUE)

ggplot(everything_else, aes(y = Taxonomy, x = prop, color = Taxonomy)) +
  geom_point() +
  geom_segment(aes(y=Taxonomy, 
                   yend=Taxonomy, 
                   x=0, 
                   xend=prop),size = 0.1, show.legend = FALSE) +
  scale_color_manual(values = colors, drop = FALSE) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"))

ggplot(u_rare, aes(y = Taxonomy, x = prop, color = Taxonomy)) +
  geom_point() +
  geom_segment(aes(y=Taxonomy, 
                   yend=Taxonomy, 
                   x=0, 
                   xend=prop),size = 0.1, show.legend = FALSE) +
  scale_color_manual(values = colors, drop = FALSE) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"))







