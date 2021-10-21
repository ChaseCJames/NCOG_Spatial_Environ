library(tidyverse)
library(vegan)
library(ragg)

load("data/18sv9_all.Rdata")
metadata <- read.csv("data/NCOG_sample_log_DNA_stvx_meta_2014-2020.csv")

samples1 <- which(duplicated(gsub("_S", "", rownames(eighteen_s)), fromLast = TRUE) == TRUE)
samples2 <- which(duplicated(gsub("_S", "", rownames(eighteen_s)), fromLast = FALSE) == TRUE)

compare <- eighteen_s[c(samples1,samples2),]

rare_rich <- rarefy(compare,5000)

r_sums <- rowSums(compare)

for (i in 1:nrow(compare)) {
  print(i)
  compare[i,] <- compare[i,]/r_sums[i]
}

compare$filter <- rep(c("GF_F", "Sterivex"), each = 312)
compare$sample <- gsub("_S", "", rownames(compare))

compare_long <- compare %>% pivot_longer(-c(filter,sample), names_to = "Hash", values_to = "rel_abun")

compare$filter <- NULL
compare$sample <- NULL

bray <- vegdist(compare)

rich_compare <- compare_long %>% group_by(sample, filter) %>% summarise(richness = length(which(rel_abun != 0)))

rich_compare <- rich_compare %>% pivot_wider(names_from = filter, values_from = richness)

rich_plot <- ggplot(rich_compare, aes(x = GF_F, y = Sterivex)) +
  geom_abline(slope = 1, intercept = 0, color = "red", lty = 2) +
  geom_point() +
  coord_fixed(xlim = c(0,2050), ylim = c(0,2050)) +
  labs(x = "GF/F filter Sample Richness", y = "Sterivex Filter Richness") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"))

agg_png("figures/sterivex_v_gf_f.png", width = 6, height = 6, units = "in", res = 400)
plot(rich_plot)
dev.off()

rare_df <- as.data.frame(rare_rich)

rare_df$filter <- rep(c("GF_F", "Sterivex"), each = 312)
rare_df$sample <- gsub("_S", "", rownames(rare_df))

rare_long <- rare_df %>% pivot_longer(-c(filter,sample), values_to = "richness")
rare_long$name <- NULL

rare_long <- rare_long %>% pivot_wider(names_from = filter, values_from = richness)

rare_rich <- ggplot(rare_long, aes(x = GF_F, y = Sterivex)) +
  geom_abline(slope = 1, intercept = 0, color = "red", lty = 2) +
  geom_point() +
  coord_fixed(xlim = c(0, 1100), ylim = c(0,1100)) +
  labs(x = "GF/F filter Rarified Richness", y = "Sterivex Filter Rarified Richness") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"))


agg_png("figures/sterivex_v_gf_f_rarefy.png", width = 6, height = 6, units = "in", res = 400)
plot(rare_rich)
dev.off()

# by groups 

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

compare_long$group <- split_taxa$group[match(compare_long$Hash, split_taxa$Feature.ID)]

rich_compare_group <- compare_long %>% group_by(sample, filter, group) %>% summarise(richness = length(which(rel_abun != 0)))

rich_compare_group <- rich_compare_group %>% pivot_wider(names_from = filter, values_from = richness)

rich_compare_group$group <- factor(rich_compare_group$group,
                                   levels = c("Chlorophytes", "Cryptophytes",
                                              "Diatoms", "Dinoflagellates",
                                              "Haptophytes", "Metazoa", "Syndiniales", "Other"))


rich_compare_group$nitracline <- metadata$NCDepth[match(rich_compare_group$sample, paste0("X",metadata$Sample.Name))]
rich_compare_group$coast <- metadata$Distance[match(rich_compare_group$sample, paste0("X",metadata$Sample.Name))]



by_group <- ggplot(rich_compare_group, aes(x = GF_F, y = Sterivex, fill = group)) +
  geom_abline(slope = 1, intercept = 0, color = "red", lty = 2) +
  geom_point(pch = 21, size = 2) +
  facet_wrap(~group, scales = "free") +
  labs(x = "GF/F Filter Richness", y = "Sterivex Filter Richness", fill = "Group") +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black"),
        legend.position=c(.85,.15),
        legend.key = element_blank(),
        strip.background = element_blank()) +
  scale_fill_manual(values = c("#c361ab",
                               "#60a75c",
                               "#8164cb",
                               "#ad963e",
                               "#6496cd",
                               "#cc693c",
                               "#cb5469", "grey60")) 
  
agg_png("figures/sterivex_v_gf_f_by_group.png", width = 10, height = 10, units = "in", res = 400)
plot(by_group)
dev.off()

mean <- compare_long %>% group_by(Hash) %>% summarise(mean = mean(rel_abun, na.rm = TRUE))

compare_filt <- compare_long %>% pivot_wider(names_from = filter, values_from = rel_abun)

compare_filt$residuals <- compare_filt$Sterivex - compare_filt$GF_F

compare_filt$ster_mean <- compare_filt$Sterivex/mean$mean[match(compare_filt$Hash, mean$Hash)]
compare_filt$gf_mean <- compare_filt$GF_F/mean$mean[match(compare_filt$Hash, mean$Hash)]

compare_filt$res_mean <- compare_filt$residuals/mean$mean[match(compare_filt$Hash, mean$Hash)]

compare_filt <- compare_filt %>% filter(Sterivex != 0 & GF_F != 0)


compare_filt$group <- factor(compare_filt$group,
                                   levels = c("Chlorophytes", "Cryptophytes",
                                              "Diatoms", "Dinoflagellates",
                                              "Haptophytes", "Metazoa", "Syndiniales", "Other"))

residuals <- ggplot(compare_filt, aes(x = Sterivex, y = Sterivex-GF_F, color = group)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       axis.text = element_text(size = 10),
                       axis.title = element_text(size = 10),
                       strip.background = element_blank(),
                       legend.position=c(.85,.15)) +
  labs(x = "GF/F Rel. Abun", y = "Sterivex Rel. Abun", color = "Group") +
  scale_y_log10() + scale_x_log10() +
  geom_hline(yintercept = 0, color = "red", lty = 2, lwd = 2) +
  facet_wrap(~group) +
  scale_color_manual(values = c("#c361ab",
                               "#60a75c",
                               "#8164cb",
                               "#ad963e",
                               "#6496cd",
                               "#cc693c",
                               "#cb5469", "grey60")) 

agg_png("figures/sterivex_v_gf_rel_abun_log_group.png", width = 10, height = 10, units = "in", res = 400)
plot(residuals)
dev.off()

ggplot(compare_filt %>% filter(sample == "X201402_086.7_033.0_10"), aes(x = gf_mean, y = ster_mean)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       axis.text = element_text(size = 10),
                       axis.title = element_text(size = 10)) +
  labs(x = "GF/F Rel. Abun/Mean Rel. Abun", y = "Sterivex Rel. Abun/Mean Rel. Abun")

save(compare, compare_long, file = "output/18sv9_filter_compare.Rdata")

