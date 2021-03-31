library(tidyverse)
library(ggplotify)

myCol <- brewer.pal(3, "Pastel2")

# load data
load("data/18sv9_all.Rdata")
load("data/18sv9_tara_oceans.Rdata")
load("data/18sv9_tara_polar.Rdata")

calcofi <- names(which(colSums(eighteen_s) != 0))
tara <- names(which(colSums(tara_dat) != 0))
polar <- names(which(colSums(polar_dat) != 0))

# tax

eight_tax_id <- eight_tax_id[,1:3]
colnames(eight_tax_id) <- c("ASV", "Taxon", "Confidence")

colnames(polar_tax) <- c("ASV", "Taxon", "Confidence")
colnames(tara_tax) <- c("ASV", "Taxon", "Confidence")

tax_full <- bind_rows(eight_tax_id, polar_tax, tara_tax)

tax_full$CalCOFI <- 0
tax_full$Tara <- 0
tax_full$Polar <- 0

tax_full <- tax_full[!duplicated(tax_full$ASV),]

tax_full$CalCOFI[which(!is.na(match(tax_full$ASV,calcofi)))] <- 1
tax_full$Tara[which(!is.na(match(tax_full$ASV,tara)))] <- 1
tax_full$Polar[which(!is.na(match(tax_full$ASV,polar)))] <- 1

split_taxa <- separate(tax_full, Taxon, sep = ";", into = c("A","B","C", "D", "E", "F", "G", "H", "I"))

split_taxa$endemic <- NA

split_taxa$endemic[which(split_taxa$CalCOFI == 1 & split_taxa$Tara == 0 & split_taxa$Polar == 0)] <- "NCOG"
split_taxa$endemic[which(split_taxa$CalCOFI == 0 & split_taxa$Tara == 1 & split_taxa$Polar == 0)] <- "TARA"
split_taxa$endemic[which(split_taxa$CalCOFI == 0 & split_taxa$Tara == 0 & split_taxa$Polar == 1)] <- "TARA Polar"

mean_groups <- split_taxa %>%
  filter(!is.na(endemic), !is.na(C)) %>%
  group_by(endemic,C) %>% summarise(count = n()) %>% mutate(prop = count/sum(count))

proportion_plot <- ggplot(mean_groups, aes(x = endemic, y = prop, fill = C)) + 
  geom_bar(stat = "identity") + theme_classic() +
  labs(x = "", y = "Proportion", fill = "Taxonomy") +
  scale_fill_manual(values = c("#dc7e83","#41c65d","#9f5dd2","#53a024","#c74eaf","#81c954","#5b6ed9",
                               "#aebe32","#8358a0","#3d933f","#da4e8d","#51c586","#cf3e5e","#54c1ae",
                               "#d34239","#4ab4dd","#df7032","#6584c7","#db9835","#cf8fd0","#809f39",
                               "#9c4969","#6ea66a","#ac5030","#30866c","#dc996b","#34753e","#986039",
                               "#52701e","#c3ab43","#576426","#abb871","#906d24","#7d8028","#807a45")) +
  theme(panel.border = element_rect(fill = NA, color = "black"))

pdf("figures/proportion_plot_TARA.pdf", width = 10, height = 8)
print(proportion_plot)
dev.off()





