library(tidyverse)
library(lubridate)
library(ragg)
library(scales)

metadata <- read.csv("data/NCOG_sample_log_DNA_stvx_meta_2014-2020.csv")

metadata$Date <- mdy(metadata$Date)

cruise_date <- metadata %>% group_by(Cruise) %>% summarise(mean_date = mean(Date, na.rm = TRUE))

filter_val <- 0.1

# 16S

six_s$sample <- rownames(six_s)

six_piv <- six_s %>% pivot_longer(-sample, names_to = "hash", values_to = "reads")

six_piv$type <- "Old Seq"

six_piv$type[which(grepl("_S", six_piv$sample, fixed = TRUE) == TRUE)] <- "New Seq"

six_piv$cruise <- metadata$Cruise[match(six_piv$sample, paste0("X",metadata$Sample.Name))]
six_piv <- six_piv %>% group_by(hash,type,cruise) %>% summarise(reads = sum(reads,na.rm = TRUE))

six_test <- six_piv %>% group_by(hash) %>% summarise(zeros  = n_distinct(cruise[which(reads == 0 & type == "New Seq")])/
                                                       n_distinct(cruise[which(type == "New Seq")]),
                                                     zeros_old  = n_distinct(cruise[which(reads == 0 & type == "Old Seq")])/
                                                       n_distinct(cruise[which(type == "Old Seq")]))  

six_test$name <- six_tax_id$silva_Taxon[match(six_test$hash, six_tax_id$Feature.ID)]

oddities <- six_test %>% filter(zeros == 1 & zeros_old < filter_val | zeros_old == 1 & zeros < filter_val )

six_cruise <- six_piv %>% filter(hash %in% c(oddities$hash)) %>%
  group_by(hash,cruise) %>%
  summarise(sum = sum(reads,na.rm=TRUE),
            samps = n_distinct(cruise)) 


six_total <- six_piv %>% group_by(cruise) %>% summarise(total_reads = sum(reads,na.rm = TRUE))

six_cruise$rel_abun <- six_cruise$sum/six_total$total_reads[match(six_cruise$cruise, six_total$cruise)]

six_cruise$mean_date <- cruise_date$mean_date[match(six_cruise$cruise,cruise_date$Cruise)]

six_cruise$name <- six_tax_id$silva_Taxon[match(six_cruise$hash, six_tax_id$Feature.ID)]

tax <- strsplit(six_cruise$name, ";")

taxa_deep <- sapply(tax, "[", 4)
taxa_deep <- sub(paste0(" o__"), "", taxa_deep) 

six_cruise$name_short <- taxa_deep

long_names <- unique(six_cruise$hash)

six_cruise$name_short <- gsub("_", " ", six_cruise$name_short)

ggplot(six_cruise, aes(x = mean_date, y = rel_abun, color = hash)) +
  geom_point(show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~name_short, scales = "free_y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank()) +
  labs(x = "Date", y = "Relative Abundance") + scale_y_continuous(labels = scales::comma_format())


# 18S

eight_s$sample <- rownames(eight_s)

eight_piv <- eight_s %>% pivot_longer(-sample, names_to = "hash", values_to = "reads")

eight_piv$type <- "Old Seq"

eight_piv$type[which(grepl("_S", eight_piv$sample, fixed = TRUE) == TRUE)] <- "New Seq"

eight_piv$cruise <- metadata$Cruise[match(eight_piv$sample, paste0("X",metadata$Sample.Name))]
eight_piv <- eight_piv %>% group_by(hash,type,cruise) %>% summarise(reads = sum(reads,na.rm = TRUE))

eight_test <- eight_piv %>% group_by(hash) %>% summarise(zeros  = n_distinct(cruise[which(reads == 0 & type == "New Seq")])/
                                                       n_distinct(cruise[which(type == "New Seq")]),
                                                     zeros_old  = n_distinct(cruise[which(reads == 0 & type == "Old Seq")])/
                                                       n_distinct(cruise[which(type == "Old Seq")]))  

eight_test$name <- eight_tax_id$pr2_Taxon[match(eight_test$hash, eight_tax_id$Feature.ID)]

oddities <- eight_test %>% filter(zeros == 1 & zeros_old < filter_val | zeros_old == 1 & zeros < filter_val )

eight_cruise <- eight_piv %>% filter(hash %in% c(oddities$hash)) %>%
  group_by(hash,cruise) %>%
  summarise(sum = sum(reads,na.rm=TRUE),
            samps = n_distinct(cruise)) 


eight_total <- eight_piv %>% group_by(cruise) %>% summarise(total_reads = sum(reads,na.rm = TRUE))

eight_cruise$rel_abun <- eight_cruise$sum/eight_total$total_reads[match(eight_cruise$cruise, eight_total$cruise)]

eight_cruise$mean_date <- cruise_date$mean_date[match(eight_cruise$cruise,cruise_date$Cruise)]

eight_cruise$name <- eight_tax_id$pr2_Taxon[match(eight_cruise$hash, eight_tax_id$Feature.ID)]

tax <- strsplit(eight_cruise$name, ";")

taxa_deep <- sapply(tax, "[", 4)

eight_cruise$name_short <- taxa_deep

long_names <- unique(eight_cruise$hash)

eight_cruise$name_short <- gsub("_", " ", eight_cruise$name_short)

ggplot(eight_cruise, aes(x = mean_date, y = rel_abun, color = hash)) +
  geom_point(show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~name_short, scales = "free_y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank()) +
  labs(x = "Date", y = "Relative Abundance")


six_cruise$name_short <- paste(six_cruise$name_short, "16S")
eight_cruise$name_short <- paste(eight_cruise$name_short, "18Sv9")
comb_cruise <- bind_rows(six_cruise, eight_cruise)

comb_cruise$name_short <- factor(comb_cruise$name_short,
                                 levels = c(unique(six_cruise$name_short),unique(eight_cruise$name_short)))

filter <- ggplot(comb_cruise, aes(x = mean_date, y = rel_abun, color = hash)) +
  geom_point(show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~name_short, scales = "free_y") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 11)) +
  labs(x = "Date", y = "Relative Abundance") +
  scale_y_continuous(labels = scales::comma_format()) 

# pdf("figures_S/suppl_fig_filter.pdf", width = 14, height = 6)
# plot(filter)
# dev.off()

six_perc <- sum(six_cruise$sum)/sum(six_piv$reads)
eight_perc <- sum(eight_cruise$sum)/sum(eight_piv$reads)

n_odd_six <- length(unique(six_cruise$hash))
n_odd_eight <- length(unique(eight_cruise$hash))

perc_asv_six <- n_odd_six/length(unique(six_piv$hash))
perc_asv_eight <- n_odd_eight/length(unique(eight_piv$hash))
