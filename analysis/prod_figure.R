library(tidyverse)
library(scales)
library(viridis)
library(patchwork)
library(standardize)
library(mgcv)

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9")

in_group_names = c("a Prochlorococcus", "b Synecococcus", "c Flavobacteriales",
                   "d Rhodobacterales", "e Sar Clade", "f Diatoms",
                   "g Dinoflagellates", "h Syndiniales", "i Haptophytes",
                   "j Chlorophytes","k Metazoans")


tsize = 12

out_plot_list <- list()

for (i in 1:length(in_group_list)) {
  
  load(paste0("output/", in_group_list[i],"_full_data.RData"))
  
  prodo_dat <- full_dat %>%
    filter(!is.na(IntC14))
  
  prodo_dat$phase <- substr(prodo_dat$Cruise,1,4)
  
  prodo_dat$phase[which(prodo_dat$phase == "2014" | prodo_dat$phase == "2015" | prodo_dat$phase == "2016")] <- "Warm"
  prodo_dat$phase[which(prodo_dat$phase == "2017" | prodo_dat$phase == "2018")] <- "Cool"
  prodo_dat$phase[which(prodo_dat$phase == "2019")] <- "2019"
  
  prodo_dat$phase <- as.factor(prodo_dat$phase)
  prodo_dat$phase <- factor(prodo_dat$phase, levels = c("Warm", "Cool", "2019"))
  
  gam_dat <- prodo_dat %>%
    filter(phase != "2019")
  
  out_plot_list[[i]] <- ggplot() +
    geom_point(data = prodo_dat, aes(x = IntC14*2, y = richness, color = phase)) + scale_x_log10() + 
    stat_smooth(data = gam_dat, aes(x = IntC14*2, y = richness, color = phase, fill = phase),
                method = "gam", se = 0.99) + labs(color = "Phase") +
    scale_color_manual(values = c("red", "blue","gold3")) +
    scale_fill_manual(values = c("red", "blue","gold3")) + ylab("Richness") +
    xlab("mgC/m3 per one light day") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.text = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          plot.title = element_text(size = tsize)) +
    ggtitle(in_group_names[i]) + guides(fill = FALSE)

  
}


out_plot <- (out_plot_list[[1]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                                        axis.ticks.x = element_blank())) +
  (out_plot_list[[2]] + theme(axis.text.x = element_blank(), axis.title = element_blank(),
                              axis.ticks.x = element_blank())) + 
  (out_plot_list[[3]] + theme(axis.text.x = element_blank(), axis.title = element_blank(),
                              axis.ticks.x = element_blank())) +
  (out_plot_list[[4]] + theme(axis.text.x = element_blank(), axis.title = element_blank(),
                              axis.ticks.x = element_blank())) +
  (out_plot_list[[5]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                              axis.ticks.x = element_blank())) +
  (out_plot_list[[6]] + theme(axis.text.x = element_blank(), axis.title = element_blank(),
                              axis.ticks.x = element_blank())) +
  (out_plot_list[[7]] + theme(axis.text.x = element_blank(), axis.title = element_blank(),
                              axis.ticks.x = element_blank())) +
  (out_plot_list[[8]] + theme(axis.title.y = element_blank())) +
  (out_plot_list[[9]] + theme()) +
  (out_plot_list[[10]] + theme(axis.title.y = element_blank())) +
  (out_plot_list[[11]] + theme(axis.title.y = element_blank())) +
     guide_area()  + plot_layout(ncol = 4, guides = "collect")

pdf("figures/prod_vs_div.pdf", width = 12, height = 9)
out_plot
dev.off()













load("output/cyano_16s_full_data.Rdata")

prodo_dat <- full_dat %>%
  filter(!is.na(IntC14))

p1 <- ggplot(prodo_dat, aes(x = IntC14*2, y = richness)) +
  stat_bin_hex(bins = 20) + scale_x_log10() + scale_fill_viridis(limits = c(0, 10), oob = scales::squish) +
  stat_smooth(method = "gam", color = "red", fill = "red", se = 0.99) +
  labs(fill = "# of Samples") + ylab("Richness") +
  xlab("mgC/m3 per one light day") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.text = element_text(size = tsize),
        legend.title = element_text(size = tsize),
        axis.text = element_text(size = tsize),
        axis.title = element_text(size = tsize),
        plot.title = element_text(size = tsize)) +
  ggtitle("Synecococcus + Prochlorococcus Richness vs Productivity")

load("output/chloro_18sv9_full_data.Rdata")

prodo_dat$cyano_rich <- full_dat$richness[match(prodo_dat$Sample.Name, full_dat$Sample.Name)]

prodo_dat$total_rich <- prodo_dat$cyano_rich + prodo_dat$richness

prodo_dat <- prodo_dat[complete.cases(prodo_dat[,c(47,65)]),]

p2 <- ggplot(prodo_dat, aes(x = IntC14*2, y = total_rich)) +
  stat_bin_hex(bins = 20) + scale_x_log10() + scale_fill_viridis(limits = c(0, 10), oob = scales::squish) +
  stat_smooth(method = "gam", color = "red", fill = "red", se = 0.99) +
  labs(fill = "# of Samples") + ylab("Richness") +
  xlab("mgC/m3 per one light day") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.text = element_text(size = tsize),
        legend.title = element_text(size = tsize),
        axis.text = element_text(size = tsize),
        axis.title = element_text(size = tsize),
        plot.title = element_text(size = tsize)) +
  ggtitle("Diatom + Chlorophyte Richness vs Productivity")

pdf(file = "figures/cyano_div_prod.pdf", width = 7.5, height = 5)
print(p1)
dev.off()

#### standardized by proportion of total asvs

load("data/18sv9_all.Rdata")

# eighteen_s[eighteen_s > 0] <- 1

load("data/16s_all.Rdata")

# sixteen_s[sixteen_s > 0] <- 1

eight_sums <- rowSums(eighteen_s)
six_sums <- rowSums(sixteen_s)


type_list = c("16s", "16s","16s", "16s", "16s", 
                  "18sv9","18sv9", "18sv9",
                  "18sv9", "18sv9", "18sv9")

in_group_names2 = c("Prochlorococcus", "Synecococcus", "Flavobacteriales",
                   "Rhodobacterales", "Sar Clade", "Diatoms",
                   "Dinoflagellates", "Syndiniales", "Haptophytes",
                   "Chlorophytes","Metazoans")

in_group_list_basic = c("16s_pro", "16s_syne","16s_flavo", "16s_rhodo", "16s_sar", 
                        "18s_diatom","18s_dino", "18s_syndin",
                        "18s_hapto", "18s_chloro", "18s_metazoa")

out_plot_list <- list()

for (i in 1:length(in_group_list)) {
  
  load(paste0("output/", in_group_list[i],"_full_data.RData"))
  load(paste0("data/",in_group_list_basic[i],".Rdata"))
  
  prodo_dat <- full_dat %>%
    filter(!is.na(IntC14))

  if(type_list[[i]] == "16s"){
    prodo_dat$tot_reads <- six_sums[match(paste0("X",prodo_dat$Sample.Name),names(six_sums))]
    group_tab <- sixteen_s[,which(!is.na(match(colnames(sixteen_s), colnames(asv_table))))]
  }
  
  if(type_list[[i]] == "18sv9"){
    prodo_dat$tot_reads <- eight_sums[match(paste0("X",prodo_dat$Sample.Name),names(eight_sums))]
    group_tab <- eighteen_s[,which(!is.na(match(colnames(eighteen_s), colnames(asv_table))))]
  }
  
  group_sum <- rowSums(group_tab)
  prodo_dat$group_reads <- group_sum[match(paste0("X",prodo_dat$Sample.Name),names(group_sum))]
  
  prodo_dat$ratio <- prodo_dat$group_reads/prodo_dat$tot_reads
  
  prodo_dat <- prodo_dat[complete.cases(prodo_dat[,c(63,66)]),]
  
  m1 <- gam(richness ~ ratio, data = prodo_dat)
  
  prodo_dat$residual_rich <- m1$residuals
  
  p2 <- ggplot(prodo_dat, aes(x = ratio, y = richness)) +
    geom_point() + 
    stat_smooth(method = "gam", fill = "red", color = "red", se = 0.99) +
    labs(fill = "# of Samples") + ylab("Richness") +
    xlab(paste0(in_group_names2[i]," Reads / Total ", type_list[i], " Reads")) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.text = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          plot.title = element_text(size = tsize)) +
    ggtitle(in_group_names[i])
  
  print(p2)

  p3 <- ggplot(prodo_dat, aes(x = IntC14*2, y = residual_rich)) +
    stat_bin_hex(bins = 20) + scale_x_log10() + scale_fill_viridis(limits = c(0, 10), oob = scales::squish) +
    stat_smooth(method = "gam", color = "red", fill = "red", se = 0.99) +
    labs(fill = "# of Samples") + ylab("Residual Richness") +
    xlab("mgC/m3 per one light day") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.text = element_text(size = tsize),
          legend.title = element_text(size = tsize),
          axis.text = element_text(size = tsize),
          axis.title = element_text(size = tsize),
          plot.title = element_text(size = tsize))

  out_plot_list[[i]] <- list(plot1 = p2, plot2 = p3)
  
}


out_plot <- out_plot_list[[1]]$plot1 + out_plot_list[[1]]$plot2 +
  out_plot_list[[2]]$plot1 + out_plot_list[[2]]$plot2 +
  out_plot_list[[3]]$plot1 + out_plot_list[[3]]$plot2 +
  out_plot_list[[4]]$plot1 + out_plot_list[[4]]$plot2 +
  out_plot_list[[5]]$plot1 + out_plot_list[[5]]$plot2 +
  out_plot_list[[6]]$plot1 + out_plot_list[[6]]$plot2 +
  out_plot_list[[7]]$plot1 + out_plot_list[[7]]$plot2 +
  out_plot_list[[8]]$plot1 + out_plot_list[[8]]$plot2 +
  out_plot_list[[9]]$plot1 + out_plot_list[[9]]$plot2 +
  out_plot_list[[10]]$plot1 + out_plot_list[[10]]$plot2 +
  out_plot_list[[11]]$plot1 + out_plot_list[[11]]$plot2 +
  plot_layout(guides= "collect", ncol = 2)

pdf("figures/adj_prod_div.pdf", width = 15, height = 40)
out_plot
dev.off()


load(paste0("output/", in_group_list[1],"_full_data.RData"))

# six and eight

prodo_dat <- full_dat %>%
  filter(!is.na(IntC14))

prodo_dat$six_rich <- six_sums[match(paste0("X",prodo_dat$Sample.Name),names(six_sums))]
prodo_dat$eight_rich <- six_sums[match(paste0("X",prodo_dat$Sample.Name),names(eight_sums))]

six_plot <- ggplot(prodo_dat, aes(x = IntC14*2, y = six_rich)) +
  stat_bin_hex(bins = 20) + scale_x_log10() + scale_fill_viridis(limits = c(0, 10), oob = scales::squish) +
  stat_smooth(method = "gam", color = "red", fill = "red", se = 0.99) +
  labs(fill = "# of Samples") + ylab("16s Richness") +
  xlab("mgC/m3 per one light day") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.text = element_text(size = tsize),
        legend.title = element_text(size = tsize),
        axis.text = element_text(size = tsize),
        axis.title = element_text(size = tsize),
        plot.title = element_text(size = tsize)) +
  ggtitle("A. 16s Richness vs Productivity")

eight_plot <- ggplot(prodo_dat, aes(x = IntC14*2, y = eight_rich)) +
  stat_bin_hex(bins = 20) + scale_x_log10() + scale_fill_viridis(limits = c(0, 10), oob = scales::squish) +
  stat_smooth(method = "lm", color = "red", fill = "red", se = 0.99) +
  labs(fill = "# of Samples") + ylab("18sv9 Richness") +
  xlab("mgC/m3 per one light day") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.text = element_text(size = tsize),
        legend.title = element_text(size = tsize),
        axis.text = element_text(size = tsize),
        axis.title = element_text(size = tsize),
        plot.title = element_text(size = tsize))+
  ggtitle("B. 18sv9 Richness vs Productivity")

pdf(file = "figures/total_pdr.pdf", width = 12, height = 6)
print(six_plot + eight_plot + plot_layout(guides = "collect"))
dev.off()
