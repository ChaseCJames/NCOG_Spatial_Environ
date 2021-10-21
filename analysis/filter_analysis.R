library(tidyverse)
library(lubridate)
library(vegan)
library(patchwork)
library(mgcv)
library(ragg)

load("output/total_dissimilar.Rdata")
load("output/total_full_data.Rdata")

dissimilar[upper.tri(dissimilar,diag = TRUE)] <- NA

dissimilar <- as.data.frame(dissimilar)

dissimilar$name <- rownames(dissimilar)

diss_long <- dissimilar %>% pivot_longer(-name, names_to = "name_2", values_to = "bc_dist")

# convert to similarity

diss_long$bc_dist <- 1 - diss_long$bc_dist

diss_long$name <- gsub("X", "", diss_long$name)
diss_long$name_2 <- gsub("X", "", diss_long$name_2)

# collect metadata 

meta_dist <- full_dat[,c(1,34,35,40,38,39,
                         42,49)]

meta_dist[,2:8] <- apply(meta_dist[,2:8], 2, function(x) (x-min(x,na.rm = TRUE))/(max(x,na.rm = TRUE)-min(x, na.rm = TRUE)))

meta_dist <- meta_dist[complete.cases(meta_dist),]

diss_long <- diss_long %>% filter(name %in% meta_dist$Sample.Name, name_2 %in% meta_dist$Sample.Name)

# filter id

full_dat$filter <- "Sterivex-GP"
full_dat$filter[which(full_dat$Date < dmy("01-01-2017"))] <- "GF/F"

diss_long$name_filt <- full_dat$filter[match(diss_long$name, full_dat$Sample.Name)]
diss_long$name_2_filt <- full_dat$filter[match(diss_long$name_2, full_dat$Sample.Name)]

# environmental distance

e_dist <- dist(meta_dist[,2:8], diag = NA, upper = TRUE) %>% as.matrix()

e_dist[upper.tri(e_dist,diag = TRUE)] <- NA

e_dist <- as.data.frame(e_dist)

colnames(e_dist) <- meta_dist$Sample.Name

e_dist$name <- meta_dist$Sample.Name

e_long <- e_dist %>% pivot_longer(-name, names_to = "name_2", values_to = "e_dist")

long_similarity <- full_join(diss_long, e_long, by = c("name", "name_2"))

long_similarity <- long_similarity[complete.cases(long_similarity),]

long_similarity$filter <- "GF/F"

long_similarity$filter[which(long_similarity$name_filt == long_similarity$name_2_filt & long_similarity$name_filt == "Sterivex-GP")] <- "Sterivex-GP"
long_similarity$filter[which(long_similarity$name_filt != long_similarity$name_2_filt)] <- "Different Filters"

long_similarity$shannon_diff <- abs(full_dat$shannon[match(long_similarity$name, full_dat$Sample.Name)] -
      full_dat$shannon[match(long_similarity$name_2, full_dat$Sample.Name)])

long_similarity$rich_diff <- abs(full_dat$richness[match(long_similarity$name, full_dat$Sample.Name)] -
                                      full_dat$richness[match(long_similarity$name_2, full_dat$Sample.Name)])


long_similarity$time_diff <- abs(difftime(full_dat$Date[match(long_similarity$name, full_dat$Sample.Name)],
                                   full_dat$Date[match(long_similarity$name_2, full_dat$Sample.Name)],units = "days")) %>% as.numeric()

dist <- ggplot(long_similarity, aes(x = e_dist, y = bc_dist, color = filter , fill = filter)) +
  stat_smooth(method = "lm", se = 0.99) +
  labs(x = "Environmental Eucledian Distance",
       y = "Bray Curtis Similarity", color = "", fill = "") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.key = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 12)) +
  scale_color_manual(values = c("#7aa457","#a46cb7","#cb6a49")) +
  scale_fill_manual(values = c("#7aa457","#a46cb7","#cb6a49"))

shan <- ggplot(long_similarity, aes(x = e_dist, y = shannon_diff, color = filter, fill = filter)) +
  stat_smooth(method = "lm", se = 0.99) +
  labs(x = "Environmental Eucledian Distance",
       y = "Difference in Shannon Index", color = "", fill = "") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.key = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 12)) +
  scale_color_manual(values = c("#7aa457","#a46cb7","#cb6a49")) +
  scale_fill_manual(values = c("#7aa457","#a46cb7","#cb6a49"))

rich <- ggplot(long_similarity, aes(x = e_dist, y = rich_diff, color = filter , fill = filter)) +
  stat_smooth(method = "lm", se = 0.99) + 
  labs(x = "Environmental Eucledian Distance",
       y = "Difference in Richness", color = "", fill = "") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.key = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 12)) +
  scale_color_manual(values = c("#7aa457","#a46cb7","#cb6a49")) +
  scale_fill_manual(values = c("#7aa457","#a46cb7","#cb6a49"))

long_bc <- long_similarity %>% filter(bc_dist > quantile(long_similarity$bc_dist, probs = 0.05),
                                      bc_dist < quantile(long_similarity$bc_dist, probs = 0.95),
                                      e_dist > quantile(long_similarity$e_dist, probs = 0.05),
                                      e_dist < quantile(long_similarity$e_dist, probs = 0.95))

long_sha <- long_similarity %>% filter(shannon_diff > quantile(long_similarity$shannon_diff, probs = 0.05),
                                       shannon_diff < quantile(long_similarity$shannon_diff, probs = 0.95),
                                       e_dist > quantile(long_similarity$e_dist, probs = 0.05),
                                       e_dist < quantile(long_similarity$e_dist, probs = 0.95))

long_rich <- long_similarity %>% filter(rich_diff > quantile(long_similarity$rich_diff, probs = 0.05),
                                        rich_diff < quantile(long_similarity$rich_diff, probs = 0.95),
                                        e_dist > quantile(long_similarity$e_dist, probs = 0.05),
                                        e_dist < quantile(long_similarity$e_dist, probs = 0.95))

# bray curtis

long_bc <- long_bc %>% mutate(bc_anom = bc_dist - mean(bc_dist))
long_bc <- long_bc %>% group_by(filter) %>% mutate(filt_bc_mean = mean(bc_dist)) %>% mutate(filt_bc_anom = bc_dist - filt_bc_mean)

out_1 <- gam(bc_dist ~ e_dist, data = long_bc) 
out_3 <- gam(bc_dist ~ e_dist + filter, data = long_bc) 

res_bc <- out_1$residuals %>% as.data.frame()
colnames(res_bc) <- "Environment"
res_bc$Filter <- long_bc$filt_bc_anom
res_bc$Both <- out_3$residuals
res_bc$Null <- long_bc$bc_anom

res_bc$type <- long_bc$filter

res_bc <- res_bc %>% pivot_longer(-type, names_to = "Model", values_to = "Residuals")
res_bc$Model[which(res_bc$Model == "Both")] <- "Environment + Filter"

res_bc$Model <- factor(res_bc$Model,
                       levels = c("Null", "Filter", "Environment", "Environment + Filter"))

# Shannon

long_sha <- long_sha %>% mutate(sha_anom = shannon_diff - mean(shannon_diff))
long_sha <- long_sha %>% group_by(filter) %>% mutate(filt_sha_mean = mean(shannon_diff)) %>% mutate(filt_sha_anom = shannon_diff - filt_sha_mean)

out_1 <- gam(shannon_diff ~ e_dist, data = long_sha) 
out_3 <- gam(shannon_diff ~ e_dist + filter, data = long_sha) 

res_sha <- out_1$residuals %>% as.data.frame()
colnames(res_sha) <- "Environment"
res_sha$Filter <- long_sha$filt_sha_anom
res_sha$Both <- out_3$residuals
res_sha$Null <- long_sha$sha_anom

res_sha <- res_sha %>% pivot_longer(everything(), names_to = "Model", values_to = "Residuals")
res_sha$Model[which(res_sha$Model == "Both")] <- "Environment + Filter"

res_sha$Model <- factor(res_sha$Model,
                        levels = c("Null", "Filter", "Environment", "Environment + Filter"))

# Richness

long_rich <- long_rich %>% mutate(rich_anom = rich_diff - mean(rich_diff))
long_rich <- long_rich %>% group_by(filter) %>% mutate(filt_rich_mean = mean(rich_diff)) %>% mutate(filt_rich_anom = rich_diff - filt_rich_mean)

out_1 <- gam(rich_diff ~ e_dist, data = long_rich) 
out_3 <- gam(rich_diff ~ e_dist + filter, data = long_rich) 

res_rich <- out_1$residuals %>% as.data.frame()
colnames(res_rich) <- "Environment"
res_rich$Filter <- long_rich$filt_rich_anom
res_rich$Both <- out_3$residuals
res_rich$Null <- long_rich$rich_anom

res_rich <- res_rich %>% pivot_longer(everything(), names_to = "Model", values_to = "Residuals")
res_rich$Model[which(res_rich$Model == "Both")] <- "Environment + Filter"

res_rich$Model <- factor(res_rich$Model,
                         levels = c("Null", "Filter", "Environment", "Environment + Filter"))


ggplot(res_bc, aes(x = Model, y = abs(Residuals))) +
  geom_violin(draw_quantiles = 0.5)

ggplot(res_sha, aes(x = Model, y = abs(Residuals))) +
  geom_violin(draw_quantiles = 0.5)

ggplot(res_rich, aes(x = Model, y = abs(Residuals))) +
  geom_violin(draw_quantiles = 0.5)



#### Specific Groups ####



dist + shan + rich + plot_layout(guides = "collect")

file_names <- c("pro_16s", "syne_16s", "flavo_16s", "rhodo_16s",
                "sar_16s", "diatom_18sv9", "dino_18sv9",
                "syndin_18sv9", "hapto_18sv9", "chloro_18sv9",
                "metazoa_18sv9")

plot_names <- c("Prochlorococcus", "Synechococcus", "Flavobacteriales", "Rhodobacterales",
                "SAR 11", "Diatoms", "Dinoflagellates",
                "Syndiniales", "Haptophytes", "Chlorophytes",
                "Metazoans")

for(i in 1:length(file_names)){
  
  print(i)
  
  load(paste0("output/",file_names[i],"_dissimilar.Rdata"))
  load(paste0("output/",file_names[i],"_full_data.Rdata"))
  
  dissimilar[upper.tri(dissimilar,diag = TRUE)] <- NA
  
  dissimilar <- as.data.frame(dissimilar)
  
  dissimilar$name <- rownames(dissimilar)
  
  diss_long <- dissimilar %>% pivot_longer(-name, names_to = "name_2", values_to = "bc_dist")
  
  # convert to similarity
  
  diss_long$bc_dist <- 1 - diss_long$bc_dist
  
  diss_long$name <- gsub("X", "", diss_long$name)
  diss_long$name_2 <- gsub("X", "", diss_long$name_2)
  
  # collect metadata 
  
  meta_dist <- full_dat[,c(1,34,35,40,38,39,
                           42,49)]
  
  meta_dist[,2:8] <- apply(meta_dist[,2:8], 2, function(x) (x-min(x,na.rm = TRUE))/(max(x,na.rm = TRUE)-min(x, na.rm = TRUE)))
  
  meta_dist <- meta_dist[complete.cases(meta_dist),]
  
  diss_long <- diss_long %>% filter(name %in% meta_dist$Sample.Name, name_2 %in% meta_dist$Sample.Name)
  
  # filter id
  
  full_dat$filter <- "Sterivex-GP"
  full_dat$filter[which(full_dat$Date < dmy("01-01-2017"))] <- "GF/F"
  
  diss_long$name_filt <- full_dat$filter[match(diss_long$name, full_dat$Sample.Name)]
  diss_long$name_2_filt <- full_dat$filter[match(diss_long$name_2, full_dat$Sample.Name)]
  
  # environmental distance
  
  e_dist <- dist(meta_dist[,2:8], diag = NA, upper = TRUE) %>% as.matrix()
  
  e_dist[upper.tri(e_dist,diag = TRUE)] <- NA
  
  e_dist <- as.data.frame(e_dist)
  
  colnames(e_dist) <- meta_dist$Sample.Name
  
  e_dist$name <- meta_dist$Sample.Name
  
  e_long <- e_dist %>% pivot_longer(-name, names_to = "name_2", values_to = "e_dist")
  
  long_similarity <- full_join(diss_long, e_long, by = c("name", "name_2"))
  
  long_similarity <- long_similarity[complete.cases(long_similarity),]
  
  long_similarity$filter <- "GF/F"
  
  long_similarity$filter[which(long_similarity$name_filt == long_similarity$name_2_filt & long_similarity$name_filt == "Sterivex-GP")] <- "Sterivex-GP"
  long_similarity$filter[which(long_similarity$name_filt != long_similarity$name_2_filt)] <- "Different Filters"
  
  long_similarity$shannon_diff <- abs(full_dat$shannon[match(long_similarity$name, full_dat$Sample.Name)] -
                                        full_dat$shannon[match(long_similarity$name_2, full_dat$Sample.Name)])
  
  long_similarity$rich_diff <- abs(full_dat$richness[match(long_similarity$name, full_dat$Sample.Name)] -
                                     full_dat$richness[match(long_similarity$name_2, full_dat$Sample.Name)])
  
  long_bc <- long_similarity %>% filter(bc_dist > quantile(long_similarity$bc_dist, probs = 0.05),
                                        bc_dist < quantile(long_similarity$bc_dist, probs = 0.95),
                                        e_dist > quantile(long_similarity$e_dist, probs = 0.05),
                                        e_dist < quantile(long_similarity$e_dist, probs = 0.95))
  
  long_sha <- long_similarity %>% filter(shannon_diff > quantile(long_similarity$shannon_diff, probs = 0.05),
                                         shannon_diff < quantile(long_similarity$shannon_diff, probs = 0.95),
                                         e_dist > quantile(long_similarity$e_dist, probs = 0.05),
                                         e_dist < quantile(long_similarity$e_dist, probs = 0.95))
  
  long_rich <- long_similarity %>% filter(rich_diff > quantile(long_similarity$rich_diff, probs = 0.05),
                                          rich_diff < quantile(long_similarity$rich_diff, probs = 0.95),
                                          e_dist > quantile(long_similarity$e_dist, probs = 0.05),
                                          e_dist < quantile(long_similarity$e_dist, probs = 0.95))
  
# bray curtis
  
  long_bc <- long_bc %>% mutate(bc_anom = bc_dist - mean(bc_dist))
  long_bc <- long_bc %>% group_by(filter) %>% mutate(filt_bc_mean = mean(bc_dist)) %>% mutate(filt_bc_anom = bc_dist - filt_bc_mean)

  out_1 <- glm(bc_dist ~ e_dist, data = long_bc) 
  out_3 <- glm(bc_dist ~ e_dist + filter, data = long_bc) 
  
  res_bc <- out_1$residuals %>% as.data.frame()
  colnames(res_bc) <- "Environment"
  res_bc$Filter <- long_bc$filt_bc_anom
  res_bc$Both <- out_3$residuals
  res_bc$Null <- long_bc$bc_anom
  
  res_bc <- res_bc %>% pivot_longer(everything(), names_to = "Model", values_to = "Residuals")
  res_bc$Model[which(res_bc$Model == "Both")] <- "Environment + Filter"
  
  res_bc$Model <- factor(res_bc$Model,
                      levels = c("Null", "Filter", "Environment", "Environment + Filter"))
  
  # Shannon
  
  long_sha <- long_sha %>% mutate(sha_anom = shannon_diff - mean(shannon_diff))
  long_sha <- long_sha %>% group_by(filter) %>% mutate(filt_sha_mean = mean(shannon_diff)) %>% mutate(filt_sha_anom = shannon_diff - filt_sha_mean)
  
  out_1 <- glm(shannon_diff ~ e_dist, data = long_sha) 
  out_3 <- glm(shannon_diff ~ e_dist + filter, data = long_sha) 
  
  res_sha <- out_1$residuals %>% as.data.frame()
  colnames(res_sha) <- "Environment"
  res_sha$Filter <- long_sha$filt_sha_anom
  res_sha$Both <- out_3$residuals
  res_sha$Null <- long_sha$sha_anom
  
  res_sha <- res_sha %>% pivot_longer(everything(), names_to = "Model", values_to = "Residuals")
  res_sha$Model[which(res_sha$Model == "Both")] <- "Environment + Filter"
  
  res_sha$Model <- factor(res_sha$Model,
                      levels = c("Null", "Filter", "Environment", "Environment + Filter"))
  
  # Richness
  
  long_rich <- long_rich %>% mutate(rich_anom = rich_diff - mean(rich_diff))
  long_rich <- long_rich %>% group_by(filter) %>% mutate(filt_rich_mean = mean(rich_diff)) %>% mutate(filt_rich_anom = rich_diff - filt_rich_mean)
  
  out_1 <- glm(rich_diff ~ e_dist, data = long_rich) 
  out_3 <- glm(rich_diff ~ e_dist + filter, data = long_rich) 
  
  res_rich <- out_1$residuals %>% as.data.frame()
  colnames(res_rich) <- "Environment"
  res_rich$Filter <- long_rich$filt_rich_anom
  res_rich$Both <- out_3$residuals
  res_rich$Null <- long_rich$rich_anom
  
  res_rich <- res_rich %>% pivot_longer(everything(), names_to = "Model", values_to = "Residuals")
  res_rich$Model[which(res_rich$Model == "Both")] <- "Environment + Filter"
  
  res_rich$Model <- factor(res_rich$Model,
                          levels = c("Null", "Filter", "Environment", "Environment + Filter"))
  

  ggplot(res_bc, aes(x = Model, y = abs(Residuals))) +
    geom_violin()
  
  ggplot(res_sha, aes(x = Model, y = abs(Residuals))) +
    geom_violin()
  
  ggplot(res_rich, aes(x = Model, y = abs(Residuals))) +
    geom_violin()
  
  
  dist <- ggplot(long_bc, aes(x = e_dist, y = bc_dist, color = filter , fill = filter)) +
    stat_smooth(method = "lm", se = 0.99) +
    labs(x = "Environmental Eucledian Distance",
         y = "Bray Curtis Similarity", color = "", fill = "") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.key = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 12)) +
    scale_color_manual(values = c("#7aa457","#a46cb7","#cb6a49")) +
    scale_fill_manual(values = c("#7aa457","#a46cb7","#cb6a49"))
  
  shan <- ggplot(long_sha, aes(x = e_dist, y = shannon_diff, color = filter, fill = filter)) +
    stat_smooth(method = "lm", se = 0.99) +
    labs(x = "Environmental Eucledian Distance",
         y = "Difference in Shannon Index", color = "", fill = "") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.key = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 12)) +
    scale_color_manual(values = c("#7aa457","#a46cb7","#cb6a49")) +
    scale_fill_manual(values = c("#7aa457","#a46cb7","#cb6a49"))
  
  rich <- ggplot(long_rich, aes(x = e_dist, y = rich_diff, color = filter , fill = filter)) +
    stat_smooth(method = "lm", se = 0.99) + 
    labs(x = "Environmental Eucledian Distance",
         y = "Difference in Richness", color = "", fill = "") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          legend.key = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = 12)) +
    scale_color_manual(values = c("#7aa457","#a46cb7","#cb6a49")) +
    scale_fill_manual(values = c("#7aa457","#a46cb7","#cb6a49"))
  
  dist + shan + rich + plot_layout(guides = "collect")
    
  
  
}





