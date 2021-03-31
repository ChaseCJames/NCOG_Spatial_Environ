library(randomForest)

load("output/total_full_data.Rdata")

full_dat$Year <- as.numeric(substr(full_dat$Cruise,1,4))

full_dat <- full_dat %>%
  filter(Depthm < 15)

year_list <- 2014:2019

aic_list <- list()

for (i in 1:length(year_list)) {
  
  partial_dat <- full_dat %>%
    filter(Year == year_list[i])
  
  out <- aic_table_func_diveristy_sign(som_maps = partial_dat, col_n = 63)
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  colnames(out) <- c("Variables","AIC", "Slope")
  out[,2] <- as.numeric(out[,2])
  out[,2] <- round(out[,2], digits = 2)
  out[,3] <- as.numeric(out[,3])
  out[,3] <- round(out[,3], digits = 3)
  
  aic_list[[i]] <- out
  
  
}

AIC_full <- full_join(aic_list[[1]], aic_list[[2]], by = "Variables")

for (i in 3:length(aic_list)) {
  
  AIC_full <- full_join(AIC_full, aic_list[[i]], by = "Variables")
  
}

AIC_full <- AIC_full[,c(1,seq(2,ncol(AIC_full),by = 2),seq(3,ncol(AIC_full),by = 2))]

colnames(AIC_full) <- c("Variables", year_list, paste0(year_list,"_slope"))

AIC_scaled <- AIC_full

for (i in 2:(length(year_list)+1)) {
  zero_one_scale <- 1-(AIC_scaled[,i]-min(AIC_scaled[,i], na.rm = TRUE))/
    abs(min(AIC_scaled[,i], na.rm = TRUE) - max(AIC_scaled[,i], na.rm = TRUE))
  AIC_scaled[,i] <- 50^zero_one_scale
  
}

plot_df <- melt(AIC_scaled[,1:(length(year_list)+1)])
plot_slope <- melt(AIC_scaled[,c(1,(length(year_list)+2):ncol(AIC_scaled))])

plot_df$slope <- plot_slope$value

colnames(plot_df) <- c("Variables", "Group", "AIC","slope")

plot_df$Variables <- as.factor(plot_df$Variables)

plot_df$Variables <- factor(plot_df$Variables, levels = c("Temp", "Salinity",
                                                          "NO3", "PO4", 
                                                          "SiO4","NCD",
                                                          "Chl-a", "O2"))

plot_df$Group <- as.factor(plot_df$Group)

plot_df$Group <- factor(plot_df$Group, levels = c("2014", "2015",
                                                  "2016", "2017", 
                                                  "2018", "2019"))


outplot <- ggplot(data = plot_df, aes(x = Variables, y = Group, size = AIC)) + 
  geom_point(color = "black", alpha = 0.6, shape = 21, fill = "red") +
  labs(fill = "Correlation") + ylab("Variable") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "right",
        panel.grid.major.x = element_line(color = "grey", linetype = 2),
        plot.title = element_text(hjust = 0, size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  scale_size_continuous(range = c(1,18), guide = FALSE) + xlab("") + ylab("") +
  ggtitle("Variable importance for total alpha diversity")

pdf(file = "figures/grant_var_import.pdf", width = 8, height = 6)
print(outplot)
dev.off()




aic_table_func_diveristy_sign <- function(som_maps = partial_dat, col_n = 63){  
  
  som_glm <- som_maps
  
  colnames(som_glm)[col_n] <- "response"
  
  som_glm$response <- as.numeric(som_glm$response)
  
  model_AIC <- matrix(NA,8,3)
  
  # temperature
  
  glm_mean_temp <- glm(response ~  T_degC, data = som_glm)
  mt_sum <- summary(glm_mean_temp)
  model_AIC[1,2] <- mt_sum$aic
  model_AIC[1,1] <- "Temp"
  model_AIC[1,3] <- cor(som_glm$T_degC, som_glm$response, use = "pairwise.complete.obs")
  
  # sea surface temperature
  
  # glm_mean_sst <- glm(response ~  sst_mean, data = som_glm)
  # mt_sum <- summary(glm_mean_sst)
  # model_AIC[2,2] <- mt_sum$aic
  # model_AIC[2,1] <- "Mean SST"
  # model_AIC[2,3] <- cor(som_glm$sst_mean, som_glm$response, use = "pairwise.complete.obs")
  
  # salinity
  
  glm_mean_sal <- glm(response ~  Salnty, data = som_glm)
  mt_sum <- summary(glm_mean_sal)
  model_AIC[2,2] <- mt_sum$aic
  model_AIC[2,1] <- "Salinity"
  model_AIC[2,3] <- cor(som_glm$Salnty, som_glm$response, use = "pairwise.complete.obs")
  
  
  # NO3
  
  glm_mean_no3 <- glm(response ~  NO3ug, data = som_glm)
  mn_sum <- summary(glm_mean_no3)
  model_AIC[3,2] <- mn_sum$aic
  model_AIC[3,1] <- "NO3"
  model_AIC[3,3] <- cor(som_glm$NO3ug, som_glm$response, use = "pairwise.complete.obs")
  
  # PO4
  
  glm_mean_po4 <- glm(response ~  PO4ug, data = som_glm)
  mp_sum <- summary(glm_mean_po4)
  model_AIC[4,2] <- mp_sum$aic
  model_AIC[4,1] <- "PO4"
  model_AIC[4,3] <- cor(som_glm$PO4ug, som_glm$response, use = "pairwise.complete.obs")
  
  # SiO3
  
  glm_mean_sio3 <- glm(response ~  SiO3ug, data = som_glm)
  ms_sum <- summary(glm_mean_sio3)
  model_AIC[5,2] <- ms_sum$aic
  model_AIC[5,1] <- "SiO4"
  model_AIC[5,3] <- cor(som_glm$SiO3ug, som_glm$response, use = "pairwise.complete.obs")
  
  # NC Depth
  
  glm_mean_nc <- glm(response ~  NCDepth, data = som_glm)
  ms_sum <- summary(glm_mean_nc)
  model_AIC[6,2] <- ms_sum$aic
  model_AIC[6,1] <- "NCD"
  model_AIC[6,3] <- cor(som_glm$NCDepth, som_glm$response, use = "pairwise.complete.obs")
  
  # chl
  
  glm_mean_mld <- glm(response ~  ChlorA, data = som_glm)
  ms_sum <- summary(glm_mean_mld)
  model_AIC[7,2] <- ms_sum$aic
  model_AIC[7,1] <- "Chl-a"
  model_AIC[7,3] <- cor(som_glm$ChlorA, som_glm$response, use = "pairwise.complete.obs")
  
  # o2
  
  glm_mean_mld <- glm(response ~  O2ml_L, data = som_glm)
  ms_sum <- summary(glm_mean_mld)
  model_AIC[8,2] <- ms_sum$aic
  model_AIC[8,1] <- "O2"
  model_AIC[8,3] <- cor(som_glm$O2ml_L, som_glm$response, use = "pairwise.complete.obs")
  
  
  return(model_AIC)
  
}


load("output/total_full_data.Rdata")

full_dat$Year <- as.numeric(substr(full_dat$Cruise,1,4))

full_dat <- full_dat %>%
  filter(Depthm < 15)
  
  partial_dat <- full_dat
  
  out <- aic_table_func_diveristy_sign(som_maps = partial_dat, col_n = 63)
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  colnames(out) <- c("Variables","AIC", "Slope")
  out[,2] <- as.numeric(out[,2])
  out[,2] <- round(out[,2], digits = 2)
  out[,3] <- as.numeric(out[,3])
  out[,3] <- round(out[,3], digits = 3)
  
  set.seed(17)
  
  rf_dat <- partial_dat[,c(63,34,35,40,38,39,49,42,37)]
  rf_dat <- rf_dat[complete.cases(rf_dat),]
  
  rf_out <- randomForest(richness ~., data = rf_dat, ntree=100, mtry=2, importance=TRUE)
  
  out$RandomForest <- rf_out$importance[,1]
  out$Slope <- NULL


AIC_full <- out

AIC_scaled <- AIC_full


  zero_one_scale <- 1-(AIC_scaled[,2]-min(AIC_scaled[,2], na.rm = TRUE))/
    abs(min(AIC_scaled[,2], na.rm = TRUE) - max(AIC_scaled[,2], na.rm = TRUE))
  AIC_scaled[,2] <- 50^zero_one_scale
  
  zero_one_scale <- 1-(max(AIC_scaled[,3], na.rm = TRUE)- AIC_scaled[,3])/
    abs(max(AIC_scaled[,3], na.rm = TRUE) - min(AIC_scaled[,3], na.rm = TRUE))
  AIC_scaled[,3] <- 50^zero_one_scale

plot_df <- melt(AIC_scaled[,1:3])

colnames(plot_df) <- c("Variables", "Metric", "Value")

plot_df$Metric <- as.character(plot_df$Metric)
plot_df$Metric[which(plot_df$Metric == "RandomForest")] <- "Random Forest\nMean Decrease Accuracy"

plot_df$Metric <- as.factor(plot_df$Metric)

plot_df$Metric <- factor(plot_df$Metric, levels = c("Random Forest\nMean Decrease Accuracy", "AIC"))

plot_df$Variables <- as.factor(plot_df$Variables)

plot_df$Variables <- factor(plot_df$Variables, levels = c("Temp", "Salinity",
                                                          "NO3", "PO4", 
                                                          "SiO4","NCD",
                                                          "Chl-a", "O2"))


outplot <- ggplot(data = plot_df, aes(x = Variables, y = Metric, size = Value)) + 
  geom_point(color = "black", alpha = 0.6, shape = 21, fill = "red") +
  labs(fill = "Correlation") + ylab("Variable") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.position = "right",
        panel.grid.major.x = element_line(color = "grey", linetype = 2),
        plot.title = element_text(hjust = 0, size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  scale_size_continuous(range = c(1,18), guide = FALSE) + xlab("") + ylab("") +
  ggtitle("Variable importance for total alpha diversity")

pdf(file = "figures/grant_var_import_glm_rf.pdf", width = 8, height = 3)
print(outplot)
dev.off()




aic_table_func_diveristy_sign <- function(som_maps = partial_dat, col_n = 63){  
  
  som_glm <- som_maps
  
  colnames(som_glm)[col_n] <- "response"
  
  som_glm$response <- as.numeric(som_glm$response)
  
  model_AIC <- matrix(NA,8,3)
  
  # temperature
  
  glm_mean_temp <- glm(response ~  T_degC, data = som_glm)
  mt_sum <- summary(glm_mean_temp)
  model_AIC[1,2] <- mt_sum$aic
  model_AIC[1,1] <- "Temp"
  model_AIC[1,3] <- cor(som_glm$T_degC, som_glm$response, use = "pairwise.complete.obs")
  
  # sea surface temperature
  
  # glm_mean_sst <- glm(response ~  sst_mean, data = som_glm)
  # mt_sum <- summary(glm_mean_sst)
  # model_AIC[2,2] <- mt_sum$aic
  # model_AIC[2,1] <- "Mean SST"
  # model_AIC[2,3] <- cor(som_glm$sst_mean, som_glm$response, use = "pairwise.complete.obs")
  
  # salinity
  
  glm_mean_sal <- glm(response ~  Salnty, data = som_glm)
  mt_sum <- summary(glm_mean_sal)
  model_AIC[2,2] <- mt_sum$aic
  model_AIC[2,1] <- "Salinity"
  model_AIC[2,3] <- cor(som_glm$Salnty, som_glm$response, use = "pairwise.complete.obs")
  
  
  # NO3
  
  glm_mean_no3 <- glm(response ~  NO3ug, data = som_glm)
  mn_sum <- summary(glm_mean_no3)
  model_AIC[3,2] <- mn_sum$aic
  model_AIC[3,1] <- "NO3"
  model_AIC[3,3] <- cor(som_glm$NO3ug, som_glm$response, use = "pairwise.complete.obs")
  
  # PO4
  
  glm_mean_po4 <- glm(response ~  PO4ug, data = som_glm)
  mp_sum <- summary(glm_mean_po4)
  model_AIC[4,2] <- mp_sum$aic
  model_AIC[4,1] <- "PO4"
  model_AIC[4,3] <- cor(som_glm$PO4ug, som_glm$response, use = "pairwise.complete.obs")
  
  # SiO3
  
  glm_mean_sio3 <- glm(response ~  SiO3ug, data = som_glm)
  ms_sum <- summary(glm_mean_sio3)
  model_AIC[5,2] <- ms_sum$aic
  model_AIC[5,1] <- "SiO4"
  model_AIC[5,3] <- cor(som_glm$SiO3ug, som_glm$response, use = "pairwise.complete.obs")
  
  # NC Depth
  
  glm_mean_nc <- glm(response ~  NCDepth, data = som_glm)
  ms_sum <- summary(glm_mean_nc)
  model_AIC[6,2] <- ms_sum$aic
  model_AIC[6,1] <- "NCD"
  model_AIC[6,3] <- cor(som_glm$NCDepth, som_glm$response, use = "pairwise.complete.obs")
  
  # chl
  
  glm_mean_mld <- glm(response ~  ChlorA, data = som_glm)
  ms_sum <- summary(glm_mean_mld)
  model_AIC[7,2] <- ms_sum$aic
  model_AIC[7,1] <- "Chl-a"
  model_AIC[7,3] <- cor(som_glm$ChlorA, som_glm$response, use = "pairwise.complete.obs")
  
  # o2
  
  glm_mean_mld <- glm(response ~  O2ml_L, data = som_glm)
  ms_sum <- summary(glm_mean_mld)
  model_AIC[8,2] <- ms_sum$aic
  model_AIC[8,1] <- "O2"
  model_AIC[8,3] <- cor(som_glm$O2ml_L, som_glm$response, use = "pairwise.complete.obs")
  
  
  return(model_AIC)
  
}



