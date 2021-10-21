n_boot <- 1000
prop <- 2/3

bc_mat <- matrix(NA,n_boot,8) %>% as.data.frame()
sha_mat <- matrix(NA,n_boot,8) %>% as.data.frame()
rich_mat <- matrix(NA,n_boot,8) %>% as.data.frame()


colnames(bc_mat) <- c("Environment", "Filter", "Environment + Filter", "Null", "Time", "Time + Filter", "Environment + Time", "All")
colnames(sha_mat) <- c("Environment", "Filter", "Environment + Filter", "Null", "Time", "Time + Filter", "Environment + Time", "All")
colnames(rich_mat) <-c("Environment", "Filter", "Environment + Filter", "Null", "Time", "Time + Filter", "Environment + Time", "All")


for (i in 1:n_boot) {
  print(i)
  
  samps <- sample(1:nrow(long_similarity), size = round(prop*nrow(long_similarity),digits = 0))
  
  sub <- long_similarity[samps,]
  
  # Bray Curtis
  
  long_bc <- sub %>% mutate(bc_anom = bc_dist - mean(bc_dist))
  long_bc <- long_bc %>% group_by(filter) %>% mutate(filt_bc_mean = mean(bc_dist)) %>% mutate(filt_bc_anom = bc_dist - filt_bc_mean)
  
  out_1 <- gam(bc_dist ~ e_dist, data = long_bc) 
  out_3 <- gam(bc_dist ~ e_dist + filter, data = long_bc) 
  out_4 <- gam(bc_dist ~ time_diff, data = long_bc) 
  out_5 <- gam(bc_dist ~ time_diff + filter, data = long_bc) 
  out_6 <- gam(bc_dist ~ e_dist + time_diff, data = long_bc) 
  out_7 <- gam(bc_dist ~ e_dist + time_diff + filter, data = long_bc) 
  
  res_bc <- out_1$residuals %>% as.data.frame()
  colnames(res_bc) <- "Environment"
  res_bc$Filter <- long_bc$filt_bc_anom
  res_bc$Both <- out_3$residuals
  res_bc$Null <- long_bc$bc_anom
  res_bc$Time <- out_4$residuals
  res_bc$Time_filt <- out_5$residuals
  res_bc$Time_env <- out_6$residuals
  res_bc$all <- out_7$residuals
  
  bc_mat[i,] <- apply(res_bc, 2, function(x) sqrt(mean((x)^2))) %>% as.numeric()
  
  
  # Shannon
  
  long_sha <- sub %>% mutate(sha_anom = shannon_diff - mean(shannon_diff))
  long_sha <- long_sha %>% group_by(filter) %>% mutate(filt_sha_mean = mean(shannon_diff)) %>% mutate(filt_sha_anom = shannon_diff - filt_sha_mean)
  
  out_1 <- gam(shannon_diff ~ e_dist, data = long_sha) 
  out_3 <- gam(shannon_diff ~ e_dist + filter, data = long_sha) 
  out_4 <- gam(shannon_diff ~ time_diff, data = long_sha) 
  out_5 <- gam(shannon_diff ~ time_diff + filter, data = long_sha) 
  out_6 <- gam(shannon_diff ~ e_dist + time_diff, data = long_sha) 
  out_7 <- gam(shannon_diff ~ e_dist + time_diff + filter, data = long_sha) 
 
  
  res_sha <- out_1$residuals %>% as.data.frame()
  colnames(res_sha) <- "Environment"
  res_sha$Filter <- long_sha$filt_sha_anom
  res_sha$Both <- out_3$residuals
  res_sha$Null <- long_sha$sha_anom
  res_sha$Time <- out_4$residuals
  res_sha$Time_filt <- out_5$residuals
  res_sha$Time_env <- out_6$residuals
  res_sha$all <- out_7$residuals

  
  sha_mat[i,] <- apply(res_sha, 2, function(x) sqrt(mean((x)^2))) %>% as.numeric()
  
  # Richness
  
  long_rich <- sub %>% mutate(rich_anom = rich_diff - mean(rich_diff))
  long_rich <- long_rich %>% group_by(filter) %>% mutate(filt_rich_mean = mean(rich_diff)) %>% mutate(filt_rich_anom = rich_diff - filt_rich_mean)
  
  out_1 <- gam(rich_diff ~ e_dist, data = long_rich) 
  out_3 <- gam(rich_diff ~ e_dist + filter, data = long_rich) 
  out_4 <- gam(rich_diff ~ time_diff, data = long_rich) 
  out_5 <- gam(rich_diff ~ time_diff + filter, data = long_rich) 
  out_6 <- gam(rich_diff ~ e_dist + time_diff, data = long_rich) 
  out_7 <- gam(rich_diff ~ e_dist + time_diff +filter, data = long_rich) 

  
  res_rich <- out_1$residuals %>% as.data.frame()
  colnames(res_rich) <- "Environment"
  res_rich$Filter <- long_rich$filt_rich_anom
  res_rich$Both <- out_3$residuals
  res_rich$Null <- long_rich$rich_anom
  res_rich$Time <- out_4$residuals
  res_rich$Time_filt <- out_5$residuals
  res_rich$Time_env <- out_6$residuals
  res_rich$all <- out_7$residuals
  
  rich_mat[i,] <- apply(res_rich, 2, function(x) sqrt(mean((x)^2))) %>% as.numeric()
}


bc_rmse <- bc_mat %>% pivot_longer(everything(), names_to = "Model", values_to = "RMSE")


bc_rmse$Model <- factor(bc_rmse$Model,
                         levels = c("Null", "Filter", "Environment", "Time",
                                    "Environment + Filter", "Time + Filter",
                                    "Environment + Time", "All"))


sha_rmse <- sha_mat %>% pivot_longer(everything(), names_to = "Model", values_to = "RMSE") 


sha_rmse$Model <- factor(sha_rmse$Model,
                         levels = c("Null", "Filter", "Environment", "Time",
                                    "Environment + Filter", "Time + Filter",
                                    "Environment + Time", "All"))


rich_rmse <- rich_mat %>% pivot_longer(everything(), names_to = "Model", values_to = "RMSE") 


rich_rmse$Model <- factor(rich_rmse$Model,
                          levels = c("Null", "Filter", "Environment", "Time",
                                     "Environment + Filter", "Time + Filter",
                                     "Environment + Time", "All"))




bc <- ggplot(bc_rmse, aes(x = Model, y = RMSE)) +
  geom_violin(draw_quantiles = 0.5, fill = "#7aa457") + 
  ggtitle("Bray-Curtis Similarity") +
  labs(x = "Model", y = "RMSE") +
  theme(panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_blank())

shan <- ggplot(sha_rmse, aes(x = Model, y = RMSE)) +
  geom_violin(draw_quantiles = 0.5, fill = "#a46cb7") + 
  ggtitle("Shannon Index") +
  labs(x = "Model", y = "RMSE") +
  theme(panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_blank())

rich <- ggplot(rich_rmse, aes(x = Model, y = RMSE)) +
  geom_violin(draw_quantiles = 0.5, fill = "#cb6a49") + 
  ggtitle("Richness") +
  labs(x = "Model", y = "RMSE") +
  theme(panel.border = element_rect(fill = NA, color = "black"),
        panel.background = element_blank())


out <- bc / shan / rich

agg_png("figures/filter_issues.png", width = 10, height = 8, units = "in", res = 400)
plot(out)
dev.off()



