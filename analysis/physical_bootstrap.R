library(tidyverse)

phys_df <- read.csv("data/NCOG_sample_log_DNA_meta_2014-2019.csv", header = TRUE)
surf <- TRUE
n_samps <- 500
set.seed(345)

if(surf == TRUE){
  phys_df <- phys_df %>%
    filter(Depthm < 15, as.numeric(substr(Sta_ID,1,3)) > 75)
}else{
  phys_df <- phys_df %>%
    filter(Depthm >= 15, as.numeric(substr(Sta_ID,1,3)) > 75)
}

phys_df <- phys_df[,c(3,34:35)]
phys_df <- phys_df[complete.cases(phys_df),]

sample_mat <- matrix(replicate(n = n_samps, sample(1:24, replace = FALSE)), nrow = 24, ncol = n_samps, byrow = FALSE)
mean_temp_diff_mat <- matrix(NA, 24, n_samps)
mean_sal_diff_mat <- matrix(NA, 24, n_samps)

cruises <- unique(phys_df$Cruise)

for (i in 1:ncol(sample_mat)) {
  print(i)
  for (j in 1:nrow(sample_mat)) {
 
    df_1 <- phys_df %>%
      filter(Cruise %in% cruises[sample_mat[1,i]])
    
    df_j <- phys_df %>%
      filter(Cruise %in% cruises[sample_mat[1:j,i]])
    
    t_similar_mat <- matrix(NA, nrow(df_1), nrow(df_j))
    s_similar_mat <- matrix(NA, nrow(df_1), nrow(df_j))
    
    for (k in 1:nrow(df_1)) {
      for (l in 1:nrow(df_j)) {
        t_similar_mat[k,l] <- abs(df_1$T_degC[k] - df_j$T_degC[l])
        s_similar_mat[k,l] <- abs(df_1$Salnty[k] - df_j$Salnty[l])
      }
    }
    
    mean_temp_diff_mat[j,i] <- mean(t_similar_mat[upper.tri(t_similar_mat, diag = FALSE)], na.rm = TRUE)
    mean_sal_diff_mat[j,i] <- mean(s_similar_mat[upper.tri(s_similar_mat, diag = FALSE)], na.rm = TRUE)
    }
       
  }


mean_temp_diff_mat <- as.data.frame(mean_temp_diff_mat)
mean_temp_diff_mat$cruise <- 1:nrow(mean_temp_diff_mat)

temp_df <- mean_temp_diff_mat %>%
  pivot_longer(-cruise, names_to = "sample", values_to = "Difference")

temp_df$Variable <- "Temperature"

mean_sal_diff_mat <- as.data.frame(mean_sal_diff_mat)
mean_sal_diff_mat$cruise <- 1:nrow(mean_sal_diff_mat)

sal_df <- mean_sal_diff_mat %>%
  pivot_longer(-cruise, names_to = "sample", values_to = "Difference")

sal_df$Variable <- "Salinity"

total_df <- bind_rows(temp_df, sal_df)

temp <- ggplot(temp_df, aes(x = cruise, y = Difference)) + 
  geom_point(alpha = 0.5, shape = 19, color = "black") +
  geom_smooth(method = "gam", level = 0.95, color = "red", fill = "red") + 
  xlab("# of Cruises Sampled") + ylab("Mean Temperature Difference (°C)") + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

sal <- ggplot(sal_df, aes(x = cruise, y = Difference)) + 
  geom_point(alpha = 0.5, shape = 19) +
  geom_smooth(method = "gam", level = 0.95, color = "blue", fill = "blue") + 
  xlab("# of Cruises Sampled") + ylab("Mean Salinity Difference") + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none")

fig_out <- temp + sal

pdf(file = "figures/temporal_bootstrap.pdf", width = 12, height = 5)
print(fig_out)
dev.off()
