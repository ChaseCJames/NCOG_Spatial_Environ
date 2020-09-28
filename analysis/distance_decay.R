library(tidyverse)
library(spatialEco)
library(sp)
library(geosphere)
library(patchwork)

in_data = "output/total_dissimilar.Rdata"
phys_dat = "data/NCOG_sample_log_DNA_meta_2014-2019.csv"
surf <- TRUE

load(in_data)
phys_df <- read.csv(phys_dat, header = TRUE)

if(surf == TRUE){
  phys_df <- phys_df %>%
    filter(Depthm < 15, as.numeric(substr(Sta_ID,1,3)) > 75)
}else{
  phys_df <- phys_df %>%
    filter(Depthm >= 15, as.numeric(substr(Sta_ID,1,3)) > 75)
}

matches <- which(!is.na(match(rownames(dissimilar), paste0("X",phys_df$Sample.Name))))
dissimilar <- dissimilar[matches,matches]

phys_matches <- which(!is.na(match(paste0("X",phys_df$Sample.Name),rownames(dissimilar))))
phys_df <- phys_df[phys_matches,]

long_lat <-  phys_df[,c(29,28)]

phys_distance <- distm(long_lat, long_lat, fun = distGeo)
phys_distance <- phys_distance/1000

rownames(phys_distance) <- paste0("X",phys_df$Sample.Name)
colnames(phys_distance) <- paste0("X",phys_df$Sample.Name)

phys_df$Scale_Temp <- scale(phys_df$T_degC)
phys_df$Scale_Sal <- scale(phys_df$Salnty)
phys_df$Scale_NO3_NH3 <- scale((phys_df$NO3ug + phys_df$NH3ug))

env_vars <- phys_df[,58:60]

env_distance <- as.matrix(dist(env_vars))

rownames(env_distance) <- paste0("X",phys_df$Sample.Name)
colnames(env_distance) <- paste0("X",phys_df$Sample.Name)

diag(dissimilar) <- NA
diag(env_distance) <- NA
diag(phys_distance) <- NA

dissimilar <- as.data.frame(dissimilar)
phys_distance <- as.data.frame(phys_distance)
env_distance <- as.data.frame(env_distance)

dissimilar$sample2 <- rownames(dissimilar)
phys_distance$sample2 <- rownames(phys_distance)
env_distance$sample2 <- rownames(env_distance)

diss_long <- pivot_longer(dissimilar, -sample2, names_to = "sample1", values_to = "Dissimilarity")
phys_long <- pivot_longer(phys_distance, -sample2, names_to = "sample1", values_to = "Physical")
env_long <- pivot_longer(env_distance, -sample2, names_to = "sample1", values_to = "Environmental")

diss_long$comb <- paste0(diss_long$sample1, diss_long$sample2)
phys_long$comb <- paste0(phys_long$sample1, phys_long$sample2)
env_long$comb <- paste0(env_long$sample1, env_long$sample2)

diss_long <- diss_long[,3:4]
phys_long <- phys_long[,3:4]
env_long <- env_long[,3:4]

data_long <- full_join(diss_long, phys_long, by = "comb")
data_long <- full_join(data_long, env_long, by = "comb")

data_long <- data_long[complete.cases(data_long),]

sub_sample <- sample(1:nrow(data_long), size = 0.01*nrow(data_long))

data_subsample <- data_long[sub_sample,]

aplot <- ggplot(data_subsample, aes(x = Physical, y = Environmental, color = Dissimilarity)) + 
  geom_point() + labs(color = "Bray-Curtis\nDissimilarity") +
  ylab("Environmental Eucledian Distance") + xlab("Distance (km)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  scale_color_gradient(low = "red", high = "blue") +
  ggtitle("A. Physical vs. Environmental Distance")

bplot <- ggplot() +
  geom_point(data = data_subsample, aes(x = Physical, y = Dissimilarity)) +
  stat_smooth(data = data_long, aes(x = Physical, y = Dissimilarity), color = "red", fill = "red", method = "gam") +
  ylab("Bray-Curtis Dissimilarity") + xlab("Distance (km)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  ggtitle("B. Phyiscal Distance vs BC Dissimilarity")

cplot <- ggplot() +
  geom_point(data = data_subsample, aes(x = Environmental, y = Dissimilarity)) +
  stat_smooth(data = data_long, aes(x = Environmental, y = Dissimilarity),color = "red", fill = "red", method = "gam") +
  ylab("Bray-Curtis Dissimilarity") + xlab("Environmental Eucledian Distance") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  ggtitle("C. Environmental Distance vs BC Dissimilarity")

out_plot <- aplot / (bplot + cplot)

pdf(file = "figures/bc_similarity_v_distance.pdf", width = 10, height = 10)
out_plot
dev.off()

