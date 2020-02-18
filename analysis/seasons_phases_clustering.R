library(tidyverse)
library(vegan)
library(SOMbrero)
library(MASS)
library(mclust)
library(sBIC)
library(cowplot)
library(cluster)

physical_dat <- read.csv(file = "data/NCOG_sample_metadata.csv", stringsAsFactors = FALSE)

full_dat <- physical_dat

physical_dat <- physical_dat[,c(3,33,34,37:39)]

winters <- unique(physical_dat$Cruise)[seq(1,20,4)]
springs <- unique(physical_dat$Cruise)[seq(2,20,4)]
summers <- unique(physical_dat$Cruise)[seq(3,20,4)]
falls <- unique(physical_dat$Cruise)[seq(4,20,4)]

early <- unique(physical_dat$Cruise)[1:12]
late <- unique(physical_dat$Cruise)[13:20]

winter_dat <- physical_dat %>%
  subset(Cruise %in% winters)

spring_dat <- physical_dat %>%
  subset(Cruise %in% springs)

summer_dat <- physical_dat %>%
  subset(Cruise %in% summers)

fall_dat <- physical_dat %>%
  subset(Cruise %in% falls)

early_dat <- physical_dat %>%
  subset(Cruise %in% early)

late_dat <- physical_dat %>%
  subset(Cruise %in% late)

winter_dat$Cruise <- NULL
spring_dat$Cruise <- NULL
summer_dat$Cruise <- NULL
fall_dat$Cruise <- NULL
early_dat$Cruise <- NULL
late_dat$Cruise <- NULL
physical_dat$Cruise <- NULL

winter_dat <- winter_dat[complete.cases(winter_dat),]
spring_dat <- spring_dat[complete.cases(spring_dat),]
summer_dat <- summer_dat[complete.cases(summer_dat),]
fall_dat <- fall_dat[complete.cases(fall_dat),]
early_dat <- early_dat[complete.cases(early_dat),]
late_dat <- late_dat[complete.cases(late_dat),]
physical_dat <- physical_dat[complete.cases(physical_dat),]


win.som <- trainSOM(x.data = winter_dat, dimension = c(5, 5), nb.save = 10, maxit = 2000, 
                    scaling = "none")

spr.som <- trainSOM(x.data = spring_dat, dimension = c(5, 5), nb.save = 10, maxit = 2000, 
                    scaling = "none")

sum.som <- trainSOM(x.data = summer_dat, dimension = c(5, 5), nb.save = 10, maxit = 2000, 
                    scaling = "none")

fal.som <- trainSOM(x.data = fall_dat, dimension = c(5, 5), nb.save = 10, maxit = 2000, 
                    scaling = "none")

ear.som <- trainSOM(x.data = early_dat, dimension = c(5, 5), nb.save = 10, maxit = 2000, 
                    scaling = "none")

lat.som <- trainSOM(x.data = late_dat, dimension = c(5, 5), nb.save = 10, maxit = 2000, 
                    scaling = "none")

phy.som <- trainSOM(x.data = physical_dat, dimension = c(5, 5), nb.save = 10, maxit = 2000, 
                    scaling = "none")


## Silohuette
silhouette_score <- function(k, df = win.som$prototypes){
  km <- kmeans(df, centers = k, nstart=25)
  ss <- silhouette(km$cluster, dist(df))
  mean(ss[, 3])
}
k <- 2:10

avg_sil <- sapply(k, silhouette_score, df = win.som$prototypes)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Winter SOM Cluster")

avg_sil <- sapply(k, silhouette_score, df = winter_dat)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Winter K-Means Cluster")

avg_sil <- sapply(k, silhouette_score, df = spr.som$prototypes)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Spring SOM Cluster")

avg_sil <- sapply(k, silhouette_score, df = spring_dat)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Spring K-Means Cluster")

avg_sil <- sapply(k, silhouette_score, df = sum.som$prototypes)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Summer SOM Cluster")

avg_sil <- sapply(k, silhouette_score, df = summer_dat)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Summer K-Means Cluster")

avg_sil <- sapply(k, silhouette_score, df = fal.som$prototypes)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Fall SOM Cluster")

avg_sil <- sapply(k, silhouette_score, df = fall_dat)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Fall K-Means Cluster")

avg_sil <- sapply(k, silhouette_score, df = ear.som$prototypes)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Early SOM Cluster")

avg_sil <- sapply(k, silhouette_score, df = early_dat)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Early K-Means Cluster")

avg_sil <- sapply(k, silhouette_score, df = lat.som$prototypes)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Late SOM Cluster")

avg_sil <- sapply(k, silhouette_score, df = late_dat)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Late K-Means Cluster")

par(mfrow = c(1,2))

avg_sil <- sapply(k, silhouette_score, df = phy.som$prototypes)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Physical SOM Cluster")

avg_sil <- sapply(k, silhouette_score, df = physical_dat)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Physical K-Means Cluster")

load("output/cyano_16s_som.Rdata")
load("data/16s_cyanos.Rdata")

avg_sil <- sapply(k, silhouette_score, df = eco.som$prototypes)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Cyano SOM Cluster")

avg_sil <- sapply(k, silhouette_score, df = scaled_inputs)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Cyano K-Means Cluster")

load("output/bacteria_m_euks_16s_som.Rdata")
load("data/16s_bacteria_m_euks.Rdata")

avg_sil <- sapply(k, silhouette_score, df = eco.som$prototypes)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Heterotrophic Bacteria/Archaea SOM Cluster")

avg_sil <- sapply(k, silhouette_score, df = scaled_inputs)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Heterotrophic Bacteria/Archaea K-Means Cluster")

load("output/euks_auto_18sv9_som.Rdata")
load("data/18s_autotrophic_euks.Rdata")

avg_sil <- sapply(k, silhouette_score, df = eco.som$prototypes)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Eukaryotic Phytoplankton SOM Cluster")

avg_sil <- sapply(k, silhouette_score, df = scaled_inputs)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Eukaryotic Phytoplankton K-Means Cluster")

load("output/euks_hetero_18sv9_som.Rdata")
load("data/18s_heterotrophic_euks.Rdata")

avg_sil <- sapply(k, silhouette_score, df = eco.som$prototypes)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Heterotrophic Eukaryotes SOM Cluster")

avg_sil <- sapply(k, silhouette_score, df = scaled_inputs)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE, main = "Heterotrophic Eukaryotes K-Means Cluster")

# winter
win.bic <- mclustBIC(winter_dat, G = 1:8)
plot(win.bic)

win.som.bic <- mclustBIC(win.som$prototypes, G = 1:8)
plot(win.som.bic)

# spring
spr.bic <- mclustBIC(spring_dat, G = 1:8)
plot(spr.bic)

spr.som.bic <- mclustBIC(spr.som$prototypes, G = 1:8)
plot(spr.som.bic)

# summer
sum.bic <- mclustBIC(summer_dat, G = 1:8)
plot(sum.bic)

sum.som.bic <- mclustBIC(sum.som$prototypes, G = 1:8)
plot(sum.som.bic)

# fall
fal.bic <- mclustBIC(fall_dat, G = 1:8)
plot(fal.bic)

fal.som.bic <- mclustBIC(fal.som$prototypes, G = 1:8)
plot(fal.som.bic)

# early
ear.bic <- mclustBIC(early_dat, G = 1:8)
plot(ear.bic)

ear.som.bic <- mclustBIC(ear.som$prototypes, G = 1:8)
plot(ear.som.bic)

# late
lat.bic <- mclustBIC(late_dat, G = 1:8)
plot(lat.bic)

lat.som.bic <- mclustBIC(lat.som$prototypes, G = 1:8)
plot(lat.som.bic)

# all
phy.bic <- mclustBIC(physical_dat, G = 1:8)
plot(phy.bic)

phy.som.bic <- mclustBIC(phy.som$prototypes, G = 1:8)
plot(phy.som.bic)

# ssi 

#winter
som_kmeans_cascade <- cascadeKM(win.som$prototypes, inf.gr = 2,
                                sup.gr = 15, 
                                iter = 10000,
                                criterion = "calinski")

summary(som_kmeans_cascade)
plot(som_kmeans_cascade, sortg = TRUE)

#spring
som_kmeans_cascade <- cascadeKM(spr.som$prototypes, inf.gr = 2,
                                sup.gr = 15, 
                                iter = 10000,
                                criterion = "calinski")

summary(som_kmeans_cascade)
plot(som_kmeans_cascade, sortg = TRUE)

#summer
som_kmeans_cascade <- cascadeKM(sum.som$prototypes, inf.gr = 2,
                                sup.gr = 15, 
                                iter = 10000,
                                criterion = "calinski")

summary(som_kmeans_cascade)
plot(som_kmeans_cascade, sortg = TRUE)

#fall
som_kmeans_cascade <- cascadeKM(fal.som$prototypes, inf.gr = 2,
                                sup.gr = 15, 
                                iter = 10000,
                                criterion = "calinski")

summary(som_kmeans_cascade)
plot(som_kmeans_cascade, sortg = TRUE)

#early
som_kmeans_cascade <- cascadeKM(ear.som$prototypes, inf.gr = 2,
                                sup.gr = 15, 
                                iter = 10000,
                                criterion = "calinski")

summary(som_kmeans_cascade)
plot(som_kmeans_cascade, sortg = TRUE)

#late
som_kmeans_cascade <- cascadeKM(lat.som$prototypes, inf.gr = 2,
                                sup.gr = 15, 
                                iter = 10000,
                                criterion = "calinski")

summary(som_kmeans_cascade)
plot(som_kmeans_cascade, sortg = TRUE)





# TS Figures

par(mfrow = c(3,2))

plot(early_dat$T_degC, early_dat$Salnty, xlim = c(10,23), ylim = c(32.6,34), col = "red", pch = 20,
     xlab = "Temperature", ylab = "Salinity", main = "Early")

plot(late_dat$T_degC, late_dat$Salnty, xlim = c(10,23), ylim = c(32.6,34), col = "blue", pch = 20,
     xlab = "Temperature", ylab = "Salinity", main = "Late")

plot(winter_dat$T_degC, winter_dat$Salnty, xlim = c(10,23), ylim = c(32.6,34), col = "blue", pch = 20,
     xlab = "Temperature", ylab = "Salinity", main = "Winter")

plot(spring_dat$T_degC, spring_dat$Salnty, xlim = c(10,23), ylim = c(32.6,34), col = "green", pch = 20,
     xlab = "Temperature", ylab = "Salinity", main = "Spring")

plot(summer_dat$T_degC, summer_dat$Salnty, xlim = c(10,23), ylim = c(32.6,34), col = "orange", pch = 20,
     xlab = "Temperature", ylab = "Salinity", main = "Summer")

plot(fall_dat$T_degC, fall_dat$Salnty, xlim = c(10,23), ylim = c(32.6,34), col = "red", pch = 20,
     xlab = "Temperature", ylab = "Salinity", main = "Fall")

# dissecting phys

full_dat$Season <- full_dat$Cruise
full_dat$Phase <- full_dat$Cruise

full_dat$Season[which(!is.na(match(full_dat$Season, winters)))] <- "Winter"
full_dat$Season[which(!is.na(match(full_dat$Season, springs)))] <- "Spring"
full_dat$Season[which(!is.na(match(full_dat$Season, summers)))] <- "Summer"
full_dat$Season[which(!is.na(match(full_dat$Season, falls)))] <- "Fall"

full_dat$Phase[which(!is.na(match(full_dat$Phase, early)))] <- "Early"
full_dat$Phase[which(!is.na(match(full_dat$Phase, late)))] <- "Late"

# TS diagrams
wint <- ggplot(full_dat %>% filter(Season == "Winter"),
       aes(x = T_degC, color = Phase, y = Salnty, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23) + ylim(32.6, 34) + ggtitle("Winter") + stat_ellipse()

sprin <- ggplot(full_dat %>% filter(Season == "Spring"),
       aes(x = T_degC, color = Phase, y = Salnty, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23) + ylim(32.6, 34) + ggtitle("Spring") + stat_ellipse()

summ <- ggplot(full_dat %>% filter(Season == "Summer"),
                aes(x = T_degC, color = Phase, y = Salnty, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23) + ylim(32.6, 34) + ggtitle("Summer") + stat_ellipse()

fall <- ggplot(full_dat %>% filter(Season == "Fall"),
                aes(x = T_degC, color = Phase, y = Salnty, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23) + ylim(32.6, 34) + ggtitle("Fall") + stat_ellipse()


plot_grid(wint, sprin, summ, fall, nrow = 2, ncol = 2)

# T-Nutrients
wint <- ggplot(full_dat %>% filter(Season == "Winter"),
               aes(x = T_degC, color = Phase, y = NO3ug, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23)  + ggtitle("Winter") + stat_ellipse() + ylim(0, 10)

sprin <- ggplot(full_dat %>% filter(Season == "Spring"),
                aes(x = T_degC, color = Phase, y = NO3ug, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23) + ylim(0, 20) + ggtitle("Spring") + stat_ellipse()

summ <- ggplot(full_dat %>% filter(Season == "Summer"),
               aes(x = T_degC, color = Phase, y = NO3ug, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23) + ylim(0, 30) + ggtitle("Summer") + stat_ellipse()

fall <- ggplot(full_dat %>% filter(Season == "Fall"),
               aes(x = T_degC, color = Phase, y = NO3ug, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23) + ylim(0, 20) + ggtitle("Fall") + stat_ellipse()


plot_grid(wint, sprin, summ, fall, nrow = 2, ncol = 2)


wint <- ggplot(full_dat %>% filter(Season == "Winter"),
               aes(x = T_degC, color = Phase, y = PO4ug, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23)  + ggtitle("Winter") + stat_ellipse() + ylim(0, 2)

sprin <- ggplot(full_dat %>% filter(Season == "Spring"),
                aes(x = T_degC, color = Phase, y = PO4ug, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23) + ylim(0, 2) + ggtitle("Spring") + stat_ellipse()

summ <- ggplot(full_dat %>% filter(Season == "Summer"),
               aes(x = T_degC, color = Phase, y = PO4ug, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23) + ylim(0, 2) + ggtitle("Summer") + stat_ellipse()

fall <- ggplot(full_dat %>% filter(Season == "Fall"),
               aes(x = T_degC, color = Phase, y = PO4ug, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23) + ylim(0, 2) + ggtitle("Fall") + stat_ellipse()


plot_grid(wint, sprin, summ, fall, nrow = 2, ncol = 2)

wint <- ggplot(full_dat %>% filter(Season == "Winter"),
               aes(x = T_degC, color = Phase, y = SiO3ug, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23)  + ggtitle("Winter") + stat_ellipse() + ylim(0, 5)

sprin <- ggplot(full_dat %>% filter(Season == "Spring"),
                aes(x = T_degC, color = Phase, y = SiO3ug, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23) + ylim(0, 5) + ggtitle("Spring") + stat_ellipse()

summ <- ggplot(full_dat %>% filter(Season == "Summer"),
               aes(x = T_degC, color = Phase, y = SiO3ug, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23) + ylim(0, 5) + ggtitle("Summer") + stat_ellipse()

fall <- ggplot(full_dat %>% filter(Season == "Fall"),
               aes(x = T_degC, color = Phase, y = SiO3ug, shape = Phase)) +
  geom_point() + theme(panel.background = element_blank(),
                       panel.border = element_rect(fill = NA, color = "black"),
                       plot.title = element_text(hjust = 0.5)) +
  xlim(10,23) + ylim(0, 5) + ggtitle("Fall") + stat_ellipse()


plot_grid(wint, sprin, summ, fall, nrow = 2, ncol = 2)
