library(SOMbrero)
library(vegan)
library(tidyverse)
library(rworldmap)
library(ggplot2)
library(ggmap)
library(scales)
library(patchwork)
library(ggmap)
library(geosphere)
library(vegan)
library(iNEXT)
library(nlstools)



som_endemics <- function(phys_dat = "data/NCOG_sample_log_DNA_meta_2014-2019.csv",
                         k.val = 2, surf = TRUE,
                         in_data = "data/16s_cyanos.Rdata", title = "Cyanobacteria",
                         out_file = "output/cyano_som_endemics.Rdata"){

in_dat <- read.csv(phys_dat)

#just south stations
in_dat <- in_dat %>% 
  filter(as.numeric(substr(Sta_ID,1,3)) > 75)

if(surf == TRUE){in_dat <- in_dat %>% filter(Depthm <= 15)}
if(surf != TRUE){in_dat <- in_dat %>% filter(Depthm > 15)}

#all phys
# in_dat <- in_dat[,c(1,34:35,37:40,42,49)]
in_dat <- in_dat[,c(1,34:35)]


in_dat <- in_dat[complete.cases(in_dat),]

scaled_dat <- in_dat[,-1]

set.seed(23)
phys.som <- trainSOM(x.data = scaled_dat, dimension = c(6, 6), nb.save = 10, maxit = 2000)

# summary(eco.som)

phys.clust <- superClass(phys.som, k = k.val)

# plot(phys.clust)

clusters <- phys.clust$cluster

ids <- phys.clust$som$clustering

som_ids <- clusters[ids]

in_dat$som_id <- som_ids

# nmds_plot <- metaMDS(scaled_dat, k = 2)
# 
# points_df <- as.data.frame(nmds_plot$points)
# points_df$cluster <- som_ids
# 
# ggplot(data = points_df, aes(x = MDS1, y = MDS2, color = as.factor(cluster))) +
#   geom_point()

load(in_data)

splits <- strsplit(rownames(asv_table), "_")

stations <- vector()
depths <- vector()

for (i in 1:length(splits)) {
  stations[i] <- paste0(splits[[i]][2]," ",splits[[i]][3])
  depths[i] <-splits[[i]][4]
}

## surface only? 
if(surf == TRUE){
  asv_table <- asv_table[which(depths <= 15),]
  scaled_inputs <- scaled_inputs[which(depths <= 15),]
  stations <- stations[which(depths <= 15)]
}else{
  asv_table <- asv_table[which(depths > 15),]
  scaled_inputs <- scaled_inputs[which(depths > 15),]
  stations <- stations[which(depths > 15)]
}

## just south stations

asv_table <- asv_table[which(as.numeric(substr(stations,1,3)) > 75),]
scaled_inputs <- scaled_inputs[which(as.numeric(substr(stations,1,3)) > 75),]
stations <- stations[which(as.numeric(substr(stations,1,3)) > 75)]

asv_rare <- asv_table
asv_rare$year_month <- paste0(substr(rownames(asv_table),6,7), "-",substr(rownames(asv_table),2,5))
asv_rare$som_id <- in_dat$som_id[match(substr(rownames(asv_table),2,22),
                                       substr(in_dat$Sample.Name,1,21))]

cruises <- unique(asv_rare$year_month)

asv_rare <- asv_rare[which(!is.na(asv_rare$som_id)),]

# endemic 

end_mat <- as.data.frame(matrix(NA,nrow = length(cruises), ncol = (ncol(asv_rare)-2)))
colnames(end_mat) <- colnames(asv_rare)[1:(ncol(asv_rare)-2)]
end_mat$Date <- as.Date(paste0("01-",cruises), format = "%d-%m-%Y")

end_station <- as.data.frame(matrix(0, nrow = (ncol(asv_rare)-2), ncol = 3))
colnames(end_station) <- c("Station", "End_Duration", "Reads")
end_station$ASV <- colnames(asv_rare)[1:(ncol(asv_rare)-2)]

for (i in 1:(ncol(asv_rare)-2)) {
  
  if(((i/250)%%1==0) == TRUE){print(paste0(i,"/",(ncol(asv_rare)-2)))}
  
  end_vect <- vector()
  stat_list <- list()
  
  for (j in 1:length(cruises)) {
    
    match_cruise <- which(!is.na(match(asv_rare$year_month, cruises[1:j])))
    
    subset <- asv_rare[match_cruise,c(i,(ncol(asv_rare)-1):ncol(asv_rare))]
    
    by_stat <- subset %>%
      group_by(som_id) %>%
      summarise_at(colnames(subset)[1], sum, na.rm = TRUE)
    
    colnames(by_stat)[2] <- "ASV"
    
    ids <- length(which(by_stat$ASV != 0))
    if(ids == 0){stat_list[[j]] <- NULL}else{stat_list[[j]] <- by_stat$som_id[which(by_stat$ASV != 0)]}
    
    
    if(ids == 1){
      if(length(which(unique(stat_list) != "NULL")) == 1){
        end_station[i,1] <- by_stat$som_id[which(by_stat$ASV != 0)]
      }
      end_station[i,3] <- end_station[i,3] + by_stat$ASV[which(by_stat$ASV != 0)]
      
    }
    
    end_mat[j,i] <- length(which(unique(stat_list) != "NULL"))
    end_vect[j] <- length(which(unique(stat_list) != "NULL"))
    
  }
  
  ends <-  which(end_vect == 1)
  if(length(ends > 0)){end_station[i,2] <- length(ends)}
  
}

end_station$Station[end_station$Station == 0] <- NA

endemic_spp <- which(!is.na(end_station$Station))
y_m <- (ncol(asv_rare)-1)

asv_ts <- asv_rare[,c(endemic_spp,y_m)]

asv_ts <- asv_ts %>%
  group_by(year_month) %>%
  summarise_all(sum)

asv_ts <- as.data.frame(asv_ts)

end_station <- end_station[complete.cases(end_station),]

decay_mat <- as.data.frame(matrix(NA,
             nrow = (max(end_station$End_Duration, na.rm = TRUE)-1), ncol = 3))
colnames(decay_mat) <- c("Cruise", "Proportion", "Count")
decay_mat$Cruise <- 1:(max(end_station$End_Duration, na.rm = TRUE)-1)

for (i in 1:(max(end_station$End_Duration, na.rm = TRUE)-1)) {
    decay_mat$Proportion[i] <- length(which(end_station$End_Duration >= i))/nrow(end_station)
    decay_mat$Count[i] <- length(which(end_station$End_Duration >= i))
}

fit <- nls(Count ~ yf + (y0 - yf) * exp(-alpha * Cruise),
           start = list(y0 = max(decay_mat$Count), yf = 5, alpha = 0.1),
           data = decay_mat, lower = c(-Inf,0,-Inf),
           algorithm = "port")

con_inter <- confint2(fit)

output <- summary(fit)
output$parameters

y0 <- output$parameters[1,1]
y0_low <- con_inter[1,1]
y0_high <- con_inter[1,2]

yf <- output$parameters[2,1]
yf_low <- con_inter[2,1]
yf_high <- con_inter[2,2]

alpha <- output$parameters[3,1]
alpha_low <- con_inter[3,1]
alpha_high <- con_inter[3,2]

timevalues <- seq(1,(max(end_station$End_Duration, na.rm = TRUE)-1), 0.1)

pred_mat <- as.data.frame(matrix(NA,length(timevalues),4))
pred_mat$V1 <- timevalues
pred_mat$V2 <- yf + (y0 - yf)*exp(-alpha*pred_mat$V1)
pred_mat$V3 <- yf_low + (y0_low - yf_low)*exp(-alpha_high*pred_mat$V1)
pred_mat$V4 <- yf_high + (y0_high - yf_high)*exp(-alpha_low*pred_mat$V1)


if(yf == 0){
  out_plot <- ggplot() +
    geom_ribbon(data = pred_mat, aes(x = V1, ymin = V3, ymax = V4),
                alpha = 0.3, color = "red", fill = "red") +
    geom_line(data = pred_mat, aes(x = V1, y = V2),
              color = "red", size  = 1) +
    geom_point(data = decay_mat, aes(x = Cruise, y = Count), size = 2) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.margin = ) +
    labs(x = "# of Cruises as Endemic (t)", y = "# of Endemics") +
    annotate("text", x=17, y=0.8*max(decay_mat$Count),
             label = paste("Proportion(t) == ",
                           round((y0-yf),2),
                           "*e^{",round(alpha,2),"*t}"),
             parse =TRUE) +
    ggtitle(paste0(title, " Endemic Decay")) +
    coord_cartesian(xlim = c(1,23), ylim = c(0,max(pred_mat)), expand = FALSE)
  
  
}else{
  out_plot <- ggplot() +
    geom_ribbon(data = pred_mat, aes(x = V1, ymin = V3, ymax = V4),
                alpha = 0.3, color = "red", fill = "red") +
    geom_line(data = pred_mat, aes(x = V1, y = V2),
              color = "red", size  = 1) +
    geom_point(data = decay_mat, aes(x = Cruise, y = Count), size = 2) + 
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          plot.margin = ) +
    labs(x = "# of Cruises as Endemic (t)", y = "# of Endemics") +
    annotate("text", x=15, y=0.8*max(decay_mat$Count),
             label = paste("Proportion(t) == ",
                           round((y0-yf),2),
                           "*e^{",round(alpha,2),"*t} + ",round(yf,2)),
             parse =TRUE) +
    ggtitle(paste0(title, " Endemic Decay")) +
    coord_cartesian(xlim = c(1,23), ylim = c(0,max(pred_mat)), expand = FALSE)
}

print(out_plot)

save(out_plot, asv_rare, decay_mat, end_mat, end_station, asv_ts, in_dat, file = out_file)

}



som_endemics(phys_dat = "data/NCOG_sample_log_DNA_meta_2014-2019.csv",
             k.val = 2, surf = TRUE,
             in_data = "data/16s_cyanos.Rdata", title = "Cyanobacteria",
             out_file = "output/cyano_som_endemics.Rdata")

som_endemics(phys_dat = "data/NCOG_sample_log_DNA_meta_2014-2019.csv",
             k.val = 2, surf = TRUE,
             in_data = "data/16s_bacteria_m_euks.Rdata", title = "Heterotrophic Bacteria",
             out_file = "output/bacter_m_euks_som_endemics.Rdata")

som_endemics(phys_dat = "data/NCOG_sample_log_DNA_meta_2014-2019.csv",
             k.val = 2, surf = TRUE,
             in_data = "data/16s_archaea.Rdata", title = "Archaea",
             out_file = "output/arch_som_endemics.Rdata")

som_endemics(phys_dat = "data/NCOG_sample_log_DNA_meta_2014-2019.csv",
             k.val = 2, surf = TRUE,
             in_data = "data/18s_autotrophic_euks.Rdata",
             title = "Photosynthetic Eukaryotic\nProtists",
             out_file = "output/euks_auto_som_endemics.Rdata")

som_endemics(phys_dat = "data/NCOG_sample_log_DNA_meta_2014-2019.csv",
             k.val = 2, surf = TRUE,
             in_data = "data/18s_heterotrophic_euks.Rdata",
             title = "Heterotrophic Eukaryotic\nProtists",
             out_file = "output/euks_hetero_som_endemics.Rdata")

# DCM

som_endemics(phys_dat = "data/NCOG_sample_log_DNA_meta_2014-2019.csv",
             k.val = 4, surf = FALSE,
             in_data = "data/16s_cyanos.Rdata", title = "Cyanobacteria",
             out_file = "output/cyano_som_endemics_dcm.Rdata")

som_endemics(phys_dat = "data/NCOG_sample_log_DNA_meta_2014-2019.csv",
             k.val = 4, surf = FALSE,
             in_data = "data/16s_bacteria_m_euks.Rdata", title = "Heterotrophic Bacteria",
             out_file = "output/bacter_m_euks_som_endemics_dcm.Rdata")

som_endemics(phys_dat = "data/NCOG_sample_log_DNA_meta_2014-2019.csv",
             k.val = 4, surf = FALSE,
             in_data = "data/16s_archaea.Rdata", title = "Archaea",
             out_file = "output/arch_som_endemics_dcm.Rdata")

som_endemics(phys_dat = "data/NCOG_sample_log_DNA_meta_2014-2019.csv",
             k.val = 4, surf = FALSE,
             in_data = "data/18s_autotrophic_euks.Rdata",
             title = "Photosynthetic Eukaryotic\nProtists",
             out_file = "output/euks_auto_som_endemics_dcm.Rdata")

som_endemics(phys_dat = "data/NCOG_sample_log_DNA_meta_2014-2019.csv",
             k.val = 4, surf = FALSE,
             in_data = "data/18s_heterotrophic_euks.Rdata",
             title = "Heterotrophic Eukaryotic\nProtists",
             out_file = "output/euks_hetero_som_endemics_dcm.Rdata")




