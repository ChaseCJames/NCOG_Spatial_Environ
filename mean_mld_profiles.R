library(tidyverse)
library(ggmap)



early_dat <- read.csv("data/ml_data_2014_2018.csv", header = TRUE, stringsAsFactors = FALSE)
late_dat <- read.csv("data/ml_data_2019.csv", header = TRUE, stringsAsFactors = FALSE)

early_dat$Date <- NULL
late_dat$Date <- NULL

late_dat$Cruise <- late_dat$CruiseAlias

combined_dat <- bind_rows(early_dat,late_dat)

combined_dat <- combined_dat[-which(substr(combined_dat$Sta_ID,2,3)<76),]

all_dat <- combined_dat %>% 
  group_by(Sta_ID) %>%
  summarise(lat = mean(DLat_Dec, na.rm = TRUE),
            lon = mean(DLon_Dec, na.rm = TRUE),
            mean_temp = mean(Temp, na.rm = TRUE),
            mean_sal = mean(Sal, na.rm = TRUE),
            mean_dens = mean(Density, na.rm = TRUE),
            mean_chl = mean(Chl, na.rm = TRUE),
            mean_po4 = mean(PO4, na.rm = TRUE),
            mean_no3 = mean(NO3, na.rm = TRUE),
            mean_sio4 = mean(SiO4, na.rm = TRUE),
            mean_mld = mean(MLD_Sigma, na.rm = TRUE),
            mean_ncd = mean(NCDepth, na.rm = TRUE), 
            coeff_temp = sd(Temp, na.rm = TRUE)/mean(Temp, na.rm = TRUE),
            coeff_sal = sd(Sal, na.rm = TRUE)/mean(Sal, na.rm = TRUE),
            coeff_dens = sd(Density, na.rm = TRUE)/mean(Density, na.rm = TRUE),
            coeff_chl = sd(Chl, na.rm = TRUE)/mean(Chl, na.rm = TRUE),
            coeff_po4 = sd(PO4, na.rm = TRUE)/mean(PO4, na.rm = TRUE),
            coeff_no3 = sd(NO3, na.rm = TRUE)/mean(NO3, na.rm = TRUE),
            coeff_sio4 = sd(SiO4, na.rm = TRUE)/mean(SiO4, na.rm = TRUE),
            coeff_mld = sd(MLD_Sigma, na.rm = TRUE)/mean(MLD_Sigma, na.rm = TRUE),
            coeff_ncd = sd(NCDepth, na.rm = TRUE)/mean(NCDepth, na.rm = TRUE))

early_phase <- combined_dat %>% 
  filter(Year < 2017) %>%
  group_by(Sta_ID) %>%
  summarise(lat = mean(DLat_Dec, na.rm = TRUE),
            lon = mean(DLon_Dec, na.rm = TRUE),
            mean_temp = mean(Temp, na.rm = TRUE),
            mean_sal = mean(Sal, na.rm = TRUE),
            mean_dens = mean(Density, na.rm = TRUE),
            mean_chl = mean(Chl, na.rm = TRUE),
            mean_po4 = mean(PO4, na.rm = TRUE),
            mean_no3 = mean(NO3, na.rm = TRUE),
            mean_sio4 = mean(SiO4, na.rm = TRUE),
            mean_mld = mean(MLD_Sigma, na.rm = TRUE),
            mean_ncd = mean(NCDepth, na.rm = TRUE), 
            coeff_temp = sd(Temp, na.rm = TRUE)/mean(Temp, na.rm = TRUE),
            coeff_sal = sd(Sal, na.rm = TRUE)/mean(Sal, na.rm = TRUE),
            coeff_dens = sd(Density, na.rm = TRUE)/mean(Density, na.rm = TRUE),
            coeff_chl = sd(Chl, na.rm = TRUE)/mean(Chl, na.rm = TRUE),
            coeff_po4 = sd(PO4, na.rm = TRUE)/mean(PO4, na.rm = TRUE),
            coeff_no3 = sd(NO3, na.rm = TRUE)/mean(NO3, na.rm = TRUE),
            coeff_sio4 = sd(SiO4, na.rm = TRUE)/mean(SiO4, na.rm = TRUE),
            coeff_mld = sd(MLD_Sigma, na.rm = TRUE)/mean(MLD_Sigma, na.rm = TRUE),
            coeff_ncd = sd(NCDepth, na.rm = TRUE)/mean(NCDepth, na.rm = TRUE))

late_phase <- combined_dat %>% 
  filter(Year > 2017 & Year < 2019) %>%
  group_by(Sta_ID) %>%
  summarise(lat = mean(DLat_Dec, na.rm = TRUE),
            lon = mean(DLon_Dec, na.rm = TRUE),
            mean_temp = mean(Temp, na.rm = TRUE),
            mean_sal = mean(Sal, na.rm = TRUE),
            mean_dens = mean(Density, na.rm = TRUE),
            mean_chl = mean(Chl, na.rm = TRUE),
            mean_po4 = mean(PO4, na.rm = TRUE),
            mean_no3 = mean(NO3, na.rm = TRUE),
            mean_sio4 = mean(SiO4, na.rm = TRUE),
            mean_mld = mean(MLD_Sigma, na.rm = TRUE),
            mean_ncd = mean(NCDepth, na.rm = TRUE), 
            coeff_temp = sd(Temp, na.rm = TRUE)/mean(Temp, na.rm = TRUE),
            coeff_sal = sd(Sal, na.rm = TRUE)/mean(Sal, na.rm = TRUE),
            coeff_dens = sd(Density, na.rm = TRUE)/mean(Density, na.rm = TRUE),
            coeff_chl = sd(Chl, na.rm = TRUE)/mean(Chl, na.rm = TRUE),
            coeff_po4 = sd(PO4, na.rm = TRUE)/mean(PO4, na.rm = TRUE),
            coeff_no3 = sd(NO3, na.rm = TRUE)/mean(NO3, na.rm = TRUE),
            coeff_sio4 = sd(SiO4, na.rm = TRUE)/mean(SiO4, na.rm = TRUE),
            coeff_mld = sd(MLD_Sigma, na.rm = TRUE)/mean(MLD_Sigma, na.rm = TRUE),
            coeff_ncd = sd(NCDepth, na.rm = TRUE)/mean(NCDepth, na.rm = TRUE))


save(all_dat, early_phase, late_phase, file = "output/mld_mean_profiles.Rdata")

