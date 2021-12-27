library(tidyverse)
library(lubridate)

meta_edit <- read.csv("data/meta_t_metadata/NCOG_sample_log_RNA_meta_2014-2018.csv")

data_2020 <- read.csv("data/meta_t_metadata/ML Data 1 Stations to CC202010.csv")
data_2019 <- read.csv("data/meta_t_metadata/ML Data 1 Stations 2019.csv")

extra_2020 <- read.csv("data/meta_t_metadata/NCOG_sample_log_2014_202010_MetaData.csv")
bottle_2020 <- read.csv("data/meta_t_metadata/TI - BottleData to CC202010.csv")

bottle_2019 <- read.csv("data/meta_t_metadata/TI - BottleData 2019.csv")

# Work on incomplete 2019 data first

meta_19 <- meta_edit %>% filter(substr(Cruise,1,4) == 2019)

# columns missing (IntChl, IntC14, MLD_Sigma, NCDepth)

meta_19 <- meta_19[,-c(45:48)]

data_2019 <- data_2019[,c("CruiseAlias","Sta_ID","IntChl", "IntC14", "MLD_Sigma", "NCDepth")]
colnames(data_2019)[1] <- "Cruise"

meta_19 <- left_join(meta_19, data_2019, by = c("Cruise", "Sta_ID"), keep = FALSE)

# get delta depths

meta_19$DelDepth[which(meta_19$Assoc_Bottle == "")] <- 0


meta_19$DelDepth[which(meta_19$Assoc_Bottle != "")] <-  meta_19$Depthm[which(meta_19$Assoc_Bottle != "")] -
  bottle_2019$Depthm[match(paste(meta_19$Cruise[which(meta_19$Assoc_Bottle != "")],
                                 meta_19$Sta_ID[which(meta_19$Assoc_Bottle != "")],
                                 meta_19$Bottle[which(meta_19$Assoc_Bottle != "")]),
                           paste(bottle_2019$Cruise, bottle_2019$Sta_ID, bottle_2019$BtlNum))]


meta_19 <- meta_19[,-c(25,26)]

bottle_2019_dt <- bottle_2019[,c(2,9,5,6)]

bottle_2019_dt <- bottle_2019_dt[!duplicated(bottle_2019_dt),]

meta_19 <- left_join(meta_19, bottle_2019_dt, by = c("Cruise", "Sta_ID"), keep = FALSE)

meta_19$NoCCSamples <- meta_19$NCOG_DNA + meta_19$NCOG_RNA

# Now work on 2020 data

meta_20 <- meta_edit %>% filter(is.na(Cruise))
meta_edit <- meta_edit %>% filter(!is.na(Cruise))
meta_edit <- meta_edit %>% filter(substr(Cruise,1,4) != 2019)

# Add data info from sample name

meta_20$Cruise <- substr(meta_20$Sample.Name,1,6)
meta_20$Sta_ID <- gsub("_"," ",substr(meta_20$Sample.Name,8,18))
meta_20$Depthm <- substr(meta_20$Sample.Name,20,26) %>% as.numeric()
meta_20$Cast_Type <- "Prodo"

meta_20$Cardinal_Sta <- meta_edit$Cardinal_Sta[match(meta_20$Sta_ID, meta_edit$Sta_ID)]
meta_20$Cardinal_Sta[which(is.na(meta_20$Cardinal_Sta))] <- FALSE

# Merge some metadata

meta_20 <- meta_20[,c(1:10,14)]

extra_2020 <- extra_2020[,c(1,8:10,12:20,26:43)]

meta_20 <- left_join(meta_20, extra_2020, by = "Sample.Name", keep = FALSE)

colnames(meta_20)[32:41] <- c("T_degC", "Salnty", "STheta", "O2ml_L",
                              "PO4ug", "SiO3ug", "NO3ug", "NH3ug",
                              "ChlorA", "Phaeop")

# add grouped data (Int, ML, NC)

data_2020 <- data_2020[,c("CruiseAlias","Sta_ID","IntChl", "IntC14", "MLD_Sigma", "NCDepth", "Distance", "Code_CCE")]
colnames(data_2020)[1] <- "Cruise"
data_2020$Cruise <- data_2020$Cruise %>% as.character()

meta_20 <- left_join(meta_20, data_2020, by = c("Cruise", "Sta_ID"), keep = FALSE)

meta_20$Spike1 <- NA
meta_20$Spike8 <- NA

bottle_2020 <- bottle_2020[,c("Cruise", "Sta_ID", "BtlNum", "RecInd")]
colnames(bottle_2020)[3] <- "Bottle"
bottle_2020$Cruise <- as.character(bottle_2020$Cruise)
bottle_2020$Bottle <- as.character(bottle_2020$Bottle)

meta_20 <- left_join(meta_20, bottle_2020, by = c("Cruise", "Sta_ID", "Bottle"), keep = FALSE)

# make columns match
meta_20$NCOG_DNA[which(meta_20$NCOG_DNA == 1)] <- TRUE
meta_20$NCOG_DNA[which(meta_20$NCOG_DNA == 0)] <- FALSE

meta_20$NCOG_RNA[which(meta_20$NCOG_RNA == 1)] <- TRUE
meta_20$NCOG_RNA[which(meta_20$NCOG_RNA == 0)] <- FALSE

meta_20$Cruise <- as.numeric(meta_20$Cruise)
meta_20$NCOG_DNA <- as.logical(meta_20$NCOG_DNA)
meta_20$NCOG_RNA <- as.logical(meta_20$NCOG_RNA)
meta_20$Filt_End <- as.numeric(meta_20$Filt_End)

meta_20$Time <- substr(meta_20$Time,12,16)

# join all data

meta_update <- bind_rows(meta_edit, meta_19)
meta_update <- bind_rows(meta_update, meta_20)

# depth categories

meta_update$depth_category <- "No Category"
meta_update$depth_category[which(meta_update$Depthm <= 13)] <- "Surface"
meta_update$depth_category[which(meta_update$Depthm > 13 & meta_update$Depthm < 150)] <- "DCM"
meta_update$depth_category[which(meta_update$Depthm > 150 )] <- "Deep"


write.csv(meta_update, file = "data/meta_t_metadata/NCOG_sample_log_RNA_meta_2014-2020.csv")

