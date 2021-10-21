library(tidyverse)

meta_pre <- read.csv("data/NCOG_sample_log_DNA_meta_2014-2020.csv")

# pull in new samples

sterivex <- read.csv("data/2014-2016-surface-final-extraction-list.csv")

sterivex <- sterivex %>% filter(Sample.Name != "blank ")

# align with old data

colnames(sterivex) <- c("New.Sample.Name", "Sample.Name", "Cruise",
                        "Order_Occ", "Sta_ID", "Bottle", "Assoc_Bottle",
                        "Depthm", "sample_num", "Plate.Position", "Sterivex_Notes")

sterivex$New.Sample.Name <- gsub("STVX","S",sterivex$New.Sample.Name)

# new ND number

max_nd <- substr(meta_pre$sample_num,3,6) %>% as.numeric() %>% max()

sterivex$sample_num <- paste0("ND",1:nrow(sterivex) + max_nd)

sterivex$Plate.Position <- NULL

# Make variables numeric if need be
sterivex$Cruise <- as.numeric(sterivex$Cruise)
sterivex$Order_Occ <- as.numeric(sterivex$Order_Occ)


# bring in filter metadata

filt <- read.csv("data/NCOG_sample_log_2014_202010_MetaData.csv")

sterivex <- left_join(sterivex, filt[,c(1:23)])

sterivex$Sample.Name <- NULL

colnames(sterivex)[1] <- "Sample.Name"

#align datasets
meta_pre$Bottle <- as.character(meta_pre$Bottle)
sterivex$Depthm <- as.numeric(sterivex$Depthm)
sterivex$Gly.Sample <- NA

new_meta <- bind_rows(meta_pre[,c(1:24)], sterivex)

# Bring in and QC Bottle and ML Data

bottle <- read.csv("data/TI-BottleData.csv")
ml_dat <- read.csv("data/ncd_data.csv")

# bottle first 

# make continous values numeric

new_meta$match_col <- paste(new_meta$Cruise, new_meta$Sta_ID, new_meta$Depthm)
bottle$match_col <- paste(bottle$Cruise, bottle$Sta_ID, bottle$Depthm)

out <- bottle[match(new_meta$match_col, bottle$match_col),c(4:6,11:13,17:28)]

new_meta <- bind_cols(new_meta, out)

bottle_na <- new_meta[which(is.na(new_meta$Distance)),]

new_meta <- new_meta %>% filter(!is.na(Distance))

new_meta$Del_Depth <- 0

new_meta$match_col <- NULL
bottle$match_col <- NULL

# go by bottle next

colnames(bottle)[15] <- "Bottle"
colnames(bottle)[16] <-"Bottle_Depth"
bottle$Bottle <- as.character(bottle$Bottle)

bottle_na$match_col <- paste(bottle_na$Cruise, bottle_na$Sta_ID, bottle_na$Bottle)
bottle$match_col <- paste(bottle$Cruise, bottle$Sta_ID, bottle$Bottle)

out <- bottle[match(bottle_na$match_col, bottle$match_col),c(4:6,11:13,16:28)]

bottle_dat <- bind_cols(bottle_na[,1:25], out)

bottle_dat$Del_Depth <- bottle_dat$Depthm - bottle_dat$Bottle_Depth

meta_join <- bottle_dat %>% filter(!is.na(Distance), abs(Del_Depth) <= 15)

new_meta <- bind_rows(new_meta, meta_join)

missing_vals <- bottle_dat %>% filter(is.na(Distance) | abs(Del_Depth) > 15)

# Now associated bottle

colnames(bottle)[15] <- "Assoc_Bottle"

missing_vals$match_col <- paste(missing_vals$Cruise, missing_vals$Sta_ID, missing_vals$Assoc_Bottle)
bottle$match_col <- paste(bottle$Cruise, bottle$Sta_ID, bottle$Assoc_Bottle)

out <- bottle[match(missing_vals$match_col, bottle$match_col),c(4:6,11:13,16:28)]

bottle_dat <- bind_cols(missing_vals[,1:25], out)

bottle_dat$Del_Depth <- bottle_dat$Depthm - bottle_dat$Bottle_Depth

meta_join <- bottle_dat %>% filter(!is.na(Distance), abs(Del_Depth) <= 15)

new_meta <- bind_rows(new_meta, meta_join)

missing_vals <- bottle_dat %>% filter(is.na(Distance) | abs(Del_Depth) > 15)

new_meta$Bottle_Depth[which(is.na(new_meta$Bottle_Depth))] <- new_meta$Depthm[which(is.na(new_meta$Bottle_Depth))]

# For now there is missing bottle data for 2019-2020, need to check with Rob about this
# will add this data from the previous metadata

# fix Del Depth for 2020

colnames(meta_pre)[32] <- "Del_Depth"

act <- meta_pre$Del_Depth[which(as.numeric(substr(meta_pre$Cruise,1,4)) == 2020)]
del <- meta_pre$CC_Depth[which(as.numeric(substr(meta_pre$Cruise,1,4)) == 2020)]

meta_pre$Del_Depth[which(as.numeric(substr(meta_pre$Cruise,1,4)) == 2020)] <- del
meta_pre$CC_Depth[which(as.numeric(substr(meta_pre$Cruise,1,4)) == 2020)] <- act


missing_vals$match_col <- paste(missing_vals$Cruise, missing_vals$Sta_ID, missing_vals$Depthm)
meta_pre$match_col <- paste(meta_pre$Cruise, meta_pre$Sta_ID, meta_pre$Depthm)

out <- meta_pre[match(missing_vals$match_col, meta_pre$match_col),c(25:30,32,34:44)]

new_year_bottle <- bind_cols(missing_vals[,1:25], out)

# remove oxy sat from meta
new_meta$O2Sat <- NULL

# add bottle depth
new_year_bottle$Bottle_Depth <- new_year_bottle$Depthm + new_year_bottle$Del_Depth

meta_join <- new_year_bottle %>% filter(as.numeric(substr(Cruise,1,4)) > 2018) 

new_meta <- bind_rows(new_meta, meta_join)

# these 11 samples have missing data, Some are ESP samples which don't appear in the bottle data
# others just don't appear in the data, particularly Station 76 80, will bind to data for now

missing_vals <- new_year_bottle %>% filter(as.numeric(substr(Cruise,1,4)) < 2019) 

new_meta <- bind_rows(new_meta, missing_vals)

new_meta <- new_meta %>% arrange(as.numeric(substr(sample_num,3,6)))

# Reorder columns

new_meta <- new_meta[,c(1:9,26:30,10:11,44,12,43,13:25,31:42)]

# add mixed layer and nitracline depth data

colnames(ml_dat)[2] <- "Cruise"

new_meta$match_col <- paste(new_meta$Cruise, new_meta$Sta_ID)
ml_dat$match_col <- paste(ml_dat$Cruise, ml_dat$Sta_ID)

out <- ml_dat[match(new_meta$match_col, ml_dat$match_col),24:25]

new_meta$match_col <- NULL

new_meta <- bind_cols(new_meta, out)

write.csv(new_meta, file = "data/NCOG_sample_log_DNA_stvx_meta_2014-2020.csv", row.names = FALSE)

