library(SOMbrero)
library(tidyverse)
library(oce)
library(rworldmap)
library(ggplot2)
library(ggmap)
library(cowplot)
library(scales)
library(leaps)
library(MASS)
library(randomForest)
library(grid)
library(gridExtra)
library(lubridate)
library(vegan)
library(spatialEco)
library(geosphere)
library(viridis)

########## Data Processing ############

excludes <- read.delim("data/exclude_samples.txt", header = FALSE)

####### Splitting 16s plastids and everything else #####

# read in 16s
sixteen_s <- read.csv("data/16s_asv_count_tax_final.csv", stringsAsFactors = FALSE)

six_id_names <- sixteen_s$Feature.ID

six_tax_id <- sixteen_s[,c(1,(ncol(sixteen_s)-3):ncol(sixteen_s))]

sixteen_s <- sixteen_s[,-c(1,(ncol(sixteen_s)-30):ncol(sixteen_s))]

sixteen_s <- apply(sixteen_s, 2, as.numeric)

sixteen_s <- t(sixteen_s)

sixteen_s <- as.data.frame(sixteen_s)

six_tp <- rownames(sixteen_s)

colnames(sixteen_s) <- six_id_names

# remove bad samples
sixteen_s <- sixteen_s[-which(!is.na(match(rownames(sixteen_s), paste0("X",excludes$V1)))),]

# split platids and rest of 16s

id_vector <- vector()

for (i in 1:nrow(six_tax_id)) {
val <- which.max(c(six_tax_id$Silva_Confidence[[i]], six_tax_id$Phytoref_Confidence[[i]]))
if(val == 1){id_vector[i] = "bacteria"}else{id_vector[i] = "plastid"}
}

six_tax_id$type <- id_vector

# split taxa

taxas <- six_tax_id[which(!is.na(match(six_tax_id$Feature.ID, colnames(sixteen_s)))),1:3]

colnames(taxas) <- c("Feature.ID", "Taxon", "Confidence")

split_taxa <- separate(taxas, Taxon, sep = ";", into = c("A","B","C", "D", "E", "F", "G", "H", "I"))

# get cyanos and everything but euks

cyanos <- which(split_taxa$B == "D_1__Cyanobacteria")

bacteria_sixteen <- sixteen_s[,which(six_tax_id$type == "bacteria")]
bacteria_copy <- bacteria_sixteen

cyano_sixteen <- sixteen_s[,which(six_tax_id$type == "bacteria" & split_taxa$B == "D_1__Cyanobacteria")]
cyano_copy <- cyano_sixteen

bacteria_m_euks_sixteen <- sixteen_s[,which(six_tax_id$type == "bacteria" & split_taxa$A != "D_0__Eukaryota" & split_taxa$B != "D_1__Cyanobacteria")]
bacteria_m_euks_copy <- bacteria_m_euks_sixteen

plastid_sixteen <- sixteen_s[,which(six_tax_id$type == "plastid")]
plastid_copy <- plastid_sixteen

bact_sums <- rowSums(bacteria_sixteen, na.rm = TRUE)
plas_sums <- rowSums(plastid_sixteen, na.rm = TRUE)
cyano_sums <- rowSums(cyano_sixteen, na.rm = TRUE)
bact_euks <- rowSums(bacteria_m_euks_sixteen, na.rm = TRUE)


bacteria_sixteen <- as.matrix(bacteria_sixteen)
plastid_sixteen <- as.matrix(plastid_sixteen)
cyano_sixteen <- as.matrix(cyano_sixteen)
bacteria_m_euks_sixteen <- as.matrix(bacteria_m_euks_sixteen)

for (i in 1:nrow(plastid_sixteen)){
  
  bacteria_sixteen[i,] <- bacteria_sixteen[i,]/bact_sums[i]
  plastid_sixteen[i,] <- plastid_sixteen[i,]/plas_sums[i]
  cyano_sixteen[i,] <- cyano_sixteen[i,]/cyano_sums[i]
  bacteria_m_euks_sixteen[i,] <- bacteria_m_euks_sixteen[i,]/bact_euks[i]
  
}

bacteria_sixteen <- as.data.frame(bacteria_sixteen)
plastid_sixteen <- as.data.frame(plastid_sixteen)
cyano_sixteen <- as.data.frame(cyano_sixteen)
bacteria_m_euks_sixteen <- as.data.frame(bacteria_m_euks_sixteen)

scaled_inputs <- bacteria_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- bacteria_copy

save(scaled_inputs, bacteria_sixteen, six_tax_id, asv_table, file = "data/16s_bacteria.Rdata")

scaled_inputs <- plastid_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- plastid_copy

save(scaled_inputs, plastid_sixteen, six_tax_id, asv_table, file = "data/16s_plastids.Rdata")

scaled_inputs <- cyano_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- cyano_copy

save(scaled_inputs, cyano_sixteen, six_tax_id, asv_table, file = "data/16s_cyanos.Rdata")

scaled_inputs <- bacteria_m_euks_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- bacteria_m_euks_copy

save(scaled_inputs, bacteria_m_euks_sixteen, six_tax_id, asv_table, file = "data/16s_bacteria_m_euks.Rdata")

###### Splitting 18sv9 between autotrophs and everything else #####

eighteen_s <- read.csv("data/18sV9_asv_count_tax_final.csv", stringsAsFactors = FALSE)

eight_id_names <- eighteen_s$Feature.ID

eight_tax_id <- eighteen_s[,c(1,(ncol(eighteen_s)-3):ncol(eighteen_s))]

eighteen_s <- eighteen_s[,-c(1,(ncol(eighteen_s)-3):ncol(eighteen_s))]

eighteen_s <- apply(eighteen_s, 2, as.numeric)

eighteen_s <- t(eighteen_s)

eighteen_s <- as.data.frame(eighteen_s)

eight_tp <- rownames(eighteen_s)

colnames(eighteen_s) <- eight_id_names

# remove mocks from data and rows with no reads
eighteen_s <- eighteen_s[-c(479:491,669:680,873:878),]

# remove bad samples
eighteen_s <- eighteen_s[-which(!is.na(match(rownames(eighteen_s), paste0("X",excludes$V1)))),]

# sort

id_vector <- vector()

for (i in 1:nrow(eight_tax_id)) {
  val <- which.max(c(eight_tax_id$Confidence_Silva[[i]], eight_tax_id$Confidence_PR2[[i]]))
  if(val == 1){a = 4}else{a = 2}
  id_vector[i] <- eight_tax_id[i,a]
}

eight_tax_id$Taxon <- id_vector

split_taxa <- separate(eight_tax_id, Taxon_PR2, sep = ";", into = c("A","B","C", "D", "E", "F", "G", "H", "I"))
split_silva <- separate(eight_tax_id, Taxon_Silva, sep = ";", into = c("A","B","C", "D", "E", "F", "G", "H", "I"))

euks <- which(split_silva$A == "D_0__Eukaryota" | split_silva$A == "Unassigned")

archaeplastids <- which(split_taxa$B == "Archaeplastida")
haptophytes <- which(split_taxa$C == "Haptophyta")
cryptophytes <- which(split_taxa$C == "Cryptophyta")
euglenids <- which(split_taxa$E == "Euglenida")
acantharea <- which(split_taxa$D == "Acantharea")
dinos_minus_syn <- which(split_taxa$C == "Dinoflagellata" &
                           split_taxa$D != "Syndiniales")
ochrophytes <- which(split_taxa$C == "Ochrophyta")

autotrophs <- c(archaeplastids, haptophytes, cryptophytes,
                euglenids, acantharea, dinos_minus_syn,
                ochrophytes)

eight_auto <- eighteen_s[,autotrophs]

euks_minus_auto <- euks[which(is.na(match(euks, autotrophs)))]

eight_hetero <- eighteen_s[,euks_minus_auto]

auto_sums <- rowSums(eight_auto, na.rm = TRUE)
hetero_sums <- rowSums(eight_hetero, na.rm = TRUE)

eight_auto <- eight_auto[-which(auto_sums == 0),]
eight_hetero <- eight_hetero[-which(hetero_sums == 0),]

auto_copy <- eight_auto
hetero_copy <- eight_hetero

auto_sums <- auto_sums[-which(auto_sums == 0)]
hetero_sums <- hetero_sums[-which(hetero_sums == 0)]

eight_auto <- as.matrix(eight_auto)
eight_hetero <- as.matrix(eight_hetero)

for (i in 1:nrow(eight_auto)){
  
  eight_auto[i,] <- eight_auto[i,]/auto_sums[i]
  eight_hetero[i,] <- eight_hetero[i,]/hetero_sums[i]
  
}

eight_auto <- as.data.frame(eight_auto)
eight_hetero <- as.data.frame(eight_hetero)

# autotrophic 18sv9

scaled_inputs <- eight_auto

scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- auto_copy

save(scaled_inputs, eight_auto, eight_tax_id, asv_table, file = "data/18s_autotrophic_euks.Rdata")

# heterotrophic 18sv9

scaled_inputs <- eight_hetero

scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- hetero_copy

save(scaled_inputs, eight_hetero, eight_tax_id, asv_table, file = "data/18s_heterotrophic_euks.Rdata")

# icthyoplankton

load("data/icthyoplankton.Rdata")

icthy_mat <- spread(calcofi_filt, key = scientific_name, value = larvae_10m2, fill = 0)

