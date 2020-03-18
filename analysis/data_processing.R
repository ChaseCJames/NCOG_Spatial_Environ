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

cyanos <- which(split_taxa$B == "D_1__Cyanobacteria" & split_taxa$D != "D_3__Chloroplast")

bacteria_sixteen <- sixteen_s[,which(six_tax_id$type == "bacteria")]
bacteria_copy <- bacteria_sixteen

cyano_sixteen <- sixteen_s[,which(split_taxa$B == "D_1__Cyanobacteria" & split_taxa$D != "D_3__Chloroplast")]
cyano_copy <- cyano_sixteen

pro_sixteen <- sixteen_s[,which(split_taxa$B == "D_1__Cyanobacteria" & substr(split_taxa$F,1,14) == "D_5__Prochloro")]
pro_copy <- pro_sixteen

syne_sixteen <- sixteen_s[,which(split_taxa$B == "D_1__Cyanobacteria" & substr(split_taxa$F,1,12) == "D_5__Synecho")]
syne_copy <- syne_sixteen

bacteria_m_euks_sixteen <- sixteen_s[,which(six_tax_id$type == "bacteria" & split_taxa$A != "D_0__Eukaryota" & split_taxa$B != "D_1__Cyanobacteria" & split_taxa$A == "D_0__Bacteria")]
bacteria_m_euks_copy <- bacteria_m_euks_sixteen

plastid_sixteen <- sixteen_s[,which(six_tax_id$type == "plastid")]
plastid_copy <- plastid_sixteen

# other groups for 16s

archaea_sixteen <- sixteen_s[,which(split_taxa$A == "D_0__Archaea")]
archaea_copy <- archaea_sixteen

flavo_sixteen <- sixteen_s[,which(split_taxa$D == "D_3__Flavobacteriales")]
flavo_copy <- flavo_sixteen

rhodo_sixteen <- sixteen_s[,which(split_taxa$D == "D_3__Rhodobacterales")]
rhodo_copy <- rhodo_sixteen

sar_sixteen <- sixteen_s[,which(split_taxa$D == "D_3__SAR11 clade")]
sar_copy <- sar_sixteen

bact_sums <- rowSums(bacteria_sixteen, na.rm = TRUE)
plas_sums <- rowSums(plastid_sixteen, na.rm = TRUE)
cyano_sums <- rowSums(cyano_sixteen, na.rm = TRUE)
bact_euks <- rowSums(bacteria_m_euks_sixteen, na.rm = TRUE)
archaea_sums <- rowSums(archaea_sixteen, na.rm = TRUE)
flavo_sums <- rowSums(flavo_sixteen, na.rm = TRUE)
rhodo_sums <- rowSums(rhodo_sixteen, na.rm = TRUE)
sar_sums <- rowSums(sar_sixteen, na.rm = TRUE)
pro_sums <- rowSums(pro_sixteen, na.rm = TRUE)
syne_sums <- rowSums(syne_sixteen, na.rm = TRUE)

pro_sums[which(pro_sums == 0)] <- 1
syne_sums[which(syne_sums == 0)] <- 1


bacteria_sixteen <- as.matrix(bacteria_sixteen)
plastid_sixteen <- as.matrix(plastid_sixteen)
cyano_sixteen <- as.matrix(cyano_sixteen)
bacteria_m_euks_sixteen <- as.matrix(bacteria_m_euks_sixteen)
archaea_sixteen <- as.matrix(archaea_sixteen)
flavo_sixteen <- as.matrix(flavo_sixteen)
rhodo_sixteen <- as.matrix(rhodo_sixteen)
sar_sixteen <- as.matrix(sar_sixteen)
pro_sixteen <- as.matrix(pro_sixteen)
syne_sixteen <- as.matrix(syne_sixteen)

for (i in 1:nrow(plastid_sixteen)){
  
  bacteria_sixteen[i,] <- bacteria_sixteen[i,]/bact_sums[i]
  plastid_sixteen[i,] <- plastid_sixteen[i,]/plas_sums[i]
  cyano_sixteen[i,] <- cyano_sixteen[i,]/cyano_sums[i]
  bacteria_m_euks_sixteen[i,] <- bacteria_m_euks_sixteen[i,]/bact_euks[i]
  archaea_sixteen[i,] <- archaea_sixteen[i,]/archaea_sums[i]
  flavo_sixteen[i,] <- flavo_sixteen[i,]/flavo_sums[i]
  rhodo_sixteen[i,] <- rhodo_sixteen[i,]/rhodo_sums[i]
  sar_sixteen[i,] <- sar_sixteen[i,]/sar_sums[i]
  pro_sixteen[i,] <- pro_sixteen[i,]/pro_sums[i]
  syne_sixteen[i,] <- syne_sixteen[i,]/syne_sums[i]
  
}

bacteria_sixteen <- as.data.frame(bacteria_sixteen)
plastid_sixteen <- as.data.frame(plastid_sixteen)
cyano_sixteen <- as.data.frame(cyano_sixteen)
bacteria_m_euks_sixteen <- as.data.frame(bacteria_m_euks_sixteen)
archaea_sixteen <- as.data.frame(archaea_sixteen)
flavo_sixteen <- as.data.frame(flavo_sixteen)
rhodo_sixteen <- as.data.frame(rhodo_sixteen)
sar_sixteen <- as.data.frame(sar_sixteen)
pro_sixteen <- as.data.frame(pro_sixteen)
syne_sixteen <- as.data.frame(syne_sixteen)

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

scaled_inputs <- archaea_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- archaea_copy

save(scaled_inputs, archaea_sixteen, six_tax_id, asv_table, file = "data/16s_archaea.Rdata")

scaled_inputs <- flavo_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- flavo_copy

save(scaled_inputs, flavo_sixteen, six_tax_id, asv_table, file = "data/16s_flavo.Rdata")

scaled_inputs <- flavo_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- flavo_copy

save(scaled_inputs, flavo_sixteen, six_tax_id, asv_table, file = "data/16s_flavo.Rdata")

scaled_inputs <- rhodo_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- rhodo_copy

save(scaled_inputs, rhodo_sixteen, six_tax_id, asv_table, file = "data/16s_rhodo.Rdata")

scaled_inputs <- sar_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- sar_copy

save(scaled_inputs, sar_sixteen, six_tax_id, asv_table, file = "data/16s_sar.Rdata")

scaled_inputs <- pro_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- pro_copy

save(scaled_inputs, pro_sixteen, six_tax_id, asv_table, file = "data/16s_pro.Rdata")

scaled_inputs <- syne_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- syne_copy

save(scaled_inputs, syne_sixteen, six_tax_id, asv_table, file = "data/16s_syne.Rdata")

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
sindins <- which(split_taxa$D == "Syndiniales")
ochrophytes <- which(split_taxa$C == "Ochrophyta")

autotrophs <- c(archaeplastids, haptophytes, cryptophytes,
                acantharea,
                ochrophytes)

eight_auto <- eighteen_s[,autotrophs]

non_auto <- 1:ncol(eighteen_s)
non_auto <- non_auto[-autotrophs]

eight_hetero <- eighteen_s[,c(non_auto, which(split_taxa$C != "Metazoa"))]

diatom_eighteen <- eighteen_s[,which(split_taxa$D == "Bacillariophyta")]

dino_eighteen <- eighteen_s[,dinos_minus_syn]

syndin_eighteen <- eighteen_s[,sindins]

hapto_eighteen <- eighteen_s[,haptophytes]

metazoa_eighteen <- eighteen_s[, which(split_taxa$C == "Metazoa")]

chloro_eighteen <- eighteen_s[,which(split_taxa$C == "Chlorophyta")]

auto_sums <- rowSums(eight_auto, na.rm = TRUE)
hetero_sums <- rowSums(eight_hetero, na.rm = TRUE)
diatom_sums <- rowSums(diatom_eighteen, na.rm = TRUE)
dino_sums <- rowSums(dino_eighteen, na.rm = TRUE)
syndin_sums <- rowSums(syndin_eighteen, na.rm = TRUE)
hapto_sums <- rowSums(hapto_eighteen, na.rm = TRUE)
metazoa_sums <- rowSums(metazoa_eighteen, na.rm = TRUE)
chloro_sums <- rowSums(chloro_eighteen, na.rm = TRUE)

diatom_sums[which(diatom_sums == 0)] <- 1
dino_sums[which(dino_sums == 0)] <- 1
syndin_sums[which(syndin_sums == 0)] <- 1
hapto_sums[which(hapto_sums == 0)] <- 1
metazoa_sums[which(metazoa_sums == 0)] <- 1
chloro_sums[which(chloro_sums == 0)] <- 1

auto_copy <- eight_auto
hetero_copy <- eight_hetero
diatom_copy <- diatom_eighteen
dino_copy <- dino_eighteen
syndin_copy <- syndin_eighteen
hapto_copy <- hapto_eighteen
metazoa_copy <- metazoa_eighteen
chloro_copy <- chloro_eighteen

eight_auto <- as.matrix(eight_auto)
eight_hetero <- as.matrix(eight_hetero)
diatom_eighteen <- as.matrix(diatom_eighteen)
dino_eighteen <- as.matrix(dino_eighteen)
syndin_eighteen <- as.matrix(syndin_eighteen)
hapto_eighteen <- as.matrix(hapto_eighteen)
metazoa_eighteen <- as.matrix(metazoa_eighteen)
chloro_eighteen <- as.matrix(chloro_eighteen)

for (i in 1:nrow(eight_auto)){
  
  eight_auto[i,] <- eight_auto[i,]/auto_sums[i]
  eight_hetero[i,] <- eight_hetero[i,]/hetero_sums[i]
  diatom_eighteen[i,] <- diatom_eighteen[i,]/diatom_sums[i]
  dino_eighteen[i,] <- dino_eighteen[i,]/dino_sums[i]
  syndin_eighteen[i,] <- syndin_eighteen[i,]/syndin_sums[i]
  hapto_eighteen[i,] <- hapto_eighteen[i,]/hapto_sums[i]
  metazoa_eighteen[i,] <- metazoa_eighteen[i,]/metazoa_sums[i]
  chloro_eighteen[i,] <- chloro_eighteen[i,]/chloro_sums[i]
  
}

eight_auto <- as.data.frame(eight_auto)
eight_hetero <- as.data.frame(eight_hetero)
diatom_eighteen <- as.data.frame(diatom_eighteen)
dino_eighteen <- as.data.frame(dino_eighteen)
syndin_eighteen <- as.data.frame(syndin_eighteen)
hapto_eighteen <- as.data.frame(hapto_eighteen)
metazoa_eighteen <- as.data.frame(metazoa_eighteen)
chloro_eighteen <- as.data.frame(chloro_eighteen)

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

# Diatom 18sv9

scaled_inputs <- diatom_eighteen
scaled_inputs <- as.matrix(scaled_inputs)
asv_table <- diatom_copy

save(scaled_inputs, diatom_eighteen, eight_tax_id, asv_table, file = "data/18s_diatom.Rdata")

# Dinoflagellates 18sv9

scaled_inputs <- dino_eighteen
scaled_inputs <- as.matrix(scaled_inputs)
asv_table <- dino_copy

save(scaled_inputs, dino_eighteen, eight_tax_id, asv_table, file = "data/18s_dino.Rdata")

# Syndiniales 18sv9

scaled_inputs <- syndin_eighteen
scaled_inputs <- as.matrix(scaled_inputs)
asv_table <- syndin_copy

save(scaled_inputs, syndin_eighteen, eight_tax_id, asv_table, file = "data/18s_syndin.Rdata")

# Haptophyte 18sv9

scaled_inputs <- hapto_eighteen
scaled_inputs <- as.matrix(scaled_inputs)
asv_table <- hapto_copy

save(scaled_inputs, hapto_eighteen, eight_tax_id, asv_table, file = "data/18s_hapto.Rdata")

# Metazoa 18sv9

scaled_inputs <- metazoa_eighteen
scaled_inputs <- as.matrix(scaled_inputs)
asv_table <- metazoa_copy

save(scaled_inputs, metazoa_eighteen, eight_tax_id, asv_table, file = "data/18s_metazoa.Rdata")

# Chlorophytes 18sv9

scaled_inputs <- chloro_eighteen
scaled_inputs <- as.matrix(scaled_inputs)
asv_table <- chloro_copy

save(scaled_inputs, chloro_eighteen, eight_tax_id, asv_table, file = "data/18s_chloro.Rdata")

# icthyoplankton

# load("data/icthyoplankton.Rdata")
# 
# icthy_mat <- spread(calcofi_filt, key = scientific_name, value = larvae_10m2, fill = 0)

