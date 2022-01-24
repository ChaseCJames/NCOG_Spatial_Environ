library(tidyverse)
library(lubridate)
library(vegan)

# load in data

eight_s <- read.csv("data/NCOG_18sV9_asv_count_tax_S.csv", header = TRUE)
six_s <- read.csv("data/NCOG_21_16S_redo2_asv_count_tax.csv", header = TRUE)
metadata <- read.csv("data/NCOG_sample_log_DNA_stvx_meta_2014-2020.csv", header = TRUE)

# remove internal spikes

# 16S

six_remove <- c("d__Bacteria; p__Deinococcota; c__Deinococci; o__Thermales; f__Thermaceae; g__Thermus; s__Thermus_thermophilus",
                "d__Bacteria; p__Deinococcota; c__Deinococci; o__Deinococcales; f__Deinococcaceae; g__Deinococcus; s__Deinococcus_radiodurans",
                "d__Bacteria; p__Firmicutes; c__Clostridia; o__Lachnospirales; f__Lachnospiraceae; g__Blautia")

six_s <- six_s %>% filter(!silva_Taxon %in% six_remove)

# 18s

eight_remove <- "Eukaryota;Opisthokonta;Fungi;Ascomycota;Taphrinomycotina;Schizosaccharomycetes;Schizosaccharomyces;Schizosaccharomyces_pombe;"

eight_s <- eight_s %>% filter(!pr2_Taxon %in% eight_remove)

# pull out taxonomy data

eight_tax_id <- eight_s[,c(1,1538:1543)]
six_tax_id <- six_s[,c(1,1071:1072)]

# remove these columns

eight_s <- eight_s[,-c(1538:1543)]
six_s <- six_s[,-c(1071:1072)]

# transform data
eight_feature <- eight_s$Feature.ID
six_feature <- six_s$Feature.ID

eight_s$Feature.ID <- NULL
six_s$Feature.ID <- NULL

eight_samples <- colnames(eight_s)
six_samples <- colnames(six_s)

# match ends for six samples

six_samples <- gsub("16S_S2", "S", six_samples)

eight_s <- eight_s %>% t() %>% as.data.frame()
six_s <- six_s %>% t() %>% as.data.frame()

# filter to samples with only metadata (only sterivex)

metadata <- metadata %>% filter(as.numeric(substr(sample_num,3,6)) > 478)

eight_s$sample <- eight_samples
six_s$sample <- six_samples

eight_s <- eight_s %>% filter(sample %in% paste0("X",metadata$Sample.Name))
six_s <- six_s %>% filter(sample %in% paste0("X",metadata$Sample.Name))

# align 18Sv9 and 16S data (missing 16s sample)
eight_s <- eight_s[-which(!eight_s$sample %in% six_s$sample),]

eight_s$sample <- NULL
six_s$sample <- NULL

# label columns

colnames(eight_s) <- eight_feature
colnames(six_s) <- six_feature

# fix six s rows

rownames(six_s) <- gsub("16S_S2", "S", rownames(six_s))

# for 18Sv9 remove prokaryotes

prok <- eight_tax_id[which(eight_tax_id$silva_Confidence > eight_tax_id$pr2_Confidence),]
prok_id <- prok$Feature.ID[which(vapply(strsplit(prok$silva_Taxon, ";", fixed = TRUE), "[", "", 1) != "d__Eukaryota")]

eight_s <- eight_s[,-c(which(!is.na(match(colnames(eight_s),prok_id))))]

# for 16S remove chloroplasts and euks

euks <- six_tax_id$Feature.ID[which(vapply(strsplit(six_tax_id$silva_Taxon, ";", fixed = TRUE), "[", "", 1) == "d__Eukaryota")]
chloro <- six_tax_id$Feature.ID[which(vapply(strsplit(six_tax_id$silva_Taxon, ";", fixed = TRUE), "[", "", 4) == " o__Chloroplast")]

six_s <- six_s[,-c(which(!is.na(match(colnames(six_s),c(euks,chloro)))))]

# remove ASVs that don't appear in this subset

eight_s <- eight_s[,-which(colSums(eight_s) == 0)]
six_s <- six_s[,-which(colSums(six_s) == 0)]

# remove asvs with batch issues
# source("analysis/filter_test.R")
# 
# six_issues <- which(colnames(six_s) %in% unique(six_cruise$hash))
# eight_issues <- which(colnames(eight_s) %in% unique(eight_cruise$hash))
# 
# six_s <- six_s[,-six_issues]
# eight_s <- eight_s[,-eight_issues]

six_s$sample <- NULL
eight_s$sample <- NULL

# set rarefy level
rare_level <- 17000

# remove samples and re-align data

eight_s <- eight_s[-which(rowSums(eight_s) < rare_level),]
six_s <- six_s[-which(rowSums(six_s) < rare_level),]

eight_s <- eight_s[-which(!rownames(eight_s) %in% rownames(six_s)),]
six_s <- six_s[-which(!rownames(six_s) %in% rownames(eight_s)),]

#rarefy to level
set.seed(961)
eight_rare <- rrarefy(eight_s, rare_level)
set.seed(81)
six_rare <- rrarefy(six_s, rare_level)

# convert to data frames

eight_rare <- eight_rare %>% as.data.frame()
six_rare <- six_rare %>% as.data.frame()

# remove asvs with no data

eight_rare <- eight_rare[,-which(colSums(eight_rare) == 0)]
six_rare <- six_rare[,-which(colSums(six_rare) == 0)]

# split taxonomy data
split_taxa_six <- separate(six_tax_id, silva_Taxon, sep = "; ", into = c("A","B","C", "D", "E", "F", "G"))
split_taxa_eight <- separate(eight_tax_id, pr2_Taxon, sep = ";", into = c("A","B","C", "D", "E", "F", "G", "H"))

split_taxa_six <- split_taxa_six %>% filter(Feature.ID %in% colnames(six_rare))
split_taxa_eight <- split_taxa_eight %>% filter(Feature.ID %in% colnames(eight_rare))

##### 16S pull out groups #####

# Major Groups
cyanobacteria <- which(split_taxa_six$B == "p__Cyanobacteria")
bacteria <- which(split_taxa_six$A == "d__Bacteria" & split_taxa_six$B != "p__Cyanobacteria")
archaea <- which(split_taxa_six$A == "d__Archaea")

# Minor Groups
flavobacteriales <- which(split_taxa_six$D == "o__Flavobacteriales")
rhodobacterales <- which(split_taxa_six$D == "o__Rhodobacterales")
sar_11 <- which(split_taxa_six$D == "o__SAR11_clade")
prochlorococcus <- which(split_taxa_six$B == "p__Cyanobacteria" & substr(split_taxa_six$F,1,12) == "g__Prochloro")
synechococcus <- which(split_taxa_six$B == "p__Cyanobacteria" & substr(split_taxa_six$F,1,10) == "g__Synecho")

# get tables
cyano_sixteen <- six_rare[,cyanobacteria]
bacteria_m_euks_sixteen <- six_rare[,bacteria]
archaea_sixteen <- six_rare[,archaea]

flavo_sixteen <- six_rare[,flavobacteriales]
rhodo_sixteen <- six_rare[,rhodobacterales]
sar_sixteen <- six_rare[,sar_11]
pro_sixteen <- six_rare[,prochlorococcus]
syne_sixteen <- six_rare[,synechococcus]

# make copies

cyano_copy <- cyano_sixteen
bacteria_m_euks_copy <- bacteria_m_euks_sixteen
archaea_copy <- archaea_sixteen

flavo_copy <- flavo_sixteen
rhodo_copy <- rhodo_sixteen
sar_copy <- sar_sixteen
pro_copy <- pro_sixteen
syne_copy <- syne_sixteen

# find library size per group

cyano_sums <- rowSums(cyano_sixteen, na.rm = TRUE)
bact_euks <- rowSums(bacteria_m_euks_sixteen, na.rm = TRUE)
archaea_sums <- rowSums(archaea_sixteen, na.rm = TRUE)
flavo_sums <- rowSums(flavo_sixteen, na.rm = TRUE)
rhodo_sums <- rowSums(rhodo_sixteen, na.rm = TRUE)
sar_sums <- rowSums(sar_sixteen, na.rm = TRUE)
pro_sums <- rowSums(pro_sixteen, na.rm = TRUE)
syne_sums <- rowSums(syne_sixteen, na.rm = TRUE)

archaea_sums[which(archaea_sums == 0)] <- 1
pro_sums[which(pro_sums == 0)] <- 1
syne_sums[which(syne_sums == 0)] <- 1
cyano_sums[which(cyano_sums == 0)] <- 1

# convert to matricies

cyano_sixteen <- as.matrix(cyano_sixteen)
bacteria_m_euks_sixteen <- as.matrix(bacteria_m_euks_sixteen)
archaea_sixteen <- as.matrix(archaea_sixteen)
flavo_sixteen <- as.matrix(flavo_sixteen)
rhodo_sixteen <- as.matrix(rhodo_sixteen)
sar_sixteen <- as.matrix(sar_sixteen)
pro_sixteen <- as.matrix(pro_sixteen)
syne_sixteen <- as.matrix(syne_sixteen)

for (i in 1:nrow(cyano_sixteen)){
  
  cyano_sixteen[i,] <- cyano_sixteen[i,]/cyano_sums[i]
  bacteria_m_euks_sixteen[i,] <- bacteria_m_euks_sixteen[i,]/bact_euks[i]
  archaea_sixteen[i,] <- archaea_sixteen[i,]/archaea_sums[i]
  flavo_sixteen[i,] <- flavo_sixteen[i,]/flavo_sums[i]
  rhodo_sixteen[i,] <- rhodo_sixteen[i,]/rhodo_sums[i]
  sar_sixteen[i,] <- sar_sixteen[i,]/sar_sums[i]
  pro_sixteen[i,] <- pro_sixteen[i,]/pro_sums[i]
  syne_sixteen[i,] <- syne_sixteen[i,]/syne_sums[i]
  
}

cyano_sixteen <- as.data.frame(cyano_sixteen)
bacteria_m_euks_sixteen <- as.data.frame(bacteria_m_euks_sixteen)
archaea_sixteen <- as.data.frame(archaea_sixteen)
flavo_sixteen <- as.data.frame(flavo_sixteen)
rhodo_sixteen <- as.data.frame(rhodo_sixteen)
sar_sixteen <- as.data.frame(sar_sixteen)
pro_sixteen <- as.data.frame(pro_sixteen)
syne_sixteen <- as.data.frame(syne_sixteen)

# save data with "_S" to denote new data

scaled_inputs <- cyano_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- cyano_copy

save(scaled_inputs, cyano_sixteen, six_tax_id, asv_table, file = "data/16s_cyanos_S.Rdata")

scaled_inputs <- bacteria_m_euks_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- bacteria_m_euks_copy

save(scaled_inputs, bacteria_m_euks_sixteen, six_tax_id, asv_table, file = "data/16s_bacteria_m_euks_S.Rdata")

scaled_inputs <- archaea_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- archaea_copy

save(scaled_inputs, archaea_sixteen, six_tax_id, asv_table, file = "data/16s_archaea_S.Rdata")

scaled_inputs <- flavo_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- flavo_copy

save(scaled_inputs, flavo_sixteen, six_tax_id, asv_table, file = "data/16s_flavo_S.Rdata")

scaled_inputs <- rhodo_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- rhodo_copy

save(scaled_inputs, rhodo_sixteen, six_tax_id, asv_table, file = "data/16s_rhodo_S.Rdata")

scaled_inputs <- sar_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- sar_copy

save(scaled_inputs, sar_sixteen, six_tax_id, asv_table, file = "data/16s_sar_S.Rdata")

scaled_inputs <- pro_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- pro_copy

save(scaled_inputs, pro_sixteen, six_tax_id, asv_table, file = "data/16s_pro_S.Rdata")

scaled_inputs <- syne_sixteen
scaled_inputs <- as.matrix(scaled_inputs)

asv_table <- syne_copy

save(scaled_inputs, syne_sixteen, six_tax_id, asv_table, file = "data/16s_syne_S.Rdata")

save(six_rare, six_tax_id, file = "data/16s_all_S.Rdata")

##### 18Sv9 pull out groups #####

# Minor Groups
diatoms <- which(split_taxa_eight$D == "Bacillariophyta")
dinos_minus_syn <- which(split_taxa_eight$C == "Dinoflagellata" &
                           split_taxa_eight$D != "Syndiniales" & !is.na(split_taxa_eight$D))
sindins <- which(split_taxa_eight$D == "Syndiniales")
haptophytes <- which(split_taxa_eight$C == "Haptophyta")
chlorophytes <- which(split_taxa_eight$C == "Chlorophyta")
metazoa <- which(split_taxa_eight$C == "Metazoa")

# groups for sorting autotrophs
cryptophytes <- which(split_taxa_eight$C == "Cryptophyta")
ochrophytes <- which(split_taxa_eight$C == "Ochrophyta")
cercozoa <- which(split_taxa_eight$C == "Cercozoa" & 
                    split_taxa_eight$D != "Filosa-Sarcomonadea")

autotrophs <- c(chlorophytes, dinos_minus_syn, cryptophytes,
                haptophytes, ochrophytes, cercozoa)

non_auto <- 1:ncol(eight_rare)
non_auto <- non_auto[-autotrophs]

non_auto <- non_auto[which(!non_auto %in% metazoa)]

# get tables

eight_auto <- eight_rare[,autotrophs]
eight_hetero <- eight_rare[,non_auto]

diatom_eighteen <- eight_rare[,diatoms]
dino_eighteen <- eight_rare[,dinos_minus_syn]
syndin_eighteen <- eight_rare[,sindins]
hapto_eighteen <- eight_rare[,haptophytes]
chloro_eighteen <- eight_rare[,chlorophytes]
metazoa_eighteen <- eight_rare[,metazoa]

# make copies

auto_copy <- eight_auto
hetero_copy <- eight_hetero

diatom_copy <- diatom_eighteen
dino_copy <- dino_eighteen
syndin_copy <- syndin_eighteen
hapto_copy <- hapto_eighteen
chloro_copy <- chloro_eighteen
metazoa_copy <- metazoa_eighteen

auto_sums <- rowSums(eight_auto, na.rm = TRUE)
hetero_sums <- rowSums(eight_hetero, na.rm = TRUE)
diatom_sums <- rowSums(diatom_eighteen, na.rm = TRUE)
dino_sums <- rowSums(dino_eighteen, na.rm = TRUE)
syndin_sums <- rowSums(syndin_eighteen, na.rm = TRUE)
hapto_sums <- rowSums(hapto_eighteen, na.rm = TRUE)
chloro_sums <- rowSums(chloro_eighteen, na.rm = TRUE)
metazoa_sums <- rowSums(metazoa_eighteen, na.rm = TRUE)

diatom_sums[which(diatom_sums == 0)] <- 1
hapto_sums[which(hapto_sums == 0)] <- 1
chloro_sums[which(chloro_sums == 0)] <- 1

# convert to matrices

eight_auto <- as.matrix(eight_auto)
eight_hetero <- as.matrix(eight_hetero)
diatom_eighteen <- as.matrix(diatom_eighteen)
dino_eighteen <- as.matrix(dino_eighteen)
syndin_eighteen <- as.matrix(syndin_eighteen)
hapto_eighteen <- as.matrix(hapto_eighteen)
chloro_eighteen <- as.matrix(chloro_eighteen)
metazoa_eighteen <- as.matrix(metazoa_eighteen)

for (i in 1:nrow(eight_auto)){

  eight_auto[i,] <- eight_auto[i,]/auto_sums[i]
  eight_hetero[i,] <- eight_hetero[i,]/hetero_sums[i]
  diatom_eighteen[i,] <- diatom_eighteen[i,]/diatom_sums[i]
  dino_eighteen[i,] <- dino_eighteen[i,]/dino_sums[i]
  syndin_eighteen[i,] <- syndin_eighteen[i,]/syndin_sums[i]
  hapto_eighteen[i,] <- hapto_eighteen[i,]/hapto_sums[i]
  chloro_eighteen[i,] <- chloro_eighteen[i,]/chloro_sums[i]
  metazoa_eighteen[i,] <- metazoa_eighteen[i,]/metazoa_sums[i]

}

eight_auto <- as.data.frame(eight_auto)
eight_hetero <- as.data.frame(eight_hetero)
diatom_eighteen <- as.data.frame(diatom_eighteen)
dino_eighteen <- as.data.frame(dino_eighteen)
syndin_eighteen <- as.data.frame(syndin_eighteen)
hapto_eighteen <- as.data.frame(hapto_eighteen)
chloro_eighteen <- as.data.frame(chloro_eighteen)
metazoa_eighteen <- as.data.frame(metazoa_eighteen)

# save data with "_S" to denote new data

# autotrophic 18sv9

scaled_inputs <- eight_auto
scaled_inputs <- as.matrix(scaled_inputs)
asv_table <- auto_copy

save(scaled_inputs, eight_auto, eight_tax_id, asv_table, file = "data/18s_autotrophic_euks_S.Rdata")

# heterotrophic 18sv9

scaled_inputs <- eight_hetero
scaled_inputs <- as.matrix(scaled_inputs)
asv_table <- hetero_copy

save(scaled_inputs, eight_hetero, eight_tax_id, asv_table, file = "data/18s_heterotrophic_euks_S.Rdata")

# Diatom 18sv9

scaled_inputs <- diatom_eighteen
scaled_inputs <- as.matrix(scaled_inputs)
asv_table <- diatom_copy

save(scaled_inputs, diatom_eighteen, eight_tax_id, asv_table, file = "data/18s_diatom_S.Rdata")

# Dinoflagellates 18sv9

scaled_inputs <- dino_eighteen
scaled_inputs <- as.matrix(scaled_inputs)
asv_table <- dino_copy

save(scaled_inputs, dino_eighteen, eight_tax_id, asv_table, file = "data/18s_dino_S.Rdata")

# Syndiniales 18sv9

scaled_inputs <- syndin_eighteen
scaled_inputs <- as.matrix(scaled_inputs)
asv_table <- syndin_copy

save(scaled_inputs, syndin_eighteen, eight_tax_id, asv_table, file = "data/18s_syndin_S.Rdata")

# Haptophyte 18sv9

scaled_inputs <- hapto_eighteen
scaled_inputs <- as.matrix(scaled_inputs)
asv_table <- hapto_copy

save(scaled_inputs, hapto_eighteen, eight_tax_id, asv_table, file = "data/18s_hapto_S.Rdata")

# Metazoa 18sv9

scaled_inputs <- metazoa_eighteen
scaled_inputs <- as.matrix(scaled_inputs)
asv_table <- metazoa_copy

save(scaled_inputs, metazoa_eighteen, eight_tax_id, asv_table, file = "data/18s_metazoa_S.Rdata")

# Chlorophytes 18sv9

scaled_inputs <- chloro_eighteen
scaled_inputs <- as.matrix(scaled_inputs)
asv_table <- chloro_copy

save(scaled_inputs, chloro_eighteen, eight_tax_id, asv_table, file = "data/18s_chloro_S.Rdata")

###### Generate Total Count Data Frames ######

six_rare <- six_rare[match(rownames(eight_rare), rownames(six_rare)),]

totals <- bind_cols(eight_rare,six_rare)

t_rows <- rownames(totals)

totals_sums <- rowSums(totals, na.rm = TRUE)
total_copy <- totals

totals <- as.matrix(totals)

for (i in 1:nrow(totals)){
  totals[i,] <- totals[i,]/totals_sums[i]
}

totals <- as.data.frame(totals)

scaled_inputs <- totals
scaled_inputs <- as.matrix(scaled_inputs)

rownames(total_copy) <- t_rows
asv_table <- total_copy

save(scaled_inputs, totals, asv_table, six_tax_id, eight_tax_id, file = "data/totals_S.Rdata")
save(eight_rare, eight_tax_id, file = "data/all_18Sv9_rare.Rdata")
save(six_rare, six_tax_id, file = "data/all_16S_rare.Rdata")
