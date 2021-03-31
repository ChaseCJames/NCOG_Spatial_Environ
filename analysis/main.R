library(SOMbrero)
library(tidyverse)
library(oce)
library(rworldmap)
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
library(purrr)
library(reshape2)
library(gganimate)
library(gifski)
library(chron)
library(magick)
library(metR)
library(ncdf4)


#### Run data processing file

source("analysis/data_processing.R")


#### Analysis

source("analysis/analysis.R")

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "bacteria_m_euks_16s", "plastid_16s", "cyano_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9", "euks_hetero_18sv9", "euks_auto_18sv9")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales", "Rhodobacterales",
                   "Sar Clade", "Archaea", "Bacteria", "Eukaryotic Phytoplankton (Plastids)",
                   "Cyanobacteria", "Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans", "Eukaryotic Protists", "Photosynthetic Eukaryotic Protists")

in_group_list_basic = c("16s_pro", "16s_syne","16s_flavo", "16s_rhodo", "16s_sar", "16s_archaea",
                        "16s_bacteria_m_euks", "16s_plastids", "16s_cyanos",
                        "18s_diatom","18s_dino", "18s_syndin",
                        "18s_hapto", "18s_chloro", "18s_metazoa", "18s_heterotrophic_euks", "18s_autotrophic_euks")

seedlist <- c(202,5,345,98,133,47,987,13,654,11,134,57,100,133,47,143,916)

for (i in 10:length(in_group_list)) {
  set.seed(seedlist[i])
  
  result_tables(input = paste0("data/", in_group_list_basic[i], ".Rdata"),
                som_file = paste0("output/", in_group_list[i], "_som.Rdata"),
                SST_data = "output/CALCOFI_temp_tables.Rdata",
                physical_data = "data/NCOG_sample_log_DNA_meta_2014-2019.csv",
                map_file = paste0("output/", in_group_list[i], "_map.Rdata"),
                regression_file = paste0("output/", in_group_list[i],"_glm.Rdata"),
                dissimmilar_matrix = paste0("output/", in_group_list[i], "_dissimilar.Rdata"),
                full_data_file = paste0("output/", in_group_list[i], "_full_data.Rdata"),
                sample_regime = "both")
  
  print(i)
  
}

set.seed(214)

result_tables(input = paste0("data/totals.Rdata"),
              som_file = paste0("output/", "total", "_som.Rdata"),
              SST_data = "output/CALCOFI_temp_tables.Rdata",
              physical_data = "data/NCOG_sample_log_DNA_meta_2014-2019.csv",
              map_file = paste0("output/", "total", "_map.Rdata"),
              regression_file = paste0("output/", "total","_glm.Rdata"),
              dissimmilar_matrix = paste0("output/", "total", "_dissimilar.Rdata"), 
              full_data_file = paste0("output/", "total", "_full_data.Rdata"),
              sample_regime = "both")

# generate all figures

# generate figures for paper
