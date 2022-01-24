source("analysis/figure_functions_S.R")

###### Outputs #####

#### SOM Maps ####

# All

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9", "bacteria_m_euks_16s",
                  "cyano_16s", "euks_hetero_18sv9", "euks_auto_18sv9")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                   "Sar 11 Clade", "Archaea","Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans", "Bacteria", "Cyanobacteria",
                   "Eukaryotic Protists", "Photosynthetic Eukaryotic Protists")


plot_list <- list()

for (i in 1:length(in_group_list)) {
  print(i)
  
  plot_list[[i]] <- som_figure(map_file = paste0("output/",in_group_list[i],"_map_S.Rdata"),
                               figure_name = paste0("figures/som_maps/", in_group_list[i],"_map_plot_S.pdf"),
                               main = in_group_names[i], cluster1 = "Nearshore", cluster2 = "Offshore", psize = 6)

  
}

#### Regression Plots #####

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9", "bacteria_m_euks_16s",
                  "cyano_16s", "euks_hetero_18sv9", "euks_auto_18sv9")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                   "Sar 11 Clade", "Archaea","Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans", "Bacteria",
                  "Cyanobacteria",
                   "Eukaryotic Protists", "Photosynthetic Eukaryotic Protists")

var_list = c("temp_mean", "sal_mean", "PO4_mean", "NO3_mean", "SiO3_mean", "NC_mean",
             "temp_coeff", "sal_coeff", "PO4_coeff", "NO3_coeff", "SiO3_coeff", "NC_coeff", "Dist_mean")

var_name_list = c("Mean Temperature (°C)", "Mean Salinity", "Mean PO4ug", "Mean NO3ug",
                  "Mean SiO3ug", "Mean Nitracline Depth (m)",
                  "Coeff. Var. Temperature", "Coeff. Var. Salinity",
                  "Coeff. Var. PO4", "Coeff. Var. NO3", "Coeff. Var. SiO3",
                  "Coeff. Var. Nitracline Depth (m)", "Distance to Coast (km)")

for (i in 1:length(in_group_list)) {
  for (j in 1:length(var_list)) {
    regression_figure(glm_file = paste0("output/",in_group_list[i], "_glm_S.Rdata"),
                      map_file = paste0("output/",in_group_list[i],"_map_S.Rdata"),   
                      figure_name = paste0("figures/glm_plots/", in_group_list[i],"_"),
                      main = in_group_names[i], cluster1 = "Nearshore", cluster2 = "Offshore",
                      var = var_list[j], var_name = var_name_list[j])
  }
}


#### Diversity Figures ####

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9", "bacteria_m_euks_16s",
                  "cyano_16s", "euks_hetero_18sv9", "euks_auto_18sv9")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                   "Sar Clade", "Archaea","Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans", "Bacteria",
                   "Cyanobacteria",
                   "Eukaryotic Protists", "Photosynthetic Eukaryotic Protists")

in_group_list_basic = c("16s_pro", "16s_syne","16s_flavo", "16s_rhodo", "16s_sar", "16s_archaea",
                        "18s_diatom","18s_dino", "18s_syndin",
                        "18s_hapto", "18s_chloro", "18s_metazoa", "16s_bacteria_m_euks",
                        "16s_cyanos", "18s_heterotrophic_euks", "18s_autotrophic_euks")


for (i in 1:length(in_group_list)) {
  
  diveristy_figure(map_file = paste0("output/", in_group_list[i], "_map_S.Rdata"),
                   full_dat = paste0("output/", in_group_list[i], "_full_data_S.Rdata"),
                   figure_start = paste0("figures/diversity/", in_group_list[i], "_"),
                   main = in_group_names[i])
  
  alpha_versus_gamma_figure(full_data_file = paste0("output/", in_group_list[i], "_full_data_S.Rdata"),
                            raw_data_file = paste0("data/", in_group_list_basic[i], "_S.Rdata"),
                            map_file = paste0("output/", in_group_list[i], "_map_S.Rdata"), minimum_tp = 8,
                            figure_name = paste0("figures/diversity/", in_group_list[i], "_alpha_gamma_S.pdf"),
                            main = in_group_names[i])
  
  # beta_diversity_figure(full_data_file = paste0("output/", in_group_list[i], "_full_data.Rdata"),
                        # bc_data_file = paste0("output/", in_group_list[i], "_dissimilar.Rdata"),
                        # raw_data_file = paste0("data/", in_group_list_basic[i], ".Rdata"),
                        # map_file = paste0("output/", in_group_list[i], "_map.Rdata"), minimum_tp = 8,
                        # figure_name = paste0("figures/diversity/", in_group_list[i],"_beta.pdf"),
                        # main = in_group_names[i])
  
  print(i)
  
}





# Diversity

# full_aic_table_figure_diversity(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
#                                                   "diatom_18sv9","dino_18sv9", "syndin_18sv9",
#                                                   "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
#                                 in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar Clade", 
#                                                    "Diatoms",
#                                                    "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"),
#                                 figure_name_2 = paste0("figures/aic_figures/small_group_aic_plot_logit_even",".pdf"),
#                                 title_name = "Variable Importance Evenness", # col 2 = even, col 3 = shan col 4 = rich
#                                 col = 2, color_fill = "purple")
# 
# full_aic_table_figure_diversity(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
#                                                   "diatom_18sv9","dino_18sv9", "syndin_18sv9",
#                                                   "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
#                                 in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar Clade", 
#                                                    "Diatoms",
#                                                    "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"),
#                                 figure_name_2 = paste0("figures/aic_figures/small_group_aic_plot_logit_shannon",".pdf"),
#                                 title_name = "Variable Importance Shannon Diversity", # col 2 = even, col 3 = shan col 4 = rich
#                                 col = 3, color_fill = "red")
# 
# full_aic_table_figure_diversity(in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", 
#                                                   "diatom_18sv9","dino_18sv9", "syndin_18sv9",
#                                                   "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9"),
#                                 in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales", "Sar Clade", 
#                                                    "Diatoms",
#                                                    "Dinoflagellates", "Syndiniales", "Haptophytes", "Chlorophytes","Metazoans"),
#                                 figure_name_2 = paste0("figures/aic_figures/small_group_aic_plot_logit_rich",".pdf"),
#                                 title_name = "Variable Importance Richness", # col 2 = even, col 3 = shan col 4 = rich
#                                 col = 4, color_fill = "blue")

# Basic Groups

# full_aic_table_figure_diversity(in_group_list = c("archaea_16s","bacteria_m_euks_16s", "cyano_16s","plastid_16s",
#                                                   "euks_hetero_18sv9"),
#                                 in_group_names = c("Archaea","Bacteria", "Cyanobacteria", "Eukaryotic Phytoplankton\n(Plastids)",
#                                                    "Eukaryotic\n Protists"),
#                                 figure_name_2 = paste0("figures/aic_figures/big_group_aic_plot_logit_even",".pdf"),
#                                 title_name = "Variable Importance Evenness", # col 2 = even, col 3 = shan col 4 = rich
#                                 col = 2, color_fill = "purple", width_plot = 8)
# 
# full_aic_table_figure_diversity(in_group_list = c("archaea_16s","bacteria_m_euks_16s", "cyano_16s","plastid_16s",
#                                                   "euks_hetero_18sv9"),
#                                 in_group_names = c("Archaea","Bacteria", "Cyanobacteria", "Eukaryotic Phytoplankton\n(Plastids)",
#                                                    "Eukaryotic\n Protists"),
#                                 figure_name_2 = paste0("figures/aic_figures/big_group_aic_plot_logit_shannon",".pdf"),
#                                 title_name = "Variable Importance\nShannon Diversity", # col 2 = even, col 3 = shan col 4 = rich
#                                 col = 3, color_fill = "red", width_plot = 8)
# 
# full_aic_table_figure_diversity(in_group_list = c("archaea_16s","bacteria_m_euks_16s", "cyano_16s","plastid_16s",
#                                                   "euks_hetero_18sv9"),
#                                 in_group_names = c("Archaea","Bacteria", "Cyanobacteria", "Eukaryotic Phytoplankton\n(Plastids)",
#                                                    "Eukaryotic\n Protists"),
#                                 figure_name_2 = paste0("figures/aic_figures/big_group_aic_plot_logit_rich",".pdf"),
#                                 title_name = "Variable Importance Richness", # col 2 = even, col 3 = shan col 4 = rich
#                                 col = 4, color_fill = "blue", width_plot = 8)

##### Community Difference Plots #####

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9", "bacteria_m_euks_16s",
                   "cyano_16s", "euks_hetero_18sv9","euks_auto_18sv9", "total")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                   "SAR 11 Clade", "Archaea","Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans", "Bacteria",
                    "Cyanobacteria",
                   "Eukaryotic Protists","Photosynthetic Eukaryotes", "All ASVs")

fig_list <- list()

for (i in 1:length(in_group_list)) {
  
  # fig_commun_map_func(in_all = paste0("output/",in_group_list[i],"_dissimilar.Rdata"),
  #                   in_dat = paste0("output/",in_group_list[i],"_full_data.Rdata"),
  #                   in_map = paste0("output/",in_group_list[i],"_map.Rdata"),
  #                   community_diff_fig = paste0("figures/community_diff/",in_group_list[i],"_surf_diff.pdf"),
  #                   tsize = 12, psize = 12, group = in_group_names[i])
  
  fig_list[[i]] <- fig_commun_map_surf_deep_func(in_all = paste0("output/",in_group_list[i],"_dissimilar_S.Rdata"),
                                in_dat = paste0("output/",in_group_list[i],"_full_data_S.Rdata"),
                                in_map = paste0("output/",in_group_list[i],"_map_S.Rdata"),
                                community_diff_fig = paste0("figures/community_diff/",
                                                            in_group_list[i],"_both_diff_S.pdf"),
                                tsize = 15, psize = 6, group = in_group_names[i])
  print(i)
  
}

###### Community Time Plots ########

# All


in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9", "bacteria_m_euks_16s",
                  "cyano_16s", "euks_hetero_18sv9","euks_auto_18sv9")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                   "Sar Clade", "Archaea","Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans", "Bacteria",
                   "Cyanobacteria",
                   "Eukaryotic Protists","Photosynthetic Eukaryotes")

for (i in 1:length(in_group_list)) {
  
  community_comparison(in_file = paste0("output/",in_group_list[i],"_full_data_S.Rdata"),
                       in_map = paste0("output/",in_group_list[i],"_map_S.Rdata"),
                       similar_mat = paste0("output/",in_group_list[i],"_dissimilar_S.Rdata"),
                       out_diff_file = paste0("output/",in_group_list[i],"_diffs_S.Rdata"),
                       title = in_group_names[i],
                       upwelling_index = "output/upwelling_indicies.Rdata",
                       index_plot = paste0("figures/",in_group_list[i],"_index_plot_S.pdf"))
  
  print(i)
  
}

# Line Specific Plots 

for (i in 1:length(in_group_list)) {
  
  community_comparison_line(in_file = paste0("output/",in_group_list[i],"_full_data.Rdata"),
                       in_map = paste0("output/",in_group_list[i],"_map.Rdata"),
                       similar_mat = paste0("output/",in_group_list[i],"_dissimilar.Rdata"),
                       out_diff_file = paste0("output/",in_group_list[i],"_line_diffs.Rdata"),
                       title = in_group_names[i])
  
  print(i)
  
}

# Station Specific

for (i in 1:length(in_group_list)) {
  
  community_comparison_station(in_file = paste0("output/",in_group_list[i],"_full_data.Rdata"),
                            in_map = paste0("output/",in_group_list[i],"_map.Rdata"),
                            similar_mat = paste0("output/",in_group_list[i],"_dissimilar.Rdata"),
                            out_diff_file = paste0("output/",in_group_list[i],"_station_diffs.Rdata"),
                            title = in_group_names[i])
  
  community_comparison_station_2(in_file = paste0("output/",in_group_list[i],"_full_data.Rdata"),
                               in_map = paste0("output/",in_group_list[i],"_map.Rdata"),
                               similar_mat = paste0("output/",in_group_list[i],"_dissimilar.Rdata"),
                               out_diff_file = paste0("output/",in_group_list[i],"_station2_diffs.Rdata"),
                               title = in_group_names[i])
  
  print(i)
  
}

###### Diversity Time Plots #####

# All

in_group_list = c("pro_16s", "syne_16s","flavo_16s", "rhodo_16s", "sar_16s", "archaea_16s",
                  "diatom_18sv9","dino_18sv9", "syndin_18sv9",
                  "hapto_18sv9", "chloro_18sv9", "metazoa_18sv9", "bacteria_m_euks_16s",
                  "cyano_16s", "euks_hetero_18sv9", "euks_auto_18sv9")

in_group_names = c("Prochlorococcus", "Synecococcus", "Flavobacteriales","Rhodobacterales",
                   "Sar Clade", "Archaea","Diatoms", "Dinoflagellates", "Syndiniales",
                   "Haptophytes", "Chlorophytes","Metazoans", "Heterotrophic Bacteria",
                    "Cyanobacteria",
                   "Heterotrophic Eukaryotic Protists", "Photosynthetic Eukaryotic Protists")

in_group_list_basic = c("16s_pro", "16s_syne","16s_flavo", "16s_rhodo", "16s_sar", "16s_archaea",
                        "18s_diatom","18s_dino", "18s_syndin",
                        "18s_hapto", "18s_chloro", "18s_metazoa", "16s_bacteria_m_euks",
                         "16s_cyanos", "18s_heterotrophic_euks", "18s_autotrophic_euks")

for (i in 1:length(in_group_list)) {
  print(i)
  
  diversity_comparison(in_file = paste0("output/",in_group_list[i],"_full_data_S.Rdata"),
                       in_map = paste0("output/",in_group_list[i],"_map_S.Rdata"),
                       in_raw = paste0("data/", in_group_list_basic[i],"_S.Rdata"),
                       out_diff_file = paste0("output/",in_group_list[i],"_diffs_div_S.Rdata"),
                       title = in_group_names[i],
                       upwelling_index = "output/upwelling_indicies.Rdata",
                       index_plot = paste0("figures/",in_group_list[i],"_index_plot_div_S.pdf"))
  
}
