# Function for running eco som analysis on groups

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

########## Data Processing ############

set.seed(532)

load("output/CALCOFI_temp_tables.Rdata")

physical_dat <- read.csv("data/NCOG_sample_metadata.csv", header = TRUE, stringsAsFactors = FALSE)

# check variables

data <- physical_dat[colnames(physical_dat)[c(1,33:41)]]

# data <- data[which(complete.cases(data)==TRUE),]

data$spice <- swSpice(data$Salnty, data$T_degC)

data$eco_name <- paste0("X20",substr(data$Sample.Name,start = 3,stop = 100))

full_dat <- physical_dat[!is.na(match(physical_dat$Sample.Name, data$Sample.Name)),]

full_dat$eco_name <- paste0("X20",substr(full_dat$Sample.Name,start = 3,stop = 100))

# read in 18s
eighteen_s <- read.csv("data/asv_count_tax_final_update.csv", stringsAsFactors = FALSE)

eight_id_names <- eighteen_s$Feature.ID

eight_tax_id <- eighteen_s[,c(1,ncol(eighteen_s)-1, ncol(eighteen_s))]

eighteen_s <- eighteen_s[,-c(1,ncol(eighteen_s)-1, ncol(eighteen_s))]

eighteen_s <- apply(eighteen_s, 2, as.numeric)

eighteen_s <- t(eighteen_s)

eighteen_s <- as.data.frame(eighteen_s)

eight_tp <- rownames(eighteen_s)

colnames(eighteen_s) <- eight_id_names

# remove mocks from data and rows with no reads
eighteen_s <- eighteen_s[-c(479:491,669:680,873:878),]

r_sums <- rowSums(eighteen_s, na.rm = TRUE)

eighteen_s <- eighteen_s[-which(r_sums == 0),]

r_sums <- r_sums[-which(r_sums == 0)]

eighteen_s <- as.matrix(eighteen_s)

# physical data match to biological data

full_dat <- full_dat[which(!is.na(match(full_dat$eco_name, rownames(eighteen_s)))),]

data <- data[which(!is.na(match(data$eco_name, rownames(eighteen_s)))),]

# split by taxa

taxas <- eight_tax_id[which(!is.na(match(eight_tax_id$Feature.ID, colnames(eighteen_s)))),]

split_taxa <- separate(taxas, Taxon, sep = ";", into = c("A","B","C", "D", "E", "F", "G", "H", "I"))

taxa_id_list <- list()
taxa_groups <- c("Diatoms", "Dinoflagellates", "Haptophytes","Crustacea",
                 "Radiolaria", "Cnidaria", "Annelida")

taxa_id_list[[1]] <- split_taxa$Feature.ID[which(split_taxa$D == "Bacillariophyta")]
taxa_id_list[[2]] <- split_taxa$Feature.ID[which(split_taxa$C == "Dinoflagellata")]
taxa_id_list[[3]] <- split_taxa$Feature.ID[which(split_taxa$C == "Haptophyta")]
taxa_id_list[[4]] <- split_taxa$Feature.ID[which(split_taxa$E == "Crustacea")]
taxa_id_list[[5]] <- split_taxa$Feature.ID[which(split_taxa$C == "Radiolaria")]
taxa_id_list[[6]] <- split_taxa$Feature.ID[which(split_taxa$D == "Cnidaria")]
taxa_id_list[[7]] <- split_taxa$Feature.ID[which(split_taxa$D == "Annelida")]

group_analysis <- function(in_amplicon = eighteen_s, group = taxa_id_list[[1]], group_name = taxa_groups[[1]],
                           environ_dat = full_dat, coeff_table = coeff_table, mean_table = mean_table){
  
  
  # Identify subset of amplicon data
  
  subset <- in_amplicon[,which(!is.na(match(colnames(in_amplicon), group)))]
  
  r_sums <- rowSums(subset, na.rm = TRUE)
  
  subset <- subset[-which(r_sums == 0),]
  in_amplicon <- in_amplicon[-which(r_sums == 0),]
  r_sums <- r_sums[-which(r_sums == 0)]
  
  for (j in 1:nrow(subset)){
    
    subset[j,] <- subset[j,]/r_sums[j]
    
  }
  
  subset <- as.data.frame(subset)
  
  # run som
  out_som <- trainSOM(x.data = subset, dimension = c(5, 5), nb.save = 10, maxit = 2000, 
                      scaling = "none")
  
  # cluster som
  out_clust <- superClass(out_som, k = 2)
  
  clusters <- out_clust$cluster
  
  ids <- out_clust$som$clustering
  
  som_ids <- clusters[ids]
  
  in_amplicon <- as.data.frame(in_amplicon)
  
  in_amplicon$som_id <- som_ids
  
  # label environmental data
  environ_dat$som_id <- in_amplicon$som_id[match(full_dat$eco_name, rownames(in_amplicon))]
  
  # diversity evenness
  
  in_amplicon$shannon_index <- diversity(subset, MARGIN = 1, index = "shannon")
  environ_dat$shannon <- in_amplicon$shannon_index[match(full_dat$eco_name, rownames(in_amplicon))]
  
  S <- apply(subset>0,1,sum)
  in_amplicon$evenness <- diversity(subset, index="shannon")/log(S)
  environ_dat$evenness <- in_amplicon$evenness[match(full_dat$eco_name, rownames(in_amplicon))]
  
  # look at stations
  som_maps <- environ_dat %>% 
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som_id == 1)/n(), som_2 = sum(som_id == 2)/n(),
              n_samps = n(), lat = mean(Lat_Dec, na.rm = TRUE), long = mean(Lon_Dec, na.rm = TRUE),
              PO4_mean = mean(PO4ug, na.rm = TRUE), NO3_mean = mean(NO3ug, na.rm = TRUE),
              SiO3_mean = mean(SiO3ug, na.rm = TRUE),
              PO4_coeff =  sd(PO4ug, na.rm = TRUE)/mean(PO4ug, na.rm = TRUE),
              NO3_coeff = sd(NO3ug, na.rm = TRUE)/mean(NO3ug, na.rm = TRUE),
              SiO3_coeff = sd(SiO3ug, na.rm = TRUE)/mean(SiO3ug, na.rm = TRUE),
              shannon = mean(shannon, na.rm = TRUE), evenness = mean(evenness, na.rm = TRUE)
    )
  
  som_plot_list <- list()
  
  names_vect <- c("Cluster 1", "Cluster 2")
  color_vect <- c("darkblue", "darkred")
  
  # FOR NOW REMOVE NORTHERN TRANSECTS
  
  som_maps <- som_maps[-c(1:12),]
  
  map <- map_data("world")    
  
  for (i in 1:2) {
    
    p <-  ggplot() + 
      geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
      coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
      xlab("Longitude") + ylab("Latitude") + 
      geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = paste0("som_",i)), color = "black", size =5, stroke = 0.1, shape = 21) +
      scale_fill_gradient(low = "white", high = color_vect[i], limits = c(0,1)) +
      ggtitle(paste0("Rel. % ",group_name," ",names_vect[i]," per Station")) +
      theme(legend.title = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA),
            title = element_text(hjust = 0.5))
    
    print(p)
    
    som_plot_list[[i]] <- p
    
  }
  
  # shannon
  shannon <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "shannon"), color = "black", size =5, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "blue", aes(min = min(shannon, na.rm = TRUE, max = max(shannon, na.rm = TRUE)), limits = c(min,max))) +
    ggtitle(paste0("Shannon Diversity")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          title = element_text(hjust = 0.5))
  
  evenness <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "evenness"), color = "black", size =5, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "red", aes(min = min(evenness, na.rm = TRUE, max = max(shannon, na.rm = TRUE)), limits = c(min,max))) +
    ggtitle(paste0("Evenness")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          title = element_text(hjust = 0.5))
  
  # Run random forests
  
  environ_dat$Date <- as.Date(environ_dat$Date, format = "%m/%d/%Y")
  environ_dat$Julian <- as.numeric(format(environ_dat$Date, "%j"))
  
  forest_dat <- environ_dat[,-which(is.na(rowSums(environ_dat[,c(27,28,33,34,37,38,39,59,60)])))]
  
  forest_dat$som_id <- as.factor(forest_dat$som_id)
  
  train <- sample(nrow(forest_dat), 0.7*nrow(forest_dat), replace = FALSE)
  
  forest_out <- randomForest(som_id ~ Lat_Dec + Lon_Dec + CC_Depth + T_degC + Salnty +
                                 PO4ug + SiO3ug + NO3ug + Julian, data = forest_dat, importance = TRUE,
                               na.action = na.omit, subset = train, mtry = 4)
  
  
  # Predicting on validation set
  predTrain <- predict(forest_out, forest_dat[train,], type = "class")
  # Checking classification accuracy
  table(predTrain, forest_dat$som_id[train]) 
  
  # Predicting on validation set
  predValid <- predict(forest_out, forest_dat[-train,], type = "class")
  # Checking classification accuracy
  table(predValid, forest_dat$som_id[-train])  

  varImpPlot(forest_out, type = 1, main = paste0("Random Forest"))
  rf_plot <- recordPlot()
  
  # GLM
  
  som_maps <- som_maps[which(som_maps$n_samps > 7),]
  
  # Identify closest grid to station 
  
  coeff_val <- vector()
  mean_val <- vector()
  
  for (i in 1:nrow(som_maps)) {
    
    lat_1 <- som_maps$lat[i]
    lon_1 <- som_maps$long[i]
    
    distance_vect <- sqrt((lat_1 - coeff_table$lat)^2 + (lon_1 - coeff_table$lon)^2)
    
    coeff_val[i] <- coeff_table$coeff_var[which.min(distance_vect)]
    mean_val[i] <- mean_table$Mean[which.min(distance_vect)]
    
    
  }
  
  som_maps$temp_mean <- mean_val
  som_maps$temp_coeff <- coeff_val
  
  # glm
  
  summary_fit <- summary(lm(som_maps$som_1~som_maps$temp_mean))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_mean <- summary_fit$r.squared
  
  summary_fit <- summary(lm(som_maps$som_1~som_maps$temp_coeff))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_coeff <- summary_fit$r.squared
  
  # mean coeff temp plots
  
  som_plots <- gather(som_maps, cluster, freq, 2:3)
  
  mean_plot <- ggplot(som_plots, aes(x = temp_mean, y = freq, color = cluster)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid")) +
    stat_smooth(method = "lm", level = 0.95) +
    xlab("Mean SST (°C)") + ylab("Frequency") +
    annotate("text", x = 15, y = 0.08, label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 2)
  
  coeff_plot <- ggplot(som_plots, aes(x = temp_coeff, y = freq, color = cluster)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none") +
    stat_smooth(method = "lm", level = 0.95) +
    xlab("SST Coeff. Var") + ylab("Frequency") +
    annotate("text", x = 0.095, y = 0.3, label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 2)
  
  theme_set(theme_cowplot(font_size=8))
  temp_plots <- plot_grid(coeff_plot, mean_plot, ncol = 2,rel_widths = c(0.8,1))
  
  # Diversity
  
  summary_fit <- summary(lm(som_maps$shannon~som_maps$temp_mean))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_mean <- summary_fit$r.squared
  
  summary_fit <- summary(lm(som_maps$shannon~som_maps$temp_coeff))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_coeff <- summary_fit$r.squared
  
  mean_plot <- ggplot(som_plots, aes(x = temp_mean, y = shannon)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid")) +
    stat_smooth(method = "lm", level = 0.95) +
    xlab("Mean SST (°C)") + ylab("Shannon Index") +
    annotate("text", x = 15, y = 1.25, label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 2)
  
  coeff_plot <- ggplot(som_plots, aes(x = temp_coeff, y = shannon)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none") +
    stat_smooth(method = "lm", level = 0.95) +
    xlab("SST Coeff. Var") + ylab("Shannon Index") +
    annotate("text", x = 0.13, y = 1.5, label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 2)
  
  theme_set(theme_cowplot(font_size=8))
  temp_diversity_plots <- plot_grid(coeff_plot, mean_plot, ncol = 1,rel_widths = c(1,1))
  
  # Evenness
  
  summary_fit <- summary(lm(som_maps$evenness~som_maps$temp_mean))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_mean <- summary_fit$r.squared
  
  summary_fit <- summary(lm(som_maps$evenness~som_maps$temp_coeff))
  
  p_val <- summary_fit$coefficients[2,4]
  r_sq_coeff <- summary_fit$r.squared
  
  mean_plot <- ggplot(som_plots, aes(x = temp_mean, y = evenness)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid")) +
    stat_smooth(method = "lm", level = 0.95) +
    xlab("Mean SST (°C)") + ylab("Evenness") +
    annotate("text", x = 15, y = 0.85, label = paste0("R Squared = ", round(r_sq_mean, digits = 3)), size = 2)
  
  coeff_plot <- ggplot(som_plots, aes(x = temp_coeff, y = evenness)) +
    geom_point(pch = 20, size = 4) +
    theme(legend.title = element_blank(),
          panel.border = element_rect(color = "black", linetype = "solid"),
          axis.line = element_blank(),
          legend.background = element_rect(color = "black", linetype = "solid"),
          legend.position = "none") +
    stat_smooth(method = "lm", level = 0.95) +
    xlab("SST Coeff. Var") + ylab("Evenness") +
    annotate("text", x = 0.1, y = 0.75, label = paste0("R Squared = ", round(r_sq_coeff, digits = 3)), size = 2)
  
  theme_set(theme_cowplot(font_size=8))
  temp_evenness_plots <- plot_grid(coeff_plot, mean_plot, ncol = 1,rel_widths = c(1,1))
  
  # subset
  
  response_df <- som_maps$som_1[which(!is.na(som_maps$som_1))]
  predictor_df <- som_maps[which(!is.na(som_maps$som_1)),5:14]
  
  regsubsetsObj <- regsubsets(x=predictor_df ,y=response_df, nbest = 2, really.big = T)
  regsubsetsObj$xnames <- c("(Intercept)", "Latitude", "Longitude", "Mean Temp", "Mean PO4",
                            "Mean NO3", "Mean SiO3", "Coeff. Var. Temp", "Coeff. Var. PO4",
                            "Coeff. Var. NO3", "Coeff. Var. SiO3")
  plot(regsubsetsObj, scale = "adjr2", main = "Backwards Subset Selection")
  sub_plot <- recordPlot()
  

  
  
  
  
  
  return(list(out_som = out_som, out_clust = out_clust,
              som_plot_list = som_plot_list, 
              rf_out = forest_out, rf_importance = rf_plot,
              temp_plots = temp_plots, glm_subset = sub_plot,
              shannon = shannon, evenness = evenness,
              temp_diversity_plots = temp_diversity_plots,
              temp_evenness_plots = temp_evenness_plots))
  
  }

# Diatoms

diatom_results <- group_analysis(in_amplicon = eighteen_s, group = taxa_id_list[[1]], group_name = taxa_groups[[1]],
                                 environ_dat = full_dat, coeff_table = coeff_table, mean_table = mean_table)

title <- ggdraw() + draw_label("Diatoms", fontface='bold')

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/diatom_results.pdf", width = 14, height = 15)
plot_grid(title, plot_grid(diatom_results$som_plot_list[[1]],diatom_results$som_plot_list[[2]]), 
          plot_grid(diatom_results$rf_importance, diatom_results$glm_subset, scale = c(0.7,0.7)),
          diatom_results$temp_plots,
          nrow = 4, ncol = 1, rel_heights = c(0.1,1,1.2,1), labels = c("", "A","B","C"))
dev.off()

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/diatom_diversity_evenness.pdf", width = 10, height = 10)
plot_grid(title, plot_grid(diatom_results$shannon,diatom_results$evenness, labels = c("A","B")), 
          plot_grid(diatom_results$temp_diversity_plots,
          diatom_results$temp_evenness_plots, ncol =  2,
          labels = c("C","D")),
          nrow = 3, ncol = 1, rel_heights = c(0.1,1,1,1))
dev.off()

#Dinoflagellates

dino_results <- group_analysis(in_amplicon = eighteen_s, group = taxa_id_list[[2]], group_name = taxa_groups[[2]],
                                 environ_dat = full_dat, coeff_table = coeff_table, mean_table = mean_table)

title <- ggdraw() + draw_label("Dinoflagellates", fontface='bold')

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/dinoflagellate_results.pdf", width = 14, height = 15)
plot_grid(title, plot_grid(dino_results$som_plot_list[[1]],dino_results$som_plot_list[[2]]), 
          plot_grid(dino_results$rf_importance, dino_results$glm_subset, scale = c(0.7,0.7)),
          dino_results$temp_plots,
          nrow = 4, ncol = 1, rel_heights = c(0.1,1,1.2,1), labels = c("", "A","B","C"))
dev.off()

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/dinoflagellate_diversity_evenness.pdf", width = 10, height = 10)
plot_grid(title, plot_grid(dino_results$shannon,dino_results$evenness, labels = c("A","B")), 
          plot_grid(dino_results$temp_diversity_plots,
                    dino_results$temp_evenness_plots, ncol =  2,
                    labels = c("C","D")),
          nrow = 3, ncol = 1, rel_heights = c(0.1,1,1,1))
dev.off()

#Haptophytes

hapto_results <- group_analysis(in_amplicon = eighteen_s, group = taxa_id_list[[3]], group_name = taxa_groups[[3]],
                               environ_dat = full_dat, coeff_table = coeff_table, mean_table = mean_table)

title <- ggdraw() + draw_label("Haptophytes", fontface='bold')

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/haptophyte_results.pdf", width = 14, height = 15)
plot_grid(title, plot_grid(hapto_results$som_plot_list[[1]],hapto_results$som_plot_list[[2]]), 
          plot_grid(hapto_results$rf_importance, hapto_results$glm_subset, scale = c(0.7,0.7)),
          hapto_results$temp_plots,
          nrow = 4, ncol = 1, rel_heights = c(0.1,1,1.2,1), labels = c("", "A","B","C"))
dev.off()

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/haptophyte_diversity_evenness.pdf", width = 10, height = 10)
plot_grid(title, plot_grid(hapto_results$shannon,hapto_results$evenness, labels = c("A","B")), 
          plot_grid(hapto_results$temp_diversity_plots,
                    hapto_results$temp_evenness_plots, ncol =  2,
                    labels = c("C","D")),
          nrow = 3, ncol = 1, rel_heights = c(0.1,1,1,1))
dev.off()

#Crustacea

crust_results <- group_analysis(in_amplicon = eighteen_s, group = taxa_id_list[[4]], group_name = taxa_groups[[4]],
                                environ_dat = full_dat, coeff_table = coeff_table, mean_table = mean_table)

title <- ggdraw() + draw_label("Crustacea", fontface='bold')

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/crustacea_results.pdf", width = 14, height = 15)
plot_grid(title, plot_grid(crust_results$som_plot_list[[1]],crust_results$som_plot_list[[2]]), 
          plot_grid(crust_results$rf_importance, crust_results$glm_subset, scale = c(0.7,0.7)),
          crust_results$temp_plots,
          nrow = 4, ncol = 1, rel_heights = c(0.1,1,1.2,1), labels = c("", "A","B","C"))
dev.off()

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/crustacea_diversity_evenness.pdf", width = 10, height = 10)
plot_grid(title, plot_grid(crust_results$shannon,crust_results$evenness, labels = c("A","B")), 
          plot_grid(crust_results$temp_diversity_plots,
                    crust_results$temp_evenness_plots, ncol =  2,
                    labels = c("C","D")),
          nrow = 3, ncol = 1, rel_heights = c(0.1,1,1,1))
dev.off()

#Radiolaria

radio_results <- group_analysis(in_amplicon = eighteen_s, group = taxa_id_list[[5]], group_name = taxa_groups[[5]],
                                environ_dat = full_dat, coeff_table = coeff_table, mean_table = mean_table)

title <- ggdraw() + draw_label("Radiolaria", fontface='bold')

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/radiolaria_results.pdf", width = 14, height = 15)
plot_grid(title, plot_grid(radio_results$som_plot_list[[1]],radio_results$som_plot_list[[2]]), 
          plot_grid(radio_results$rf_importance, radio_results$glm_subset, scale = c(0.7,0.7)),
          radio_results$temp_plots,
          nrow = 4, ncol = 1, rel_heights = c(0.1,1,1.2,1), labels = c("", "A","B","C"))
dev.off()

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/radiolaria_diversity_evenness.pdf", width = 10, height = 10)
plot_grid(title, plot_grid(radio_results$shannon,radio_results$evenness, labels = c("A","B")), 
          plot_grid(radio_results$temp_diversity_plots,
                    radio_results$temp_evenness_plots, ncol =  2,
                    labels = c("C","D")),
          nrow = 3, ncol = 1, rel_heights = c(0.1,1,1,1))
dev.off()

#Cnidaria

cnidaria_results <- group_analysis(in_amplicon = eighteen_s, group = taxa_id_list[[6]], group_name = taxa_groups[[6]],
                                environ_dat = full_dat, coeff_table = coeff_table, mean_table = mean_table)

title <- ggdraw() + draw_label("Cnidaria", fontface='bold')

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/cnidaria_results.pdf", width = 14, height = 15)
plot_grid(title, plot_grid(cnidaria_results$som_plot_list[[1]],cnidaria_results$som_plot_list[[2]]), 
          plot_grid(cnidaria_results$rf_importance, cnidaria_results$glm_subset, scale = c(0.7,0.7)),
          cnidaria_results$temp_plots,
          nrow = 4, ncol = 1, rel_heights = c(0.1,1,1.2,1), labels = c("", "A","B","C"))
dev.off()

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/cnidaria_diversity_evenness.pdf", width = 10, height = 10)
plot_grid(title, plot_grid(cnidaria_results$shannon,cnidaria_results$evenness, labels = c("A","B")), 
          plot_grid(cnidaria_results$temp_diversity_plots,
                    cnidaria_results$temp_evenness_plots, ncol =  2,
                    labels = c("C","D")),
          nrow = 3, ncol = 1, rel_heights = c(0.1,1,1,1))
dev.off()


#Annelida

annelida_results <- group_analysis(in_amplicon = eighteen_s, group = taxa_id_list[[6]], group_name = taxa_groups[[6]],
                                   environ_dat = full_dat, coeff_table = coeff_table, mean_table = mean_table)

title <- ggdraw() + draw_label("Annelida", fontface='bold')

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/annelida_results.pdf", width = 14, height = 15)
plot_grid(title, plot_grid(annelida_results$som_plot_list[[1]],annelida_results$som_plot_list[[2]]), 
          plot_grid(annelida_results$rf_importance, annelida_results$glm_subset, scale = c(0.7,0.7)),
          annelida_results$temp_plots,
          nrow = 4, ncol = 1, rel_heights = c(0.1,1,1.2,1), labels = c("", "A","B","C"))
dev.off()

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/annelida_diversity_evenness.pdf", width = 10, height = 10)
plot_grid(title, plot_grid(annelida_results$shannon,annelida_results$evenness, labels = c("A","B")), 
          plot_grid(annelida_results$temp_diversity_plots,
                    annelida_results$temp_evenness_plots, ncol =  2,
                    labels = c("C","D")),
          nrow = 3, ncol = 1, rel_heights = c(0.1,1,1,1))
dev.off()
