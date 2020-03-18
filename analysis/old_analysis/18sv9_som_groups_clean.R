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
library(spatialEco)
library(geosphere)
library(viridis)
library(RColorBrewer)

########## Data Processing ############

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
  out_som <- trainSOM(x.data = subset, dimension = c(5, 5), maxit = 5000, 
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
  
  if(length(which(is.na(environ_dat$som_id)))>0){environ_dat <- environ_dat[-which(is.na(environ_dat$som_id)),]}
  
  # look at stations
  som_maps <- environ_dat %>% 
    group_by(Sta_ID) %>%
    summarise(som_1 = sum(som_id == 1,na.rm = TRUE)/n(), som_2 = sum(som_id == 2,na.rm = TRUE)/n(),
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
  centroid_vect <- c("red","blue")
  
  # FOR NOW REMOVE NORTHERN TRANSECTS
  
  som_maps <- som_maps[-c(1:12),]
  
  map <- map_data("world")    
  
  # find centroids
  clean_som <- som_maps[,c(6,5,2,3,13,14)]
  clean_som <- clean_som[which(!is.na(rowSums(clean_som))),]
  
  centroid_df <- SpatialPointsDataFrame(coords = clean_som[,1:2], data = clean_som)
  
  centroid_list <- list()
  
  centroid_list[[1]] <- wt.centroid(x = centroid_df , p = 3)
  centroid_list[[2]] <- wt.centroid(x = centroid_df , p = 4)
  centroid_list[[3]] <- wt.centroid(x = centroid_df , p = 5)
  centroid_list[[4]] <- wt.centroid(x = centroid_df , p = 6)
  


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
  
  som_plot_list[[1]] <- som_plot_list[[1]] + geom_point(aes(x = centroid_list[[1]]@coords[1], y = centroid_list[[1]]@coords[2]),
                      color = centroid_vect[1], size = 5, pch = 10)
  
  som_plot_list[[2]] <- som_plot_list[[2]] + geom_point(aes(x = centroid_list[[2]]@coords[1], y = centroid_list[[2]]@coords[2]),
                                                        color = centroid_vect[2], size = 5, pch = 10)
  
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

  shannon <- shannon + geom_point(aes(x = centroid_list[[3]]@coords[1], y = centroid_list[[3]]@coords[2]), color = "red", size = 5, pch = 10)
  
  evenness <- ggplot() + 
    geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
    coord_fixed(xlim = c(-125.3206, -116.3055),ylim= c(28.84998,38.08734), 1.3) +
    xlab("Longitude") + ylab("Latitude") + 
    geom_point(data = som_maps, aes_string(x = "long", y = "lat", fill = "evenness"), color = "black", size =5, stroke = 0.1, shape = 21) +
    scale_fill_gradient(low = "white", high = "red", aes(min = min(evenness, na.rm = TRUE, max = max(evenness, na.rm = TRUE)), limits = c(min,max))) +
    ggtitle(paste0("Evenness")) +
    theme(legend.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          title = element_text(hjust = 0.5))

  evenness <- evenness + geom_point(aes(x = centroid_list[[4]]@coords[1], y = centroid_list[[4]]@coords[2]), color = "blue", size = 5, pch = 10)
  
  # distance to coast
  
  coast_calc <- vector()
  
  for (i in 1:nrow(environ_dat)) {
    
    long_lat <-  c(environ_dat$Lon_Dec[i], environ_dat$Lat_Dec[i])
    
    if (length(which(!is.na(long_lat))) == 2) {
      
      map_dist <- map[,1:2]
      
      distances <- distm(long_lat, map_dist, fun = distGeo)
      
      coast_calc[i] <- min(distances, na.rm = TRUE)
      
    }
    
    else{coast_calc[i] <- NA}
    
  }
  
  environ_dat$dist_to_coast <- coast_calc/1000
  
  # Run random forests
  
  environ_dat$Date <- as.Date(environ_dat$Date, format = "%m/%d/%Y")
  environ_dat$Julian <- as.numeric(format(environ_dat$Date, "%j"))
  
  forest_dat <- environ_dat[-which(is.na(rowSums(environ_dat[,c(33,34,37,38,39,59,62,63)]))),]
  
  forest_dat$som_id <- as.factor(forest_dat$som_id)
  
  train <- sample(nrow(forest_dat), 0.7*nrow(forest_dat), replace = FALSE)
  
  forest_out <- randomForest(som_id ~  CC_Depth + T_degC +
                                 PO4ug + SiO3ug + NO3ug, data = forest_dat, importance = TRUE,
                               na.action = na.omit, subset = train, mtry = 4)
  
  
  # Predicting on validation set
  predTrain <- predict(forest_out, forest_dat[train,], type = "class")
  # Checking classification accuracy
  table(predTrain, forest_dat$som_id[train]) 
  
  # Predicting on validation set
  predValid <- predict(forest_out, forest_dat[-train,], type = "class")
  # Checking classification accuracy
  table(predValid, forest_dat$som_id[-train])  

  varImpPlot(forest_out, type = 2, main = paste0("Random Forest"))
  rf_plot <- recordPlot()
  
  # GLM
  
  som_maps <- som_maps[which(som_maps$n_samps > 7),]
  
  # Identify closest grid to station 
  
  coeff_val <- vector()
  mean_val <- vector()
  
  coast_distance <- vector()
  
  for (i in 1:nrow(som_maps)) {
    
    long_lat <-  c(som_maps$long[i], som_maps$lat[i])
    sst_dat <- coeff_table[,1:2]  
    map_dat <- map[,1:2]
    
    sst_dist <- distm(long_lat, sst_dat, fun = distGeo)
    map_dist <- distm(long_lat, map_dat, fun = distGeo)
    
    coeff_val[i] <- coeff_table$coeff_var[which.min(sst_dist)]
    mean_val[i] <- mean_table$Mean[which.min(sst_dist)]
    
    # find distance between coast and location
    coast_distance[i] <- min(map_dist, na.rm = TRUE)
    
    
  }
  
  som_maps$temp_mean <- mean_val
  som_maps$temp_coeff <- coeff_val
  som_maps$dist_to_coast <- coast_distance/1000
  
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
  

# extra plots
  
  plot(environ_dat$dist_to_coast, environ_dat$shannon, xlab = "Distance to Coast (km)",
       ylab = "Diversity (Shannon)")
  
  abline(lm(environ_dat$shannon~environ_dat$dist_to_coast), col = "red", lty = 1, lwd = 2)
  
  plot(environ_dat$dist_to_coast, environ_dat$evenness, xlab = "Distance to Coast (km)",
       ylab = "Evenness")
  
  abline(lm(environ_dat$evenness~environ_dat$dist_to_coast), col = "red", lty = 1, lwd = 2)
  
  shannon_depth_dist <- ggplot(environ_dat, aes(x = dist_to_coast, y = CC_Depth, fill = shannon)) + 
    geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 2, name = "Shannon Diversity") +
    ylab("Depth (m)") + xlab("Distance to Coast (km)") + ggtitle(group_name)
  
  even_depth_dist <- ggplot(environ_dat, aes(x = dist_to_coast, y = CC_Depth, fill = evenness)) + 
    geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred",
                         name = "Evenness", midpoint = 0.7) +
    ylab("Depth (m)") + xlab("Distance to Coast (km)") + ggtitle(group_name)
  
  environ_dat$Date <- as.Date(environ_dat$Date, format = "%m/%d/%Y")
  
  shannon_timeplot <- ggplot(environ_dat, aes(x = Date, y = dist_to_coast, fill = shannon)) + 
    geom_jitter(size = 4, pch = 21, color = "black",height = 0, width = 10) + 
    scale_fill_viridis(name = "Shannon", option = "D") +
    ylab("Distance to Coast") + xlab("Time") 
  
  even_timeplot <- ggplot(environ_dat, aes(x = Date, y = dist_to_coast, fill = evenness)) + 
    geom_jitter(size = 4, pch = 21, color = "black",height = 0, width = 10) + 
    scale_fill_viridis(name = "Evenness", option = "D") +
    ylab("Distance to Coast") + xlab("Time") 

  
  return(list(out_som = out_som, out_clust = out_clust,
              som_plot_list = som_plot_list, 
              rf_out = forest_out, rf_importance = rf_plot,
              temp_plots = temp_plots, glm_subset = sub_plot,
              shannon = shannon, evenness = evenness,
              temp_diversity_plots = temp_diversity_plots,
              temp_evenness_plots = temp_evenness_plots,
              environ_dat = environ_dat, som_maps = som_maps,
              shannon_depth_dist = shannon_depth_dist,
              even_depth_dist = even_depth_dist,
              centroid_list = centroid_list,
              shannon_timeplot = shannon_timeplot,
              even_timeplot = even_timeplot,
              asv_list = subset))
  
  }

# Diatoms

set.seed(15)

diatom_results <- group_analysis(in_amplicon = eighteen_s, group = taxa_id_list[[1]], group_name = taxa_groups[[1]],
                                 environ_dat = full_dat, coeff_table = coeff_table, mean_table = mean_table)

title <- ggdraw() + draw_label("Diatoms", fontface='bold')

theme_set(theme_cowplot(font_size=10))
# pdf(file = "figures/diatom_results.pdf", width = 14, height = 15)
# plot_grid(title, plot_grid(diatom_results$som_plot_list[[1]],diatom_results$som_plot_list[[2]]), 
#           plot_grid(diatom_results$rf_importance, diatom_results$glm_subset, scale = c(0.7,0.7)),
#           diatom_results$temp_plots,
#           nrow = 4, ncol = 1, rel_heights = c(0.1,1,1.2,1), labels = c("", "A","B","C"))
# dev.off()

# additional diversity evenness plots

shannon_depth_dist <- ggplot(diatom_results$environ_dat, aes(x = dist_to_coast, y = CC_Depth, fill = shannon)) + 
  geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
  scale_fill_viridis(name = "Shannon Diversity", option = "D") +
  ylab("Depth (m)") + xlab("Distance to Coast (km)") 

even_depth_dist <- ggplot(diatom_results$environ_dat, aes(x = dist_to_coast, y = CC_Depth, fill = evenness)) + 
  geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
  scale_fill_viridis(name = "Evenness", option = "D") +
  ylab("Depth (m)") + xlab("Distance to Coast (km)") 

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/diatom_diversity_evenness.pdf", width = 10, height = 13)
plot_grid(title, plot_grid(diatom_results$shannon,diatom_results$evenness, labels = c("A","B")),
          plot_grid(diatom_results$temp_diversity_plots,
          diatom_results$temp_evenness_plots, ncol =  2,
          labels = c("C","D")), plot_grid(shannon_depth_dist, even_depth_dist, labels = c("E","F")),
          plot_grid(diatom_results$shannon_timeplot, diatom_results$even_timeplot, ncol = 2, labels = c("G", "H")),
          nrow = 5, ncol = 1, rel_heights = c(0.1,1,1.2,0.8,0.8))
dev.off()

#Dinoflagellates

set.seed(123)

dino_results <- group_analysis(in_amplicon = eighteen_s, group = taxa_id_list[[2]], group_name = taxa_groups[[2]],
                                 environ_dat = full_dat, coeff_table = coeff_table, mean_table = mean_table)

title <- ggdraw() + draw_label("Dinoflagellates", fontface='bold')

theme_set(theme_cowplot(font_size=10))
# pdf(file = "figures/dinoflagellate_results.pdf", width = 14, height = 15)
# plot_grid(title, plot_grid(dino_results$som_plot_list[[1]],dino_results$som_plot_list[[2]]), 
#           plot_grid(dino_results$rf_importance, dino_results$glm_subset, scale = c(0.7,0.7)),
#           dino_results$temp_plots,
#           nrow = 4, ncol = 1, rel_heights = c(0.1,1,1.2,1), labels = c("", "A","B","C"))
# dev.off()

shannon_depth_dist <- ggplot(dino_results$environ_dat, aes(x = dist_to_coast, y = CC_Depth, fill = shannon)) + 
  geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
  scale_fill_viridis(name = "Shannon Diversity", option = "D") +
  ylab("Depth (m)") + xlab("Distance to Coast (km)") 

even_depth_dist <- ggplot(dino_results$environ_dat, aes(x = dist_to_coast, y = CC_Depth, fill = evenness)) + 
  geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
  scale_fill_viridis(name = "Evenness", option = "D") +
  ylab("Depth (m)") + xlab("Distance to Coast (km)") 

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/dinoflagellate_diversity_evenness.pdf", width = 10, height = 13)
plot_grid(title, plot_grid(dino_results$shannon,dino_results$evenness, labels = c("A","B")),
          plot_grid(dino_results$temp_diversity_plots,
                    dino_results$temp_evenness_plots, ncol =  2,
                    labels = c("C","D")), plot_grid(shannon_depth_dist, even_depth_dist, labels = c("E","F")),
          plot_grid(dino_results$shannon_timeplot, dino_results$even_timeplot, ncol = 2, labels = c("G", "H")),
          nrow = 5, ncol = 1, rel_heights = c(0.1,1,1,1.2,0.8,0.8))
dev.off()

#Haptophytes

set.seed(4)

hapto_results <- group_analysis(in_amplicon = eighteen_s, group = taxa_id_list[[3]], group_name = taxa_groups[[3]],
                               environ_dat = full_dat, coeff_table = coeff_table, mean_table = mean_table)

title <- ggdraw() + draw_label("Haptophytes", fontface='bold')

theme_set(theme_cowplot(font_size=10))
# pdf(file = "figures/haptophyte_results.pdf", width = 14, height = 15)
# plot_grid(title, plot_grid(hapto_results$som_plot_list[[1]],hapto_results$som_plot_list[[2]]),
#           plot_grid(hapto_results$rf_importance, hapto_results$glm_subset, scale = c(0.7,0.7)),
#           hapto_results$temp_plots,
#           nrow = 4, ncol = 1, rel_heights = c(0.1,1,1.2,1), labels = c("", "A","B","C"))
# dev.off()

shannon_depth_dist <- ggplot(hapto_results$environ_dat, aes(x = dist_to_coast, y = CC_Depth, fill = shannon)) + 
  geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
  scale_fill_viridis(name = "Shannon Diversity", option = "D") +
  ylab("Depth (m)") + xlab("Distance to Coast (km)") 

even_depth_dist <- ggplot(hapto_results$environ_dat, aes(x = dist_to_coast, y = CC_Depth, fill = evenness)) + 
  geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
  scale_fill_viridis(name = "Evenness", option = "D") +
  ylab("Depth (m)") + xlab("Distance to Coast (km)") 

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/haptophyte_diversity_evenness.pdf", width = 10, height = 13)
plot_grid(title, plot_grid(hapto_results$shannon,hapto_results$evenness, labels = c("A","B")),
          plot_grid(hapto_results$temp_diversity_plots,
                    hapto_results$temp_evenness_plots, ncol =  2,
                    labels = c("C","D")), plot_grid(shannon_depth_dist, even_depth_dist, labels = c("E","F")),
          plot_grid(hapto_results$shannon_timeplot, hapto_results$even_timeplot, ncol = 2, labels = c("G", "H")),
          nrow = 5, ncol = 1, rel_heights = c(0.1,1,1,1.2,0.8,0.8))
dev.off()

#Crustacea

set.seed(452)

crust_results <- group_analysis(in_amplicon = eighteen_s, group = taxa_id_list[[4]], group_name = taxa_groups[[4]],
                                environ_dat = full_dat, coeff_table = coeff_table, mean_table = mean_table)

title <- ggdraw() + draw_label("Crustacea", fontface='bold')

theme_set(theme_cowplot(font_size=10))
# pdf(file = "figures/crustacea_results.pdf", width = 14, height = 15)
# plot_grid(title, plot_grid(crust_results$som_plot_list[[1]],crust_results$som_plot_list[[2]]), 
#           plot_grid(crust_results$rf_importance, crust_results$glm_subset, scale = c(0.7,0.7)),
#           crust_results$temp_plots,
#           nrow = 4, ncol = 1, rel_heights = c(0.1,1,1.2,1), labels = c("", "A","B","C"))
# dev.off()

shannon_depth_dist <- ggplot(crust_results$environ_dat, aes(x = dist_to_coast, y = CC_Depth, fill = shannon)) + 
  geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
  scale_fill_viridis(name = "Shannon Diversity", option = "D") +
  ylab("Depth (m)") + xlab("Distance to Coast (km)") 

even_depth_dist <- ggplot(crust_results$environ_dat, aes(x = dist_to_coast, y = CC_Depth, fill = evenness)) + 
  geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
  scale_fill_viridis(name = "Evenness", option = "D") +
  ylab("Depth (m)") + xlab("Distance to Coast (km)") 


theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/crustacea_diversity_evenness.pdf", width = 10, height = 13)
plot_grid(title, plot_grid(crust_results$shannon,crust_results$evenness, labels = c("A","B")),
          plot_grid(crust_results$temp_diversity_plots,
                    crust_results$temp_evenness_plots, ncol =  2,
                    labels = c("C","D")), plot_grid(shannon_depth_dist, even_depth_dist, labels = c("E","F")),
          plot_grid(crust_results$shannon_timeplot, crust_results$even_timeplot, ncol = 2, labels = c("G", "H")),
          nrow = 5, ncol = 1, rel_heights = c(0.1,1,1,1.2,0.8,0.8))
dev.off()

#Radiolaria

set.seed(56)

radio_results <- group_analysis(in_amplicon = eighteen_s, group = taxa_id_list[[5]], group_name = taxa_groups[[5]],
                                environ_dat = full_dat, coeff_table = coeff_table, mean_table = mean_table)

title <- ggdraw() + draw_label("Radiolaria", fontface='bold')

theme_set(theme_cowplot(font_size=10))
# pdf(file = "figures/radiolaria_results.pdf", width = 14, height = 15)
# plot_grid(title, plot_grid(radio_results$som_plot_list[[1]],radio_results$som_plot_list[[2]]), 
#           plot_grid(radio_results$rf_importance, radio_results$glm_subset, scale = c(0.7,0.7)),
#           radio_results$temp_plots,
#           nrow = 4, ncol = 1, rel_heights = c(0.1,1,1.2,1), labels = c("", "A","B","C"))
# dev.off()

shannon_depth_dist <- ggplot(radio_results$environ_dat, aes(x = dist_to_coast, y = CC_Depth, fill = shannon)) + 
  geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
  scale_fill_viridis(name = "Shannon Diversity", option = "D") +
  ylab("Depth (m)") + xlab("Distance to Coast (km)") 

even_depth_dist <- ggplot(radio_results$environ_dat, aes(x = dist_to_coast, y = CC_Depth, fill = evenness)) + 
  geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
  scale_fill_viridis(name = "Evenness", option = "D") +
  ylab("Depth (m)") + xlab("Distance to Coast (km)") 

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/radiolaria_diversity_evenness.pdf", width = 10, height = 13)
plot_grid(title, plot_grid(radio_results$shannon,radio_results$evenness, labels = c("A","B")),
          plot_grid(radio_results$temp_diversity_plots,
                    radio_results$temp_evenness_plots, ncol =  2,
                    labels = c("C","D")), plot_grid(shannon_depth_dist, even_depth_dist, labels = c("E","F")),
          plot_grid(radio_results$shannon_timeplot, radio_results$even_timeplot, ncol = 2, labels = c("G", "H")),
          nrow = 5, ncol = 1, rel_heights = c(0.1,1,1,1.2,0.8,0.8))
dev.off()

#Cnidaria

set.seed(981)

cnidaria_results <- group_analysis(in_amplicon = eighteen_s, group = taxa_id_list[[6]], group_name = taxa_groups[[6]],
                                environ_dat = full_dat, coeff_table = coeff_table, mean_table = mean_table)

title <- ggdraw() + draw_label("Cnidaria", fontface='bold')

theme_set(theme_cowplot(font_size=10))
# pdf(file = "figures/cnidaria_results.pdf", width = 14, height = 15)
# plot_grid(title, plot_grid(cnidaria_results$som_plot_list[[1]],cnidaria_results$som_plot_list[[2]]), 
#           plot_grid(cnidaria_results$rf_importance, cnidaria_results$glm_subset, scale = c(0.7,0.7)),
#           cnidaria_results$temp_plots,
#           nrow = 4, ncol = 1, rel_heights = c(0.1,1,1.2,1), labels = c("", "A","B","C"))
# dev.off()

shannon_depth_dist <- ggplot(cnidaria_results$environ_dat, aes(x = dist_to_coast, y = CC_Depth, fill = shannon)) + 
  geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
  scale_fill_viridis(name = "Shannon Diversity", option = "D") +
  ylab("Depth (m)") + xlab("Distance to Coast (km)") 

even_depth_dist <- ggplot(cnidaria_results$environ_dat, aes(x = dist_to_coast, y = CC_Depth, fill = evenness)) + 
  geom_point(size = 4, pch = 21, color = "black") + ylim(150,0) + 
  scale_fill_viridis(name = "Evenness", option = "D") +
  ylab("Depth (m)") + xlab("Distance to Coast (km)") 

theme_set(theme_cowplot(font_size=10))
pdf(file = "figures/cnidaria_diversity_evenness.pdf", width = 10, height = 10)
plot_grid(title, plot_grid(cnidaria_results$shannon,cnidaria_results$evenness, labels = c("A","B")),
          plot_grid(cnidaria_results$temp_diversity_plots,
                    cnidaria_results$temp_evenness_plots, ncol =  2,
                    labels = c("C","D")), plot_grid(shannon_depth_dist, even_depth_dist, labels = c("E","F")),
          plot_grid(cnidaria_results$shannon_timeplot, cnidaria_results$even_timeplot, ncol = 2, labels = c("G", "H")),
          nrow = 5, ncol = 1, rel_heights = c(0.1,1,1,1.2,0.8,0.8))
dev.off()

# Summary Plots

map <- map_data("world")  
map_dist <- map[,1:2]

summary_stats <- function(in_dat = diatom_results){

lat_lon_1 <- c(in_dat$centroid_list[[1]]@coords[1], in_dat$centroid_list[[1]]@coords[2])
lat_lon_2 <- c(in_dat$centroid_list[[2]]@coords[1], in_dat$centroid_list[[2]]@coords[2])

distances <- distm(lat_lon_1, lat_lon_2, fun = distGeo)
separation_clust <- round(distances/1000, digits = 0)

ordered <- in_dat$rf_out$importance[order(in_dat$rf_out$importance[,4], decreasing = TRUE),]

var_name <- paste0(rownames(ordered)[1],", ",rownames(ordered)[2])

mean_lm <- summary(lm(som_1 ~ temp_mean, data = in_dat$som_maps))
coeff_lm <- summary(lm(som_1 ~ temp_coeff, data = in_dat$som_maps))

lm_results <- c(round(mean_lm$r.squared,digits = 2), round(mean_lm$coefficients[2,4], digits = 4),
                round(coeff_lm$r.squared, digits = 2), round(coeff_lm$coefficients[2,4], digits = 4))

total_results <- c(separation_clust,lm_results, var_name)

return(total_results)

}

result_table <- as.data.frame(matrix(NA,6,6))

colnames(result_table) <- c("Distance between\n Centriod of Clusters (km)", "Mean SST\n R-Sq", "Mean SST\n p-val", "Coeff. Var.\n SST R-Sq", "Coeff. Var.\n SST p-val", "Most Important\n Random Forest Variables")

rownames(result_table) <- c("Diatoms", "Dinoflagellates", "Haptophytes", "Crustacea", "Radiolaria", "Cnidaria")

result_table[1,] <- summary_stats(in_dat = diatom_results)
result_table[2,] <- summary_stats(in_dat = dino_results)
result_table[3,] <- summary_stats(in_dat = hapto_results)
result_table[4,] <- summary_stats(in_dat = crust_results)
result_table[5,] <- summary_stats(in_dat = radio_results)
result_table[6,] <- summary_stats(in_dat = cnidaria_results)


g <- tableGrob(result_table)

# pdf(file = "figures/18sv9_summary_table_0723.pdf", width = 10, height = 4)
# plot.new()
# grid.draw(g)
# dev.off()


########## Relative Abundance Plots ##########

rel_abun_plots <- function(in_asv = diatom_results$asv_list, title = "Diatoms"){
  
  mean_som <- in_asv %>%
    group_by(som_id) %>%
    summarise_all(mean)
  
  diff <- abs(mean_som[1,]-mean_som[2,])
  diff <- diff[,-1]
  
  diff <- diff[,order(diff[1,], decreasing = TRUE)]
  
  mean_som <- gather(mean_som, "ASV", "mean_rel_abun", -som_id)
  
  mean_som$som_id[mean_som$som_id == 1] <- "som_1"
  mean_som$som_id[mean_som$som_id == 2] <- "som_2"
  
  colourCount = length(unique(mean_som$ASV))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  rel_plot <- ggplot(mean_som, aes(x = som_id, y = mean_rel_abun, fill = ASV)) + geom_col() + 
    theme(legend.position = "none") + scale_fill_manual(values = getPalette(colourCount)) +
    xlab("SOM Cluster") + ylab("Mean Relative Abundance") + ggtitle(title)
  
  print(rel_plot)
  
  # top taxa
  
  som1 <- mean_som[which(mean_som$som_id == "som_1"),]
  som2 <- mean_som[which(mean_som$som_id == "som_2"),]
  
  som1 <- som1[order(som1$mean_rel_abun, decreasing = TRUE),]
  som2 <- som2[order(som2$mean_rel_abun, decreasing = TRUE),]
  
  som1_taxa <- taxas[match(som1$ASV[1:20],taxas$Feature.ID),]
  som2_taxa <- taxas[match(som2$ASV[1:20],taxas$Feature.ID),]
  
  diff_taxa <- taxas[match(names(diff)[1:20],taxas$Feature.ID),]
  
  diff_taxa$diff <- as.numeric(diff[1:20])
  
  return(list(plot = rel_plot,
              mean_som = mean_som,
              diff = diff,
              som1_taxa = som1_taxa,
              som2_taxa = som2_taxa,
              diff_taxa = diff_taxa))
  
}

diatom_rel <- rel_abun_plots(in_asv = diatom_results$asv_list, title = "Diatoms")
dino_rel <- rel_abun_plots(in_asv = dino_results$asv_list, title = "Dinoflagellates")
hapto_rel <- rel_abun_plots(in_asv = hapto_results$asv_list, title = "Haptophytes")
crust_rel <- rel_abun_plots(in_asv = crust_results$asv_list, title = "Crustacea")
radio_rel <- rel_abun_plots(in_asv = radio_results$asv_list, title = "Radiolaria")
cnidaria_rel <- rel_abun_plots(in_asv = cnidaria_results$asv_list, title = "Cnidaria")

# summary plot

ASV_plot <- function(in_rel = diatom_rel, title = "Diatoms", n = 5){
  
  colourCount = length(unique(in_rel$mean_som$ASV))
  getPalette = colorRampPalette(brewer.pal(9, "Greys"))
  colors <- getPalette(colourCount)
  ASVs <- unique(in_rel$mean_som$ASV)
  
  #som 1
  som_1_palette <- colorRampPalette(c("deepskyblue", "darkblue"))
  in_rel$som1_taxa$color <- som_1_palette(n)
  som1_locations <- which(!is.na(match(ASVs, in_rel$som1_taxa$Feature.ID[1:n])))
  colors[som1_locations] <- in_rel$som1_taxa$color[match(ASVs, in_rel$som1_taxa$Feature.ID[1:n])[som1_locations]]
  
  #som 2
  som_2_palette <- colorRampPalette(c("indianred1", "red4"))
  in_rel$som2_taxa$color <- som_2_palette(n)
  som2_locations <- which(!is.na(match(ASVs, in_rel$som2_taxa$Feature.ID[1:n])))
  colors[som2_locations] <- in_rel$som2_taxa$color[match(ASVs, in_rel$som2_taxa$Feature.ID[1:n])[som2_locations]]
  
  both <- which(!is.na(match(in_rel$som1_taxa$Feature.ID[1:n], in_rel$som2_taxa$Feature.ID[1:n])))
  both2 <- which(!is.na(match(in_rel$som2_taxa$Feature.ID[1:n], in_rel$som1_taxa$Feature.ID[1:n])))
  
  if(length(both)>0){
    som_both_palette <- colorRampPalette(c("seagreen1", "darkgreen"))
    in_rel$som1_taxa$color[both] <- som_both_palette(length(both))
    in_rel$som2_taxa$color[both2] <- in_rel$som1_taxa$color[match(in_rel$som2_taxa$Feature.ID[both2], in_rel$som1_taxa$Feature.ID)]
    
    som_both_locations <- which(!is.na(match(ASVs, in_rel$som1_taxa$Feature.ID[1:n])))
    colors[som_both_locations] <- in_rel$som1_taxa$color[match(ASVs, in_rel$som1_taxa$Feature.ID[1:n])[som1_locations]]
  }
  
  
  in_rel_df <- in_rel$mean_som %>% 
    # Get frequency of "Swelling 1" within each level of Genotype
    # Order by frequency of "Swelling 1"
    arrange(mean_rel_abun)
  
  colors <- colors[match(unique(in_rel_df$ASV),unique(in_rel$mean_som$ASV))]
  
  in_rel_df$ASV <- factor(in_rel_df$ASV, levels = unique(in_rel_df$ASV))
  
  rel_plot <- ggplot(in_rel_df, aes(x = som_id, y = mean_rel_abun, fill = ASV)) + geom_col() + theme(legend.position = "none")  + scale_fill_manual(values = colors) +
    xlab("SOM Cluster") + ylab("Mean Relative Abundance") + 
    scale_x_discrete(labels = c("SOM 1", "SOM 2")) 
  
  som_table <- matrix("",5,2)
  colnames(som_table) <- c("SOM 1", "SOM 2")
  cols <- matrix("black", nrow(som_table), ncol(som_table))
  
  
  for (i in 1:n) {
    som1_str <- strsplit(in_rel$som1_taxa$Taxon[i], ";")
    som2_str <- strsplit(in_rel$som2_taxa$Taxon[i], ";")
    
    som_table[i,1] <- som1_str[[1]][length(som1_str[[1]])]
    som_table[i,2] <- som2_str[[1]][length(som2_str[[1]])]
    
    cols[i,1] <- in_rel$som1_taxa$color[i]
    cols[i,2] <- in_rel$som2_taxa$color[i]
    
  }
  
  bgs <- matrix("white", nrow(som_table), ncol(som_table))
  
  tt <- ttheme_default(core=list(fg_params = list(col = cols),
                                 bg_params = list(fill=bgs)),
                       rowhead=list(bg_params = list(col=NA)),
                       colhead=list(bg_params = list(col=NA)))
  
  table_plot <- tableGrob(som_table, theme = tt)
  
  title1 <- ggdraw() + draw_label(title, fontface='bold')
  
  plot_grid(title1,
            plot_grid(rel_plot, table_plot, ncol = 2, labels = c("A","B")),
            nrow = 2, rel_heights = c(0.1,1))
  
}


pdf(file = "figures/diatom_asvs.pdf", width = 10, height = 6)
ASV_plot(in_rel = diatom_rel, title = "Diatoms", n = 5)
dev.off()

pdf(file = "figures/dino_asvs.pdf", width = 10, height = 6)
ASV_plot(in_rel = dino_rel, title = "Dinoflagellates", n = 5)
dev.off()

pdf(file = "figures/hapto_asvs.pdf", width = 12, height = 6)
ASV_plot(in_rel = hapto_rel, title = "Haptophytes", n = 5)
dev.off()

pdf(file = "figures/crust_asvs.pdf", width = 10, height = 6)
ASV_plot(in_rel = crust_rel, title = "Crustacea", n = 5)
dev.off()

pdf(file = "figures/radio_asvs.pdf", width = 13, height = 6)
ASV_plot(in_rel = radio_rel, title = "Radiolaria", n = 5)
dev.off()

pdf(file = "figures/cnidaria_asvs.pdf", width = 10, height = 6)
ASV_plot(in_rel = cnidaria_rel, title = "Cnidaria", n = 5)
dev.off()

# Save output files

save(diatom_results, dino_results, hapto_results,
     crust_results, radio_results, cnidaria_results,
     diatom_rel, dino_rel, hapto_rel,
     crust_rel, radio_rel, cnidaria_rel,
     file = "output/18sv9_group_results.Rdata")
