library(tidyverse)
library(VennDiagram)
library(RColorBrewer)
library(ggplotify)

                          surf = "dcm"
                          in_data = "data/18sv9_all.Rdata"
                          phys_dat = "data/NCOG_sample_log_DNA_meta_2014-2019.csv"
                          file_name = "figures/tara_ultra_overlap_dcm.pdf"
                          
  
  load("data/18sv9_tara_oceans.Rdata")
  load("data/18sv9_tara_polar.Rdata")
  load(in_data)
  phys_df <- read.csv(phys_dat, header = TRUE)
  
  asv_table <- eighteen_s
  
  splits <- strsplit(rownames(asv_table), "_")
  
  stations <- vector()
  depths <- vector()
  
  for (i in 1:length(splits)) {
    stations[i] <- paste0(splits[[i]][2]," ",splits[[i]][3])
    depths[i] <-splits[[i]][4]
  }
  
  
  if(surf == "surf"){
    phys_df <- phys_df %>%
      filter(Depthm < 15, as.numeric(substr(Sta_ID,1,3)) > 75)
    asv_table <- asv_table[which(depths < 15),]
    stations <- stations[which(depths < 15)]
  }
  if(surf == "dcm"){
    phys_df <- phys_df %>%
      filter(Depthm >= 15, as.numeric(substr(Sta_ID,1,3)) > 75)
    asv_table <- asv_table[which(depths > 15),]
    stations <- stations[which(depths > 15)]
  }
  
  if(surf == "both"){
    phys_df <- phys_df %>%
      filter(as.numeric(substr(Sta_ID,1,3)) > 75)
  }
  
 


  asv_table <- asv_table[which(as.numeric(substr(stations,1,3)) > 75),]
  stations <- stations[which(as.numeric(substr(stations,1,3)) > 75)]
  
  asv_table <- asv_table[,-which(colSums(asv_table, na.rm = TRUE) == 0)]
  
  matches <- which(!is.na(match(rownames(asv_table), paste0("X",phys_df$Sample.Name))))
  asv_table <- asv_table[matches,]
  
  rank_abun <- asv_table
  rank_abun[rank_abun > 0] <- 1
  ultra_rare <- names(which(colSums(rank_abun) == 1))
  
  tara <- names(which(colSums(tara_dat) != 0))
  polar <- names(which(colSums(polar_dat) != 0))
  
  venn_plot <- venn.diagram(x = list(tara, polar, ultra_rare),
                            category.names = c("Tara Oceans" , "Tara Polar" , "CalCOFI (Ultra Rare Taxa)"),
                            filename = NULL,
                            output=FALSE,
                            units = "in",
                            height = 6 , 
                            width = 6, 
                            resolution = 400,
                            compression = "lzw",
                            
                            # Circles
                            lwd = 1,
                            col=c("#440154ff", '#21908dff', '#fde725ff'),
                            fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
                            
                            
                            # Numbers
                            cex = 1.5,
                            fontface = "bold",
                            fontfamily = "sans",
                            
                            # Set names
                            cat.cex = 1.5,
                            cat.default.pos = "outer",
                            cat.pos = c(-20, 20, 130),
                            cat.dist = c(0.055, 0.055, 0.055),
                            cat.fontfamily = "sans",
                            cat.col = c("black", "black", "black"),
                            rotation = 1
  )
  
  pdf(file = file_name, width = 8, height = 8)
  grid.draw(venn_plot)
  dev.off()
