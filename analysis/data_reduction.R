# Data Reduction Methods

library(tidyverse)
library(ggbiplot)
library(phateR)



data_reduction_compare <- function(in_file = "data/16s_cyanos.Rdata"){
  
  load(in_file)
  
  scaled_inputs <- scaled_inputs[,apply(scaled_inputs, 2, var, na.rm=TRUE) != 0]
  
  pca.cyano <- prcomp(scaled_inputs, center = TRUE, scale. = FALSE)
  
  ggbiplot(pca.cyano, var.axes = FALSE)
  
}