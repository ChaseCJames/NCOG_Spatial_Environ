library(tidyverse)
library(VennDiagram)
library(RColorBrewer)
library(ggplotify)

myCol <- brewer.pal(3, "Pastel2")

# load data

load("data/18sv9_all.Rdata")
load("data/18sv9_tara_oceans.Rdata")
load("data/18sv9_tara_polar.Rdata")

eighteen_s <- eighteen_s[-which(substr(rownames(eighteen_s),10,11) < 76),]


calcofi <- names(which(colSums(eighteen_s) != 0))
tara <- names(which(colSums(tara_dat) != 0))
polar <- names(which(colSums(polar_dat) != 0))

venn_plot <- venn.diagram(x = list(tara, polar, calcofi),
  category.names = c("Tara Oceans" , "Tara Polar" , "CalCOFI"),
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

save(venn_plot, file = "output/venn_diagram.Rdata")
