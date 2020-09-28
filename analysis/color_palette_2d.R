## colorPalette.R
## builds an (m × n) 2D palette
## by mixing 2 hues (col1, col2)
## and across two luminosities (lum1,lum2)
## returns a matrix of the hex RGB values

rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)

makePalette <- function(col1 = RGB(0,0,0),
                        col2 = RGB(0,0,255),
                        col1.1 = RGB(0,0,0),
                        col2.2 = RGB(255,0,0),m,n,...) {
  C <- matrix(data=NA,ncol=m,nrow=n)
  alpha <- seq(0,1,length.out=m)
  
    for (j in 1:m) {
      for (i in 1:n) {
        
        color_1 <- mixcolor(alpha[i], col1, col2)
        color_2 <- mixcolor(alpha[j], col1.1, col2.2)
        
        c1_cord <- coords(color_1)
        c2_cord <- coords(color_2)
        color_3 <- RGB(c1_cord+c2_cord)

        out <- coords(color_3)
        
        hexc <- rgb2hex(out[1,1], out[1,2], out[1,3])
        C[i,j] <- hexc
      }
    }


  return(C)
}

## plot a vector or matrix of RGB colors
plotPalette <- function(C) {

  pal_df <- gather(as.data.frame(C))
  pal_df$x <- rep(1:10, each = 10)
  pal_df$y <- rep(1:10, times = 10)
  
  out_plot <- ggplot(pal_df, aes(x = x, y = y, fill = value)) +
    geom_tile(show.legend = FALSE) +
    scale_fill_manual(values = pal_df$value) + coord_fixed() +
    theme(panel.background = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    ylab("PCA 2") + xlab("PCA 1")
  
  return(out_plot)
  
}

## Let's put it all together.
## We make two colors in the LAB space, and then plot a 2D palette
## going from 60 to 25 luminosity values.
library(colorspace)
C <- makePalette(m=10, n=10)
out_pal <- plotPalette(C)
print(out_pal)

save(out_pal, file = "output/2dcolorpal.Rdata")
