library(tidyverse)
library(zoo)
library(lubridate)
library(broom)
library(patchwork)
library(ragg)
library(ggmap)

load("output/upwelling_indicies.Rdata")
load("output/total_full_data.Rdata")
ncd_data <- read.csv("data/ncd_data.csv")


full_dat <- full_dat %>% filter(substr(Sta_ID,2,3) > 75)

soms <- full_dat %>% group_by(som_id) %>% summarise(dist = mean(dist_to_coast, na.rm = TRUE)) 
near <- soms$som_id[which.min(soms$dist)]

cruise <- full_dat %>% group_by(Cruise) %>% summarise(mean_date = mean(Date, na.rm = TRUE),
                                                      min_date = min(Date, na.rm = TRUE),
                                                      max_date = max(Date, na.rm = TRUE),
                                                      prop_som = sum(som_id == near)/n())


ncd_data <- ncd_data %>% filter(CruiseAlias %in% cruise$Cruise, !is.na(NCDepth),
                                St_Line > 75)


slopes <- ncd_data %>% group_by(CruiseAlias) %>%
  do(fit_ncslope = tidy(lm(NCDepth ~ Distance.from.Shore, data = .))) %>% 
  unnest(fit_ncslope)

slopes <- slopes[seq(2,nrow(slopes), by = 2),]

cruise$nc_slope <- slopes$estimate[match(cruise$Cruise, slopes$CruiseAlias)]

ncd_data$Date <- dmy(paste0("01-",substr(ncd_data$CruiseAlias,5,6),"-",
                            substr(ncd_data$CruiseAlias,1,4)))


flat <- "201501"
steep <- "201907"

range_nc <- c(max(ncd_data$NCDepth, na.rm = TRUE), min(ncd_data$NCDepth, na.rm = TRUE))

map <- map_data("world") 

a <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = ncd_data %>% filter(CruiseAlias == flat), aes(x = DLon_Dec, y = DLat_Dec, fill = NCDepth),
             color = "black", size =2, stroke = 0.1, shape = 21) +
  scale_fill_gradient(name = "Mean Nitracline\nDepth (m)", low = "darkblue", high = "cyan",
                      trans = 'reverse', limits = c(100,10), oob = scales::squish,
                      breaks = c(10,25,50,75,100), labels = c("< 10","25","50","75","> 100")) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA,colour = "black", linetype = "solid")) +
  annotate(geom = "text", label = "Winter Cruise 2015", x = -121, y = 28.2, size = 2.5) +
  annotate(geom = "text", x = -126.95, y = 36.95, label = "a", fontface = "bold", size = 4) +
  scale_x_continuous(breaks = c(-125,-122,-119))
  

b <- ggplot() + 
  geom_polygon(data = map, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = c(-127, -116),ylim= c(28,37), 1.3) +
  xlab("Longitude") + ylab("Latitude") + 
  geom_point(data = ncd_data %>% filter(CruiseAlias == steep), aes(x = DLon_Dec, y = DLat_Dec, fill = NCDepth),
             color = "black", size =2, stroke = 0.1, shape = 21) +
  scale_fill_gradient(name = "Mean Nitracline\nDepth (m)", low = "darkblue", high = "cyan",
                      trans = 'reverse', limits = c(100,10), oob = scales::squish,
                      breaks = c(10,25,50,75,100), labels = c("< 10","25","50","75","> 100")) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA,colour = "black", linetype = "solid")) +
  annotate(geom = "text", label = "Summer Cruise 2019", x = -121, y = 28.2, size = 2.5) +
  annotate(geom = "text", x = -126.95, y = 36.95, label = "b", fontface = "bold", size = 4) +
  scale_x_continuous(breaks = c(-125,-122,-119))


yr_plot <- ggplot(ncd_data, aes(x = -Distance.from.Shore, y = NCDepth,
                     color = Date, group = Date)) +
  stat_smooth(method = "lm", fill = NA) + geom_point() +
  scale_y_reverse() +
  labs(x = "Distance from Shore (km)", y = "Nitracline Depth (m)") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  scale_color_date(low = "lightblue3", high = "darkblue") +
  annotate(geom = "text", x = -10, y = 0, label = "c", fontface = "bold", size = 4)


plot_nc <- (a + b + plot_layout(guides = "collect")) / yr_plot 

agg_png("figures/nitracline_examples.png", width = 8, height = 8, units = "in", res = 400)
plot(plot_nc)
dev.off()


# CUTI BEUTI NITRATE

up_mat <- matrix(NA,nrow(cruise),10) %>% as.data.frame()

colnames(up_mat) <- c("Cruise", "CUTI_mean", "BEUTI_mean", "Nitrate_mean",
                      "CUTI_tot", "BEUTI_tot", "Nitrate_tot",
                      "CUTI_2wk", "BEUTI_2wk", "Nitrate_2wk")

up_mat$Cruise <- cruise$Cruise

for (i in 1:nrow(up_mat)) {
  
  
  mean_dat <- cruise$mean_date[match(up_mat$Cruise[i], cruise$Cruise)]
  
  mean_val <- index_vals_mat %>% filter(Date >= (mean_dat %m+% months(-1)),
                            Date <= (mean_dat %m+% months(1))) 
  
  up_mat[i,2:4] <- colMeans(mean_val[,2:4], na.rm = TRUE)
  

  tot_val <- index_vals_mat %>% filter(Date >= cruise$min_date[match(up_mat$Cruise[i], cruise$Cruise)],
                                        Date <= cruise$max_date[match(up_mat$Cruise[i], cruise$Cruise)]) 
  
  up_mat[i,5:7] <- colMeans(tot_val[,2:4], na.rm = TRUE)
  
  wk2_val <- index_vals_mat %>% filter(Date >= (cruise$min_date[match(up_mat$Cruise[i], cruise$Cruise)] - 14),
                                       Date <= cruise$max_date[match(up_mat$Cruise[i], cruise$Cruise)]) 
  
  up_mat[i,8:10] <- colMeans(wk2_val[,2:4], na.rm = TRUE)
  
}

cruise_dat <- full_join(cruise, up_mat, by = "Cruise")



index_vals_mat$CUTI_mean <- rollapply(index_vals_mat$CUTI, width = 30, fill = NA, FUN = function(x) mean(x, na.rm = TRUE))
index_vals_mat$BEUTI_mean <- rollapply(index_vals_mat$BEUTI, width = 30, fill = NA, FUN = function(x) mean(x, na.rm = TRUE))
index_vals_mat$Nitrate_mean <- rollapply(index_vals_mat$Nitrate, width = 30, fill = NA, FUN = function(x) mean(x, na.rm = TRUE))

index_pivot <- index_vals_mat[,c(1,5:7)] %>% pivot_longer(-Date, names_to = "index", values_to = "value", values_drop_na = TRUE)

cruise_pivot <- cruise_dat[,c(2,6:15)] %>% pivot_longer(-c(mean_date, nc_slope))

 nc_slope_fit <- -mean(cruise_dat$BEUTI_mean)/mean(cruise_dat$nc_slope)
 
 cols <- c("Daily BEUTI"="grey80","Mean BEUTI"="#462796","Nitracline Slope"="#b02812")
 
p1 <- ggplot(data = cruise_dat, aes(x = mean_date)) +
  geom_line(data = index_pivot %>% filter(index == "BEUTI_mean", Date < dmy("01-01-2020")),
            aes(x = Date, y = value, color = "Daily BEUTI"), size = 1) +
  geom_point(aes(y = BEUTI_mean, color = "Mean BEUTI"), size = 5) +
  geom_line(aes(y = BEUTI_mean, color = "Mean BEUTI"), size = 1) +
  geom_point(aes(y = -nc_slope*nc_slope_fit, color = "Nitracline Slope"), size = 5) +
  geom_line(aes(y = -nc_slope*nc_slope_fit, color = "Nitracline Slope"), size = 1) +
  scale_y_continuous(name = "BEUTI\n(Biologically Effective Upwelling Trasport Index)",
                     sec.axis = sec_axis(~.*-1/nc_slope_fit,
                                         name = "Nitracline Slope (m/km)",
                                         breaks = c(0,-0.1,-0.2,-0.3, -0.4, -0.5),
                                         labels = c("0","-0.1","-0.2","-0.3", "-0.4", "<-0.5"))) +
  labs(x = "Date", color = "") + 
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        panel.grid.major = element_line(color = "grey92"),
        panel.grid.minor = element_line(color = "grey92"),
        legend.key = element_blank()) +
  scale_color_manual(values = cols) +
  annotate(geom = "text", x = dmy("01-01-2014"), y = 10.3, label = "d",
           fontface = "bold", size = 4)

 
p2 <- ggplot(data = cruise_dat, aes(x = BEUTI_mean, y = nc_slope)) +
  stat_smooth(method = "lm", color = "grey60") + geom_point(size = 3) +
  labs(x = "Mean BEUTI",
       y= "Nitracline Slope") +
  theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.title.x = element_text(color = "#462796"),
        axis.title.y = element_text(color = "#b02812"),
        plot.background = element_blank())

out_plot <- p1 + inset_element(p2, 0.01,0.6,0.52,0.9)

comb <- plot_nc | out_plot

pdf("figures_S/supp_fig_14_S.pdf", width = 12, height = 5)
plot(comb)
dev.off()











