library(tidyverse)
library(ncdf4)
library(chron)
library(lattice)
library(RColorBrewer)
library(rworldmap)
library(ggplot2)
library(ggmap)
library(cowplot)
library(reshape2)
library(metR)
library(oce)
library(zoo)
library(lubridate)
library(fpp2)

# GLOBAL OCEAN GRIDDED L4 SEA SURFACE HEIGHTS AND DERIVED VARIABLES REPROCESSED (1993-ONGOING)
# http://marine.copernicus.eu/services-portfolio/access-to-products/

# set path and filename

ncpath <- "data/"
ncname_CUTI <- "CUTI_daily" # use full dataset
ncname_BEUTI <- "BEUTI_daily" # use full dataset
ncfname_CUTI <- paste(ncpath, ncname_CUTI, ".nc", sep="")
ncfname_BEUTI <- paste(ncpath, ncname_BEUTI, ".nc", sep="")

dname_CUTI <- "CUTI"  
dname_BEUTI <- "BEUTI" 

ncin_CUTI <- nc_open(ncfname_CUTI)
ncin_BEUTI <- nc_open(ncfname_BEUTI)


# get lat and long

lat <- ncin_CUTI$var$CUTI$dim[[1]]$vals


# get time

year <- ncvar_get(ncin_CUTI,"year")
month <- ncvar_get(ncin_CUTI,"month")
month[which(month < 10)] <- paste0("0",month[which(month < 10)])
day <- ncvar_get(ncin_CUTI,"day")
day[which(day < 10)] <- paste0("0",day[which(day < 10)])

date <- paste0(day,"/",month,"/",year)

date <- as.Date(date, format = "%d/%m/%Y")

# Get vars

CUTI_array <- ncvar_get(ncin_CUTI,dname_CUTI)
BEUTI_array <- ncvar_get(ncin_BEUTI,dname_BEUTI)

dunits_CUTI <- ncatt_get(ncin_CUTI,dname_CUTI,"units")
fillvalue_CUTI <- ncatt_get(ncin_CUTI,dname_CUTI,"_FillValue")

dunits_BEUTI <- ncatt_get(ncin_BEUTI,dname_BEUTI,"units")
fillvalue_BEUTI <- ncatt_get(ncin_BEUTI,dname_BEUTI,"_FillValue")

# Convert Time

ts <- date

CUTI_array[CUTI_array==fillvalue_CUTI$value] <- NA
BEUTI_array[BEUTI_array==fillvalue_BEUTI$value] <- NA

socal_upwell <- which(lat < 36)
time_period <- which(year > 2013 & year < 2019)

CUTI_trim <- CUTI_array[socal_upwell,time_period]
BEUTI_trim <- BEUTI_array[socal_upwell,time_period]

# set negative BEUTI to 0 <- only interested in upwelling

BEUTI_trim[which(BEUTI_trim < 0)] <- 0

regional_CUTI <- colMeans(CUTI_trim)
regional_BEUTI <- colMeans(BEUTI_trim)
regional_nitrate <- colMeans(BEUTI_trim / CUTI_trim)

index_vals <- as.data.frame(matrix(NA,length(regional_BEUTI)*3,3))
colnames(index_vals) <- c("Date","Index","Value")

index_vals$Date <- rep(date[time_period],3)

index_vals$Value <- c(regional_CUTI, regional_BEUTI, regional_nitrate)

index_vals$Index <- c(rep("CUTI", length(regional_CUTI)), rep("BEUTI", length(regional_BEUTI)),
                      rep("Nitrate", length(regional_nitrate)))

cuti_df <- filter(index_vals, Index == "CUTI")

cuti_df$month2_MA <- rollapply(cuti_df$Value, width = 60, fill = NA, FUN = function(x) mean(x, na.rm = TRUE))
cuti_df$Year <- substr(cuti_df$Date,1,4)
cuti_df$Year_Day <- yday(cuti_df$Date)
cuti_df$Phase <- rep(NA,nrow(cuti_df))
cuti_df$Phase[which(!is.na(match(cuti_df$Year, c(2014,2015,2016))))] <- "Early"
cuti_df$Phase[which(is.na(cuti_df$Phase))] <- "Late"

cuti_plot <- ggplot(cuti_df, aes(x = Date, y = Value, color = Index)) + geom_line(color = "#F8766D") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none") +
  ylab("Coastal Upwelling Transport Index\n(CUTI)")

yearly_cuti <- ggplot(cuti_df, aes(x = Year_Day, y = month2_MA, group = interaction(Year, Phase), color = Phase)) +
  geom_line(size = 1) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none") +
  ylab("Coastal Upwelling Transport Index\n(CUTI)") + xlab("Year Day")

beuti_df <- filter(index_vals, Index == "BEUTI")

beuti_df$month2_MA <- rollapply(beuti_df$Value, width = 60, fill = NA, FUN = function(x) mean(x, na.rm = TRUE))
beuti_df$Year <- substr(beuti_df$Date,1,4)
beuti_df$Year_Day <- yday(beuti_df$Date)
beuti_df$Phase <- rep(NA,nrow(beuti_df))
beuti_df$Phase[which(!is.na(match(beuti_df$Year, c(2014,2015,2016))))] <- "Early"
beuti_df$Phase[which(is.na(beuti_df$Phase))] <- "Late"

beuti_plot <- ggplot(beuti_df, aes(x = Date, y = Value, color = Index)) + geom_line(color = "#00BFC4") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  ylab("Biologically Effective Upwelling Transport Index\n(BEUTI)")

yearly_beuti <- ggplot(beuti_df, aes(x = Year_Day, y = month2_MA, group = interaction(Year, Phase), color = Phase)) +
  geom_line(size = 1) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none") +
  ylab("Biologically Effective Upwelling Transport Index\n(BEUTI)") + xlab("Year Day")

nitrate_df <- filter(index_vals, Index == "Nitrate")

nitrate_df$month2_MA <- rollapply(nitrate_df$Value, width = 60, fill = NA, FUN = function(x) mean(x, na.rm = TRUE))
nitrate_df$Year <- substr(nitrate_df$Date,1,4)
nitrate_df$Year_Day <- yday(nitrate_df$Date)
nitrate_df$Phase <- rep(NA,nrow(nitrate_df))
nitrate_df$Phase[which(!is.na(match(nitrate_df$Year, c(2014,2015,2016))))] <- "Early"
nitrate_df$Phase[which(is.na(nitrate_df$Phase))] <- "Late"

nitrate_plot <- ggplot(nitrate_df, aes(x = Date, y = Value, color = Index)) + geom_line(color = "#7CAE00") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  ylab("Regionally Availible Nitrate\n(BEUTI/CUTI)")

yearly_nitrate <- ggplot(nitrate_df, aes(x = Year_Day, y = month2_MA, group = interaction(Year, Phase), color = Phase)) +
  geom_line(size = 1) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  ylab("Regionally Availible Nitrate\n(BEUTI/CUTI)") + xlab("Year Day")

yearly_nitrate1 <- ggplot(nitrate_df, aes(x = Year_Day, y = month2_MA, group = interaction(Year, Phase), color = Phase)) +
  geom_line(size = 1) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "bottom") +
  ylab("Regionally Availible Nitrate\n(BEUTI/CUTI)") + xlab("Year Day")

legend_1 <- get_legend(yearly_nitrate)
legend_2 <- get_legend(yearly_nitrate1)

yearly_nitrate <- yearly_nitrate + theme(legend.position = "none")

# pdf(file = "figures/coastal_upwelling_indicies.pdf", width = 12, height = 4)
# print(plot_grid(cuti_plot, beuti_plot, nitrate_plot, ncol = 3))
# dev.off()

pdf("figures/yearly_upwelling_index.pdf", width = 12, height = 4)
plot_grid(yearly_cuti,yearly_beuti,yearly_nitrate,legend_1,
          ncol = 4, rel_widths = c(1,1,1,0.2))
dev.off()

index_vals_mat <- as.data.frame(matrix(NA,length(regional_BEUTI),4))

colnames(index_vals_mat) <- c("Date", "CUTI", "BEUTI", "Nitrate")

index_vals_mat$Date <- date[time_period]
index_vals_mat$CUTI <- regional_CUTI
index_vals_mat$BEUTI <- regional_BEUTI
index_vals_mat$Nitrate <- regional_nitrate

# overall time series plots

cuti_ma <- ggplot() + 
  geom_point(data = cuti_df, aes(x = Date, y = Value), color = "grey", alpha = 0.5) +
  geom_line(data = cuti_df, aes(x = Date, y = month2_MA), color = "#F8766D", size = 1) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none") +
  ylab("Coastal Upwelling Transport Index\n(CUTI)") +
  scale_x_date(breaks = pretty_breaks(6)) + xlab("")

beuti_ma <- ggplot() +
  geom_point(data = beuti_df, aes(x = Date, y = Value), color = "grey", alpha = 0.5) +
  geom_line(data = beuti_df, aes(x = Date, y = month2_MA), color = "#00BFC4", size = 1) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none") +
  ylab("Biologically Effective Upwelling Transport Index\n(BEUTI)") + 
  scale_x_date(breaks = pretty_breaks(6)) + xlab("")

nitrate_ma <- ggplot() +
  geom_point(data = nitrate_df, aes(x = Date, y = Value), color = "grey", alpha = 0.5) +
  geom_line(data = nitrate_df, aes(x = Date, y = month2_MA), color = "#7CAE00", size = 1) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none") +
  ylab("Regionally Availible Nitrate\n(BEUTI/CUTI)") + 
  scale_x_date(breaks = pretty_breaks(6))

title_plot <- ggdraw() + draw_label("2-Month Moving Averages\nfor Upwelling Indicies", fontface='bold')

pdf("figures/upwelling_indicies_time.pdf", width = 5, height = 12)
plot_grid(title_plot, cuti_ma, beuti_ma, nitrate_ma, nrow = 4,
          rel_heights = c(0.2,1,1,1))
dev.off()

cuti_ma_bw <- ggplot() + 
  geom_point(data = cuti_df, aes(x = Date, y = Value), color = "grey", alpha = 0.5) +
  geom_line(data = cuti_df, aes(x = Date, y = month2_MA), color = "black", size = 1) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none") +
  ylab("Coastal Upwelling Transport Index\n(CUTI)") +
  scale_x_date(breaks = pretty_breaks(6)) + xlab("")

beuti_ma_bw <- ggplot() +
  geom_point(data = beuti_df, aes(x = Date, y = Value), color = "grey", alpha = 0.5) +
  geom_line(data = beuti_df, aes(x = Date, y = month2_MA), color = "black", size = 1) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none") +
  ylab("Biologically Effective Upwelling Transport Index\n(BEUTI)") + 
  scale_x_date(breaks = pretty_breaks(6)) + xlab("")

nitrate_ma_bw <- ggplot() +
  geom_point(data = nitrate_df, aes(x = Date, y = Value), color = "grey", alpha = 0.5) +
  geom_line(data = nitrate_df, aes(x = Date, y = month2_MA), color = "black", size = 1) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        legend.position = "none") +
  ylab("Regionally Availible Nitrate\n(BEUTI/CUTI)") + 
  scale_x_date(breaks = pretty_breaks(6))

side_panel_a <- plot_grid(title_plot, cuti_ma_bw, beuti_ma_bw, nitrate_ma_bw, nrow = 5,
                        rel_heights = c(0.2,1,1,1,0.1))

title_plot2 <- ggdraw() + draw_label("Yearly Upwelling Indicies\nby Phase", fontface='bold')

yearly_cuti <- yearly_cuti + xlab("") + ylab("")
yearly_beuti <- yearly_beuti + xlab("") + ylab("")
yearly_nitrate <- yearly_nitrate  + ylab("")

side_panel_b <- plot_grid(title_plot2,yearly_cuti,yearly_beuti,yearly_nitrate,legend_2,
          nrow = 5, rel_heights = c(0.2,1,1,1,0.1))

pdf("figures/full_upwelling_figure.pdf", width = 12, height = 12)
plot_grid(side_panel_a, side_panel_b, labels = c("A", "B"), ncol = 2)
dev.off()

save(index_vals_mat, file = "output/upwelling_indicies.Rdata")
