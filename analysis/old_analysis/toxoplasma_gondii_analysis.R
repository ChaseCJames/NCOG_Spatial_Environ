# Apicomplexa (Toxoplama gondii)

library(smooth)

apicomplexa <- split_taxa$ï..Feature.ID[which(split_taxa$C == "Apicomplexa")]

api_vals <- scaled_inputs[,match(apicomplexa, colnames(scaled_inputs))]

total_apis <- rowSums(api_vals, na.rm = TRUE)

full_dat$apicomplexa_totals <- total_apis[match(full_dat$eco_name, names(total_apis))]

# covariates / random forest

tox_glm <- full_dat[,c(64,33,37,38,39,61,62)]

tox_glm <- na.omit(tox_glm)

lm_tox <- lm(apicomplexa_totals ~ T_degC + PO4ug + SiO3ug + NO3ug + evenness + richness, data = tox_glm)

response_df <- tox_glm$apicomplexa_totals
predictor_df <- tox_glm[,2:7]

regsubsetsObj <- regsubsets(x=predictor_df ,y=response_df, nbest = 2, really.big = T)
regsubsetsObj$xnames <- c("(Intercept)", "Temp (C)", "PO4", "SiO3", "NO3", "Eveness", "Richness")

# pdf(file = "figures/18sv9_stepwise_regression_leaps_0723.pdf", width = 6, height = 5)
plot(regsubsetsObj, scale = "adjr2", main = "Backwards Subset Selection")
# dev.off()

# 
# glm_mean_var <- glm(apicomplexa_totals ~ lat + long + PO4_mean + NO3_mean + SiO3_mean + temp_mean + 
#       PO4_coeff + NO3_coeff + SiO3_coeff + temp_coeff, data = som_maps)
# 
# stepAIC(glm_mean_var)

model_AIC <- as.data.frame(matrix(NA,8,2))
colnames(model_AIC) <- c("Model","AIC")

# Temp

glm_mean_temp <- glm(apicomplexa_totals ~  T_degC, data = tox_glm)
mt_sum <- summary(glm_mean_temp)
model_AIC[1,2] <- mt_sum$aic
model_AIC[1,1] <- "Temp (C)"

# NO3

glm_mean_no3 <- glm(apicomplexa_totals ~  NO3ug, data = tox_glm)
mn_sum <- summary(glm_mean_no3)
model_AIC[2,2] <- mn_sum$aic
model_AIC[2,1] <- "NO3"

# PO4

glm_mean_po4 <- glm(apicomplexa_totals ~  PO4ug, data = tox_glm)
mp_sum <- summary(glm_mean_po4)
model_AIC[3,2] <- mp_sum$aic
model_AIC[3,1] <- "PO4"

# SiO3

glm_mean_sio3 <- glm(apicomplexa_totals ~  SiO3ug, data = tox_glm)
ms_sum <- summary(glm_mean_sio3)
model_AIC[4,2] <- ms_sum$aic
model_AIC[4,1] <- "SiO3"

# evenness

glm_mean_evenness <- glm(apicomplexa_totals ~  evenness, data = tox_glm)
ms_sum <- summary(glm_mean_evenness)
model_AIC[5,2] <- ms_sum$aic
model_AIC[5,1] <- "evenness"

# richness

glm_mean_richness <- glm(apicomplexa_totals ~  richness, data = tox_glm)
ms_sum <- summary(glm_mean_richness)
model_AIC[6,2] <- ms_sum$aic
model_AIC[6,1] <- "richness"

# Everything

glm_mean_var <- glm(apicomplexa_totals ~ T_degC + PO4ug + SiO3ug + NO3ug + evenness + richness, data = tox_glm)

all_sum <- summary(glm_mean_var)
model_AIC[7,2] <- all_sum$aic
model_AIC[7,1] <- "Full Model"

# best fit

stepAIC(glm_mean_var)

glm_simple <- glm(apicomplexa_totals ~ T_degC + SiO3ug + NO3ug +
                    richness, data = tox_glm)

simple_sum <- summary(glm_simple)
model_AIC[8,2] <- simple_sum$aic
model_AIC[8,1] <- "Temp + SiO3 +\n NO3 + richness"

AIC_table <- model_AIC[order(model_AIC$AIC, decreasing = TRUE),]

AIC_table$AIC <- round(AIC_table$AIC, 3)

pdf(file = "figures/toxoplasma_glm.pdf", height = 5, width = 7)
grid.table(AIC_table, rows = NULL)
dev.off()

####### Random Forests ##########

train <- sample(nrow(tox_glm), 0.7*nrow(tox_glm), replace = FALSE)

forest_out <- randomForest(apicomplexa_totals ~ T_degC + PO4ug + SiO3ug + NO3ug, data = tox_glm, importance = TRUE,
                           na.action = na.omit, subset = train, mtry = 3)


# Predicting on train set
predTrain <- predict(forest_out, tox_glm[train,])
# Checking classification accuracy
plot(predTrain, tox_glm$apicomplexa_totals[train], xlab = "Prediction", ylab = "Observation", xlim = c(0,0.1), ylim = c(0,0.1))
abline(a = 0, b = 1, col = "red", lty = 2)

# Predicting on validation set
predValid <- predict(forest_out, tox_glm[-train,])
# Checking classification accuracy
plot(predValid, tox_glm$apicomplexa_totals[-train], xlab = "Prediction", ylab = "Observation", xlim = c(0,0.05), ylim = c(0,0.05))
abline(a = 0, b = 1, col = "red", lty = 2)

varImpPlot(forest_out, type = 2, main = paste0("18sv9 Random Forest"))

pdf(file = "figures/toxoplasma_rf.pdf", width = 5, height = 4)
varImpPlot(forest_out, type = 2, main = paste0("Toxoplasma Random Forest"))
dev.off()

full_dat$Date <- as.Date(full_dat$Date, format = "%m/%d/%Y")

# plots

rainfall <- read.csv("data/rainfall.csv", header = FALSE, stringsAsFactors = FALSE)

means <- full_dat %>%
  group_by(Cruise) %>%
  summarise(mean_date = mean.Date(Date, na.rm = TRUE), mean_T = mean(T_degC, na.rm = TRUE),
            mean_NO3 = mean(NO3ug, na.rm = TRUE), mean_SiO3 = mean(SiO3ug, na.rm = TRUE),
            mean_PO4 = mean(PO4ug, na.rm = TRUE), mean_apicomplexa = mean(apicomplexa_totals, na.rm = TRUE))
  

pdf(file = "figures/apicomplexa_v_temp.pdf", width = 8, height = 6)
par(mar = c(5,5,2,5))
plot(full_dat$Date, full_dat$apicomplexa_totals, xlab = "Date", ylab = "Apicomplexa Relative Abundance", pch = 20, col = "black")

points(means$mean_date, means$mean_apicomplexa, type = "l", col = "black", lty = 2, lwd = 2)

par(new = T)
plot(means$mean_date, means$mean_T, type = "l", col = "blue", axes=F, xlab=NA, ylab=NA, lty = 2, ylim = c(10,22), lwd = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Temperature')

legend("topright",
       legend=c("Apicomplexa","Temperature"),
       lty=c(2,2), pch=c(NA, NA), col=c("black", "blue"))
dev.off()


# rainfall

pdf(file = "figures/apicomplexa_v_rain.pdf", width = 8, height = 6)
colnames(rainfall) <- c("Date", "Rainfall")
rainfall$Date <- as.Date(rainfall$Date, format = "%m/%d/%Y")

par(mar = c(5,5,2,5))
plot(full_dat$Date, full_dat$apicomplexa_totals, xlab = "Date", ylab = "Apicomplexa Relative Abundance", pch = 20, col = "black")

points(means$mean_date, means$mean_apicomplexa, type = "l", col = "black", lty = 2, lwd = 2)

par(new = T)
plot(rainfall$Date, rainfall$Rainfall/100, type = "l", col = "purple", axes=F, xlab=NA, ylab=NA, lty = 2, lwd = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Rainfall (inches)')

legend("topright",
       legend=c("Apicomplexa","Rainfall"),
       lty=c(2,2), lwd=c(2, 2), col=c("black", "purple"))
dev.off()
