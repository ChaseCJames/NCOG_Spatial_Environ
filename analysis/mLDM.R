library(mLDM)
library(corrplot)


load("data/16s_cyanos.Rdata")
load("output/cyano_16s_full_data.Rdata")

asv_short <- asv_table[,which(colSums(asv_table) > 20000)]

full_dat <- full_dat[match(rownames(asv_short), full_dat$eco_name),]

full_dat <- full_dat[,c(33,34,36,37,38,39,47,48)]

full_dat <- as.matrix(full_dat)

asv_short <- asv_short[complete.cases(full_dat),]
full_dat <- full_dat[complete.cases(full_dat),]

asv_short <- as.matrix(asv_short)

asv_short <- apply(asv_short, 2, as.numeric)
asv_short <- apply(asv_short, 2, as.integer)

#small test
# asv_short <- asv_short[1:600,]
# full_dat <- full_dat[1:600,]

test <- mLDM(X = asv_short, M = full_dat, Z_mean = 1, max_iteration = 10, max_iteration_B = 10)

resultOptimal <- test$optimal
B <- resultOptimal[[1]]
B0 <- resultOptimal[[2]]

Associations_OTU_OTU <- resultOptimal[[9]]
Associations_OTU_OTU <- as.data.frame(Associations_OTU_OTU)
colnames(Associations_OTU_OTU) <- c(1:48)
rownames(Associations_OTU_OTU) <- c(1:48)

Associations_EF_OTU <- resultOptimal[[8]]
Associations_EF_OTU <- as.data.frame(Associations_EF_OTU)
colnames(Associations_EF_OTU) <- c(1:48)
rownames(Associations_EF_OTU) <- colnames(full_dat)

diag(Associations_OTU_OTU) <- 0

corrplot(as.matrix(Associations_OTU_OTU), method = "color")
corrplot(as.matrix(Associations_EF_OTU), method = "color")
