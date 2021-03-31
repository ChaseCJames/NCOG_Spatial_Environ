load("output/cyano_16s_full_data.Rdata")

full_dat <- full_dat[which(!is.na(full_dat$dist_to_coast)),]

full_dat <- full_dat %>%
  filter(as.numeric(substr(Sta_ID,1,3)) >= 76)

full_dat <- full_dat %>%
  filter(Depthm < 15)

coords_dat <- full_dat[,c(6,29,28)]
out <- t(as.matrix(as.data.frame(strsplit(coords_dat$Sta_ID, " "))))
out <- as.data.frame(apply(out, 2, as.numeric))
coords_dat <- bind_cols(out,coords_dat)

lines <- unique(coords_dat$V1)

boundary_mat <- coords_dat %>% filter(V1 == lines[1])
boundary_mat <- boundary_mat[c(which.min(boundary_mat$Lon_Dec), which.max(boundary_mat$Lon_Dec)),]

for (i in 2:length(lines)) {
  line_mat <- coords_dat %>% filter(V1 == lines[i])
  line_mat <- line_mat[c(which.min(line_mat$Lon_Dec), which.max(line_mat$Lon_Dec)),]
  boundary_mat <- bind_rows(boundary_mat, line_mat)
}

boundary_mat <- boundary_mat[c(13,14,12,10,8,6,4,2,1,3,5,7,11,13),]