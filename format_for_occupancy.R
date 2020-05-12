# Format for occupancy modeling
warbler_2018_array <- readRDS('warbler_2018_array.RDS')
detection_array <- warbler_2018_array$detection_array
species <- warbler_2018_array$species
species$SCINAME <- paste(species$Genus, species$Species, sep = ' ')
species$SCINAME <- gsub('Oreothlypis', 'Leiothlypis', species$SCINAME)
species$code <- c("OVEN", "WEWA", "LOWA", "NOWA", "GWWA", "BWWA", "BAWW", "PROW", "SWWA", "TEWA",
                  "OCWA", "LUWA", "NAWA", "VIWA", "CONW", "MGWA", "MOWA", "KEWA", "COYE", "HOWA",
                  "AMRE", "KIWA", "CMWA", "CERW", "NOPA", "TRPA", "MAWA", "BBWA", "BLBW", "YEWA",
                  "CSWA", "BLPW", "BTBW", "PAWA", "PIWA", "MYWA", "AUWA", "YTWA", "PRAW", "GRWA",
                  "BTYW", "TOWA", "HEWA", "GCWA", "BTNW", "CAWA", "WIWA", "RFWA", "PARE")
sites_prelim <-  st_as_sf(warbler_2018_array$sites, coords = c('Longitude', 'Latitude'), crs = 4326)
sites <- st_transform(sites_prelim, AEAstring)



detection_slice <- list()
for(k in 1:49){
  detection_slice[[k]] <- detection_array[,,k]
}
flattened_data <- as.data.frame(do.call(rbind, detection_slice))
names(flattened_data) <- c('v1', 'v2', 'v3', 'v4', 'v5')
flattened_data$species <- rep(species$code, each = nrow(sites))
flattened_data$sp_id <- rep(1:49, each = nrow(sites))
flattened_data$distance <- as.vector(distances_allpoints)
flattened_data$distance_updated <- as.vector(distances_allpoints_updated)
flattened_data$Q <- rowSums(flattened_data[,c(1,2,4,5)]) > 0