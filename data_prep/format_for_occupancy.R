# Format for occupancy modeling
library(sf)
setwd('/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM')

AEAstring <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


warbler_2018_array <- readRDS('warbler_2018_array.RDS')
detection_array <- warbler_2018_array$detection_array
species <- warbler_2018_array$species
sites_prelim <-  st_as_sf(warbler_2018_array$sites, coords = c('Longitude', 'Latitude'), crs = 4326)
sites <- st_transform(sites_prelim, AEAstring)

rangemap_distances <- readRDS('rangemap_distances.RDS')
distances <- rangemap_distances$distances
distances_updated <- rangemap_distances$distances_updated

detection_slice <- list()
for(k in 1:nrow(species)){
  detection_slice[[k]] <- detection_array[,,k]
}
flattened_data <- as.data.frame(do.call(rbind, detection_slice))
names(flattened_data) <- c('v1', 'v2', 'v3', 'v4', 'v5')
flattened_data$species <- rep(species$code, each = nrow(sites))
flattened_data$sp_id <- rep(1:nrow(species), each = nrow(sites))
flattened_data$site <- rep(sites$routeID, nrow(species))
flattened_data$distance <- as.vector(distances)
flattened_data$distance_updated <- as.vector(distances_updated)
flattened_data$Q <- rowSums(flattened_data[,c(1,2,4,5)]) > 0

saveRDS(flattened_data, file = 'flattened_data.RDS')
