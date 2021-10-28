# Format for occupancy modeling
library(sf)
library(ggplot2)
library(reticulate)
setwd('/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM')

year <- 2018

AEAstring <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


warbler_array <- readRDS(paste0('warbler_', year, '_array.RDS'))
detection_array <- warbler_array$detection_array
species <- warbler_array$species
sites_prelim <-  st_as_sf(warbler_array$sites, coords = c('Longitude', 'Latitude'), crs = 4326)
sites_prelim$lon <- warbler_array$sites$Longitude
sites_prelim$lat <- warbler_array$sites$Latitude

# get elevations for all sites
# get ALOS elevations
use_condaenv('gee_interface', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication

ALOS <- ee$Image('JAXA/ALOS/AW3D30/V2_2')
ALOS_elev <- ALOS$select('AVE_DSM')
# Featurecollection of point geometries
geompts <- sapply(1:nrow(sites_prelim),function(x)ee$Geometry$Point(c(sites_prelim$lon[x],sites_prelim$lat[x])))
geompts <- ee$FeatureCollection(c(unlist(geompts)))
# Extract ALOS elevations for all points - combine into dataframe
pts_elev <- ALOS$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
ALOSelev <- sapply(c(1:length(pts_elev$features)),function(x)pts_elev$features[[x]]$properties$AVE_DSM)

sites_prelim$elev_ALOS <- ALOSelev
# Check that this looks right
ggplot(sites_prelim) + geom_sf(aes(col=elev_ALOS))

sites <- st_transform(sites_prelim, AEAstring)

rangemap_distances <- readRDS(paste0('rangemap_distances_2way_coastclip', year, '.RDS'))
distances <- rangemap_distances
# distances <- rangemap_distances$distances
# distances_updated <- rangemap_distances$distances_updated

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
# flattened_data$distance_updated <- as.vector(distances_updated)
flattened_data$elev <- rep(sites$elev_ALOS, nrow(species))
flattened_data$Q <- rowSums(flattened_data[,c(1:5)]) > 0
flattened_data$N <- rowSums(flattened_data[,c(1:5)])
flattened_data$distance_scaled <- flattened_data$distance/400000
p <- dd <- vector()
dev.off()
for(sp in unique(flattened_data$species)){
  fd_sp <- flattened_data#[flattened_data$species == sp, ]
  for(i in 1:50){
    j <- -2 + 7*i/50
    j0 <- j - 7/50
    p[i] <- mean(fd_sp$Q[fd_sp$distance_scaled>j0 & fd_sp$distance_scaled < j])
    dd[i] <- j
  }
  dd2 <- boot::inv.logit(dd*3)
  plot(boot::logit(p) ~ dd2, main = sp)
}
flattened_data$distance_transformed <- boot::inv.logit(flattened_data$distance_scaled*3)


saveRDS(flattened_data, file = paste0('flattened_data_', year, '.RDS'))


# Check distances
fdq <- flattened_data[flattened_data$Q == 1,]
hist(fdq$distance)
max(fdq$distance)
nrow(fdq)
fdq[fdq$distance > 50000, ]
