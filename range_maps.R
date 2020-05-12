# This script processes BBS and range map data in order to understand the extent of extralimital BBS records
# Loads object created by bbs_import.R

library('sf')
`%ni%` <- Negate(`%in%`)
setwd('/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM')
AEAstring <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# Load warbler BBS data, update taxonomy, and create sf object for sites
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

# Load recast_range_maps (see about.txt in enclosing folder for metadata); extract and project breeding ranges
load('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/birdlife_maps/recast_range_maps.Rdata')
warbler_ranges <- recast_range_maps[recast_range_maps$SCINAME %in% species$SCINAME, ]
warbler_breeding_prelim <- warbler_ranges[warbler_ranges$SEASONAL %in% c(1,2), ]
warbler_breeding <- st_transform(warbler_breeding_prelim, AEAstring)

# Get distances from each point to the range of each species (based on BirdLife maps)
distances_allpoints <- matrix(data = NA, nrow = nrow(sites), ncol = nrow(species))
for(k in 1:nrow(species)){
  species_range <- warbler_breeding[warbler_breeding$SCINAME == species$SCINAME[k], ]
  distance_matrix <- st_distance(sites, species_range)
  distances_allpoints[, k] <- apply(distance_matrix, 1, min)
}

# distances between points with detections and range of each species
naive_Z <- apply(detection_array, MARGIN = c(1,3), FUN = function(x){return(sum(x) > 0)})

distances <- list()
for(k in 1:nrow(species)){
  species_sites <- sites[naive_Z[,k], ]
  species_range <- warbler_breeding[warbler_breeding$SCINAME == species$SCINAME[k], ]
  distance_matrix <- st_distance(species_sites, species_range)
  distances[[k]] <- apply(distance_matrix, 1, min)
}

distances2 <- unlist(distances)
sum(distances2>0)/length(distances2) # 7 percent of all site-records are outside of the mapped range
hist(distances2[distances2>0], breaks = 20)
hist(log(distances2[distances2>0]), breaks = 20) # note that the clustering of the logarithms is due to much larger areas available in larger bins

sum(distances2>50000)/length(distances2) # 2 percent
sum(distances2>100000)/length(distances2) # 0.5 percent
sum(distances2>150000)/length(distances2) # 0.2 percent

max_dists <- data.frame(species = species$English, m5 = NA, m4 = NA, m3 = NA, m2 = NA, m1 = NA)
for(k in 1:nrow(species)){
  if(length(distances[[k]]) > 4){
    max_dists[k, 2:6] <- distances[[k]][order(distances[[k]])][(length(distances[[k]])-4):length(distances[[k]])]
  }else if(length(distances[[k]]) == 4){
    max_dists[k, 3:6] <- distances[[k]][order(distances[[k]])]
  }else if(length(distances[[k]]) == 3){
    max_dists[k, 4:6] <- distances[[k]][order(distances[[k]])]
  }else if(length(distances[[k]]) == 2){
    max_dists[k, 5:6] <- distances[[k]][order(distances[[k]])]
  }else if(length(distances[[k]]) == 1){
    max_dists[k, 6] <- distances[[k]]
  }
}

#View(max_dists)

## Plot extralimital detections
# states_prelim <- USAboundaries::us_boundaries(type = "state")
# states_prelim2 <- states_prelim[states_prelim$state_abbr %ni% c('HI', 'PR', 'DC', 'AK'), ]
# states <- st_transform(states_prelim2, AEAstring)
# colors <- c('gray90', 'red3')
# sites_test <- sites
# k <- 25
# species$English[k]
# sites_test$problem <- 0
# sites_test$problem[which(naive_Z[,k] == 1)[distances[[k]] > 150000]] <- 1
# plot(st_geometry(warbler_breeding[warbler_breeding$SCINAME == species$SCINAME[k], ]), col = 'gray95')
# plot(st_geometry(states), add = T) 
# plot(st_geometry(sites_test), add = T, col = colors[sites_test$problem + 1])


##### Update range maps #######
wb <- warbler_breeding
for(k in 1:nrow(species)){
  update_dir <- paste0("Range_updates/", species$code[k])
  kmz_files <- list.files(update_dir)
  if(length(kmz_files) > 0){
    for(m in 1:length(kmz_files)){
      unzip(paste0(update_dir, '/', kmz_files[m]))
      kml <- st_read('doc.kml')
      file.remove('doc.kml')
      newshape <- st_cast(st_transform(st_zm(kml), AEAstring), to = 'MULTIPOLYGON')
      test <- st_sf(as.data.frame(matrix(rep(NA, 17), nrow = 1)), newshape$geometry)
      names(test) <- names(wb)
      st_geometry(test) <- 'Shape'
      test$SCINAME <- species$SCINAME[k]
      wb <- rbind(wb, test)
    }
  }
}
warbler_breeding_updated <- wb

distances_allpoints_updated <- matrix(data = NA, nrow = nrow(sites), ncol = nrow(species))
for(k in 1:nrow(species)){
  species_range <- warbler_breeding_updated[warbler_breeding_updated$SCINAME == species$SCINAME[k], ]
  distance_matrix <- st_distance(sites, species_range)
  distances_allpoints_updated[, k] <- apply(distance_matrix, 1, min)
}


distances_updated <- list()
for(k in 1:nrow(species)){
  species_sites <- sites[naive_Z[,k], ]
  species_range <- warbler_breeding_updated[warbler_breeding_updated$SCINAME == species$SCINAME[k], ]
  distance_matrix <- st_distance(species_sites, species_range)
  distances_updated[[k]] <- apply(distance_matrix, 1, min)
}

distances2_updated <- unlist(distances_updated)
sum(distances2_updated > 0)/length(distances2_updated) # 3.5 percent of all site-records are outside of the updated range
hist(distances2_updated[distances2_updated > 0], breaks = 20)
hist(log(distances2_updated[distances2_updated > 0]), breaks = 20) # note that the clustering of the logarithms is due to much larger areas available in larger bins

sum(distances2_updated > 50000)/length(distances2_updated) # 0.6 percent
sum(distances2_updated > 100000)/length(distances2_updated) # 0.1 percent
sum(distances2_updated > 150000)/length(distances2_updated) # 0.03 percent (3 records)

max_dists_updated <- data.frame(species = species$English, m5 = NA, m4 = NA, m3 = NA, m2 = NA, m1 = NA)
for(k in 1:nrow(species)){
  if(length(distances_updated[[k]]) > 4){
    max_dists_updated[k, 2:6] <- distances_updated[[k]][order(distances_updated[[k]])][(length(distances_updated[[k]])-4):length(distances_updated[[k]])]
  }else if(length(distances_updated[[k]]) == 4){
    max_dists_updated[k, 3:6] <- distances_updated[[k]][order(distances_updated[[k]])]
  }else if(length(distances_updated[[k]]) == 3){
    max_dists_updated[k, 4:6] <- distances_updated[[k]][order(distances_updated[[k]])]
  }else if(length(distances_updated[[k]]) == 2){
    max_dists_updated[k, 5:6] <- distances_updated[[k]][order(distances_updated[[k]])]
  }else if(length(distances_updated[[k]]) == 1){
    max_dists_updated[k, 6] <- distances_updated[[k]]
  }
}

# View(max_dists_updated)

## Plot extralimital detections
# states_prelim <- USAboundaries::us_boundaries(type = "state")
# states_prelim2 <- states_prelim[states_prelim$state_abbr %ni% c('HI', 'PR', 'DC', 'AK'), ]
# states <- st_transform(states_prelim2, AEAstring)
# colors <- c('gray90', 'red3')
# sites_test <- sites
# k <- 25
# species$English[k]
# sites_test$problem <- 0
# sites_test$problem[which(naive_Z[,k] == 1)[distances_updated[[k]] > 150000]] <- 1
# plot(st_geometry(warbler_breeding_updated[warbler_breeding_updated$SCINAME == species$SCINAME[k], ]), col = 'gray95')
# plot(st_geometry(states), add = T) 
# plot(st_geometry(sites_test), add = T, col = colors[sites_test$problem + 1])

