# Import range maps, get distances of BBS records to species ranges, and explore result

# Loads object created by bbs_import.R
year <- 2018

library('sf')
library('ggplot2')
`%ni%` <- Negate(`%in%`)
setwd('/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM')
AEAstring <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

states <- raster::getData(country="USA", level=1)
provinces <- raster::getData(country="Canada", level=1)
estados <- raster::getData(country="Mexico", level=1)

states2 <- states[states$NAME_1 != "Hawaii", ]

states_sf <- st_as_sf(states2)
us_sf <- st_union(states_sf)

provinces_sf <- st_as_sf(provinces)
ca_sf <- st_union(provinces_sf)

estados_sf <- st_as_sf(estados)
mx_sf <- st_union(estados_sf)

n_am_sf <- st_union(us_sf, mx_sf)
n_am_sf <- st_union(n_am_sf, ca_sf)

n_am_AEA <- st_transform(n_am_sf, AEAstring)

n_am_buffin <- st_buffer(n_am_AEA, -30000)

saveRDS(n_am_buffin, "/Users/JacobSocolar/Desktop/n_am_buffin.Rdata")
n_am_buffin <- readRDS("/Users/JacobSocolar/Desktop/n_am_buffin.Rdata")

lakes <- st_read("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/great_lakes/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp")
big_lakes <- lakes[lakes$Continent == "North America" & lakes$Lake_area > 1000, ]
big_lakes <- st_union(st_make_valid(big_lakes))
big_lakes <- st_transform(big_lakes, AEAstring)

big_lakes_buffout <- st_buffer(big_lakes, 30000)

bc <- st_read("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/great_lakes/coastal_bc.kml")
la <- st_read("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/great_lakes/coastal_la.kml")
sl <- st_read("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/great_lakes/st_lawrence_snippet.kml")
extra_clip <- st_union(bc, la)
extra_clip <- st_union(extra_clip, sl)
extra_clip <- st_transform(extra_clip, AEAstring)

n_am_map <- st_difference(n_am_buffin, big_lakes_buffout)
n_am_map <- st_difference(n_am_map, extra_clip)

plot(n_am_map)
saveRDS(n_am_map, "/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/n_am_map_buffin.Rdata")



###########
# Load warbler BBS data, update taxonomy, and create sf object for sites
n_am_map <- readRDS("/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/n_am_map_buffin.Rdata")
warbler_array <- readRDS(paste0('warbler_', year, '_array.RDS'))
detection_array <- warbler_array$detection_array
species <- warbler_array$species

sites_prelim <-  st_as_sf(warbler_array$sites, 
                          coords = c('Longitude', 'Latitude'), crs = 4326)
sites <- st_transform(sites_prelim, AEAstring)

# Load recast_range_maps (see about.txt in enclosing folder for metadata); extract and project breeding ranges
load('/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/birdlife_maps/recast_range_maps.Rdata')
warbler_ranges <- recast_range_maps[recast_range_maps$SCINAME %in% species$SCINAME, ]
warbler_breeding_prelim <- warbler_ranges[warbler_ranges$SEASONAL %in% c(1,2), ]
warbler_breeding <- st_transform(warbler_breeding_prelim, AEAstring)

# Get distances from each point to the range of each species (based on BirdLife maps)
distances_allpoints <- matrix(data = NA, nrow = nrow(sites), ncol = nrow(species))
for(k in 1:nrow(species)){
  print(k)
  species_range <- st_union(st_make_valid(warbler_breeding[warbler_breeding$SCINAME == species$SCINAME[k], ]))
  species_border <- st_cast(species_range, to = "MULTILINESTRING")
  species_border_crop <- st_intersection(species_border, n_am_map)
  distances <- (as.numeric(st_distance(sites, species_range)) > 0) * as.numeric(st_distance(sites, species_border)) -
    (as.numeric(st_distance(sites, species_range)) <= 0) * as.numeric(st_distance(sites, species_border_crop))
    
# Positive distances for outside-of-range and negative distances for in-range.  Turns 
# out that as.numeric(st_distance(sites, species_range))>0) is 
# much faster than st_within(sites, species_range)
  distances_allpoints[, k] <- distances
}

# 
# # # Get distances between points with detections and range of each species (useful for data exploration; not neccesary for analysis)
# # # This could be sped up by subsetting distances_allpoints, but the script works fine as written and doesn't take to long, so I haven't messed with it
# # naive_Z <- apply(detection_array, MARGIN = c(1,3), 
# #                  FUN = function(x){return(sum(x) > 0)}) # Get detection/nondetection matrix (naive Z-matrix)
# # 
# # distances <- list()
# # for(k in 1:nrow(species)){
# #   species_sites <- sites[naive_Z[,k], ]
# #   species_range <- warbler_breeding[warbler_breeding$SCINAME == species$SCINAME[k], ]
# #   distance_matrix <- st_distance(species_sites, species_range)
# #   distances[[k]] <- apply(distance_matrix, 1, min)
# # }
# # distances2 <- unlist(distances)
# 
# # # Explore extralimital distances
# # sum(distances2>0)/length(distances2) # 7 percent of all site-records are outside of the mapped range
# # hist(distances2[distances2>0], breaks = 20)
# # hist(log(distances2[distances2>0]), breaks = 20) # note that the clustering of the logarithms is due to much larger areas available in larger bins
# # 
# # sum(distances2>50000)/length(distances2) # 2 percent
# # sum(distances2>100000)/length(distances2) # 0.5 percent
# # sum(distances2>150000)/length(distances2) # 0.2 percent
# 
# # max_dists <- data.frame(species = species$English, m5 = NA, m4 = NA, m3 = NA, m2 = NA, m1 = NA)
# # for(k in 1:nrow(species)){
# #   if(length(distances[[k]]) > 4){
# #     max_dists[k, 2:6] <- distances[[k]][order(distances[[k]])][(length(distances[[k]])-4):length(distances[[k]])]
# #   }else if(length(distances[[k]]) == 4){
# #     max_dists[k, 3:6] <- distances[[k]][order(distances[[k]])]
# #   }else if(length(distances[[k]]) == 3){
# #     max_dists[k, 4:6] <- distances[[k]][order(distances[[k]])]
# #   }else if(length(distances[[k]]) == 2){
# #     max_dists[k, 5:6] <- distances[[k]][order(distances[[k]])]
# #   }else if(length(distances[[k]]) == 1){
# #     max_dists[k, 6] <- distances[[k]]
# #   }
# # }
# #View(max_dists)
# 
# # # Plot extralimital detections
# # states_prelim <- USAboundaries::us_boundaries(type = "state")
# # states_prelim2 <- states_prelim[states_prelim$state_abbr %ni% c('HI', 'PR', 'DC', 'AK'), ]
# # states <- st_transform(states_prelim2, AEAstring)
# # colors <- c('gray90', 'red3')
# # sites_test <- sites
# # k <- 25
# # species$English[k]
# # sites_test$problem <- 0
# # sites_test$problem[which(naive_Z[,k] == 1)[distances[[k]] > 150000]] <- 1
# # plot(st_geometry(warbler_breeding[warbler_breeding$SCINAME == species$SCINAME[k], ]), col = 'gray95')
# # plot(st_geometry(states), add = T) 
# # plot(st_geometry(sites_test), add = T, col = colors[sites_test$problem + 1])
# 
# 
# ##### Update range maps #######
# wb <- warbler_breeding
# for(k in 1:nrow(species)){
#   update_dir <- paste0("Range_updates/", species$code[k])
#   kmz_files <- list.files(update_dir)
#   if(length(kmz_files) > 0){
#     for(m in 1:length(kmz_files)){
#       unzip(paste0(update_dir, '/', kmz_files[m]))
#       kml <- st_read('doc.kml')
#       file.remove('doc.kml')
#       newshape <- st_cast(st_transform(st_zm(kml), AEAstring), to = 'MULTIPOLYGON')
#       test <- st_sf(as.data.frame(matrix(rep(NA, 17), nrow = 1)), newshape$geometry)
#       names(test) <- names(wb)
#       st_geometry(test) <- 'Shape'
#       test$SCINAME <- species$SCINAME[k]
#       wb <- rbind(wb, test)
#     }
#   }
# }
# warbler_breeding_updated <- wb
# 
# distances_allpoints_updated <- matrix(data = NA, nrow = nrow(sites), ncol = nrow(species))
# for(k in 1:nrow(species)){
#   species_range <- warbler_breeding_updated[warbler_breeding_updated$SCINAME == species$SCINAME[k], ]
#   distance_matrix <- st_distance(sites, species_range)
#   distances_allpoints_updated[, k] <- apply(distance_matrix, 1, min)
# }
# 

# distances_updated <- list()
# for(k in 1:nrow(species)){
#   species_sites <- sites[naive_Z[,k], ]
#   species_range <- warbler_breeding_updated[warbler_breeding_updated$SCINAME == species$SCINAME[k], ]
#   distance_matrix <- st_distance(species_sites, species_range)
#   distances_updated[[k]] <- apply(distance_matrix, 1, min)
# }
# 
# distances2_updated <- unlist(distances_updated)

# # Explore extralimital distances
# sum(distances2_updated > 0)/length(distances2_updated) # 3.5 percent of all site-records are outside of the updated range
# hist(distances2_updated[distances2_updated > 0], breaks = 20)
# hist(log(distances2_updated[distances2_updated > 0]), breaks = 20) # note that the clustering of the logarithms is due to much larger areas available in larger bins
# 
# sum(distances2_updated > 50000)/length(distances2_updated) # 0.6 percent
# sum(distances2_updated > 100000)/length(distances2_updated) # 0.1 percent
# sum(distances2_updated > 150000)/length(distances2_updated) # 0.03 percent (3 records)

# max_dists_updated <- data.frame(species = species$English, m5 = NA, m4 = NA, m3 = NA, m2 = NA, m1 = NA)
# for(k in 1:nrow(species)){
#   if(length(distances_updated[[k]]) > 4){
#     max_dists_updated[k, 2:6] <- distances_updated[[k]][order(distances_updated[[k]])][(length(distances_updated[[k]])-4):length(distances_updated[[k]])]
#   }else if(length(distances_updated[[k]]) == 4){
#     max_dists_updated[k, 3:6] <- distances_updated[[k]][order(distances_updated[[k]])]
#   }else if(length(distances_updated[[k]]) == 3){
#     max_dists_updated[k, 4:6] <- distances_updated[[k]][order(distances_updated[[k]])]
#   }else if(length(distances_updated[[k]]) == 2){
#     max_dists_updated[k, 5:6] <- distances_updated[[k]][order(distances_updated[[k]])]
#   }else if(length(distances_updated[[k]]) == 1){
#     max_dists_updated[k, 6] <- distances_updated[[k]]
#   }
# }
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
# 
# rangemap_distances <- list(distances = distances_allpoints, 
#                            distances_updated = distances_allpoints_updated)
rangemap_distances <- distances_allpoints
saveRDS(rangemap_distances, file = paste0('rangemap_distances_2way_coastclip', year, '.RDS'))





