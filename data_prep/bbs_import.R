# Import and explore BBS data for 2018

setwd('/Users/JacobSocolar/Desktop/useful_datasets/BBS/BBS_50_stop')
`%ni%` <- Negate(`%in%`)

##### 50-stop data #####
filenames <- paste0('Fifty', 1:10, '.csv')

bbs2018_list <- list()
for(i in 1:10){
  bbsImport <- read.csv(filenames[i])
  bbs2018_list[[i]] <- bbsImport[bbsImport$Year == 2018, ]
  rm(bbsImport)
}

colnames(bbs2018_list[[8]])[3] <- 'StateNum'  # typo in column names from BBS, with third column given as 'statenum' instead of 'StateNum'
bbs2018 <- do.call(rbind, bbs2018_list)

# restrict to continental US (this gives us more data-poor species by excluding Canadian detections)
bbs2018 <- bbs2018[(bbs2018$CountryNum == 840) & (bbs2018$StateNum != 3), ] 

# Aggregate data to five ten-stop blocks per route
bbs2018$b1 <- as.numeric(rowSums(bbs2018[,8:17]) > 0)
bbs2018$b2 <- as.numeric(rowSums(bbs2018[,18:27]) > 0)
bbs2018$b3 <- as.numeric(rowSums(bbs2018[,28:37]) > 0)
bbs2018$b4 <- as.numeric(rowSums(bbs2018[,38:47]) > 0)
bbs2018$b5 <- as.numeric(rowSums(bbs2018[,48:57]) > 0)
bbs2018 <- bbs2018[,c(1:7, 58:62)]

# Some species with just a few detections:
# KIWA: 1 detection (visit 2) on one route
# TRPA: 1 detection (visit 3) on one route
# BLPW: 1 detection (visit 3) on one route
# PAWA: 4 routes
# GCWA: 4 routes
# BBWA: 5 routes
# PARE: 7 routes
# CONW: 7 routes
# RFWA: 8 routes
# TEWA: 9 routes

##### Route data #####
bbs2018$routeID <- paste(bbs2018$CountryNum, bbs2018$StateNum, bbs2018$Route, 
                         sep = '_') # create unique route IDs
rte_list_2018 <- unique(bbs2018$routeID) # unique routes run in 2018

routes <- read.csv('routes.csv')
routes$routeID <- paste(routes$CountryNum, routes$StateNum, routes$Route, 
                        sep = '_') # unique route IDs
routes2018 <- routes[routes$routeID %in% rte_list_2018, ]
routes2018 <- routes2018[routes2018$RouteTypeDetailID == 1, ] # just use the random 50-stop routes

bbs2018 <- bbs2018[bbs2018$routeID %in% routes2018$routeID,]

##### Species data #####
species <- read.fwf('SpeciesList.txt', skip = 10, strip.white = TRUE, 
                    header = FALSE,
                    colClasses = c("integer", "integer", rep("character", 7)),
                    widths = c(6, -1, 5, -1, 50, -1, 50, -1, 50, -1, 50, -1, 
                               50, -1, 50, -1, 50),
                    fileEncoding = "iso-8859-1")
species <- species[,c(2,3,6,7,8,9)]
colnames(species) <- c('AOU', 'English', 'Order', 'Family', 'Genus', 'Species')

##### Get species-site-visit array for Parulids #####
warblers <- species[species$Family == "Parulidae", ]
sum(bbs2018$AOU == 06556) # All YRWA from 2018 are classified as MYWA or AUWA, which is convenient for HBW matching
warblers <- warblers[warblers$AOU %ni% c(6412, 6413, 6556, 6686, 6685), ] # Remove slashes and hybrids
warblers$Species[warblers$Species == 'coronata coronata'] <- 'coronata'
warblers$Species[warblers$Species == 'coronata audoboni'] <- 'auduboni'

# Add two species that regularly breed in lower 48 but never sampled in BBS
extra_warblers <- data.frame(AOU = NA, 
                             English = c("Colima Warbler", 
                                         "Rufous-capped Warbler"), 
                             Order = "Passeriformes",
                             Family = "Parulidae", 
                             Genus = c("Leiothlypis", "Basileuterus"), 
                             Species = c("crissalis", "rufifrons"))
warblers <- rbind(warblers, extra_warblers)
warblers$Genus <- gsub('Oreothlypis', 'Leiothlypis', warblers$Genus) # Update harmonize to range map taxonomy (BirdLife)
warblers$SCINAME <- paste(warblers$Genus, warblers$Species, sep = ' ')
warblers$code <- c("OVEN", "WEWA", "LOWA", "NOWA", "GWWA", "BWWA", "BAWW", 
                   "PROW", "SWWA", "TEWA", "OCWA", "LUWA", "NAWA", "VIWA", 
                   "CONW", "MGWA", "MOWA", "KEWA", "COYE", "HOWA", "AMRE", 
                   "KIWA", "CMWA", "CERW", "NOPA", "TRPA", "MAWA", "BBWA", 
                   "BLBW", "YEWA", "CSWA", "BLPW", "BTBW", "PAWA", "PIWA", 
                   "MYWA", "AUWA", "YTWA", "PRAW", "GRWA", "BTYW", "TOWA", 
                   "HEWA", "GCWA", "BTNW", "CAWA", "WIWA", "RFWA", "PARE", 
                   "COWA", "RCWA")
bbs2018_warblers <- bbs2018[bbs2018$AOU %in% warblers$AOU, ] # COWA and RCWA were never detected in BBS, so we don't need to worry about their codes.

detection_array <- array(data = 0, dim = c(nrow(routes2018), 5, nrow(warblers)))
for(k in 1:nrow(warblers)){
  spdata <- bbs2018_warblers[bbs2018_warblers$AOU == warblers$AOU[k], ]
  for(i in 1:nrow(routes2018)){
    if(routes2018$routeID[i] %in% spdata$routeID){
      detection_array[i, ,k] <- as.numeric(spdata[spdata$routeID == routes2018$routeID[i], 8:12])
    }
  }
}

warbler_2018_array <- list(detection_array = detection_array, sites = routes2018, species = warblers)
saveRDS(warbler_2018_array, file = '/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/warbler_2018_array.RDS')
