# Import and explore BBS data for 2018

setwd('/Users/JacobSocolar/Desktop/useful_datasets/BBS/BBS_50_stop')
`%ni%` <- Negate(`%in%`)

##### 50-stop data #####
filenames <- paste0('Fifty', 1:10, '.csv')

year <- 2018

bbs_list <- list()
for(i in 1:10){
  bbsImport <- read.csv(filenames[i])
  bbs_list[[i]] <- bbsImport[bbsImport$Year == year, ]
  rm(bbsImport)
}

#2001 one strikeable TRPA
#2011 no TRPA, one strikeable BLPW
#2012 no TRPA, one strikeable BLPW, no GCWA
#2014 no GCWA, one strikeable TRPA

colnames(bbs_list[[8]])[3] <- 'StateNum'  # typo in column names from BBS, with third column given as 'statenum' instead of 'StateNum'
bbs <- do.call(rbind, bbs_list)

# restrict to continental US (this gives us more data-poor species by excluding Canadian detections)
bbs <- bbs[(bbs$CountryNum == 840) & (bbs$StateNum != 3), ] 

# Aggregate data to five ten-stop blocks per route
bbs$b1 <- as.numeric(rowSums(bbs[,8:17]) > 0)
bbs$b2 <- as.numeric(rowSums(bbs[,18:27]) > 0)
bbs$b3 <- as.numeric(rowSums(bbs[,28:37]) > 0)
bbs$b4 <- as.numeric(rowSums(bbs[,38:47]) > 0)
bbs$b5 <- as.numeric(rowSums(bbs[,48:57]) > 0)
bbs <- bbs[,c(1:7, 58:62)]

# Some species with just a few detections:
# KIWA:
bbs[bbs$AOU == 6700,]
# TRPA: 
bbs[bbs$AOU == 6490,]
# BLPW: 
bbs[bbs$AOU == 6610,]
# PAWA: 
bbs[bbs$AOU == 6720,]
# GCWA: 
bbs[bbs$AOU == 6660,]
# BBWA: 
bbs[bbs$AOU == 6600,]
# PARE: 
bbs[bbs$AOU == 6880,]
# CONW: 
bbs[bbs$AOU == 6780,]
# RFWA: 
bbs[bbs$AOU == 6900,]
# TEWA: 
bbs[bbs$AOU == 6470,]
# GOLW
bbs[bbs$AOU == 6920,]
# RCWA
bbs[bbs$AOU == 6921,]

##### Route data #####
bbs$routeID <- paste(bbs$CountryNum, bbs$StateNum, bbs$Route, 
                         sep = '_') # create unique route IDs
rte_list <- unique(bbs$routeID) # unique routes run in 2018

routes <- read.csv('routes.csv')
routes$routeID <- paste(routes$CountryNum, routes$StateNum, routes$Route, 
                        sep = '_') # unique route IDs
routes_year <- routes[routes$routeID %in% rte_list, ]
routes_year <- routes_year[routes_year$RouteTypeDetailID == 1, ] # just use the random 50-stop routes

bbs <- bbs[bbs$routeID %in% routes_year$routeID,]

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
sum(bbs$AOU == 06556) # All YRWA from 2018 are classified as MYWA or AUWA, which is convenient for HBW matching
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
bbs_warblers <- bbs[bbs$AOU %in% warblers$AOU, ] # COWA and RCWA were never detected in BBS, so we don't need to worry about their codes.

detection_array <- array(data = 0, dim = c(nrow(routes_year), 5, nrow(warblers)))
for(k in 1:nrow(warblers)){
  spdata <- bbs_warblers[bbs_warblers$AOU == warblers$AOU[k], ]
  spdata <- spdata[!duplicated(spdata$routeID),]
  for(i in 1:nrow(routes_year)){
    if(routes_year$routeID[i] %in% spdata$routeID){
      detection_array[i, ,k] <- as.numeric(spdata[spdata$routeID == routes_year$routeID[i], 8:12])
    }
  }
}

warbler_array <- list(detection_array = detection_array, sites = routes_year, species = warblers)
saveRDS(warbler_array, file = paste0('/Users/jacobsocolar/Dropbox/Work/Occupancy/biogeographicMSOM/warbler_', year, '_array.RDS'))
