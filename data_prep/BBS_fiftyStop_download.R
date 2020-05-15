# Download 50-stop BBS data, species file, and routes file

setwd('/Users/JacobSocolar/Desktop/useful_datasets/BBS/BBS_50_stop')

filenames <- paste0('Fifty', 1:10, '.zip')

for(i in 1:10){
  download.file(
    paste0('ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/50-StopData/1997ToPresent_SurveyWide/', filenames[i]),
    destfile = filenames[i])
  unzip(filenames[i])
  file.remove(filenames[i])
}

download.file('ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/SpeciesList.txt', 'SpeciesList.txt')
download.file('ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/routes.zip', 'routes.zip')
unzip('routes.zip')
file.remove('routes.zip')