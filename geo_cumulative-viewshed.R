# abstract = Generates raster boolean layer based on cumulative viewshed of multi-point layer (i.e. union of each instance)
# Inp: Layer of points
#      DEM
# Out: Raster of cumulative viewshed of all points
# Author: Julian Rosser
###############################################################################

# wps.des: geo.cumulative-viewshed;
# wps.in: data, application/x-zipped-shp;
# wps.in: dem, type = tiff;

library(rgeos); library(maptools); library(rgdal); library(spatstat); library(spgrass6); library(raster); library(sp);

#Global variables for type of processing
outputRasterIsViewshedCount = 1; #0 = simple logical union of viewsheds. 1 = cumulative viewshed (viewshed count)
maxDistance = 500; #maximum viewing distance in metres for each geo-point.
writeViewshedsToDisk = 0; #for debugging after each viewshed calcuation by GRASS



#function to add up the viewshed layers
mapCalc <- function(los, productRaster){
  #these conditions seem to be needed for lidar data
  if (any(is.na(los[[1]])) == TRUE) {
    losValues <- ifelse(is.na(los@data[[1]]), 0, 1) #if na values, set to 0, else set to 1
  } else {
    losValues = los@data[[1]]
  }
  sum <- losValues + productRaster
  return(sum)
}


#Read in inputs from WPS execute request, this must come before any wd changes!
ViewPoints = readShapePoints(data)
baseElevationModel = readGDAL(dem)



#GRASS initialisation
#setwd("C:/tmp_r_grass") #manual setting R and GRASS processing directory, don't use this for deployment
#Remember when called in WPS execution mode getwd() will return the temporary WPSR dir.
loc <- initGRASS("C:/programme/GRASS7.0.svn",home=getwd(), gisDbase="GRASS_TEMP", override=TRUE )


initialWd = getwd()



# wps.off; #Hard-coded test data

inputPointsDir = "C:/Users/ezzjfr/post_doc/data/spatial_accuracy/"
inputPointsFilename ="Flickr_for_spatial_accuracy"

baseDemFilename ="base_dem.tif"
initialWdBaseDem= paste0(initialWd,"/", baseDemFilename)

setwd(inputPointsDir)

#shpFiles = list.files(inputPointsDir,full.names=FALSE, glob2rx(paste0(inputPointsFilename,".*")))
shpFiles = list.files(inputPointsDir,full.names=FALSE, glob2rx(paste0(inputPointsFilename,".*")))

#file.copy(from = shpFiles, to = initialWd,  overwrite=TRUE)
ViewPointsDataFrame = readOGR(dsn=paste0(inputPointsFilename, ".shp"), layer=inputPointsFilename)

# wps.on;

setwd(initialWd )

#Read in the raster data
# wps.off;
file.copy(from = "C:/Users/ezzjfr/post_doc/data/spatial_accuracy/Lidar_DTM_downsampled_to_5m_spatial_accuracy.tif", to = initialWdBaseDem  ,overwrite = TRUE )


baseElevationModel = readGDAL(baseDemFilename)
# wps.on;

#get projection
demProjection = proj4string(baseElevationModel)

execGRASS("r.in.gdal", flags="o", parameters=list(input=baseDemFilename, output="DEM"))
execGRASS("g.region", parameters=list(raster="DEM"))


# GRID CREATING & GETTING THE NUMBER OF CELLS(for further programming)
grd <- gmeta2grd()
ncells <- grd@cells.dim[1]*grd@cells.dim[2]

bottomLeftX = grd@cellcentre.offset[1]
bottomLeftY = grd@cellcentre.offset[2]

demBbox = bbox(baseElevationModel)

demExtentPoly = as(extent(baseElevationModel), "SpatialPolygons")
proj4string(demExtentPoly) = proj4string(ViewPointsDataFrame)
#crop the points to the dem extent
ViewPoints = crop(ViewPointsDataFrame,demExtentPoly)


#ViewPoints prep
aFrame = data.frame(ViewPoints) # create data frame
#coords = aFrame[2:3] # take the coordinate cols
coordsX = aFrame[["coords.x1"]]
coordsY = aFrame[["coords.x2"]]
coords = cbind(coordsX,coordsY)

# GRID VALUES INITIALIZATION before LOOPING
sumV <- rep(0, ncells)

uniqueOwners = unique(aFrame$ownerNum)


#First viewshed
execGRASS("r.viewshed", parameters = list(input = "DEM", output = "cumulativeViewshed", max_distance= maxDistance, coordinates = as.integer(coords[1,])), flags = c("overwrite" , "b"))
los <- readRAST6("cumulativeViewshed")
#add this firest layer into the cumulative raster
sumV = mapCalc(los,sumV)

# Now loop over any remaining points, creating viewsheds and adding into cumulative viewshed with logical union mapcalc
for (i in 2:nrow(coords)) {
#for (i in 2:4) {
  losLayerName = toString(i)
  #initial cleaning
  #execGRASS("g.remove", parameters = list(type = "raster",name=losLayerName), flags="f")
  #rm(losInR)

  #Compute this shed instance
  execGRASS("r.viewshed", parameters = list(input = "DEM", output = losLayerName, max_distance=maxDistance, coordinates = as.integer(coords[i,])), flags = c("overwrite" , "b"))

  los <- readRAST6(losLayerName)

  #Write this viewshedinstance to disk
  if (writeViewshedsToDisk == 1)
    writeGDAL(los, paste(i,".tiff"), drivername="GTiff", type="Float32")
  end

  #add it to the cumulative viewshed with map calc
  sumV = mapCalc(los,sumV)
}

# 0 VALUES TO NA
#sumV[sumV==0]<-NA #ARC doesn't like this, probably encodes nulls weirdly
#save(sumV, file="sumV.RData")


#############################
### MAPPING DATA TO GRID  ###
#############################
sgdf <- SpatialGridDataFrame(grd, data = data.frame(sum=sumV))
proj4string(sgdf) = demProjection

# EXPORT DATA FILLED GRID TO TIFF
outFile = "cumulativeViewshed.tif";
out = writeGDAL(sgdf["sum"], outFile, drivername="GTiff", type="Float32")


#clean up GRASS junk
unlink(paste(getwd(), "GRASS_TEMP", sep="/"), recursive=TRUE)
file.remove(paste(getwd(), ".grassrc6", sep="/"))

#list the layers in grass workspace
#execGRASS("g.list", parameters=list(type="rast"), intern=TRUE)

#Export 52N WPS using annotations
# wps.out: id = out, type = tiff, title = a cumulative viewshed;

